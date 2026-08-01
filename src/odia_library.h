// odia_library.h — the assay library as a compact object graph.
//
// WHAT THE CALLER SEES: proteins, peptides and transitions, each addressed by an opaque handle,
// with accessors that return ordinary `std::string_view`s. No offsets, no chunk indices, no spans,
// no interning. Handles are strongly typed, so a Peptide cannot be passed where a Protein belongs.
//
// WHAT IS UNDERNEATH, and why it is worth hiding: measured on the Astral benchmark library
// (2026-08-01, phase-resolved probe), materialising 7,149,966 compounds and 78,569,077 transitions
// as ordinary objects cost 38.79 GB of RSS, of which 11.29 GB was memory glibc had freed but could
// not return -- ~471 million small string allocations had fragmented the arena around the live
// ones. Three properties of proteomics data remove nearly all of it:
//
//   1. A tryptic peptide is a SUBSTRING of its protein. Store the ~20k FASTA sequences once
//      (~11 MB for the human proteome) and a peptide costs (protein, offset, length) -- no
//      characters of its own. Peptides that are not substrings (variants, semi-tryptic products,
//      anything synthesised) are stored explicitly; the caller cannot tell the difference and does
//      not need to.
//   2. Fragment annotations repeat massively -- "y7", "b3" and so on across millions of
//      transitions. Interned, they cost 4 bytes per reference.
//   3. `peptide_ref` was 78.6M references to ~7.1M distinct values. As a handle it is 4 bytes and
//      no allocation.
//
// The complexity is real, so it lives in one place with one test, behind an interface that reads
// like the domain. `stats()` exists so the saving is a measurement rather than a claim.

#ifndef ODIA_LIBRARY_H
#define ODIA_LIBRARY_H

#include "odia_seqstore.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace odia
{

class CompactLibrary
{
public:
  /// Opaque handles. `enum class` rather than a typedef so the compiler rejects
  /// `sequence(some_protein)` where a peptide was meant -- a mistake that would otherwise return a
  /// plausible wrong sequence rather than failing.
  enum class Protein : std::uint32_t {};
  enum class Peptide : std::uint32_t {};
  enum class Transition : std::uint32_t {};

  static constexpr Protein no_protein = Protein(~std::uint32_t(0));

  // ---- building ---------------------------------------------------------------------------
  /// @param accession e.g. "sp|P02768|ALBU_HUMAN"
  /// @param sequence  the FASTA amino-acid sequence; may be empty if unknown, in which case its
  ///                  peptides simply store themselves.
  Protein addProtein(std::string_view accession, std::string_view sequence = {})
  {
    ProteinRec r;
    r.accession = store_.intern(accession);
    r.sequence = sequence.empty() ? SequenceStore::npos : store_.append(sequence);
    proteins_.push_back(r);
    return Protein(std::uint32_t(proteins_.size() - 1));
  }

  /// Add a peptide belonging to `parent`. Whether it ends up as a substring reference or as its own
  /// stored sequence is decided here and never surfaces to the caller.
  Peptide addPeptide(Protein parent, std::string_view sequence, int charge = 0)
  {
    PeptideRec r;
    r.parent = parent;
    r.charge = std::int16_t(charge);
    const SequenceStore::Id parent_seq =
      (parent == no_protein) ? SequenceStore::npos : proteins_[idx(parent)].sequence;
    bool derived = false;
    r.span = store_.spanOfOrIntern(parent_seq, sequence, &derived);
    derived_ += derived ? 1 : 0;
    peptides_.push_back(r);
    return Peptide(std::uint32_t(peptides_.size() - 1));
  }

  /// A fragment of `p`. `annotation` ("y7", "b3+1") is interned: it repeats across the library.
  Transition addTransition(Peptide p, double product_mz, float intensity,
                           std::string_view annotation = {})
  {
    TransitionRec r;
    r.peptide = p;
    r.product_mz = product_mz;
    r.intensity = intensity;
    r.annotation = annotation.empty() ? SequenceStore::npos : store_.intern(annotation);
    transitions_.push_back(r);
    return Transition(std::uint32_t(transitions_.size() - 1));
  }

  /// Pre-size everything. Parquet's footer gives exact row counts and total column bytes, so a
  /// caller loading from parquet never has to guess -- and the vectors then never reallocate, which
  /// on the benchmark library is ~24 GB of avoidable memmove and ~6 GB of transient peak.
  void reserve(std::size_t n_proteins, std::size_t n_peptides, std::size_t n_transitions,
               std::size_t total_sequence_chars)
  {
    proteins_.reserve(n_proteins);
    peptides_.reserve(n_peptides);
    transitions_.reserve(n_transitions);
    store_.reserve(n_proteins + n_peptides / 8, total_sequence_chars);
  }

  /// Call when loading is finished: releases the interning index, which exists only for building.
  void finalize() { store_.freezeIndex(); }

  // ---- reading ----------------------------------------------------------------------------
  std::string_view accession(Protein p) const { return store_.view(proteins_[idx(p)].accession); }
  std::string_view proteinSequence(Protein p) const
  {
    const auto id = proteins_[idx(p)].sequence;
    return id == SequenceStore::npos ? std::string_view{} : store_.view(id);
  }
  std::string_view sequence(Peptide p) const { return store_.view(peptides_[idx(p)].span); }
  Protein parent(Peptide p) const { return peptides_[idx(p)].parent; }
  int charge(Peptide p) const { return peptides_[idx(p)].charge; }

  Peptide peptideOf(Transition t) const { return transitions_[idx(t)].peptide; }
  double productMz(Transition t) const { return transitions_[idx(t)].product_mz; }
  float intensity(Transition t) const { return transitions_[idx(t)].intensity; }
  std::string_view annotation(Transition t) const
  {
    const auto id = transitions_[idx(t)].annotation;
    return id == SequenceStore::npos ? std::string_view{} : store_.view(id);
  }
  /// Convenience: the sequence a transition ultimately belongs to, without the caller walking there.
  std::string_view sequenceOf(Transition t) const { return sequence(peptideOf(t)); }

  std::size_t proteinCount() const { return proteins_.size(); }
  std::size_t peptideCount() const { return peptides_.size(); }
  std::size_t transitionCount() const { return transitions_.size(); }

  // ---- measurement ------------------------------------------------------------------------
  struct Stats
  {
    std::size_t proteins, peptides, transitions;
    std::size_t peptides_from_fasta;   ///< stored as a substring of their protein: zero characters
    std::size_t peptides_standalone;   ///< variants etc., stored explicitly
    std::size_t bytes;                 ///< what this actually holds
    std::size_t bytes_as_objects;      ///< the same data as std::string-bearing structs
  };

  Stats stats() const
  {
    Stats s{};
    s.proteins = proteins_.size();
    s.peptides = peptides_.size();
    s.transitions = transitions_.size();
    s.peptides_from_fasta = derived_;
    s.peptides_standalone = peptides_.size() - derived_;
    s.bytes = sizeof(*this) + store_.bytes()
            + proteins_.capacity() * sizeof(ProteinRec)
            + peptides_.capacity() * sizeof(PeptideRec)
            + transitions_.capacity() * sizeof(TransitionRec);
    // The comparison layout: OpenMS's LightCompound is 376 B and LightTransition 128 B, both
    // measured, and each owning std::strings whose heap blocks are counted per REFERENCE because
    // copying a std::string always allocates.
    std::size_t heap = 0;
    for (const auto& p : peptides_)
    {
      const std::size_t n = store_.view(p.span).size();
      if (n > 15) { heap += ((n + 1 + 15) / 16) * 16 + 16; }
    }
    s.bytes_as_objects = peptides_.size() * 376 + transitions_.size() * 128 + heap
                       + transitions_.size() * (32 + 32);   // peptide_ref + annotation per transition
    return s;
  }

private:
  struct ProteinRec { SequenceStore::Id accession, sequence; };
  struct PeptideRec { SequenceStore::Span span; Protein parent; std::int16_t charge; };
  struct TransitionRec { double product_mz; float intensity; Peptide peptide; SequenceStore::Id annotation; };

  static std::size_t idx(Protein p) { return std::size_t(static_cast<std::uint32_t>(p)); }
  static std::size_t idx(Peptide p) { return std::size_t(static_cast<std::uint32_t>(p)); }
  static std::size_t idx(Transition t) { return std::size_t(static_cast<std::uint32_t>(t)); }

  SequenceStore store_;
  std::vector<ProteinRec> proteins_;
  std::vector<PeptideRec> peptides_;
  std::vector<TransitionRec> transitions_;
  std::size_t derived_ = 0;
};

} // namespace odia

#endif // ODIA_LIBRARY_H

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
#include <limits>

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
  /// Opaque handle for an interned fragment annotation.
  using AnnotationId = SequenceStore::Id;

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

  /// Call once after all proteins are added and before any peptide. Indexes every position of
  /// every protein sequence so a peptide can be located in the PROTEOME rather than in whichever
  /// protein its accession happens to name.
  void indexProteins(std::size_t k = 5) { store_.buildIndex(k); }

  /// Fold everything interned since the last index call INTO the index -- decoys above all.
  ///
  /// Decoys are shuffled, so none is a substring of a protein and each is interned on first sight.
  /// But they are not unrelated to each other: OpenSWATH shuffles within a peptide, so decoys of
  /// homologous or repeated targets share subsequences, and a decoy seen once can host later ones.
  /// Leaving them out of the index guarantees every one of the ~3.5M pays for its own characters.
  ///
  /// Call periodically during a read, not per peptide: each call sorts the new tail and merges it,
  /// so the cost amortises over a batch instead of being paid per insertion.
  void indexInterned() { store_.indexAppended(); }

  std::size_t indexBytes() const { return store_.indexBytes(); }

  /// Add a peptide. Resolution order, and the order matters:
  ///   1. look the sequence up in the indexed proteome -- if it occurs anywhere, it costs
  ///      (protein, offset, length) and NO characters;
  ///   2. only if it occurs nowhere, intern it in the hashed unique-string set.
  ///
  /// Searching the proteome rather than the accession's protein is the point. An accession can be
  /// absent, renamed, versioned differently, or the peptide can be shared between proteins -- in
  /// every one of those cases an accession-keyed lookup fails and stores characters for a sequence
  /// that is right there in the FASTA. `parent` is retained for protein rollup, not for lookup.
  Peptide addPeptide(Protein parent, std::string_view sequence, int charge = 0)
  {
    PeptideRec r;
    r.parent = parent;
    r.charge = std::int16_t(charge);
    r.span = store_.findAnywhere(sequence);              // 1. the proteome
    if (r.span.valid()) { ++derived_; }
    else { r.span = store_.internAsSpan(sequence); }     // 2. hashed unique set, deduplicated
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


  // ---- parallel bulk fill ------------------------------------------------------------------
  // The row-by-row readers run at 1.0 core (measured: 4.8 s for 7.1M precursors, 20.6 s for 78.6M
  // transitions). The rows are independent, so the only obstacles are the shared vectors and the
  // interning map. Both are removed by pre-sizing and pre-interning, after which threads write to
  // disjoint indices and touch nothing shared.
  //
  // Contract, deliberately narrow: resize first, intern every annotation BEFORE the parallel
  // region, then fill each index exactly once. Filling an index twice or leaving one unset is a
  // caller error this cannot detect cheaply -- which is why it is a separate, explicitly named API
  // rather than an overload of addTransition().

  /// Pre-size for index-addressed filling. Parquet's footer gives the exact count.
  void resizeTransitions(std::size_t n) { transitions_.assign(n, TransitionRec{}); }
  void resizePeptides(std::size_t n) { peptides_.assign(n, PeptideRec{}); }

  /// Intern an annotation up front, OUTSIDE any parallel region, and reuse the handle inside.
  /// Fragment annotations are a tiny distinct set ("y7", "b3", ...) repeated across millions of
  /// transitions, so this is a few hundred calls, not 78 million.
  AnnotationId internAnnotation(std::string_view a)
  {
    return a.empty() ? SequenceStore::npos : store_.intern(a);
  }

  /// Thread-safe when every index is written exactly once by one thread.
  void setTransition(std::size_t i, Peptide p, double product_mz, float intensity, AnnotationId ann)
  {
    TransitionRec& r = transitions_[i];
    r.peptide = p;
    r.product_mz = product_mz;
    r.intensity = intensity;
    r.annotation = ann;
  }

  /// Peptides cannot be filled in parallel the same way: locating a peptide inside its protein
  /// mutates nothing, but a peptide that is NOT a substring must be appended to the store, which
  /// does. So the parallel form takes a pre-resolved span, and the caller decides how to obtain it
  /// -- see locateOrNull()/internPeptideSequence() below.
  void setPeptide(std::size_t i, Protein parent, SequenceStore::Span span, int charge, bool derived)
  {
    PeptideRec& r = peptides_[i];
    r.parent = parent;
    r.span = span;
    r.charge = std::int16_t(charge);
    if (derived) { ++derived_; }
  }

  /// Proteome-wide lookup, read-only and therefore safe from many threads at once. Returns an
  /// invalid span when the sequence occurs in no protein -- which the caller then interns.
  SequenceStore::Span locateOrNull(std::string_view seq) const { return store_.findAnywhere(seq); }

  /// Mutating: store a sequence that is not a substring. Call from ONE thread (or serially after
  /// a parallel locate pass has identified the misses).
  SequenceStore::Span internPeptideSequence(std::string_view seq) { return store_.internAsSpan(seq); }

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



  // ---- fields the materialised LightTargetedExperiment needs --------------------------------
  void setPrecursor(Peptide p, double mz, double rt, double drift)
  {
    auto& r = peptides_[idx(p)];
    r.precursor_mz = mz;
    r.rt = float(rt);
    r.drift_time = float(drift);
  }
  double precursorMz(Peptide p) const { return peptides_[idx(p)].precursor_mz; }
  double rt(Peptide p) const { return peptides_[idx(p)].rt; }
  double driftTime(Peptide p) const { return peptides_[idx(p)].drift_time; }

  enum Flag : std::uint8_t { Decoy = 1, Detecting = 2, Identifying = 4, Quantifying = 8 };
  void setTransitionFlags(std::size_t i, std::int8_t frag_charge, std::uint8_t flags)
  {
    transitions_[i].fragment_charge = frag_charge;
    transitions_[i].flags = flags;
  }
  std::int8_t fragmentCharge(Transition t) const { return transitions_[idx(t)].fragment_charge; }
  std::uint8_t transitionFlags(Transition t) const { return transitions_[idx(t)].flags; }

  // ---- original identifiers -----------------------------------------------------------------
  // The library's own ids ("DECOY_PEPTIDEK_2") are 16+ characters and therefore heap-allocated:
  // 78.6M of them is where the fragmentation comes from. Materialising a LightTargetedExperiment
  // with SHORT synthetic ids keeps every string inside libstdc++'s 15-char SSO buffer, so the
  // materialised library allocates nothing per row -- but the original id must survive for output
  // and for anything that pairs on it. It is stored here, once, interned.
  void setOriginalId(Peptide p, std::string_view id) { peptides_[idx(p)].original_id = store_.intern(id); }
  std::string_view originalId(Peptide p) const
  {
    const auto id = peptides_[idx(p)].original_id;
    return id == SequenceStore::npos ? std::string_view{} : store_.view(id);
  }
  void setDecoy(Peptide p, bool d) { peptides_[idx(p)].decoy = d; }
  bool isDecoy(Peptide p) const { return peptides_[idx(p)].decoy; }

  /// Synthetic id for the materialised view: "p<n>" for targets, "DECOY_p<n>" for decoys.
  /// Both stay within the SSO buffer (max "DECOY_p4294967295" is 17 -- so the index is emitted in
  /// base-36, bounding it at "DECOY_p1z141z3" = 14 characters for any uint32).
  static std::string syntheticId(std::uint32_t index, bool decoy)
  {
    static const char* D = "0123456789abcdefghijklmnopqrstuvwxyz";
    char buf[8];
    int n = 0;
    std::uint32_t v = index;
    do { buf[n++] = D[v % 36]; v /= 36; } while (v && n < 7);
    std::string out;
    out.reserve(15);
    if (decoy) { out += "DECOY_"; }
    out += 'p';
    while (n) { out += buf[--n]; }
    return out;                                  // <= 6+1+7 = 14 chars: always SSO
  }

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
  // 48 B. Everything a LightCompound needs EXCEPT the strings, which are a span and two handles.
  struct PeptideRec
  {
    SequenceStore::Span span;                              // sequence, as a slice of its protein
    double precursor_mz = 0.0;
    float rt = std::numeric_limits<float>::quiet_NaN();    // library RT
    float drift_time = -1.0f;
    Protein parent = no_protein;
    SequenceStore::Id original_id = SequenceStore::npos;
    std::int16_t charge = 0;
    bool decoy = false;
  };
  // 24 B. The flags are the four booleans LightTransition carries.
  struct TransitionRec
  {
    double product_mz = 0.0;
    float intensity = 0.0f;
    Peptide peptide = Peptide(0);
    SequenceStore::Id annotation = SequenceStore::npos;
    std::int8_t fragment_charge = 0;
    std::uint8_t flags = 0;      // bit0 decoy, bit1 detecting, bit2 identifying, bit3 quantifying
  };

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

// Self-check for odia_library.h. Standalone: no OpenMS, no framework.
#include "odia_library.h"
#include <cassert>
#include <cstdio>
#include <random>
#include <string>
#include <utility>

using odia::CompactLibrary;

int main()
{
  std::printf("odia_library_test\n");

  // T1 -- the interface reads like the domain; nothing about spans or interning is visible.
  {
    CompactLibrary lib;
    const auto alb = lib.addProtein("sp|P02768|ALBU_HUMAN", "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFK");
    lib.indexProteins(6);
    const auto pep = lib.addPeptide(alb, "DAHKSEVAHRFK", 2);
    const auto tr  = lib.addTransition(pep, 512.27, 0.83f, "y7");
    assert(lib.sequence(pep) == "DAHKSEVAHRFK");
    assert(lib.accession(alb) == "sp|P02768|ALBU_HUMAN");
    assert(lib.parent(pep) == alb);
    assert(lib.charge(pep) == 2);
    assert(lib.annotation(tr) == "y7");
    assert(lib.sequenceOf(tr) == "DAHKSEVAHRFK");
    std::printf("  T1 domain accessors return plain string_views: OK\n");
  }

  // T2 -- a peptide that IS a substring costs no characters; one that is not still works.
  {
    CompactLibrary lib;
    const auto p = lib.addProtein("P1", "MKWVTFISLLFLFSSAYSR");
    lib.indexProteins(6);
    lib.addPeptide(p, "FISLLFLFSSAYSR");          // substring
    lib.addPeptide(p, "FISLLFLFSSAYSK");          // variant: K for R, not present
    const auto s = lib.stats();
    assert(s.peptides_from_fasta == 1);
    assert(s.peptides_standalone == 1);
    assert(lib.sequence(CompactLibrary::Peptide(0)) == "FISLLFLFSSAYSR");
    assert(lib.sequence(CompactLibrary::Peptide(1)) == "FISLLFLFSSAYSK");   // must NOT silently
                                                                            // resolve to the R form
    std::printf("  T2 substring vs variant both correct, and distinguishable: OK\n");
  }

  // T3 -- a peptide with no parent protein sequence still round-trips.
  {
    CompactLibrary lib;
    const auto p = lib.addProtein("NOSEQ");        // accession only
    lib.indexProteins(6);
    const auto pep = lib.addPeptide(p, "SYNTHETICPEPTIDER");
    assert(lib.sequence(pep) == "SYNTHETICPEPTIDER");
    assert(lib.proteinSequence(p).empty());
    std::printf("  T3 protein without a sequence: peptides store themselves: OK\n");
  }

  // T4 -- finalize() drops the build-time index; reads must be unaffected.
  {
    CompactLibrary lib;
    const auto p = lib.addProtein("P", "AAAKBBBKCCCK");
    lib.indexProteins(4);
    const auto a = lib.addPeptide(p, "BBBK");
    lib.finalize();
    assert(lib.sequence(a) == "BBBK");
    std::printf("  T4 reads survive finalize(): OK\n");
  }

  // T5 -- the measurement that justifies the design, at proteome scale and realistic ratios.
  {
    CompactLibrary lib;
    std::mt19937 rng(7);
    const int n_prot = 2000, pep_per_prot = 40, tr_per_pep = 6;
    const char* aa = "ACDEFGHIKLMNPQRSTVWY";
    std::vector<std::pair<CompactLibrary::Protein, std::string>> prots;
    for (int i = 0; i < n_prot; ++i)
    {
      std::string seq;
      seq.reserve(600);
      for (int j = 0; j < 600; ++j) { seq += aa[rng() % 20]; }
      prots.push_back({lib.addProtein("sp|ACC" + std::to_string(i) + "|X_HUMAN", seq), seq});
    }
    lib.indexProteins();                       // index the whole proteome ONCE, then add peptides
    for (auto& pp : prots)
    {
      const auto pr = pp.first; const std::string& seq = pp.second;
      for (int k = 0; k < pep_per_prot; ++k)
      {
        const std::size_t off = rng() % 570;
        const std::size_t len = 12 + rng() % 15;
        const auto pep = lib.addPeptide(pr, std::string_view(seq).substr(off, len), 2);
        for (int t = 0; t < tr_per_pep; ++t)
        {
          lib.addTransition(pep, 400.0 + t, 1.0f, (t % 2 ? "y" : "b") + std::to_string(t + 3));
        }
      }
    }
    lib.finalize();
    const auto s = lib.stats();
    std::printf("  T5 %zu proteins, %zu peptides (%zu from FASTA, %zu standalone), %zu transitions\n",
                s.proteins, s.peptides, s.peptides_from_fasta, s.peptides_standalone, s.transitions);
    std::printf("     compact %.1f MB  vs  string-bearing objects %.1f MB  -> %.1fx smaller\n",
                s.bytes / 1048576.0, s.bytes_as_objects / 1048576.0,
                double(s.bytes_as_objects) / double(s.bytes));
    assert(s.bytes < s.bytes_as_objects);
    assert(s.peptides_from_fasta > s.peptides * 9 / 10);   // nearly all should be substrings
  }

  // T6 -- synthetic ids must ALWAYS fit libstdc++'s 15-char SSO buffer, or the whole point
  // (materialising a library with zero per-row heap allocations) is lost. Check the extremes.
  {
    for (std::uint32_t v : {0u, 1u, 999u, 3603424u, 78569076u, 4294967295u})
    {
      const auto t = CompactLibrary::syntheticId(v, false);
      const auto d = CompactLibrary::syntheticId(v, true);
      assert(t.size() <= 15 && d.size() <= 15);
      assert(d.rfind("DECOY_", 0) == 0);
      assert(d.substr(6) == t);              // pairing is by DECOY_ + target id
    }
    std::printf("  T6 synthetic ids <=15 chars (SSO) incl. UINT32_MAX, decoy = DECOY_+target: OK\n");
  }

  // T7 -- original ids survive alongside the synthetic ones.
  {
    CompactLibrary lib;
    const auto pr = lib.addProtein("P", "AAAKBBBKCCCK");
    lib.indexProteins(4);
    const auto a = lib.addPeptide(pr, "BBBK", 2);
    lib.setOriginalId(a, "DECOY_PEPTIDEK_2");
    lib.setDecoy(a, true);
    assert(lib.originalId(a) == "DECOY_PEPTIDEK_2");
    assert(lib.isDecoy(a));
    assert(lib.sequence(a) == "BBBK");
    std::printf("  T7 original id + decoy flag retained: OK\n");
  }

  // T8 -- a peptide SHORTER than the k-mer must still be found. Returning "not found" for it is a
  // silent false negative: on the real library it sent 127,853 peptides to standalone storage.
  {
    CompactLibrary lib;
    lib.addProtein("P", "MKWVTFISLLFLFSSAYSRGVFRR");
    lib.indexProteins(8);                       // deliberately longer than the peptides below
    const auto short5 = lib.addPeptide(CompactLibrary::no_protein, "GVFRR");   // 5 < k
    const auto short3 = lib.addPeptide(CompactLibrary::no_protein, "SSA");     // 3 < k
    const auto absent = lib.addPeptide(CompactLibrary::no_protein, "WWWW");    // truly absent
    assert(lib.sequence(short5) == "GVFRR");
    assert(lib.sequence(short3) == "SSA");
    assert(lib.sequence(absent) == "WWWW");
    const auto st = lib.stats();
    assert(st.peptides_from_fasta == 2);        // both short ones found despite being under k
    assert(st.peptides_standalone == 1);        // only the genuinely absent one interned
    std::printf("  T8 needles shorter than k are found, not silently missed: OK\n");
  }

  // T9 -- decoys must become searchable once indexed, so a later decoy sharing a subsequence with
  // an earlier one costs a span rather than its own characters.
  {
    CompactLibrary lib;
    lib.addProtein("P", "MKWVTFISLLFLFSSAYSR");
    lib.indexProteins(5);
    // A shuffled decoy: not in the protein, so it is interned.
    const auto d1 = lib.addPeptide(CompactLibrary::no_protein, "RSYASSFLFLLSIFTVWKM");
    auto st1 = lib.stats();
    assert(st1.peptides_standalone == 1);
    lib.indexInterned();                       // fold the decoy into the index
    // A second decoy that is a SUBSTRING of the first must now be found, not re-stored.
    const auto d2 = lib.addPeptide(CompactLibrary::no_protein, "ASSFLFLLSIF");
    auto st2 = lib.stats();
    assert(lib.sequence(d1) == "RSYASSFLFLLSIFTVWKM");
    assert(lib.sequence(d2) == "ASSFLFLLSIF");
    assert(st2.peptides_from_fasta == st1.peptides_from_fasta + 1);   // resolved against the decoy
    assert(st2.peptides_standalone == 1);                            // nothing new interned
    std::printf("  T9 interned decoys become searchable after indexInterned(): OK\n");
  }

  std::printf("odia_library_test OK\n");
  return 0;
}

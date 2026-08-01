// Self-check for odia_library.h. Standalone: no OpenMS, no framework.
#include "odia_library.h"
#include <cassert>
#include <cstdio>
#include <random>
#include <string>

using odia::CompactLibrary;

int main()
{
  std::printf("odia_library_test\n");

  // T1 -- the interface reads like the domain; nothing about spans or interning is visible.
  {
    CompactLibrary lib;
    const auto alb = lib.addProtein("sp|P02768|ALBU_HUMAN", "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFK");
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
    const auto pep = lib.addPeptide(p, "SYNTHETICPEPTIDER");
    assert(lib.sequence(pep) == "SYNTHETICPEPTIDER");
    assert(lib.proteinSequence(p).empty());
    std::printf("  T3 protein without a sequence: peptides store themselves: OK\n");
  }

  // T4 -- finalize() drops the build-time index; reads must be unaffected.
  {
    CompactLibrary lib;
    const auto p = lib.addProtein("P", "AAAKBBBKCCCK");
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
    for (int i = 0; i < n_prot; ++i)
    {
      std::string seq;
      seq.reserve(600);
      for (int j = 0; j < 600; ++j) { seq += aa[rng() % 20]; }
      const auto pr = lib.addProtein("sp|ACC" + std::to_string(i) + "|X_HUMAN", seq);
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

  std::printf("odia_library_test OK\n");
  return 0;
}

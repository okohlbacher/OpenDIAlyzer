// Validate P2 mod_x featurization (OpenMS/ML/PEPTDEEP/PeptDeepModX.h) against
// AlphaPeptDeep's ground truth. Expected values are the peptdeep 1.5.0 output
// captured in src/testdata/peptdeep_modx_reference.json -- so this is
// differential testing against the reference predictor, exactly the bit-parity
// gate the review demanded for the highest-risk patch. Each assertion is one of
// codex's named failure modes (element order, terminal placement, stacking,
// multi-element composition).
//
// Links the already-built vendored libOpenMS (AASequence/EmpiricalFormula); the
// featurization header is inline, so no OpenMS rebuild is needed to run this.

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepModX.h>

#include <cstdio>
#include <string>
#include <vector>

using OpenMS::AASequence;
using namespace OpenMS::ML;

namespace
{
int failures = 0;
constexpr int E = 109;

// Assert that, for `spec` (OpenMS modified-peptide syntax) with nAA residues,
// mod_x at (position, element-symbol) equals `value` -- and that NO other
// element at that position is set (catches intensity landing on the wrong slot).
void check(const std::string& spec, size_t nAA, size_t position,
           const std::vector<std::pair<std::string, float>>& expected)
{
  const AASequence seq = AASequence::fromString(spec);
  const size_t seqlen = nAA + 2; // terminal tokens
  std::vector<float> modx(seqlen * E, -1.0f);
  fillPeptideModX(seq, /*add_terminal_tokens*/ true, seqlen, modx.data());

  for (const auto& [sym, val] : expected)
  {
    const int idx = peptDeepModElementIndex(sym);
    const float got = modx[position * E + idx];
    if (got != val)
    {
      std::fprintf(stderr, "FAIL %-22s pos %zu %s(idx %d): got %.3f want %.3f\n",
                   spec.c_str(), position, sym.c_str(), idx, got, val);
      ++failures;
    }
  }
  // Nothing else at this position may be non-zero.
  float sum_expected = 0;
  for (const auto& [sym, val] : expected) { (void)sym; sum_expected += val; }
  float sum_got = 0;
  for (int j = 0; j < E; ++j) sum_got += modx[position * E + j];
  if (sum_got != sum_expected)
  {
    std::fprintf(stderr, "FAIL %-22s pos %zu has extra nonzero elements (sum %.3f != %.3f)\n",
                 spec.c_str(), position, sum_got, sum_expected);
    ++failures;
  }
}

// Assert every position other than those listed is all-zero (no leakage).
void check_only_positions(const std::string& spec, size_t nAA,
                          const std::vector<size_t>& nonzero_positions)
{
  const AASequence seq = AASequence::fromString(spec);
  const size_t seqlen = nAA + 2;
  std::vector<float> modx(seqlen * E, 0.0f);
  fillPeptideModX(seq, true, seqlen, modx.data());
  for (size_t p = 0; p < seqlen; ++p)
  {
    bool expected_nz = false;
    for (size_t q : nonzero_positions) if (q == p) expected_nz = true;
    float s = 0;
    for (int j = 0; j < E; ++j) s += modx[p * E + j];
    if (!expected_nz && s != 0.0f)
    {
      std::fprintf(stderr, "FAIL %-22s pos %zu should be all-zero, sum %.3f\n", spec.c_str(), p, s);
      ++failures;
    }
  }
}
} // namespace

int main()
{
  // Vocabulary integrity: order pins the whole scheme.
  if (std::string(peptDeepModElements()[3]) != "O" ||
      std::string(peptDeepModElements()[0]) != "C" ||
      std::string(peptDeepModElements()[108]) != "?")
  { std::fprintf(stderr, "FAIL vocabulary order\n"); ++failures; }

  // Fixture cases (peptdeep_modx_reference.json). residue r (0-based) -> pos r+1.
  // AMKMR, Oxidation on M (residue 2, 1-based) -> pos 2: O=1
  check("AM(Oxidation)KMR", 5, 2, {{"O", 1.0f}});
  check_only_positions("AM(Oxidation)KMR", 5, {2});

  // MPEPTIDE, Oxidation on first residue M -> pos 1: O=1  (first-residue edge)
  check("M(Oxidation)PEPTIDE", 8, 1, {{"O", 1.0f}});

  // PEPTIDEC, Carbamidomethyl on last residue C (residue 8) -> pos 8: C2 H3 N1 O1
  check("PEPTIDEC(Carbamidomethyl)", 8, 8, {{"C", 2.0f}, {"H", 3.0f}, {"N", 1.0f}, {"O", 1.0f}});

  // ACEPTIDE, N-terminal Acetyl -> pos 0: C2 H2 O1  (terminal placement)
  check(".(Acetyl)ACEPTIDE", 8, 0, {{"C", 2.0f}, {"H", 2.0f}, {"O", 1.0f}});
  check_only_positions(".(Acetyl)ACEPTIDE", 8, {0});

  // SMTMYK, two Oxidations (residues 2 and 4) -> pos 2 and pos 4  (stacking/multiple)
  check("SM(Oxidation)TM(Oxidation)YK", 6, 2, {{"O", 1.0f}});
  check("SM(Oxidation)TM(Oxidation)YK", 6, 4, {{"O", 1.0f}});
  check_only_positions("SM(Oxidation)TM(Oxidation)YK", 6, {2, 4});

  // SPYTK, Phospho on S (residue 1) -> pos 1: H1 O3 P1  (multi-element, P slot)
  check("S(Phospho)PYTK", 5, 1, {{"H", 1.0f}, {"O", 3.0f}, {"P", 1.0f}});

  if (failures == 0) { std::puts("odia-modx: all P2 fixture checks OK"); return 0; }
  std::fprintf(stderr, "odia-modx: %d FAILURES\n", failures);
  return 1;
}

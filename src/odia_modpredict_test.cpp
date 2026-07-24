// Integration gate for P2: prove mod_x flows through inference and changes the
// prediction. The featurization is already bit-parity validated (odia-modx-test)
// and the unmodified ONNX path is 1e-7 vs AlphaPeptDeep (PeptDeepInference_test);
// this checks the wiring between them -- that predicting a MODIFIED peptide
// differs from its unmodified base, in a way only a populated mod_x could cause.

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace OpenMS;

int main(int argc, char** argv)
{
  if (argc < 2) { std::fprintf(stderr, "usage: odia-modpredict-test <model_dir>\n"); return 2; }
  const std::string dir = argv[1];
  int failures = 0;

  PeptDeepMS2Inference ms2(dir + "/peptdeep_ms2_dynamic.onnx", 1, 64);
  PeptDeepRTInference rt(dir + "/peptdeep_rt_dynamic.onnx", 1, 64);

  // Same base sequence, with and without Oxidation on the M.
  const std::vector<AASequence> unmod{AASequence::fromString("AMKMR")};
  const std::vector<AASequence> mod{AASequence::fromString("AM(Oxidation)KMR")};
  const std::vector<float> ch{2.0f}, nce{30.0f};
  const std::vector<int64_t> instr{0};

  auto a = ms2.predictMS2(unmod, ch, nce, instr);
  auto b = ms2.predictMS2(mod, ch, nce, instr);

  // Same shape (same base sequence length), but the intensities must differ --
  // a zero mod_x (the old behaviour) would make them identical.
  if (a.size() != 1 || b.size() != 1 || a[0].size() != b[0].size())
  { std::fprintf(stderr, "FAIL shape mismatch\n"); ++failures; }
  else
  {
    double maxdiff = 0.0, sumabs = 0.0;
    for (size_t i = 0; i < a[0].size(); ++i)
    { const double d = std::fabs(a[0][i] - b[0][i]); maxdiff = std::max(maxdiff, d); sumabs += d; }
    std::printf("MS2 unmod-vs-Ox(M): max|d|=%.4g  sum|d|=%.4g  (n=%zu)\n", maxdiff, sumabs, a[0].size());
    if (maxdiff < 1e-4) { std::fprintf(stderr, "FAIL: modified MS2 identical to unmodified (mod_x not applied)\n"); ++failures; }
  }

  // RT: oxidation shifts hydrophobicity, so predicted iRT should move too.
  auto ru = rt.predictRT(unmod), rm = rt.predictRT(mod);
  std::printf("RT unmod=%.4f  Ox(M)=%.4f  |d|=%.4g\n", ru[0], rm[0], std::fabs(ru[0] - rm[0]));
  if (std::fabs(ru[0] - rm[0]) < 1e-5) { std::fprintf(stderr, "FAIL: modified RT identical to unmodified\n"); ++failures; }

  // Sanity: the unmodified AASequence path must equal the string path exactly
  // (delegation must be transparent).
  auto s = ms2.predictMS2(std::vector<std::string>{"AMKMR"}, ch, nce, instr);
  for (size_t i = 0; i < a[0].size(); ++i)
    if (a[0][i] != s[0][i]) { std::fprintf(stderr, "FAIL: AASequence path != string path for unmodified\n"); ++failures; break; }

  if (failures == 0) { std::puts("odia-modpredict: P2 inference integration OK"); return 0; }
  std::fprintf(stderr, "odia-modpredict: %d FAILURES\n", failures);
  return 1;
}

// odia_chromatogram_test.cpp — self-check for the compact XIC store.
//
// The point of this file is that every claim odia_chromatogram.h makes is MEASURED here, not
// asserted in a comment:
//
//   T1 round-trip fidelity: RTs are exact (they are the shared grid, not a copy), intensities are
//      within float32's relative epsilon of the doubles that went in.
//   T2 the float32 narrowing does not move any score that matters -- Pearson correlation and dot
//      product between two realistic XICs, computed from float32 vs from double, agree to ~1e-7.
//   T3 the memory claim, on the real benchmark shape: 4x, and the exact GB figure.
//   T4 the shared-grid invariant holds for chromatograms with DIFFERENT sub-ranges of one grid --
//      this is the assumption the design rests on, so a sub-range must map to the right RTs.
//   T5 out-of-range spans THROW rather than clamp (a clamped span pairs intensities with the wrong
//      retention times, which is a wrong answer, not a smaller one).
//   T6 degenerate shapes: empty chromatogram, single point, full-grid.
//
// Build: c++ -std=c++17 -O2 src/odia_chromatogram_test.cpp -o t && ./t

#include "odia_chromatogram.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

static int failures = 0;
static void check(bool ok, const char* what)
{
  if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++failures; }
}

/// A realistic XIC: gaussian peak on a noisy baseline, intensities spanning ~5 decades.
static std::vector<double> makeXIC(std::mt19937& rng, std::size_t n, double apex, double sigma,
                                   double height)
{
  std::uniform_real_distribution<double> noise(0.0, height * 1e-3);
  std::vector<double> v(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    const double d = (double(i) - apex) / sigma;
    v[i] = height * std::exp(-0.5 * d * d) + noise(rng);
  }
  return v;
}

static double pearson(const std::vector<double>& a, const std::vector<double>& b)
{
  const std::size_t n = a.size();
  double ma = 0, mb = 0;
  for (std::size_t i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  double num = 0, da = 0, db = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    num += (a[i] - ma) * (b[i] - mb);
    da += (a[i] - ma) * (a[i] - ma);
    db += (b[i] - mb) * (b[i] - mb);
  }
  return num / std::sqrt(da * db);
}

int main()
{
  std::mt19937 rng(1234);

  // A realistic grid: 2333 s gradient, 0.6 s cycle -> 3889 points, as on the Astral benchmark.
  const std::size_t GRID = 3889;
  std::vector<double> grid(GRID);
  for (std::size_t i = 0; i < GRID; ++i) { grid[i] = 2333.4 * double(i) / double(GRID - 1); }

  // ---- T1 round-trip ------------------------------------------------------------------
  {
    odia::ChromatogramStore s;
    s.setGrid(grid);
    const std::size_t first = 500, n = 2392;
    auto xic = makeXIC(rng, n, 1200.0, 40.0, 5.0e8);
    s.add(/*transition*/ 42, first, xic);

    check(s.size() == 1 && s.points(0) == n, "T1 one chromatogram of the right length");
    check(s.entry(0).transition_index == 42, "T1 transition index preserved");

    bool rt_exact = true, int_ok = true;
    double worst_rel = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
      if (s.rt(0, j) != grid[first + j]) { rt_exact = false; }
      const double rel = std::abs(double(s.intensity(0, j)) - xic[j]) / std::abs(xic[j]);
      worst_rel = std::max(worst_rel, rel);
      if (rel > 1e-6) { int_ok = false; }
    }
    check(rt_exact, "T1 RTs are EXACT (read from the shared grid, never copied or narrowed)");
    check(int_ok, "T1 intensities within float32 relative epsilon");
    std::fprintf(stderr, "  T1 worst relative intensity error: %.3e (float32 eps = %.3e)\n",
                 worst_rel, double(std::numeric_limits<float>::epsilon()));
  }

  // ---- T2 the narrowing does not move the scores --------------------------------------
  // This is the claim that actually matters: nobody consumes an intensity, they consume a
  // correlation or a dot product between two XICs. Compute both ways and compare.
  {
    const std::size_t n = 2392;
    auto a = makeXIC(rng, n, 1200.0, 40.0, 5.0e8);
    auto b = makeXIC(rng, n, 1203.0, 42.0, 1.7e7);   // different apex and 30x lower -- a real pair
    std::vector<double> af(n), bf(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      af[i] = double(float(a[i]));
      bf[i] = double(float(b[i]));
    }
    const double r_d = pearson(a, b), r_f = pearson(af, bf);
    double dot_d = 0, dot_f = 0;
    for (std::size_t i = 0; i < n; ++i) { dot_d += a[i] * b[i]; dot_f += af[i] * bf[i]; }
    const double dr = std::abs(r_d - r_f);
    const double dd = std::abs(dot_d - dot_f) / std::abs(dot_d);
    std::fprintf(stderr, "  T2 correlation %.12f (double) vs %.12f (float32), |delta| = %.3e\n",
                 r_d, r_f, dr);
    std::fprintf(stderr, "  T2 dot product relative delta = %.3e\n", dd);
    check(dr < 1e-6, "T2 float32 does not move the XIC correlation");
    check(dd < 1e-6, "T2 float32 does not move the XIC dot product");
  }

  // ---- T3 the memory claim, on the real benchmark shape --------------------------------
  {
    odia::ChromatogramStore s;
    s.setGrid(grid);
    const std::size_t N = 20000, n = 2392;          // sample; scale the answer to the full run
    s.reserve(N, N * n);
    std::vector<float> flat(n, 1.0f);
    for (std::size_t i = 0; i < N; ++i) { s.add(std::uint32_t(i), 500, flat.data(), n); }

    const double got = double(s.bytes());
    const double ref = double(s.bytesAsMSChromatogram());
    const double ratio = ref / got;
    std::fprintf(stderr, "  T3 %zu chromatograms x %zu points:\n", N, n);
    std::fprintf(stderr, "       MSChromatogram %.2f GB   compact %.2f GB   ratio %.2fx\n",
                 ref / 1073741824.0, got / 1073741824.0, ratio);
    // Scale to the measured benchmark: 4,463,919 MS2 + 423,079*4 MS1 chromatograms
    const double scale = (4463919.0 + 423079.0 * 4.0) / double(N);
    std::fprintf(stderr, "       full run (6,156,235 chrom @ 1435 s window): %.0f GB -> %.0f GB, "
                         "saves %.0f GB\n",
                 ref * scale / 1073741824.0, got * scale / 1073741824.0,
                 (ref - got) * scale / 1073741824.0);
    check(ratio > 3.9, "T3 compact store is >3.9x smaller than vector<MSChromatogram>");
    // Per-point cost must be exactly 4 bytes; anything else means padding crept in.
    const double per_point = (got - double(sizeof(odia::ChromatogramStore)) -
                              double(grid.capacity() * sizeof(double)) -
                              double(N * sizeof(odia::ChromatogramStore::Entry))) / double(N * n);
    std::fprintf(stderr, "  T3 marginal bytes per point: %.3f\n", per_point);
    check(per_point > 3.99 && per_point < 4.01, "T3 exactly 4 bytes per point (no per-point RT)");
    check(sizeof(odia::ChromatogramStore::Entry) == 24, "T3 Entry is 24 bytes");
  }

  // ---- T4 different sub-ranges of the SAME grid ----------------------------------------
  // The design rests on chromatograms occupying different contiguous slices of one grid, which is
  // what the extractor produces when precursors have different RT windows. A sub-range must read
  // back the RIGHT retention times, not the grid's leading points.
  {
    odia::ChromatogramStore s;
    s.setGrid(grid);
    struct Case { std::uint32_t first, n; };
    const Case cases[] = {{0, 100}, {1500, 800}, {3000, 889}, {2000, 1}};
    std::vector<float> buf(1000, 2.0f);
    for (const Case& c : cases) { s.add(0, c.first, buf.data(), c.n); }
    bool all_ok = true;
    for (std::size_t i = 0; i < 4; ++i)
    {
      std::vector<double> expect(grid.begin() + cases[i].first,
                                 grid.begin() + cases[i].first + cases[i].n);
      if (!s.gridMatches(i, expect)) { all_ok = false; }
      if (s.rt(i, 0) != grid[cases[i].first]) { all_ok = false; }
    }
    check(all_ok, "T4 every sub-range reads back its own RTs, exactly");
    check(s.gridIndexAtOrAfter(grid[1500]) == 1500, "T4 gridIndexAtOrAfter finds an exact point");
    check(s.gridIndexAtOrAfter(-1.0) == 0, "T4 gridIndexAtOrAfter below the grid -> 0");
    check(s.gridIndexAtOrAfter(1e9) == GRID, "T4 gridIndexAtOrAfter above the grid -> size");
  }

  // ---- T5 out-of-range spans must THROW ------------------------------------------------
  {
    odia::ChromatogramStore s;
    s.setGrid(grid);
    std::vector<float> buf(100, 1.0f);
    bool threw = false;
    try { s.add(0, std::uint32_t(GRID - 10), buf.data(), 100); } catch (const std::out_of_range&) { threw = true; }
    check(threw, "T5 a span running off the end of the grid throws (never silently clamps)");

    bool threw2 = false;
    odia::ChromatogramStore s2;
    try { s2.setGrid({5.0, 4.0, 3.0}); } catch (const std::invalid_argument&) { threw2 = true; }
    check(threw2, "T5 an unsorted grid is rejected (would mean mixed SWATH maps)");

    bool threw3 = false;
    odia::ChromatogramStore s3;
    s3.setGrid(grid);
    s3.add(0, 0, buf.data(), 100);
    try { s3.setGrid(grid); } catch (const std::logic_error&) { threw3 = true; }
    check(threw3, "T5 changing the grid after adding data throws");
  }

  // ---- T6 degenerate shapes -------------------------------------------------------------
  {
    odia::ChromatogramStore s;
    s.setGrid(grid);
    std::vector<float> one(1, 7.0f);
    s.add(3, 0, nullptr, 0);                       // empty
    s.add(4, 100, one.data(), 1);                  // single point
    std::vector<float> full(GRID, 1.0f);
    s.add(5, 0, full.data(), GRID);                // whole grid
    check(s.points(0) == 0, "T6 empty chromatogram stored");
    check(s.points(1) == 1 && s.rt(1, 0) == grid[100] && s.intensity(1, 0) == 7.0f,
          "T6 single-point chromatogram reads back correctly");
    check(s.points(2) == GRID && s.rt(2, GRID - 1) == grid[GRID - 1],
          "T6 full-grid chromatogram reads back correctly");
    bool threw = false;
    try { (void)s.rt(0, 0); } catch (const std::out_of_range&) { threw = true; }
    check(threw, "T6 indexing into an empty chromatogram throws");
  }

  if (failures) { std::fprintf(stderr, "odia_chromatogram_test FAILED (%d)\n", failures); return 1; }
  std::fprintf(stderr, "odia_chromatogram_test OK\n");
  return 0;
}

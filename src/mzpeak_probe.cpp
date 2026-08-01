// M0 gate: prove the mzPeak reader streams a real diaPASEF run reading metadata
// ONLY (RT / isolation window / ion mobility) without decoding peaks -- the
// bounded-memory access the streaming extractor is built on. Also closes the
// open mzPeak handoff item (ion-mobility values were UNVERIFIED on real IM data).
//
//   mzpeak_probe <file.mzpeak>
#include <mzpeak/open.h>
#include <mzpeak/index.h>
#include <mzpeak/spectra.h>
#include <mzpeak/spectrum.h>
#include <mzpeak/spectrum_metadata.h>

#include <sys/resource.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <optional>

static long rss_kb() { struct rusage r; getrusage(RUSAGE_SELF, &r); return r.ru_maxrss; }

int main(int argc, char** argv)
{
  if (argc < 2) { std::fprintf(stderr, "usage: mzpeak_probe <file.mzpeak>\n"); return 2; }
  auto index = MzPeak::open(argv[1]);
  auto spectra = index.spectra();

  size_t n = 0, ms1 = 0, ms2 = 0, im_nonnull = 0, iso_present = 0;
  double rt_min = 1e30, rt_max = -1e30;
  std::map<long, size_t> windows;   // rounded isolation-window centre -> count
  double im_min = 1e30, im_max = -1e30;

  for (auto sp : spectra)            // lazy per-spectrum; we never call mz()
  {
    ++n;
    if (auto rt = sp.retention_time()) { rt_min = std::min(rt_min, *rt); rt_max = std::max(rt_max, *rt); }
    const int level = sp.ms_level();
    if (level == 1) ++ms1;
    else if (level == 2)
    {
      ++ms2;
      if (auto im = sp.ion_mobility()) { ++im_nonnull; im_min = std::min(im_min, *im); im_max = std::max(im_max, *im); }
      for (const auto& p : sp.precursors())
        if (p.isolation_window.target_mz)
        { ++iso_present; windows[std::lround(*p.isolation_window.target_mz)]++; }
    }
    if (n % 5000 == 0)
      std::fprintf(stderr, "  ...%zu spectra, RSS=%ld MB\r", n, rss_kb() / 1024);
  }

  std::printf("\nspectra          %zu  (MS1 %zu, MS2 %zu)\n", n, ms1, ms2);
  std::printf("RT range (s)     %.2f .. %.2f\n", rt_min, rt_max);
  std::printf("isolation windows %zu distinct centres (diaPASEF expects ~24-25)\n", windows.size());
  std::printf("MS2 with isolation %zu / %zu\n", iso_present, ms2);
  std::printf("MS2 with 1/K0     %zu / %zu  (%.1f%%)  range %.3f .. %.3f\n",
              im_nonnull, ms2, ms2 ? 100.0 * im_nonnull / ms2 : 0.0,
              im_nonnull ? im_min : 0.0, im_nonnull ? im_max : 0.0);
  std::printf("peak RSS         %ld MB   <- metadata-only pass, should be small\n", rss_kb() / 1024);

  // acceptance: diaPASEF window scheme + IM present + RT in seconds
  const bool ok = windows.size() >= 20 && windows.size() <= 40
               && im_nonnull > ms2 * 0.9 && rt_max > 60.0;
  std::printf("%s\n", ok ? "M0 OK: streaming metadata read verified" : "M0 CHECK: see numbers above");
  return ok ? 0 : 1;
}

// Degenerate-input check for the extraction-split derivation. Mirrors the shipped arithmetic.
#include "odia_split.h"

#include <cstdio>
static bool check(int threads, std::size_t n_ms2, std::size_t compounds, const char* label)
{
  const odia::ExtractionSplit s = odia::deriveExtractionSplit(threads, n_ms2, compounds);
  const int outer = s.outer_threads, inner = s.inner_threads, batch = s.batch_size;
  const bool ok = outer >= 1 && outer <= static_cast<int>(n_ms2) && inner >= 1 && batch >= 50;
  std::printf("  %-28s threads=%-4d windows=%-5zu -> outer=%-4d inner=%-4d batch=%-6d %s\n",
              label, threads, n_ms2, outer, inner, batch, ok ? "ok" : "FAIL");
  return ok;
}
int main()
{
  std::printf("extraction-split derivation, degenerate inputs\n");
  bool all = true;
  all &= check(224, 150, 423079, "the benchmark");
  all &= check(1,   150, 423079, "single thread");
  all &= check(8,   150, 423079, "few threads");
  all &= check(224, 1,   423079, "ONE window");
  all &= check(224, 150, 0,      "zero compounds");
  all &= check(224, 500, 7149966,"large library, many windows");
  all &= check(13,  150, 423079, "threads < 14");
  all &= check(448, 2,   100,    "more threads than work");
  std::printf("%s\n", all ? "derivetest OK" : "derivetest FAILED");
  return all ? 0 : 1;
}

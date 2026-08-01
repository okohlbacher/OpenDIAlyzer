// odia_lda_realdata_test.cpp — regression guard: the LDA must NOT degenerate to constant
// output on REAL extracted features. codex's FDR-honesty rewrite passed the synthetic
// adversarial suite but produced an ALL-EQUAL d-score (~1.16) on real OpenSWATH sub-scores
// (-> 0 IDs), while a reference sklearn LDA on the same features separates the classes.
// Fixture (testdata/lda_fixture.txt): first line "nrows ncols", then "group label f0..fN"
// (label 1=target, 0=decoy; missing = "nan"). ~44.5k real rows / 2500 precursor groups.
//
// Build: c++ -std=c++17 -O2 src/odia_lda_realdata_test.cpp -o t && ./t

#include "odia_lda.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  const char* path = argc > 1 ? argv[1] : "testdata/lda_fixture.txt";
  std::ifstream in(path);
  if (!in) { std::fprintf(stderr, "cannot open %s\n", path); return 2; }
  size_t n = 0, m = 0;
  in >> n >> m;
  std::string rest;
  std::getline(in, rest);
  std::vector<std::vector<double>> X;
  std::vector<int> lab;
  std::vector<long long> grp;
  X.reserve(n); lab.reserve(n); grp.reserve(n);
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty()) { continue; }
    std::istringstream ss(line);
    long long g; int l;
    ss >> g >> l;
    std::vector<double> x(m);
    for (size_t j = 0; j < m; ++j)
    {
      std::string tok; ss >> tok;
      x[j] = (tok == "nan" || tok == "-nan") ? std::nan("") : std::stod(tok);
    }
    grp.push_back(g); lab.push_back(l); X.push_back(std::move(x));
  }
  std::fprintf(stderr, "loaded %zu rows x %zu cols\n", X.size(), m);

  odia::ScoredGroups s = odia::scoreSemiSupervisedLDA(X, lab, grp);

  std::set<long long> distinct;
  for (double d : s.dscore) { distinct.insert(std::llround(d * 1e6)); }

  std::map<long long, std::pair<double, int>> best;
  for (size_t i = 0; i < X.size(); ++i)
  {
    auto it = best.find(grp[i]);
    if (it == best.end() || s.dscore[i] > it->second.first) { best[grp[i]] = {s.dscore[i], lab[i]}; }
  }
  double tsum = 0, dsum = 0; size_t tn = 0, dn = 0;
  for (const auto& kv : best) { if (kv.second.second == 1) { tsum += kv.second.first; ++tn; } else { dsum += kv.second.first; ++dn; } }
  const double tmean = tn ? tsum / tn : 0, dmean = dn ? dsum / dn : 0;
  std::fprintf(stderr, "distinct d-scores: %zu ; target groups %zu (mean %.3f) decoy %zu (mean %.3f)\n",
               distinct.size(), tn, tmean, dn, dmean);

  int fails = 0;
  if (distinct.size() < 100) { std::fprintf(stderr, "FAIL: d-scores nearly constant (%zu distinct) -> degenerate LDA on real data\n", distinct.size()); ++fails; }
  if (!(tmean > dmean))      { std::fprintf(stderr, "FAIL: targets do not outscore decoys on real data (%.3f vs %.3f)\n", tmean, dmean); ++fails; }
  if (fails) { std::fprintf(stderr, "realdata test FAILED\n"); return 1; }
  std::fprintf(stderr, "realdata test OK\n");
  return 0;
}

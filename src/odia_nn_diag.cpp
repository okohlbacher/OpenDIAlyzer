// Why does the network's seed fit produce no confident positives on real sub-scores?
//
// The ablation says nn = 0 IDs with 0/9 iterations trained, meaning the seed model ranked fewer
// than 26 target groups at q <= 0.15 out of ~100k. The GBT, given the identical seed rows and the
// identical standardised matrix, reaches 449. So this is not the pipeline: it is the network, or
// how it is being asked to learn.
//
// This reproduces JUST the seed step -- standardise, pick each group's best row by the one-feature
// bootstrap, fit, rank -- and reports group-level AUC, which is what the q-values consume. Sweeping
// hyperparameters here costs seconds against a subsampled fixture instead of an hour against the
// full one.
//
// usage: odia-nn-diag <fixture> [lr] [epochs] [batch] [hidden_csv]

#include "odia_gbt.h"
#include "odia_nn.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{

double aucOf(const std::vector<double>& pos, std::vector<double> neg)
{
  if (pos.empty() || neg.empty()) { return 0.5; }
  std::sort(neg.begin(), neg.end());
  double a = 0.0;
  for (const double s : pos)
  {
    const auto lo = std::lower_bound(neg.begin(), neg.end(), s);
    const auto hi = std::upper_bound(neg.begin(), neg.end(), s);
    // Ties count a half, or a constant scorer would read as AUC 0 or 1 instead of 0.5.
    a += (static_cast<double>(lo - neg.begin()) +
          0.5 * static_cast<double>(hi - lo)) / static_cast<double>(neg.size());
  }
  return a / static_cast<double>(pos.size());
}

} // namespace

int main(int argc, char** argv)
{
  if (argc < 2) { std::fprintf(stderr, "usage: odia-nn-diag <fixture> [lr] [epochs] [batch] [h,h,h]\n"); return 2; }

  std::ifstream in(argv[1]);
  if (!in) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
  std::size_t n = 0, m = 0;
  std::string head;
  std::getline(in, head);
  { std::istringstream hs(head); hs >> n >> m; }

  std::vector<std::vector<double>> X;
  std::vector<int> lab;
  std::vector<long long> grp;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty()) { continue; }
    std::istringstream ss(line);
    long long g; int l; std::string tid;
    ss >> g >> l >> tid;
    std::vector<double> x(m);
    for (std::size_t j = 0; j < m; ++j)
    {
      std::string t; ss >> t;
      x[j] = (t == "nan" || t == "-nan") ? std::nan("") : std::strtod(t.c_str(), nullptr);
    }
    grp.push_back(g); lab.push_back(l); X.push_back(std::move(x));
  }
  std::printf("loaded %zu rows x %zu cols\n", X.size(), m);

  // Standardise exactly as odia_lda.h does, including imputing missing to the column mean (0 after
  // centring) -- otherwise this would be diagnosing a different matrix than the one that failed.
  std::vector<double> mean(m, 0.0), sd(m, 0.0);
  for (std::size_t j = 0; j < m; ++j)
  {
    double s = 0.0; std::size_t k = 0;
    for (const auto& r : X) { if (std::isfinite(r[j])) { s += r[j]; ++k; } }
    mean[j] = k ? s / static_cast<double>(k) : 0.0;
    double v = 0.0;
    for (const auto& r : X) { if (std::isfinite(r[j])) { v += (r[j] - mean[j]) * (r[j] - mean[j]); } }
    sd[j] = (k > 1) ? std::sqrt(v / static_cast<double>(k - 1)) : 0.0;
    if (!(sd[j] > 1e-12) || !std::isfinite(sd[j])) { sd[j] = 1.0; }
  }
  std::vector<std::vector<double>> Z(X.size(), std::vector<double>(m, 0.0));
  for (std::size_t i = 0; i < X.size(); ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      Z[i][j] = std::isfinite(X[i][j]) ? (X[i][j] - mean[j]) / sd[j] : 0.0;
    }
  }

  // The one-feature bootstrap: largest |Welch t| between target and decoy rows, then each group's
  // best row by that single feature. This is what selects the seed rows in production.
  std::size_t best_j = 0; double best_t = -1.0, sign = 1.0;
  for (std::size_t j = 0; j < m; ++j)
  {
    double s[2] = {0, 0}, q[2] = {0, 0}; std::size_t c[2] = {0, 0};
    for (std::size_t i = 0; i < Z.size(); ++i)
    { const int L = lab[i]; s[L] += Z[i][j]; q[L] += Z[i][j] * Z[i][j]; ++c[L]; }
    if (c[0] < 2 || c[1] < 2) { continue; }
    const double m0 = s[0] / c[0], m1 = s[1] / c[1];
    const double v0 = q[0] / c[0] - m0 * m0, v1 = q[1] / c[1] - m1 * m1;
    const double se = std::sqrt(std::max(1e-30, v0 / c[0] + v1 / c[1]));
    const double t = (m1 - m0) / se;
    if (std::fabs(t) > best_t) { best_t = std::fabs(t); best_j = j; sign = (t < 0 ? -1.0 : 1.0); }
  }
  std::printf("bootstrap feature: col %zu, |t| = %.1f, sign %+.0f\n", best_j, best_t, sign);

  std::map<long long, std::size_t> best;
  for (std::size_t i = 0; i < Z.size(); ++i)
  {
    auto it = best.find(grp[i]);
    if (it == best.end() || sign * Z[i][best_j] > sign * Z[it->second][best_j]) { best[grp[i]] = i; }
  }
  // HELD-OUT SPLIT BY GROUP. Without it this tool measures training-set fit, which is not a
  // comparison of anything: the first version reported GBT AUC 0.5518 untransformed and 0.9811
  // rank-transformed, and a monotone per-column transform CANNOT change a tree's ranking that much.
  // What it changed was how much the tree could MEMORISE -- 64 equal-width histogram bins put
  // almost every row of a 27-sigma-tailed column into one bin, and ranking spreads them out. Both
  // numbers were training accuracy and neither was evidence.
  std::vector<std::size_t> pos, neg, hpos, hneg;
  for (const auto& kv : best)
  {
    const bool held = (kv.first % 2) != 0;                 // by GROUP id, deterministic
    if (lab[kv.second] == 1) { (held ? hpos : pos).push_back(kv.second); }
    else                     { (held ? hneg : neg).push_back(kv.second); }
  }
  std::printf("seed rows: train %zu T / %zu D, held-out %zu T / %zu D\n", pos.size(), neg.size(),
              hpos.size(), hneg.size());

  auto report = [&](const char* what, auto&& score_fn) {
    std::vector<double> sp, sn;
    for (const std::size_t r : hpos) { sp.push_back(score_fn(Z[r])); }   // HELD OUT
    for (const std::size_t r : hneg) { sn.push_back(score_fn(Z[r])); }
    auto stat = [](std::vector<double> v) {
      std::sort(v.begin(), v.end());
      return std::make_pair(v.front(), v.back());
    };
    const auto ps = stat(sp), ns = stat(sn);
    // Targets ranking above the 99th percentile of decoys: the regime a 1% FDR actually lives in.
    std::vector<double> sorted_n = sn;
    std::sort(sorted_n.begin(), sorted_n.end());
    const double cut = sorted_n[static_cast<std::size_t>(0.99 * sorted_n.size())];
    std::size_t above = 0;
    for (const double s : sp) { if (s > cut) { ++above; } }
    std::printf("  %-22s heldAUC %.4f  above-99th-pct %6zu/%zu (%.4f)  range T[%.3g,%.3g] D[%.3g,%.3g]\n",
                what, aucOf(sp, sn), above, sp.size(),
                static_cast<double>(above) / static_cast<double>(sp.size()),
                ps.first, ps.second, ns.first, ns.second);
  };

  report("bootstrap (1 feature)", [&](const std::vector<double>& x) { return sign * x[best_j]; });

  // INPUT TRANSFORMS. The standardised features reach 27 sigma -- these sub-scores are extremely
  // heavy-tailed. A tree is invariant to any monotone transform of a feature, so the GBT does not
  // care; a tanh network does, because one 27-sigma input saturates the first layer and the
  // gradient through it dies. If that is the difference, a squashing or rank transform closes it
  // and no amount of learning-rate tuning will.
  const char* tf_name = std::getenv("ODIA_NN_TF");
  const std::string tf = tf_name ? tf_name : "none";
  if (tf == "clip3")
  {
    for (auto& r : Z) { for (double& v : r) { v = std::max(-3.0, std::min(3.0, v)); } }
  }
  else if (tf == "tanh")
  {
    for (auto& r : Z) { for (double& v : r) { v = std::tanh(v); } }
  }
  else if (tf == "asinh")
  {
    for (auto& r : Z) { for (double& v : r) { v = std::asinh(v); } }
  }
  else if (tf == "rank")
  {
    // Per-column rank -> centred uniform. Fully monotone, so it discards nothing a tree could use,
    // and bounds every input to [-1, 1] whatever the tail does.
    const std::size_t nrows = Z.size();
    std::vector<std::size_t> ord(nrows);
    for (std::size_t j = 0; j < m; ++j)
    {
      for (std::size_t i = 0; i < nrows; ++i) { ord[i] = i; }
      std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return Z[a][j] < Z[b][j]; });
      // AVERAGE RANKS FOR TIES -- not the sort position. This matters enormously here: 37.5% of
      // the rows share the value 0 in some columns, and the fixture is written CLASS-ORDERED (its
      // first 200k rows are 100% decoy, its last 200k 100% target). Breaking ties by sort position
      // therefore assigns tied rows a rank that encodes WHICH BLOCK OF THE FILE they came from,
      // i.e. the label. That is exactly what the first version of this did, and it reported a
      // held-out GBT AUC of 0.9811 against 0.5250 untransformed -- a number that survived a proper
      // group-wise held-out split and was still pure leakage, because the leak was inside the
      // feature rather than across the split.
      std::vector<double> tmp(nrows);
      std::size_t k = 0;
      while (k < nrows)
      {
        std::size_t e = k;
        while (e + 1 < nrows && Z[ord[e + 1]][j] == Z[ord[k]][j]) { ++e; }
        const double avg = 0.5 * (static_cast<double>(k) + static_cast<double>(e));
        const double val = 2.0 * (avg / static_cast<double>(nrows - 1)) - 1.0;
        for (std::size_t t = k; t <= e; ++t) { tmp[ord[t]] = val; }
        k = e + 1;
      }
      for (std::size_t i = 0; i < nrows; ++i) { Z[i][j] = tmp[i]; }
    }
  }
  std::printf("input transform: %s\n", tf.c_str());

  {
    odia::GBT g;
    odia::GBTParams gp;
    if (g.fit(Z, pos, neg, gp)) { report("GBT", [&](const std::vector<double>& x) { return g.score(x); }); }
    else { std::printf("  GBT                    FIT FAILED\n"); }
  }

  odia::NNParams np;
  if (argc > 2) { np.lr = std::atof(argv[2]); }
  if (argc > 3) { np.epochs = std::atoi(argv[3]); }
  if (argc > 4) { np.batch_size = std::atoi(argv[4]); }
  if (argc > 5)
  {
    np.hidden.clear();
    std::stringstream hs(argv[5]);
    std::string tok;
    while (std::getline(hs, tok, ',')) { np.hidden.push_back(std::atoi(tok.c_str())); }
  }
  char tag[128];
  std::string hid;
  for (std::size_t i = 0; i < np.hidden.size(); ++i)
  { hid += (i ? "-" : "") + std::to_string(np.hidden[i]); }
  std::snprintf(tag, sizeof(tag), "NN lr%.3g e%d b%d h%s", np.lr, np.epochs, np.batch_size, hid.c_str());

  odia::NNEnsemble e;
  if (e.fit(Z, pos, neg, np)) { report(tag, [&](const std::vector<double>& x) { return e.score(x); }); }
  else { std::printf("  %-22s FIT FAILED\n", tag); }
  return 0;
}

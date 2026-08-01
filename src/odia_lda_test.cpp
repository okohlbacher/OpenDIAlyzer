// odia_lda_test.cpp — oracle test for the in-process semi-supervised LDA.
// Builds synthetic target/decoy peak groups where a known subset of target
// precursors are "true hits" with elevated (multivariate) sub-scores; the rest of
// the target rows and all decoy rows are null. A correct semi-supervised LDA must
// (a) rank true-hit rows above null rows, (b) recover a good fraction of true-hit
// precursors at q<0.01, and (c) keep the empirical target-decoy FDR controlled.
//
// Build (standalone, no OpenMS/Eigen):  c++ -std=c++17 -O2 odia_lda_test.cpp -o t && ./t

#include "odia_lda.h"

#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>
#include <cstdio>
#include <random>
#include <vector>

int main()
{
  std::mt19937 rng(123);
  std::normal_distribution<double> noise(0.0, 1.0);

  const int M = 8;                 // sub-scores
  const int n_target_prec = 2000;  // target precursors
  const int n_decoy_prec = 2000;   // decoy precursors
  const int groups_per_prec = 4;   // candidate peak groups per precursor
  const double true_frac = 0.6;    // fraction of target precursors that are real hits
  const double mu = 2.5;           // signal on the informative dims (strong: honest q<0.01
                                   // recall must clear 0.35 under the target-decoy FDR rule)
  const int informative[] = {0, 2, 5};  // only some dims carry signal

  std::vector<std::vector<double>> feats;
  std::vector<int> labels;
  std::vector<long long> group;
  std::vector<char> prec_is_true;   // per precursor id: is it a real hit
  long long gid = 0;

  auto add_prec = [&](bool is_target, bool is_true) {
    const long long id = gid++;
    prec_is_true.push_back(is_true ? 1 : 0);
    for (int g = 0; g < groups_per_prec; ++g)
    {
      std::vector<double> x(M);
      for (int j = 0; j < M; ++j) { x[j] = noise(rng); }
      // exactly one peak group of a true-hit target carries the signal
      if (is_true && g == 0) { for (int k : informative) { x[k] += mu; } }
      feats.push_back(std::move(x));
      labels.push_back(is_target ? 1 : 0);
      group.push_back(id);
    }
  };

  std::uniform_real_distribution<double> u(0.0, 1.0);
  int n_true = 0;
  for (int i = 0; i < n_target_prec; ++i) { bool t = u(rng) < true_frac; n_true += t; add_prec(true, t); }
  for (int i = 0; i < n_decoy_prec; ++i)  { add_prec(false, false); }

  odia::LDAParams p;  // defaults
  odia::ScoredGroups s = odia::scoreSemiSupervisedLDA(feats, labels, group, p);

  // ---- row order must not change the result ----
  // Fold assignment used a group's first-occurrence POSITION, and row order is whatever the
  // parallel extraction produced -- so identical input gave 6487 / 6565 / 6433 IDs across three
  // real runs. Permute the rows, score again, and require the per-row d-scores to match exactly.
  {
    std::vector<std::size_t> perm(feats.size());
    for (std::size_t i = 0; i < perm.size(); ++i) { perm[i] = i; }
    std::mt19937 prng(999);
    std::shuffle(perm.begin(), perm.end(), prng);

    std::vector<std::vector<double>> f2(feats.size());
    std::vector<int> l2(feats.size());
    std::vector<long long> g2_ids(feats.size());
    for (std::size_t i = 0; i < perm.size(); ++i)
    {
      f2[i] = feats[perm[i]]; l2[i] = labels[perm[i]]; g2_ids[i] = group[perm[i]];
    }
    const odia::ScoredGroups s2 = odia::scoreSemiSupervisedLDA(f2, l2, g2_ids, p);
    double worst = 0.0;
    for (std::size_t i = 0; i < perm.size(); ++i)
    {
      worst = std::max(worst, std::fabs(s2.dscore[i] - s.dscore[perm[i]]));
    }
    std::fprintf(stderr, "permuted-row d-score max |delta| = %.3e\n", worst);
    if (worst > 1e-9)
    {
      std::fprintf(stderr, "FAIL: row order changed the scores (max delta %.3e)\n", worst);
      return 1;
    }

    // GBT is the classifier the benchmark actually runs (-classifier gbt), and it is NOT
    // order-invariant on raw input: its histogram reduction partitions rows into chunks by INDEX
    // (lo = n_rows*c/nchunk), so permuting rows regroups the partial sums and the floating-point
    // difference flips split points. Measured max |delta| 1.965 -- and only with OpenMP compiled
    // in, identical at 1/8/64 threads, so it reads as deterministic under a serial test build.
    //
    // Chunking by index is deliberate: that loop is ~93% of a fit and contiguous chunks are what
    // make it cache-friendly. So determinism is supplied UPSTREAM -- OswRows::canonicalize() sorts
    // by (precursor, apex RT, feature id) before scoring. This asserts that contract: canonically
    // ordered input gives identical scores no matter what order the rows arrived in.
    odia::LDAParams pg = p;
    pg.classifier = odia::LDAParams::Classifier::GBT;
    const auto canon = [&](std::vector<std::vector<double>> F, std::vector<int> L,
                           std::vector<long long> G) {
      std::vector<std::size_t> ix(F.size());
      for (std::size_t i = 0; i < ix.size(); ++i) { ix[i] = i; }
      // same key shape as OswRows::canonicalize(): group first, then a per-row physical value
      std::stable_sort(ix.begin(), ix.end(), [&](std::size_t a, std::size_t b) {
        if (G[a] != G[b]) { return G[a] < G[b]; }
        return F[a][0] < F[b][0];
      });
      std::vector<std::vector<double>> F2; std::vector<int> L2; std::vector<long long> G2;
      for (std::size_t i : ix) { F2.push_back(F[i]); L2.push_back(L[i]); G2.push_back(G[i]); }
      return std::make_tuple(F2, L2, G2);
    };
    const auto [fa, la, ga] = canon(feats, labels, group);
    const auto [fb, lb, gb] = canon(f2, l2, g2_ids);
    const odia::ScoredGroups g1 = odia::scoreSemiSupervisedLDA(fa, la, ga, pg);
    const odia::ScoredGroups gg2 = odia::scoreSemiSupervisedLDA(fb, lb, gb, pg);
    double gworst = 0.0;
    for (std::size_t i = 0; i < g1.dscore.size(); ++i)
    {
      gworst = std::max(gworst, std::fabs(gg2.dscore[i] - g1.dscore[i]));
    }
    std::fprintf(stderr, "canonicalised-then-permuted GBT d-score max |delta| = %.3e\n", gworst);
    if (gworst > 1e-9)
    {
      std::fprintf(stderr, "FAIL: canonical ordering did not make GBT reproducible (%.3e)\n", gworst);
      return 1;
    }
  }

  // ---- basic contract ----
  if (s.dscore.size() != feats.size() || s.qvalue.size() != feats.size())
  {
    std::fprintf(stderr, "FAIL: output size mismatch\n");
    return 1;
  }

  // ---- separation: true-hit rows should outscore null rows on average ----
  double sum_true = 0, sum_null = 0; int n_t = 0, n_n = 0;
  for (size_t i = 0; i < feats.size(); ++i)
  {
    const bool true_hit = (labels[i] == 1 && prec_is_true[group[i]] && (i % groups_per_prec) == 0);
    if (true_hit) { sum_true += s.dscore[i]; ++n_t; }
    else          { sum_null += s.dscore[i]; ++n_n; }
  }
  const double mean_true = sum_true / n_t, mean_null = sum_null / n_n;
  if (!(mean_true > mean_null))
  {
    std::fprintf(stderr, "FAIL: true-hit d-score %.3f not above null %.3f\n", mean_true, mean_null);
    return 1;
  }

  // ---- q<0.01: recover a good fraction of true targets, controlled FDR ----
  // Count distinct precursors whose best row is at q<0.01, split target/decoy.
  std::vector<double> best_q(gid, 2.0);
  std::vector<int> best_lab(gid, -1);
  std::vector<double> best_d(gid, -1e300);
  for (size_t i = 0; i < feats.size(); ++i)
  {
    if (s.dscore[i] > best_d[group[i]]) { best_d[group[i]] = s.dscore[i]; best_q[group[i]] = s.qvalue[i]; best_lab[group[i]] = labels[i]; }
  }
  int rec_target = 0, rec_decoy = 0, rec_true = 0;
  for (long long id = 0; id < gid; ++id)
  {
    if (best_q[id] < 0.01)
    {
      if (best_lab[id] == 1) { ++rec_target; if (prec_is_true[id]) ++rec_true; }
      else                   { ++rec_decoy; }
    }
  }
  const double emp_fdr = rec_target > 0 ? (double)rec_decoy / (double)rec_target : 1.0;
  const double recall = (double)rec_true / (double)n_true;
  std::fprintf(stderr, "recovered target-precursors@q<0.01=%d (true=%d/%d, recall=%.2f), decoys=%d, emp_FDR=%.3f\n",
               rec_target, rec_true, n_true, recall, rec_decoy, emp_fdr);

  if (recall < 0.35)   { std::fprintf(stderr, "FAIL: recall %.2f < 0.35\n", recall); return 1; }
  if (emp_fdr > 0.05)  { std::fprintf(stderr, "FAIL: empirical FDR %.3f > 0.05\n", emp_fdr); return 1; }

  // p-value / q-value / PEP must be three DISTINCT statistics. Every .osw written before this had
  // all three bound to the q-value, which is invisible unless something checks. Guard the two
  // properties that make them distinct and the one that makes PEP usable:
  //   (a) they are not all identical;
  //   (b) PEP >= q at the same peak group -- LOCAL FDR at a threshold cannot be below the AVERAGE
  //       FDR of everything above it (that is what the two mean);
  //   (c) PEP is monotone non-increasing in d-score.
  {
    std::size_t n_all_equal = 0, n_checked = 0, n_pep_below_q = 0;
    for (std::size_t i = 0; i < s.qvalue.size(); ++i)
    {
      if (s.pvalue.empty() || s.pep.empty()) { break; }
      ++n_checked;
      if (s.pvalue[i] == s.qvalue[i] && s.qvalue[i] == s.pep[i]) { ++n_all_equal; }
      if (s.pep[i] < s.qvalue[i] - 1e-9) { ++n_pep_below_q; }
    }
    if (n_checked == 0)
    {
      std::fprintf(stderr, "FAIL: no p-value/PEP produced\n"); return 1;
    }
    if (n_all_equal == n_checked)
    {
      std::fprintf(stderr, "FAIL: p-value, q-value and PEP are identical on all %zu rows\n", n_checked);
      return 1;
    }
    // A handful of ties at the extremes is fine; a systematic violation is not.
    if (n_pep_below_q > n_checked / 20)
    {
      std::fprintf(stderr, "FAIL: PEP < q-value on %zu/%zu rows (local FDR below average FDR)\n",
                   n_pep_below_q, n_checked);
      return 1;
    }
    // Monotonicity of PEP in d-score.
    std::vector<std::pair<double, double>> ds;
    ds.reserve(n_checked);
    for (std::size_t i = 0; i < n_checked; ++i) { ds.emplace_back(s.dscore[i], s.pep[i]); }
    std::sort(ds.begin(), ds.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    for (std::size_t i = 1; i < ds.size(); ++i)
    {
      if (ds[i].second < ds[i - 1].second - 1e-9)
      {
        std::fprintf(stderr, "FAIL: PEP increases with d-score at rank %zu (%.4f -> %.4f)\n",
                     i, ds[i - 1].second, ds[i].second);
        return 1;
      }
    }
    std::fprintf(stderr, "p/q/PEP distinct OK (all-equal rows: %zu/%zu)\n", n_all_equal, n_checked);
  }

  std::fprintf(stderr, "odia_lda_test OK\n");
  return 0;
}

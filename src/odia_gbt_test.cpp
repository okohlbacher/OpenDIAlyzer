// odia_gbt_test.cpp — self-check for the histogram GBT.
//
// The claim being tested is NOT "boosting is better" in the abstract. It is the specific reason to
// prefer it here: OpenSWATH sub-scores interact non-linearly, and a single hyperplane cannot express
// that. So the decisive test is a dataset where the classes are NOT linearly separable but ARE
// separable by an interaction — the GBT must win there, and must NOT lose on the linear case.
//
//   T1 determinism: identical output across repeated fits, and ties broken reproducibly.
//   T2 linear data: GBT is competitive with LDA (it must not PAY for the extra capacity).
//   T3 XOR-like interaction: GBT wins decisively where LDA is at chance. This is the whole point.
//   T4 the model actually learns: training AUC improves over the boosting rounds.
//   T5 missing values are a splittable category, not silently imputed.
//   T6 degenerate inputs are refused rather than producing a garbage model.
//
// Build: c++ -std=c++17 -O2 src/odia_gbt_test.cpp -o t && ./t

#include "odia_gbt.h"
#include "odia_lda.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

static int failures = 0;
static void check(bool ok, const char* what)
{
  if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++failures; }
}

/// Area under the ROC curve from scores and labels (rank-based, ties averaged).
static double auc(const std::vector<double>& s, const std::vector<int>& y)
{
  std::vector<std::size_t> idx(s.size());
  for (std::size_t i = 0; i < idx.size(); ++i) { idx[i] = i; }
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return s[a] < s[b]; });
  double rank_sum = 0.0;
  std::size_t n_pos = 0, n_neg = 0;
  for (std::size_t k = 0; k < idx.size();)
  {
    std::size_t j = k;
    while (j < idx.size() && s[idx[j]] == s[idx[k]]) { ++j; }
    const double avg_rank = 0.5 * (double(k + 1) + double(j)); // 1-based, ties averaged
    for (std::size_t t = k; t < j; ++t)
    {
      if (y[idx[t]] == 1) { rank_sum += avg_rank; ++n_pos; } else { ++n_neg; }
    }
    k = j;
  }
  if (n_pos == 0 || n_neg == 0) { return 0.5; }
  return (rank_sum - double(n_pos) * (n_pos + 1) / 2.0) / (double(n_pos) * double(n_neg));
}

/// Fit an LDA on the same two row sets, for a like-for-like comparison. Fisher's discriminant with
/// a small ridge -- the same estimator odia_lda.h uses inside its loop.
static std::vector<double> ldaScores(const std::vector<std::vector<double>>& X,
                                     const std::vector<std::size_t>& pos,
                                     const std::vector<std::size_t>& neg)
{
  const std::size_t m = X[0].size();
  std::vector<double> mp(m, 0.0), mn(m, 0.0);
  for (auto r : pos) { for (std::size_t j = 0; j < m; ++j) { mp[j] += X[r][j]; } }
  for (auto r : neg) { for (std::size_t j = 0; j < m; ++j) { mn[j] += X[r][j]; } }
  for (std::size_t j = 0; j < m; ++j) { mp[j] /= pos.size(); mn[j] /= neg.size(); }
  std::vector<double> S(m * m, 0.0);
  auto acc = [&](const std::vector<std::size_t>& rows, const std::vector<double>& mu) {
    for (auto r : rows)
      for (std::size_t j = 0; j < m; ++j)
        for (std::size_t k = 0; k <= j; ++k)
          S[j * m + k] += (X[r][j] - mu[j]) * (X[r][k] - mu[k]);
  };
  acc(pos, mp); acc(neg, mn);
  const double dof = double(pos.size() + neg.size() - 2);
  for (std::size_t j = 0; j < m; ++j)
    for (std::size_t k = 0; k <= j; ++k) { S[j*m+k] /= dof; S[k*m+j] = S[j*m+k]; S[j*m+j] += 0.0; }
  for (std::size_t j = 0; j < m; ++j) { S[j * m + j] += 1e-6; }
  std::vector<double> d(m), w;
  for (std::size_t j = 0; j < m; ++j) { d[j] = mp[j] - mn[j]; }
  if (!odia::lda_detail::choleskySolve(S, d, w)) { w.assign(m, 0.0); }
  std::vector<double> out(X.size(), 0.0);
  for (std::size_t i = 0; i < X.size(); ++i) { out[i] = odia::lda_detail::dot(w, X[i]); }
  return out;
}

int main()
{
  // ---- T2 / T3: build both a LINEAR and an INTERACTION dataset ------------------------
  // Same size, same noise, same class balance; the only difference is how the label is generated.
  std::mt19937 rng(20260731);
  std::normal_distribution<double> nd(0.0, 1.0);
  const std::size_t N = 6000, M = 8;

  auto build = [&](bool interaction, std::vector<std::vector<double>>& X, std::vector<int>& y,
                   std::vector<std::size_t>& pos, std::vector<std::size_t>& neg) {
    X.assign(N, std::vector<double>(M, 0.0));
    y.assign(N, 0);
    pos.clear(); neg.clear();
    for (std::size_t i = 0; i < N; ++i)
    {
      for (std::size_t j = 0; j < M; ++j) { X[i][j] = nd(rng); }
      // Linear: a plane through features 0 and 1.
      // Interaction: the SIGN of the product -- an XOR, which no hyperplane separates. The
      // remaining 6 features are pure noise in both cases, so the classifier must also ignore them.
      const double signal = interaction ? (X[i][0] * X[i][1]) : (X[i][0] + X[i][1]);
      y[i] = (signal + 0.35 * nd(rng) > 0.0) ? 1 : 0;
      (y[i] == 1 ? pos : neg).push_back(i);
    }
  };

  odia::GBTParams gp;   // defaults

  for (int interaction = 0; interaction <= 1; ++interaction)
  {
    std::vector<std::vector<double>> X; std::vector<int> y;
    std::vector<std::size_t> pos, neg;
    build(interaction != 0, X, y, pos, neg);

    odia::GBT gbt;
    const bool ok = gbt.fit(X, pos, neg, gp);
    check(ok, interaction ? "T3 GBT fits the interaction data" : "T2 GBT fits the linear data");
    std::vector<double> gs(N);
    for (std::size_t i = 0; i < N; ++i) { gs[i] = gbt.score(X[i]); }
    const double a_gbt = auc(gs, y);
    const double a_lda = auc(ldaScores(X, pos, neg), y);

    std::fprintf(stderr, "  %s data: LDA AUC %.4f   GBT AUC %.4f   (%zu trees)\n",
                 interaction ? "INTERACTION" : "linear     ", a_lda, a_gbt, gbt.nTrees());

    if (!interaction)
    {
      // The linear case is the LDA's home ground. The GBT must not be badly worse -- if it is, the
      // extra capacity is costing more than it earns and the default params are wrong.
      check(a_gbt > a_lda - 0.03, "T2 GBT stays within 0.03 AUC of LDA on linear data");
      check(a_lda > 0.85, "T2 sanity: LDA does well on linear data");
    }
    else
    {
      // The whole reason to add boosting. LDA should be near chance; GBT should be decisively better.
      check(a_lda < 0.60, "T3 sanity: LDA is near chance on an interaction (no hyperplane exists)");
      check(a_gbt > 0.80, "T3 GBT learns the interaction");
      check(a_gbt > a_lda + 0.20, "T3 GBT beats LDA by a wide margin on interaction data");
    }
  }

  // ---- T1 determinism -----------------------------------------------------------------
  {
    std::vector<std::vector<double>> X; std::vector<int> y;
    std::vector<std::size_t> pos, neg;
    build(true, X, y, pos, neg);
    odia::GBT a, b;
    a.fit(X, pos, neg, gp);
    b.fit(X, pos, neg, gp);
    double worst = 0.0;
    for (std::size_t i = 0; i < N; ++i) { worst = std::max(worst, std::abs(a.score(X[i]) - b.score(X[i]))); }
    check(worst == 0.0, "T1 two fits on identical data are BIT-identical (no RNG anywhere)");
    // And row order must not matter beyond the deliberate tie-break rule: reversing the negative
    // set changes accumulation order but not the histogram sums, so the model must be unchanged.
    std::vector<std::size_t> neg_rev(neg.rbegin(), neg.rend());
    odia::GBT c;
    c.fit(X, pos, neg_rev, gp);
    double worst2 = 0.0;
    for (std::size_t i = 0; i < N; ++i) { worst2 = std::max(worst2, std::abs(a.score(X[i]) - c.score(X[i]))); }
    std::fprintf(stderr, "  T1 max |delta| across fits: %.3e (identical), %.3e (row order reversed)\n",
                 worst, worst2);
    check(worst2 < 1e-9, "T1 reversing row order does not change the model");
  }

  // ---- T4 the model actually learns as rounds accumulate --------------------------------
  {
    std::vector<std::vector<double>> X; std::vector<int> y;
    std::vector<std::size_t> pos, neg;
    build(true, X, y, pos, neg);
    double prev = 0.5;
    bool monotone_ish = true;
    for (int n : {1, 5, 20, 120})
    {
      odia::GBTParams p = gp; p.n_trees = n;
      odia::GBT g; g.fit(X, pos, neg, p);
      std::vector<double> s(N);
      for (std::size_t i = 0; i < N; ++i) { s[i] = g.score(X[i]); }
      const double a = auc(s, y);
      std::fprintf(stderr, "  T4 %3d trees -> AUC %.4f\n", n, a);
      if (a < prev - 0.02) { monotone_ish = false; }
      prev = a;
    }
    check(monotone_ish, "T4 AUC does not degrade as trees are added");
    check(prev > 0.80, "T4 the full model is well above chance");
  }

  // ---- T5 missing values are a category, not an imputation -------------------------------
  {
    // Feature 0 is missing for exactly the positives. A model that imputes to a mean learns nothing;
    // one that treats missing as its own bin separates perfectly.
    const std::size_t n = 2000;
    std::vector<std::vector<double>> X(n, std::vector<double>(2, 0.0));
    std::vector<std::size_t> pos, neg;
    for (std::size_t i = 0; i < n; ++i)
    {
      const bool is_pos = (i % 2 == 0);
      X[i][0] = is_pos ? std::nan("") : nd(rng);
      X[i][1] = nd(rng);
      (is_pos ? pos : neg).push_back(i);
    }
    odia::GBT g;
    check(g.fit(X, pos, neg, gp), "T5 fits with missing values present");
    std::vector<double> s(n); std::vector<int> yy(n);
    for (std::size_t i = 0; i < n; ++i) { s[i] = g.score(X[i]); yy[i] = (i % 2 == 0) ? 1 : 0; }
    const double a = auc(s, yy);
    std::fprintf(stderr, "  T5 missing-as-category AUC: %.4f\n", a);
    check(a > 0.95, "T5 missingness is learnable (treated as its own bin, not imputed)");
  }

  // ---- T6 degenerate inputs --------------------------------------------------------------
  {
    std::vector<std::vector<double>> X(10, std::vector<double>(3, 1.0));
    std::vector<std::size_t> one{0}, many{1, 2, 3};
    odia::GBT g;
    check(!g.fit(X, one, many, gp), "T6 refuses a single-row positive class");
    check(!g.fit(X, many, one, gp), "T6 refuses a single-row negative class");
    std::vector<std::vector<double>> empty;
    check(!g.fit(empty, many, many, gp), "T6 refuses empty data");
    check(g.score(std::vector<double>(3, 0.0)) == 0.0, "T6 an untrained model scores 0, not garbage");
  }

  // ---- T7 end-to-end through the REAL semi-supervised loop -------------------------------
  // Everything above tests the learner in isolation. This drives scoreSemiSupervisedLDA() itself
  // with both classifiers on identical input, which is what actually ships: same folds, same
  // training-set selection, same fold normalisation, same q-values -- only the learner differs.
  {
    std::mt19937 r2(99);
    std::normal_distribution<double> g(0.0, 1.0);
    const std::size_t n_groups = 3000, per_group = 4, M2 = 6;
    std::vector<std::vector<double>> X;
    std::vector<int> lab;
    std::vector<long long> grp;
    std::size_t n_true = 0;
    for (std::size_t gi = 0; gi < n_groups; ++gi)
    {
      const bool is_target = (gi % 2 == 0);
      // A quarter of target precursors are real hits, and their signal is an INTERACTION between
      // two sub-scores -- the case a hyperplane cannot see.
      const bool is_hit = is_target && (gi % 4 == 0);
      if (is_hit) { ++n_true; }
      for (std::size_t k = 0; k < per_group; ++k)
      {
        std::vector<double> row(M2);
        for (std::size_t j = 0; j < M2; ++j) { row[j] = g(r2); }
        if (is_hit && k == 0)
        {
          // push both features to the same side so their PRODUCT is large and positive
          const double sgn = (gi % 16 == 0) ? 1.0 : -1.0;
          row[0] = sgn * (2.2 + 0.3 * std::abs(g(r2)));
          row[1] = sgn * (2.2 + 0.3 * std::abs(g(r2)));
        }
        X.push_back(std::move(row));
        lab.push_back(is_target ? 1 : 0);
        grp.push_back(static_cast<long long>(gi));
      }
    }
    auto run = [&](odia::LDAParams::Classifier c) {
      odia::LDAParams p;
      p.classifier = c;
      auto s = odia::scoreSemiSupervisedLDA(X, lab, grp, p);
      std::size_t ids = 0, dec = 0;
      std::vector<char> seen(n_groups, 0);
      for (std::size_t i = 0; i < X.size(); ++i)
      {
        const std::size_t gi = static_cast<std::size_t>(grp[i]);
        if (seen[gi] || s.qvalue[i] > 0.01) { continue; }
        seen[gi] = 1;
        (lab[i] == 1 ? ids : dec)++;
      }
      // AUC of the cross-validated d-score over top rows: says whether the loop learned anything
      // at all, independently of whether the q-value cutoff was reached.
      std::vector<double> best(n_groups, -1e300);
      std::vector<int> blab(n_groups, 0);
      for (std::size_t i = 0; i < X.size(); ++i)
      {
        const std::size_t gi = static_cast<std::size_t>(grp[i]);
        if (s.dscore[i] > best[gi]) { best[gi] = s.dscore[i]; blab[gi] = lab[i]; }
      }
      std::fprintf(stderr, "         (d-score AUC target-vs-decoy = %.4f)\n", auc(best, blab));
      return std::make_pair(ids, dec);
    };
    const auto a = run(odia::LDAParams::Classifier::LDA);
    const auto b = run(odia::LDAParams::Classifier::GBT);
    std::fprintf(stderr, "  T7 full loop, interaction signal (%zu true hits):\n", n_true);
    std::fprintf(stderr, "       LDA  %zu IDs @q<0.01 (%zu decoys)\n", a.first, a.second);
    std::fprintf(stderr, "       GBT  %zu IDs @q<0.01 (%zu decoys)\n", b.first, b.second);
    check(b.first > a.first, "T7 GBT recovers more true precursors than LDA through the full loop");
    check(b.second <= b.first / 20 + 5, "T7 GBT keeps the empirical FDR controlled");
  }

  // ---- T8 thread-count independence -------------------------------------------------------
  // The histogram is a floating-point reduction, so the number of partial sums decides the
  // rounding. Partitioning by THREAD would make the model depend on OMP_NUM_THREADS; the chunk
  // count depends only on the row count, so it must not. This is the claim that makes the
  // parallelisation safe for an FDR-producing score, so it is tested rather than asserted.
  {
    std::vector<std::vector<double>> X; std::vector<int> y;
    std::vector<std::size_t> pos, neg;
    build(true, X, y, pos, neg);
    odia::GBTParams p1 = gp; p1.n_threads = 1;
    odia::GBTParams p8 = gp; p8.n_threads = 8;
    odia::GBTParams p64 = gp; p64.n_threads = 64;
    odia::GBT g1, g8, g64;
    g1.fit(X, pos, neg, p1);
    g8.fit(X, pos, neg, p8);
    g64.fit(X, pos, neg, p64);
    double d8 = 0.0, d64 = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
      d8  = std::max(d8,  std::abs(g1.score(X[i]) - g8.score(X[i])));
      d64 = std::max(d64, std::abs(g1.score(X[i]) - g64.score(X[i])));
    }
    std::fprintf(stderr, "  T8 max |delta| vs 1 thread:  8 threads %.3e   64 threads %.3e\n", d8, d64);
    check(d8 == 0.0 && d64 == 0.0,
          "T8 the model is BIT-identical at 1, 8 and 64 threads (chunked, not per-thread, reduction)");
  }

  if (failures) { std::fprintf(stderr, "odia_gbt_test FAILED (%d)\n", failures); return 1; }
  std::fprintf(stderr, "odia_gbt_test OK\n");
  return 0;
}

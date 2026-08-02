// odia_lda.h — in-process semi-supervised LDA scorer for OpenDIAlyzer.
//
// A pyprophet/mProphet-style confidence model, reimplemented in-process so
// OpenDIAlyzer needs no external pyprophet. Pure C++ (no OpenMS, no Eigen), so it
// compiles and self-tests standalone. Used in two places:
//   (1) recalibration: pick RT anchors as targets at q < threshold (a real
//       semi-supervised d-score, not raw VAR_XCORR_SHAPE/VAR_LIBRARY_CORR cutoffs);
//   (2) final FDR: q-values on the last pass.
//
// Start with LDA (deterministic, dependency-light). Gradient-boosted trees are the
// backlog upgrade (see docs/OpenDIAlyzer-scoring-fdr-backlog.md).
//
// ALGORITHM (target for implementation):
//   Inputs: features (N rows x M sub-scores), labels (1=target, 0=decoy), group
//   (precursor id; a precursor has several candidate peak groups = rows). FDR is a
//   per-precursor quantity, so it is computed on the best-scoring row per group.
//   1. z-standardize each feature column (mean/sd over all rows).
//   2. k-fold cross-validation BY GROUP (all rows of a precursor fall in one fold),
//      so a row is never scored by a model trained on its own precursor.
//   3. For each fold, train on the other folds by semi-supervised iteration:
//        - initialise the discriminant with the single most target/decoy-separating
//          feature (max |two-sample t|);
//        - pick a confident-target training set: best row per target group whose
//          current score beats the decoy-derived cutoff at a lenient train FDR
//          (e.g. 0.05); all decoy rows are negatives;
//        - Fisher LDA: w = Sw^{-1} (mu_pos - mu_neg), Sw = pooled within-class
//          covariance (solve SPD system by Cholesky; ridge-regularise the diagonal
//          if needed); score = w . x; repeat n_iter times.
//      Apply the fold's final w to its held-out rows -> cross-validated d-score.
//   4. q-values: take the best d-score per precursor (targets and decoys), sort
//      descending; at each cut FDR = #decoys_above / #targets_above (target-decoy),
//      monotonise to q-values; broadcast each precursor's q to its rows.
//   Deterministic given `seed` (only the fold assignment is randomised).
//
#ifndef ODIA_LDA_H
#define ODIA_LDA_H

#include "odia_gbt.h"
#include "odia_anchor_training.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

namespace odia
{

struct ScoredGroups
{
  std::vector<double> dscore;   ///< per input row (peak group)
  std::vector<double> qvalue;   ///< per input row (its precursor's q-value)
  std::vector<double> pvalue;   ///< per input row: tail probability under the decoy null
  std::vector<double> pep;      ///< per input row: local FDR (posterior error probability)
  /// Diagnostics: how many semi-supervised iterations actually FITTED a discriminant, and how many
  /// bailed out for want of confident positives. If trained==0 the returned scores come from the
  /// single-feature initialisation, not from a learned model -- which looks like a working LDA and
  /// is not one.
  int n_iterations_trained = 0;
  int n_iterations_skipped = 0;
};

struct LDAParams
{
  int n_folds = 3;          ///< cross-validation folds (by group)
  int n_iter = 3;           ///< semi-supervised iterations
  double train_fdr_initial = 0.15;  ///< FDR for the FIRST training-set selection. Deliberately more
                            ///< lenient than the rest: the first pass is seeded by a single feature,
                            ///< so a strict cut can select too few positives to fit anything and the
                            ///< model never bootstraps. pyprophet uses 0.15 here for the same reason.
  double train_fdr = 0.05;  ///< FDR for subsequent iterations, once a real discriminant exists.
  double ridge = 1e-6;      ///< diagonal regularisation for the within-class covariance solve
  unsigned seed = 42;       ///< RNG seed (fold assignment only) — determinism
  bool use_pi0 = false;     ///< Storey pi0 correction. false = HONEST/conservative (true 1% FDR,
                            ///< fewer IDs); true = pyprophet/DIA-NN parity (more IDs, but a nominal
                            ///< 1% is ~2% actual — matches their calibration, incl. its optimism).
  bool top_decoys_only = true;  ///< Train the negative class on each decoy precursor's BEST peak
                            ///< group only, not on all of its candidate peak groups. See the note at
                            ///< the negative-class construction: the positive class is top-peaks-only
                            ///< by definition, so using all decoy rows makes the two classes differ
                            ///< in peak RANK as well as in label, and the discriminant partly learns
                            ///< "is this the best peak of its group" instead of "is this a target".
                            ///< pyprophet trains on get_top_decoy_peaks() for exactly this reason.
  /// Which learner fills the semi-supervised loop. LDA fits one hyperplane; GBT fits an additive
  /// ensemble of depth-limited trees and can express interactions between sub-scores that no
  /// hyperplane can. Measured on synthetic data with the SAME loop around it: on linearly separable
  /// data LDA 0.981 vs GBT 0.986 AUC (no penalty for the extra capacity), on an interaction LDA
  /// 0.514 -- chance -- vs GBT 0.926. Everything else (folds, training-set selection, fold
  /// normalisation, q-values) is identical, so the two are directly comparable.
  ///
  /// NN adds a third: an ensemble of small tanh MLPs, the shape DIA-NN uses. It can express the
  /// same interactions the GBT can, plus smooth ones a tree approximates in steps.
  enum class Classifier { LDA, GBT, NN };
  Classifier classifier = Classifier::LDA;
  GBTParams gbt;            ///< only consulted when classifier == GBT
  NNParams nn;              ///< only consulted when classifier == NN

  // ---- the anti-circularity mechanisms of odia_anchor_training.h, each toggleable ALONE ---------
  // They are separable on purpose. Bundled, an improvement and a regression cancel and the net
  // result is "no effect"; separately, each one's contribution is attributable.

  /// Mechanism 1. Feature mask for the SEED fit only (empty = no exclusion). The seed selects the
  /// anchors; if it selects them using the very features they will be used to calibrate, the
  /// correction is biased toward the uncorrected state.
  std::vector<char> seed_mask;
  /// Mechanism 3. Share of precursor GROUPS each ensemble member trains on (1.0 = off, which is
  /// DIA-NN's arrangement: members differ only in weight init, so they share one seed bias).
  double bag_fraction = 1.0;
  /// Mechanism 5. Stop when the positive SET stops changing, or when it SHRINKS. Counting
  /// identifications cannot see a collapse -- they rise throughout one.
  bool stop_on_composition = false;
  double stop_jaccard = 0.98;
  /// DIAGNOSTIC ONLY, and FDR-INVALID when true: every group trains the model that scores it.
  ///
  /// It exists because `n_folds = 1` cannot express this -- the fold count is clamped to >= 2 a few
  /// lines into the routine, so an ablation arm that set n_folds=1 silently ran 2-fold CV and
  /// "demonstrated" nothing while claiming to demonstrate that cross-validation is load-bearing.
  /// A flag that says what it does cannot be defeated by a clamp.
  bool disable_cv = false;
  bool normalize_folds = true;  ///< Rescale each fold's held-out scores to its own decoy null
                            ///< (mean 0, sd 1) before pooling. Each fold has its OWN weight vector,
                            ///< with its own arbitrary scale and offset, so the raw scores are not
                            ///< comparable across folds; pooling them into one ranking without this
                            ///< mixes incommensurable scales.
};

/// Semi-supervised LDA scoring with cross-validation and target-decoy q-values.
/// features[i] has the same length for all i (M sub-scores). labels[i] in {0,1}.
/// group[i] is the precursor id shared by that precursor's candidate peak groups.
/// Returns a d-score and q-value per input row. Deterministic given params.seed.
ScoredGroups scoreSemiSupervisedLDA(const std::vector<std::vector<double>>& features,
                                    const std::vector<int>& labels,
                                    const std::vector<long long>& group,
                                    const LDAParams& params = LDAParams());

namespace lda_detail
{

struct RankedGroup
{
  std::size_t group_index;
  std::size_t best_row;
  int label;
  double score;
  double qvalue;
  double pvalue = 1.0;   ///< empirical p from the decoy null (TAIL probability at this score)
  double pep = 1.0;      ///< posterior error probability = LOCAL false-discovery rate at this score
};

inline double dot(const std::vector<double>& a, const std::vector<double>& b)
{
  double result = 0.0;
  for (std::size_t j = 0; j < a.size(); ++j) { result += a[j] * b[j]; }
  return result;
}

// Assign target-decoy q-values to a list containing one best score per group.
// Equal scores are treated as one threshold, avoiding order-dependent q-values.
inline void assignQValues(std::vector<RankedGroup>& ranked, bool use_pi0)
{
  std::sort(ranked.begin(), ranked.end(), [](const RankedGroup& a, const RankedGroup& b) {
    if (a.score != b.score) { return a.score > b.score; }
    return a.group_index < b.group_index;
  });

  std::size_t Ntar = 0, Ndec = 0;
  for (const auto& r : ranked) { if (r.label == 1) { ++Ntar; } else { ++Ndec; } }

  // Storey pi0 = estimated fraction of TRUE-NULL targets. Walking high->low score, a
  // target's empirical p-value from the decoy null is (decoys_seen_so_far / Ndec); null
  // targets have ~uniform p, so the mass with p>lambda estimates pi0. Without pi0 the
  // target-decoy estimator assumes pi0=1 (every target could be false) and is
  // over-conservative vs pyprophet/DIA-NN (measured: 22,959 vs 37,539 IDs on the same osw,
  // despite our discriminant separating MORE clean targets). pi0 recovers the honest count.
  double pi0 = 1.0;
  if (use_pi0 && Ntar > 0 && Ndec > 0)
  {
    const double lambda = 0.5;
    const std::size_t dthr = static_cast<std::size_t>(lambda * static_cast<double>(Ndec));
    std::size_t dec_seen = 0, tar_hi = 0;
    for (const auto& r : ranked)
    {
      if (r.label == 0) { ++dec_seen; }
      else if (dec_seen > dthr) { ++tar_hi; }
    }
    pi0 = static_cast<double>(tar_hi) / ((1.0 - lambda) * static_cast<double>(Ntar));
    if (!(pi0 > 0.0)) { pi0 = 1.0 / static_cast<double>(Ntar); }   // never 0
    if (pi0 > 1.0) { pi0 = 1.0; }                                   // never > 1
  }

  std::size_t targets = 0;
  std::size_t decoys = 0;
  for (std::size_t begin = 0; begin < ranked.size();)
  {
    std::size_t end = begin + 1;
    while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
    for (std::size_t i = begin; i < end; ++i)
    {
      if (ranked[i].label == 1) { ++targets; }
      else                      { ++decoys; }
    }
    // FDR(t) = pi0 * (decoys_above/Ndec) / (targets_above/Ntar)  -- p-value based, ratio-corrected.
    double fdr;
    if (targets == 0) { fdr = std::numeric_limits<double>::infinity(); }
    else
    {
      // +1 finite-sample correction on the decoy count (Käll): keeps the estimator
      // honest at the tail where single decoys otherwise make it anti-conservative.
      const double dr = (Ndec > 0) ? (static_cast<double>(decoys) + 1.0) / static_cast<double>(Ndec) : 0.0;
      const double tr = static_cast<double>(targets) / static_cast<double>(Ntar);
      fdr = std::min(1.0, pi0 * dr / tr);
    }
    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = fdr; }
    begin = end;
  }

  double running_min = 1.0;
  for (std::size_t end = ranked.size(); end > 0;)
  {
    std::size_t begin = end - 1;
    while (begin > 0 && ranked[begin - 1].score == ranked[end - 1].score) { --begin; }
    running_min = std::min(running_min, ranked[begin].qvalue);
    for (std::size_t i = begin; i < end; ++i) { ranked[i].qvalue = running_min; }
    end = begin;
  }

  // --- p-value and PEP -------------------------------------------------------------------
  // These were previously never computed: the caller wrote the q-value into the PVALUE, QVALUE and
  // PEP columns alike, so 100% of rows had all three identical and anything downstream reading PEP
  // (IPF, protein-level inference) was silently consuming a q-value. They are different statistics:
  //
  //   p-value : TAIL probability under the null -- P(a null score >= this one), from the decoys.
  //   q-value : minimum FDR of the SET selected at this threshold (already computed above).
  //   PEP     : LOCAL FDR -- the probability that THIS peak group specifically is null. A group at
  //             q = 0.01 sitting right at the threshold can easily have PEP ~ 0.3.
  //
  // p is the standard conservative empirical estimate (Käll's +1 on both counts).
  {
    std::size_t dec_at_or_above = 0;
    for (std::size_t begin = 0; begin < ranked.size();)
    {
      std::size_t end = begin + 1;
      while (end < ranked.size() && ranked[end].score == ranked[begin].score) { ++end; }
      for (std::size_t i = begin; i < end; ++i) { if (ranked[i].label == 0) { ++dec_at_or_above; } }
      const double p = (Ndec > 0)
                         ? (static_cast<double>(dec_at_or_above) + 1.0) / (static_cast<double>(Ndec) + 1.0)
                         : 1.0;
      for (std::size_t i = begin; i < end; ++i) { ranked[i].pvalue = std::min(1.0, p); }
      begin = end;
    }
  }

  // PEP by the LOCAL analogue of the global estimator used above: in a score neighbourhood holding
  // t targets and d decoys, the decoys estimate the null target density (scaled by Ntar/Ndec and
  // pi0), so the expected null-target count is pi0*d*(Ntar/Ndec) and PEP ~ that divided by t.
  // Window is a fixed fraction of the list so it adapts to size; a raw local ratio is very noisy,
  // hence the monotonicity pass afterwards.
  {
    const std::size_t n = ranked.size();
    const std::size_t half = std::max<std::size_t>(50, n / 200);   // ~0.5% of the list, >=100 wide
    const double scale = (Ndec > 0) ? (static_cast<double>(Ntar) / static_cast<double>(Ndec)) : 0.0;
    // SLIDING counts, not a recount per element. The window is n/200 wide, so recomputing it inside
    // the loop is O(n^2/100): at n = 2M peak groups that is ~4e10 operations and turns a seconds-long
    // step into a multi-minute one. Maintaining running counts makes the whole pass O(n) -- each
    // element enters and leaves the window exactly once.
    std::size_t t = 0, d = 0, lo = 0, hi = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      const std::size_t new_lo = (i > half) ? i - half : 0;
      const std::size_t new_hi = std::min(n, i + half + 1);
      while (hi < new_hi) { if (ranked[hi].label == 1) { ++t; } else { ++d; } ++hi; }
      while (lo < new_lo) { if (ranked[lo].label == 1) { --t; } else { --d; } ++lo; }
      double pep = 1.0;
      if (t > 0) { pep = pi0 * static_cast<double>(d) * scale / static_cast<double>(t); }
      ranked[i].pep = std::min(1.0, std::max(0.0, pep));
    }
    // PEP must not increase with score. The list is sorted high->low, so sweep from the low-score
    // end keeping a running max: this is the isotonic projection under the ordering constraint and
    // removes the local-ratio noise without smoothing away the trend.
    double running_max = 0.0;
    for (std::size_t end = n; end > 0; --end)
    {
      running_max = std::max(running_max, ranked[end - 1].pep);
      ranked[end - 1].pep = running_max;
    }
  }
}

// Solve A x = b for symmetric positive-definite A. The input matrix is
// row-major and is replaced by its lower-triangular Cholesky factor.
inline bool choleskySolve(std::vector<double> a,
                          const std::vector<double>& b,
                          std::vector<double>& x)
{
  const std::size_t m = b.size();
  for (std::size_t i = 0; i < m; ++i)
  {
    for (std::size_t j = 0; j <= i; ++j)
    {
      double value = a[i * m + j];
      for (std::size_t k = 0; k < j; ++k) { value -= a[i * m + k] * a[j * m + k]; }
      if (i == j)
      {
        if (!(value > 0.0) || !std::isfinite(value)) { return false; }
        a[i * m + i] = std::sqrt(value);
      }
      else
      {
        a[i * m + j] = value / a[j * m + j];
      }
    }
  }

  std::vector<double> y(m, 0.0);
  for (std::size_t i = 0; i < m; ++i)
  {
    double value = b[i];
    for (std::size_t j = 0; j < i; ++j) { value -= a[i * m + j] * y[j]; }
    y[i] = value / a[i * m + i];
  }

  x.assign(m, 0.0);
  for (std::size_t ii = m; ii > 0; --ii)
  {
    const std::size_t i = ii - 1;
    double value = y[i];
    for (std::size_t j = i + 1; j < m; ++j) { value -= a[j * m + i] * x[j]; }
    x[i] = value / a[i * m + i];
    if (!std::isfinite(x[i])) { return false; }
  }
  return true;
}

} // namespace lda_detail

inline ScoredGroups scoreSemiSupervisedLDA(
  const std::vector<std::vector<double>>& features,
  const std::vector<int>& labels,
  const std::vector<long long>& group,
  const LDAParams& params)
{
  const std::size_t n = features.size();
  ScoredGroups result;
  result.dscore.assign(n, 0.0);
  result.qvalue.assign(n, 1.0);
  result.pvalue.assign(n, 1.0);
  result.pep.assign(n, 1.0);
  if (n == 0 || labels.size() != n || group.size() != n) { return result; }

  const std::size_t m = features.front().size();
  if (m == 0) { return result; }
  for (const auto& row : features)
  {
    if (row.size() != m) { return result; }
  }

  // Global z-standardisation is unsupervised. Non-finite cells are treated as
  // missing and become zero (the column mean) after standardisation.
  std::vector<double> mean(m, 0.0);
  std::vector<double> count(m, 0.0);
  for (const auto& row : features)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      if (std::isfinite(row[j]))
      {
        mean[j] += row[j];
        count[j] += 1.0;
      }
    }
  }
  for (std::size_t j = 0; j < m; ++j)
  {
    if (count[j] > 0.0) { mean[j] /= count[j]; }
  }

  std::vector<double> sum_squared_deviation(m, 0.0);
  std::vector<double> sd(m, 1.0);
  std::vector<std::vector<double>> z(n, std::vector<double>(m, 0.0));
  for (const auto& row : features)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      if (std::isfinite(row[j]))
      {
        const double d = row[j] - mean[j];
        sum_squared_deviation[j] += d * d;
      }
    }
  }
  for (std::size_t j = 0; j < m; ++j)
  {
    sd[j] = count[j] > 1.0
              ? std::sqrt(sum_squared_deviation[j] / (count[j] - 1.0))
              : 0.0;
    if (!(sd[j] > std::numeric_limits<double>::epsilon()) || !std::isfinite(sd[j]))
    {
      sd[j] = 1.0;
    }
  }
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      z[i][j] = std::isfinite(features[i][j]) ? (features[i][j] - mean[j]) / sd[j] : 0.0;
    }
  }

  // Build a stable, first-occurrence ordering of precursor groups.
  std::unordered_map<long long, std::size_t> group_lookup;
  group_lookup.reserve(n);
  std::vector<std::vector<std::size_t>> group_rows;
  std::vector<int> group_label;
  std::vector<long long> group_id;                 // the precursor id behind each group index
  for (std::size_t i = 0; i < n; ++i)
  {
    auto inserted = group_lookup.emplace(group[i], group_rows.size());
    if (inserted.second)
    {
      group_rows.emplace_back();
      group_label.push_back(labels[i] == 1 ? 1 : 0);
      group_id.push_back(group[i]);
    }
    const std::size_t g = inserted.first->second;
    group_rows[g].push_back(i);
  }

  const std::size_t group_count = group_rows.size();
  if (group_count < 2) { return result; } // leakage-free training is impossible

  int folds = params.n_folds;
  if (folds < 2) { folds = 2; }
  if (folds > static_cast<int>(group_count)) { folds = static_cast<int>(group_count); }

  std::vector<std::size_t> target_groups;
  std::vector<std::size_t> decoy_groups;
  for (std::size_t g = 0; g < group_count; ++g)
  {
    (group_label[g] == 1 ? target_groups : decoy_groups).push_back(g);
  }
  // Order these by the precursor's IDENTITY before shuffling. A group's index is its
  // first-occurrence position in row order, and row order is whatever the parallel extraction
  // happened to produce -- so a seeded shuffle of indices still put the same precursor in a
  // different fold on every run, training a different model and reporting a different ID count.
  // Measured spread on byte-identical input: 6487 / 6565 / 6433 (+/-1%), which is larger than
  // most of the effects being A/B tested. Sorting by id first costs one sort and makes the fold
  // assignment a function of the data alone; the shuffle still balances fold sizes exactly.
  const auto by_id = [&](std::size_t a, std::size_t b) { return group_id[a] < group_id[b]; };
  std::sort(target_groups.begin(), target_groups.end(), by_id);
  std::sort(decoy_groups.begin(), decoy_groups.end(), by_id);
  std::mt19937 rng(params.seed);
  std::shuffle(target_groups.begin(), target_groups.end(), rng);
  std::shuffle(decoy_groups.begin(), decoy_groups.end(), rng);
  std::vector<int> group_fold(group_count, 0);
  for (std::size_t i = 0; i < target_groups.size(); ++i)
  {
    group_fold[target_groups[i]] =
      static_cast<int>(i % static_cast<std::size_t>(folds));
  }
  for (std::size_t i = 0; i < decoy_groups.size(); ++i)
  {
    group_fold[decoy_groups[i]] =
      static_cast<int>(i % static_cast<std::size_t>(folds));
  }

  // Folds are independent BY CONSTRUCTION: fold f trains on the groups not assigned to f and writes
  // result.dscore only for the rows of groups that ARE assigned to f. So no two iterations read or
  // write the same dscore element, and none of them touch shared state except the two skip/train
  // counters, which are reduced. Everything else (`z`, `group_rows`, `group_fold`) is read-only here.
  // This loop was serial, and on the benchmark feature table the whole scoring step is 142 s of the
  // 1,284 s serial tail that caps the run's speedup at 1.90x -- see
  // docs/OpenDIAlyzer-parallel-efficiency.md. Cheap to fix, so fixed.
  //
  // Determinism is preserved: fold assignment comes from the seeded RNG above, each fold's model
  // depends only on its own training set, and the q-values are computed after the loop from the
  // pooled scores. The result does not depend on completion order.
  // The GBT parallelises internally over rows, but fit() is called from INSIDE the fold loop below.
  // A nested parallel region defaults to a team of one, so without raising the active-level limit
  // the inner pragmas would be dead code. Splitting the available threads as
  // folds x (threads/folds) uses the machine without oversubscribing it. The GBT's result does not
  // depend on this number (odia_gbt_test T8), so it is purely a speed knob.
  GBTParams gbt_params = params.gbt;
#ifdef _OPENMP
  if (params.classifier == LDAParams::Classifier::GBT)
  {
    omp_set_max_active_levels(2);
    gbt_params.n_threads = std::max(1, omp_get_max_threads() / std::max(1, folds));
  }
#endif

  int n_trained = 0, n_skipped = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : n_trained, n_skipped)
#endif
  for (int fold = 0; fold < folds; ++fold)
  {
    std::vector<std::size_t> train_rows;
    std::vector<std::size_t> train_groups;
    train_rows.reserve(n);
    train_groups.reserve(group_count);
    for (std::size_t g = 0; g < group_count; ++g)
    {
      if (!params.disable_cv && group_fold[g] == fold) { continue; }
      train_groups.push_back(g);
      train_rows.insert(train_rows.end(), group_rows[g].begin(), group_rows[g].end());
    }
    if (train_rows.empty()) { continue; }

    // w = S_W^-1 (mu_pos - mu_neg): Fisher's discriminant on two explicit row sets.
    auto fit_lda = [&](const std::vector<std::size_t>& positive_rows,
                       const std::vector<std::size_t>& negative_rows,
                       std::vector<double>& out) -> bool {
      if (positive_rows.size() < 2 || negative_rows.size() < 2) { return false; }
      std::vector<double> positive_mean(m, 0.0);
      std::vector<double> negative_mean(m, 0.0);
      for (const std::size_t row : positive_rows)
      {
        for (std::size_t j = 0; j < m; ++j) { positive_mean[j] += z[row][j]; }
      }
      for (const std::size_t row : negative_rows)
      {
        for (std::size_t j = 0; j < m; ++j) { negative_mean[j] += z[row][j]; }
      }
      for (std::size_t j = 0; j < m; ++j)
      {
        positive_mean[j] /= static_cast<double>(positive_rows.size());
        negative_mean[j] /= static_cast<double>(negative_rows.size());
      }

      std::vector<double> covariance(m * m, 0.0);
      auto add_within_class = [&](const std::vector<std::size_t>& rows,
                                  const std::vector<double>& class_mean) {
        for (const std::size_t row : rows)
        {
          for (std::size_t j = 0; j < m; ++j)
          {
            const double dj = z[row][j] - class_mean[j];
            for (std::size_t k = 0; k <= j; ++k)
            {
              covariance[j * m + k] += dj * (z[row][k] - class_mean[k]);
            }
          }
        }
      };
      add_within_class(positive_rows, positive_mean);
      add_within_class(negative_rows, negative_mean);
      const double dof =
        static_cast<double>(positive_rows.size() + negative_rows.size() - 2);
      for (std::size_t j = 0; j < m; ++j)
      {
        for (std::size_t k = 0; k <= j; ++k)
        {
          covariance[j * m + k] /= dof;
          covariance[k * m + j] = covariance[j * m + k];
        }
      }

      std::vector<double> difference(m);
      for (std::size_t j = 0; j < m; ++j)
      {
        difference[j] = positive_mean[j] - negative_mean[j];
      }

      // A positive user ridge is tried exactly; a zero/negative ridge starts
      // with a tiny numerical ridge. Increase it geometrically on failure.
      double ridge = params.ridge > 0.0 && std::isfinite(params.ridge)
                       ? params.ridge
                       : 1e-12;
      bool solved = false;
      for (int attempt = 0; attempt < 10 && !solved; ++attempt)
      {
        std::vector<double> regularised = covariance;
        for (std::size_t j = 0; j < m; ++j) { regularised[j * m + j] += ridge; }
        solved = lda_detail::choleskySolve(regularised, difference, out);
        ridge *= 10.0;
      }
      double norm_squared = 0.0;
      for (const double value : out) { norm_squared += value * value; }
      return solved && norm_squared > std::numeric_limits<double>::epsilon();
    };

    // The learner for this fold. LDA carries a weight vector; GBT carries a tree ensemble. Every
    // score in this fold goes through score_row(), so the choice is made in exactly ONE place and
    // the surrounding machinery -- training-set selection, fold normalisation, q-values -- is
    // literally the same code for both.
    GBT gbt;
    std::vector<NNEnsemble> nn_members;      // one per bag; size 1 when bagging is off
    const bool use_gbt = (params.classifier == LDAParams::Classifier::GBT);
    const bool use_nn = (params.classifier == LDAParams::Classifier::NN);


    // Initial direction: the signed feature with the largest absolute Welch
    // two-sample t statistic between all target and decoy training rows.
    std::vector<double> w(m, 0.0);
    auto score_row = [&](std::size_t row) -> double {
      if (use_nn && !nn_members.empty()) { return baggedScore(nn_members, z[row]); }
      return (use_gbt && gbt.trained()) ? gbt.score(z[row]) : lda_detail::dot(w, z[row]);
    };
    // Fit whichever learner this run selected, on the given rows. `mask` is honoured by the NN
    // only -- for a tree, a feature the fit never split on is already inert at scoring time, and
    // for LDA a zero weight is likewise inert, so neither needs one.
    auto fit_learner = [&](const std::vector<std::size_t>& pos, const std::vector<std::size_t>& neg,
                           const std::vector<char>& mask) -> bool {
      if (use_nn)
      {
        NNParams np = params.nn;
        np.mask = mask;
        if (params.bag_fraction >= 1.0)
        {
          NNEnsemble e;
          if (!e.fit(z, pos, neg, np)) { return false; }
          nn_members.assign(1, std::move(e));
          return true;
        }
        AnchorTrainingParams ap;
        ap.nn = np;
        ap.bag_fraction = params.bag_fraction;
        ap.min_anchors = 1;                  // the outer loop already refuses to fit on nothing
        ap.seed = params.seed;
        auto members = trainBaggedOnAnchors(z, pos, neg, group, ap);
        if (members.empty()) { return false; }
        nn_members = std::move(members);
        return true;
      }
      if (use_gbt) { GBT g; if (!g.fit(z, pos, neg, gbt_params)) { return false; }
                     gbt = std::move(g); return true; }
      std::vector<double> next_w;
      if (!fit_lda(pos, neg, next_w)) { return false; }
      w.swap(next_w);
      return true;
    };

    std::size_t best_feature = 0;
    double best_abs_t = -1.0;
    double best_difference = 1.0;
    for (std::size_t j = 0; j < m; ++j)
    {
      double sum[2] = {0.0, 0.0};
      double sum_sq[2] = {0.0, 0.0};
      std::size_t class_n[2] = {0, 0};
      for (const std::size_t row : train_rows)
      {
        const int cls = labels[row] == 1 ? 1 : 0;
        sum[cls] += z[row][j];
        sum_sq[cls] += z[row][j] * z[row][j];
        ++class_n[cls];
      }
      if (class_n[0] == 0 || class_n[1] == 0) { continue; }
      const double class_mean0 = sum[0] / static_cast<double>(class_n[0]);
      const double class_mean1 = sum[1] / static_cast<double>(class_n[1]);
      const double var0 = class_n[0] > 1
                            ? std::max(0.0, (sum_sq[0] - sum[0] * class_mean0) /
                                                static_cast<double>(class_n[0] - 1))
                            : 0.0;
      const double var1 = class_n[1] > 1
                            ? std::max(0.0, (sum_sq[1] - sum[1] * class_mean1) /
                                                static_cast<double>(class_n[1] - 1))
                            : 0.0;
      const double difference = class_mean1 - class_mean0;
      const double standard_error =
        std::sqrt(var0 / static_cast<double>(class_n[0]) +
                  var1 / static_cast<double>(class_n[1]));
      const double abs_t = standard_error > 0.0
                             ? std::abs(difference) / standard_error
                             : (difference == 0.0 ? 0.0
                                                  : std::numeric_limits<double>::infinity());
      if (abs_t > best_abs_t)
      {
        best_abs_t = abs_t;
        best_feature = j;
        best_difference = difference;
      }
    }
    w[best_feature] = best_difference < 0.0 ? -1.0 : 1.0;

    // ...but a SINGLE feature is a far weaker seed than it looks, and when it is too weak the
    // semi-supervised loop never ignites: iteration 0 selects positives with this w, finds none at
    // q<=train_fdr_initial, skips the fit, so w is unchanged and every later iteration skips too.
    // Measured on testdata/lda_fixture.txt (44,568 real OpenSWATH rows / 2,500 precursors): the best
    // single feature reaches |t|=13.7, which sounds decisive but over 29,482 rows is a 0.16 sd
    // per-row effect -- and with ~18 candidate peak groups per precursor, the max-over-group that
    // actually does the ranking is dominated by extreme-value noise. Result: a FLAT q of 0.965 for
    // every target group, 0/9 iterations trained, and 0 IDs.
    //
    // The fix is to seed with a real multivariate direction. The target/decoy labels are known
    // outright -- no FDR estimate needed -- so an LDA of all top-target rows against all top-decoy
    // rows is available for free and is enormously stronger than one column. The target class is
    // contaminated (most target precursors are false), which shrinks the fitted direction toward
    // zero but does not rotate it systematically: the contaminant IS the decoy distribution, so it
    // biases mu_pos toward mu_neg and costs magnitude, not orientation. It only has to be good
    // enough to ignite the loop; iteration 0 then re-selects positives properly.
    // pyprophet has no equivalent problem because OpenSWATH hands it a composite `main_score` to
    // rank by (semi_supervised.py: train.rank_by("main_score")). ODIA has no such column.
    {
      std::vector<std::size_t> seed_pos, seed_neg;
      seed_pos.reserve(train_groups.size());
      seed_neg.reserve(train_groups.size());
      for (const std::size_t g : train_groups)
      {
        std::size_t best_row = group_rows[g].front();
        double best = score_row(best_row);
        for (const std::size_t row : group_rows[g])
        {
          const double score = score_row(row);
          if (score > best) { best = score; best_row = row; }
        }
        (group_label[g] == 1 ? seed_pos : seed_neg).push_back(best_row);
      }
      // The seed must use the SAME learner as the loop. Seeding a GBT run with an LDA fit
      // reintroduces exactly the cold start this block exists to prevent, one level up: on data
      // whose signal is an interaction, the LDA seed is at chance by construction, so iteration 0
      // selects no positives, the GBT never fits, and the run silently falls back to ranking by a
      // hyperplane that cannot see the signal. Caught by odia_gbt_test T7, which reported 0 IDs
      // for BOTH classifiers before this.
      //
      // params.seed_mask applies HERE and only here (mechanism 1) -- the seed picks the anchors,
      // so it is the seed that must not see the calibrated features.
      fit_learner(seed_pos, seed_neg, params.seed_mask);
    }

    std::vector<std::size_t> prev_positives;   // mechanism 5's state, sorted
    for (int iteration = 0; iteration < std::max(0, params.n_iter); ++iteration)
    {
      // Reduce training scores to the best candidate row per precursor, then
      // estimate group-level q-values for confident positive selection.
      std::vector<lda_detail::RankedGroup> ranked;
      ranked.reserve(train_groups.size());
      for (const std::size_t g : train_groups)
      {
        std::size_t best_row = group_rows[g].front();
        double best_score = score_row(best_row);
        for (const std::size_t row : group_rows[g])
        {
          const double score = score_row(row);
          if (score > best_score)
          {
            best_score = score;
            best_row = row;
          }
        }
        ranked.push_back({g, best_row, group_label[g], best_score, 1.0});
      }
      lda_detail::assignQValues(ranked, params.use_pi0);

      // First pass is seeded by a SINGLE feature, so a strict cut can select too few positives to
      // fit anything and the model never bootstraps -- the failure this split exists to prevent.
      const double raw_fdr = (iteration == 0) ? params.train_fdr_initial : params.train_fdr;
      const double train_fdr = std::max(0.0, std::min(1.0, raw_fdr));
      std::vector<std::size_t> positive_rows;
      for (const auto& candidate : ranked)
      {
        if (candidate.label == 1 && candidate.qvalue <= train_fdr)
        {
          positive_rows.push_back(candidate.best_row);
        }
      }
      // Too few confident positives to fit an m-dimensional discriminant. Skipping is right, but it
      // used to be SILENT -- and silence here is dangerous: if every iteration skips, `w` stays at
      // its initialisation (a single feature, weight +/-1), so the "LDA" degenerates to ranking by
      // one sub-score and nothing says so. Record it; the caller reports it.
      if (positive_rows.size() < m + 2)
      {
        ++n_skipped;
        continue;
      }
      ++n_trained;

      // The positive class above is one row per precursor -- `candidate.best_row`, the highest
      // scoring peak group. Taking ALL decoy rows as negatives would therefore make the two classes
      // differ in two ways at once: target vs decoy (wanted) AND rank-1 vs runner-up (not wanted).
      // With ~4-6 candidate peak groups per precursor the negative class is then dominated by
      // runner-ups, so the fitted direction partly separates "best peak in its group" from "not the
      // best peak" -- a real, learnable axis that carries no target/decoy information. Both the
      // negative mean and the within-class covariance are pulled by it. Matching the positive side
      // (top peak per precursor) removes the confound; this is what pyprophet does
      // (semi_supervised.py: td_peaks = train.get_top_decoy_peaks()).
      std::vector<std::size_t> negative_rows;
      for (const std::size_t g : train_groups)
      {
        if (group_label[g] == 1) { continue; }
        if (!params.top_decoys_only)
        {
          negative_rows.insert(negative_rows.end(), group_rows[g].begin(), group_rows[g].end());
          continue;
        }
        std::size_t best_row = group_rows[g].front();
        double best_score = score_row(best_row);
        for (const std::size_t row : group_rows[g])
        {
          const double score = score_row(row);
          if (score > best_score) { best_score = score; best_row = row; }
        }
        negative_rows.push_back(best_row);
      }
      if (negative_rows.size() < 2) { continue; }

      // A failed fit leaves the previous model in place, exactly as a failed Cholesky leaves the
      // previous w -- the iteration is skipped, not replaced with something degenerate.
      fit_learner(positive_rows, negative_rows, {});

      // Mechanism 5. Watch the positive SET, not the score. A model collapsing onto a subset of
      // its seed shows a SHRINKING positive set while its identification count rises, so a rule
      // reading the score is blind to exactly the failure worth catching.
      if (params.stop_on_composition && iteration >= 1)
      {
        // FROM ITERATION 1 ONWARD, NOT FROM 0. Iteration 0 selects positives at
        // train_fdr_initial (0.15) and every later iteration at train_fdr (0.05) -- a 3x stricter
        // cut. Comparing across that change makes a HEALTHY run look collapsed at the first
        // opportunity: the set is smaller because the threshold moved, not because the model
        // narrowed. With the rule armed that way it broke out after a single iteration every time,
        // so the mechanism meant to DETECT a collapse instead silently truncated training.
        // Only sets selected at the SAME threshold are comparable.
        std::vector<std::size_t> curr = positive_rows;
        std::sort(curr.begin(), curr.end());
        AnchorTrainingParams ap;
        ap.stop_jaccard = params.stop_jaccard;
        ap.max_iterations = params.n_iter;
        AnchorTrainingReport rep;
        rep.iterations_run = iteration;
        // The production path CALLS the tested helper rather than reimplementing it. The inline
        // copy that used to live here was the untested twin of a tested function, which is how the
        // threshold bug above survived its own unit test.
        const bool go = anchorIterationShouldContinue(prev_positives, curr, ap, rep);
        prev_positives.swap(curr);
        if (!go) { break; }
      }
      else if (params.stop_on_composition)
      {
        prev_positives = positive_rows;
        std::sort(prev_positives.begin(), prev_positives.end());
      }
    }

    // This model has seen no row from the groups scored.
    for (std::size_t g = 0; g < group_count; ++g)
    {
      if (group_fold[g] != fold) { continue; }
      for (const std::size_t row : group_rows[g])
      {
        result.dscore[row] = score_row(row);
      }
    }

    // Every fold has its OWN weight vector, and an LDA direction is defined only up to scale (and
    // the mean-offset depends on that fold's training set), so fold A's score of 4.0 and fold B's
    // of 4.0 mean different things. The final q-values below pool all folds into ONE ranking, which
    // silently assumes they are commensurable. They are not, and the pooled ranking is then partly
    // sorted by which fold a precursor happened to land in. Rescaling each fold's held-out scores to
    // that fold's own DECOY null (mean 0, sd 1) puts every fold on the one scale target-decoy FDR
    // actually cares about -- "how many decoy sds above the null" -- and makes pooling valid.
    // pyprophet sidesteps the problem differently: it averages the fold weight vectors into a single
    // model and rescores everything with it (at the cost of the leakage-free property kept here).
    if (params.normalize_folds)
    {
      double sum = 0.0, sum_sq = 0.0;
      std::size_t decoy_n = 0;
      for (std::size_t g = 0; g < group_count; ++g)
      {
        if (group_fold[g] != fold || group_label[g] == 1) { continue; }
        std::size_t best_row = group_rows[g].front();
        for (const std::size_t row : group_rows[g])
        {
          if (result.dscore[row] > result.dscore[best_row]) { best_row = row; }
        }
        sum += result.dscore[best_row];
        sum_sq += result.dscore[best_row] * result.dscore[best_row];
        ++decoy_n;
      }
      if (decoy_n >= 2)
      {
        const double mu = sum / static_cast<double>(decoy_n);
        const double var = std::max(0.0, (sum_sq - sum * mu) / static_cast<double>(decoy_n - 1));
        const double sigma = std::sqrt(var);
        // A degenerate null (all decoys identical) carries no scale information; leaving those
        // scores unscaled is the only honest option, and shifting them alone would be worse.
        if (sigma > std::numeric_limits<double>::epsilon() && std::isfinite(sigma))
        {
          for (std::size_t g = 0; g < group_count; ++g)
          {
            if (group_fold[g] != fold) { continue; }
            for (const std::size_t row : group_rows[g])
            {
              result.dscore[row] = (result.dscore[row] - mu) / sigma;
            }
          }
        }
      }
    }
  }

  result.n_iterations_trained = n_trained;
  result.n_iterations_skipped = n_skipped;

  std::vector<lda_detail::RankedGroup> final_ranked;
  final_ranked.reserve(group_count);
  for (std::size_t g = 0; g < group_count; ++g)
  {
    std::size_t best_row = group_rows[g].front();
    for (const std::size_t row : group_rows[g])
    {
      if (result.dscore[row] > result.dscore[best_row]) { best_row = row; }
    }
    final_ranked.push_back(
      {g, best_row, group_label[g], result.dscore[best_row], 1.0});
  }
  lda_detail::assignQValues(final_ranked, params.use_pi0);
  for (const auto& ranked_group : final_ranked)
  {
    for (const std::size_t row : group_rows[ranked_group.group_index])
    {
      result.qvalue[row] = ranked_group.qvalue;
      result.pvalue[row] = ranked_group.pvalue;
      result.pep[row]    = ranked_group.pep;
    }
  }
  return result;
}

} // namespace odia

#endif // ODIA_LDA_H

#pragma once
// A LEARNED library pre-filter, replacing the fixed "N of the top-6 fragments in one spectrum" rule.
//
// WHY. The fixed rule thresholds ONE feature at an integer. Measured against a reference ID list
// over 3,603,425 target candidates, the enrichment of true identifications by match depth is:
//
//     depth 6 : 3479 refs /   4,768 targets = 729.7 per 1000
//     depth 5 : 1399 refs /   8,059 targets = 173.6
//     depth 4 : 1175 refs / 107,686 targets =  10.9
//     depth 3 : 1084 refs / 914,675 targets =   1.19
//     depth 2 :  636 refs / 1,645,520       =   0.39
//
// A 600x dynamic range collapsed into one bit. Depth 6 is 67x more informative than depth 4 and the
// threshold treats them identically. Worse, the threshold sits on a cliff: 4 -> 3 is a 9x drop in
// enrichment, so relaxing it by one admits 913,591 near-certain noise candidates for 1,084 real
// ones and measured 6930 -> 5580 identifications.
//
// The evidence the filter already computes is richer than the bit it emits: match depth, total
// matches, how many spectra qualified, MS1 hits, and matched intensities. This scores a precursor
// on all of it and keeps the best N, instead of keeping everything above an integer.
//
// ---------------------------------------------------------------------------------------------
// VALIDITY. This is selection on data, upstream of a target-decoy FDR, so it can invalidate the
// null if done carelessly. Three rules, all load-bearing:
//
// 1. LABEL SYMMETRY. The model must be applied to targets and decoys through identical code, and
//    the same number retained from each class. The existing threshold already has this property --
//    the surrounding code scans a decoy view with the same criterion and keeps a pair when either
//    member passes -- and it exists because a target set that cleared an evidence bar, paired with
//    a decoy set that did not, makes retained decoys score systematically lower than retained
//    false targets. That is an anti-conservative FDR.
//
// 2. NO SEQUENCE-DERIVED FEATURES. Decoys are shuffled sequences, so anything reading composition
//    separates the classes for a reason unrelated to presence in the sample, and the model would
//    learn "is this a decoy" rather than "is this present". Precursor m/z and charge ARE safe: a
//    decoy preserves its target's mass by construction, so they carry no label information.
//
// 3. THE SELECTION IS CONSERVATIVE BY CONSTRUCTION, and it is worth understanding why. Training on
//    target-vs-decoy and keeping the top-N of each class retains the most TARGET-LIKE decoys --
//    the hardest negatives. That makes the downstream null harder, not easier. The failure mode to
//    watch is the opposite one: keeping only hard cases could compress the score range the
//    downstream classifier sees. Measure identifications, not just the retained-set purity.
//
// This is deliberately NOT trained on identification outcomes. Training on which precursors were
// identified would select on the downstream test statistic and is how a pre-filter turns into an
// FDR violation.

#include "odia_gbt.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

namespace odia
{

/// The evidence a pre-filter scan produces for ONE precursor candidate. Mirrors the fields an
/// evidence filter already computes, so no new measurement is required.
struct PrefilterEvidence
{
  double precursor_mz = 0.0;
  int charge = 0;
  std::size_t ms2_best_fragment_hits = 0;   ///< deepest co-occurrence in a single spectrum
  std::size_t ms2_hit_count = 0;            ///< total distinct matches across spectra
  std::size_t ms2_qualifying_spectra = 0;   ///< how many spectra met the depth bar
  double ms2_max_intensity = 0.0;
  double ms2_sum_intensity = 0.0;
  std::size_t ms1_hit_count = 0;
  double ms1_max_intensity = 0.0;
  double ms1_sum_intensity = 0.0;
  std::size_t n_indexed_transitions = 6;    ///< denominator for the depth fraction
};

/// Feature names, in the order featurise() emits them. Exposed so a caller can log which inputs the
/// model actually had -- a silently-missing column is otherwise invisible.
inline std::vector<std::string> prefilterFeatureNames()
{
  return {"log_best_hits", "depth_fraction", "log_hit_count", "log_qualifying_spectra",
          "log_ms2_max_int", "log_ms2_sum_int", "ms2_int_concentration",
          "log_ms1_hits", "log_ms1_max_int", "has_ms1",
          "precursor_mz_scaled", "charge"};
}

/// Map one candidate's evidence to a fixed-width feature row.
///
/// Counts are log1p-compressed: the enrichment table above spans 600x, and a tree learns a
/// monotone split either way, but the compression keeps the histogram binning from spending all
/// 64 bins on the tail. Intensities likewise.
inline void featurise(const PrefilterEvidence& e, std::vector<double>& out)
{
  const double denom = static_cast<double>(std::max<std::size_t>(1, e.n_indexed_transitions));
  out.clear();
  out.push_back(std::log1p(static_cast<double>(e.ms2_best_fragment_hits)));
  out.push_back(static_cast<double>(e.ms2_best_fragment_hits) / denom);
  out.push_back(std::log1p(static_cast<double>(e.ms2_hit_count)));
  out.push_back(std::log1p(static_cast<double>(e.ms2_qualifying_spectra)));
  out.push_back(std::log1p(e.ms2_max_intensity));
  out.push_back(std::log1p(e.ms2_sum_intensity));
  // How concentrated the matched signal is in its strongest peak. A real fragment set spreads
  // intensity across several channels; a single coincident spike does not.
  out.push_back(e.ms2_sum_intensity > 0.0 ? e.ms2_max_intensity / e.ms2_sum_intensity : 0.0);
  out.push_back(std::log1p(static_cast<double>(e.ms1_hit_count)));
  out.push_back(std::log1p(e.ms1_max_intensity));
  out.push_back(e.ms1_hit_count > 0 ? 1.0 : 0.0);
  out.push_back(e.precursor_mz / 1000.0);
  out.push_back(static_cast<double>(e.charge));
}

struct PrefilterModelParams
{
  GBTParams gbt;                    ///< the boosted-tree settings; deterministic, no subsampling
  double keep_fraction = 0.06;      ///< share of each class retained; 0.06 ~ the fixed rule's 5.9%
  std::size_t min_train_rows = 2000;///< below this, fall back rather than fit noise
};

/// A fitted model. Selection is a SET operation (see selectPrefilter), not a per-row threshold,
/// so no cut is stored: GBT scores are heavily tied -- integer match depths dominate the features
/// and thousands of candidates share a score -- and a `>= cut` rule keeps every tied row. Asking
/// for 6% that way measured 28% of targets and 15% of decoys retained, which is both wrong and
/// ASYMMETRIC, the one thing the null cannot tolerate.
struct PrefilterModel
{
  GBT gbt;
  bool trained = false;
  std::vector<std::string> feature_names;

  double score(const PrefilterEvidence& e) const
  {
    std::vector<double> row;
    featurise(e, row);
    return gbt.score(row);
  }
};

/// Train on evidence with known labels (decoy = known negative).
///
/// Returns trained=false if there is too little data, in which case the caller must fall back to
/// the fixed rule rather than proceed with an unfitted model.
inline PrefilterModel fitPrefilterModel(const std::vector<PrefilterEvidence>& ev,
                                        const std::vector<bool>& is_decoy,
                                        const PrefilterModelParams& p = {})
{
  PrefilterModel m;
  m.feature_names = prefilterFeatureNames();
  if (ev.size() != is_decoy.size() || ev.size() < p.min_train_rows) { return m; }

  std::vector<std::vector<double>> X;
  X.reserve(ev.size());
  std::vector<double> row;
  for (const auto& e : ev) { featurise(e, row); X.push_back(row); }

  std::vector<std::size_t> pos, neg;
  for (std::size_t i = 0; i < is_decoy.size(); ++i) { (is_decoy[i] ? neg : pos).push_back(i); }
  if (pos.size() < p.min_train_rows / 2 || neg.size() < p.min_train_rows / 2) { return m; }
  if (!m.gbt.fit(X, pos, neg, p.gbt)) { return m; }
  m.trained = true;
  return m;
}

/// Choose which candidates to keep: the top `keep_fraction` of EACH class by model score.
///
/// Per-class ranking, not one shared threshold. Targets score higher than decoys by construction,
/// so a shared cut retains more targets than decoys and hands the downstream FDR an unbalanced
/// null. Ranking each class separately makes the retained COUNTS proportional, which is the
/// property that matters.
///
/// Ties are broken by candidate INDEX -- library order, which is a property of the input file and
/// not of the run. Breaking them by score alone is impossible (they are equal) and breaking them
/// randomly would make the pre-filter non-reproducible.
///
/// Returns an all-true mask when the model is untrained: failing OPEN is the only safe direction,
/// since a filter that silently rejects everything looks exactly like one that finds nothing.
inline std::vector<char> selectPrefilter(const PrefilterModel& m,
                                         const std::vector<PrefilterEvidence>& ev,
                                         const std::vector<bool>& is_decoy,
                                         double keep_fraction)
{
  std::vector<char> keep(ev.size(), 1);
  if (!m.trained || ev.size() != is_decoy.size()) { return keep; }
  std::fill(keep.begin(), keep.end(), 0);

  for (int cls = 0; cls < 2; ++cls)
  {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < ev.size(); ++i)
    {
      if (static_cast<int>(is_decoy[i]) == cls) { idx.push_back(i); }
    }
    if (idx.empty()) { continue; }
    std::vector<double> sc(idx.size());
    for (std::size_t j = 0; j < idx.size(); ++j) { sc[j] = m.score(ev[idx[j]]); }

    std::vector<std::size_t> order(idx.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      if (sc[a] != sc[b]) { return sc[a] > sc[b]; }
      return idx[a] < idx[b];                       // library order: deterministic tie-break
    });
    const std::size_t n_keep = static_cast<std::size_t>(
      std::llround(static_cast<double>(idx.size()) * keep_fraction));
    for (std::size_t j = 0; j < n_keep && j < order.size(); ++j) { keep[idx[order[j]]] = 1; }
  }
  return keep;
}

} // namespace odia

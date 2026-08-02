#pragma once
// Training the peak-group classifier on the run's own confident anchors, robustly.
//
// THE PROBLEM. A semi-supervised classifier needs positives, and the only positives available are
// the ones it identified itself. That is circular by construction, and the circularity has teeth:
// the seed model is weak, so whatever it happens to prefer becomes the training set, and every
// subsequent iteration reinforces that preference. A model that overfits its seed looks *better*
// by every internal measure -- target/decoy separation improves, q-values shrink, identification
// counts rise -- while identifying a narrower and more self-similar set of peptides.
//
// The anchors used for RT recalibration are the natural seed: they are the highest-confidence
// identifications the run produces, they already exist, and they cost nothing extra. But they are
// selected BY a d-score, so using them naively closes exactly the loop described above.
//
// This is the protocol that breaks it. Five mechanisms, each addressing a distinct failure.
//
// ---------------------------------------------------------------------------------------------
// 1. FEATURE EXCLUSION AT SEED TIME (anti-circularity for the recalibrated quantity)
//
// The anchors exist to correct RT. If they are selected using RT-agreement sub-scores, they are
// selected by the very agreement they are meant to measure, and the correction they produce is
// biased toward the uncorrected state. So the SEED scoring pass excludes the RT-derived features,
// and only the full feature set is used once the transform is fitted.
//
// This project already applies exactly this rule to its LDA anchor selection. The same rule must
// carry over to the network, and it generalises: exclude, at seed time, any feature the anchors
// will be used to calibrate.
//
// 2. LENIENT POSITIVE SELECTION (keeps the seed broad)
//
// Take the confident set at q <= 0.10, not 0.01. DIA-NN uses 0.10 for exactly this purpose. A
// strict cutoff yields a small, homogeneous seed -- the easiest peptides, which are the least
// informative about the boundary -- and the model learns the easy cases and nothing else. The
// looser threshold admits some false positives, which is the price for a seed that spans the
// decision boundary rather than hugging one side of it.
//
// 3. BAGGING OVER ANCHORS (the ensemble must disagree about the seed, not just about init)
//
// DIA-NN's ensemble members differ only in weight initialisation, so they see identical training
// data and share identical seed bias -- averaging them reduces variance but not that bias. Here
// each member trains on a deterministic 70% subsample of the anchors. Members then disagree about
// which anchors mattered, and the disagreement is measurable: `anchorDisagreement()` returns the
// spread of member scores, and a spread that collapses toward zero means the ensemble has
// converged onto the seed rather than generalising from it.
//
// 4. GROUP-WISE CROSS-VALIDATION (unchanged, and still load-bearing)
//
// Folds are assigned by precursor group, so a peak group is never scored by a model trained on its
// own precursor. Without this the circularity is total: a row's own presence in the training set
// guarantees its score.
//
// 5. A STOPPING RULE THAT WATCHES THE SEED, NOT THE SCORE
//
// Iterate while the positive set is still CHANGING in composition, and stop when it stabilises or
// starts shrinking. A growing-then-shrinking positive set is the signature of a model collapsing
// onto a subset. Counting identifications cannot detect this -- they rise throughout.
//
// ---------------------------------------------------------------------------------------------
// WHAT THIS SCHEME CANNOT DO. None of the above detects a seed that is *uniformly* biased -- if
// the initial scoring systematically prefers, say, high-intensity precursors, every mechanism here
// preserves that preference while making it more confident. Only an external standard detects
// that: an entrapment set, or a second run. Treat the internal diagnostics as necessary and not
// sufficient.

#include "odia_nn.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

namespace odia
{

struct AnchorTrainingParams
{
  NNParams nn;
  double seed_qvalue = 0.10;        ///< lenient, per mechanism 2
  double bag_fraction = 0.70;       ///< share of anchors each ensemble member sees (mechanism 3)
  int max_iterations = 4;
  double stop_jaccard = 0.98;       ///< stop when consecutive positive sets overlap this much
  std::size_t min_anchors = 200;    ///< below this, do not train -- fall back to the linear path
  std::uint64_t seed = 42;
};

struct AnchorTrainingReport
{
  int iterations_run = 0;
  std::vector<std::size_t> positive_counts;   ///< per iteration; a fall is the collapse signature
  std::vector<double> jaccard;                ///< overlap with the previous iteration's positives
  double final_disagreement = 0.0;            ///< ensemble spread on the anchors
  bool converged = false;
  bool collapsed = false;                     ///< positive set shrank -- see mechanism 5
  std::string note;
};

/// Deterministic membership test for bagging. A hash of (member, group id), never an RNG: the bag
/// must be identical on every machine and every thread count, and it must be keyed on the PRECURSOR
/// GROUP rather than the row, so all rows of a precursor land in the same bag together.
inline bool inBag(std::uint64_t member, std::int64_t group, std::uint64_t seed, double fraction)
{
  std::uint64_t h = seed ^ (member * 0x9E3779B97F4A7C15ull);
  h ^= static_cast<std::uint64_t>(group) + 0xBF58476D1CE4E5B9ull + (h << 6) + (h >> 2);
  h *= 0x94D049BB133111EBull;
  h ^= h >> 31;
  return static_cast<double>(h >> 11) / static_cast<double>(1ull << 53) < fraction;
}

/// Jaccard overlap between two sorted index sets. Used to detect a stabilised (or collapsing)
/// positive set, which is the quantity mechanism 5 watches.
inline double jaccardOverlap(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b)
{
  if (a.empty() && b.empty()) { return 1.0; }
  std::size_t inter = 0, i = 0, j = 0;
  while (i < a.size() && j < b.size())
  {
    if (a[i] == b[j]) { ++inter; ++i; ++j; }
    else if (a[i] < b[j]) { ++i; }
    else { ++j; }
  }
  const std::size_t uni = a.size() + b.size() - inter;
  return uni ? static_cast<double>(inter) / static_cast<double>(uni) : 1.0;
}

/// Spread of ensemble member scores on a sample of rows. Near zero means the members agree
/// completely, which after bagging means they have converged onto the seed rather than
/// generalising from it -- mechanism 3's diagnostic.
inline double anchorDisagreement(const std::vector<NNEnsemble>& members,
                                 const std::vector<std::vector<double>>& X,
                                 const std::vector<std::size_t>& rows)
{
  if (members.size() < 2 || rows.empty()) { return 0.0; }
  double acc = 0.0;
  for (std::size_t r : rows)
  {
    double m = 0.0, m2 = 0.0;
    for (const auto& e : members)
    {
      const double s = e.score(X[r]);
      m += s;
      m2 += s * s;
    }
    const double n = static_cast<double>(members.size());
    const double var = std::max(0.0, m2 / n - (m / n) * (m / n));
    acc += std::sqrt(var);
  }
  return acc / static_cast<double>(rows.size());
}

/// Train an ensemble on anchors, with bagging over precursor groups.
///
/// `positives` are row indices of the current confident set; `negatives` are decoy rows; `group`
/// gives each row's precursor id so a bag takes whole precursors rather than splitting them.
///
/// Returns one trained ensemble per member. They are NOT averaged here -- the caller averages, and
/// can first inspect their disagreement, which is the point of bagging them separately.
template <typename GroupT>
inline std::vector<NNEnsemble> trainBaggedOnAnchors(const std::vector<std::vector<double>>& X,
                                                    const std::vector<std::size_t>& positives,
                                                    const std::vector<std::size_t>& negatives,
                                                    const std::vector<GroupT>& group,
                                                    const AnchorTrainingParams& p)
{
  std::vector<NNEnsemble> members;
  if (positives.size() < p.min_anchors) { return members; }

  const int n_members = std::max(1, p.nn.n_nets);
  members.reserve(static_cast<std::size_t>(n_members));
  for (int m = 0; m < n_members; ++m)
  {
    std::vector<std::size_t> pos_bag, neg_bag;
    pos_bag.reserve(positives.size());
    neg_bag.reserve(negatives.size());
    for (std::size_t r : positives)
    {
      if (inBag(static_cast<std::uint64_t>(m), static_cast<std::int64_t>(group[r]), p.seed, p.bag_fraction))
      {
        pos_bag.push_back(r);
      }
    }
    // Negatives are bagged too, at the same rate, so the class ratio each member sees matches the
    // full set. Bagging only the positives would make every member's null larger than its signal
    // by a different factor.
    for (std::size_t r : negatives)
    {
      if (inBag(static_cast<std::uint64_t>(m), static_cast<std::int64_t>(group[r]),
                p.seed ^ 0xABCDEFull, p.bag_fraction))
      {
        neg_bag.push_back(r);
      }
    }
    if (pos_bag.empty() || neg_bag.empty()) { continue; }

    NNParams np = p.nn;
    np.n_nets = 1;                      // one net per bag; the bagging IS the ensemble
    np.seed = p.nn.seed + static_cast<std::uint64_t>(m);
    NNEnsemble e;
    if (e.fit(X, pos_bag, neg_bag, np)) { members.push_back(std::move(e)); }
  }
  return members;
}

/// Mean score across bagged members.
inline double baggedScore(const std::vector<NNEnsemble>& members, const std::vector<double>& x)
{
  if (members.empty()) { return 0.0; }
  double s = 0.0;
  for (const auto& e : members) { s += e.score(x); }
  return s / static_cast<double>(members.size());
}

/// Decide whether to continue iterating, given the positive sets from the last two rounds.
/// Fills `rep` and returns true to continue.
///
/// Stops on stabilisation (Jaccard >= stop_jaccard) or on COLLAPSE (the positive set shrank). The
/// second is the one worth having: identification counts rise throughout a collapse, so a stopping
/// rule that watches the score cannot see it.
inline bool anchorIterationShouldContinue(const std::vector<std::size_t>& prev,
                                          const std::vector<std::size_t>& curr,
                                          const AnchorTrainingParams& p,
                                          AnchorTrainingReport& rep)
{
  // Incremented HERE, by the function that reads it. It was previously only ever set by callers,
  // so `iterations_run < max_iterations` read 0 < 4 forever and the iteration cap never fired --
  // the one stopping condition that does not depend on the data was dead.
  ++rep.iterations_run;
  rep.positive_counts.push_back(curr.size());
  if (prev.empty()) { rep.jaccard.push_back(0.0); return rep.iterations_run < p.max_iterations; }

  const double j = jaccardOverlap(prev, curr);
  rep.jaccard.push_back(j);
  if (curr.size() < prev.size())
  {
    rep.collapsed = true;
    rep.note = "positive set shrank (" + std::to_string(prev.size()) + " -> " +
               std::to_string(curr.size()) + "); model is converging onto a subset of its seed";
    return false;
  }
  if (j >= p.stop_jaccard)
  {
    rep.converged = true;
    rep.note = "positive set stabilised at Jaccard " + std::to_string(j);
    return false;
  }
  return rep.iterations_run < p.max_iterations;
}

} // namespace odia

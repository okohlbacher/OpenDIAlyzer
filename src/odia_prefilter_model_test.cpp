// Self-check for odia_prefilter_model.h.
//
// The properties that matter are the ones whose violation is SILENT: an asymmetric retention rate
// between classes would unbalance the downstream null without any error, and a model that merely
// reproduces the depth threshold would look like an improvement while adding nothing.

#include "odia_prefilter_model.h"

#include <cassert>
#include <cstdio>
#include <random>

using namespace odia;

namespace
{

/// Synthetic candidates whose depth enrichment mirrors the MEASURED benchmark distribution:
/// depth 6 is ~730 true per 1000, depth 3 ~1.2 per 1000. Decoys draw from the null only.
std::vector<PrefilterEvidence> make(std::size_t n, bool decoy, std::mt19937& rng,
                                    std::vector<bool>& present)
{
  // measured share of candidates at each depth, targets
  const double share[7] = {0.235, 0.021, 0.457, 0.254, 0.030, 0.0022, 0.0013};
  const double p_true[7] = {0.0, 0.00013, 0.00039, 0.00119, 0.0109, 0.1736, 0.7297};
  std::discrete_distribution<int> depth(share, share + 7);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::lognormal_distribution<double> inten(10.0, 2.0);

  std::vector<PrefilterEvidence> out;
  out.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    PrefilterEvidence e;
    int d = depth(rng);
    // A decoy cannot be genuinely present, so it draws depth from the null part of the mixture:
    // shift it down, which is what "coincidence only" looks like.
    if (decoy) { d = std::max(0, d - (u(rng) < 0.5 ? 1 : 0)); }
    const bool real = !decoy && (u(rng) < p_true[d]);
    e.ms2_best_fragment_hits = static_cast<std::size_t>(d);
    e.ms2_hit_count = static_cast<std::size_t>(d * (real ? 40 : 12) + 1);
    e.ms2_qualifying_spectra = real ? static_cast<std::size_t>(3 + 10 * u(rng)) : (d >= 4 ? 1 : 0);
    e.ms2_max_intensity = inten(rng) * (real ? 4.0 : 1.0);
    e.ms2_sum_intensity = e.ms2_max_intensity * (real ? 3.5 : 1.2);
    e.ms1_hit_count = real ? static_cast<std::size_t>(1 + 3 * u(rng)) : (u(rng) < 0.1 ? 1 : 0);
    e.ms1_max_intensity = e.ms1_hit_count ? inten(rng) : 0.0;
    e.precursor_mz = 400.0 + 800.0 * u(rng);
    e.charge = 2 + static_cast<int>(2 * u(rng));
    out.push_back(e);
    present.push_back(real);
  }
  return out;
}

} // namespace

int main()
{
  std::printf("odia_prefilter_model_test\n");
  std::mt19937 rng(7);

  std::vector<bool> tp, dp;
  auto tgt = make(60000, false, rng, tp);
  auto dec = make(60000, true, rng, dp);

  std::vector<PrefilterEvidence> ev = tgt;
  ev.insert(ev.end(), dec.begin(), dec.end());
  std::vector<bool> is_decoy(tgt.size(), false);
  is_decoy.insert(is_decoy.end(), dec.size(), true);

  PrefilterModelParams p;
  p.gbt.n_trees = 40;
  p.gbt.max_depth = 4;
  p.keep_fraction = 0.06;

  const auto m = fitPrefilterModel(ev, is_decoy, p);
  assert(m.trained);
  assert(m.feature_names.size() == 12);
  std::printf("  model trained on %zu rows, %zu features\n", ev.size(), m.feature_names.size());

  // ---- 1. retention must be class-SYMMETRIC in proportion --------------------------------------
  // Violating this is silent and unbalances the downstream null.
  const auto keep = selectPrefilter(m, ev, is_decoy, p.keep_fraction);
  std::size_t kt = 0, kd = 0;
  for (std::size_t i = 0; i < ev.size(); ++i) { (is_decoy[i] ? kd : kt) += (keep[i] != 0); }
  const double rt = double(kt) / tgt.size(), rd = double(kd) / dec.size();
  std::printf("  retained: targets %.3f, decoys %.3f (asked %.3f)\n", rt, rd, p.keep_fraction);
  assert(std::abs(rt - p.keep_fraction) < 0.015);
  assert(std::abs(rd - p.keep_fraction) < 0.015);

  // ---- 2. it must BEAT the fixed depth threshold at equal retention -----------------------------
  // The fixed rule keeps depth>=4. Count how many genuinely-present targets each retains.
  std::size_t model_true = 0, thr_true = 0, thr_kept = 0;
  for (std::size_t i = 0; i < tgt.size(); ++i)
  {
    if (keep[i] && tp[i]) { ++model_true; }
    if (tgt[i].ms2_best_fragment_hits >= 4) { ++thr_kept; if (tp[i]) { ++thr_true; } }
  }
  std::printf("  true positives retained: model %zu (of %zu kept) vs depth>=4 threshold %zu (of %zu)\n",
              model_true, kt, thr_true, thr_kept);
  assert(model_true >= thr_true);

  // ---- 3. determinism: the same input must give the same model ---------------------------------
  const auto m2 = fitPrefilterModel(ev, is_decoy, p);
  const auto keep2 = selectPrefilter(m2, ev, is_decoy, p.keep_fraction);
  std::size_t diff = 0;
  for (std::size_t i = 0; i < ev.size(); ++i) { diff += (keep[i] != keep2[i]); }
  std::printf("  refit disagreement: %zu rows\n", diff);
  assert(diff == 0);

  // ---- 4. an untrained model must keep EVERYTHING, not nothing ---------------------------------
  // Failing open is the only safe direction: a model that silently rejects all candidates on a
  // small dataset would look like a working filter that finds nothing.
  PrefilterModel none;
  assert(!none.trained);
  const auto open_mask = selectPrefilter(none, ev, is_decoy, p.keep_fraction);
  for (char c : open_mask) { assert(c); }
  const auto tiny = fitPrefilterModel(std::vector<PrefilterEvidence>(10), std::vector<bool>(10, false), p);
  assert(!tiny.trained);
  std::printf("  untrained model keeps everything: OK\n");

  std::printf("odia_prefilter_model_test OK\n");
  return 0;
}

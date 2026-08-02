// Ablation harness for the classifier and the anti-circularity mechanisms of
// odia_anchor_training.h.
//
// WHY OFFLINE. Extraction takes ~15 minutes and is bit-identical across classifier variants; the
// classifier takes seconds. Running the whole pipeline per arm would spend two hours measuring
// something that changes in the last thirty seconds, and would let node contention and extraction
// noise into a comparison that has nothing to do with either. One `-score_fixture` dump, N arms
// off the same matrix: the only thing that differs between arms is the thing under test.
//
// WHY EACH MECHANISM SEPARATELY. Bundled, an improvement and a regression cancel and the result
// reads as "no effect" -- which is the same number a mechanism that does nothing produces. The
// arms below turn exactly one thing on at a time against a stated baseline.
//
// WHAT THE NUMBERS MEAN, AND DO NOT. Identification count at 1% FDR is the headline, and it is
// self-reported: every arm computes its own q-values from its own decoy null, so an arm that
// breaks the null reports MORE identifications for being more wrong. That is not a hypothetical --
// the no-cross-validation arm exists in this table precisely to show it. Read the count together
// with the diagnostics printed beside it, and treat any arm that wins while its null degrades as
// a failure, not a result.
//
// usage: odia-ablate <fixture> [arm ...]

#include "odia_lda.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{

struct Fixture
{
  std::vector<std::vector<double>> X;
  std::vector<int> label;
  std::vector<long long> group;
  std::vector<std::string> names;
  std::vector<std::string> traml;    ///< precursor identity, for comparison against a reference list
};

bool load(const std::string& path, Fixture& f)
{
  std::ifstream in(path);
  if (!in) { return false; }
  std::size_t n = 0, m = 0;
  std::string head;
  if (!std::getline(in, head)) { return false; }
  {
    std::istringstream hs(head);
    hs >> n >> m;
    std::string nm;
    while (hs >> nm) { f.names.push_back(nm); }
  }
  f.X.reserve(n);
  f.label.reserve(n);
  f.group.reserve(n);

  // Two fixture generations exist and both are legitimate input: the original
  // testdata/lda_fixture.txt is `group label f0..fN`, and -score_fixture writes
  // `group label traml_id f0..fN`. Detect by COUNTING TOKENS on the first data row against the
  // declared width rather than by guessing from the header -- a guess that is wrong shifts every
  // feature by one column and still parses, which is the worst possible outcome.
  bool has_traml = false;
  {
    const std::streampos mark = in.tellg();
    std::string probe;
    while (std::getline(in, probe) && probe.empty()) {}
    std::istringstream ps(probe);
    std::size_t tok = 0;
    std::string t;
    while (ps >> t) { ++tok; }
    if (tok == m + 3) { has_traml = true; }
    else if (tok != m + 2)
    {
      std::fprintf(stderr, "odia-ablate: header declares %zu features but the first row has %zu "
                           "tokens (expected %zu or %zu)\n", m, tok, m + 2, m + 3);
      return false;
    }
    in.clear();
    in.seekg(mark);
  }

  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty()) { continue; }
    std::istringstream ss(line);
    long long g = 0;
    int l = 0;
    std::string tid = "-";
    ss >> g >> l;
    if (has_traml) { ss >> tid; }
    // A short row used to zero-fill in silence, which desyncs every column after it AND makes the
    // name-derived feature mask (mechanism 1) point at the wrong columns. Fail on the first one.
    std::vector<double> x(m);
    for (std::size_t j = 0; j < m; ++j)
    {
      std::string tok;
      if (!(ss >> tok))
      {
        std::fprintf(stderr, "odia-ablate: row %zu has %zu of %zu columns -- truncated fixture\n",
                     f.X.size() + 1, j, m);
        return false;
      }
      x[j] = (tok == "nan" || tok == "-nan") ? std::nan("") : std::strtod(tok.c_str(), nullptr);
    }
    f.group.push_back(g);
    f.label.push_back(l);
    f.traml.push_back(tid);
    f.X.push_back(std::move(x));
  }
  return !f.X.empty();
}

/// Identifications at a q-value cut: TARGET precursors, counted once each, best peak group only.
/// Counting rows would count a precursor as many times as it had candidate peaks.
std::size_t idsAt(const Fixture& f, const odia::ScoredGroups& s, double q)
{
  std::map<long long, std::size_t> best;
  for (std::size_t i = 0; i < f.X.size(); ++i)
  {
    auto it = best.find(f.group[i]);
    if (it == best.end() || s.dscore[i] > s.dscore[it->second]) { best[f.group[i]] = i; }
  }
  std::size_t n = 0;
  for (const auto& kv : best)
  {
    const std::size_t i = kv.second;
    if (f.label[i] == 1 && s.qvalue[i] < q) { ++n; }
  }
  return n;
}

/// Standardised separation between the class medians, in decoy standard deviations.
///
/// WHAT IT IS NOT: a test of null validity. Genuine discrimination and label leakage BOTH increase
/// it, so it cannot tell them apart, and an earlier version of this comment claimed it could.
/// Read it as a scale check on the score -- an arm whose gap collapses to ~0 has stopped
/// discriminating, and an arm whose gap explodes while its identification count also jumps is
/// worth suspecting -- but the only real test of the null is an EXTERNAL standard: an entrapment
/// set, or a second run. NaN is returned for a degenerate distribution rather than 0.0, because a
/// 0.0 there is indistinguishable from a healthy-but-uninformative arm.
double medianGap(const Fixture& f, const odia::ScoredGroups& s)
{
  std::vector<double> t, d;
  std::map<long long, std::size_t> best;
  for (std::size_t i = 0; i < f.X.size(); ++i)
  {
    auto it = best.find(f.group[i]);
    if (it == best.end() || s.dscore[i] > s.dscore[it->second]) { best[f.group[i]] = i; }
  }
  for (const auto& kv : best) { (f.label[kv.second] == 1 ? t : d).push_back(s.dscore[kv.second]); }
  if (t.size() < 3 || d.size() < 3) { return std::nan(""); }
  auto med = [](std::vector<double>& v) {
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
  };
  double mean = 0.0;
  for (double v : d) { mean += v; }
  mean /= static_cast<double>(d.size());
  double var = 0.0;
  for (double v : d) { var += (v - mean) * (v - mean); }
  const double sd = std::sqrt(var / static_cast<double>(d.size()));
  const double mt = med(t), md = med(d);
  return (sd > 0.0 && std::isfinite(sd)) ? (mt - md) / sd : std::nan("");
}

/// Mask hiding every RT-derived sub-score: mechanism 1's input.
std::vector<char> rtMask(const std::vector<std::string>& names, std::string& listed)
{
  std::vector<char> mask(names.size(), 1);
  for (std::size_t j = 0; j < names.size(); ++j)
  {
    std::string up = names[j];
    std::transform(up.begin(), up.end(), up.begin(), ::toupper);
    if (up.find("_RT") != std::string::npos || up.find("RT_") != std::string::npos)
    {
      mask[j] = 0;
      listed += (listed.empty() ? "" : ",") + names[j];
    }
  }
  return mask;
}

struct Arm
{
  std::string name;
  std::string what;
  odia::LDAParams p;
};

} // namespace

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::fprintf(stderr, "usage: odia-ablate <fixture> [arm ...]\n");
    return 2;
  }
  Fixture f;
  if (!load(argv[1], f))
  {
    std::fprintf(stderr, "odia-ablate: cannot read fixture %s\n", argv[1]);
    return 2;
  }
  std::set<long long> groups(f.group.begin(), f.group.end());
  std::size_t n_t = 0;
  for (int l : f.label) { n_t += (l == 1); }
  std::string rt_listed;
  const std::vector<char> mask = rtMask(f.names, rt_listed);

  // Group-level class balance, and the FINITE-DECOY FLOOR it implies.
  //
  // An earlier version of this said that more decoy groups than target ones makes 1% unreachable.
  // That is backwards, and adversarial review caught it. The estimator is
  //     q = (decoys_above + 1)/Ndec  /  (targets_above/Ntar)
  // which ratio-normalises both classes, so at perfect separation the q-value of the k-th target
  // is Ntar / (Ndec * k). Reaching q < 0.01 therefore needs k > 100 * Ntar/Ndec targets ranked
  // above every decoy -- and MORE decoys make that floor LOWER, not higher. What the balance
  // actually tells you is how many clean targets an arm must find before any of them can be
  // reported at all, which is worth printing; "unreachable" was simply wrong.
  std::map<long long, int> glab;
  for (std::size_t i = 0; i < f.X.size(); ++i) { glab[f.group[i]] = f.label[i]; }
  std::size_t gt = 0, gd = 0;
  for (const auto& kv : glab) { (kv.second == 1 ? gt : gd)++; }
  std::printf("fixture %s: %zu rows, %zu cols, %zu groups, %zu target rows\n", argv[1], f.X.size(),
              f.X.empty() ? 0 : f.X[0].size(), groups.size(), n_t);
  const long long floor_k =
    gd ? static_cast<long long>(std::ceil(100.0 * static_cast<double>(gt) / static_cast<double>(gd)))
       : -1;
  std::printf("group balance: %zu target, %zu decoy (%.2f decoy per target)\n", gt, gd,
              gt ? static_cast<double>(gd) / static_cast<double>(gt) : 0.0);
  if (floor_k > 0)
  {
    std::printf("finite-decoy floor: at perfect separation the %lldth-ranked target is the first "
                "to reach q<0.01\n", floor_k);
  }
  else
  {
    std::printf("finite-decoy floor: NO DECOYS -- q-values are meaningless here, every arm below "
                "is uninterpretable\n");
  }
  std::printf("RT-derived columns hidden by mechanism 1: %s\n",
              rt_listed.empty() ? "(none found -- mechanism 1 is a no-op on this fixture)"
                                : rt_listed.c_str());

  using C = odia::LDAParams::Classifier;
  std::vector<Arm> arms;
  auto add = [&](const std::string& n, const std::string& w, C c,
                 void (*tweak)(odia::LDAParams&) = nullptr) {
    Arm a;
    a.name = n;
    a.what = w;
    a.p.classifier = c;
    if (tweak) { tweak(a.p); }
    arms.push_back(std::move(a));
  };

  // ---- references ------------------------------------------------------------------------------
  add("lda", "linear discriminant (the original)", C::LDA);
  add("gbt", "boosted trees (the CURRENT default -- this is 'what we have')", C::GBT);

  // ---- the NN baseline: DIA-NN's arrangement, no mechanism on ----------------------------------
  // 12 nets differing ONLY in weight init, fixed iteration count, whole feature set at seed time.
  // Every mechanism below is measured against THIS, not against the GBT, or the comparison would
  // confound "does the mechanism help" with "is a network better than trees here".
  add("nn", "MLP ensemble, DIA-NN arrangement (no mechanism)", C::NN);

  // ---- one mechanism at a time -----------------------------------------------------------------
  add("nn_m1_seedmask", "M1: RT sub-scores hidden from the SEED fit", C::NN);
  arms.back().p.seed_mask = mask;

  add("nn_m2_strictseed", "M2 REMOVED: strict seed (0.01) instead of the lenient 0.15", C::NN,
      [](odia::LDAParams& p) { p.train_fdr_initial = 0.01; });

  add("nn_m3_bag70", "M3: each member trains on 70% of precursor groups", C::NN,
      [](odia::LDAParams& p) { p.bag_fraction = 0.70; });

  // NOT n_folds=1: that is clamped to 2 inside the scorer, so the arm ran ordinary 2-fold CV and
  // reported it as "no cross-validation". Found by adversarial review; the flag below cannot be
  // silently neutralised the same way.
  add("nn_m4_nocv", "M4 REMOVED: every group trains the model that scores it (FDR-INVALID)", C::NN,
      [](odia::LDAParams& p) { p.disable_cv = true; });

  // MATCHED PAIR. Both M5 arms raise n_iter from 3 to 8 so the rule has something to stop, which
  // means an M5-vs-default comparison confounds "the rule fired" with "five more iterations ran".
  // The control below is n_iter=8 with the rule OFF, so the difference between the two is the rule
  // alone. Without it the arm measures two changes and attributes them to one.
  add("nn_iter8_ctrl", "control for M5: n_iter=8, rule OFF (isolates the extra iterations)", C::NN,
      [](odia::LDAParams& p) { p.n_iter = 8; });

  add("nn_m5_compstop", "M5: stop when the positive SET stabilises or shrinks (vs nn_iter8_ctrl)",
      C::NN, [](odia::LDAParams& p) { p.stop_on_composition = true; p.n_iter = 8; });

  // ---- everything that is meant to be on -------------------------------------------------------
  add("nn_all", "M1+M3+M5 together (M2 and M4 are on by default)", C::NN,
      [](odia::LDAParams& p) {
        p.bag_fraction = 0.70;
        p.stop_on_composition = true;
        p.n_iter = 8;
      });
  arms.back().p.seed_mask = mask;

  // The stopping rule is learner-agnostic, so it gets a GBT arm too: if it helps the network and
  // not the trees, that is a fact about the network, and this is the only arm that can tell them
  // apart.
  add("gbt_m5_compstop", "M5 on the CURRENT default learner, to separate rule from learner", C::GBT,
      [](odia::LDAParams& p) { p.stop_on_composition = true; p.n_iter = 8; });

  if (argc > 2)
  {
    std::set<std::string> want(argv + 2, argv + argc);
    std::set<std::string> known;
    for (const Arm& a : arms) { known.insert(a.name); }
    for (const std::string& w : want)
    {
      if (!known.count(w))
      {
        // A typo used to erase every arm and exit 0 -- an empty table that reads as "nothing to
        // report" rather than "you asked for an arm that does not exist".
        std::fprintf(stderr, "odia-ablate: unknown arm '%s'. Known:", w.c_str());
        for (const std::string& k : known) { std::fprintf(stderr, " %s", k.c_str()); }
        std::fprintf(stderr, "\n");
        return 2;
      }
    }
    arms.erase(std::remove_if(arms.begin(), arms.end(),
                              [&](const Arm& a) { return !want.count(a.name); }),
               arms.end());
  }

  std::printf("\n%-18s %8s %8s %8s %7s %6s  %s\n", "arm", "IDs@1%", "IDs@5%", "vs gbt", "gap/sd",
              "iters", "what");
  std::printf("%s\n", std::string(112, '-').c_str());

  // Arms run one after another, each internally parallel. Running them concurrently instead would
  // make every timing a measure of contention, and the classifier's own threading is what the
  // production path uses.
  long long ref = -1;
  bool all_zero = true;
  for (Arm& a : arms)
  {
    const odia::ScoredGroups s = odia::scoreSemiSupervisedLDA(f.X, f.label, f.group, a.p);
    const std::size_t i1 = idsAt(f, s, 0.01), i5 = idsAt(f, s, 0.05);
    if (a.name == "gbt") { ref = static_cast<long long>(i1); }
    if (i1 > 0) { all_zero = false; }
    char delta[32] = "     ---";
    if (ref >= 0 && a.name != "gbt")
    {
      std::snprintf(delta, sizeof(delta), "%+8lld", static_cast<long long>(i1) - ref);
    }
    const double gap = medianGap(f, s);
    char gaps[16];
    if (std::isfinite(gap)) { std::snprintf(gaps, sizeof(gaps), "%7.2f", gap); }
    else                    { std::snprintf(gaps, sizeof(gaps), "%7s", "degen"); }
    std::printf("%-18s %8zu %8zu %8s %7s %3d/%-2d  %s\n", a.name.c_str(), i1, i5, delta, gaps,
                s.n_iterations_trained, s.n_iterations_skipped, a.what.c_str());
    std::fflush(stdout);
  }
  if (all_zero)
  {
    std::printf("\nEVERY arm returned 0. That is a property of the FIXTURE, not of the arms:\n"
                "  - %zu target vs %zu decoy groups%s\n"
                "  - no arm fitted a discriminant at all if the iters column reads 0/N\n"
                "Use a fixture from a full run (-score_fixture) before drawing any conclusion.\n",
                gt, gd, gd > gt ? " -- inverted, so no q-value can fall below 1%" : "");
  }
  std::printf("\niters = fitted/skipped, summed over folds. A skipped iteration means the arm found\n"
              "too few confident positives to fit -- its scores then come from an earlier model.\n");
  return 0;
}

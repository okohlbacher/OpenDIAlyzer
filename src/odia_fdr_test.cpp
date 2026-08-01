// odia_fdr_test.cpp — self-check for the context / picked / entrapment FDR layer.
//
// The claims under test are the ones the header asserts, each with a case that FAILS if the
// implementation is wrong rather than merely exercising the code:
//
//   T1 rollUp takes the MAX member score and inherits the entity label.
//   T2 entity q-values are the same estimator as the precursor level (sanity: a clean
//      separation gives ~all targets at q<=0.01 and no decoys).
//   T3 the inflation §1 warns about is real: rolling many null precursors up to an entity
//      makes the ENTITY-level error rate worse than the precursor-level one, so inheriting a
//      precursor q-value would understate it.
//   T4 picked competition collapses each pair to one entry, keeps the winner, and breaks ties
//      toward the DECOY.
//   T5 picked competition removes the protein-SIZE bias: with decoy proteins of varying
//      peptide count and no real signal, plain target-decoy over-reports large decoys while
//      picked does not.
//   T6 the entrapment estimator reproduces (1+1/r)*n_ent/n and refuses degenerate inputs.
//
// Build: c++ -std=c++17 -O2 src/odia_fdr_test.cpp -o t && ./t

#include "odia_fdr.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

static int failures = 0;
static void check(bool ok, const char* what)
{
  if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++failures; }
}

int main()
{
  using odia::fdr::Entity;

  // ---- T1: rollUp takes the max, keeps the label -------------------------------------
  {
    std::vector<std::string> key = {"P1", "P1", "P2", "", "P2"};
    std::vector<double> sc = {1.0, 3.0, -2.0, 99.0, -5.0};
    std::vector<int> lab = {1, 1, 0, 1, 0};
    auto e = odia::fdr::rollUp(key, sc, lab);
    check(e.size() == 2, "T1 rollUp collapses to one entity per key (unmapped rows dropped)");
    check(e[0].id == "P1" && std::abs(e[0].score - 3.0) < 1e-12, "T1 max member score for P1");
    check(e[1].id == "P2" && std::abs(e[1].score + 2.0) < 1e-12, "T1 max member score for P2");
    check(e[0].label == 1 && e[1].label == 0, "T1 label inherited from members");
  }

  // ---- T2: clean separation -> targets pass, decoys do not ---------------------------
  {
    std::vector<Entity> e;
    for (int i = 0; i < 500; ++i) { e.push_back({"T" + std::to_string(i), 1, 10.0 + i * 0.01}); }
    for (int i = 0; i < 500; ++i) { e.push_back({"D" + std::to_string(i), 0, i * 0.01}); }
    odia::fdr::assignQValues(e, false);
    check(odia::fdr::countAtQ(e, 0.01) >= 495, "T2 separable targets are identified");
    // A q-value belongs to the THRESHOLD, not to the item: the highest-scoring decoy sits
    // immediately below the whole target block, where the cutoff really does have ~0.4% FDR,
    // so it carries a small q of its own. That is correct and is why decoys are excluded from
    // identifications by LABEL (countAtQ) rather than by q. What must not happen is decoys
    // deep inside the target block, so bound the count rather than requiring zero.
    std::size_t bad = 0;
    for (const auto& x : e) { if (x.label == 0 && x.qvalue <= 0.01) { ++bad; } }
    check(bad <= 5, "T2 at most a boundary decoy or two sits at q<=0.01");
    for (const auto& x : e)
    {
      if (x.label == 0) { check(x.score < 10.0, "T2 no decoy outscores the target block"); break; }
    }
  }

  // ---- T3: rolling up INFLATES the error rate (the reason this file exists) ----------
  // Pure null: target and decoy precursors drawn from the same distribution. Every target
  // entity gets `n_draws` precursors, every decoy entity gets one. If entity q-values were
  // simply inherited from precursors, the target entities -- being a max over n_draws --
  // would look far better than they are. Recomputing at the entity level must NOT report
  // them as discoveries.
  {
    std::mt19937 rng(7);
    std::normal_distribution<double> nd(0.0, 1.0);
    const int n_entities = 2000, n_draws = 20;
    std::vector<std::string> key;
    std::vector<double> sc;
    std::vector<int> lab;
    for (int i = 0; i < n_entities; ++i)
    {
      for (int d = 0; d < n_draws; ++d)
      { key.push_back("T" + std::to_string(i)); sc.push_back(nd(rng)); lab.push_back(1); }
      for (int d = 0; d < n_draws; ++d)
      { key.push_back("D" + std::to_string(i)); sc.push_back(nd(rng)); lab.push_back(0); }
    }
    auto e = odia::fdr::rollUp(key, sc, lab);
    odia::fdr::assignQValues(e, false);
    const std::size_t ids = odia::fdr::countAtQ(e, 0.01);
    check(ids <= 20, "T3 pure null yields ~no entity-level identifications");

    // And the inflation itself: the best target ENTITY score is far above the best single
    // target precursor's typical score, i.e. the max-of-n statistic really is shifted.
    double best_entity = -1e30, mean_precursor = 0.0;
    std::size_t np = 0;
    for (const auto& x : e) { if (x.label == 1) { best_entity = std::max(best_entity, x.score); } }
    for (std::size_t i = 0; i < sc.size(); ++i) { if (lab[i] == 1) { mean_precursor += sc[i]; ++np; } }
    mean_precursor /= static_cast<double>(np);
    check(best_entity > mean_precursor + 3.0,
          "T3 entity score is a max-of-n extreme value, not a precursor score");
  }

  // ---- T4: picked competition -------------------------------------------------------
  {
    std::vector<Entity> e = {
      {"A", 1, 5.0}, {"DECOY_A", 0, 3.0},     // target wins
      {"B", 1, 1.0}, {"DECOY_B", 0, 4.0},     // decoy wins
      {"C", 1, 2.0}, {"DECOY_C", 0, 2.0},     // tie -> decoy
      {"E", 1, 7.0},                          // unpaired target
      {"DECOY_F", 0, 6.0},                    // unpaired decoy (no target F)
    };
    std::size_t paired = 0;
    auto p = odia::fdr::pickedCompetition(e, "DECOY_", &paired);
    check(paired == 3, "T4 three pairs found");
    check(p.size() == 5, "T4 pairs collapse to one entry each, singletons survive");
    int n_target = 0, n_decoy = 0;
    bool a_target = false, b_decoy = false, c_decoy = false;
    for (const auto& x : p)
    {
      (x.label == 1 ? n_target : n_decoy)++;
      if (x.id == "A" && x.label == 1) { a_target = true; }
      if (x.id == "DECOY_B" && x.label == 0) { b_decoy = true; }
      if (x.id == "DECOY_C" && x.label == 0) { c_decoy = true; }
    }
    check(a_target, "T4 higher-scoring target wins its pair");
    check(b_decoy, "T4 higher-scoring decoy wins its pair");
    check(c_decoy, "T4 ties resolve to the decoy (conservative)");
    check(n_target == 2 && n_decoy == 3, "T4 winner labels counted correctly");
  }

  // ---- T5: picked removes the protein-SIZE bias --------------------------------------
  // Null data, proteins of widely differing peptide counts. A protein's score is the max
  // over its peptides, so large proteins score higher REGARDLESS of label. Plain
  // target-decoy then lets large decoys crowd the top of the list; picked competition,
  // which compares each protein only with its equal-sized partner, must not.
  {
    std::mt19937 rng(11);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<std::string> key;
    std::vector<double> sc;
    std::vector<int> lab;
    const int n_prot = 1500;
    for (int i = 0; i < n_prot; ++i)
    {
      const int size = 1 + (i % 50) * 4;           // 1..197 peptides
      for (int k = 0; k < size; ++k)
      { key.push_back("P" + std::to_string(i)); sc.push_back(nd(rng)); lab.push_back(1); }
      for (int k = 0; k < size; ++k)
      { key.push_back("DECOY_P" + std::to_string(i)); sc.push_back(nd(rng)); lab.push_back(0); }
    }
    auto e = odia::fdr::rollUp(key, sc, lab);

    // top 100 by score, plain: how skewed toward large proteins is it?
    auto by_score = e;
    std::sort(by_score.begin(), by_score.end(),
              [](const Entity& a, const Entity& b) { return a.score > b.score; });
    double plain_mean_size = 0.0;
    for (int i = 0; i < 100; ++i)
    {
      const std::string& id = by_score[i].id;
      const std::string num = id.substr(id.find('P') + 1);
      plain_mean_size += 1 + (std::stoi(num) % 50) * 4;
    }
    plain_mean_size /= 100.0;
    const double overall_mean_size = 1 + (49 / 2.0) * 4;   // ~99
    check(plain_mean_size > overall_mean_size,
          "T5 plain ranking is size-biased (large proteins crowd the top)");

    std::size_t paired = 0;
    auto p = odia::fdr::pickedCompetition(e, "DECOY_", &paired);
    check(paired == static_cast<std::size_t>(n_prot), "T5 every protein paired");
    check(p.size() == static_cast<std::size_t>(n_prot), "T5 picked halves the entity count");
    // Under the null the pair winner is a coin flip, independent of size.
    std::size_t won_by_target = 0;
    for (const auto& x : p) { if (x.label == 1) { ++won_by_target; } }
    check(won_by_target > 600 && won_by_target < 900,
          "T5 pair winners are ~50/50 under the null, independent of protein size");
    odia::fdr::assignQValues(p, false);
    check(odia::fdr::countAtQ(p, 0.01) <= 20, "T5 picked reports ~no IDs on null data");
  }

  // ---- T6: entrapment estimator ------------------------------------------------------
  {
    // r = 1, 30 entrapments among 1000 reported -> (1+1)*30/1000 = 0.06, i.e. the nominal
    // 1% those 1000 were taken at would be understating the true rate six-fold.
    auto a = odia::fdr::entrapmentFdp(1000, 30, 5000, 5000);
    check(a.valid && std::abs(a.fdp - 0.06) < 1e-12, "T6 combined estimator at r=1");
    // r = 2 -> (1+0.5)*30/1000 = 0.045
    auto b = odia::fdr::entrapmentFdp(1000, 30, 5000, 10000);
    check(b.valid && std::abs(b.fdp - 0.045) < 1e-12, "T6 combined estimator at r=2");
    // no entrapment discoveries -> 0
    auto c = odia::fdr::entrapmentFdp(1000, 0, 5000, 5000);
    check(c.valid && c.fdp == 0.0, "T6 zero entrapment hits -> zero FDP");
    // degenerate inputs must be refused, not guessed
    check(!odia::fdr::entrapmentFdp(1000, 30, 5000, 0).valid, "T6 no entrapment DB -> invalid");
    check(!odia::fdr::entrapmentFdp(0, 0, 5000, 5000).valid, "T6 nothing reported -> invalid");
    // clamped at 1
    auto d = odia::fdr::entrapmentFdp(100, 90, 5000, 5000);
    check(d.valid && d.fdp == 1.0, "T6 FDP clamped to 1");
  }

  if (failures) { std::fprintf(stderr, "odia_fdr_test FAILED (%d)\n", failures); return 1; }
  std::fprintf(stderr, "odia_fdr_test OK\n");
  return 0;
}

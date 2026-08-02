// Self-check for odia_anchor_training.h.
//
// The failure this file exists to prevent is silent and looks like success: a model that overfits
// its own seed improves every internal metric while identifying a narrower set. So the tests are
// about the DETECTORS, not about accuracy.

#include "odia_anchor_training.h"

#include <cassert>
#include <cstdio>
#include <random>

using namespace odia;

int main()
{
  std::printf("odia_anchor_training_test\n");

  // ---- 1. bagging is deterministic and keyed on the GROUP, not the row -------------------------
  // All rows of one precursor must land in the same bag, or a precursor is half-trained-on and
  // half-held-out, which is the leak group-wise CV exists to prevent.
  std::size_t in = 0;
  for (std::int64_t g = 0; g < 10000; ++g) { in += inBag(3, g, 42, 0.7); }
  const double frac = double(in) / 10000.0;
  std::printf("  bag fraction for member 3: %.3f (asked 0.700)\n", frac);
  assert(std::fabs(frac - 0.7) < 0.02);
  assert(inBag(3, 12345, 42, 0.7) == inBag(3, 12345, 42, 0.7));      // reproducible
  std::size_t differ = 0;
  for (std::int64_t g = 0; g < 1000; ++g) { differ += (inBag(0, g, 42, 0.7) != inBag(1, g, 42, 0.7)); }
  std::printf("  members 0 and 1 disagree on %zu/1000 groups\n", differ);
  assert(differ > 200);   // members must actually see different data, else bagging is decorative

  // ---- 2. Jaccard detects stabilisation --------------------------------------------------------
  std::vector<std::size_t> a = {1,2,3,4,5,6,7,8,9,10}, b = a, c = {1,2,3,4,5};
  assert(jaccardOverlap(a, b) == 1.0);
  assert(jaccardOverlap(a, c) == 0.5);
  std::printf("  jaccard identical=1.0, half-subset=0.5: OK\n");

  // ---- 3. THE COLLAPSE DETECTOR ----------------------------------------------------------------
  // A shrinking positive set is the signature of a model converging onto a subset of its seed.
  // Identification counts RISE during this, so a rule watching the score cannot see it.
  AnchorTrainingParams p;
  AnchorTrainingReport rep;
  rep.iterations_run = 1;
  std::vector<std::size_t> prev(1000), curr(700);
  std::iota(prev.begin(), prev.end(), 0);
  std::iota(curr.begin(), curr.end(), 0);
  const bool cont = anchorIterationShouldContinue(prev, curr, p, rep);
  std::printf("  shrink 1000 -> 700: continue=%d collapsed=%d\n", (int)cont, (int)rep.collapsed);
  assert(!cont && rep.collapsed);
  assert(!rep.note.empty());

  // growth that stabilises must stop WITHOUT flagging collapse
  AnchorTrainingReport rep2;
  rep2.iterations_run = 1;
  std::vector<std::size_t> grew(1002);
  std::iota(grew.begin(), grew.end(), 0);
  const bool cont2 = anchorIterationShouldContinue(prev, grew, p, rep2);
  std::printf("  grow 1000 -> 1002 (jaccard %.3f): continue=%d converged=%d collapsed=%d\n",
              rep2.jaccard.back(), (int)cont2, (int)rep2.converged, (int)rep2.collapsed);
  assert(!rep2.collapsed);
  assert(rep2.converged && !cont2);

  // ---- 4. too few anchors must REFUSE to train, not train badly --------------------------------
  std::vector<std::vector<double>> X(50, std::vector<double>(4, 0.0));
  std::vector<std::size_t> few = {0,1,2}, negs = {10,11,12};
  std::vector<std::int64_t> grp(50, 0);
  for (std::size_t i = 0; i < grp.size(); ++i) { grp[i] = static_cast<std::int64_t>(i); }
  const auto none = trainBaggedOnAnchors(X, few, negs, grp, p);
  assert(none.empty());
  assert(baggedScore(none, X[0]) == 0.0);
  std::printf("  %zu anchors < min %zu: refuses to train, scores 0: OK\n", few.size(), p.min_anchors);

  // ---- 5. disagreement is measurable and non-zero on a real fit ---------------------------------
  std::mt19937 rng(3);
  std::normal_distribution<double> nz(0.0, 1.0);
  std::vector<std::vector<double>> XX;
  std::vector<std::size_t> pos, neg;
  std::vector<std::int64_t> group;
  for (int i = 0; i < 3000; ++i)
  {
    const bool t = (i % 2 == 0);
    XX.push_back({(t ? 1.2 : -1.2) + nz(rng), nz(rng), nz(rng), nz(rng)});
    group.push_back(i / 3);                       // 3 rows per precursor
    (t ? pos : neg).push_back(XX.size() - 1);
  }
  AnchorTrainingParams pp;
  pp.nn.hidden = {8, 4};
  pp.nn.n_nets = 6;
  pp.nn.epochs = 60;
  pp.nn.lr = 0.4;
  pp.min_anchors = 100;
  const auto members = trainBaggedOnAnchors(XX, pos, neg, group, pp);
  std::printf("  trained %zu bagged members\n", members.size());
  assert(members.size() >= 4);
  std::vector<std::size_t> sample(pos.begin(), pos.begin() + 200);
  const double dis = anchorDisagreement(members, XX, sample);
  std::printf("  ensemble disagreement on anchors: %.4f\n", dis);
  assert(dis > 0.0);        // zero would mean bagging changed nothing

  std::size_t right = 0;
  for (std::size_t i = 0; i < XX.size(); ++i)
  {
    right += ((baggedScore(members, XX[i]) > 0.0) == (i % 2 == 0));
  }
  std::printf("  bagged accuracy: %.3f\n", double(right) / XX.size());
  assert(double(right) / XX.size() > 0.75);

  std::printf("odia_anchor_training_test OK\n");
  return 0;
}

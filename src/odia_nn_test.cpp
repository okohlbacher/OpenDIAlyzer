// Self-check for odia_nn.h.
//
// The properties whose violation is SILENT: a network that is not bit-reproducible reintroduces the
// exact class of bug that produced a phantom +/-1% identification spread in this project, and an
// IM-feature layout that imputes a missing score as zero makes "no mobility measured" look
// identical to "a perfect mobility match".

#include "odia_nn.h"

#include <cassert>
#include <cstdio>
#include <random>

using namespace odia;

int main()
{
  std::printf("odia_nn_test\n");

  // ---- 1. IM layout: absent means NARROWER, present means paired with an indicator -------------
  const std::vector<std::string> ms2 = {"a", "b", "c"}, ms1 = {"m1", "m2"};
  const auto im = imFeatureNames();
  FeatureBlocks nb; nb.ms2 = true; nb.ms1 = false; nb.im = false;
  FeatureBlocks wb; wb.ms2 = true; wb.ms1 = true;  wb.im = true;
  const auto lay_no = nnFeatureLayout(ms2, ms1, im, nb);
  const auto lay_yes = nnFeatureLayout(ms2, ms1, im, wb);
  assert(lay_no.size() == 3);
  assert(lay_yes.size() == 3 + 2 + im.size() * 2);          // each IM score + its presence flag
  std::size_t flags = 0;
  for (const auto& n : lay_yes) { flags += (n.size() > 9 && n.rfind("__present") == n.size() - 9); }
  assert(flags == im.size());
  std::printf("  layout: no-IM %zu features, with-IM %zu (%zu presence flags): OK\n",
              lay_no.size(), lay_yes.size(), flags);

  // ---- 2. init is order-independent ------------------------------------------------------------
  // Same coordinates -> same weight, regardless of the order they are requested in.
  const double w1 = nnInitWeight(3, 1, 7, 2, 42, 0.5);
  const double w0 = nnInitWeight(0, 0, 0, 0, 42, 0.5);
  assert(nnInitWeight(3, 1, 7, 2, 42, 0.5) == w1);
  assert(w0 != w1);
  std::printf("  counter-based init reproducible and coordinate-dependent: OK\n");

  // ---- 3. it must LEARN a signal a linear model cannot: an XOR-like interaction ----------------
  // Two features where the label depends on their INTERACTION, not either alone. This is the whole
  // justification for a network over the LDA path.
  std::mt19937 rng(1);
  std::normal_distribution<double> nz(0.0, 0.35);
  std::vector<std::vector<double>> X;
  std::vector<std::size_t> pos, neg;
  for (int i = 0; i < 4000; ++i)
  {
    const double u = (i % 2) ? 1.0 : -1.0;
    const double v = ((i / 2) % 2) ? 1.0 : -1.0;
    X.push_back({u + nz(rng), v + nz(rng), nz(rng)});
    ((u * v > 0) ? pos : neg).push_back(X.size() - 1);
  }
  NNParams p;
  p.hidden = {16, 8};
  p.n_nets = 6;
  p.epochs = 400;                 // the real pipeline uses 1 per semi-supervised iteration; here we
  p.lr = 0.5;                     // train to convergence to prove the capacity exists at all
  NNEnsemble net;
  assert(net.fit(X, pos, neg, p));
  std::size_t right = 0;
  for (std::size_t i = 0; i < X.size(); ++i)
  {
    const bool is_pos = (X[i][0] * X[i][1] > 0);
    right += ((net.score(X[i]) > 0.0) == is_pos);
  }
  const double acc = double(right) / X.size();
  std::printf("  XOR-interaction accuracy: %.3f (a linear discriminant scores ~0.5)\n", acc);
  assert(acc > 0.85);

  // ---- 4. bit-reproducibility ------------------------------------------------------------------
  NNEnsemble net2;
  assert(net2.fit(X, pos, neg, p));
  double worst = 0.0;
  for (std::size_t i = 0; i < X.size(); ++i)
  {
    worst = std::max(worst, std::fabs(net.score(X[i]) - net2.score(X[i])));
  }
  std::printf("  refit max |delta|: %.3e\n", worst);
  assert(worst == 0.0);

  // ---- 5. an untrained ensemble must be inert, not random --------------------------------------
  NNEnsemble none;
  assert(!none.trained());
  assert(none.score(X[0]) == 0.0);
  std::printf("  untrained ensemble scores 0: OK\n");

  // ---- thread-count invariance -----------------------------------------------------------------
  // The reason trainOne partitions into a FIXED 512 chunks rather than one per thread. Floating-
  // point addition is not associative, so a gradient summed over a thread-dependent partition is a
  // thread-dependent model. This is the test that would have caught the same bug in the GBT
  // histogram reduction, and it must be BIT-exact, not approximate: a tolerance here would pass
  // while the model quietly depended on the machine it ran on.
  {
    NNParams tp = p;
    tp.n_threads = 1;
    NNEnsemble one;
    assert(one.fit(X, pos, neg, tp));
    double worst = 0.0;
    for (int t : {2, 7, 16, 64})
    {
      NNParams mp = p;
      mp.n_threads = t;
      NNEnsemble many;
      assert(many.fit(X, pos, neg, mp));
      for (const auto& row : X)
      {
        worst = std::max(worst, std::fabs(one.score(row) - many.score(row)));
      }
    }
    std::printf("  thread-count invariance (1 vs 2/7/16/64): max |delta| %.3e\n", worst);
    assert(worst == 0.0);
  }

  // ---- THE DEFAULTS MUST TRAIN --------------------------------------------------------------
  // The test above uses epochs=400, lr=0.5 to show the CAPACITY exists. That is not the same claim
  // as "the configuration we ship learns anything", and the gap between the two is precisely how a
  // full-batch implementation reached the benchmark: at epochs=1, full-batch is ONE gradient step,
  // the output stayed near-constant, and the run reported 0 identifications with all 9
  // fold-iterations skipped for want of confident positives -- while this file printed OK.
  //
  // So: default NNParams, a linearly separable problem, and an accuracy floor. Anything that makes
  // the shipped defaults stop learning now fails here instead of on the cluster 40 minutes in.
  {
    std::mt19937 rng(99);
    std::normal_distribution<double> g(0.0, 1.0);
    std::vector<std::vector<double>> D;
    std::vector<std::size_t> dpos, dneg;
    const std::size_t kN = 20000;
    for (std::size_t i = 0; i < kN; ++i)
    {
      const bool positive = (i % 2 == 0);
      std::vector<double> row(8);
      for (double& v : row) { v = g(rng); }
      if (positive) { row[0] += 1.5; row[3] += 1.0; }
      D.push_back(row);
      (positive ? dpos : dneg).push_back(i);
    }
    NNParams dp;                                   // DEFAULTS. Deliberately not tuned here.
    NNEnsemble e;
    assert(e.fit(D, dpos, dneg, dp));
    std::size_t correct = 0;
    for (std::size_t i = 0; i < kN; ++i)
    {
      const bool predicted_positive = e.score(D[i]) > 0.0;
      if (predicted_positive == (i % 2 == 0)) { ++correct; }
    }
    const double acc = static_cast<double>(correct) / static_cast<double>(kN);
    std::printf("  DEFAULT params (epochs=%d, batch=%d, lr=%.3g) accuracy: %.3f\n", dp.epochs,
                dp.batch_size, dp.lr, acc);
    // 0.70 is far below what this problem allows; it is a floor against "did not train at all",
    // not a quality bar. One full-batch step scores ~0.5 here.
    assert(acc > 0.70);
  }

  std::printf("odia_nn_test OK\n");
  return 0;
}

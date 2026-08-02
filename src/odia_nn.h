#pragma once
// A deterministic feed-forward network ensemble for peak-group scoring, shaped after DIA-NN's.
//
// WHY A NETWORK AT ALL. The current classifier consumes 23 MS2-scope sub-scores (35 with
// -ms1_scores). DIA-NN trains an ensemble of feed-forward networks on 73. Measured on this
// benchmark, target and decoy discriminant distributions nearly coincide -- medians -0.205 vs
// -0.295, with 46.8% of targets inside the decoy interquartile range -- so the binding constraint
// is discrimination, not extraction: 1,261 of 1,458 recoverable misses are extracted and scored but
// not separated.
//
// A linear discriminant cannot represent an interaction, and the interactions are the point. The
// clearest measured example: RT deviation means different things at different points in the
// gradient -- median |delta_rt| runs 6.8 s in the first RT octile to 95.9 s in the last, a 14x
// swing. A linear model must pick one slope; a network can learn "90 s is unremarkable late and
// damning early". Trees can too, which is why the existing GBT is the fallback and the honest
// baseline to beat.
//
// WHAT IS COPIED FROM DIA-NN, AND WHAT IS NOT
//   copied : ensemble of small nets, tanh hidden layers, ONE epoch per semi-supervised iteration,
//            retrained per run rather than shipped pretrained.
//   not    : 73 hand-built sub-scores we do not have; a softmax head (for two classes it is a
//            reparameterisation of a single logit, so we use the logit and keep the score on the
//            same scale as the LDA/GBT path).
//
// DETERMINISM IS A HARD REQUIREMENT, NOT A PREFERENCE. Two order-dependency bugs in this project
// produced a +/-1% identification spread that was mistaken for a noise floor and used to size
// experiments for weeks. Both were invisible to the obvious tests: one was identical at 1/8/64
// threads, and neither appeared without OpenMP compiled in. Everything below is written so that
// the output is bit-identical regardless of thread count:
//   - weights initialise from a counter-based hash of (net index, layer, row, col), never an RNG
//     whose draw order depends on scheduling;
//   - gradients accumulate into a FIXED number of chunk buffers determined by row count alone,
//     never by thread count, and merge in chunk-index order;
//   - no dropout, no batch shuffling, no atomics on floats.

// TRAINING DATA. There is none to collect and none to ship -- the model is fitted per run from the
// run's own decoys, which is what makes it robust to instrument, gradient and library changing
// underneath it.
//
//   negatives : every decoy peak group. Known-absent by construction, ~103k rank-1 rows on this
//               benchmark against ~104k targets, so the classes are naturally balanced at rank 1.
//   positives : NOT "all targets". The confident subset -- best row per target group clearing a
//               LENIENT training-FDR cutoff (0.15 on the first semi-supervised iteration, 0.05
//               after). Most library targets are genuinely absent, so treating them all as
//               positives trains on a label that is ~99% wrong.
//   held out  : nothing, at fit time. Generalisation is enforced by the surrounding k-fold split
//               BY PRECURSOR GROUP, so a row is never scored by a model trained on its own
//               precursor.
//
// What must NEVER become training data: identification outcomes from a previous pass of the same
// run, or anything derived from the q-values this model's output will be used to compute. That is
// selection on the test statistic and it inflates FDR without any visible symptom.
//
// An entrapment library validates this; it does not train it.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odia
{

// ---------------------------------------------------------------------------------------------
// Feature assembly, including ion mobility when the acquisition has it
// ---------------------------------------------------------------------------------------------

/// Which optional feature blocks a run actually provides. IM sub-scores exist only for
/// mobility-separated acquisitions (diaPASEF and similar); on an Orbitrap/Astral run the columns
/// are absent entirely.
struct FeatureBlocks
{
  bool ms2 = true;      ///< always present
  bool ms1 = false;     ///< -ms1_scores
  bool im = false;      ///< ion mobility present in the acquisition
};

/// Assemble the input layout for a run.
///
/// IM HANDLING. A missing IM score must NOT be imputed as zero and fed in silently: zero is a
/// legitimate value for several of these (a delta of 0 is a perfect match), so an absent column
/// would masquerade as a perfect one. Two rules:
///   1. when the run has no IM, the IM columns are not in the layout at all -- the network is
///      simply narrower, rather than being fed constants that it must learn to ignore;
///   2. when the run HAS IM but an individual row lacks a value (no mobilogram extracted for that
///      peak group), the value is imputed to the column mean AND a companion presence indicator is
///      set to 0. The network can then represent "IM says nothing here" separately from "IM says
///      the delta is zero".
inline std::vector<std::string> nnFeatureLayout(const std::vector<std::string>& ms2_names,
                                                const std::vector<std::string>& ms1_names,
                                                const std::vector<std::string>& im_names,
                                                const FeatureBlocks& blocks)
{
  std::vector<std::string> out;
  if (blocks.ms2) { out.insert(out.end(), ms2_names.begin(), ms2_names.end()); }
  if (blocks.ms1) { out.insert(out.end(), ms1_names.begin(), ms1_names.end()); }
  if (blocks.im)
  {
    out.insert(out.end(), im_names.begin(), im_names.end());
    for (const auto& n : im_names) { out.push_back(n + "__present"); }
  }
  return out;
}

/// The IM sub-scores this project's extraction emits, when the acquisition is mobility-separated.
/// Kept as a named list so a run can be checked against it and a missing column reported rather
/// than silently dropped.
inline std::vector<std::string> imFeatureNames()
{
  return {"var_ms2_im_xcorr_shape", "var_ms2_im_xcorr_coelution", "var_ms2_im_delta_score",
          "var_ms2_im_log_intensity", "var_ms1_im_ms1_delta_score",
          "var_im_xcorr_coelution_contrast", "var_im_xcorr_shape_contrast",
          "var_im_xcorr_coelution_combined", "var_im_xcorr_shape_combined", "delta_im"};
}

// ---------------------------------------------------------------------------------------------
// The network
// ---------------------------------------------------------------------------------------------

struct NNParams
{
  std::vector<int> hidden = {32, 16, 16, 8, 8};  ///< DIA-NN's shape is 5 tanh layers; widths here
                                                 ///< are scaled to a 35-60 input rather than 73
  int n_nets = 12;            ///< ensemble size, as DIA-NN
  int epochs = 1;             ///< ONE epoch per semi-supervised iteration, as DIA-NN. The ensemble
                              ///< and the outer iteration do the work, not long training.
  double lr = 0.05;
  double l2 = 1e-5;
  std::uint64_t seed = 42;
  int n_threads = 0;
  /// Optional input mask, one entry per feature; empty means "use everything". A 0 entry zeroes
  /// that input in the FORWARD pass, so the feature can neither influence the output nor receive
  /// gradient (dW = delta * x = 0). This is mechanism 1 of odia_anchor_training.h: the seed fit
  /// must not see the features the anchors will be used to calibrate.
  ///
  /// A mask rather than a narrower matrix because the mask travels WITH the fitted model -- a
  /// model trained without a feature must also be SCORED without it, and an ensemble whose
  /// members were fitted under different masks still shares one X.
  std::vector<char> mask;
};

/// Counter-based weight init. A hash of the coordinates, NOT a sequential RNG: a draw order that
/// depends on thread scheduling is exactly the class of bug that cost this project weeks.
inline double nnInitWeight(std::uint64_t net, std::uint64_t layer, std::uint64_t i,
                           std::uint64_t j, std::uint64_t seed, double scale)
{
  std::uint64_t h = seed;
  for (std::uint64_t v : {net, layer, i, j})
  {
    h ^= v + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
    h *= 0xBF58476D1CE4E5B9ull;
    h ^= h >> 31;
  }
  // uniform in [-1,1), then scaled -- deterministic for a given coordinate, order-independent
  const double u = static_cast<double>(h >> 11) / static_cast<double>(1ull << 53);
  return (2.0 * u - 1.0) * scale;
}

/// One multilayer perceptron with tanh hidden activations and a single linear output (the logit).
class MLP
{
public:
  void init(int n_in, const std::vector<int>& hidden, std::uint64_t net, std::uint64_t seed,
            const std::vector<char>& mask = {})
  {
    mask_ = (static_cast<int>(mask.size()) == n_in) ? mask : std::vector<char>();
    dims_.clear();
    dims_.push_back(n_in);
    for (int h : hidden) { dims_.push_back(h); }
    dims_.push_back(1);
    W_.assign(dims_.size() - 1, {});
    b_.assign(dims_.size() - 1, {});
    for (std::size_t l = 0; l + 1 < dims_.size(); ++l)
    {
      const int fan_in = dims_[l], fan_out = dims_[l + 1];
      // Xavier: keeps tanh in its responsive range instead of saturating on the first pass
      const double scale = std::sqrt(6.0 / static_cast<double>(fan_in + fan_out));
      W_[l].resize(static_cast<std::size_t>(fan_in) * fan_out);
      b_[l].assign(fan_out, 0.0);
      for (int i = 0; i < fan_in; ++i)
      {
        for (int j = 0; j < fan_out; ++j)
        {
          W_[l][static_cast<std::size_t>(i) * fan_out + j] =
            nnInitWeight(net, l, static_cast<std::uint64_t>(i), static_cast<std::uint64_t>(j),
                         seed, scale);
        }
      }
    }
  }

  /// Forward pass. `act` is filled with each layer's activations so backward can reuse them.
  double forward(const double* x, std::vector<std::vector<double>>& act) const
  {
    act.resize(dims_.size());
    act[0].assign(x, x + dims_[0]);
    if (!mask_.empty())
    {
      for (int i = 0; i < dims_[0]; ++i) { if (!mask_[static_cast<std::size_t>(i)]) { act[0][i] = 0.0; } }
    }
    for (std::size_t l = 0; l + 1 < dims_.size(); ++l)
    {
      const int fi = dims_[l], fo = dims_[l + 1];
      act[l + 1].assign(static_cast<std::size_t>(fo), 0.0);
      for (int j = 0; j < fo; ++j)
      {
        double s = b_[l][j];
        for (int i = 0; i < fi; ++i) { s += act[l][i] * W_[l][static_cast<std::size_t>(i) * fo + j]; }
        // hidden layers use tanh; the output layer is linear (the logit)
        act[l + 1][j] = (l + 2 < dims_.size()) ? std::tanh(s) : s;
      }
    }
    return act.back()[0];
  }

  double score(const std::vector<double>& x) const
  {
    std::vector<std::vector<double>> act;
    return forward(x.data(), act);
  }

  const std::vector<int>& dims() const { return dims_; }
  std::vector<std::vector<double>>& weights() { return W_; }
  std::vector<std::vector<double>>& biases() { return b_; }
  const std::vector<std::vector<double>>& weights() const { return W_; }
  const std::vector<std::vector<double>>& biases() const { return b_; }

private:
  std::vector<int> dims_;
  std::vector<std::vector<double>> W_, b_;
  std::vector<char> mask_;
};

/// Ensemble of MLPs, averaged. Each member differs only in its init hash, so the ensemble is
/// reproducible and its members are genuinely different.
class NNEnsemble
{
public:
  /// Train on standardised features. `pos` / `neg` are row indices; the caller owns the
  /// semi-supervised selection of `pos` exactly as for the LDA and GBT paths.
  bool fit(const std::vector<std::vector<double>>& X,
           const std::vector<std::size_t>& pos,
           const std::vector<std::size_t>& neg,
           const NNParams& p)
  {
    if (X.empty() || pos.empty() || neg.empty()) { return false; }
    n_in_ = static_cast<int>(X[0].size());
    params_ = p;
    nets_.assign(static_cast<std::size_t>(p.n_nets), MLP());

    // Class weighting: the negative set is typically far larger, and an unweighted fit collapses
    // to predicting "decoy" for everything.
    const double w_pos = 1.0;
    const double w_neg = static_cast<double>(pos.size()) / static_cast<double>(neg.size());

    std::vector<std::size_t> rows;
    rows.reserve(pos.size() + neg.size());
    rows.insert(rows.end(), pos.begin(), pos.end());
    rows.insert(rows.end(), neg.begin(), neg.end());
    std::vector<double> y(rows.size(), 1.0), w(rows.size(), w_pos);
    for (std::size_t k = pos.size(); k < rows.size(); ++k) { y[k] = 0.0; w[k] = w_neg; }
    // Fixed order, independent of how the caller assembled pos/neg.
    std::sort(rows.begin(), rows.end());

    for (int n = 0; n < p.n_nets; ++n)
    {
      nets_[static_cast<std::size_t>(n)].init(n_in_, p.hidden, static_cast<std::uint64_t>(n),
                                              p.seed, p.mask);
      trainOne(nets_[static_cast<std::size_t>(n)], X, pos, neg, w_neg, p);
    }
    trained_ = true;
    return true;
  }

  /// Mean logit over the ensemble.
  double score(const std::vector<double>& x) const
  {
    if (!trained_) { return 0.0; }
    double s = 0.0;
    for (const auto& net : nets_) { s += net.score(x); }
    return s / static_cast<double>(nets_.size());
  }

  bool trained() const { return trained_; }
  int inputWidth() const { return n_in_; }

private:
  /// One net, `epochs` passes of full-batch gradient descent on the logistic loss.
  ///
  /// Full-batch, not minibatch: a minibatch order is one more thing that has to be made
  /// order-independent, and at one epoch the ensemble supplies the variance that minibatching
  /// would.
  static void trainOne(MLP& net, const std::vector<std::vector<double>>& X,
                       const std::vector<std::size_t>& pos, const std::vector<std::size_t>& neg,
                       double w_neg, const NNParams& p)
  {
    auto& W = net.weights();
    auto& B = net.biases();
    const auto& dims = net.dims();

    // Flat parameter layout, so one gradient buffer per chunk is ONE allocation rather than a
    // vector-of-vectors per layer. At 512 chunks the difference is 3072 allocations an epoch.
    std::vector<std::size_t> wo(W.size()), bo(B.size());
    std::size_t np = 0;
    for (std::size_t l = 0; l < W.size(); ++l) { wo[l] = np; np += W[l].size(); }
    for (std::size_t l = 0; l < B.size(); ++l) { bo[l] = np; np += B[l].size(); }

    // One flat row list: positives first, then negatives, each in the caller's sorted order.
    std::vector<std::size_t> rows;
    rows.reserve(pos.size() + neg.size());
    rows.insert(rows.end(), pos.begin(), pos.end());
    rows.insert(rows.end(), neg.begin(), neg.end());
    const std::size_t n_pos = pos.size(), n_rows = rows.size();
    if (n_rows == 0) { return; }

    // A FIXED chunk count, deliberately independent of the thread count.
    //
    // Floating-point addition is not associative, so a partial sum depends on WHICH rows went into
    // it. Partitioning by thread count would therefore make the trained weights a function of how
    // many threads happened to be available -- the exact class of bug that cost this project weeks
    // in the GBT histogram reduction, where `lo = n_rows*c/nchunk` had the same shape. Fixing the
    // partition at 512 and summing the chunks in INDEX order makes the result bit-identical on 1
    // thread and on 224, which is the property the determinism tests check.
    //
    // 512 rather than 64: with 224 threads a chunk count near the thread count leaves the tail
    // badly balanced, and the buffer is only 512 * ~1.8k doubles = 7 MB.
    constexpr std::size_t kChunks = 512;
    const std::size_t n_chunks = std::min<std::size_t>(kChunks, n_rows);

    for (int ep = 0; ep < p.epochs; ++ep)
    {
      std::vector<double> part(n_chunks * np, 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(p.n_threads > 0 ? p.n_threads : omp_get_max_threads())
#endif
      for (long long c = 0; c < static_cast<long long>(n_chunks); ++c)
      {
        const std::size_t lo = n_rows * static_cast<std::size_t>(c) / n_chunks;
        const std::size_t hi = n_rows * static_cast<std::size_t>(c + 1) / n_chunks;
        double* g = part.data() + static_cast<std::size_t>(c) * np;
        std::vector<std::vector<double>> act, delta;
        for (std::size_t k = lo; k < hi; ++k)
        {
          const std::size_t r = rows[k];
          const double label = (k < n_pos) ? 1.0 : 0.0;
          const double weight = (k < n_pos) ? 1.0 : w_neg;
          const double logit = net.forward(X[r].data(), act);
          const double pr = 1.0 / (1.0 + std::exp(-logit));
          delta.assign(dims.size(), {});
          delta.back().assign(1, weight * (pr - label));       // dL/dlogit
          for (std::size_t l = W.size(); l-- > 0;)
          {
            const int fi = dims[l], fo = dims[l + 1];
            for (int j = 0; j < fo; ++j)
            {
              const double d = delta[l + 1][static_cast<std::size_t>(j)];
              g[bo[l] + static_cast<std::size_t>(j)] += d;
              for (int i = 0; i < fi; ++i)
              {
                g[wo[l] + static_cast<std::size_t>(i) * fo + j] +=
                  d * act[l][static_cast<std::size_t>(i)];
              }
            }
            if (l > 0)
            {
              delta[l].assign(static_cast<std::size_t>(fi), 0.0);
              for (int i = 0; i < fi; ++i)
              {
                double s = 0.0;
                for (int j = 0; j < fo; ++j)
                {
                  s += delta[l + 1][static_cast<std::size_t>(j)] *
                       W[l][static_cast<std::size_t>(i) * fo + j];
                }
                const double a = act[l][static_cast<std::size_t>(i)];
                delta[l][static_cast<std::size_t>(i)] = s * (1.0 - a * a);   // tanh'
              }
            }
          }
        }
      }

      // Chunk reduction in INDEX order -- serial and cheap (512 * 1.8k adds), and the thing that
      // makes the whole loop thread-count invariant.
      std::vector<double> gsum(np, 0.0);
      for (std::size_t c = 0; c < n_chunks; ++c)
      {
        const double* g = part.data() + c * np;
        for (std::size_t k = 0; k < np; ++k) { gsum[k] += g[k]; }
      }

      const double n = static_cast<double>(n_rows);
      for (std::size_t l = 0; l < W.size(); ++l)
      {
        for (std::size_t k = 0; k < W[l].size(); ++k)
        {
          W[l][k] -= p.lr * (gsum[wo[l] + k] / n + p.l2 * W[l][k]);
        }
        for (std::size_t k = 0; k < B[l].size(); ++k) { B[l][k] -= p.lr * (gsum[bo[l] + k] / n); }
      }
    }
  }


  std::vector<MLP> nets_;
  NNParams params_;
  int n_in_ = 0;
  bool trained_ = false;
};

} // namespace odia

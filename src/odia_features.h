#pragma once
// Flat, ODIA-owned scored features.
//
// WHY THIS EXISTS. Measured on the Astral benchmark, per precursor, from the component walk:
//
//   4.89 OpenMS Feature objects        x 296 B          =  1,448 B
//   71.3 subordinate Features          x 296 B          = 21,111 B
//   269 meta values in string-keyed flat_maps x ~24 B   =  6,459 B
//   76 separate MetaInfo heap allocations               =  2,439 B
//                                                        --------
//                                                         30.7 KB   and 76 allocations
//
// DIA-NN's counterpart, verified in diann_1.8.cpp:2924-2934, is
//
//   class PrecursorEntry { float scores[pN], decoy_scores[pN]; };   // pN = 110
//
// 880 B, contiguous, inline, with ZERO heap allocations and ZERO string keys -- a 35.7x difference
// in size and 76-to-0 in allocation count.
//
// That gap is a representation difference, not a tuning one, and it explains the fragmentation as
// well as the size: ~32M allocations in the feature map alone (on top of ~471M at library load),
// against live data that never exceeds 33 GB while the glibc arena reaches 186 GB.
//
// This is the columnar equivalent: one contiguous float block for all sub-scores, the score NAMES
// stored once rather than per feature, and no nested objects. It is deliberately NOT an OpenMS type
// -- CLEAN-ROOM keeps OpenMS unvendored, so ODIA owns its own scored-feature representation.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace odia
{

/// Sub-scores for many features, laid out row-major in ONE allocation.
/// Column names live in `names`, once, instead of being a string key on every value.
class FeatureTable
{
public:
  using Index = std::size_t;

  void setColumns(std::vector<std::string> names)
  {
    names_ = std::move(names);
  }
  const std::vector<std::string>& columns() const { return names_; }
  std::size_t columnCount() const { return names_.size(); }
  std::size_t size() const { return feature_id_.size(); }

  void reserve(std::size_t n)
  {
    feature_id_.reserve(n); precursor_id_.reserve(n);
    exp_rt_.reserve(n); delta_rt_.reserve(n); norm_rt_.reserve(n);
    left_width_.reserve(n); right_width_.reserve(n); exp_im_.reserve(n);
    scores_.reserve(n * names_.size());
  }

  /// Append one feature. `row` must hold exactly columnCount() values, in `columns()` order.
  /// Missing sub-scores are the caller's problem -- pass a sentinel; there is no per-row key set to
  /// consult, which is the point.
  void append(std::int64_t feature_id, std::int64_t precursor_id,
              float exp_rt, float delta_rt, float norm_rt,
              float left_width, float right_width, float exp_im,
              const float* row)
  {
    feature_id_.push_back(feature_id);
    precursor_id_.push_back(precursor_id);
    exp_rt_.push_back(exp_rt);
    delta_rt_.push_back(delta_rt);
    norm_rt_.push_back(norm_rt);
    left_width_.push_back(left_width);
    right_width_.push_back(right_width);
    exp_im_.push_back(exp_im);
    scores_.insert(scores_.end(), row, row + names_.size());
  }

  const float* row(Index i) const { return scores_.data() + i * names_.size(); }
  std::int64_t featureId(Index i) const { return feature_id_[i]; }
  std::int64_t precursorId(Index i) const { return precursor_id_[i]; }
  float expRt(Index i) const { return exp_rt_[i]; }
  float deltaRt(Index i) const { return delta_rt_[i]; }
  float normRt(Index i) const { return norm_rt_[i]; }

  const std::vector<std::int64_t>& featureIds() const { return feature_id_; }
  const std::vector<std::int64_t>& precursorIds() const { return precursor_id_; }

  /// Per-transition rows, kept separately because they are ~15x more numerous than features and
  /// carry no sub-scores -- the subordinate Features that hold this today are 21.1 KB of the
  /// 30.7 KB per precursor.
  void appendTransition(std::int64_t feature_id, std::int64_t transition_id,
                        float intensity, float apex_intensity)
  {
    t_feature_id_.push_back(feature_id);
    t_transition_id_.push_back(transition_id);
    t_intensity_.push_back(intensity);
    t_apex_.push_back(apex_intensity);
  }
  std::size_t transitionCount() const { return t_feature_id_.size(); }

  /// Bytes actually held, ignoring reserve slack. Comparable to the component walk's figures.
  std::size_t bytes() const
  {
    std::size_t b = 0;
    b += feature_id_.size() * sizeof(std::int64_t) * 2;                 // feature + precursor id
    b += exp_rt_.size() * sizeof(float) * 6;                            // the six scalar columns
    b += scores_.size() * sizeof(float);
    b += t_feature_id_.size() * (sizeof(std::int64_t) * 2 + sizeof(float) * 2);
    for (const auto& n : names_) { b += n.capacity() + sizeof(std::string); }
    return b;
  }

  /// What the same content costs as OpenMS Features with string-keyed MetaInfo, using the measured
  /// per-object sizes. For reporting the ratio, not for allocation.
  std::size_t bytesAsFeatureMap() const
  {
    constexpr std::size_t kFeature = 296;      // sizeof(OpenMS::Feature)
    constexpr std::size_t kMetaPair = 24;      // flat_map<UInt,DataValue> element, padded
    constexpr std::size_t kAllocHdr = 32;      // per separate MetaInfo allocation
    const std::size_t feats = feature_id_.size();
    const std::size_t subs = t_feature_id_.size();
    return feats * kFeature + subs * kFeature
         + feats * names_.size() * kMetaPair
         + (feats + subs) * kAllocHdr;
  }

  /// Allocation count: one block per column vector, versus one MetaInfo per feature AND per
  /// subordinate in the OpenMS representation.
  std::size_t allocationCount() const { return 12 + names_.size(); }
  std::size_t allocationCountAsFeatureMap() const
  {
    return feature_id_.size() + t_feature_id_.size();
  }

private:
  std::vector<std::string> names_;
  std::vector<std::int64_t> feature_id_, precursor_id_;
  std::vector<float> exp_rt_, delta_rt_, norm_rt_, left_width_, right_width_, exp_im_;
  std::vector<float> scores_;                 // size() * columnCount(), row-major, ONE allocation
  std::vector<std::int64_t> t_feature_id_, t_transition_id_;
  std::vector<float> t_intensity_, t_apex_;
};

} // namespace odia

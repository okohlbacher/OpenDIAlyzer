#pragma once
// Scored features, transitions and precursors, sized from what the columns ACTUALLY contain.
//
// WHY, measured. Releasing pass-1's FeatureMap drops in_use from 29.64 to 9.09 GB, so the OpenMS
// representation of these three tables is 20.55 GB live. Written out they are 0.95 GB of parquet.
// The 21x penalty is not mysterious: per precursor, 4.89 Feature objects at 296 B, 71.3
// subordinates at 296 B, 269 meta values in string-keyed flat_maps and 76 separate MetaInfo heap
// allocations. MetaInfoInterface is an 8-byte pointer to a heap MetaInfo, so the meta values are
// invisible to sizeof() and are the larger term.
//
// THREE THINGS THE COLUMN CONTENT SAYS, none of which a schema shows (measured over 400k rows of
// a real bundle -- see docs/feature-table-memory-plan.md):
//
//   1. MOST COLUMNS ARE EMPTY. feature_transition declares 44 columns and 35 of them are ALL NULL
//      on this run; features declares 65 and 16 are all-null with 2 more constant. They are ion
//      mobility and peak-shape columns that only PASEF or an enabled elution-model scorer fills.
//      Storing a column that is entirely null costs nothing here and is re-emitted as null by the
//      writer, so the schema is preserved without the bytes.
//   2. run_id IS CONSTANT. One value per run, stored once rather than 25.2M times.
//   3. THE JOIN KEYS ARE NOT DATA. feature_id repeats across every subordinate row purely to
//      re-associate rows that were contiguous to begin with. Grouped CSR-style (first, count) per
//      feature, it disappears entirely -- 21.9M int64 keys, 175 MB, for information already
//      implied by position.
//
// Per row, naive 8 B/column against the smallest LOSSLESS width for the values actually present:
//
//     features             520 -> 357 B/row
//     feature_transition   352 ->  53 B/row     (44 B once feature_id is implied by grouping)
//     feature_precursor     40 ->  22 B/row     (13 B likewise)
//
// which is ~1.74 GB against 20.55 GB, an 11.8x reduction WITH NO PRECISION CHANGE. Narrowing the
// remaining f64 score columns to f32 would reach ~1.05 GB, and is deliberately NOT done here:
// the spectrum store's float32 m/z has 0.059 ppm error against a 10 ppm window -- two orders of
// margin -- and still came back 22 peptides short on a full run, unexplained. Precision changes
// get their own change and their own measurement, never a ride on a representation swap.

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace odia
{

/// Which of a table's declared columns actually carry data in this run.
///
/// The writer must still emit every declared column so the output schema does not change between
/// runs; it emits null for anything absent here. Presence is decided ONCE, from the data, rather
/// than per row -- a per-row null mask on 21.9M rows would cost more than the columns it elides.
class ColumnSet
{
public:
  void declare(std::vector<std::string> all) { all_ = std::move(all); present_.assign(all_.size(), false); slot_.assign(all_.size(), kAbsent); }
  /// Mark a declared column as carrying data, and give it a storage slot.
  void present(const std::string& name)
  {
    for (std::size_t i = 0; i < all_.size(); ++i)
    {
      if (all_[i] == name)
      {
        if (!present_[i]) { present_[i] = true; slot_[i] = n_present_++; }
        return;
      }
    }
    throw std::invalid_argument("odia::ColumnSet: '" + name + "' is not a declared column");
  }
  bool isPresent(std::size_t i) const { return i < present_.size() && present_[i]; }
  /// Storage slot of declared column @p i, or kAbsent.
  std::uint32_t slot(std::size_t i) const { return i < slot_.size() ? slot_[i] : kAbsent; }
  std::size_t declared() const { return all_.size(); }
  std::size_t stored() const { return n_present_; }
  const std::vector<std::string>& names() const { return all_; }
  static constexpr std::uint32_t kAbsent = 0xFFFFFFFFu;

private:
  std::vector<std::string> all_;
  std::vector<bool> present_;
  std::vector<std::uint32_t> slot_;
  std::uint32_t n_present_ = 0;
};

/// One scored feature per row; the sub-scores in one contiguous block.
class FeatureRows
{
public:
  void setColumns(ColumnSet cs) { cols_ = std::move(cs); }
  const ColumnSet& columns() const { return cols_; }
  void setRunId(std::int64_t r) { run_id_ = r; }          ///< constant: stored once, not per row
  std::int64_t runId() const { return run_id_; }

  void reserve(std::size_t n)
  {
    feature_id_.reserve(n); precursor_id_.reserve(n);
    exp_rt_.reserve(n); delta_rt_.reserve(n); norm_rt_.reserve(n);
    left_.reserve(n); right_.reserve(n);
    scores_.reserve(n * cols_.stored());
  }

  /// `row` holds exactly columns().stored() values, in slot order.
  std::size_t append(std::int64_t feature_id, std::int64_t precursor_id,
                     float exp_rt, float delta_rt, float norm_rt,
                     float left_width, float right_width, const double* row)
  {
    feature_id_.push_back(feature_id);
    precursor_id_.push_back(precursor_id);
    exp_rt_.push_back(exp_rt); delta_rt_.push_back(delta_rt); norm_rt_.push_back(norm_rt);
    left_.push_back(left_width); right_.push_back(right_width);
    scores_.insert(scores_.end(), row, row + cols_.stored());
    return feature_id_.size() - 1;
  }

  std::size_t size() const { return feature_id_.size(); }
  std::int64_t featureId(std::size_t i) const { return feature_id_.at(i); }
  std::int64_t precursorId(std::size_t i) const { return precursor_id_.at(i); }
  float expRt(std::size_t i) const { return exp_rt_.at(i); }
  float deltaRt(std::size_t i) const { return delta_rt_.at(i); }
  float normRt(std::size_t i) const { return norm_rt_.at(i); }
  float leftWidth(std::size_t i) const { return left_.at(i); }
  float rightWidth(std::size_t i) const { return right_.at(i); }
  const double* scores(std::size_t i) const { return scores_.data() + i * cols_.stored(); }

  std::size_t bytes() const
  {
    return (feature_id_.capacity() + precursor_id_.capacity()) * sizeof(std::int64_t)
         + (exp_rt_.capacity() + delta_rt_.capacity() + norm_rt_.capacity()
            + left_.capacity() + right_.capacity()) * sizeof(float)
         + scores_.capacity() * sizeof(double);
  }
  double bytesPerRow() const { return size() ? double(bytes()) / double(size()) : 0.0; }
  void compact()
  {
    feature_id_.shrink_to_fit(); precursor_id_.shrink_to_fit();
    exp_rt_.shrink_to_fit(); delta_rt_.shrink_to_fit(); norm_rt_.shrink_to_fit();
    left_.shrink_to_fit(); right_.shrink_to_fit(); scores_.shrink_to_fit();
  }

private:
  ColumnSet cols_;
  std::int64_t run_id_ = 0;
  std::vector<std::int64_t> feature_id_, precursor_id_;
  std::vector<float> exp_rt_, delta_rt_, norm_rt_, left_, right_;
  std::vector<double> scores_;
};

/// Subordinate rows, grouped by their parent feature so the parent key is never stored.
///
/// This is the table that dominates: 21.9M rows against 2.07M features. Of its 44 declared
/// columns, 35 are all-null on this run and one (run_id) is constant, so nine remain -- and one of
/// those nine is the parent key, which grouping removes.
class TransitionRows
{
public:
  void reserve(std::size_t n_features, std::size_t n_rows)
  {
    first_.reserve(n_features + 1); first_.push_back(0);
    transition_id_.reserve(n_rows);
    area_.reserve(n_rows); total_area_.reserve(n_rows); apex_.reserve(n_rows);
    apex_rt_.reserve(n_rows); fwhm_.reserve(n_rows); masserr_.reserve(n_rows);
  }

  /// Append one subordinate to the feature currently being built.
  void append(std::uint32_t transition_id, float area, double total_area, double apex,
              float apex_rt, double fwhm, double masserror_ppm)
  {
    transition_id_.push_back(transition_id);
    area_.push_back(area); total_area_.push_back(total_area); apex_.push_back(apex);
    apex_rt_.push_back(apex_rt); fwhm_.push_back(fwhm); masserr_.push_back(masserror_ppm);
  }
  /// Close the current feature's run. Must be called once per feature, in feature order.
  void endFeature()
  {
    if (first_.empty()) { first_.push_back(0); }
    first_.push_back(static_cast<std::uint64_t>(transition_id_.size()));
  }

  std::size_t features() const { return first_.empty() ? 0 : first_.size() - 1; }
  std::size_t size() const { return transition_id_.size(); }
  std::uint64_t begin(std::size_t f) const { return first_.at(f); }
  std::uint64_t end(std::size_t f) const { return first_.at(f + 1); }
  std::uint32_t count(std::size_t f) const { return static_cast<std::uint32_t>(end(f) - begin(f)); }

  std::uint32_t transitionId(std::uint64_t r) const { return transition_id_.at(r); }
  float area(std::uint64_t r) const { return area_.at(r); }
  double totalArea(std::uint64_t r) const { return total_area_.at(r); }
  double apex(std::uint64_t r) const { return apex_.at(r); }
  float apexRt(std::uint64_t r) const { return apex_rt_.at(r); }
  double fwhm(std::uint64_t r) const { return fwhm_.at(r); }
  double massErrorPpm(std::uint64_t r) const { return masserr_.at(r); }

  std::size_t bytes() const
  {
    return first_.capacity() * sizeof(std::uint64_t)
         + transition_id_.capacity() * sizeof(std::uint32_t)
         + (area_.capacity() + apex_rt_.capacity()) * sizeof(float)
         + (total_area_.capacity() + apex_.capacity() + fwhm_.capacity() + masserr_.capacity())
           * sizeof(double);
  }
  double bytesPerRow() const { return size() ? double(bytes()) / double(size()) : 0.0; }
  void compact()
  {
    first_.shrink_to_fit(); transition_id_.shrink_to_fit();
    area_.shrink_to_fit(); total_area_.shrink_to_fit(); apex_.shrink_to_fit();
    apex_rt_.shrink_to_fit(); fwhm_.shrink_to_fit(); masserr_.shrink_to_fit();
  }

private:
  std::vector<std::uint64_t> first_;                     ///< CSR offsets: feature f owns [f, f+1)
  std::vector<std::uint32_t> transition_id_;
  std::vector<float> area_, apex_rt_;
  std::vector<double> total_area_, apex_, fwhm_, masserr_;
};

/// Precursor isotope rows, grouped by feature for the same reason.
class PrecursorRows
{
public:
  void reserve(std::size_t n_features, std::size_t n_rows)
  {
    first_.reserve(n_features + 1); first_.push_back(0);
    isotope_.reserve(n_rows); area_.reserve(n_rows); apex_.reserve(n_rows);
  }
  /// precursor_isotope is 0..3 on this data, so a byte holds it.
  void append(std::uint8_t isotope, float area, double apex)
  {
    isotope_.push_back(isotope); area_.push_back(area); apex_.push_back(apex);
  }
  void endFeature()
  {
    if (first_.empty()) { first_.push_back(0); }
    first_.push_back(static_cast<std::uint64_t>(isotope_.size()));
  }
  std::size_t features() const { return first_.empty() ? 0 : first_.size() - 1; }
  std::size_t size() const { return isotope_.size(); }
  std::uint64_t begin(std::size_t f) const { return first_.at(f); }
  std::uint64_t end(std::size_t f) const { return first_.at(f + 1); }
  std::uint8_t isotope(std::uint64_t r) const { return isotope_.at(r); }
  float area(std::uint64_t r) const { return area_.at(r); }
  double apex(std::uint64_t r) const { return apex_.at(r); }

  std::size_t bytes() const
  {
    return first_.capacity() * sizeof(std::uint64_t) + isotope_.capacity()
         + area_.capacity() * sizeof(float) + apex_.capacity() * sizeof(double);
  }
  double bytesPerRow() const { return size() ? double(bytes()) / double(size()) : 0.0; }
  void compact()
  {
    first_.shrink_to_fit(); isotope_.shrink_to_fit(); area_.shrink_to_fit(); apex_.shrink_to_fit();
  }

private:
  std::vector<std::uint64_t> first_;
  std::vector<std::uint8_t> isotope_;
  std::vector<float> area_;
  std::vector<double> apex_;
};

} // namespace odia

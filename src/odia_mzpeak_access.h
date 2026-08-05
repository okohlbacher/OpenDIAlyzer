// odia_mzpeak_access.h -- stream DIA data from mzPeak into OpenSWATH.
//
// WHY: the current input path (OpenMS SwathFile, readOptions "normal") materialises the whole
// run in RAM -- measured 500 GB - 2.0 TB peak RSS on a 1/8 library. mzPeak is a columnar,
// lazily-decoded store, so an adapter that only decodes what is asked for turns that into a
// bounded working set. The seam is exactly one interface: OpenSwathWorkflow consumes
// std::vector<OpenSwath::SwathMap>, each holding a shared_ptr<OpenSwath::ISpectrumAccess>.
//
// THE STRUCTURAL MISMATCH THIS FILE SOLVES
// mzPeak stores a diaPASEF FRAME as ONE spectrum carrying N isolation windows over disjoint
// ion-mobility ranges, plus a PER-PEAK mobility array. (mzPeakConverter does this deliberately:
// mzML has nowhere to put the mobility dimension and so splits a frame into N spectra, but
// mzPeak does, so the frame is kept whole with N precursors attached.) OpenSWATH instead
// expects ONE spectrum per (frame, isolation window) with parallel mz / intensity / drift-time
// ARRAYS. So the adapter emits one group per (frame, window) and, when a frame is shared,
// gives each window only the peaks inside its own mobility band -- otherwise the windows
// bleed into one another and, if only the first precursor is read, N-1 of them vanish.
//
// Metadata is read once up front WITHOUT decoding peaks (proven in src/mzpeak_probe.cpp);
// peak arrays are decoded only when a group is requested.
//
// Memory: O(one requested group), not O(run). Measured 156 MB peak RSS on a 13.7 GB run,
// against 500 GB - 2.0 TB for the SwathFile path.
//
// The older per-IM-slice layout (one spectrum per slice with a SCALAR ion_mobility()) is
// still handled: a group may hold several slice indices, and the scalar fills the drift array.
#pragma once

#include <mzpeak/open.h>

#include "odia_mmap_archive.h"
#include "odia_spectrumstore.h"
#include <mzpeak/index.h>
#include <mzpeak/spectra.h>
#include <mzpeak/spectrum.h>
#include <mzpeak/spectrum_metadata.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>   // OpenSwath::SwathMap

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <limits>
#include <mutex>
#include <optional>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace odia
{

/// One (frame, isolation-window) group: the mzPeak spectrum indices that make it up.
struct MzPeakGroup
{
  double rt = 0.0;                       ///< retention time (seconds)
  std::vector<uint64_t> slices;          ///< mzPeak spectrum indices, one per IM slice
  /// Mobility band (1/K0) owned by THIS isolation window inside its frame. A diaPASEF
  /// frame carries several windows over disjoint mobility ranges in ONE mzPeak spectrum,
  /// so a group's peaks must be restricted to its own band or the windows bleed into each
  /// other. Unbounded when the frame has a single window.
  double im_lo = -std::numeric_limits<double>::infinity();
  double im_hi =  std::numeric_limits<double>::infinity();
};

/// One thread's decoder: an Index and the Spectra over it, constructed in place. Neither type
/// is movable (Spectra binds its fetch callback to `this`), so they are built directly into
/// this aggregate and only ever accessed by reference.
/// THE ONE INDEX FOR THE PROCESS. Metadata is parsed once and the archive is opened once.
///
/// Measured on the Astral benchmark before this existed: one Index per worker thread meant 1,060
/// open file descriptors (each Reader opens the archive plus several parquet members) and 224
/// metadata-footer parses, and the run peaked at 127.6 GB -- against 74.5 GB for the SAME analysis
/// reading mzML. The columnar input, which should be the cheap path, was 71% more expensive than
/// XML because of how it was opened.
///
/// A function-local static is the right shape here and not merely convenient: MzPeak::Index has a
/// user-declared destructor and holds a unique_ptr, so it is neither movable nor copyable and
/// cannot be put in an optional, a vector, or on the heap from a returned prvalue. Copy-
/// initialising a static from the prvalue elides the move, and C++11 magic statics make the
/// initialisation thread-safe and exactly-once without a lock of our own.
///
/// s_path_ must be set before any worker calls this -- there is one input per run, so it is.
/// One mapping and one parsed central directory for the whole process.
inline const odia::MmapArchive& sharedArchive(const std::string& path)
{
  static odia::MmapArchive a(path);          // magic static: once, thread-safe, no lock of ours
  return a;
}

/// A fresh Archive over the SHARED mapping. Index takes ownership of an Archive, so each Index
/// needs its own object -- but the object is now a pair of shared_ptrs, not a file open and not a
/// directory scan.
inline std::unique_ptr<MzPeak::IO::Archive> archiveView(const std::string& path)
{
  const auto& a = sharedArchive(path);
  return std::make_unique<odia::MmapArchive>(a.mapping(), a.directory());
}

/// THE one Index. Metadata parsed once, over the shared mapping.
///
/// Centralising this ALONE was measured and was not enough: descriptors went 1,060 -> 911 because
/// the handles come from the per-thread Spectra, not from the Index. What fixes it is that every
/// File handed out below is a span of the mapping, so no Spectra can open anything.
inline const MzPeak::Index& sharedIndex(const std::string& path)
{
  static MzPeak::Index idx{archiveView(path)};
  return idx;
}

struct Reader
{
  explicit Reader(const std::string& path) : index(MzPeak::open(path)), spectra(index.spectra()) {}
  Reader(const Reader&) = delete;
  Reader& operator=(const Reader&) = delete;
  MzPeak::Index index;
  MzPeak::Spectra spectra;
};

/// Shared, immutable index over an mzPeak file: built once, shared by every clone and every
/// SWATH. Holds only metadata (no peak arrays).
class MzPeakIndex
{
public:
  explicit MzPeakIndex(const std::string& path)
    : path_(path), index_(MzPeak::open(path)), spectra_(index_.spectra())
  {
    s_path_ = path;                              // must precede any worker thread

    auto& spectra = spectra_;
    uint64_t i = 0;
    // key: isolation centre rounded to 1e-3 (MS1 -> a sentinel), then RT rounded to 1e-4
    std::map<long, std::map<long, MzPeakGroup>> by_window;
    for (auto sp : spectra)              // metadata-only sweep: never touches mz()/intensity()
    {
      const double rt = sp.retention_time() ? *sp.retention_time() : 0.0;
      const int level = sp.ms_level();
      if (level == 1)
      {
        auto& grp = by_window[kMs1Key][std::lround(rt * 10000.0)];
        grp.rt = rt;
        grp.slices.push_back(i);
        ++i;
        continue;
      }

      // A diaPASEF frame is ONE mzPeak spectrum carrying N isolation windows over
      // disjoint mobility ranges (mzPeakConverter keeps the frame whole because mzPeak,
      // unlike mzML, has somewhere to put the mobility dimension). Taking only the first
      // precursor here silently drops the other N-1 windows and misattributes their peaks
      // -- on S08 that is 12 of 24 windows. Collect them ALL.
      struct Win { double centre, lower, upper, im, im_lo, im_hi; };
      std::vector<Win> wins;
      for (const auto& p : sp.precursors())
      {
        if (!p.isolation_window.target_mz) { continue; }
        const double centre = *p.isolation_window.target_mz;
        // window half-widths are optional in the format; fall back to a symmetric guess
        const double lo = p.isolation_window.lower_offset ? *p.isolation_window.lower_offset : 0.0;
        const double up = p.isolation_window.upper_offset ? *p.isolation_window.upper_offset : 0.0;
        // Per-window mobility: the true [lower, upper] band when the writer records it,
        // plus the MIDPOINT (MS:1002815) that older writers emit alone. NaN when absent.
        double im = std::numeric_limits<double>::quiet_NaN();
        double im_lo = std::numeric_limits<double>::quiet_NaN();
        double im_hi = std::numeric_limits<double>::quiet_NaN();
        for (const auto& si : p.selected_ions)
        {
          if (si.ion_mobility_value) { im = *si.ion_mobility_value; }
          if (si.ion_mobility_lower_limit) { im_lo = *si.ion_mobility_lower_limit; }
          if (si.ion_mobility_upper_limit) { im_hi = *si.ion_mobility_upper_limit; }
          if (si.ion_mobility_value || si.ion_mobility_lower_limit) { break; }
        }
        wins.push_back({centre, centre - lo, centre + up, im, im_lo, im_hi});
      }
      if (wins.empty()) { ++i; continue; }             // MS2 without isolation info -> skip

      // TRUE bounds if the writer recorded them. They are required, not a nicety: a frame's
      // windows cover ASYMMETRIC scan ranges, so splitting at the midpoint between adjacent
      // window centres puts the boundary in the wrong place -- measured against the vendor
      // analysis.tdf on S08, 3.5% of the mobility axis on average and 10.3% at worst lands in
      // a window that never fragmented it. The vendor scan->mobility calibration is also
      // non-linear, so averaging in mobility space is wrong even for equal-width windows.
      const bool have_bounds = std::all_of(wins.begin(), wins.end(), [](const Win& w) {
        return std::isfinite(w.im_lo) && std::isfinite(w.im_hi) && w.im_lo < w.im_hi;
      });
      const bool have_im = std::all_of(wins.begin(), wins.end(),
                                       [](const Win& w) { return std::isfinite(w.im); });
      std::vector<std::size_t> order(wins.size());
      for (std::size_t k = 0; k < wins.size(); ++k) { order[k] = k; }
      if (have_im)
      {
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return wins[a].im < wins[b].im; });
      }
      if (!have_bounds && wins.size() > 1 && !im_split_warned_)
      {
        im_split_warned_ = true;
        std::fprintf(stderr,
                     "OpenDIAlyzer/mzPeak: frames carry %zu isolation windows but the file has no "
                     "per-window ion-mobility BOUNDS%s. Re-convert with a writer that records "
                     "them; results from this file are approximate.\n",
                     wins.size(),
                     have_im ? " (falling back to a midpoint split, which misplaces the boundary "
                               "between asymmetric windows)"
                             : " and no midpoints either, so windows cannot be separated at all");
      }
      // Ambiguous geometry: two windows claiming the same midpoint cannot be told apart, and
      // silently handing one an empty band would look like a clean run.
      if (!have_bounds && have_im && wins.size() > 1)
      {
        for (std::size_t k = 1; k < order.size(); ++k)
        {
          if (wins[order[k]].im == wins[order[k - 1]].im && !im_tie_warned_)
          {
            im_tie_warned_ = true;
            std::fprintf(stderr, "OpenDIAlyzer/mzPeak: two isolation windows share an ion-mobility "
                                 "midpoint; their peaks cannot be separated.\n");
          }
        }
      }

      for (std::size_t k = 0; k < order.size(); ++k)
      {
        const Win& w = wins[order[k]];
        const long wkey = std::lround(w.centre * 1000.0);
        auto& grp = by_window[wkey][std::lround(rt * 10000.0)];
        grp.rt = rt;
        grp.slices.push_back(i);
        if (have_bounds)
        {
          grp.im_lo = w.im_lo; grp.im_hi = w.im_hi;
        }
        else if (have_im && order.size() > 1)
        {
          if (k > 0) { grp.im_lo = 0.5 * (wins[order[k - 1]].im + w.im); }
          if (k + 1 < order.size()) { grp.im_hi = 0.5 * (w.im + wins[order[k + 1]].im); }
        }
        if (windows_.find(wkey) == windows_.end()) { windows_[wkey] = {w.centre, w.lower, w.upper}; }
        // SwathMap mobility limits. OpenSWATH assigns transitions with STRICT
        // imLower < precursorIM < imUpper, so these must be the window's real span. Seeding
        // from a default-constructed pair would give [0, midpoint] and reject the upper half
        // of every window, so track the range explicitly.
        if (have_bounds)
        {
          auto [it, fresh] = im_range_.try_emplace(wkey, w.im_lo, w.im_hi);
          if (!fresh)
          {
            it->second.first = std::min(it->second.first, w.im_lo);
            it->second.second = std::max(it->second.second, w.im_hi);
          }
        }
      }
      ++i;
    }
    n_spectra_total_ = i;
    // flatten each window's RT map into an RT-ordered vector (OpenSWATH assumes RT order)
    for (auto& w : by_window)
    {
      auto& v = groups_[w.first];
      v.reserve(w.second.size());
      for (auto& g : w.second) { v.push_back(std::move(g.second)); }
      std::sort(v.begin(), v.end(), [](const MzPeakGroup& a, const MzPeakGroup& b) { return a.rt < b.rt; });
    }
  }

  static constexpr long kMs1Key = -1;

  struct Window { double centre, lower, upper; };

  const std::vector<MzPeakGroup>& groups(long wkey) const
  {
    static const std::vector<MzPeakGroup> kEmpty;
    auto it = groups_.find(wkey);
    return it == groups_.end() ? kEmpty : it->second;
  }
  std::vector<long> windowKeys() const
  {
    std::vector<long> k;
    for (const auto& w : groups_) { if (w.first != kMs1Key) { k.push_back(w.first); } }
    return k;
  }
  bool hasMs1() const { return groups_.count(kMs1Key) > 0; }
  Window window(long wkey) const { auto it = windows_.find(wkey); return it == windows_.end() ? Window{-1,-1,-1} : it->second; }
  std::pair<double,double> imRange(long wkey) const
  {
    auto it = im_range_.find(wkey);
    if (it == im_range_.end() || it->second.first > it->second.second) { return {-1.0, -1.0}; }
    return it->second;
  }
  uint64_t nSpectraTotal() const { return n_spectra_total_; }
  const std::string& path() const { return path_; }

  /// Decode one group into a fresh OpenSWATH spectrum (mz / intensity / drift arrays).
  /// Thread-safe: mzPeak decoding is serialised here because the reader is not documented
  /// as concurrently reentrant. Chunked callers therefore parallelise over *work*, not IO.
  /// Decode EVERY group of every window once, into per-window SpectrumStores.
  ///
  /// The run walks the spectra five times (prefilter targets, prefilter decoys, calibration,
  /// pass 1, pass 2). mzML parsed once and left them resident, so passes 2-5 were walks over RAM;
  /// mzPeak streams, so each pass genuinely re-decoded parquet -- the same analysis went from
  /// 20:16 to over 83 minutes still inside the prefilter. Decoding once here restores the property
  /// mzML had by accident, at 8 B/peak instead of 16.
  ///
  /// Parallel over WINDOWS, serial within one: each store is filled by exactly one thread, so no
  /// locking and no shared append. ~150 windows over the available workers.
  void populate(int threads)
  {
    std::vector<long> keys;
    keys.reserve(windows_.size() + 1);
    for (const auto& kv : groups_) { keys.push_back(kv.first); }
    stores_.clear();
    for (long k : keys) { stores_.emplace(k, odia::SpectrumStore{}); }

    std::size_t n_spec = 0, n_peak = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::max(1, threads)) \
        reduction(+ : n_spec, n_peak)
#endif
    for (std::size_t ki = 0; ki < keys.size(); ++ki)
    {
      const long wkey = keys[ki];
      auto& store = stores_.at(wkey);
      const auto& gs = groups_.at(wkey);
      store.reserve(gs.size(), 0);
      for (std::size_t gi = 0; gi < gs.size(); ++gi)
      {
        OpenSwath::SpectrumPtr sp = decodeRaw(wkey, gi);
        odia::SpectrumStore::Meta m;
        m.rt = gs[gi].rt;
        m.ms_level = (wkey == kMs1Key) ? 1 : 2;
        const auto wit = windows_.find(wkey);
        if (wit != windows_.end())
        {
          m.iso_lower = wit->second.lower;
          m.iso_upper = wit->second.upper;
          m.precursor_mz = 0.5 * (wit->second.lower + wit->second.upper);
        }
        const auto& mzv = sp->getMZArray()->data;
        const auto& inv = sp->getIntensityArray()->data;
        std::vector<float> fi(inv.size());
        for (std::size_t k = 0; k < inv.size(); ++k) { fi[k] = static_cast<float>(inv[k]); }
        const auto dtp = sp->getDriftTimeArray();
        const bool has_dt = dtp && dtp->data.size() == mzv.size() &&
                            std::any_of(dtp->data.begin(), dtp->data.end(),
                                        [](double x) { return x >= 0.0; });
        store.add(mzv.data(), fi.data(), static_cast<std::uint32_t>(mzv.size()), m,
                  has_dt ? dtp->data.data() : nullptr);
        n_spec += 1;
        n_peak += mzv.size();
      }
    }
    populated_ = true;
    std::size_t b = 0;
    for (const auto& kv : stores_) { b += kv.second.bytes(); }
    std::fprintf(stderr,
                 "OpenDIAlyzer[mzpeak] decoded %zu spectra, %zu peaks ONCE into %.2f GB "
                 "(%.2f GB as double pairs, %.2fx); later passes read memory\n",
                 n_spec, n_peak, b / 1073741824.0, n_peak * 16.0 / 1073741824.0,
                 (n_peak * 16.0) / std::max<double>(1.0, double(b)));
  }

  bool populated() const { return populated_; }

  /// Serve from the store when it is populated; otherwise decode.
  OpenSwath::SpectrumPtr decode(long wkey, std::size_t gi)
  {
    if (populated_)
    {
      const auto it = stores_.find(wkey);
      if (it != stores_.end() && gi < it->second.size())
      {
        OpenSwath::SpectrumPtr out(new OpenSwath::Spectrum);
        OpenSwath::BinaryDataArrayPtr mz(new OpenSwath::BinaryDataArray);
        OpenSwath::BinaryDataArrayPtr in(new OpenSwath::BinaryDataArray);
        OpenSwath::BinaryDataArrayPtr dt(new OpenSwath::BinaryDataArray);
        it->second.widen(gi, mz->data, in->data, &dt->data);
        out->setMZArray(mz);
        out->setIntensityArray(in);
        out->setDriftTimeArray(dt);
        return out;
      }
    }
    return decodeRaw(wkey, gi);
  }

  OpenSwath::SpectrumPtr decodeRaw(long wkey, std::size_t gi)
  {
    OpenSwath::SpectrumPtr out(new OpenSwath::Spectrum);
    OpenSwath::BinaryDataArrayPtr mz(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr in(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr dt(new OpenSwath::BinaryDataArray);
    const auto& gs = groups(wkey);
    if (gi >= gs.size())
    {
      out->setMZArray(mz); out->setIntensityArray(in);
      return out;
    }
    {
      // PER-THREAD DECODER. The index (groups_/windows_) is immutable and shared, but decoding
      // used to run under one global mutex -- measured at exactly 1 core of 224 while the .d
      // path reached ~14, making the streaming input 12x slower despite using 250x less memory.
      //
      // Each thread instead opens its own Index+Spectra over the same file: mzPeak decode state
      // is per-Spectra, so nothing is shared and no lock is needed. Cost is one file handle and
      // one metadata footer per worker.
      //
      // These MUST be thread_local OBJECTS, not unique_ptrs: MzPeak::Spectra binds its fetch
      // callback to `this`, so heap-allocating it from a returned temporary would move it and
      // leave the callback pointing at the moved-from object (this exact mistake produced a
      // segfault earlier). Copy-initialising a thread_local from a prvalue elides the move.
      // MzPeak::Index is neither copyable nor movable, so it cannot go in an optional or on
      // the heap -- it must be COPY-INITIALISED from the prvalue, where the move is elided.
      // The path therefore comes from a static set before any worker starts (one input per run).
      // One Reader per thread, for the life of the thread.
      //
      // TESTED AND REJECTED: recycling this Reader every N decodes, on the theory that the
      // reader accumulates decode state. It does not help -- RSS still climbed 1.7 -> 9.8 GB
      // over 17 min with recycling on, and SUPERLINEARLY, which is the shape of work completed
      // rather than of reader lifetime. The growth is OpenSwathWorkflow accumulating extracted
      // chromatograms and features in memory; it is independent of the input backend and is not
      // something this adapter can fix. Do not re-add recycling without first measuring where
      // the memory actually goes (massif/heaptrack), or it is complexity for nothing.
      // CENTRALISED METADATA, PER-THREAD STREAMING. Index::spectra() is const, so one shared Index
      // serves every worker: the archive is opened once and the footer parsed once, while each
      // thread still gets its own Spectra and therefore its own decode state, so no lock is needed
      // and the 1-core-of-224 regression that per-thread Readers were introduced to avoid does not
      // come back. What goes away is 1,060 file descriptors and 223 redundant metadata parses.
      //
      // Spectra is non-copyable AND non-movable (it binds its fetch callback to `this`), so it is
      // COPY-INITIALISED from the prvalue, where the move is elided -- the same reason the Index
      // above is a static rather than an optional. It must not be wrapped in optional/unique_ptr.
      thread_local MzPeak::Spectra t_spectra = sharedIndex(s_path_).spectra();
      // Use the LONG-LIVED Spectra (spectra_) rather than a fresh index_.spectra() per call.
      // Rationale: MzPeak::Spectra binds its fetch callback to `this`
      // (`std::bind(std::mem_fn(&Spectra::fetch), this, _1)`, src/spectra.cpp), so any Spectra
      // that is copied or moved carries a callback aimed at the *source* object. Guaranteed
      // copy-elision may make the naive form work by accident, which is exactly why it should
      // not be relied on. UNVERIFIED as the cause of the observed "0 peaks decoded" -- the
      // iterator-vs-operator[] isolation probe (scratch: mzp_iso.cpp) settles that; this change
      // is sound either way and also avoids reconstructing Spectra on every call.
      auto& spectra = t_spectra;   // thread_local object, not a pointer
      for (uint64_t si : gs[gi].slices)
      {
        // Spectra::fetch() is private; EnumerableProxy::operator[] is the public random-access
        // entry point (it forwards to the same lazy fetch).
        auto sp = spectra[static_cast<std::size_t>(si)];
        // UPSTREAM LIMITATION: Spectrum::intensity() is typed `const std::vector<float>&` and
        // throws MzPeak::TypeError ("Expected float or double but got: int32") when the file
        // stores intensities as integers -- which mzML permits and real files use. There is no
        // public type-agnostic accessor, so we can only fail that spectrum loudly-once rather
        // than crash the whole run. Fix belongs in mzPeak (a variant/lifted accessor).
        try
        {
          const std::vector<double>& m = sp.mz();
          const std::vector<float>& y = sp.intensity();
          const std::size_t n = std::min(m.size(), y.size());
          // PER-PEAK ion mobility where the file has it. Both writers seen so far emit one
          // spectrum per FRAME carrying a full mobility ARRAY -- collapsing that to the
          // spectrum's scalar ion_mobility() would give every peak the same drift time and
          // throw away the IM separation diaPASEF exists for. The scalar path remains for
          // per-slice writers, where it is the correct value for the whole slice.
          const std::vector<double>& imv = sp.ion_mobility_array();
          const bool per_peak_im = imv.size() >= n;
          const double lo = gs[gi].im_lo, hi = gs[gi].im_hi;
          // A frame shared by several isolation windows must contribute only the peaks in
          // THIS window's mobility band; the rest belong to its siblings. Half-open [lo, hi)
          // so adjacent bands partition the axis without double-counting a boundary peak.
          const bool banded = per_peak_im && (std::isfinite(lo) || std::isfinite(hi));
          for (std::size_t k = 0; k < n; ++k)
          {
            // NaN fails every comparison, so a naive range test would KEEP a NaN-mobility
            // peak in every sibling window and break the partition. Require finiteness.
            if (banded && !(imv[k] >= lo && imv[k] < hi)) { continue; }
            mz->data.push_back(m[k]);
            in->data.push_back(static_cast<double>(y[k]));
            dt->data.push_back(per_peak_im ? imv[k]
                                           : (sp.ion_mobility() ? *sp.ion_mobility() : -1.0));
          }
        }
        catch (const std::exception& e)
        {
          if (!decode_error_reported_)
          {
            decode_error_reported_ = true;
            std::fprintf(stderr, "OpenDIAlyzer/mzPeak: cannot decode spectrum %llu (%s). "
                                 "Unsupported intensity encoding -- needs a type-agnostic "
                                 "accessor in mzPeak.\n", (unsigned long long) si, e.what());
          }
        }
      }
    }
    // OpenSWATH REQUIRES m/z-sorted peaks (extractChromatograms throws otherwise). Concatenated
    // IM slices are not globally sorted, so sort the three arrays together by m/z.
    const std::size_t n = mz->data.size();
    std::vector<std::size_t> ord(n);
    for (std::size_t k = 0; k < n; ++k) { ord[k] = k; }
    std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return mz->data[a] < mz->data[b]; });
    OpenSwath::BinaryDataArrayPtr mz2(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr in2(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr dt2(new OpenSwath::BinaryDataArray);
    mz2->data.reserve(n); in2->data.reserve(n); dt2->data.reserve(n);
    for (std::size_t k : ord)
    {
      mz2->data.push_back(mz->data[k]);
      in2->data.push_back(in->data[k]);
      dt2->data.push_back(dt->data[k]);
    }
    out->setMZArray(mz2);
    out->setIntensityArray(in2);
    out->setDriftTimeArray(dt2);
    return out;
  }

private:
  std::string path_;
  // Read by the thread_local decoder initialisers, which run on first decode in each worker.
  // Set once in the constructor, before any extraction thread exists.
  static inline std::string s_path_;
  MzPeak::Index index_;
  MzPeak::Spectra spectra_;      ///< index-build only; decode() uses a thread_local Spectra
  std::map<long, odia::SpectrumStore> stores_;   ///< per window; filled once by populate()
  bool populated_ = false;
  std::map<long, std::vector<MzPeakGroup>> groups_;
  std::map<long, Window> windows_;
  std::map<long, std::pair<double,double>> im_range_{};
  uint64_t n_spectra_total_ = 0;
  bool im_split_warned_ = false;
  bool im_tie_warned_ = false;
  // (io_mutex_ removed: decoding is per-thread, so there is nothing left to serialise)
  bool decode_error_reported_ = false;
};

using MzPeakIndexPtr = std::shared_ptr<MzPeakIndex>;

/// ISpectrumAccess over ONE isolation window of an mzPeak file.
/// lightClone() shares the index (cheap, as OpenSWATH clones this per scoring job).
class MzPeakSpectrumAccess : public OpenSwath::ISpectrumAccess
{
public:
  MzPeakSpectrumAccess(MzPeakIndexPtr idx, long wkey) : idx_(std::move(idx)), wkey_(wkey) {}
  ~MzPeakSpectrumAccess() override = default;

  std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override
  {
    return std::make_shared<MzPeakSpectrumAccess>(idx_, wkey_);
  }

  OpenSwath::SpectrumPtr getSpectrumById(int id) override
  {
    if (id < 0) { return OpenSwath::SpectrumPtr(new OpenSwath::Spectrum); }
    return idx_->decode(wkey_, static_cast<std::size_t>(id));
  }

  OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override
  {
    OpenSwath::SpectrumMeta m;
    const auto& g = idx_->groups(wkey_);
    if (id >= 0 && static_cast<std::size_t>(id) < g.size())
    {
      m.RT = g[id].rt;
      m.ms_level = (wkey_ == MzPeakIndex::kMs1Key) ? 1 : 2;
      m.id = std::to_string(id);
    }
    return m;
  }

  std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override
  {
    const auto& g = idx_->groups(wkey_);
    std::vector<std::size_t> out;
    if (g.empty()) { return out; }
    if (deltaRT <= 0.0)                                  // nearest single spectrum
    {
      std::size_t best = 0; double bd = std::abs(g[0].rt - RT);
      for (std::size_t i = 1; i < g.size(); ++i)
      {
        const double d = std::abs(g[i].rt - RT);
        if (d < bd) { bd = d; best = i; }
      }
      out.push_back(best);
      return out;
    }
    for (std::size_t i = 0; i < g.size(); ++i)           // groups are RT-sorted
    {
      if (g[i].rt < RT - deltaRT) { continue; }
      if (g[i].rt > RT + deltaRT) { break; }
      out.push_back(i);
    }
    return out;
  }

  std::size_t getNrSpectra() const override { return idx_->groups(wkey_).size(); }

  // DIA input carries no chromatograms; these exist only to satisfy the interface.
  OpenSwath::ChromatogramPtr getChromatogramById(int /*id*/) override
  { return OpenSwath::ChromatogramPtr(new OpenSwath::Chromatogram); }
  std::size_t getNrChromatograms() const override { return 0; }
  std::string getChromatogramNativeID(int /*id*/) const override { return ""; }

private:
  MzPeakIndexPtr idx_;
  long wkey_;
};

/// Build the SwathMap vector OpenSwathWorkflow consumes, straight from an mzPeak file.
/// Mirrors what SwathFile::loadBrukerTdf produces, but streaming.
inline std::vector<OpenSwath::SwathMap> loadMzPeakSwathMaps(const std::string& path,
                                                            MzPeakIndexPtr& idx_out)
{
  idx_out = std::make_shared<MzPeakIndex>(path);
  std::vector<OpenSwath::SwathMap> maps;
  if (idx_out->hasMs1())
  {
    OpenSwath::SwathMap m;
    m.sptr = std::make_shared<MzPeakSpectrumAccess>(idx_out, MzPeakIndex::kMs1Key);
    m.ms1 = true;
    m.lower = m.upper = m.center = 0.0;
    m.imLower = m.imUpper = -1.0;
    maps.push_back(m);
  }
  for (long wkey : idx_out->windowKeys())
  {
    const auto w = idx_out->window(wkey);
    const auto im = idx_out->imRange(wkey);
    OpenSwath::SwathMap m;
    m.sptr = std::make_shared<MzPeakSpectrumAccess>(idx_out, wkey);
    m.ms1 = false;
    m.center = w.centre;
    m.lower = (w.lower >= 0) ? w.lower : w.centre;
    m.upper = (w.upper >= 0) ? w.upper : w.centre;
    m.imLower = im.first;
    m.imUpper = im.second;
    maps.push_back(m);
  }
  return maps;
}

} // namespace odia

// OpenDIALibGen — Stage 1: peptide list -> predicted OpenSWATH transition TSV.
//
// One job only: read a TSV of (PeptideSequence, PrecursorCharge), predict RT
// and MS2 fragment intensities with the vendored PeptDeep ONNX models, pair
// each predicted intensity to its theoretical b/y fragment m/z, and write the
// OpenSWATH transition TSV that experiments/diann_lib_to_openswath.py emits
// for the DIA-NN arm (same column set, so both arms feed OpenSWATH the same
// way).
//
// Explicitly NOT here (later stages): FASTA digest, modifications, CCS,
// decoys (feed this output to OpenSwathDecoyGenerator), assay refinement.
// Prediction parity vs AlphaPeptDeep Python is inherited from OpenMS's own
// PeptDeepInference_test (RT max err 5.96e-08, MS2 7.15e-07); this file owns
// only the assembly logic, which --selftest covers as a shape/sanity check.

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>
// Stage 2: FASTA digest + in-process assay refinement + decoy generation,
// folding in what OpenSwathAssayGenerator / OpenSwathDecoyGenerator do via the
// same library classes (MRMAssay / MRMDecoy).
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
// --- Fixed prediction conditions -------------------------------------------
// AlphaPeptDeep's pretrained MS2 model is conditioned on (charge, NCE,
// instrument); wrong values give plausible-but-wrong spectra.
// NCE 30: AlphaPeptDeep's default for the pretrained ms2 model
// (peptdeep pretrained_models, default nce=30, instrument='QE'), and the
// central value in OpenMS's parity reference data
// (proteomicsml_test_data_ms2_spectra.csv, nce range 25.03-32.5).
// Instrument index 0 = QE: that reference data pairs instrument_index 0 with
// instrument "QE" and PASSES PeptDeepInference_test against AlphaPeptDeep
// Python ONNX output, so 0->QE is the empirically verified encoding.
// NOTE: the docstring in PeptDeepMS2Inference.h claims 0=Lumos; the parity
// data contradicts it and is the only source backed by a passing test.
// REVIEWER TODO: cross-check against AlphaPeptDeep's instrument_dict for the
// pinned checkpoint before Stage 2.
constexpr float kNCE = 30.0f;
constexpr int64_t kInstrument = 0; // QE, per the OpenMS parity reference data

// PeptDeep MS2 output layout (from PeptDeepInference_test and the reference
// CSV's fragment_position/ion_index/ion_type columns): per peptide a flat
// array of (nAA-1)*8 floats, fragment-major — index = fragment_position*8 +
// ion_index, fragment_position = 0-based peptide-bond index, ion_index:
//   0=b_z1  1=b_z2  2=y_z1  3=y_z2  4..7=b/y modloss (zero for unmodified)
// Pairing rule to TheoreticalSpectrumGenerator annotations:
//   b-series number s  <-> fragment_position s-1
//   y-series number s  <-> fragment_position nAA-1-s
// predictMS2() already base-peak-normalizes to 1.0 and zeroes negatives,
// values < 1e-4, and b1_z1 (flat index 0).
constexpr int kIonChannels = 8;

// Column set: exactly what diann_lib_to_openswath.py emits (its RENAME
// targets + synthesized columns), minus GeneName (no proteins on a bare
// peptide list) and PrecursorIonMobility (CCS is out of Stage 1). Covers the
// converter's REQUIRED list.
const char* kHeader =
    "PrecursorMz\tProductMz\tPrecursorCharge\tProductCharge\tLibraryIntensity\t"
    "NormalizedRetentionTime\tPeptideSequence\tModifiedPeptideSequence\tProteinId\t"
    "FragmentType\tFragmentSeriesNumber\tTransitionGroupId\tTransitionId\tDecoy\t"
    "DetectingTransition\tIdentifyingTransition\tQuantifyingTransition";

struct Precursor
{
  std::string sequence; // unmodified, uppercase
  int charge;
  std::string protein = "UNKNOWN"; // FASTA identifier in Stage 2; UNKNOWN for bare peptide lists
};

// Stage 1 accepts only what the pinned unmodified checkpoint was validated
// on: canonical uppercase AAs, length >= 2 (MS2 has nAA-1 fragment bonds),
// charge >= 1.
bool valid_precursor(const Precursor& p)
{
  if (p.sequence.size() < 2 || p.charge < 1) return false;
  return p.sequence.find_first_not_of("ACDEFGHIKLMNPQRSTVWY") == std::string::npos;
}

// The whole pipeline, into a string so --selftest can diff two runs for the
// byte-identical-output requirement without touching the filesystem.
std::string generate(const std::vector<Precursor>& precursors,
                     const std::string& model_dir, int threads)
{
  std::vector<std::string> peptides;
  std::vector<float> charges, nces;
  std::vector<int64_t> instruments;
  peptides.reserve(precursors.size());
  for (const auto& p : precursors)
  {
    peptides.push_back(p.sequence);
    charges.push_back(static_cast<float>(p.charge));
    nces.push_back(kNCE);
    instruments.push_back(kInstrument);
  }

  OpenMS::PeptDeepRTInference rt_model(model_dir + "/peptdeep_rt_dynamic.onnx", threads);
  OpenMS::PeptDeepMS2Inference ms2_model(model_dir + "/peptdeep_ms2_dynamic.onnx", threads);
  const std::vector<float> rts = rt_model.predictRT(peptides);
  const std::vector<std::vector<float>> ms2 =
      ms2_model.predictMS2(peptides, charges, nces, instruments);

  // Defaults already give plain b/y ions only (no losses, isotopes, precursor
  // or immonium peaks). Two overrides: add_metainfo is what makes the
  // annotation pairing possible — peak i carries ion name "b3+"/"y5++" in
  // StringDataArray[0]; add_first_prefix_ion adds b1, so EVERY predicted
  // channel (including b1_z2) has a theoretical m/z partner. (b1_z1 is
  // generated but never emitted: predictMS2 forces its intensity to 0.)
  OpenMS::TheoreticalSpectrumGenerator tsg;
  OpenMS::Param tsg_param = tsg.getParameters();
  tsg_param.setValue("add_metainfo", "true");
  tsg_param.setValue("add_first_prefix_ion", "true");
  tsg.setParameters(tsg_param);

  std::string out = kHeader;
  out += '\n';
  char line[512];

  for (size_t i = 0; i < precursors.size(); ++i)
  {
    const std::string& seq = precursors[i].sequence;
    const int z = precursors[i].charge;
    const int nAA = static_cast<int>(seq.size());
    const OpenMS::AASequence aaseq = OpenMS::AASequence::fromString(seq);
    const double prec_mz = aaseq.getMonoWeight(OpenMS::Residue::Full, z) / z;

    // PeptDeep predicts z1/z2 channels only; a fragment charge above the
    // precursor charge is physically implausible (AlphaPeptDeep's library
    // builder caps fragment charge at min(precursor charge, 2) the same way).
    const int frag_zmax = std::min(z, 2);

    OpenMS::PeakSpectrum spec;
    tsg.getSpectrum(spec, aaseq, 1, frag_zmax, z);
    if (spec.getStringDataArrays().empty())
    {
      throw std::runtime_error("TheoreticalSpectrumGenerator returned no ion annotations");
    }
    const auto& names = spec.getStringDataArrays()[0];

    const std::string tg = seq + "_" + std::to_string(z);

    // Collect (ProductMz, intensity, ion) rows, then emit sorted by ProductMz
    // — fixed formatting + fixed order is what makes the output byte-stable.
    struct Row
    {
      double mz;
      float intensity;
      char type; // 'b' or 'y'
      int series;
      int fz;
      bool operator<(const Row& o) const
      {
        if (mz != o.mz) return mz < o.mz;
        if (type != o.type) return type < o.type;
        if (series != o.series) return series < o.series;
        return fz < o.fz;
      }
    };
    std::vector<Row> rows;

    for (size_t k = 0; k < spec.size(); ++k)
    {
      const std::string& name = names[k];
      const size_t plus = name.find('+');
      if (plus == std::string::npos || plus < 2) continue; // not a plain b/y annotation
      const char type = name[0];
      if (type != 'b' && type != 'y') continue;
      const int series = std::stoi(name.substr(1, plus - 1));
      const int fz = static_cast<int>(name.size() - plus); // count of '+'
      if (fz < 1 || fz > 2) continue;

      const int frag_pos = (type == 'b') ? series - 1 : nAA - 1 - series;
      const int channel = (type == 'b') ? fz - 1 : 2 + fz - 1;
      const size_t idx = static_cast<size_t>(frag_pos) * kIonChannels + channel;
      if (frag_pos < 0 || idx >= ms2[i].size()) continue;

      const float intensity = ms2[i][idx];
      if (intensity <= 0.0f) continue; // zeroed by predictMS2 (floor 1e-4)
      rows.push_back({spec[k].getMZ(), intensity, type, series, fz});
    }
    std::sort(rows.begin(), rows.end());

    // TransitionId mirrors the converter's groupby cumcount: tg + "_" + a
    // 0-based running index in emission order.
    for (size_t t = 0; t < rows.size(); ++t)
    {
      const Row& r = rows[t];
      const std::string tid = tg + "_" + std::to_string(t);
      std::snprintf(line, sizeof(line),
                    "%.6f\t%.6f\t%d\t%d\t%.6f\t%.4f\t%s\t%s\t%s\t%c\t%d\t%s\t%s\t0\t1\t0\t1\n",
                    prec_mz, r.mz, z, r.fz, static_cast<double>(r.intensity),
                    static_cast<double>(rts[i]), seq.c_str(), seq.c_str(),
                    precursors[i].protein.c_str(),
                    r.type, r.series, tg.c_str(), tid.c_str());
      out += line;
    }
  }
  return out;
}

std::vector<std::string> split_tab(const std::string& s)
{
  std::vector<std::string> f;
  size_t start = 0;
  for (;;)
  {
    const size_t tab = s.find('\t', start);
    if (tab == std::string::npos) { f.push_back(s.substr(start)); break; }
    f.push_back(s.substr(start, tab - start));
    start = tab + 1;
  }
  return f;
}

int read_precursors(const char* path, std::vector<Precursor>& out)
{
  std::ifstream in(path);
  if (!in) { std::fprintf(stderr, "error: cannot open %s\n", path); return 1; }
  std::string header;
  if (!std::getline(in, header))
  {
    std::fprintf(stderr, "error: %s is empty\n", path);
    return 1;
  }
  const std::vector<std::string> cols = split_tab(header);
  int iseq = -1, ichg = -1;
  for (size_t i = 0; i < cols.size(); ++i)
  {
    if (cols[i] == "PeptideSequence") iseq = static_cast<int>(i);
    if (cols[i] == "PrecursorCharge") ichg = static_cast<int>(i);
  }
  if (iseq < 0 || ichg < 0)
  {
    std::fprintf(stderr,
                 "error: %s needs columns PeptideSequence and PrecursorCharge\n", path);
    return 1;
  }

  std::string row;
  size_t lineno = 1;
  while (std::getline(in, row))
  {
    ++lineno;
    if (row.empty()) continue;
    const std::vector<std::string> f = split_tab(row);
    if (static_cast<int>(f.size()) <= std::max(iseq, ichg))
    {
      std::fprintf(stderr, "error: %s:%zu: too few columns\n", path, lineno);
      return 1;
    }
    char* end = nullptr;
    errno = 0;
    const long z = std::strtol(f[ichg].c_str(), &end, 10);
    if (errno != 0 || end == f[ichg].c_str() || *end != '\0')
    {
      std::fprintf(stderr, "error: %s:%zu: bad PrecursorCharge '%s'\n",
                   path, lineno, f[ichg].c_str());
      return 1;
    }
    const Precursor p{f[iseq], static_cast<int>(z)};
    if (!valid_precursor(p))
    {
      std::fprintf(stderr,
                   "error: %s:%zu: '%s' charge %ld unsupported "
                   "(need unmodified uppercase AAs, length >= 2, charge >= 1)\n",
                   path, lineno, f[iseq].c_str(), z);
      return 1;
    }
    out.push_back(p);
  }
  if (out.empty())
  {
    std::fprintf(stderr, "error: %s has no precursor rows\n", path);
    return 1;
  }
  return 0;
}

#define CHECK(cond)                                                                                \
  if (!(cond))                                                                                     \
  {                                                                                                \
    std::fprintf(stderr, "selftest FAILED at line %d: %s\n", __LINE__, #cond);                     \
    return 1;                                                                                      \
  }

// Stage 2 digest options; defaults match a standard tryptic search.
struct DigestOpts
{
  int missed = 1, min_len = 7, max_len = 30, min_charge = 2, max_charge = 3;
};

// FASTA -> tryptic peptides -> Precursors (one per peptide x charge), keeping
// the source protein identifier. Deduplicated on (sequence, charge); the first
// protein seen wins, so output is order-deterministic in FASTA order.
int digest_fasta(const char* path, const DigestOpts& d, std::vector<Precursor>& out)
{
  std::vector<OpenMS::FASTAFile::FASTAEntry> entries;
  try { OpenMS::FASTAFile().load(path, entries); }
  catch (const std::exception& e) { std::fprintf(stderr, "error: FASTA %s: %s\n", path, e.what()); return 1; }

  OpenMS::ProteaseDigestion dig;
  dig.setEnzyme("Trypsin");
  dig.setMissedCleavages(d.missed);           // specificity defaults to full-tryptic

  std::vector<std::pair<std::string, int>> seen; // insertion-ordered dedup key
  for (const auto& e : entries)
  {
    std::vector<OpenMS::AASequence> peps;
    dig.digest(OpenMS::AASequence::fromString(e.sequence), peps, d.min_len, d.max_len);
    for (const auto& pep : peps)
    {
      const std::string s = pep.toUnmodifiedString();
      if (!valid_precursor({s, d.min_charge})) continue; // canonical AAs only (drops X/B/Z/U)
      for (int z = d.min_charge; z <= d.max_charge; ++z)
      {
        const std::pair<std::string, int> key{s, z};
        // ponytail: O(n^2) dedup, fine for Stage-2 test proteomes; swap for an
        // unordered_set if this ever runs whole-proteome.
        bool dup = false;
        for (const auto& k : seen) if (k == key) { dup = true; break; }
        if (dup) continue;
        seen.push_back(key);
        out.push_back({s, z, e.identifier});
      }
    }
  }
  return 0;
}

int selftest(const std::string& model_dir, int threads)
{
  // ponytail: shape/sanity only. Parity vs AlphaPeptDeep Python is OpenMS's
  // PeptDeepInference_test, not this.
  const std::vector<Precursor> fx = {{"PEPTIDEK", 2}, {"SLYNTVATL", 3}, {"AGHCEWQMKYR", 2}};
  std::string tsv;
  try
  {
    tsv = generate(fx, model_dir, threads);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "selftest FAILED: %s\n", e.what());
    return 1;
  }

  std::istringstream in(tsv);
  std::string line;
  CHECK(std::getline(in, line));
  CHECK(line == kHeader);
  const std::vector<std::string> cols = split_tab(line);
  auto col_of = [&](const char* name) -> int {
    for (size_t i = 0; i < cols.size(); ++i)
      if (cols[i] == name) return static_cast<int>(i);
    return -1;
  };
  const int c_pmz = col_of("PrecursorMz"), c_mz = col_of("ProductMz");
  const int c_int = col_of("LibraryIntensity"), c_rt = col_of("NormalizedRetentionTime");
  const int c_seq = col_of("PeptideSequence");
  CHECK(c_pmz >= 0 && c_mz >= 0 && c_int >= 0 && c_rt >= 0 && c_seq >= 0);

  size_t rows = 0, seen[3] = {0, 0, 0};
  double max_int[3] = {0, 0, 0};
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    const std::vector<std::string> f = split_tab(line);
    CHECK(f.size() == cols.size());
    const double pmz = std::strtod(f[c_pmz].c_str(), nullptr);
    const double mz = std::strtod(f[c_mz].c_str(), nullptr);
    const double inten = std::strtod(f[c_int].c_str(), nullptr);
    const double rt = std::strtod(f[c_rt].c_str(), nullptr);
    CHECK(std::isfinite(pmz) && pmz > 0.0);
    CHECK(std::isfinite(mz) && mz > 0.0);
    CHECK(std::isfinite(inten) && inten > 0.0 && inten <= 1.0 + 1e-5);
    CHECK(std::isfinite(rt));
    size_t w = 3;
    if (f[c_seq] == "PEPTIDEK") w = 0;
    else if (f[c_seq] == "SLYNTVATL") w = 1;
    else if (f[c_seq] == "AGHCEWQMKYR") w = 2;
    CHECK(w < 3);
    ++seen[w];
    if (inten > max_int[w]) max_int[w] = inten;
    ++rows;
  }
  CHECK(rows >= 12); // three peptides of length 8/9/11, several b/y rows each
  for (size_t w = 0; w < 3; ++w)
  {
    CHECK(seen[w] > 0);
    // predictMS2 base-peak-normalizes: each peptide's max must be exactly 1.
    CHECK(std::fabs(max_int[w] - 1.0) < 1e-4);
  }

  // Byte-identical output for identical input is a hard requirement.
  std::string again;
  try
  {
    again = generate(fx, model_dir, threads);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "selftest FAILED: %s\n", e.what());
    return 1;
  }
  CHECK(tsv == again);

  // Stage 2: FASTA digest branch. A temp FASTA exercises the real FASTAFile +
  // ProteaseDigestion path. Sequence is proline-free so the K/R-before-P
  // exception doesn't complicate the expectation. ELVISK|LIVESR|NOWK.
  {
    const char* fasta = ">TESTPROT test\nELVISKLIVESRNOWK\n";
    const std::string ftmp = "/tmp/odialibgen_selftest.fasta";
    auto write_fasta = [&] { std::ofstream o(ftmp); o << fasta; };
    DigestOpts d;
    d.missed = 0; d.min_len = 4; d.max_len = 30; d.min_charge = 2; d.max_charge = 2;
    write_fasta();
    std::vector<Precursor> dp;
    CHECK(digest_fasta(ftmp.c_str(), d, dp) == 0);
    // Fully-tryptic peptides, each at charge 2, protein set. Assertions are
    // chemically robust (exact count depends on OpenMS terminus handling):
    // within length bounds, internal peptides end in K/R.
    CHECK(dp.size() >= 2);
    size_t kr_terminal = 0;
    for (const auto& p : dp)
    {
      CHECK(p.charge == 2);
      CHECK(p.protein == "TESTPROT");
      CHECK(p.sequence.size() >= 4 && p.sequence.size() <= 30);
      const char last = p.sequence.back();
      if (last == 'K' || last == 'R') ++kr_terminal;
    }
    CHECK(kr_terminal >= 2); // at least the two internal peptides end in K/R
    // Determinism: identical FASTA -> identical precursor list.
    write_fasta();
    std::vector<Precursor> dp2;
    CHECK(digest_fasta(ftmp.c_str(), d, dp2) == 0);
    std::remove(ftmp.c_str());
    CHECK(dp2.size() == dp.size());
    for (size_t k = 0; k < dp.size(); ++k)
      CHECK(dp[k].sequence == dp2[k].sequence && dp[k].charge == dp2[k].charge);
  }

  std::puts("selftest OK");
  return 0;
}
#undef CHECK

} // namespace

// Load the raw target TSV into a TargetedExperiment, run the OpenSwathAssayGenerator
// refinement (reannotate -> restrict -> select detecting), optionally append
// OpenSwathDecoyGenerator shuffle decoys, and write the final TSV. Parameter
// values are copied verbatim from the two TOPP tools so output matches them;
// the equivalence gate checks this empirically. Uses temp files because the
// TransitionTSVFile / MRM* APIs are file- and TargetedExperiment-based.
//
// NOTE (vendored patch P3): MRMDecoy shuffle seeds with time(nullptr) in a
// release build, so the DECOY rows are NOT reproducible -- identical to the
// standalone OpenSwathDecoyGenerator. The target library is deterministic.
int refine_and_decoy(const std::string& raw_tsv, const std::string& out_path, bool decoys)
{
  using namespace OpenMS;
  const std::string tmp_raw = out_path + ".raw.tmp.tsv";
  { std::ofstream o(tmp_raw, std::ios::binary); if (!o) { std::fprintf(stderr, "error: temp %s\n", tmp_raw.c_str()); return 1; } o << raw_tsv; }

  const std::vector<std::string> frag_types{"b", "y"};
  const std::vector<size_t> frag_charges{1, 2, 3, 4};        // AssayGenerator default
  const std::vector<std::pair<double, double>> no_swathes;   // no SWATH restriction
  try
  {
    TargetedExperiment exp;
    TransitionTSVFile().convertTSVToTargetedExperiment(tmp_raw.c_str(), FileTypes::TSV, exp);

    // Values copied verbatim from OpenSwathAssayGenerator's registered defaults
    // (product m/z limits 350-2000 Th are the ones that actually change the
    // transition set; using 0/inf kept small fragments the tool drops).
    MRMAssay assay;
    assay.reannotateTransitions(exp, 0.025, 0.025, frag_types, frag_charges, false, false);
    assay.restrictTransitions(exp, 350.0, 2000.0, no_swathes);
    assay.detectingTransitions(exp, 6, 6);                   // min=max=6

    if (decoys)
    {
      TargetedExperiment dec;
      MRMDecoy().generateDecoys(exp, dec, "shuffle", 1.0, /*switchKR*/ true, "DECOY_",
                                /*max_attempts*/ 10, /*identity_threshold*/ 0.7,
                                /*precursor_mz_shift*/ 0.0, /*product_mz_shift*/ 20.0,
                                /*product_mz_threshold*/ 0.025, frag_types, frag_charges,
                                false, false);
      exp += dec; // TargetedExperiment::operator+= merges proteins/peptides/transitions
    }

    TransitionTSVFile().convertTargetedExperimentToTSV(out_path.c_str(), exp);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: refine/decoy: %s\n", e.what());
    std::remove(tmp_raw.c_str());
    return 1;
  }
  std::remove(tmp_raw.c_str());
  return 0;
}

int main(int argc, char** argv)
{
  std::string in_path, fasta_path, out_path, model_dir;
  int threads = 1; // explicit, never omp_get_max_threads(); 1 = deterministic default
  bool do_selftest = false, raw = false, decoys = true;
  DigestOpts dig;

  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "--selftest") do_selftest = true;
    else if (a == "-in" && i + 1 < argc) in_path = argv[++i];
    else if (a == "-fasta" && i + 1 < argc) fasta_path = argv[++i];
    else if (a == "-out" && i + 1 < argc) out_path = argv[++i];
    else if (a == "-model_dir" && i + 1 < argc) model_dir = argv[++i];
    else if (a == "-threads" && i + 1 < argc) threads = std::atoi(argv[++i]);
    else if (a == "-missed_cleavages" && i + 1 < argc) dig.missed = std::atoi(argv[++i]);
    else if (a == "-min_len" && i + 1 < argc) dig.min_len = std::atoi(argv[++i]);
    else if (a == "-max_len" && i + 1 < argc) dig.max_len = std::atoi(argv[++i]);
    else if (a == "-min_charge" && i + 1 < argc) dig.min_charge = std::atoi(argv[++i]);
    else if (a == "-max_charge" && i + 1 < argc) dig.max_charge = std::atoi(argv[++i]);
    else if (a == "-no_decoys") decoys = false;
    else if (a == "-raw") raw = true; // emit pre-refinement, pre-decoy transitions
    else
    {
      std::fprintf(stderr, "error: unknown or incomplete argument '%s'\n", a.c_str());
      return 2;
    }
  }

  const bool have_input = !in_path.empty() ^ !fasta_path.empty(); // exactly one
  if (threads < 1 || model_dir.empty() || (!do_selftest && (!have_input || out_path.empty())))
  {
    std::fprintf(stderr,
                 "usage (Stage 1): OpenDIALibGen -in peptides.tsv -out lib.tsv -model_dir <dir> [-threads N] [-raw]\n"
                 "usage (Stage 2): OpenDIALibGen -fasta db.fasta -out lib.tsv -model_dir <dir>\n"
                 "                   [-missed_cleavages 1] [-min_len 7] [-max_len 30]\n"
                 "                   [-min_charge 2] [-max_charge 3] [-no_decoys] [-raw] [-threads N]\n"
                 "         selftest: OpenDIALibGen --selftest -model_dir <dir> [-threads N]\n"
                 "model_dir must contain peptdeep_rt_dynamic.onnx and peptdeep_ms2_dynamic.onnx\n");
    return 2;
  }

  if (do_selftest) return selftest(model_dir, threads);

  std::vector<Precursor> precursors;
  if (!fasta_path.empty())
  {
    if (digest_fasta(fasta_path.c_str(), dig, precursors) != 0) return 1;
  }
  else if (read_precursors(in_path.c_str(), precursors) != 0) return 1;

  std::string tsv;
  try { tsv = generate(precursors, model_dir, threads); }
  catch (const std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 1; }

  if (raw || fasta_path.empty())
  {
    // Stage 1 semantics, or -raw: write the raw target transitions as-is.
    std::ofstream out(out_path, std::ios::binary);
    if (!out) { std::fprintf(stderr, "error: cannot write %s\n", out_path.c_str()); return 1; }
    out << tsv;
    out.close();
    std::printf("%zu precursors -> %s\n", precursors.size(), out_path.c_str());
    std::puts(raw ? "raw target transitions (pre-refine, pre-decoy)"
                  : "targets only - generate decoys with OpenSwathDecoyGenerator");
    return 0;
  }

  // Stage 2: fold in assay refinement + decoys in-process.
  if (refine_and_decoy(tsv, out_path, decoys) != 0) return 1;
  std::printf("%zu precursors digested+refined%s -> %s\n", precursors.size(),
              decoys ? "+decoyed" : "", out_path.c_str());
  return 0;
}

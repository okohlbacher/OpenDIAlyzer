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
                    "%.6f\t%.6f\t%d\t%d\t%.6f\t%.4f\t%s\t%s\tUNKNOWN\t%c\t%d\t%s\t%s\t0\t1\t0\t1\n",
                    prec_mz, r.mz, z, r.fz, static_cast<double>(r.intensity),
                    static_cast<double>(rts[i]), seq.c_str(), seq.c_str(),
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

  std::puts("selftest OK");
  return 0;
}
#undef CHECK

} // namespace

int main(int argc, char** argv)
{
  std::string in_path, out_path, model_dir;
  int threads = 1; // explicit, never omp_get_max_threads(); 1 = deterministic default
  bool do_selftest = false;

  for (int i = 1; i < argc; ++i)
  {
    const std::string a = argv[i];
    if (a == "--selftest") do_selftest = true;
    else if (a == "-in" && i + 1 < argc) in_path = argv[++i];
    else if (a == "-out" && i + 1 < argc) out_path = argv[++i];
    else if (a == "-model_dir" && i + 1 < argc) model_dir = argv[++i];
    else if (a == "-threads" && i + 1 < argc) threads = std::atoi(argv[++i]);
    else
    {
      std::fprintf(stderr, "error: unknown or incomplete argument '%s'\n", a.c_str());
      return 2;
    }
  }

  if (threads < 1 || model_dir.empty() || (!do_selftest && (in_path.empty() || out_path.empty())))
  {
    std::fprintf(stderr,
                 "usage: OpenDIALibGen -in peptides.tsv -out library.tsv -model_dir <dir> [-threads N]\n"
                 "       OpenDIALibGen --selftest -model_dir <dir> [-threads N]\n"
                 "model_dir must contain peptdeep_rt_dynamic.onnx and peptdeep_ms2_dynamic.onnx\n"
                 "(OpenMS build tree: share/OpenMS/models/)\n");
    return 2;
  }

  if (do_selftest) return selftest(model_dir, threads);

  std::vector<Precursor> precursors;
  if (read_precursors(in_path.c_str(), precursors) != 0) return 1;

  std::string tsv;
  try
  {
    tsv = generate(precursors, model_dir, threads);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }

  std::ofstream out(out_path, std::ios::binary);
  if (!out) { std::fprintf(stderr, "error: cannot write %s\n", out_path.c_str()); return 1; }
  out << tsv;
  out.close();

  std::printf("%zu precursors -> %s\n", precursors.size(), out_path.c_str());
  std::puts("targets only - generate decoys with OpenSwathDecoyGenerator");
  return 0;
}

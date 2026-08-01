// Self-check for odia_seqstore.h. Standalone: no OpenMS, no framework.
//
// The properties that matter are the ones whose violation is SILENT: a dangling view after growth
// would read freed memory and produce a plausible-looking wrong sequence, and a broken dedup would
// quietly reintroduce the duplication the store exists to remove.

#include "odia_seqstore.h"

#include <cassert>
#include <cstdio>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

using odia::SequenceStore;

static void t1_roundtrip()
{
  SequenceStore s;
  const std::vector<std::string> in = {"", "A", "PEPTIDEK_2", "DECOY_PEPTIDEK_2",
                                       std::string(1000, 'x'), "y7^1"};
  std::vector<SequenceStore::Id> ids;
  for (const auto& v : in) { ids.push_back(s.append(v)); }
  for (std::size_t i = 0; i < in.size(); ++i) { assert(s.view(ids[i]) == in[i]); }
  std::printf("  T1 round-trip incl. empty and 1000-char string: OK\n");
}

static void t2_dedup()
{
  SequenceStore s;
  const auto a = s.intern("sp|P02768|ALBU_HUMAN");
  const auto b = s.intern("sp|P02768|ALBU_HUMAN");
  const auto c = s.intern("sp|P01023|A2MG_HUMAN");
  assert(a == b);
  assert(a != c);
  assert(s.size() == 2);                       // stored twice would defeat the whole point
  assert(s.find("sp|P02768|ALBU_HUMAN") == a);
  assert(s.find("not present") == SequenceStore::npos);
  std::printf("  T2 intern dedups, find works, absent -> npos: OK\n");
}

/// THE hazard this design exists to avoid. A single growing buffer would reallocate and invalidate
/// every earlier view; chunks must not. Take views BEFORE appending a great deal more, then verify
/// them afterwards.
static void t3_views_survive_growth()
{
  SequenceStore s(4096);                       // tiny chunks: force many of them
  std::vector<SequenceStore::Id> ids;
  std::vector<std::string> expect;
  for (int i = 0; i < 50; ++i)
  {
    expect.push_back("early_" + std::to_string(i) + std::string(60, 'a'));
    ids.push_back(s.append(expect.back()));
  }
  std::vector<std::string_view> early;
  for (auto id : ids) { early.push_back(s.view(id)); }

  for (int i = 0; i < 20000; ++i) { s.append("filler_" + std::to_string(i) + std::string(80, 'b')); }

  for (std::size_t i = 0; i < early.size(); ++i)
  {
    assert(early[i] == expect[i]);             // would read freed memory if chunks reallocated
    assert(s.view(ids[i]) == expect[i]);
  }
  std::printf("  T3 views taken before 20k appends still valid (%zu chunks): OK\n",
              std::size_t(s.bytes() / 4096));
}

static void t4_oversized_string()
{
  SequenceStore s(1024);
  const std::string big(5000, 'z');            // larger than a whole chunk
  const auto id = s.append(big);
  const auto after = s.append("small");
  assert(s.view(id) == big);                   // must not be split or truncated
  assert(s.view(after) == "small");
  std::printf("  T4 string larger than a chunk stored whole: OK\n");
}

static void t5_freeze_index()
{
  SequenceStore s;
  const auto a = s.intern("KEEP_ME");
  s.intern("AND_ME");
  const std::size_t before = s.bytes();
  s.freezeIndex();
  assert(s.frozen());
  assert(s.view(a) == "KEEP_ME");              // reads must survive dropping the dedup map
  assert(s.bytes() < before);
  std::printf("  T5 freezeIndex frees the map, reads still work (%zu -> %zu B): OK\n",
              before, s.bytes());
}

/// The measurement that justifies the class: peptide_ref on the real library is 78.6M references
/// to ~7.1M distinct values. Model that ratio and report both costs.
static void t6_duplication_saving()
{
  SequenceStore s;
  const std::size_t distinct = 200000;
  const std::size_t refs = distinct * 11;      // the measured 78.6M / 7.1M ratio
  std::mt19937 rng(1234);
  std::vector<SequenceStore::Id> held;
  held.reserve(refs);
  for (std::size_t i = 0; i < distinct; ++i)
  {
    s.intern("DECOY_PEPTIDESEQ" + std::to_string(i) + "_2");
  }
  for (std::size_t i = 0; i < refs; ++i)
  {
    held.push_back(SequenceStore::Id(rng() % distinct));
  }
  const std::size_t store_b = s.bytes() + held.capacity() * SequenceStore::bytesPerRef();
  const std::size_t str_b   = s.bytesAsStdString(refs);
  std::printf("  T6 %zu refs to %zu distinct: std::string %.1f MB vs store %.1f MB (%.1fx)\n",
              refs, distinct, str_b / 1048576.0, store_b / 1048576.0,
              double(str_b) / double(store_b));
  assert(store_b < str_b);
}

int main()
{
  std::printf("odia_seqstore_test\n");
  t1_roundtrip();
  t2_dedup();
  t3_views_survive_growth();
  t4_oversized_string();
  t5_freeze_index();
  t6_duplication_saving();
  std::printf("odia_seqstore_test OK\n");
  return 0;
}

// odia_fdr.h — error control BEYOND the precursor-level target-decoy q-value.
//
// Three things the scoring/FDR backlog (docs/OpenDIAlyzer-scoring-fdr-backlog.md, Area 2)
// lists as "prototype-next", all algorithmic, no new dependencies:
//
//   1. CONTEXT FDR — peptide- and protein-level q-values by rolling precursor scores up to
//      the entity and re-running target-decoy on that set (Rosenberger 2017, PyProphet).
//      Controlling FDR at 1% on precursors does NOT give 1% on peptides or proteins: a
//      protein collects the best of many independent precursor draws, so its score is an
//      extreme-value statistic and its error rate is strictly worse than the precursor one.
//      Reporting a precursor q-value as if it were a protein q-value is the single most
//      common way a DIA result overstates its confidence.
//
//   2. PICKED protein-group FDR (Savitski 2015; The 2022) — instead of ranking all target
//      and all decoy proteins together, let each target compete with ITS OWN decoy partner
//      and enter the ranking once, as the winner. Plain target-decoy at the protein level is
//      biased by protein size: a large protein offers more peptides for a random decoy hit
//      to land on, so decoy proteins are systematically over-represented among large entries
//      and the estimate is inflated in a size-dependent way. Competition within the pair
//      removes that, because the pair members have the same size by construction.
//
//   3. ENTRAPMENT FDP (Wen/Freestone/Keich/Noble 2025) — a validation rig, not a model.
//      Target-decoy FDR in DIA is repeatedly found OPTIMISTIC; the way to know whether a
//      reported 1% is real is to search a database salted with sequences that cannot be
//      present and count how many come back. See entrapmentFdp() for exactly which estimator
//      this is and what it does not do.
//
// The q-value / PEP mathematics is NOT reimplemented here — it is the same step-down
// target-decoy estimator the precursor level uses (lda_detail::assignQValues), applied to a
// different set of items. One implementation, three levels.

#ifndef ODIA_FDR_H
#define ODIA_FDR_H

#include "odia_lda.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace odia
{
namespace fdr
{

/// One item at whatever level is being controlled: a peptide, or a protein (group).
struct Entity
{
  std::string id;            ///< modified sequence, or protein accession / group key
  int label = 1;             ///< 1 = target, 0 = decoy
  double score = 0.0;        ///< best member score (see rollUp)
  double qvalue = 1.0;
  double pvalue = 1.0;
  double pep = 1.0;
};

/// Roll member scores up to entities: each entity takes the score of its BEST member.
///
/// `member_key[i]` names the entity row i belongs to, `member_score[i]` is its score and
/// `member_label[i]` its class. Rows whose key is empty are skipped (unmapped precursors).
/// The maximum is the right statistic here and is what makes §1 above true: it is exactly
/// the "best of n draws" that inflates the entity-level error rate relative to the member
/// one, which is why the q-values must be recomputed at this level rather than inherited.
inline std::vector<Entity> rollUp(const std::vector<std::string>& member_key,
                                  const std::vector<double>& member_score,
                                  const std::vector<int>& member_label)
{
  std::vector<Entity> out;
  const std::size_t n = member_key.size();
  if (member_score.size() != n || member_label.size() != n) { return out; }

  std::unordered_map<std::string, std::size_t> index;
  index.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    if (member_key[i].empty()) { continue; }
    auto inserted = index.emplace(member_key[i], out.size());
    if (inserted.second)
    {
      Entity e;
      e.id = member_key[i];
      e.label = member_label[i] == 1 ? 1 : 0;
      e.score = member_score[i];
      out.push_back(std::move(e));
    }
    else if (member_score[i] > out[inserted.first->second].score)
    {
      out[inserted.first->second].score = member_score[i];
    }
  }
  return out;
}

/// Target-decoy q-value / p-value / PEP over an entity set, identical estimator to the
/// precursor level. Entities are updated in place; their order is preserved.
inline void assignQValues(std::vector<Entity>& entities, bool use_pi0)
{
  std::vector<lda_detail::RankedGroup> ranked;
  ranked.reserve(entities.size());
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    lda_detail::RankedGroup r;
    r.group_index = i;
    r.best_row = i;
    r.label = entities[i].label;
    r.score = entities[i].score;
    r.qvalue = 1.0;
    ranked.push_back(r);
  }
  lda_detail::assignQValues(ranked, use_pi0);
  for (const auto& r : ranked)
  {
    entities[r.group_index].qvalue = r.qvalue;
    entities[r.group_index].pvalue = r.pvalue;
    entities[r.group_index].pep = r.pep;
  }
}

/// True if `id` is `decoy_tag` + something; writes the partner (target) id to `target_id`.
inline bool decoyPartnerId(const std::string& id, const std::string& decoy_tag,
                           std::string& target_id)
{
  if (decoy_tag.empty() || id.size() <= decoy_tag.size()) { return false; }
  if (id.compare(0, decoy_tag.size(), decoy_tag) != 0) { return false; }
  target_id = id.substr(decoy_tag.size());
  return true;
}

/// Picked target-decoy competition (Savitski 2015; The 2022).
///
/// Each target/decoy pair is collapsed to ONE entry carrying the better of the two scores
/// and the label of whichever won; entities with no partner pass through unchanged. The
/// caller then runs assignQValues() on the result.
///
/// Why this is not merely a refinement: a decoy protein's score is the max over its decoy
/// peptides, so it grows with the protein's peptide count. Ranking all targets against all
/// decoys therefore mixes a size-dependent null into one list, and the resulting estimate is
/// wrong in opposite directions for large and small proteins. A target and its own decoy
/// have the same size by construction, so competing them cancels that dependence exactly.
///
/// Assumes the decoy id is `decoy_tag` + the target id -- the MRMDecoy convention (see
/// MRMDecoy.cpp:842 for compounds, :507 for proteins). Entities whose ids do not follow it
/// simply never pair, so the function degrades to plain target-decoy rather than misbehaving.
inline std::vector<Entity> pickedCompetition(const std::vector<Entity>& entities,
                                             const std::string& decoy_tag,
                                             std::size_t* n_paired = nullptr)
{
  std::unordered_map<std::string, std::size_t> target_index;
  target_index.reserve(entities.size());
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    if (entities[i].label == 1) { target_index.emplace(entities[i].id, i); }
  }

  std::vector<Entity> out;
  out.reserve(entities.size());
  std::vector<char> consumed(entities.size(), 0);
  std::size_t paired = 0;

  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    if (entities[i].label == 1) { continue; }
    std::string partner;
    if (!decoyPartnerId(entities[i].id, decoy_tag, partner)) { continue; }
    const auto it = target_index.find(partner);
    if (it == target_index.end()) { continue; }

    const Entity& decoy = entities[i];
    const Entity& target = entities[it->second];
    // Ties go to the DECOY. A tie carries no evidence either way, and resolving it in the
    // target's favour would silently bias the estimate downward -- the exact direction the
    // 2025 entrapment critiques find DIA FDR already errs in.
    out.push_back(target.score > decoy.score ? target : decoy);
    consumed[i] = 1;
    consumed[it->second] = 1;
    ++paired;
  }
  for (std::size_t i = 0; i < entities.size(); ++i)
  {
    if (!consumed[i]) { out.push_back(entities[i]); }
  }
  if (n_paired) { *n_paired = paired; }
  return out;
}

/// Result of the entrapment check. `valid` is false when the inputs cannot support an
/// estimate (no entrapment sequences in the database, or nothing reported).
struct EntrapmentEstimate
{
  std::size_t n_reported = 0;      ///< discoveries at the threshold (targets + entrapments)
  std::size_t n_entrapment = 0;    ///< of those, how many were entrapment sequences
  std::size_t db_target = 0;       ///< real target sequences searched
  std::size_t db_entrapment = 0;   ///< entrapment sequences searched
  double ratio = 0.0;              ///< db_entrapment / db_target
  double fdp = 1.0;                ///< estimated false discovery PROPORTION
  bool valid = false;
};

/// Combined entrapment estimator (Wen/Freestone/Keich/Noble, Nat Methods 2025).
///
///     FDP = (1 + 1/r) * n_entrapment / n_reported ,   r = db_entrapment / db_target
///
/// Reading: entrapment sequences cannot be in the sample, so every entrapment discovery is
/// a false one. They sample the false-discovery process at rate r relative to the real
/// targets, so the n_entrapment false discoveries seen imply about n_entrapment/r further
/// false discoveries hiding among the reported TARGETS. Total false = n_entrapment(1 + 1/r).
///
/// This is a measurement of the reported FDR, not a replacement for it: run it on a search
/// against an entrapment-salted library and compare `fdp` to the nominal q-value cutoff the
/// discoveries were taken at. They should agree; the 2025 literature finds they frequently
/// do not, and that DIA tools tend to be optimistic.
///
/// NOT the paired estimator. The paired variant of that paper is more powerful but needs an
/// explicit target<->entrapment pairing in the library and a different formula; implementing
/// it from a half-remembered description would be worse than not having it. Combined is the
/// conservative one and is well defined from the counts alone.
inline EntrapmentEstimate entrapmentFdp(std::size_t n_reported, std::size_t n_entrapment,
                                        std::size_t db_target, std::size_t db_entrapment)
{
  EntrapmentEstimate e;
  e.n_reported = n_reported;
  e.n_entrapment = n_entrapment;
  e.db_target = db_target;
  e.db_entrapment = db_entrapment;
  if (db_target == 0 || db_entrapment == 0 || n_reported == 0) { return e; }
  e.ratio = static_cast<double>(db_entrapment) / static_cast<double>(db_target);
  e.fdp = std::min(1.0, (1.0 + 1.0 / e.ratio) * static_cast<double>(n_entrapment) /
                          static_cast<double>(n_reported));
  e.valid = true;
  return e;
}

/// Targets at or below `q` (decoys never count as identifications).
inline std::size_t countAtQ(const std::vector<Entity>& entities, double q)
{
  return static_cast<std::size_t>(
    std::count_if(entities.begin(), entities.end(),
                  [&](const Entity& e) { return e.label == 1 && e.qvalue <= q; }));
}

} // namespace fdr
} // namespace odia

#endif // ODIA_FDR_H

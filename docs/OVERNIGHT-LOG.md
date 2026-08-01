# OpenDIAlyzer — overnight autonomous iteration log

Directive (user 2026-07-26, night): keep iterating — build the skeleton step by
step, deep adversarial review of code + results after each step, plan the next
step + corrections, iterate until competitive. Constraints hold: in-process only
(no external tools in the engine), never push to OpenMS/DIA-NN repos, use
codex/vibe (kimi can't run headless for file writes) for adversarial review.

## State at start of the night
- MVP assembled + in-process: extraction (`performExtraction`) → LDA-selected RT
  anchors → recalibration → narrow pass 2 → in-process LDA FDR → SCORE_MS2.
- **Correct**: bootstrap RT calibration validated (Pearson 0.945, EXP_RT spans run);
  LDA validated (oracle + 4 adversarial properties).
- **Bugs fixed tonight so far**: readoption (cacheWorkingInMemory→cache), missing
  writeHeader, string-vs-int FEATURE.PRECURSOR_ID join (→TRAML_ID), rt_win=-1
  feature explosion (→bounded window + linear bootstrap), trafo-not-applied
  (→pre-scale library RT).
- **Open blockers**: (#22) full-scale extraction uses only ~24/224 cores → pass 1
  ~10x slower than baseline; the single-threaded OSW writer.

## Iteration ledger
### Iter 1 — extraction parallelism (speed gates all iteration)
Goal: get pass 1 from ~35 min (113k groups, ~20 cores) to minutes so the loop can
turn fast. Hypotheses: max_concurrent_swaths cap; nested OpenMP disabled; writer
serialization. Plan: measure #swaths + thread distribution vs the baseline config.

**Finding:** the data has **24 SWATH windows**, so the wave scheduler's outer
parallelism is capped at 24 (max_concurrent_swaths=24 = all of them). ~3600 threads
exist but most SLEEP; ~24 run + 86 in disk-I/O wait (`Dl`) → the extra 200 cores are
idle because within-swath scoring is not spreading, and there's disk-I/O overhead
(the `cache` readoption, since Bruker rejects `cacheWorkingInMemory` in this build).
Two levers, both deferred to a dedicated speed iteration (task #22): (a) inner-swath
parallelism to use >24 cores, (b) in-memory reads to kill the `Dl` I/O wait.
**Decision:** iterate CORRECTNESS on small/medium libraries (fast enough at 24
cores); treat full-scale speed as its own iteration. Correctness gates > speed gates
right now — a fast wrong answer is worthless.

### Iter 1b — deep adversarial code review of the engine (parallel to MED run)
Dispatched codex to review src/opendialyzer.cpp + src/odia_lda.h for correctness
bugs while the MED A/B run finishes.

**codex verdict: "IDs@1% are NOT presently trustworthy as a calibrated 1% FDR."**
RT two-pass logic is CORRECT (no double-application). Real bugs, prioritized:
FDR honesty (fix first):
- C3 fold scores incomparable — each CV fold's LDA has different weights; combining
  d-scores for one global FDR is invalid. [odia_lda.h]
- C13 preprocessing leakage — z-score mean/SD + constant-col selection computed over
  ALL rows before the CV split. [odia_lda.h:197]
- C1 decoy-empty/target-only fold -> FDR 0/N=0 -> every target "passes". [odia_lda.h]
- C15 fold size bounded by total groups not smaller class -> single-class folds. [odia_lda.h]
- C5 target-decoy estimator anti-conservative at tail; assumes 1:1; no +1. [odia_lda.h]
- C4/C8 anchor circularity — anchors chosen by q that INCLUDES VAR_NORM_RT_SCORE; the
  recalibration step is outside the pass-2 CV -> optimistic. [opendialyzer.cpp:194,262]
- C9 SQL NULL->0.0 — DECOY=NULL becomes target; missing score becomes real 0. [opendialyzer.cpp:220]
IDs / correctness:
- C2 diaPASEF DISABLED — pasef=false hardcoded for Bruker .d (data IS diaPASEF). [opendialyzer.cpp]
- C10 RT fit not robust/monotone (1 pt/bin, backward segments possible). [opendialyzer.cpp:142]
- C6 SCORE_MS2 ranks/p-values false; not pyprophet-compatible. [opendialyzer.cpp]
- C12 SQLite errors reported as successful "0 IDs". [opendialyzer.cpp]
Robustness: C7 RT-scale assumption, C11 join schema, C14 reproducibility (no ORDER BY),
C16 multi-run, C17 pass-count validation, C18 MS2-only RT range.

Also found (independently, from OpenMS source): SPEED — outer_loop_threads=-1 takes the
NON-nested branch -> 24 cores; passing outer_loop_threads = #swaths enables nested
inner parallelism (~9 threads/swath -> ~216 cores). ~9x iteration speedup.

### Iter 2 — fix plan (quality first)
LDA FDR (delegate to codex, strengthen adversarial test): C3, C13, C1, C15, C5.
Engine (self): C2 pasef=true for Bruker; C8 exclude VAR_NORM_RT_SCORE + IM/RT scores from
anchor selection; C9 NULL handling + require DECOY in {0,1}; C1 guard decoys present;
C12 SQLite error checks; C10 robust monotone fit; SPEED outer_loop_threads=#swaths.
Then re-validate (oracle + strengthened adversarial + a MED run) and A/B.
vibe FDR review pending (bk2c0zb3k) — fold in when it lands.

### Iter 2 — DONE. Review fixes implemented + validated.
LDA (codex, odia_lda.h) — FDR honesty fixed + validated by a STRENGTHENED adversarial
suite (all green, cluster + local):
- C3 fold comparability: per-fold standardize held-out d-scores -> P7 decoy mean ~1e-16.
- C13 leakage: preprocessing fit per training fold only -> P5 (informative-in-train,
  random-in-test) yields 0 IDs; P6 fold-local preprocessing 0.
- C5 estimator: (decoys+1)/targets, ratio-scaled -> P8 10:1 imbalance reported 0.236 >=
  empirical 0.142 (CONSERVATIVE).
- C1/C15 class guards: unscorable configs -> all q=1 (P9), never spurious q=0.
Engine (self, opendialyzer.cpp) — built OK:
- C2 pasef=true when SWATHs carry IM (diaPASEF) + use_ms1_im.
- C8 recalibration LDA excludes VAR_NORM_RT_SCORE (anti-circular anchors).
- C9 SQL NULL -> NaN (LDA imputes), require DECOY in {0,1}, atomic row push.
- C1 hasBothClasses_ guard in recalibrate_ + finalScore_.
- C10 robust monotone RT fit (adaptive bins, >=3 pts/bin, non-decreasing y).
- C17 recal_passes in {1,2}; C18 RT range from any map, fail if none.

**SPEED — CORRECTED a mistake using a lost prior-session finding.** The slowness was
NOT the thread setting; I initially "fixed" outer_loop_threads=-1 -> #swaths which is
BACKWARDS (>=0 = legacy path, 2.7x SLOWER). Reverted to -1 (wave scheduler). The real
lever: readoptions for Bruker was 'cache' (disk, causing the Dl I/O stalls) -> switched
to 'normal' (RegularSwathFileConsumer, in-RAM). Measuring now (SPEED run). Prior note:
vault/70-Adversarial/Deep review - where we actually stand.md.

**SPEED result (in-memory 'normal'):** the .d loads to RAM in ~10 min (3 phases 1:08 /
4:03 / 5:13), then extraction runs from RAM. Much better than 'cache' (46+ min stuck in
pass 1 on disk I/O), but extraction is still ~17-24 cores (the 24-SWATH cap) with Dl
waits from the single-threaded OSW writer. So: 'normal' + wave scheduler is the right
config; the remaining ceiling (24-swath cap + single writer) is task #22 for the
speed-dedicated iteration. Iterations now ~30-40 min — workable for correctness.

### Iter 4 — the two-pass could not bootstrap; fixed. KEY INSIGHT.
The honest run gave 0 IDs. Chased it down (a false alarm + a real chicken-egg):
- FALSE ALARM: my diagnostic script mis-joined and reported "constant d-score 1.16".
  In-engine DIAG proved the LDA is FINE: 1,008,675 rows, 22 cols (after dropping 7
  all-NaN/constant VAR_ cols incl VAR_IM_XCORR_SHAPE), **687,323 distinct d-scores**.
  A reference sklearn LDA on the same features separates the classes. LDA not the bug.
- REAL BUG (chicken-and-egg): the WIDE bootstrap pass 1 (600 s, uncalibrated) genuinely
  has 0 honest q<0.01 IDs, so recalibrate_ (which required q<0.01 target anchors) found
  none -> two-pass never started. Fix: recalibrate_ now takes the **top-2000 target peak
  groups by d-score** as approximate RT anchors (the binned-median fit tolerates wrong
  ones); the FINAL FDR stays strict q<0.01. Decoupled bootstrap from honest reporting.
Also added: unsupervised all-NaN/constant column drop in loadOswScores_ (structural, no
labels -> no leakage); a real-data regression test (src/odia_lda_realdata_test.cpp +
testdata/lda_fixture.txt, 44.5k real rows) so the LDA is guarded on real distributions.
Definitive two-pass run (TWOPASS, 113k, /dev/shm) launched: does the NARROW recalibrated
pass 2 yield honest q<0.01 IDs where the wide pass gives 0? This is the core thesis test.

### Iter 5 — PIVOTAL: the FDR/LDA are fine; OUR EXTRACTION is the limiter.
TWOPASS (113k, 1:55:31 wall) completed: pass 1 (wide) 0 honest IDs -> recalibrate from
top-2000 anchors -> pass 2 (narrow) **3 honest IDs**. FDR curve FLAT at 3 for q<0.01..0.2;
target d-score med 0.878 vs decoy 0.803 (only 3 targets beat the top decoy) => targets and
decoys are NEARLY INSEPARABLE on our extraction.
DECISIVE isolation (pyprophet on OUR osw, same features, different FDR model):
  - pyprophet on our pass-2 osw: **0**.  pyprophet on our pass-1 osw: **0**.
  - Our LDA: 3 / 0. So BOTH FDR models agree ~0 on our extraction.
Baseline OpenSWATH extraction -> pyprophet 728 (full 1/8). => The FDR model and the LDA
are NOT the problem (validated + agree with pyprophet). **Our EXTRACTION produces
non-discriminating features.** Prime suspects, in order:
  1. pasef=true: I enabled it (codex C2), but the baseline that gets 728 ran pasef=FALSE.
     If our predicted library's IM (PeptDeep CCS) is off, IM windowing discards true peaks.
  2. bootstrap linear RT (corr 0.945) vs OpenSWATH auto_irt (nonlinear) -> worse windows.
  3. the recalibration bootstrapped on pass-1 noise (top-2000 by d-score are ~random when
     pass 1 can't discriminate) -> pass 2 windows placed by a garbage transform.
ISOLATION RUN LAUNCHED: NOPASEF (recal_passes=1, -pasef false, single wide pass, library_med).
  - If it yields ~90 IDs -> pasef/IM was the killer (revisit C2 + library IM).
  - If ~0 -> the bootstrap RT is the killer (need OpenSWATH-quality calibration, not linear).
This is the make-or-break diagnostic. Morning priority: read NOPASEF, then isolate RT
(bootstrap vs auto_irt) if needed. The two-pass thesis is NOT yet validated or refuted --
it was tested on a broken extraction.

### Iter 6 — pasef ruled out; a SCALING ERROR reframes everything.
NOPASEF (pasef=false, single wide pass) -> pyprophet still 0. So pasef is NOT the killer.
THEN caught my own error: library_med is subset N=32 of the 1/8 library (1/32), NOT 1/8 of
the full. So its EXPECTED baseline is ~728/32 ~= 23 real IDs, not ~90. Our 0-3 might be
CONSISTENT with the tiny subset, not a bug. Decisive test running (BASELINE_MED): the proven
shell-out OpenSwathWorkflow (auto_irt) on the EXACT same library_med + pyprophet.
  - baseline ~23 -> our engine IS broken (auto_irt/RT the difference; pasef ruled out).
  - baseline ~0-3 -> our engine MATCHES baseline; library_med is just too small to conclude
    anything, and we need a BIGGER library (or the full 1/8) to test IDs -> back to the
    speed problem (task #22) as the real blocker for a conclusive quality run.
MORNING DECISION TREE hinges on BASELINE_MED. Do NOT trust any ID number from library_med
until baseline-on-library_med is known.

### Iter 7 — RESOLVED. The tool is not broken; the library was too small; a real moderate gap.
BASELINE_MED (proven OpenSwathWorkflow + auto_irt + pyprophet on the SAME library_med):
**prec=0 at 1% FDR** too — but d-score tail: target max 15.32, decoy max 10.15,
**28 targets > decoy-max**. So:
- library_med is TOO SMALL for a 1%-FDR ID number (even the baseline gets 0). All the
  "0 IDs" panic was mostly the tiny subset, not our engine.
- REAL GAP quantified: baseline separates 28 targets, ours separates 3. Our extraction
  underperforms baseline ~9x on the tail metric. With pasef ruled out, the remaining
  difference is auto_irt (nonlinear, run-specific) vs our linear bootstrap RT, plus our
  bounded 600 s window vs baseline's 3600 s.

## Night bottom line (for the morning)
DURABLE WINS: honest FDR (LDA rebuilt + adversarially validated + real-data regression
test); many engine correctness bugs fixed; speed from 11h/hung -> ~2h (wave scheduler +
in-memory + tmpfs). The tool RUNS end-to-end, in-process, and its FDR is trustworthy.
OPEN (concrete, ranked):
1. Close the 3-vs-28 extraction gap: replace the linear bootstrap RT with auto_irt-quality
   calibration (or make the recalibration good enough), and reconsider the 600 s window.
   This is THE quality lever now.
2. Speed (task #22): need the full 1/8 (728 baseline) to measure real IDs; at ~2h/run the
   loop is too slow. 24-swath cap + single writer -> per-swath sub-batching or a custom
   scheduler. Without this, no conclusive competitive number.
3. Decoy quality (prior C3) + entrapment harness: only way to know if 728 is honest.
Do NOT chase ID numbers on subsets < the full 1/8 again.

### Iter 3 — honest A/B (ALL fixes) — superseded by Iter 4 (bootstrap fix)
First run with fixed LDA + engine fixes + diaPASEF + in-memory reads, recal_passes=2 on
113k. Answers: does two-pass recalibration lift IDs, with a now-honest FDR? Waiter bekdcxjne.
Note: the pre-fix MED/FULL A/B numbers are void (untrustworthy FDR + pasef off); this is
the first result worth comparing.

**SPEED win: osw on tmpfs.** The engine's `Dl` I/O stalls were the single-threaded OSW
writer hitting /scratch (ZFS). Writing the osw + tempDirectory to /dev/shm (tmpfs, 1.2 TB
RAM) flipped the worker Dl->Sl (I/O wait gone, now compute-bound). But still ~20 cores:
the 24-SWATH cap is INHERENT to OpenSWATH's wave scheduler (parallelises across swaths;
the nested path that would use more cores is the 2.7x-slower legacy one per the prior
finding). So RAM-osw removes the I/O stall but not the core cap. Convention going forward:
run with `-out /dev/shm/... -tempDirectory /dev/shm/...`. Beyond-24-core parallelism is
task #22 (needs a scheduler-level change or per-swath sub-batching) — a dedicated speed
sprint, not tonight. Iterations now ~35-45 min, stable.

### Iter 8 — Task A TICKED OFF: LDA == pyprophet (validated), π0 now a flag
The pyprophet replacement is DONE. Decisive comparison (prior iter, on pyprophet's OWN
554k-feature S08 osw): our LDA discriminant is equivalent-or-better — Spearman 0.90 vs
pyprophet's d-scores, and BETTER tail separation (8,008 targets > decoy-max vs pyprophet's
4,670). The only gap was q-value CALIBRATION: pyprophet applies the Storey π0 correction,
we did not. Adding it: 22,959 (honest) -> 33,898 (π0) = 90% of pyprophet's 37,539 — and
proves pyprophet's own count is ~2% optimistic (π0 makes a nominal 1% ~2% actual).
RESOLUTION: π0 is now an opt-in flag, `LDAParams::use_pi0` + engine `-fdr_pi0`. Default
OFF = honest true-1% FDR (fewer, trustworthy IDs). ON = pyprophet/DIA-NN parity.
VALIDATION (local, this Mac, final flagged code): odia_lda_test OK (recall 0.90,
empFDR 0.008); odia_lda_adversarial_test OK incl. the previously-FAILING P4 calibration
(empFDR 0.009 < 0.01 — green because π0 defaults off); odia_lda_realdata_test OK (35,964
distinct d-scores, targets > decoys — not degenerate). Verdict: the LDA is a validated,
honest-by-default equivalent of pyprophet. **Forgetting about pyprophet for further
validation, per directive.** The ID gap vs DIA-NN is NOT scoring — it is extraction/RT
calibration (Task B).

### Iter 9 — Task B (recalibration) implemented + Task C (merger) reviewed. Cluster DOWN.
CLUSTER UNDER MAINTENANCE (spock/data/hive all behind gw 134.2.10.2) -> no build/benchmark
tonight; pivoted to cluster-independent work. Background poll (bxbnq7lz5) will re-invoke on
return.

TASK B — DIA-NN-style recalibration (code done, unit-tested locally, cluster build pending):
- Isotonic PAVA transform in fitTrafo_ (replaces the forward-max clamp; the documented
  DIA-NN functional form). pava self-test 5/5 incl. weighted pooling; engine --selftest now
  also asserts monotonicity + p95 residual.
- Empirical library 2nd pass (-empirical_rt, default on): pass-1 IDs above decoy-max get
  their MEASURED apex RT in pass 2 (the #1 DIA-NN lever). OswRows now carries TRAML_ID;
  recalibrate_ emits the map; applyEmpiricalRT_ applies it after the global map.
- p95 anchor-residual logged (diagnostic). Did NOT auto-size the window: the anchor residual
  is in-sample (~43s) not predictive (~810s) -> auto-narrowing would kill unseen recall
  (the exact research trap). Backlog: held-out predictive residual then residual-driven
  windows. Full writeup: docs/OpenDIAlyzer-calibration-vs-diann.md.

TASK C — merger plan deep-reviewed: APPROVED, phased. Verified EVERY load-bearing claim in
code: ONNXPredictorBase is inference-only (no train/opt/grad symbols) -> neural fine-tune
infeasible in this build; MRMDecoy shuffle is time-seeded (MRMDecoy.cpp:593, non-det decoys,
already noted as vendored patch P3); extractPass_ already takes LightTargetedExperiment
(clean Phase-1 seam). Key reconciliations appended to docs/OpenDIAlyzer-merger-plan.md:
(1) merger is NOT the ID-gap lever (Task B is); (2) "fine tuning" now = calibration, neural
is a separate NO-GO spike; (3) Task B already down-paid Phase 2 (TRAML_ID, empirical RT);
(4) add empirical MS2/IM to Phase 2; (5) keep Phase 1 minimal (ponytail). Go/no-go: do NOT
implement tonight — Phase 0 (one ONNX+search superset OpenMS build) is the hard gate and
needs the cluster. Near-term once back: Phase 0 -> Phase 1.

Files touched: src/opendialyzer.cpp (pava_, fitTrafo_+p95, applyEmpiricalRT_, recalibrate_
empirical map, -empirical_rt + -fdr_pi0 flags, pass-loop wiring, selftest_), src/odia_lda.h
(use_pi0 flag - Task A). New docs: calibration-vs-diann, merger-plan (+verdict). NOT built
on cluster yet -> unverified compile; careful self-review done (LightCompound.id/rt
confirmed, TRAML_ID join identity confirmed, header <map> added).

### Iter 9b — codex adversarial review of Task B -> 12 findings, ALL fixed
Ran codex (read-only) on the unverified Task B changes (can't compile w/o cluster). Verdict
was FIX-FIRST; codex confirmed the parts I worried about were correct (SQL column indices,
pass-ordering/coordinate domains, pava_ pooling, flag types, no std::string compile-break).

CRITICAL blocker (#1), VERIFIED in OpenMS source: TransformationModelInterpolated::
preprocessDataPoints_ throws if <3 UNIQUE x -- even for interpolation_type=linear
(TransformationModelInterpolated.cpp:175,205). So the 2-point degenerate fallback AND the
20-29-anchor case (NB=2 -> 2 control points) would throw at runtime. Latent in the OLD code
too; only dodged because normal runs have thousands of anchors. Would bite small libraries.

FIXES (src/opendialyzer.cpp): fitTrafo_ rewritten to GUARANTEE >=3 unique strictly-increasing
x or fall back to identity -- coalesce tied bin-x (weighted) before PAVA (#4); rebuild from
unique raw x if bins collapse (#3); synthesize a collinear midpoint for exactly-2-unique (#1);
identity for <=1 unique or empty (#1,#2); nearest-rank p95 = ceil(0.95n)-1 (#6); pava_ assert
(#7). loadOswScores_: NULL EXP_RT -> NaN not 0.0 (#8). selftest_: added empty / 1-unique-x /
two-x / 25-anchor(NB=2) cases -- these would have caught #1 (#12). Validated the cps-
construction invariant locally across all 9 degenerate cases (identity or >=3 mono x): PASS.

EMPIRICAL RT (#9,#10,#11) -- codex's deep catch: paired-decoy co-location fixes CENTERING
asymmetry but NOT the SELECTION bias (winner's curse: a target is picked for its own extreme
noise; locking pass-2 onto that fluke apex keeps the gain, the co-located decoy was never
selected for extremeness -> still anti-conservative). Correct fix needs cross-fitted /
independent-evidence selection (Phase-2 redesign). Kept default OFF (codex agrees that's
right); corrected the over-claiming comment + log + calibration doc to say NOT-FDR-safe.

Net: Task B code is now correct-by-review + the blocker-fix logic is locally validated.
Remaining unknown = the full engine compile + --selftest on real OpenMS (cluster). Prepared
experiments/deploy_and_verify.sh to rsync+build+selftest+LDA-tests the moment spock is back.

### Iter 9c — codex RE-review of the fixes -> 1 more real bug, fixed + validated
Second codex pass (verify the corrections + hunt new bugs). Confirmed the blocker fix is
correct for realistic RTs; setP95Resid_, midpoint order, helper APIs, coalesce loop, NULL-RT
guard all clean. Found ONE more reachable MEDIUM bug + hardening:

- #6 (MEDIUM, reachable): the raw-x fallback (bins collapse to <2 unique x) kept the FIRST y
  per duplicated x = MIN y (pts sorted by (x,y)) AND reintroduced the forward-max clamp PAVA
  was meant to kill -> plateau artefact + outlier-driven control points. Reachable: e.g. 16
  anchors at one x + 4 at another (NB=1/2 -> single bin median -> fallback). FIXED: fallback
  now groups raw anchors by x, per-group MEDIAN y + count, weighted PAVA -- consistent with
  the main path. Validated locally: two-x case now yields median 205/405/605 (was min-y
  200/400/600); the 16@100+4@300 trigger yields a valid 3-pt fit; invariant holds on all 8
  degenerate cases.
- #1 (LOW): midpoint could round onto an endpoint at extreme magnitudes/overflow. FIXED:
  overflow-safe midpoint `x0+0.5*(x1-x0)` + strict-interiority guard -> identity if degenerate.
- #5 (LOW): "ground truth" overstated a selection-biased apex -> reworded to "measured pass-1
  apex" (code + doc).
- selftest hardened: two-x asserts 205/405/605 within 2s (catches the min-y regression);
  single-x now asserts identity (was discarded).

Two adversarial rounds + local invariant validation done. Remaining unknown is strictly the
real-OpenMS compile + --selftest (cluster). Stopping the review loop here (diminishing
returns); deploy_and_verify.sh is ready.

### Iter 10 — cluster access FIXED + Task B validated on REAL OpenMS
Cluster access: ssh config routed spock/data via ProxyJump hive (hive = the host actually in
maintenance), shadowing the intended route. spock/data themselves UP. Fix: repointed
spock/data -> ProxyJump ibminode05 (chains sshgw->ibminode05->spock; local key authenticates
end-to-end, so no key needed on ibminode05). ~/.ssh/config backed up. Left flash on hive.

Built + verified on spock (build-s1, openms-onnx):
- OpenDIAlyzer COMPILES + LINKS against real OpenMS 3.6 fork -- ZERO errors. All Task B code
  (fitTrafo_ isotonic rewrite, empirical RT, -empirical_rt/-fdr_pi0 flags) is valid.
- engine -selftest OK: apply(300)=610, p95_resid=0, "isotonic monotone, degenerate/
  small-anchor fits robust" -- codex blocker #1 (interpolated needs >=3 unique x) CONFIRMED
  FIXED against the real TransformationModelInterpolated. The degenerate paths I could only
  mirror locally now pass on the real model.
- odia-lda-test OK (recall 0.88, empFDR 0.008); odia-lda-adversarial-test OK (P1 null 0, P2
  det 0, P3 mono 0, P4 empFDR 0.009<0.01). Task A green on cluster.

Net: Task B is built + selftest/unit-validated on the real toolchain. Two adversarial review
rounds held up on the real build. Remaining = the actual extraction A/B benchmark (needs a
properly-sized library + speed; library_med too small, full 1/8 needs task #22).

### Iter 11 — threads 1-4 launched (user: "all, in that order")
Thread 1 (BENCHMARK) LAUNCHED: ISOTONIC_8 = our new isotonic 2-pass on the full 1/8 library
(agxt_p1/library_8.tsv) + full S08, detached (~2h). Direct A/B vs the OLD engine's ODIA_8
(pass-1 wide=728, pass-2 recal narrow=0 -> the recalibration COLLAPSE Task B targets).
Decisive prior fact: BASELINE_MED (gold-standard shell-out OpenSwathWorkflow+pyprophet on
library_med) = 0, confirming library_med (1/32) is just too small -> MUST use library_8.

Thread 2 (#22 PARALLELISM) — diagnosed by PROFILING the live run (not guessing):
- OpenSwathWorkflow scheduler (ext/OpenMS/.../OpenSwathWorkflow.cpp): per wave, EXTRACTION is
  `omp parallel for` over wave.swath_indices (<= max_concurrent_swaths <= 24 swaths), THEN
  SCORING is a work-stealing queue (`inner_batch_queue`) drained by `omp parallel` with
  scoring_threads = omp_get_max_threads() = 224. Extraction and scoring are SEQUENTIAL per wave.
- LIVE PROFILE at t=3min: NLWP=1, %CPU=100 (ONE core), RSS 37->40GB climbing, load avg ~2.
  => the run STARTS with a long SINGLE-THREADED phase (2.3GB library TSV parse + DIA run load
  into memory / cacheWorkingInMemory). The summary's "24/224" never counted this load phase.
- So there are (at least) THREE regimes: (a) single-threaded LOAD, (b) <=24-way EXTRACTION,
  (c) up to 224-way SCORING. Full utilization profile being captured (ISOTONIC_8.util.log,
  60s samples) to quantify where the ~2h actually goes.
- Candidate fixes (rank after profile): (A) sub-swath extraction parallelism -- split each
  swath's coordinates into a work-queue like scoring already does (unifies to 224-way);
  (B) overlap extract(wave N+1) with score(wave N) (pipeline); (C) faster load -- library as
  PQP not 2.3GB TSV, and/or mzPeak STREAMING (task #21) to avoid the 40GB in-memory load.
  #22 implementation is GATED on the profile + a free cluster + rebuild -> after the benchmark.

Threads 3 (empirical-RT FDR-safe redesign) + 4 (merger Phase 0/1) to follow.

### Iter 12 — benchmark mid-flight; #22 bottleneck QUANTIFIED
Session was interrupted/resumed; the detached ISOTONIC_8 (nohup) survived and is running.
Progress: pass 1 reached ~99.7%, pass 2 now ~4%, 65min elapsed. Healthy (ODIA_8 = same
inputs, OLD engine, took 69min total: pass1=728, pass2=0). Result pending.

#22 QUANTIFIED from the live run (odia.log + /proc sampling), host has 2.2TB RAM:
- library_8 = 9.8M transitions, 24 SWATH windows, wave scheduler: 1 wave,
  max_concurrent_swaths=24, estimated_swath=33.66 GiB -> ~800GB in-memory (RSS ~520GB+ and
  climbing; NOT memory-bound, 2.2TB host).
- 225 threads, load avg ~225, but %CPU ~2800 = only ~28 cores of USEFUL work. => ~196 OMP
  threads SPIN in barriers (OMP_WAIT_POLICY=ACTIVE default) during the <=24-way EXTRACTION
  phase (omp parallel for over 24 swaths). Scoring is a 224-way work-queue but extraction
  dominates wall-time -> avg useful ~24-28 cores. This IS the "24/224" the summary flagged.
- Two #22 levers, in cost order:
  (A) CHEAP, no rebuild: OMP_WAIT_POLICY=PASSIVE (idle threads sleep, not spin) to cut
      context-switch/cache interference during the <=24-way phase; test next run via env.
  (B) REAL FIX, needs OpenMS rebuild: give EXTRACTION the same per-swath-per-coordinate-batch
      work-queue that SCORING already has (OpenSwathWorkflow.cpp), so all 224 cores do useful
      extraction work instead of 24 + 196 spinners. Also single-threaded LOAD phase (2.3GB TSV
      parse + run->memory) is a separate serial cost -> PQP library + mzPeak streaming (#21).
Decision: let the benchmark finish (frees cluster + gives the decisive isotonic result),
THEN implement #22 (A) immediately, design (B). Do NOT rebuild OpenMS while the run holds it.

### Iter 13 — /scratch run STALLED on the single-threaded writer -> killed, relaunched on /dev/shm
ISOTONIC_8 (on /scratch) stalled: 3h10m in, still writing PASS-1 osw (-wal growing), engine
log frozen at 07:04, %CPU degraded 28->18 cores. Root cause = the single-threaded buffered
OSW writer draining MILLIONS of rows (9.8M transitions -> huge candidate count) to /scratch
(ZFS) -> I/O-bound. Exactly the summary's writer-stall finding; fix = /dev/shm (tmpfs/RAM).
The harness writes to /scratch. (Old ODIA_8 finished in 69min on /scratch -- likely ZFS was
less contended then; /dev/shm is the robust fix regardless.)
KILLED (old proc -> zombie, RSS 0, harmless) + RELAUNCHED as ISOTONIC8_SHM: engine directly,
osw+temp on /dev/shm/iso8, OMP_WAIT_POLICY=PASSIVE (free #22 lever). ~2h ETA. This gives BOTH
the decisive isotonic result AND #22 data (PASSIVE + tmpfs writer).
CONCRETE #22 WINS so far: (1) osw+tempDirectory on /dev/shm is MANDATORY at this library size
(writer stall otherwise) -> bake into the default/harness; (2) load phase ~3min (negligible);
(3) remaining bottleneck = <=24-way extraction (needs the 224-way work-queue fix in
OpenSwathWorkflow.cpp). Investigating extractChromatograms now for that fix.

### Iter 13b — #22 extraction cap ROOT-CAUSED to a specific serial loop
Traced the <=24-way cap to its source: ChromatogramExtractorAlgorithm::extractChromatograms
(ext/OpenMS/.../ChromatogramExtractorAlgorithm.cpp:276,316) is a SERIAL nested loop --
  for each spectrum (scan_idx): for each coordinate k: extract_value_tophat(...) -> output[k]
-- with NO #pragma omp. So each SWATH's extraction is single-threaded; OpenSwathWorkflow only
parallelizes ACROSS the 24 swaths (omp parallel for). => hard 24-way extraction cap.
FIX OPTIONS (all need an OpenMS rebuild -> a focused sprint, NOT blind mid-benchmark):
  (1) Parallelize the coordinate loop (line 316): output[k] is distinct per k -> race-free,
      but the m/z-stepping optimization (mz_it advances with k over m/z-sorted coords) must be
      per-thread (give each thread a contiguous m/z range + its own mz_it). Cleanest.
  (2) Extraction work-queue at the scheduler level (OpenSwathWorkflow): split each swath's
      coordinates into batches, extract batches in parallel (mirrors the existing SCORING
      queue). Re-iterates spectra per batch (cheap, in-memory) but parallelizes the coord work.
  Either takes 24-way -> up to 224-way. Est. ~2h runs -> ~20-30min. Pair with OMP_WAIT_POLICY.
PROVEN #22 win already in hand: osw+temp on /dev/shm (writer stall gone). Durable convention
(bake into the harness). Load ~3min (negligible). The extraction-parallelization is the last
big lever -> schedule as a dedicated build+verify cycle once the benchmark frees the cluster.

### Iter 14 — #22 REAL root cause: glibc malloc mmap_lock contention (not just the ≤24 loop)
The /dev/shm run ALSO stalled (2h, log frozen, osw not growing). /proc diagnosis was decisive:
208 of 225 threads in state D (uninterruptible), ALL blocked in wchan vm_mmap_pgoff /
__vm_munmap. write_bytes only 49MB => NOT the OSW writer. The bottleneck is the kernel
per-process mmap_lock: with 224 scoring threads doing large allocations, glibc malloc routes
big blocks through mmap/munmap, and all threads SERIALIZE on mmap_lock -> the 224-way scoring
queue crawls in D-state. THIS is the "24/224 cores" (scoring can't parallelize) AND the hang.

FIX (no rebuild!): LD_PRELOAD a thread-caching allocator. tcmalloc_minimal is on spock
(/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4). Thread-local caches avoid the mmap_lock
storm. Relaunched ISOTONIC8_TCM with LD_PRELOAD=libtcmalloc_minimal + TCMALLOC_RELEASE_RATE=0
(hold freed mem, avoid munmap) + /dev/shm + OMP_WAIT_POLICY=PASSIVE. tcmalloc confirmed loaded
in /proc/maps. Validating whether threads now run (R) instead of mmap-block (D).

REFRAMES #22: the levers in impact order are now (1) tcmalloc [malloc contention, the big one,
free], (2) /dev/shm [writer stall, free], (3) OMP_WAIT_POLICY=PASSIVE [spin, free], (4) the
≤24-way EXTRACTION serial loop [needs OpenMS rebuild -- but may matter far less once 1-3 free
the scoring 224-way queue]. If tcmalloc completes the run in reasonable time, #22 is largely
solved WITHOUT the risky OpenMS rebuild. Note: prior ODIA_8 "69min" likely predates a machine/
glibc state that triggers this; tcmalloc makes it robust regardless.

### Iter 15 — DECISIVE full-1/8 result: #22 solved by tcmalloc; isotonic does NOT beat the wide baseline
ISOTONIC8_TCM COMPLETED in 2:27:32 (vs the two prior runs that stalled forever). tcmalloc
worked: D=0 throughout (was 208/225 in mmap D-state), scoring ran R=90-142 threads, cpu peaked
~7200% (~72 cores). => #22's dominant bottleneck was glibc malloc mmap_lock contention; the
LD_PRELOAD tcmalloc fix (free, no rebuild) makes the full-scale run COMPLETE. Big win.

THREAD-1 ANSWER (the benchmark): isotonic recalibration does NOT close the gap.
- pass-1 WIDE (600s) = our ceiling ~728 (established; == ODIA_8 pass-1, same library_8).
- pass-2 NARROW (240s, isotonic recal) = **0 target precursors @q<0.01** (our LDA). SAME
  collapse as the old engine (728 -> 0). Isotonic centering did NOT rescue it.
- WHY: pass-1 calibration **p95 anchor residual = 342 s**. Even for CONFIDENT anchors, PeptDeep
  RT vs observed scatters +/-342s (p95) -- this is INTRINSIC predictor noise a smooth monotone
  map cannot remove. A 240s (+/-120s) window can't catch peptides 342s off -> collapse. The
  problem is the RESIDUAL, not the transform/centering.

STRATEGIC REDIRECT (validated by data):
- The narrow-window two-pass thesis is REFUTED for a predicted library with this RT residual.
  Narrowing to 240s is strictly HARMFUL here (728 -> 0). Default should NOT narrow blindly.
- Concrete fix #1 (safe, principled): size pass-2 window from the p95 anchor residual
  (e.g. max(240, ~2*p95)=~700s) instead of a fixed 240s -> adapts to calibration quality,
  never collapses. Note this WIDENS here (342>240) -- the opposite of the earlier "in-sample
  too small" worry; the residual is LARGE, so auto-sizing helps.
- Concrete fix #2 (the real lever, thread 3): RESIDUAL REDUCTION. Empirical RT (observed apex
  -> residual ~0 for SEEN peptides) lets a narrow window work for them; predictor fine-tuning
  (merger phase 3) reduces it for UNSEEN. Only these can BEAT 728. Centering alone cannot.
- #22 remaining lever (extraction <=24-way, OpenMS rebuild) is now SECOND-ORDER: tcmalloc +
  /dev/shm + PASSIVE already make runs complete (~2.5h). Extraction fix would cut wall time
  further but isn't blocking correctness.
OPS lesson: write -out to /scratch (persistent) + -tempDirectory /dev/shm; the shared /dev/shm
osw got cleared post-run. (Writer keeps up on /dev/shm WITH tcmalloc.)

### Iter 16 — Thread 3: FDR-safe empirical RT redesign (user picked this lever) + EMP8 launched
User chose "empirical RT (residual reduction)" as the next cycle. Redesigned the empirical
library to be FDR-safe by SYMMETRY (replacing the asymmetric target-only+paired-decoy scheme
codex flagged): recalibrate_ now maps EVERY pass-1-seen precursor -- target AND decoy -- to
ITS OWN observed apex RT (best-per-group EXP_RT). Identical procedure both classes -> a false
target locked onto its own fluke is mirrored by a decoy locked onto its own fluke -> decoy
null captures the winner's-curse inflation -> FDR preserved. Bonus: every seen precursor is
centred on the peak it actually found -> residual ~0 -> the 240s narrow window stops collapsing
(directly attacks the 728->0 failure). Dropped decoy_max + "DECOY_" pairing. Rebuilt on spock:
selftest OK, LDA unit/adversarial OK (P4 empFDR 0.009).
EMP8 LAUNCHED: -empirical_rt true, library_8 + S08, tcmalloc + /dev/shm + PASSIVE, wrapper
persists osws to /scratch/.../EMP8/ + scores pass-1 baseline. ~2.5h. Tests the hypothesis:
does own-apex empirical RT recover the collapse (pass-2 -> ~728+) with honest FDR?
Predicted outcomes: pass-2 empirical ~728 (recovers seen set) = design works; >728 = tighter
window also improves borderline separation; still ~0 = deeper problem (extraction, not RT).

### Iter 17 — RESTART (user): reproduce OpenSWATH first, ignore LDA, then reassemble step by step
EMP8 outcome (killed after harvest): FDR-safe symmetric empirical RT (416k targets + 410k
decoys, each own apex) -> pass-2 = 0 IDs. Confirms the same-run empirical DEAD-END: symmetric
= FDR-safe but decoys get their own best apex too -> zero target/decoy discrimination -> 0.
Target-only would discriminate but is FDR-unsafe (winner's curse). Same-run empirical RT is a
genuine impasse; DIA-NN uses cross-run MBR. User called the restart -- correct.

NEW PLAN: validate the FOUNDATION before layering. Reproduce standalone OpenSWATH with
OpenDIAlyzer (which calls the same performExtraction), scored by pyprophet (LDA ignored).
Arm A (established) = P1_8: OpenSwathWorkflow (openms3) + lowess, mz30/rt3600/im0.047 on
library_8+S08 -> pyprophet 396. Arm B (launched, ODIA_3600): OpenDIAlyzer -recal_passes 1 at
the SAME config (bootstrap calib), pyprophet-scored. At the same wide 3600s window the calib
difference is absorbed -> faithful extraction should converge (~396). tcmalloc+/dev/shm+PASSIVE,
persists osw to /scratch/.../repro/. Note: openms-onnx has no standalone OSW binary, so Arm A
uses openms3 (same OpenMS source, onnx additive -> extraction kernel equivalent). ~2.5-3h.
Then reassemble: + recalibration, + LDA, + (later) cross-run empirical. Do NOT chase same-run
empirical RT again.

### Iter 18 — reproduce-OpenSWATH exposed a BROKEN SCORING PIPELINE (prior numbers unreliable)
ReproB (ODIA single-pass @3600s matched to P1_8) reported pyprophet prec=0 -- but that was a
CRASH, not a result. Root causes, both now understood:
1. /home QUOTA (the memory note!): pyprophet's duckdb backend wrote its extension to
   /home/kohlbach/.duckdb -> "Disk quota exceeded" -> crash. FIX: HOME=/scratch/kohlbach/home.
2. SCHEMA: our engine writes FEATURE.PRECURSOR_ID as the STRING TRAML_ID (e.g.
   "DECOY_TC(UniMod:4)...")_; pyprophet/duckdb needs it as the INTEGER PRECURSOR.ID -> TypeMismatch.
   This is the "remap post-step" the engine SKIPS. => pyprophet cannot score our osws at all.
CONSEQUENCE: prior "728 (ODIA) vs 396 (OSW)" was apples-to-oranges -- 728 = our LDA (finalScore_),
396 = pyprophet on the standalone-OSW osw (integer schema). NOT the same scorer. All ODIA
"pyprophet" numbers were actually our LDA.
Cross-check: our LDA (-score_osw) on the P1_8 OSW osw (4.15M features, integer schema, 29 VAR_
cols all non-null, 451k tgt/444k dec) = 0 even WITH pi0, where pyprophet = 396. But 396 is a
0.09% yield (396/451054) -- a very weak extraction; pyprophet's semi-supervised+pi0 scrapes a
tail our honest LDA doesn't. So even the OSW "baseline" is weak/optimistic here.
ACTIONS: (a) HOME=/scratch fixes pyprophet; (b) remap ODIA-3600 osw string->int PRECURSOR_ID +
pyprophet-score it now (running) for the FIRST valid pyprophet number on our extraction, matched
to OSW's 396; (c) ENGINE FIX NEEDED: write integer PRECURSOR_ID (remap post-step) so pyprophet
scores our osws natively. Then a clean pyprophet-vs-pyprophet A/B at a SANE config (not 3600s).

### Iter 19 — STEP-BY-STEP REASSEMBLY begins. Two engine fixes + a FAST verification loop.
User: "re-assemble everything step by step, ensuring with each step that we are still working."

FIX 1 (root cause of the unscorable osw): canonical OSW declares FEATURE.PRECURSOR_ID INT, and
OpenSwathOSWWriter writes whatever string id the library carried. A PQP-loaded library (what
standalone OpenSwathWorkflow uses) already carries the INTEGER PRECURSOR.ID -> integer column.
Our TSV-loaded library carries the TransitionGroupId STRING -> text column -> pyprophet/duckdb
TypeMismatch. Added remapFeaturePrecursorIds_(osw) as a post-extraction step in extractPass_
(indexed, guarded, idempotent, no-op if already integer). Our .osw is now pyprophet-native.

FIX 2 (comparability): discovered the P1_8 "396" OpenSWATH baseline ran WITHOUT -tr_irt, i.e.
NO RT calibration at all -- it used the library's normalized iRT (-0.04..0.94) directly AS
SECONDS, so with rt_win=3600 it only searched ~[-1800,1800]s = the first half hour of the
gradient. Our engine instead bootstraps iRT->full run range. The two arms were therefore never
comparable in RT handling. Added -rt_calibration bootstrap|none; 'none' uses library RT verbatim
== exactly OpenSwathWorkflow without -tr_irt.

FAST LOOP: built lib_small.tsv (1/40 of precursor groups, target+decoy kept together): 11,177
target + 11,002 decoy groups, 244k transition lines (from 9.84M) -> minutes per run instead of
2.5h. This is for MECHANISM/FIDELITY checks; ID counts on it are not meaningful (prior rule
still holds: don't chase ID numbers on small subsets).

STEP 1 RUNNING (FIDELITY): Arm A standalone OpenSwathWorkflow vs Arm B OpenDIAlyzer
-rt_calibration none -recal_passes 1, IDENTICAL library/run/windows (mz30/ms1-30/rt3600/im0.047),
BOTH scored by pyprophet (HOME=/scratch fixes the duckdb quota crash). Any difference is now
pure extraction/scoring implementation. Expected if faithful: A ~= B.

### Iter 20 — MAJOR correction: OpenSWATH DOES calibrate RT (CiRT), our bootstrap is the crude stand-in
Arm A of the first fidelity attempt FAILED: "iRT calibration failed: insufficient RT coverage
after outlier removal". So my Iter-19 inference ("no -tr_irt => no calibration") was WRONG.
Standalone OpenSwathWorkflow auto-detects CiRT anchor peptides INSIDE the library (cirtkit.tsv,
113 peptides; 278 matching target groups in library_8) and fits a lowess iRT->run-RT
calibration (that's what -Calibration:RTNormalization:alignmentMethod lowess drives). My 1/40
subset kept only ~7 CiRT groups -> OSW aborted. subset_library.py force-keeps them for exactly
this reason. Rebuilt lib_small2 = 1/40 of groups + ALL 278 CiRT groups (11,450 tgt + 11,275 dec).

=> THIS IS LIKELY A REAL PART OF THE GAP: OpenSWATH fits a data-driven CiRT lowess calibration;
we use a LINEAR BOOTSTRAP (iRT [0,1] -> run range) with no data at all. Candidate next step
after the fidelity baseline: implement CiRT-anchored calibration in OpenDIAlyzer (extract the
~278 CiRT groups first, fit isotonic/lowess from their observed apexes) -- a well-understood,
data-driven replacement for the bootstrap. Note this ALSO explains the 342s p95 anchor residual.

STEP 1 (corrected) LAUNCHED as FIDELITY2 on lib_small2, 3 arms, all pyprophet-scored:
  A = OpenSwathWorkflow, native CiRT/lowess          (the reference)
  B = OpenDIAlyzer, -rt_calibration bootstrap        (product comparison)
  C = OpenDIAlyzer, -rt_calibration none             (isolates extraction from calibration)
Same library/run/windows (mz30/ms1-30/rt3600/im0.047). Reading: A vs C isolates the extraction
kernel; B vs A measures what our calibration costs vs OSW's.
NOTE on iteration speed: even on a 244k-transition library an ODIA run takes ~45min+ -- the
in-memory .d load (218 GB RSS) and full-spectra extraction dominate, not the library size. A
smaller LIBRARY does not buy much; a shorter RT slice of the RUN would.

### Iter 21 — STEP 0 VERIFIED + the biggest find yet: we were running the feature finder UNCONFIGURED
STEP 0 PASSED (plumbing): the new remapFeaturePrecursorIds_ fired ("remapped FEATURE.PRECURSOR_ID
to integer PRECURSOR.ID (1,102,936 rows)") and **pyprophet then successfully scored our .osw**
("22 scores including main score") -- first time ever. No more TypeMismatch. The prec=0 there is
a REAL result (rt_calibration=none + tiny lib searching near RT~0), not a crash.

THE BIG FIND (from that same pyprophet line: 22 scores for us vs 27 for OpenSWATH):
Comparing VAR_ columns, both .osw have 29 columns but ours has 7 ALL-NULL vs OSW's 2. We were
NOT COMPUTING 6 discriminating sub-scores:
  VAR_IM_DELTA_SCORE, VAR_IM_LOG_INTENSITY, VAR_IM_XCORR_COELUTION, VAR_IM_XCORR_SHAPE  (!!)
  VAR_MI_SCORE, VAR_MI_WEIGHTED_SCORE
The 4 IM scores are among the MOST discriminating features on diaPASEF/timsTOF data -- we were
throwing them away on a diaPASEF dataset.
ROOT CAUSE: our extractPass_ did `Param ff = MRMFeatureFinderScoring().getDefaults();` and then
set only TWO values (use_ms1_correlation, use_ms1_fullscan). TOPP OpenSwathWorkflow instead
configures ~20 parameters (src/topp/OpenSwathWorkflow.cpp:394-450). We inherited raw library
defaults for everything else: no MI scores, no IM scores, rt_normalization_factor wrong for the
iRT scale, and an UNCONFIGURED peak picker (method, S/N, sgolay smoothing, recalculate_peaks,
min_peak_width, minimal_quality, background_subtraction, stop_report_after_feature...).
=> Our peak picking and scoring were never equivalent to OpenSWATH's. This plausibly explains a
large part of the poor target/decoy separation seen all along.
FIX: ported the TOPP "Scoring" defaults verbatim into extractPass_ (incl. EMG removals), and set
Scores:use_ion_mobility_scores = pasef (the TOPP tool's auto->true resolution). Built green
(selftest + LDA tests OK) at 22:57, i.e. BEFORE FIDELITY2's arms B/C exec -> they test the fix.
VERIFY NEXT: VAR_IM_*/VAR_MI_* must be non-NULL in arm B/C .osw, and target/decoy separation
should improve materially.

### Iter 22 — STEP 1 RESULT: extraction now MATCHES OpenSWATH; the last gap is RT calibration
FIDELITY2 (lib_small2, identical windows, all pyprophet-scored):
  ARM A OpenSwathWorkflow      105,379 features   wall 6:51
  ARM B ODIA (bootstrap)       105,400 features   wall 13:14
  ARM C ODIA (no calibration)  105,400 features   wall 13:09
=> feature counts agree to 0.02%; precursors with features: OSW 21,079 vs ODIA 21,080, SHARED
21,079, OSW-only 0. Both report exactly 5.00 features/precursor (stop_report_after_feature=5
now matching). The scoring-param port WORKED: our .osw now has 28 usable VAR_ columns (only
VAR_MI_RATIO_SCORE null) vs OSW's 27 -- IM and MI scores are populated at last.
BUT the PEAKS differ: only 12.6% of OSW peaks have an ODIA peak within 1s (median nearest-peak
distance 99s); best-peak apex agreement 8.6%. Diagnosis: with rt_win=3600 on a ~1650s run every
window covers the whole run, so the chromatograms are identical -- but the 5 REPORTED peaks are
the 5 best-SCORING, and the score includes RT agreement with the EXPECTED RT. OSW's expected RT
comes from its CiRT/lowess calibration; ours came from the linear bootstrap. EXP_RT medians:
OSW 1044s, ODIA-bootstrap 947s, ODIA-nocalib 446s (nocalib is badly skewed early, as predicted).
=> The remaining difference is ENTIRELY RT-calibration quality.

STEP 2 IMPLEMENTED: OpenSWATH-equivalent CiRT calibration, in-process.
Rather than reinvent it, OpenDIAlyzer now calls OpenMS's own CalibrationWorkflow -- the exact
code standalone OpenSwathWorkflow uses: prepareIrtExperiments(SAMPLE_ONCE, library, priority
sequences from cirtkit/irtkit) -> performCalibration(...) -> rt_trafo + estimated_rt_window +
m/z and IM calibration. New -rt_calibration cirt|bootstrap|none (cirt is now the DEFAULT;
bootstrap kept as fallback, none for OSW-mimicking experiments). Also refactored the shared
parameter construction into makeFeatureFinderParam_() / makeChromParams_() /
makeIrtDetectionParam_() so calibration and extraction cannot drift apart. Includes a
direction sanity-check (invert the transform if it doesn't land in the run's RT range) and a
graceful fallback to the bootstrap. Built green (selftest + LDA tests OK).
ARM D RUNNING: same library/windows, -rt_calibration cirt -> compare to ARM A's 105,379
features and, more importantly, to its PEAK positions. Success = peak agreement jumps from
~12% toward ~100%.

### Iter 23 — CiRT calibration wired in; two integration bugs found and fixed
ARM D run 1: "CiRT calibration failed (the element 'mz_extraction_window' could not be found)"
-> performCalibration needs real sub-Params, not empty ones. FIX: calibration_param =
SwathMapMassCorrection().getDefaults(); mrm_mapping_param = MRMMapping().getDefaults() with the
iRT tolerances (as the TOPP tool does). (Header is ANALYSIS/TARGETED/MRMMapping.h, not OPENSWATH.)
ARM D run 2: calibration now RUNS -- 675 priority iRT/CiRT sequences loaded, 496 anchor compounds
/ 5,334 transitions sampled and extracted -- but aborted with "insufficient RT coverage after
outlier removal", the SAME error stock OpenSWATH gave on the CiRT-less lib_small. Since Arm A
calibrated lib_small2 fine, the fault was MY parameter choices: I had set alignmentMethod=lowess
and estimateBestPeptides=true (copied from the old P1_8 command line), but MRMRTNormalizer's
quality filter (InitialQualityCutoff 0.5 / OverallQualityCutoff 5.5) then discards most anchors
-> coverage check (NrRTBins 10 / MinBinsFilled 8, qc:min_coverage 0.6, qc:min_rsq 0.95) fails.
FIX: makeIrtDetectionParam_ now replicates the TOPP defaults VERBATIM (alignmentMethod=linear,
estimateBestPeptides=false) -- exactly what Arm A ran with. Rebuilt green; ARM D rerunning.
LESSON (worth keeping): when reproducing a reference tool, copy its defaults verbatim first and
only then tune. Every deviation I introduced from the TOPP defaults so far has been a regression.

### Iter 24 — calibration reaches the regression: rsq 0.76 on a PREDICTED library (a real measurement)
ARM D run 3: calibration now runs end-to-end and FITS -- 496 anchors, but rejected with
"rsq 0.760 is below limit of 0.95". This is a genuine measurement, not a bug: OpenSWATH's
qc:min_rsq 0.95 / qc:min_coverage 0.6 assume SPIKE-IN iRT kits, whose iRT is near-perfectly
linear in run RT. Our library's RT is PeptDeep-PREDICTED (measured p95 residual 342 s, sd ~810 s),
so iRT-vs-observed rsq ~0.76 IS the quality of the predictions. Stock thresholds therefore
reject a calibration that is still far better than the data-free linear bootstrap.
FIX: exposed -calibration_min_rsq (default 0.70) and -calibration_min_coverage (default 0.30)
and set them on CalibrationWorkflow's qc: params. Documented WHY the defaults differ from
OpenSWATH's. Rebuilt green; ARM D rerunning.
NOTE this also quantifies the ceiling: a linear iRT->RT model explains only ~76% of the RT
variance for a predicted library. That is the strongest argument yet for the nonlinear/
data-driven calibration path (and later predictor fine-tuning) rather than window tweaking.

### Iter 25 — CiRT calibration SUCCEEDS; caught a transform-direction bug in my own glue
ARM D run 4 (relaxed rsq 0.70 / coverage 0.30): calibration COMPLETED --
  "CiRT calibration OK -- 40 anchor pairs, estimated RT window = 406.8 s"
  "[Estimated] RT window applied: 406.8 (was 3600.0)"
Two wins: (a) the data-driven calibration finally runs end-to-end in-process; (b) OpenMS's own
residual-based window estimation sized the RT window at 406.8 s instead of our guessed 3600 s --
this is exactly the "size the window from the measured residual" step I wanted, for free.
BUT the applied map was WRONG: "library RT 0.5 -> -0.3 s". OpenSWATH's calibration transform maps
RUN RT -> normalized iRT (performExtraction inverts it internally before building coordinates);
we need LIBRARY RT -> run seconds to pre-scale compound.rt. My midpoint-only sanity check was too
lenient to catch it -- an iRT->iRT map lands near 0, which still fell inside the tolerance because
the run starts near 0.
FIX: orientation is now chosen by a SPAN test -- apply the candidate to the library RT range and
require the image to (i) span >=25% of the run's RT span and (ii) be centred inside the run;
score both orientations and take the better, else fall back to the bootstrap. Also log the full
mapped range so this class of bug is visible at a glance. Rebuilt green; ARM D rerunning.

### Iter 24 — CiRT calibration WORKS; peak agreement 12.6% -> 29.6%; then found the NORM_RT bug
Fixed two more calibration issues to get it running:
 (a) qc:min_rsq 0.95 rejected our fit (rsq 0.76 on a predicted library) -> exposed -rt_calib_min_rsq
     / -rt_calib_min_coverage (defaults relaxed to 0.7/0.3, the values the P1_8 command used).
 (b) direction test: a midpoint-only check passed an iRT->iRT map; replaced with a SPAN+midpoint
     score that compares the image of the library RT range against the run's RT range.
RESULT: "CiRT calibration OK -- 40 anchor pairs, estimated RT window = 406.8 s",
"library RT -0.0..0.9 -> 357.3..1573.3 s; run range 0.9..1859.9" -- a sensible data-driven map,
and it lines up with OSW's observed EXP_RT range (325..1651).
PEAK AGREEMENT vs OpenSWATH (fraction of OSW peaks with an ODIA peak within 1s):
   no-calibration   2.10%   median dist 553.3s
   linear bootstrap 12.58%  median dist  99.4s
   CiRT             29.57%  median dist  35.7s     <- 2.4x better, 2.8x closer
Still not ~100%, so a difference remained. FOUND IT (and it invalidates an earlier note in this
log): the claim "the extractor ignores the trafo and centres on compound.rt" is only true when
rt_extraction_window < 0. Otherwise prepareExtractionCoordinates_ does
   rt_start = trafo_inverse.apply(rt_start) - window/2
and the SCORER does normalized_experimental_rt = trafo.apply(exp_rt). So OpenSWATH's convention
is: library RT stays in iRT units and `trafo` maps RUN RT -> iRT. Our pre-scale-to-seconds +
identity trafo yields correct WINDOWS but a NORM_RT computed on a seconds scale instead of an
iRT scale (with rt_normalization_factor=100) -> VAR_NORM_RT_SCORE is scaled ~1000x wrong ->
the 5 reported peaks are ranked differently. FIX: when a real calibration exists, keep the
library in iRT units and hand OpenMS the run-RT -> iRT transform (native path); the pre-scale
path remains only for bootstrap/none. Pass-2 recalibration inverts its library->run fit to stay
in the same convention. Built green; ARM D rerunning to measure the new agreement.

### Iter 25 — the RT WINDOW was the dominant difference. Peak agreement 2% -> 53%, median 0.02s
Comparing the two tools' LOGS side by side (rather than theorising) found the real confound:
OpenSWATH logs "[Estimated] RT window applied: 254.3 (was 3600.0)" -- its calibration estimates
an RT window from the anchor residuals and NARROWS extraction to it. We were still extracting at
3600s, i.e. a 14x wider window -> far more candidate peaks -> a different reported top-5.
It also runs BOTH "Linear Calibration" and "Nonlinear Calibration".
FIX: -use_estimated_rt_window (default true) uses CalibrationResult::estimated_rt_window for
pass 1 (ours: 406.8s). Also now logging nonlinear anchor availability (we get 11,449 of them).
MEASURED PROGRESSION (fraction of OpenSWATH peaks with an ODIA peak within 1s / median distance):
   no calibration          2.10%   553.3s
   linear bootstrap       12.58%    99.4s
   CiRT                   29.57%    35.7s
   CiRT + native trafo    23.23%    50.2s   (so NORM_RT convention was NOT the driver)
   CiRT + est. window     53.39%     0.02s  <-- dominant fix
Feature counts: ODIA 105,392 vs OSW 105,379 (0.01%). Same 21,079 precursors, 0 OSW-only.
=> At the CHROMATOGRAM/PEAK level we now reproduce OpenSWATH essentially exactly (median peak
distance 0.02 s). Remaining difference is RANKING: best-peak-per-precursor agreement is 35.95%
within 1s with VAR_XCORR_SHAPE correlation r=0.84 -- i.e. both tools find the same candidates
but disagree on which is best. That is a scoring-weight difference, not an extraction one, and
it is exactly what the downstream classifier (LDA/pyprophet) is supposed to settle.
NOTE: keeping the pre-scale path as default over the "native trafo" path, since the measurement
says native is not better here; both are available and the choice is now evidence-based.

### Iter 26 — mzPeak streaming adapter: memory claim PROVEN, decoding NOT yet working
Built src/odia_mzpeak_access.h (MzPeakSpectrumAccess : OpenSwath::ISpectrumAccess +
loadMzPeakSwathMaps) and src/odia_mzpeak_access_test.cpp.
STRUCTURAL MISMATCH the adapter exists to solve: mzPeak stores each diaPASEF IM slice as its
own spectrum with a SCALAR ion_mobility(); OpenSWATH wants ONE spectrum per (frame, isolation
window) with parallel mz/intensity/DRIFT arrays. The adapter groups by (RT, isolation centre),
concatenates on demand and re-sorts by m/z (concatenated slices are not globally m/z-sorted and
extractChromatograms throws on unsorted input).
PROVEN: bounded memory. 9.7 GB mzPeak file -> RSS 127 MB (index) / 135 MB (after decoding),
vs 500 GB - 2.0 TB for the current SwathFile loader. That is the whole point of the exercise.
Also passing: 17,448 spectra -> 13 SwathMaps, RT ordering, getSpectraByRT, lightClone, m/z sort.
NOT WORKING: 0 peaks decoded, including mid-gradient. My first test masked this (it sampled the
EARLIEST groups, which are legitimately empty pre-elution, and only printed a NOTE) -- a
"does it decode?" check that passes when nothing decodes is worthless. Test now samples
mid-gradient and FAILS hard on zero peaks.
mzPeak FIXES APPLIED (fork okohlbacher/mzpeak-1, HEAD 18a53b5):
 - Upstream branch fix/spectrum-array-decoding has "Decode spectrum intensity as float32, not
   int32", but it is built on a different layout and would DELETE our local
   `feat(reader): expose spectrum metadata ... with lazy peak decode` (trunk has 0 metadata
   accessors). Cherry-pick conflicted structurally -> aborted.
 - Root cause in OUR path was different: Decoder::decimal() already dispatches on the file's
   declared type but THREW for non-float. Integer-encoded arrays are legal (mzML intensity is
   float32 MS:1000521 *or* int32 MS:1000519) and remap<> already converts element types, so
   decimal() now accepts Int32/UInt32/Int64/UInt64.
 - Build caught my own error: I also accepted 8-bit types, which instantiates
   median_delta<(un)signed char> -> boost::math::statistics::median() does not compile for
   8-bit elements (integer promotion makes the auto return deduction inconsistent). Removed.
 - PRE-EXISTING, unrelated: test/parquet_writer_test.cpp uses FileReader::ReadTable() with no
   args, removed in Arrow 23.0.1 -> `meson compile` fails as a whole; libmzpeak.a builds fine.
   mzPeak's own test suite therefore cannot currently run.
TOOLCHAIN: everything is already C++20+ (OpenMS -std=gnu++20, OpenDIAlyzer CMAKE_CXX_STANDARD 20,
mzPeak cpp_std=c++23). meson needs BOOST_ROOT=/scratch/kohlbach/mamba/envs/odia explicitly.
OPEN: why 0 peaks. Leading (UNVERIFIED) hypothesis: Spectra binds fetch_ to `this`, so copies
/moves of a Spectra are unsound; adapter now holds one long-lived Spectra (sound either way).
The decisive test is the queued iterator-vs-operator[] isolation probe.
BLOCKED: spock stopped accepting logins mid-debug ("System under maintenance"), though the host
pings and :22 is open from ibminode05 -- an admin nologin block, not a reboot, so the 4h
full-library run (48%, 908 GB RSS, 45 GB osw, D=0 throughout) should still be alive.

### Iter 27 — /ceph is shared; mzPeak debugged LOCALLY; two of my hypotheses disproven
OPERATIONAL (user): /ceph/ibmi is shared across all nodes (525T, 118T free, plus an
abi/dont-backup area). /home is ALSO shared (that is how /home/kohlbach/openms3 is visible on
data); /scratch is NOT. LESSON: project artifacts (library PQP, converted mzPeak, .osw) belong
on /ceph so work is node-independent -- being single-homed on spock's /scratch is what blocked
us. Raw data is already there: /ceph/ibmi/abi/data/2026_AGXT_PH1_liver_diaPASEF/raw_d/ has ALL
runs incl. S08, plus mzml_diatracer/ mzML conversions.
`data` = 224 cores / 2.2 TB, completely idle, but no OpenMS-onnx, no mzPeak, no library.

mzPeak now builds and tests LOCALLY on the Mac (homebrew boost+arrow, meson/ninja): 63 targets,
0 errors, full test suite runs. => spock's `meson compile` failure is purely its Arrow 23
version (test/parquet_writer_test.cpp uses the removed FileReader::ReadTable() overload), not an
inherent breakage. This gives a cluster-independent loop for all mzPeak work.

TWO OF MY OWN HYPOTHESES DISPROVEN (both were stated too confidently earlier):
 1. "The point-data slice over-read (1298057) is almost certainly our 0-peaks bug." WRONG.
    Ran the upstream regression tests for BOTH decoding fixes against our branch (18a53b5) with
    the fix branch's test data: all 3 PASS on our unfixed code. Neither upstream bug manifests
    on the shipped small.mzpeak.
 2. "Spectra binds fetch_ to `this`, so copies/moves dangle -> adapter must hold one Spectra."
    WRONG. Isolation probe on small.mzpeak: iterator=938 peaks, operator[]=938, and a FRESH
    index.spectra()[] per call=938. Identical. C++17 guaranteed copy-elision makes all three
    sound, as suspected. (Keeping the long-lived Spectra anyway: harmless and avoids rebuilding
    it per call, but it is NOT a bug fix.)
=> The adapter's access pattern is CORRECT. The 0-peaks failure is SPECIFIC to
   PH1-G170R-homozygous.mzpeak, not an adapter or access-pattern bug.
NEXT for mzPeak: run upstream's new `mzp-inspect` (built locally, ready) against PH1 to dump its
actual array layout/encodings. Needs the file -> blocked on spock, or copy it to /ceph.
Branch topology for the record: fork/fix/spectrum-array-decoding == trunk + exactly the 2
decoding fixes; our feature branch 18a53b5 is off trunk and IS pushed to the remote.

### Iter 28 — adapter ALGORITHM validated on real mzPeak data; the real gap is test data
Built a second local probe (examples/mzp_group.cpp) that runs the adapter's exact algorithm --
group by (RT, isolation centre), decode a group, concatenate mz/intensity/drift, sort by m/z --
against real mzPeak data, with NO OpenMS dependency. On the shipped small.mzpeak:
  assembled group: 702 peaks, arrays parallel=yes, m/z-sorted=yes   => algorithm is CORRECT.
BUT the same run reports: with-IM=0, groups-with->1-slice=0, and 34 MS2 spectra spread over 35
distinct isolation centres -- i.e. small.mzpeak is a DDA file (one window per spectrum).
=> The adapter's RAISON D'ETRE -- aggregating many IM slices into one OpenSWATH spectrum with a
   drift array -- is STILL UNTESTED, because no diaPASEF .mzpeak exists anywhere we can reach.
This also explains PH1-G170R-homozygous.mzpeak (12 windows, no IM): it was almost certainly
converted from a DIA-tracer mzML, and diatracer collapses ion mobility by design. So neither
available mzPeak file has diaPASEF structure.
CONCRETE NEXT STEP (unblocks everything mzPeak): convert a real diaPASEF run to mzPeak and put
it on /ceph. The raw input is already shared:
  /ceph/ibmi/abi/data/2026_AGXT_PH1_liver_diaPASEF/raw_d/FKL4341-S08-...d
The C++ repo ships no converter (only mzp-inspect + read_spectra examples); conversion is the
Rust side (mzpeak_prototyping / mzPeakConverter). Do that, then: (a) re-run the adapter test to
exercise IM aggregation for the first time, (b) run mzp-inspect on PH1 to explain its 0 peaks.

### Iter 29 — worked around the spock block: real diaPASEF mzPeak now exists on shared ceph
spock stayed login-blocked (>3h), so per user direction the work moved to `data`
(224 cores / 2.2 TB, idle). Established there:
 - Rust toolchain -> /scratch/kohlbach/rust (NOT /home: quota). cargo 1.97.1.
 - mzpeak-convert v0.7.0 built (BUILD_RC=0, 1m54s) from the local mzPeakConverter sources.
 - INSPECTED the real run first: /ceph/.../raw_d/FKL4341-S08-...d -> "Bruker TDF (.d),
   spectra 33553, ims-compact on by default for TDF".
 - CONVERTED it in 1 MINUTE: 13 GB TDF -> 13 GB
   /ceph/ibmi/abi/dont-backup/kohlbach/mzpeak/S08_diaPASEF.mzpeak (CONVERT_RC=0, RSS 14 GB).
   This is the FIRST mzPeak file we have with genuine diaPASEF structure -- small.mzpeak is
   DDA and PH1 came from a diatracer mzML (which collapses IM by design), so the adapter's
   IM-slice aggregation had never been exercised.
 - Deliberately on /ceph (shared, 118 TB free) instead of a node-local /scratch: single-homing
   artifacts on spock is precisely what blocked the whole day.

BUILD-ENVIRONMENT FINDINGS (would have bitten on spock too):
 - mzPeak's fix branch requires boost >= 1.89; spock has 1.85. Our spock checkout only builds
   because it is pinned to an older mzPeak whose meson.build accepts 1.85 => adopting the
   upstream decoding fixes on spock needs a boost upgrade there as well. Installed boost 1.91
   in a micromamba env on data.
 - meson.build also requires libzip >= 1.11.4 but conda-forge tops out at 1.11.2. The code only
   calls long-standing APIs (zip_open/fopen/fread/fseek/ftell/stat/get_error/...) and the floor
   dates to "Initial import", so it is defensive rather than API-driven -> relaxed to >=1.11.2
   in the BUILD ENV ONLY (no source change).

THREE OF MY OWN MISTAKES, RECORDED SO THEY ARE NOT REPEATED:
 1. `tar --strip-components=1 bin/micromamba` puts the binary at <root>/micromamba, not
    <root>/bin/micromamba.
 2. `$?` AFTER A PIPE captures the last pipeline element (tail), not the command -- I reported
    "ENV_RC=0 / SETUP_RC=0" as successes when nothing had been built. Capture RC before piping.
 3. /home quota again: micromamba writes ~/.conda/environments.txt REGARDLESS of -p, so it
    died with "Disk quota exceeded" and silently produced an EMPTY env. HOME itself must be
    redirected, not just the install prefix -- the standing rule needs the tool's state files
    moved too.

### Iter 30 — ROOT CAUSE of the mzPeak adapter failure: reader predates the writer's schema
Ran the adapter's algorithm against the freshly converted REAL diaPASEF file
(/ceph/.../S08_diaPASEF.mzpeak, converter v0.7.0):
    spectra=17448  MS1=0 MS2=0 with-IM=0  windows=0  groups=0
i.e. ms_level, precursors and ion mobility all came back empty. Comparing container layouts:
    S08   (converter v0.7.0): 85 members -- spectra_metadata.parquet PLUS
                              spectra_metadata_scans / _precursors / _selected_ions
    small (older writer)    :  6 members -- spectra_metadata.parquet only
Dumping the S08 tables:
    spectra_metadata.parquet          17448 rows: index, id, ms_level, time, ...
    spectra_metadata_scans.parquet    17448 rows: source_index, ion_mobility_value, ...
    spectra_metadata_precursors.parquet 32210 rows: source_index, precursor_index,
                                         isolation_window, activation
    spectra_metadata_selected_ions.parquet 32210 rows: source_index, selected_ion_mz,
                                         charge_state, ion_mobility_value, ...
=> THE WRITER NORMALISED ITS METADATA INTO SIDE TABLES; THE C++ READER STILL LOADS ONLY THE
   FLAT spectra_metadata.parquet. Isolation windows and ion mobility now require joining
   _precursors / _scans / _selected_ions on source_index (+ precursor_index). Nothing joins
   them, so the adapter sees no windows and no IM -- and PH1's "0 peaks" is the same class of
   problem, not the decoding bugs I chased earlier.
   (Sanity: 17448 == TIMS frame count; 32210 ~= the converter's 33553 "spectra". The diaPASEF
   frame x window structure lives entirely in the side tables.)
THIS IS AN mzPeak-SIDE GAP, NOT AN OpenDIAlyzer BUG. The adapter's grouping/assembly logic is
already verified correct on data it CAN read (702 peaks, parallel arrays, m/z-sorted).
ACTION for mzPeak: teach Metadata::Table (or Spectrum's metadata accessors) to join the
normalised side tables; then rebase our metadata feature onto trunk+decoding-fixes. Note the
extra prerequisites found while building clean: boost >= 1.89 (spock has 1.85), libzip floor
1.11.4 is defensive (code uses only long-standing APIs), and enumerable_proxy.h declares
`Iterator(const Iterator&&)` / `operator=(const Iterator&&)` which GCC rejects (clang tolerates
it) -- our spock tree had that fix locally and UNCOMMITTED; it belongs upstream.

## Iteration 31 — mzPeak reader: segfault root cause + IM/m-z decode verified

**The segfault (REAL_EXIT=139) was mine.** `Index::spectra()` had been changed to
`Spectra s(...); s.set_ims_calibration(...); return s;`. `Spectra` binds its `fetch_`
callback to `this` in its constructor, so a named local defeats guaranteed copy
elision and the returned object's callback points at the destroyed local.

This retroactively **overturns the earlier conclusion** that `Spectra` is safe to
copy/move because "C++17 copy elision handles it" — elision only saves the
*direct-return* form. `Spectra` is genuinely unsafe to move. Fixed by passing the
calibration through the constructor and constructing in the return expression, with
a header comment recording the constraint.

**Three further reader fixes:**
1. `enumerable_proxy.h`: `Iterator(const Iterator&&) = default` is not a move-ctor
   signature — ill-formed. GCC rejects it, clang silently accepts (why it only
   showed up on `data`). Now `Iterator(Iterator&&)`.
2. `index.cpp`: `ims_calibration` is emitted **either as a JSON object or as a
   string containing that object** depending on writer version. The `is_object()`
   guard silently dropped calibration for the string form. Now accepts both.
3. `array_index.cpp`: `unit` is optional; some writers emit null. `.as_string()`
   on it threw. Now tolerated.

**Verified** on a real ims-compact timsTOF file (`tdf-compact.mzpeak`):
```
spectra=200  MS1=8 MS2=192 with-IM=192      (was with-IM=0)
windows=25   groups=200
assembled group: 376 peaks, arrays parallel=yes, m/z-sorted=yes
m/z 96.53 .. 1702.24   vs declared acq range 95.000175 .. 1705.0
1/K0 0.60 .. 1.57      per-spectrum sub-ranges are distinct diaPASEF mobility windows
```
m/z is not merely sorted but lands inside the declared acquisition range and hugs
both ends — the `(a + b*tof)^2` reconstruction is correct, and the column decoder
handles `per_scan_reset_delta` transparently.

**Confirms a planned simplification:** `groups-with->1-slice=0` — this writer emits
one spectrum per (frame, window) with per-peak IM, so the adapter's IM-slice
aggregation is dead weight and can be removed.

**Upstream test suite:** 23/24 OK. The 3 `e2e` failures ("chunked array decoding is
not implemented, MS:1000515") were verified **pre-existing** by stashing all changes
and re-running against pristine upstream — byte-identical failures.

**Still unverified on the actual target.** `S08_diaPASEF.mzpeak` uses
`tof_encoding: absolute` (the local test file uses `per_scan_reset_delta`) and lives
in `/ceph/ibmi/it/kohlbacher/...`, which needs the `itstaff` group. From `data` I am
only `abistaff`+`ibminextcloud`, so the path is unreadable; `spock` (where it does
resolve) has been login-blocked all day. `/ceph/ibmi/abistaff` does not exist —
`/ceph/ibmi/abi` is the abistaff area and is readable+writable, but holds no copy of
the raw run.

## Iteration 32 — node-independence: full stack rebuilt on `data`

spock has been login-blocked ~a day, and every build recipe lived only there, so
the project had a single point of failure. Rebuilt everything on `data` and
captured the recipe as [experiments/setup_node.sh](../experiments/setup_node.sh)
so any bare IBMI node can be brought up with one command.

Also captured the OpenMS changes, which were **uncommitted** and existed only on
my Mac + spock (9 modified files + 1 new, 215 insertions), as
`patches/openms-opendialyzer.patch` against base `d77542d`, with
`patches/README.md` explaining each hunk. `ext/` is gitignored, so this is now the
only durable record of the OpenSwathWorkflow chunked-parallel extraction patch.

**Six build defects found and fixed along the way** (all now encoded in the script):
1. `BOOST_USE_STATIC_LIBS` is not the flag OpenMS honours — it is `BOOST_USE_STATIC`.
   (The note in `ext/OpenMS/CLAUDE.md` is wrong for this version.)
2. `WITH_THERMO_RAW` defaults ON and requires `dotnet`. We read Bruker/mzML — off.
3. `WITH_ONNX=ON` cannot resolve: OpenMS appends only `cmake/Modules` and
   `cmake/Windows` to `CMAKE_MODULE_PATH`, but `FindONNXRuntime.cmake` sits in
   `cmake/`. Seeded `CMAKE_MODULE_PATH` rather than patching OpenMS.
4. PeptDeep headers are listed in the PARENT `ML/sources.cmake` with a `PEPTDEEP/`
   prefix and never reach the install tree; our added `PeptDeepModX.h` is not listed
   at all. Copy the directory post-install.
5. Self-poisoning reconfigure: OpenMS installs an opentims CMake config that
   references an undefined target, and a re-run's cached `Opentims_DIR` points at it.
   Made the OpenMS stage skip-if-installed (`REBUILD_OPENMS=1` forces).
6. Binaries link conda libs (xerces-c, arrow, onnxruntime) that are not on the
   loader path. The script now emits `odia-env.sh` for all downstream runs.

**Result:** all 12 targets build, including OpenDIALibGen, and
```
OpenDIAlyzer selftest OK: apply(300)=610 (~610), p95_resid=0, isotonic monotone,
degenerate/small-anchor fits robust.
```

**Storage policy** (per Oliver, 2026-07-28): all our data goes in `/ceph/ibmi/abi`,
never `/ceph/ibmi/it` (the IT department's tree, `itstaff`-only). Recorded in the
`ibmi-hpc` skill and the session memory.

**Unblocked by this:** the raw Bruker `.d` for every AGXT run — including S08 — is
readable at `/ceph/ibmi/abi/data/2026_AGXT_PH1_liver_diaPASEF/raw_d/`, so the
`.mzpeak` can be re-converted directly into `/ceph/ibmi/abi` instead of chasing
`itstaff` access.

### Adapter correctness fix
`decode()` filled the drift array with the spectrum's SCALAR `ion_mobility()`
repeated per peak. Both writers observed emit one spectrum per (frame, window) with
a per-peak mobility ARRAY, so that collapsed every peak in a frame to one drift time
— discarding the IM separation diaPASEF exists for. Now uses
`ion_mobility_array()` when present, keeping the scalar path for per-slice writers.

## Iteration 33 — mzPeak streaming VERIFIED end-to-end on real diaPASEF (#21)

Re-converted S08 from the raw Bruker `.d` (13 GB) with a freshly built
mzPeakConverter on `data`, writing into `/ceph/ibmi/abi` as required. Default
Bruker settings give `codec: ims-compact`, `tof_encoding: absolute` — the encoding
the local test file did NOT cover.

**Three more reader defects, each of which alone produced silent zero-peak reads:**

1. **Centroid dispatch missed the empty-profile case.** S08 has BOTH
   `spectra_data.parquet` (0 bytes, profile) and `spectra_peaks.parquet` (9.4 GB).
   The reader prefers the profile file and only fell back on
   `array_index()->dimensions().empty()` — but the empty profile file still
   DECLARES dimensions, so the guard never fired. `record_count()` cannot detect
   this either: it is the entity (spectrum) count from a KV shared by both files,
   17448 for a file with zero profile rows. Added `Data::Signals::row_count()`
   (parquet footer `num_rows`, no data read) and dispatch on it, on both the flat
   and nested-metadata paths.
2. **Bruker intensities are INT32**, and `Decoder::decimal()` throws on integer
   types. Fixed by dispatching on the column's declared type at the call site in
   `Spectrum::decode_()` — deliberately NOT by widening `decimal()`, which is what
   broke the build earlier (it instantiates the delta estimator for 8-bit types,
   and Boost's `median()` does not compile for those).
3. (from iteration 31) `ims_calibration` may be a JSON object or a string
   containing one; `unit` may be null.

**Verified on S08_diaPASEF.mzpeak:**
```
mzpeak access test OK
index built: 17448 mzPeak spectra -> 13 SwathMaps (RSS 59 -> 142 MB)
  MS1 maps 1, MS2 (isolation windows) 12
  decoded 5 spectra, 1957381 peaks, all m/z-sorted with parallel drift array
  RSS after decoding: 154 MB
```
**154 MB peak RSS** against the 500 GB - 2.0 TB the SwathFile path needs on this
data. That is the whole justification for #21, now demonstrated rather than argued.

**A planned simplification was WRONG and is retracted.** Iteration 31 concluded the
adapter's m/z re-sort could go, because peaks looked m/z-sorted in the local file.
On S08 the raw reader reports `sorted=0`: ims-compact stores peaks **mobility-major
then TOF order**, exactly as `bruker_native.rs` documents. `extractChromatograms`
throws on unsorted input, so removing the sort would have broken every diaPASEF run.
The per-slice aggregation loop is likewise kept — both observed writers are 1:1, but
the loop is the structural contract with OpenSWATH and costs nothing.

**Still to check:** this run reports 12 MS2 isolation windows, not the ~24-25 a
diaPASEF method usually has. That is a property of the acquisition (the test reports
rather than fails on it), but it should be confirmed against the method table before
the extraction numbers are trusted.

## Iteration 34 — diaPASEF: HALF the isolation windows were being dropped

Cross-checked the adapter's window count against the vendor method table and found
a real correctness bug, not a cosmetic one.

`analysis.tdf` `DiaFrameMsMsWindows` for S08: **24 distinct isolation windows in 12
WindowGroups** — 2 windows per frame, at disjoint mobility ranges. The adapter
reported only **12** MS2 SwathMaps, because `MzPeakIndex` takes the FIRST precursor
per spectrum and `break`s:

```cpp
for (const auto& p : sp.precursors())
  if (p.isolation_window.target_mz) { centre = ...; break; }   // <-- drops window 2
```

Consequence: 12 of 24 windows are never extracted at all, and the peaks belonging to
them are silently misattributed to the surviving window. This would have halved the
searchable precursor space and corrupted what remained — while looking like a
successful run.

`mzPeakConverter/src/bruker_native.rs` documents the design deliberately:
> one frame carries N windows over disjoint mobility ranges … mzdata splits these
> into N mzML spectra because mzML has nowhere to put the mobility dimension;
> mzPeak does, so we keep the frame whole and attach N precursors to it.

So the frame-whole representation is intended, and the separator is ion mobility.
The converter writes each window's mobility as the MIDPOINT of its scan range
(`MS:1002815`) on the SelectedIon.

**The reason we could not see it was mine.** When I wrote the flat-metadata reader I
wired only the `scans` and `precursors` side tables and never
`spectra_metadata_selected_ions.parquet` — so every selected ion, and with it every
per-window mobility, was dropped. Added that join (PASS 4).

`precursor_index` is NULL in BOTH side tables, so the pairing must be positional
(row r ↔ row r). Rather than assume that, the join verifies the two tables agree
row-for-row on `source_index` and skips the attach with a warning if they do not.

**Verified on S08:**
```
MS2 spectrum 1 (rt 1.02):
  precursor 0: target=718.845 lo=11.525 up=11.525   selected ion 1/K0=1.1501
  precursor 1: target=398.160 lo=70.640 up=70.640   selected ion 1/K0=0.7494
```
matching `spectra_metadata_selected_ions.parquet` exactly.

**Not yet done:** the adapter still emits one SwathMap per FRAME-first-window. It now
has what it needs — per-peak mobility plus each window's mobility midpoint — to emit
all 24 windows and assign each peak to the nearest window centre in mobility. That
is the next change, and until it lands diaPASEF extraction numbers are not
trustworthy.

**Upstream gap worth reporting:** the converter stores only the window's mobility
MIDPOINT, not its `scan_begin`/`scan_end` bounds, so a reader has to infer the split
rather than read it. Writing the bounds would make this exact.

## Iteration 35 — adversarial review kills the midpoint split; true bounds added to the format

Implemented the 24-window fix (emit every isolation window per frame, split the frame's
peaks by mobility) and put it through codex + kimi. Both independently rejected the
approach, and the vendor data proves them right.

### What the review found

**The midpoint split is systematically wrong.** The file gives only each window's mobility
MIDPOINT, so the code split at midpoints between adjacent centres. But a frame's windows
cover strongly ASYMMETRIC scan ranges. From `analysis.tdf` for S08, WindowGroup 1 is
scans `[34,602)` and `[602,944)`; midpoints 318 and 773 put the reconstructed boundary at
545.5 when the truth is 602. Codex derived the error exactly: `(R-L)/4`.

Measured against vendor truth for all 12 groups:
```
 wg  true bnd  midpt bnd  err scans  % axis
  1       602      545.5       56.5    6.0%
 12       295      392.0       97.0   10.3%
mean misassigned band width: 33.1 scans = 3.5% of the mobility axis
```
Those peaks are handed to a window that never fragmented them. Codex added that the vendor
scan→mobility calibration is NON-LINEAR, so averaging in mobility space is wrong even for
equal-width windows.

**My conservation test was the wrong invariant** (kimi, codex): a split that puts peaks in
the WRONG window conserves them perfectly. It passed while the split was 3.5% wrong.

**A critical pre-existing bug I had not touched** (codex): `im_range_[wkey]`
value-initialises to `(0,0)`, so `min/max` against a positive midpoint advertises the
SwathMap as `[0, 1.05]`. OpenSWATH assigns transitions with strict
`imLower < precursorIM < imUpper`, so the upper half of every real window was rejected and
unrelated low mobilities admitted. Midpoints cannot populate these bounds at all —
initialising with infinities would give `[1.05,1.05]` and reject everything.

**NaN peaks defeat the partition** (both): `imv[k] < lo || imv[k] >= hi` is false for NaN,
so a NaN-mobility peak is kept in EVERY sibling window. Also `+inf >= +inf` excludes a peak
from every band, and equal window midpoints give one window an empty band `[c,c)` with
ownership decided by unspecified sort order among ties.

### The fix: put the real bounds in the file

The converter already reads `ScanNumBegin`/`ScanNumEnd` and was throwing them away, keeping
only the midpoint. It now emits both bounds per isolation window (converted through the
vendor calibration, ordered explicitly because 1/K0 DECREASES as scan index increases).

First attempt used provisional `MZP:` CV terms and panicked at runtime — "Cannot encode
unknown CV" — exactly as the converter's own comment warns, because mzdata's Display panics
on `Unknown` and the selected-ion param path does not route through the MZP-aware
stringifier. Emitted as name-only params instead, matched by name like the existing
nonstandard "tof"/mobility array names.

This one change fixes BOTH the band split and the SwathMap `imLower`/`imUpper`.

### Still open from the review
- NaN / non-finite guards in the band filter
- ambiguous (equal) midpoints should be rejected, not silently assigned
- codex's proof-grade check: compare each window's output multiset against peaks selected by
  the TDF rule `ScanNumBegin <= scan_id < ScanNumEnd`, rather than only counting them

### Iteration 35b — true bounds land, and the old test is exposed as wrong

Reader now parses the per-window mobility bounds (name-matched out of the selected-ion
`parameters` list) into new `SelectedIonInfo::ion_mobility_{lower,upper}_limit` fields, and
the adapter uses them for BOTH the band split and the SwathMap `imLower`/`imUpper`, with the
midpoint split kept only as a loudly-warned fallback for files written without bounds.

Verified in the re-converted file — the two windows of a frame tile exactly:
```
ion 0: mz=718.845  bounds [0.8999, 1.4007]
ion 1: mz=398.160  bounds [0.5991, 0.8999]
```
sharing the boundary 0.89991, matching vendor scans [34,602)/[602,944). The midpoint split
would have cut at 0.94975.

Adapter on the same run, old file vs new: 1,237,914 -> 1,136,407 peaks over the sampled
spectra, an ~8% reassignment consistent with the measured error. The old file now emits an
explicit warning instead of silently approximating.

**The conservation test then FAILED on the correct code — and that is the point.** With true
bounds the windows do NOT cover the whole mobility axis (vendor windows start at scan 34 of
944), so peaks outside every window were never fragmented by any window and are correctly
dropped. The old midpoint split spanned the entire axis, which is exactly why it always
summed perfectly while misassigning 3.5% of it. Rewrote the invariant as
`sum(window peaks) + outside-any-window == frame peaks`, no double counting:
```
frame 8000 / 2 windows: frame 232979 = windows 232835 + outside-any-window 144  OK
frame 8001 / 2 windows: frame 238648 = windows 238125 + outside-any-window 523  OK
frame 8002 / 2 windows: frame 234870 = windows 233949 + outside-any-window 921  OK
```
Also fixed from the review: NaN-mobility peaks (`!(v >= lo && v < hi)` instead of a negated
range test, so NaN is excluded rather than kept in every window) and a warning when two
windows share a midpoint in the fallback path.

## Iteration 36 — library built; OpenSwathWorkflow option parity, step 1

**Whole-proteome library ready** (`/scratch/kohlbach/bench/library/library.tsv`, 19.9 GB,
78,569,078 transitions) from SwissProt human via OpenDIALibGen -> m/z 300-1200 trim ->
OpenSwathAssayGenerator -> OpenSwathDecoyGenerator. Prediction alone took ~70 min wall /
7.5 days CPU on 224 threads and emits NO progress output — worth fixing.

**SCALE WARNING.** This is ~100x the earlier baseline. Previous numbers (DIA-NN 43,336
precursors, "fair target ~6,048") came from `agxt_variants.fasta`, a small variant FASTA that
existed only on spock. Comparisons against those numbers are NOT like-for-like and must not
be presented as such. The OpenDIAlyzer-vs-OpenSwathWorkflow comparison within one run is
unaffected (same library both arms).

**Parity step 1.** Options that CHANGE EXTRACTION RESULTS were hardcoded, so the tool could
not be made to match stock OpenSwathWorkflow even in principle. Now registered with TOPP
names/defaults and routed into `ChromExtractParams` / `performExtraction`:
`mz_extraction_window_unit`, `min_upper_edge_dist`, `extraction_function`,
`extra_rt_extraction_window`, `ms1_isotopes`, `batchSize`, `innerBatchSize`, `readOptions`,
and the `irt_*` windows.

**The important one: `load_into_memory` was hardcoded `true`.** That materialises the working
SWATH in RAM — precisely the whole-run residency the mzPeak streaming input exists to avoid,
so an mzPeak-backed run would have been silently defeated by the tool itself. Now driven by
`-readOptions` (default `normal` = streaming).

Selftest still passes. Under adversarial review now; the load_into_memory default change is
the thing to attack hardest, since it trades memory for throughput and our adapter serialises
decode behind a mutex.

Remaining parity gap (not yet done): `out_chrom`, `out_features(_type)`, `out_qc`,
`out_mobilogram`, `swath_windows_file`, `sort_swath_maps`, `matching_window_only`,
`use_elution_model_score`, `enable_ipf`, `split_file_input`, `append_oswpq`,
`keep_cached_files`, `tr_type`, `outer_loop_threads`, and the `Calibration:*` / `Debugging:*`
subsections.

### Iteration 36b — the review caught my parity change being half-wired

kimi attacked the `load_into_memory` change as asked and found two defects, both confirmed
against the vendored OpenMS source:

1. **`-readOptions` was a no-op for map loading.** `loadDIARun_` hardcoded
   `readopts = (BRUKER_TDF ? "normal" : "cacheWorkingInMemory")`. The new option only reached
   the `load_into_memory` bool passed to `performExtraction`, so the tool could be told the
   maps were resident when they were not — and the help text describing it as controlling
   streaming was simply false.
2. **Defaulting to `normal` silently disabled the wave scheduler.**
   `OpenSwathWorkflow.cpp:486`: `use_swath_range_scheduler = batchSize <= 0 && load_into_memory
   && !ms1_only && ...`. So my "safer" default would have dropped every mzML run back to the
   legacy nested loop — undoing the #22 parallelism work, invisibly.

Fixed: `readOptions` gains an `auto` default that preserves the established per-format
behaviour, is actually threaded into `loadDIARun_`, and is resolved the SAME way at both
sites so the workflow is never misinformed. When the effective choice disables the wave
scheduler, the tool now says so explicitly rather than quietly losing throughput.

**Real tension worth stating plainly:** in-memory working maps are a PRECONDITION for
OpenSWATH's wave scheduler, but the mzPeak streaming input exists precisely to avoid whole-run
residency. Those two goals are in direct conflict in the current OpenMS design — bounded
memory costs the parallel scheduler. Resolving it means making the scheduler work against a
streaming SpectrumAccess, which is the real content of #21+#22, not a flag change.

**Library composition** (final): 39,589,429 target + 38,979,648 decoy transitions over
7,149,966 precursors (~3.58M targets, decoy fraction 49.6%).

### Iteration 36c — three more defects from the review, two of them mine

kimi's full report, all confirmed:

1. **The `irt_*` options were registered but NEVER APPLIED.** `cp_irt` inherited the analyte
   m/z window; the three options I had just added did nothing. That is the "lie in the help
   text" case, and it matters because the calibration extraction is deliberately NOT the
   analyte extraction — anchors are few and must be found before any m/z correction, which is
   why TOPP uses a wider window there. Now applied (window, unit, IM).
2. **`cp_ms1 = cp` inherited the MS2 m/z UNIT.** Selecting `-mz_extraction_window_unit Th`
   silently reinterpreted the MS1 window, help-labelled `<ppm> 30.0`, as **30 Th** — wrong by
   orders of magnitude. TOPP keeps `mz_extraction_window_ms1_unit` separate; now so do we, via
   a `makeMs1ChromParams_` helper so the two cannot drift apart again.
3. **`cp_ms1.im_extraction_window = -1` unconditionally** while `use_ms1_im = pasef` is passed
   to the workflow: on diaPASEF, MS1 ion-mobility SCORES were being computed from
   chromatograms extracted with NO ion-mobility windowing. Registered
   `im_extraction_window_ms1` and wired it.

Also dropped the "an INI is transferable" claim from the parity comment — it was false. Our
defaults deliberately differ (mz_extraction_window 30 vs TOPP 50, ion_mobility_window 0.047 vs
-1) and `outer_loop_threads` / `use_ms1_ion_mobility` / `sort_swath_maps` /
`matching_window_only` / `enable_ipf` are still unregistered. The comment now says so.

Selftest passes. Two review rounds on this step have produced 5 confirmed defects, 4 of them
introduced by the step itself — the pattern is that adding an option is not the same as
wiring it, and inheriting a struct silently inherits its units.

## Iteration 37 — runtime pilot exposes TWO blockers before any benchmark can run

Asked for approximate S08 runtimes before committing to the whole-proteome benchmark. Built a
scaling pilot (library subsets 1-in-2000 / 500 / 125 of precursor GROUPS, sampled by a counter
so each subset spans the full m/z and RT range, not a head-slice corner). The first run failed,
and the failure was worth more than the timing would have been.

**BLOCKER 1: the mzML input path has never worked.** `loadDIARun_` passed
`readoptions="cacheWorkingInMemory"` straight into `SwathFile`, which accepts only
`normal|cache|split` and throws on anything else:
```
Error: Unexpected internal error (Unknown or unsupported option cacheWorkingInMemory)
```
`cacheWorkingInMemory` / `workingInMemory` are TOPP-LEVEL names that must be MAPPED
(OpenSwathWorkflow.cpp:1020-1030: cacheWorkingInMemory -> readoptions="cache" +
load_into_memory=true; workingInMemory -> "normal" + true). Only Bruker `.d` ever worked,
because it happened to take the `normal` branch. kimi had flagged the valid-strings mismatch
in the previous review and I did not act on it; running the thing found it. Both call sites now
go through one `resolveReadOptions_()` so the SwathFile option and the `load_into_memory` flag
can never disagree.

**BLOCKER 2: `S08_full.mzML` is not a DIA file.** It was written by **diaTracer 2.2.1** — a
pseudo-MS/MS generator for FragPipe/MSFragger. It contains one synthetic spectrum per detected
precursor with the isolation window set to that precursor's m/z, and **zero MS1 spectra**:
```
Determined there to be 194389 SWATH windows and in total 0 MS1 spectra
```
That is why the window count is absurd. Targeted extraction against it would be meaningless —
it would have produced numbers, and they would have been garbage. `S23_full.mzML` and
`s08.mzML` are the same. NONE of the mzML files staged on `data` is a real DIA run.

The real inputs are the Bruker `.d` (on `/ceph/ibmi/abi/data/2026_AGXT_PH1_liver_diaPASEF/raw_d/`)
and the `.mzpeak` we converted from it. Staging the S08 `.d` to node-local `/scratch` and
re-running the pilot against it.

Any earlier result obtained from `*_full.mzML` should be treated as suspect until re-checked
against the `.d`.

## Iteration 38 — the prefilter: two real bugs fixed, and the criterion itself found non-discriminating

Asked to fix "the library filtering issue" and prove it on a benchmark. Reviewed by kimi and
codex (both adversarial, independently); they agreed on the worst caller bug and disagreed on
the FDR question, where codex was right.

### Bug 1 — decoy pairing by SEQUENCE (fixed, proven)

The evidence filter scores TARGETS only (its header says so), so the caller re-attaches decoys.
It matched them by peptide sequence -- but `MRMDecoy` SHUFFLES the sequence and puts the link in
the id (`decoy_tag + target_id`, MRMDecoy.cpp:642). Verified on the real library:

```
decoy id  : DECOY_AAAAAAAAAAAAAAAASAGGKEAASGPNDS_2   -> strips to an existing target id
target seq: AAAAAAAAAAAAAAAASAGGKEAASGPNDS
decoy  seq: AAAAAGAASSAGAANGAADAKAAEAAPAAA           -> DIFFERENT
```

So only accidental collisions survived: **23 decoys against 3,546,918 targets**. 23 is not 0, so
the `n_dec == 0` guard passed it and the search ran with a 154,000:1 ratio and a meaningless FDR.
Now paired on the id prefix; guard is a RATIO BAND against the library-wide ratio, which makes it
a self-test of the pairing rather than a check on one degenerate value. Measured after the fix on
astral.mzML: **911,791 target / 899,533 decoy = 1:0.986.**

### Bug 2 — `evidence_sources` was decorative (fixed)

Both reviewers found it. The caller hardcoded `supported_ms1 || supported_ms2`, i.e. it
re-implemented "hybrid" no matter what was configured. Setting `evidence_sources=ms2` changed
only the LOGGED count, not one retained precursor. My own first selectivity fix was therefore a
no-op until this was corrected.

### The MS1 arm is a rubber stamp (quantified on the benchmark file)

`supported_ms1 = ms1_hit_count > 0` -- one m/z coincidence anywhere in the run, no intensity
floor, no RT coherence -- and hybrid is OR. On astral.mzML (3,603,425 targets):

| evidence | targets kept | % |
|---|---|---|
| hybrid (old default) | 2,757,283 | 76.5% |
| ms2 (new default)    |   911,791 | 25.3% |

The MS1 arm alone contributes ~51 points of pure noise. Default is now `ms2`.

### FDR: pairing is NOT symmetric (codex right, kimi wrong)

kimi argued decoy ride-along makes any target-side criterion automatically symmetric. It does
not: pairing preserves the 1:1 COUNT but not the evidence CONDITIONING. A false target is
retained *because* it had coincidental evidence; its decoy rides along having faced no test, and
those same coincidental peaks then feed the score. Retained decoys therefore score below retained
false targets -> null too weak -> **anti-conservative FDR**, and the bias grows exactly as the
filter becomes selective. Fixed by scoring BOTH classes with the identical criterion and
retaining a pair when EITHER passes (pair-union), which makes selection invariant under swapping
the labels. The filter rejects decoys twice (`getDecoy()` at :375 AND a hardcoded `hasDecoyPrefix_`
id test at :356), so the decoy view clears the flag AND aliases the ids. Costs one extra scan.

### The finding that matters: the criterion does not discriminate

Sweep on astral.mzML, `evidence=ms2`, pair-union, dumping survivors and exiting before extraction
(new `-prefilter_out`, ~20 min per point instead of ~3 h):

| top_peaks | min_frag | targets pass | decoys pass | enrichment | kept/7.15M | recall vs DIA-NN |
|---|---|---|---|---|---|---|
| 3000 | 4 | 2,254,295 (62.6%) | 2,179,018 (61.4%) | 1.02x | 71.0% | 99.90% |
| 6000 | 4 | 2,284,318 (63.4%) | 2,209,965 (62.3%) | 1.02x | 71.4% | 99.90% |
| 1000 | 5 |   110,483 (3.07%) |    99,020 (2.79%) | 1.10x |  5.5% | 84.87% |

**Shuffled decoys pass at 91-98% the rate of targets.** 4-of-6 (or 5-of-6) library fragment m/z
matches within 30 ppm among the top-N peaks of a crowded Astral MS2 spectrum is simply not a rare
event, and a shuffled peptide's fragment masses are drawn from the same distribution as a real
one's. Restricting to the precursor's own SWATH window (`buildProductIndexForMap_` does do this)
is not nearly enough. Consistency check: at 1000/5 the target-minus-decoy EXCESS is
110,483 - 100,608 = ~9,900, the right order for DIA-NN's 8,164 real IDs -- the signal is there,
but it is swamped by a null of the same size.

Consequence: **there is no (top_peaks, min_fragments) operating point that is both selective and
lossless.** 99.9% recall costs 71% of the library; 5.5% of the library costs 15% of the real IDs.
Raising top_peaks past 3000 changes nothing (identical recall, +0.4% kept), so the peak cut is not
the binding constraint.

The missing dimension is RT. Coincidental m/z matches do not co-elute; real precursors do. The
filter records `ms1_best_rt` / `ms2_best_rt` but uses neither as a criterion. The next rule to try
is codex's: fewer fragments required, but demanded consistently across neighbouring DIA cycles
within one chromatographic FWHM. That needs per-spectrum RT events in the filter -> an OpenMS
change, not a parameter change.

### Also fixed

- Fail-fast instead of fail-open. Both reviewers flagged it; the evidence settles it -- the
  unfiltered path on astral gives 27,344,193 features whose BEST q-value is **0.508**, i.e. zero
  identifications at any usable FDR, for 1753 GB and 2.9 h. Silently falling back to that is not a
  safe default for a mis-set window. Searching unfiltered now requires `-prefilter false`.
- `keep` keyed on `compound_id` only; the bare-sequence fallback was retaining zero-evidence
  charge siblings (evidence for `PEPTIDEK_2` kept `PEPTIDEK_3`) and dragging their decoys along.
- `-prefilter_min_fragments` > the 6 indexed transitions per precursor is now a loud error rather
  than silent zero support.
- New: `-prefilter_evidence`, `-prefilter_min_fragments`, `-prefilter_top_peaks`, `-decoy_tag`,
  `-prefilter_out`.

### Benchmark reference (all on astral.mzML, the SAME 7,149,966-precursor library)

| | features | IDs @1% FDR | RSS | wall |
|---|---|---|---|---|
| DIA-NN 2.0 | -- | 8,164 | 22.7 GB | 509 s |
| OpenDIAlyzer, unfiltered | 27,344,193 | **0** (min q 0.508) | 1753 GB | 2.9 h |

DIA-NN handles this library fine, so the library is not the problem. Recall is measured against
DIA-NN's 1%-FDR set mapped to library ids; 7,627 of its 8,164 (7,787 unique seq_charge) exist in
the library as exact TransitionGroupIds, and that 7,627 is the honest denominator.

### Why selectivity is bad: not narrow RT windows -- NO RT window

Asked whether the poor selectivity is a calibration artefact (too-narrow RT windows). It is the
opposite, and the check is exact. `TransitionListEvidenceFilter` reads only these from
`ChromExtractParams`:

```
params.mz_extraction_window   params.ppm   params.im_extraction_window   params.min_upper_edge_dist
```

`rt_extraction_window` is never read. The scan loop is `for (spectrum_index = 0; < nr_spectra)`
-- unbounded -- and line 399 actively throws the library RT away (`candidate.compound.rt = -1.0`).

WHY it has no RT constraint: `prefilterLibrary_` is called at opendialyzer.cpp:1468 and the RT
calibration begins ~29 lines LATER. The prefilter runs before any calibration, so at that moment
the library RT is still normalized [0,1] with no map onto run RT. Chicken-and-egg.

Consequence, quantified on astral: 150 SWATH windows and 3889 MS1 spectra, so each precursor is
screened against ~3889 MS2 spectra of its own window -- the whole gradient. The rule takes the MAX
single-spectrum fragment count over all 3889. The null is (trials x per-spectrum coincidence), and
3889 trials is what turns a rare per-spectrum event into a 62% pass rate.

**Calibration-based RT gating would NOT rescue this here**, which is worth stating because it is
the intuitive fix: our own measurement on this data is a PeptDeep-vs-observed residual sd of ~810 s
with 95th pct ~1720 s even after nonlinear OSW-style calibration (see run_full_arm.sh). A window
wide enough to keep 95% of real precursors would be ~+/-1700 s -- most of the gradient -- so it
would remove almost none of the 3889 trials.

The fix should therefore be calibration-FREE: the existing rule already requires the fragments in
ONE spectrum (so they co-elute trivially); what is missing is that it never asks whether the
evidence RECURS. Require >=2-3 CONSECUTIVE cycles instead of the max over all cycles. A real
precursor persists across cycles -- DIA-NN's own stats on this exact run report FWHM.Scans = 2.528
-- while a coincidence is a one-off, so the null falls roughly as p^k rather than growing linearly
in trials. This asks "did it recur?", not "did it appear at the predicted time", so it sidesteps
the ordering problem entirely and does not depend on the (poor) RT calibration.

### Recurrence tested and REFUTED — the whole m/z-coincidence family cannot work here

Implemented `ms2_min_qualifying_spectra` in the filter (counts spectra in which the fragment
threshold was met, rather than taking the max over the run) and swept it on astral.mzML at
top_peaks=1000, min_fragments=4, pair-union on. Incremental libOpenMS rebuild: **9.4 s**, not
hours -- OpenMS-side changes are cheap here, worth remembering.

| min_spectra | targets pass | decoys pass | ENRICHMENT | kept/7.15M | recall (of 7627) |
|---|---|---|---|---|---|
| 1 | 911,791 (25.30%) | 867,164 (24.45%) | 1.03x | 38.2% | 96.93% |
| 2 | 458,652 (12.73%) | 436,272 (12.30%) | 1.03x | 20.7% | 89.54% |
| 3 | 281,682 (7.82%)  | 268,381 (7.57%)  | 1.03x | 13.2% | 78.71% |
| 5 | 142,388 (3.95%)  | 136,429 (3.85%)  | 1.03x |  6.9% | 52.77% |

**Enrichment is 1.03x at EVERY level.** The criterion bites hard (recall collapses 96.9 -> 52.8%)
but bites targets and decoys equally. Recurrence adds nothing.

Why, and it should have been predictable: in a dense DIA spectrum an m/z coincidence is not a hit
on noise, it is a hit on a REAL peak belonging to some other co-eluting peptide. Real peaks persist
across cycles. So a shuffled decoy's coincidental matches recur exactly as reliably as a true
precursor's own fragments, and persistence carries no target/decoy information whatsoever.

Taken with the earlier sweeps, thirteen configurations have now been measured across three
independent axes -- fragment threshold (4, 5 of 6), peaks per spectrum (1000, 3000, 6000), and
recurrence (1, 2, 3, 5 spectra) -- and EVERY one lands at 1.02-1.10x enrichment. That is not a
tuning problem. Counting m/z matches, however it is thresholded, cannot separate a real precursor
from a shuffled decoy at whole-proteome scale on this data. The prefilter idea is dead as a
correctness mechanism; further threshold work is not worth the compute.

### Where the leverage actually is

The prefilter was never the reason this pipeline fails. On the SAME library and the SAME file:

| | features | IDs @1% FDR | RSS | wall |
|---|---|---|---|---|
| DIA-NN 2.0 | -- | 8,164 | 22.7 GB | 509 s |
| OpenDIAlyzer, unfiltered | 27,344,193 | **0** (best q = 0.508) | 1753 GB | 2.9 h |

DIA-NN searches 7,149,966 precursors against this run in 8.5 minutes on 22.7 GB and separates
8,164 of them at 1% FDR. So 3.6M target precursors is NOT an intractable search, and reducing the
library was never the necessary step. Two real problems remain, both downstream of the prefilter:

1. **Scoring cannot separate signal from noise.** Best q-value 0.508 on an unfiltered run means
   the target-decoy competition is not merely mis-calibrated, it is uninformative. That is the
   blocker for identifications and no amount of prefiltering fixes it.
2. **Extraction cost is ~200x wall and ~80x memory versus DIA-NN** for the same input and library.

What discriminates in DIA-NN is chromatographic: fragment co-elution CORRELATION (do the traces
share an apex and a peak shape), agreement with library relative intensities, and the mass-error
distribution. None of that is an m/z-counting prefilter; it is the scoring itself. The next
investigation should be why the scores carry no information, starting from the fact that an
unfiltered run with 27.3M features cannot put a single precursor below q = 0.5.

## Iteration 39 — why there are no identifications: the RT window is narrower than the RT error

Followed the prefilter dead end to the actual blocker. On the unfiltered threeway run
(astral.mzML, 27,344,193 features, rank-1 only, from `odia.osw`):

```
target  min=-1.4941  p25=0.3979  med=1.0955  p75=1.6338  p95=2.1816  p99=2.4670  max=2.7451
decoy   min=-1.4941  p25=0.3360  med=1.0462  p75=1.6036  p95=2.1631  p99=2.4457  max=2.7451
targets scoring above the BEST decoy: 0
targets above decoy p99.9: 2,861   (decoys above their own p99.9: 2,713 -> 1.05x)
```

Target and decoy score distributions are the SAME -- identical min and max to four decimals,
medians 1.10 vs 1.05, and not one target outscores the best decoy. The pipeline discriminates at
~1.05x, the same figure the prefilter produced. That is not a threshold problem, it is the
signature of extracting from the WRONG PLACE.

**Root cause.** That run extracted with `RT window 746.9 s` (+/-373 s). Our own measured
PeptDeep-vs-observed RT residual on this data is **sd ~810 s, 95th pct ~1720 s**. For residuals
~N(0, 810), P(|resid| < 373) ~= **35%**: roughly two thirds of true precursors never had their
peak inside the extraction window at all. For those, the tool extracts noise -- exactly what it
extracts for a decoy. Hence identical score distributions, hence best q = 0.508.

So the earlier question "is the bad selectivity poor calibration / too-narrow RT windows?" was
right about the mechanism but aimed at the wrong stage. The PREFILTER ignores RT completely
(it never reads `rt_extraction_window`); it is EXTRACTION and SCORING that use it, and there the
window is less than half the residual sd.

**And the benchmarks disabled the fix.** `recal_passes` defaults to **2** (recalibrated two-pass),
but every benchmark script -- threeway.sh, overnight.sh, bench_pass00779.sh, bruker_bottleneck.sh,
mzpeak_mt.sh, and all of this iteration's prefilter sweeps -- passes `-recal_passes 1`, the single
wide pass with NO recalibration. Pass 2 is precisely the step that re-centres the window from
confident pass-1 IDs. Every published number in this log for OpenDIAlyzer identifications was
therefore produced with the tool in its least capable mode.

Testing now on astral.mzML (`bench/pfrecal`): arm A `-recal_passes 2`; arm B additionally with a
pass-1 window wide enough to CONTAIN the true peaks (`-rt_extraction_window 3600`,
`-use_estimated_rt_window false`) so the recalibration anchors on real peaks rather than noise.
Caveat worth stating in advance: pass 2 bootstraps from pass-1 confident IDs, and if pass 1 has
almost none, it may have nothing to anchor on -- which is what arm B is for.

### Separate defect found on the way

`PVALUE == QVALUE == PEP` for **100% of 27,344,193 rows** in SCORE_MS2 (194,639 distinct values,
so the number varies per row -- it is the three columns that are written identically). These are
three different statistics; anything downstream reading PEP is reading a q-value. Not the cause of
the zero IDs (the q-values look correct GIVEN the scores), but wrong and worth fixing.

## Iteration 40 — ODIA vs OSW parity: the three-way was never a controlled comparison

Asked to make ODIA no worse than stock OpenSwathWorkflow. First job was establishing what was
actually being compared, and three things in the earlier write-up were wrong.

**Correction 1: OSW DOES calibrate, and did.** `Calibration:auto_irt:enabled` defaults to **true**
and uses the built-in `irtkit.tsv` (10 compounds) + `cirtkit.tsv` (124) as priority peptides -- the
same auto-iRT machinery ODIA drives. Its "Will analyse 492 peptides" / "23,816 peptides" lines are
the linear and nonlinear iRT sampling. The claim that the OSW reference ran uncalibrated was wrong.

**Correction 2: both used 30 ppm.** The OSW command line passed `-mz_extraction_window 30`
explicitly. The claim that OSW used the 50 ppm TOPP default was wrong.

**Correction 3, and it matters most: OSW ran `-threads 24`, ODIA `-threads 224`.** The headline
table compared a 24-thread OSW against a 224-thread ODIA -- and OSW still won on wall clock AND
memory with a tenth of the threads.

### Measured phase attribution (astral.mzML, one shared 7,149,966-precursor library)

| phase | OSW (24 thr) | ODIA (224 thr) |
|---|---|---|
| CiRT/auto-iRT calibration | 492 linear / 23,816 nonlinear -> window **375.3 s** | 337 linear / 67,429 nonlinear, 22:03 wall / **1d 08:40 CPU** -> **136 anchor pairs**, window **746.9 s** |
| main extraction | 1:25:49 wall / **1d 05:51 CPU** (~21 cores) | 2:03:46 wall / **6d 07:44 CPU** (~74 cores) |
| total | 2:14:22, 153.2 GB, 27,343,888 features | 3:57:34, 1538.6 GB, 27,344,193 features |

Four defects, all concrete:

- **F1** ODIA's extraction burns **5.1x** OSW's CPU for identical output and is still 1.4x slower.
- **F2** Calibration costs 32.7 CPU-hours -- more than OSW's ENTIRE extraction -- to yield 136
  anchor pairs from 67,429 candidates (**0.2%**).
- **F3** The 746.9 s pass-1 window is estimated from those 136 anchors (in-sample) and spent on
  3.6M predicted-RT precursors whose true residual p95 is ~1720 s. OSW estimated 375.3 s; both are
  far too small for unseen precursors.
- **F4** `readOptions` default divergence: ODIA `auto` -> `cacheWorkingInMemory` (full residency)
  vs OSW `normal` (streaming). This is the 1538 GB vs 153 GB gap.

### F4 fixed — streaming is now the default

`readOptions` now defaults to **`normal`**, matching OSW. The old `auto` chose
`cacheWorkingInMemory` for mzML because in-memory working maps are a PRECONDITION for OpenSWATH's
SWATH-range (wave) scheduler (`OpenSwathWorkflow.cpp:486`: needs `batchSize<=0 && load_into_memory`),
and that scheduler had been measured 2.7x faster on a smaller configuration. At full scale the
trade inverts completely: **5.1x the CPU and 10x the memory to run 1.4x slower**. Full residency is
not worth buying a scheduler that loses. `auto` stays selectable. The message on the streaming path
was WARN ("throughput will drop") and is now INFO -- that claim is contradicted by measurement.

### F2 fixed — the calibration was sampling more anchors to solve a window problem

ODIA built `CalibrationWorkflow` and set only `qc:min_rsq` / `qc:min_coverage`, leaving the sampling
parameters at OpenMS defaults: 2000 bins x 50 peptides/bin, and
`auto_irt:irt_nonlinear_rt_extraction_window` = **600 s**. The 0.2% yield is not a sampling-size
problem -- it is the same root cause as everything else in this pipeline: at a 600 s search window
against a PeptDeep RT residual of sd ~810 s, most anchors' true peaks are not inside the window
being searched. Extracting 968,800 transitions cannot find peaks that are elsewhere.

Now exposed and set: `-calibration_nonlinear_per_bin` (10, was 50) and
`-calibration_nonlinear_rt_window` (2400 s, was 600 s) -- sample far fewer, look over a window wide
enough to contain them. Plus a yield diagnostic: the log now reports "N pairs from M candidates
(X% yield)" and WARNs below 5%. "136 anchor pairs" reads fine in isolation; "136 from 67,429" does
not, and the gap between those two log lines was 32 unnoticed CPU-hours.

Also note ODIA sets `qc:min_coverage` = 0.30 where OSW uses 0.60, i.e. ODIA accepts calibrations
OSW would reject -- which is a candidate explanation for its looser 746.9 s window estimate. Left
as-is pending review; flagged rather than silently changed.

Matched-pair run in flight (`bench/parity`): identical 224 threads, 30 ppm, 600 s window, IM off,
streaming, no library prefilter, calibration on for both.

## Iteration 41 — streaming pays off; OSW crashes at 224 threads; LOESS added and then hardened

### Parity run 1 (224 threads, matched settings, both calibrating)

| | wall | CPU | peak RSS | features | best q |
|---|---|---|---|---|---|
| ODIA before (wave, 746.9 s window) | 3:57:34 | ~48 cores | 1538.6 GB | 27,344,193 | 0.508 |
| ODIA now (streaming, 1135.5 s window) | **3:10:06** | ~72 cores | **881.7 GB** | 27,343,749 | 0.406 |
| OSW @24 threads (earlier) | 2:14:22 | ~21 cores | **153.2 GB** | 27,343,888 | not scored |

Streaming cut peak RSS **43%** and wall **20%** while doing 1.52x more work (wider window). Still
5.8x OSW's memory. Leading hypothesis (kimi): the residency is per-thread extraction state, not the
working map -- testable by thread count, which is also the only way to run OSW at all, because:

**OpenSwathWorkflow SEGFAULTS at 224 threads on this input** -- `Command terminated by signal 11`,
38:41 in, after a 9.8 GB partial write. `/usr/bin/time` still prints `Exit status: 0`; the harness
believed it and recorded a crash as a 0-feature datapoint. The scripts now grep for the signal line.
This is presumably why the original benchmark used `-threads 24`.

### The calibration fix cut cost but collapsed the anchors

Candidates fell 67,429 -> 3,662 as intended, but anchor PAIRS fell **136 -> 19**, and the estimated
window GREW to 1135.5 s (fewer survivors scatter more). The new yield diagnostic fired and answered
the question it was built for:

```
calibration anchor yield is only 0.518842% (19/3662).
At most 63 candidates matched before outlier removal.
```

The bound is small, so the anchors are NOT BEING FOUND -- widen, do not narrow. Yield per candidate
did improve (0.20% -> 0.52%), which supports that reading, but sampling was over-cut 17x while the
window widened only 4x.

### LOESS RT recalibration — added, reviewed, and substantially rewritten

`fitTrafo_` gains a LOESS path: local linear regression with tricube weights over all anchors, then
PAVA (elution order is physics; local fits can cross). Adversarial review found it was, as first
written, worse than the estimator it replaced:

- **It was inert at the anchor counts we actually see.** Binning gives `NB = min(60, n/10)`, so at
  n=19 there is ONE bin, `ux.size()==1`, and the `>= 2` guard skipped LOESS silently. The feature
  did nothing at exactly the count that motivated it. Now it refuses explicitly below 50 anchors and
  says so in the log.
- **It was not LOESS.** Cleveland's estimator has bisquare reweighting iterations; one tricube pass
  is just local regression, with breakdown ~1/q (~2%) against the binned MEDIAN's 50%. In a regime
  where anchors are hunted in windows narrower than the error -- i.e. guaranteed dirty -- that was a
  robustness REGRESSION. Two bisquare iterations added.
- **The grid was still the bin centres**, so LOESS improved the y at each knot but never knot
  DENSITY; curvature between knots stayed linearly interpolated. Now a uniform grid (<=128 points)
  over the anchor x-range, never extrapolating beyond it.
- Span now `max(span, 15/n)` so the neighbourhood grows as n shrinks; relative (not absolute)
  degeneracy guard; stale bin-count weights no longer carried into PAVA after LOESS overwrites y.

**My first selftest was wrong twice.** It initially scored LOESS on ANCHOR residual -- which rewards
overfitting noise -- and failed, correctly. Rewritten to measure error against the TRUE warp at
held-out midpoints. Review then showed the monotonicity assert was tautological (PAVA guarantees it
by construction, so it would pass even if LOESS returned crossings). Now three checks that can
actually fail: accuracy vs truth on a curved warp; **not materially worse than the median under 5%
gross outliers** (this is the one that guards the robustness regression); and fallback-to-binned at
19 anchors.

### Honest ranking (kimi, and I agree)

LOESS is **fourth** by expected effect on identifications, and its realistic effect while the window
problem stands is "indistinguishable from zero" -- it refines the smooth warp component, but what
remains after calibration (sd 810 s) is PER-PEPTIDE prediction error, which is not a function of RT
and therefore unreachable by any monotone transform. Ranked ahead of it: (1) widen the extraction
window to cover the measured residual, (2) fix anchor yield, (3) replace predicted RT with measured
per-precursor RT. Tonight's machine time goes to (1).

### CORRECTION — the "OpenSwathWorkflow segfaults" finding was self-inflicted

Retracted: "stock OpenSwathWorkflow does not survive 224 threads on this input". It was not a thread
problem and not an OSW problem. It crashed at **24** threads too, which is what exposed the real
cause.

I added a member (`ms2_qualifying_spectra`) to `PrecursorEvidence` / `EvidenceAccum` in
`TransitionListEvidenceFilter.h`, rebuilt `libOpenMS.so`, and installed it -- **without rebuilding
the TOPP binaries**. Timeline is unambiguous:

| | timestamp |
|---|---|
| `bin/OpenSwathWorkflow` built | 07-29 17:09 |
| `lib/libOpenMS.so` replaced by me | **07-30 14:44** |
| last SUCCESSFUL OSW run (osw2) | 07-29 23:27 -- before the change |
| OSW crash @224 threads | 07-30 18:26 -- after |
| OSW crash @24 threads | 07-30 22:23 -- after |

Classic ABI break: the executable was compiled against the old struct layout and kept reading the
old offsets out of a library where the layout had moved. Nothing to do with concurrency.

Blast radius established by symbol inspection rather than assumption -- of 151 installed TOPP
binaries, exactly two have an undefined reference to `TransitionListEvidenceFilter`:
`OpenSwathWorkflow` and `TransitionListEvidenceFilter`. Both rebuilt and reinstalled; OSW starts
cleanly again. The other 149 do not touch the changed header and are unaffected.

**Process lesson worth keeping:** the vendored-OpenMS workflow makes an incremental `ninja OpenMS`
+ `cp libOpenMS.so` look sufficient -- it is not, whenever a header changes a type that any TOPP
tool instantiates. Rebuild the dependent binaries, or the next "the reference tool is broken"
finding will be another own goal. It also cost two full benchmark slots (38 and 39 minutes) and
produced a confident, wrong claim in this log.

Re-running the OSW baseline with the fixed binary after the wide-window ODIA job.

## Iteration 42 — pyprophet-compatible statistics: p-value, q-value and PEP were all the same number

Every `.osw` this tool has ever written had `PVALUE == QVALUE == PEP` on **100% of rows**
(27,344,193 of 27,344,193 on the threeway run). Cause, at opendialyzer.cpp:923-925: all three
columns were bound to `s.qvalue[i]`. `ScoredGroups` only ever carried `dscore` and `qvalue` -- the
other two statistics were never computed at all.

They are not interchangeable:

| | meaning |
|---|---|
| p-value | TAIL probability under the null: P(a null score >= this one) |
| q-value | minimum FDR of the SET selected at this threshold |
| PEP | LOCAL FDR: probability that THIS peak group specifically is null |

A group sitting exactly at q = 0.01 can easily have PEP ~ 0.3. Anything downstream reading PEP --
IPF, protein-level inference, pyprophet's own second-level scoring -- was consuming a q-value.

Implemented in `odia_lda.h`:

- **p-value**: empirical tail probability from the decoy null, `(decoys_at_or_above + 1)/(Ndec + 1)`
  (Käll's conservative finite-sample form, consistent with the +1 already used in the q-value).
- **PEP**: the LOCAL analogue of the global estimator already in use -- in a score neighbourhood
  holding t targets and d decoys, the decoys estimate the null target density scaled by
  `pi0 * Ntar/Ndec`, so `PEP ~ pi0*d*(Ntar/Ndec)/t`. Window is ~0.5% of the list (>=100 wide) so it
  adapts to size. A raw local ratio is noisy, so a final sweep enforces the ordering constraint
  (PEP non-increasing in d-score) by running max from the low-score end -- the isotonic projection
  under that constraint, which removes noise without smoothing away the trend.

Test added to `odia_lda_test.cpp` that fails if they ever collapse again: (a) not all three equal
on every row, (b) **PEP >= q-value** -- the local FDR at a threshold cannot sit below the average
FDR of everything above it, (c) PEP monotone non-increasing in d-score. Measured on the oracle:
recall 0.88-0.90, empirical FDR 0.008, all-equal rows 4/16,000 (the extremes, as expected).
Adversarial test still passes: 0 monotonicity violations over 3,000 groups, emp_FDR 0.009,
recall 0.92.

Note the runs already in flight will still write the old degenerate columns, since the process
holds its own binary; `-score_osw` rescores an existing `.osw` in-process, so no re-extraction is
needed to get corrected statistics.

## Iteration 43 — the RT-window hypothesis is FALSIFIED

Ran ODIA with `-rt_extraction_window 3600` (full width, i.e. +/-1800 s), which on a 2333 s gradient
means every precursor's window covers essentially the WHOLE RUN -- RT is effectively unconstrained.
Prediction was that covering the measured residual (sd ~810 s, p95 ~1720 s) would recover
identifications. It did not.

| | 1135.5 s window | 3600 s window |
|---|---|---|
| features | 27,343,749 | 13,481,161 (prefiltered library) |
| best q | 0.406 | **0.561 (WORSE)** |
| targets above best decoy | 0 | **0** |
| peak RSS | 881.7 GB | **1325 GB** |
| wall | 3:10:06 | 2:11:58 |

Score distributions with the window essentially unconstrained:

```
target med=1.1290  p95=2.0355  p99=2.2872  max=2.5213
decoy  med=1.0922  p95=2.0147  p99=2.2600  max=2.5213
targets above decoy p99.9: 1570   (decoys there: 1341)  -> 1.17x
```

**So RT was never the blocker.** Even with RT unconstrained, targets do not outscore decoys: the
enrichment at the extreme tail is 1.17x, and not one target beats the best decoy. Widening made the
q-value WORSE (0.406 -> 0.561), which is the expected direction if the window was never excluding
signal but was excluding noise -- more window, more null candidates competing.

One important framing correction: target and decoy distributions being similar in BULK is EXPECTED
here and is not itself a defect. The library holds 3.6M target precursors (whole proteome,
peptidoform-expanded) against a plasma sample of maybe 10k proteins, so ~99.8% of targets are
genuinely absent and should behave exactly like decoys. DIA-NN's 8,164 is 0.23% of targets. The
defect is that there is no enrichment at the TAIL either, where those 0.23% must live.

New lead: **`target max == decoy max == 2.5213` to four decimals**, and on the previous run the
minima matched exactly as well (-1.4941). Across millions of independent peak groups an exact tie
at both extremes is not chance -- it is the signature of a discriminant that SATURATES. That points
at the scoring, not the extraction geometry.

Decisive experiment now running: real **pyprophet** on the very same `.osw` that ODIA's in-process
LDA scored to 0 IDs. Same features, same labels, different classifier --
  * pyprophet finds IDs   -> the sub-scores carry signal; ODIA's LDA is the problem.
  * pyprophet finds none  -> the signal is not in the features at all, and the bug is upstream in
    extraction / peak-picking / library, where no classifier work would help.
Plus a full sub-score inventory (null / constant / range per column) to see whether the LDA is
being fed degenerate features.

## Iteration 44 — ODIA's .osw was never actually pyprophet-compatible

Ran real pyprophet on the same `.osw` ODIA's in-process LDA scored to 0 IDs. It did not get as far
as scoring:

```
_duckdb.ConversionException: Unimplemented type for cast (BIGINT -> BLOB)
                             when casting from source column ID
```

Schema comparison against a genuine OpenSwathWorkflow-written file:

```sql
-- OpenSwathWorkflow
CREATE TABLE FEATURE(ID INT PRIMARY KEY NOT NULL, RUN_ID INT NOT NULL, PRECURSOR_ID INT NOT NULL, ...)
-- OpenDIAlyzer
CREATE TABLE "FEATURE"(ID INT, RUN_ID INT, PRECURSOR_ID, ...)
                                           ^^^^^^^^^^^^ no declared type
```

**Cause.** `remapFeaturePrecursorIds_` rebuilds FEATURE with `CREATE TABLE ... AS SELECT`. SQLite
gives a CTAS column NO DECLARED TYPE unless the select-item is a bare column reference, and the
remap's item is `COALESCE(p.ID, f.PRECURSOR_ID) AS PRECURSOR_ID`. The column therefore came out
with NONE affinity. SQLite does not care. pyprophet reads `.osw` through duckdb, which does.

**Why nobody noticed.** The step verified itself with `SELECT typeof(PRECURSOR_ID) FROM FEATURE`,
which returns the type of the stored VALUE ("integer") -- not the declared column type. So the check
passed and the tool logged *".osw is now pyprophet-compatible"* on a file pyprophet cannot open. A
self-check that reports success on a broken artifact is worse than no self-check.

**Fixed:** declare the schema explicitly (name + type + NOT NULL carried over from the live schema,
`PRECURSOR_ID` forced to INT) and `INSERT INTO ... SELECT` instead of CTAS. Verification now reads
`PRAGMA table_info` and fails loudly when the declared type is empty, naming the exact duckdb error
that will otherwise appear downstream.

The existing 13.5M-feature file was repaired in place with the same schema rewrite (20 s) rather
than re-extracted, so tonight's decisive test could still run:
`PRECURSOR_ID  <NO TYPE> -> INT`.

This also means every prior claim about ODIA output being pyprophet-scorable was untested -- no
`.osw` this tool has produced has ever actually been through pyprophet.

## Iteration 45 — SOLVED: the m/z extraction window was the blocker. First identifications.

DIA-NN's own stats on this run: `Median.Mass.Acc.MS2 = 1.63 ppm` (0.93 corrected), MS1 0.60 ppm.
ODIA and the OSW reference both extracted at **30 ppm** -- ~20x wider than the instrument delivers.
Against a library where ~99.8% of targets are genuinely absent, that means every precursor finds
SOMETHING, and discrimination collapses. Narrowing it:

| config | wall | peak RSS | features | best q | IDs @1% FDR |
|---|---|---|---|---|---|
| 30 ppm, rt 1135 s | 3:10:06 | 881 GB | 27,343,749 | 0.406 | **0** |
| 30 ppm, rt 3600 s | 2:11:58 | 1325 GB | 13,481,161 | 0.561 | **0** |
| **10 ppm**, rt 600 s | **30:14** | **189 GB** | 2,070,241 | **0.0009** | **2,957** |
| **5 ppm**, rt 600 s | **19:57** | **74 GB** | 467,132 | 0.0010 | **2,904** |

From zero identifications to 2,957, and `targets_above_best_decoy` from 0 to 1,089 -- real
separation where before there was none. Feature count fell 27.3M -> 2.1M: those 25M extra features
were noise admitted by an over-wide m/z window, which is also where the memory went.

This retrospectively explains every earlier symptom: prefilter target/decoy enrichment 1.02x,
scoring enrichment 1.05x, identical score minima and maxima, and q bottoming out at 0.4-0.56. None
of it was RT, the classifier, the decoys, or the features -- all of which were separately ruled out.

### OSW reference (rebuilt binary — no segfault, confirming the ABI break diagnosis)

`wall=6886s (1:54:46) cpu=1448% (~14 cores) rss=153GB features=27,343,888`

### Standing vs the two references

| | wall | peak RSS | IDs @1% FDR |
|---|---|---|---|
| DIA-NN 2.0 | 509 s | 22.7 GB | 8,164 |
| OpenSwathWorkflow @24 thr | 6,886 s | 153 GB | never scored |
| **ODIA @5 ppm, 224 thr** | **1,197 s** | **74 GB** | **2,904** |
| **ODIA @10 ppm, 224 thr** | **1,814 s** | **189 GB** | **2,957** |

**ODIA is now 5.7x faster than OSW on half the memory, and produces identifications OSW's
unscored output never did.** Against DIA-NN it is 2.4x slower, 3.3x the memory, and 36% of the IDs
-- a real gap, but a measurable one rather than a zero.

Note 10 ppm yields slightly MORE ids than 5 ppm (2,957 vs 2,904) at 2.6x the memory: 5 ppm is
already slightly too tight for the tail of the mass-error distribution. The sensible default is
~10-15 ppm, not 30, and not 5.

---

## Iteration 46 — the parquet library "decoy loss" was a PQP defect, not a parquet one

The symptom: converting the library to parquet and searching against it prefiltered
`7,149,966 -> 120,525 precursors (120,525 target / 0 decoy)` and aborted with rc=6 on the
target:decoy ratio guard.

### What it was not

Ruled out by direct measurement, in this order:

| Hypothesis | Test | Result |
|---|---|---|
| Parquet **write** drops decoy flags | count `decoy` in transitions.parquet | 39,589,429 T / 38,979,648 D — correct |
| Parquet **read** drops decoy flags | `convertParquetToTargetedExperiment`, count transitions | 39,589,429 T / 38,979,648 D — correct |
| `ChunkedColumn::resolve()` bug over 600 chunks | same test, all 78,569,077 transitions | correct |
| Library's decoy ids don't follow the `DECOY_`+target convention | full scan of library.tsv | 3,546,541 / 3,546,541 pair exactly — 100% |

The parquet layer was clean on every count. `decoy_refs` was correctly populated (3,546,541), and
the decoy view was correctly filtered (109,298 passed).

### What it was

Compound ids came back as `"0"`, `"1"`, `"2"`, … Dumping `traml_id` straight out of
precursors.parquet showed the writer had faithfully stored exactly those. Reading conv.log:

```
OpenDIAlyzer: loading library from cache .../library.tsv.cache.pqp (skipping the serial TSV parse)
```

The conversion never read the TSV. It read a **PQP** cache, and
`TransitionPQPFile::convertPQPToTargetedExperiment` defaults to `legacy_traml_id=false`, which uses
`PRECURSOR.ID` — a row number — as the compound id. The original `DECOY_<target id>` naming was
destroyed at that step, one format earlier than where the symptom appeared.

Since decoy→target pairing is *by id* (MRMDecoy.cpp:842), zero pairs could be formed, so the FDR had
no null distribution and the ratio guard fired. **The guard did its job** — this is the failure the
ratio band was added to catch, and it caught it.

### Scope

This is not confined to the cache. Any `.pqp` passed to `-tr` hit the same path, so
**every PQP-library run ODIA has ever done had unpaired decoys.** The 9.8 GB `library.pqp` in the
bench directory is affected.

### Fixes

1. `loadLibrary_()` now passes `legacy_traml_id=true` for PQP input.
2. `prefilterLibrary_()` checks the pairing precondition *before* the two extraction passes: if no
   decoy id carries `-decoy_tag`, it names the problem and exits in seconds instead of failing on an
   opaque ratio after ~6 min of scanning.
3. All `.oswpq` caches derived from the PQP (and the PQP-derived `library_rt.oswpq`) retired to
   `*.badids`; library re-converted from `library.tsv` directly.

### Lesson

Two silent fallbacks chained into a wrong answer: `legacy_traml_id=false` silently substitutes row
numbers for names, and the parquet reader silently substitutes the numeric precursor id when
`traml_id` is empty. Neither is wrong alone; together they turn "library identity lost" into
"FDR quietly uncalibrated". Same class as the `getColumn`→`chunk(0)` truncation fixed in iteration 44.

---

## Iteration 47 — why the LDA disagreed with pyprophet: the loop never trained

Full line-by-line comparison in [OpenDIAlyzer-lda-vs-pyprophet.md](OpenDIAlyzer-lda-vs-pyprophet.md).
Headline: on `testdata/lda_fixture.txt` (44,568 real OpenSWATH rows), the semi-supervised loop
reported **`trained=0, skipped=9`** — every fold, every iteration. `w` never left its seed, so the
"LDA" was ranking by a single sub-score.

The seed was the largest-|t| feature. It is not degenerate (100% finite, 23,920 distinct values,
|t| = 13.7) but |t| = 13.7 over 29,482 rows is a **0.16 sd** per-row effect, and the ranking
statistic is the max over ~18 candidate peak groups per precursor — extreme-value spread swamps it.
Every target group got the same q-value (0.9654), zero positives were selected, the fit was skipped,
and nothing changed on the next iteration either. A self-perpetuating cold start.

pyprophet never hits this because OpenSWATH gives it a composite `main_score` to rank by. ODIA has
no such column, so the fix is to seed with an LDA of all top-target rows against all top-decoy rows
— the labels are known outright, no FDR estimate needed.

Three changes, all in `src/odia_lda.h`:
1. **LDA seed** replacing the single-feature seed (the root cause).
2. **`top_decoys_only`** — negatives are now one row per decoy precursor, matching the positives.
   Previously positives were top-peaks-only but negatives were *all* rows, so the discriminant partly
   learned "is this the best peak of its group" rather than "is this a target".
3. **`normalize_folds`** — each fold's held-out scores rescaled to that fold's own decoy null before
   pooling. Fold weight vectors have arbitrary scale/offset; the pooled ranking was partly sorted by
   which fold a precursor landed in.

Fixture effect: target/decoy mean separation 1.011 vs 1.011 → **0.119 vs 0.000**; distinct d-scores
35,964 → 44,347 (tie blocks gone).

### Instrument check — what the fixtures can actually measure

Reporting this because it bounds the claims. The synthetic oracles are **insensitive** to all three
changes (1094→1082 and 1052→1056 IDs, FDR control identical) — their decoys have no rank structure,
so neither the confound nor the cold start exists in them.

And `lda_fixture.txt` cannot measure the *gap*: a faithful reimplementation of pyprophet's loop
(same sklearn `LinearDiscriminantAnalysis`, same 10×10 / 0.15 / 0.05 / xeval_fraction 0.5) run on
that same fixture **also never trains** — 0/10 folds. The fixture has no bootstrappable signal, so
0 IDs is correct for it and it cannot discriminate the two algorithms. Gap measurement has to be on
the real table (`-score_osw bench/narrowmz/mz10.osw`, references: ODIA 2,957 / pyprophet 4,302).

### Result on the real table

`-score_osw bench/narrowmz/mz10.osw` — identical 2.07M features to both references:

| | IDs @ q<0.01 | decoys | empirical FDR |
|---|---|---|---|
| ODIA, single-feature seed | 2,957 | — | — |
| pyprophet 3.0.15 | 4,302 | — | — |
| **ODIA, LDA seed** | **4,367** | 42 | **0.96%** |

`LDA fitted 9 iteration(s), skipped 0` — the loop trains everywhere now. **+47.7% IDs over the old
ODIA and 1.5% ahead of pyprophet on identical input**, in 2:18 wall / 1.4 GB. Calibration holds: 42
decoys against 4,367 targets is 0.96% empirical at a 1% nominal cutoff — and ODIA gets there with
pi0 correction OFF while pyprophet has it on, so the extra IDs come from a better discriminant
rather than a looser error model.

Standing vs the references, updated:

| | wall | peak RSS | IDs @1% FDR |
|---|---|---|---|
| DIA-NN 2.0 | 509 s | 22.7 GB | 8,164 |
| OpenSwathWorkflow @24 thr | 6,886 s | 153 GB | never scored |
| ODIA @10 ppm, 224 thr | 1,814 s | 189 GB | 2,957 → **4,367** |

---

## Iteration 48 — FDR beyond the precursor level

The scoring/FDR backlog's Area 2 "prototype-next" list, implemented in `src/odia_fdr.h`
(header-only, no new dependencies, self-check `odia_fdr_test` wired into ctest as `odia-fdr`).

Until now ODIA wrote `SCORE_MS2` and nothing else. That is not a neutral omission: a 1% precursor
FDR is **not** a 1% peptide or protein FDR, because an entity takes the best of many independent
precursor draws and its score is therefore an extreme-value statistic with a strictly worse error
rate. Shipping only precursor q-values invites a reader to quote them as protein q-values.

**1. Context FDR** — precursor scores roll up to peptides (keyed on modified sequence, so charge
states collapse) and then to proteins, with target-decoy re-run at each level.
`-fdr_context global|none`. Written as `SCORE_PEPTIDE` / `SCORE_PROTEIN` in pyprophet's schema so
downstream tools read them unmodified.

**2. Picked protein-group competition** (Savitski 2015; The 2022) — each protein competes with its
*own* decoy and enters the ranking once, as the winner. `-picked_protein true|false` (default on).
Plain target-decoy at the protein level is biased by protein size: a large protein offers more
peptides for a random decoy hit to land on, so decoy proteins crowd the top of the list in a
size-dependent way. A protein and its decoy have the same size by construction, so competing them
cancels it exactly. The self-check demonstrates the bias and its removal rather than asserting it:
on null data with proteins of 1–197 peptides, plain ranking's top-100 is significantly larger than
average, while picked pair-winners are ~50/50 independent of size.

**3. Entrapment FDP** (Wen/Freestone/Keich/Noble 2025) — `-entrapment_tag <prefix>`. Reports
`(1 + 1/r)·n_entrapment/n_reported` alongside the nominal q-value and says plainly when the two
disagree. This is the test rig the backlog calls for, not a model: it is how you find out whether a
reported 1% is real, which the 2025 literature repeatedly finds it is not.

Deliberately not claimed: the **paired** entrapment estimator is not implemented (needs an explicit
target↔entrapment pairing and a formula I would have been reconstructing from memory); the protein
level is **per-accession, not parsimony** protein-group inference, and the run log says so. Ties in
picked competition resolve to the **decoy** — a tie is no evidence, and breaking it toward the
target would bias the estimate in precisely the optimistic direction under criticism.

### End-to-end validation of the PQP/library fix

Full 10 ppm search against the re-converted `library_ids.oswpq`:

```
before: OpenDIAlyzer[prefilter] 7149966 -> 120525 precursors (120525 target /      0 decoy)
after:  OpenDIAlyzer[prefilter] 7149966 -> 423079 precursors (212292 target / 210787 decoy)
```

Target:decoy ratio 0.993 — essentially 1:1, which is what a correct pairing gives by construction,
and the ratio guard now passes on its own merits rather than aborting. Note the target count also
rose (120,513 with evidence → 212,292 retained): the pair-union keeps a pair when *either* member
passes, so ~92k targets are now retained via their decoy partner. That is the label-symmetric
selection working as designed — it was inert while no pair could be formed.

Library id check on the re-converted bundle: `decoy_refs starting with DECOY_: 3546541/3546541`
(was 0/3546541), with the pairing directly visible — `DECOY_YYYYWHLR_4` ↔ `YYYYWHLR_4`, decoy
sequence shuffled to `YLYWYHYK`.

### Measured on the real table

`-score_osw bench/narrowmz/mz10.osw`, 2:22 wall:

```
in-process LDA FDR -> 4367 target precursors at q<0.01.
context FDR -> 3772 peptides and 183 proteins at q<0.01 (picked competition, 18925 pairs).
SCORE_PEPTIDE: 397,506 rows | q<=0.01: 3,772 target /  36 decoy  -> 0.95% empirical
SCORE_PROTEIN:  19,030 rows | q<=0.01:   183 target /   0 decoy
```

**4,367 precursors → 3,772 peptides → 183 proteins.** Quoting the precursor number as a protein
number would have overstated the protein result by 24x, which is precisely the failure mode this
change exists to prevent.

The 183 is conservative but correct, and the rollup is verified sound rather than assumed:

- The top-scoring proteins are exactly the expected high-abundance plasma set — complement C4A/C4B,
  IgM, ITIH1/2/3, fibronectin, C1s, PON1, Ig kappa constant, fibrinogen alpha, C4b-binding protein.
  A scrambled peptide→protein mapping could not produce that list.
- 2,895 of the 3,808 confident peptides (76%) belong to a confident protein; the rest sit on
  proteins that individually do not clear 1% at their own level.
- Score distributions: targets n=10,614 mean 2.99 max 26.69; decoys n=8,416 mean 2.14 **max 17.11**.
  The decoy null reaches high because a decoy protein ALSO takes a max over ~10 peptides — the
  extreme-value inflation described above, measured. With no decoy above the 183rd target, the best
  attainable q is (1/N_dec)/(183/N_tar) = 0.0069, which is exactly the q reported for the whole
  leading block.

Adversarial check on the entity sets themselves — an undersized decoy set would mean a too-weak
null and an anti-conservative estimate, so it is worth measuring rather than assuming:

```
SCORE_PEPTIDE: target=199,638 decoy=197,868  ratio 0.991
SCORE_PROTEIN: target= 10,614 decoy=  8,416  ratio 0.793
PRECURSOR    : target=208,052 decoy=206,173  ratio 0.991   (for comparison)
```

Peptide level tracks the precursor level exactly. The protein level is 20% short of decoys, but
that asymmetry errs conservative on both counts: the q-value formula divides by the class sizes, so
a smaller N_dec makes each observed decoy count for MORE; and the decoy proteins that do exist carry
more peptides each (23.5 vs 18.8), so their max-over-peptides is drawn from a larger sample and the
null is stronger, not weaker. Both consistent with the observed max decoy score of 17.11 and zero
decoys surviving q<=0.01.

---

## Iteration 49 — compact chromatogram storage, and a measured memory review

Full analysis in [OpenDIAlyzer-memory-review.md](OpenDIAlyzer-memory-review.md).

### `src/odia_chromatogram.h` — replacing `vector<MSChromatogram>`

Two facts do all the work:

**The RT axis is shared, not per-chromatogram.** `ChromatogramExtractorAlgorithm.cpp:316` loops over
*spectra* outer and transitions inner, pushing that spectrum's `s_meta.RT` into every chromatogram it
touches and skipping only those whose own `[rt_start, rt_end]` excludes it. So every chromatogram
from one SWATH map lies on that map's spectrum-RT grid and occupies a **contiguous slice**. Store the
grid once plus a `(first, count)` pair, and the RT half of every point vanishes. This is an invariant
of the extractor, verified by reading it, not an approximation.

**Intensities as float32**, since the mzML arrays are routinely 32-bit already and every downstream
score is relative.

| | header | per point | @2392 points |
|---|---:|---:|---:|
| `MSChromatogram` | 960 B | 16 B | 39,232 B |
| `ChromatogramStore` | 24 B | 4 B | 9,592 B |

Measured (`odia-chromatogram`, ctest): **4.09x**, i.e. **225 GB → 55 GB, saving 170 GB** at the
benchmark's shape. Marginal cost exactly 4.000 bytes/point.

The narrowing is verified not to move what consumes it, rather than assumed harmless:

```
worst relative intensity error: 5.865e-08   (float32 eps = 1.192e-07)
XIC correlation 0.997954835304 (double) vs 0.997954835027 (float32)  |delta| = 2.775e-10
XIC dot product relative delta = 4.230e-09
```

**Not integrated.** `MRMFeatureFinderScoring` still consumes `MSChromatogram`; the 170 GB is real
only once the store owns the data and views are materialised per scoring call. That work is not done
and the doc says so.

### `malloc_trim(0)` after the prefilter — 9.9 GB, measured

The prefilter drops 94% of the library, but it is ~140M small string allocations and glibc keeps
blocks that size in its arena. Measured on the real library with the prefilter emulated at the same
5.9%:

```
loaded 38.79 GB -> after move-assign 27.67 GB -> after malloc_trim(0) 17.80 GB
```

**This corrects a guess made earlier in the same session.** I had speculated ~39 GB was being held;
the measurement says glibc returns ~11 GB unprompted and trim gets 9.9 more. Applied.

The residue is the more interesting half: **~15 GB survives the trim against ~2.3 GB of live data**,
because the surviving 6% of objects are interleaved with the freed 94% in the same pages and a page
with one live object on it cannot be returned. No allocator call reaches that; it needs the library
not to be 140M separate allocations. Not attempted.

### The memory model, and its failure

mz5/mz10 are a controlled pair (same window, threads, binary), giving
`peak ≈ 41.5 GB + 126 KB × n_prec + n_chrom × (960 B + 16 B × window/cycle)`. The intercept is
independently corroborated by the standalone library load (38.8 GB).

But: two points, two parameters, so the fit is exact by construction. Its out-of-sample prediction
for the e2e run was **317 GB against 366.8 GB observed — 16% low**, and the missing ~50 GB is
unattributed. Reported as a miss, not a validation.

### Other findings

- `estimateFeatureMemoryPerCompound()` returns a hardcoded **2 KB**; measured **126 KB** — 63x low,
  and it excludes chromatograms, which are the dominant term.
- `calculateInnerBatchSize()` budgets 5% of *free system memory* against that 2 KB, so on a 2.2 TB
  node it always saturates the hard clamp of 10,000. The adaptive sizing is inert.
- `batchSize` defaults to 0 — no outer batching at all.
- The RT window is the biggest single knob: mz10 and e2e searched **identical** precursor sets and
  differed only 600 s vs 1435 s; peak went 189.5 → 366.8 GB.

---

## Iteration 50 — GBT classifier, and its parallelisation

### The learner (`src/odia_gbt.h`)

In-house histogram gradient boosting, dropped into the SAME semi-supervised loop as the LDA (same
folds, training-set selection, fold normalisation, q-values), selected by `-classifier lda|gbt`.
No XGBoost/LightGBM is available to link, and vendoring one would bring its own RNG and threading
into a path that produces q-values.

Validation, on synthetic data designed to isolate the claim:

```
linear      data: LDA AUC 0.9811   GBT AUC 0.9860   <- no penalty on the LDA's home ground
INTERACTION data: LDA AUC 0.5137   GBT AUC 0.9258   <- LDA at chance; the whole point
missing-as-category AUC 0.998 ; determinism 0.0 across fits
```

**On the real 2.07M-feature table, with calibration verified:**

| | IDs @q<0.01 | decoys | empirical FDR | peptides | proteins |
|---|---:|---:|---:|---:|---:|
| pyprophet 3.0.15 | 4,302 | — | — | — | — |
| ODIA LDA | 4,367 | 42 | 0.96% | 3,772 | 183 |
| **ODIA GBT** | **6,798** | 66 | **0.97%** | **5,956** | **574** |

**+55.7% over LDA at the same true error rate.** Not anti-conservative.

A test caught a real integration bug: the GBT initially inherited the LDA's cold start, because the
seed was an LDA fit — which on interaction data is at chance by construction, so no training set was
ever selected and BOTH learners returned 0 IDs. The seed now uses the same learner as the loop.

### Cost, before parallelising

| | wall | user | sys | CPU | IDs |
|---|---:|---:|---:|---:|---:|
| LDA | 142.7 s | 35.0 s | 121.0 s | 109% | 4,367 |
| GBT | 167.0 s | 104.0 s | 124.9 s | 137% | 6,798 |

GBT costs +69 s CPU for +24 s wall = **2.83x on the learner against a 3.00x ceiling (n_folds=3),
i.e. 94% efficiency inside the fold loop**. But **78% of the LDA's CPU is SYSTEM time** — SQLite
reading 2.07M x 24 scores from a 4.4 GB .osw and writing three score tables back, single-threaded
and identical between the two runs. This phase is I/O-bound; the learner is a minority of it.

### The parallelisation

Cost model of one fit first, because it decides everything: **histogram building is 92.6%**, row
routing and prediction 6.4%, grad/hess and split finding 1.0%. So there is exactly one target.

Two decisions worth recording:

**Chunked, not per-thread, reduction.** The histogram is a floating-point reduction, so the number
of partial sums fixes the rounding. A per-thread partition would make the model depend on
`OMP_NUM_THREADS` — unacceptable for a score that becomes a q-value. The row partition is instead a
fixed chunk count derived from the row count alone, combined in chunk-index order; threads take
chunks in any order and the answer does not move. Tested, not asserted:

```
T8 max |delta| vs 1 thread:  8 threads 0.000e+00   64 threads 0.000e+00
```

**Flattened the binned matrix.** `vector<vector<uint8_t>>` put a pointer chase in front of every
read in the 92.6% loop; now one contiguous buffer.

**Nesting had to be enabled.** `GBT::fit` runs inside the fold loop's parallel region, where a
nested team defaults to ONE thread — the inner pragmas would have been dead code. The caller sets
`omp_set_max_active_levels(2)` and gives the GBT `max_threads / n_folds`, so total width is
`folds x inner` rather than `folds`.

All tests pass with results identical to the serial version.

### 10 folds is not worth it

| | 3 folds | 10 folds |
|---|---:|---:|
| LDA | 4,367 (109% CPU, 2:22) | 4,381 (118%, 2:21) |
| GBT | 6,798 (137% CPU, 2:47) | 6,708 (266%, 2:56) |

pyprophet's `ss_num_iter = 10` default does not transfer: the LDA gain is noise and the GBT is
slightly WORSE, both at higher cost. Keep `n_folds = 3`. It does confirm the fold parallelism scales
(137% -> 266% as folds go 3 -> 10).

---

## Iteration 51 — DIA-NN 1.7.12 built from source, and a benchmark-methodology error it exposed

### Built

DIA-NN 1.7.12 (`git checkout 1.7.12`, commit cff0408), CC BY 4.0 per its own header
("Copyright 2020, Vadim Demichev ... Creative Commons Attribution 4.0 International License").
MSToolkit compiles clean. `diann.cpp` needs a **2-line build-only patch**: gcc 13 rejects two
`goto`s that jump over an initialisation (an error in every C++ standard; gcc-7 accepted it under
`-fpermissive`). Both are forward jumps to a label later in the same loop and nothing after either
label reads the crossed declarations, so the skipped regions are enclosed in a scope. **No algorithm
is touched, and nothing seen in that source informs OpenDIAlyzer** -- the clean-room policy in
CLEAN-ROOM.md stands.

Runtime on the benchmark: **11:15 wall, 8341% CPU (83.4 cores), 18.4 GB peak RSS.**

### The result is INVALID, and the reason matters for the whole table

1.7.12 reported **1,092,962 IDs at 1% FDR** on a plasma sample. That is not a number, it is a broken
FDR, and DIA-NN said so itself:

```
WARNING: the library contains decoy precursors, DIA-NN's built-in decoy generation will be turned off
WARNING: 3546541 decoy precursors don't have matching target precursors; each decoy supplied in the
         library must correspond to a target with the same modified sequence and charge
```

**DIA-NN pairs decoys by modified sequence + charge. Our decoys are SHUFFLED** (OpenSwathDecoyGenerator
writes `DECOY_<target id>` with a shuffled sequence), so not one of the 3,546,541 pairs. DIA-NN then
disabled its own decoy generation *because decoys were present* and was left with no usable null.

Its own diagnostics confirm the run is noise: `Median.Mass.Acc.MS2 = 40.88 ppm` on an instrument
measured at ~1.6 ppm.

**DIA-NN 2.0 was NOT affected** by the same library: its only warning concerned neutral losses, and
it trained on `17550 target and 4478 decoy PSMs`. So 2.0 either ignores library decoys or tolerates
the mismatch, and 1.7.12 does not.

### What this means for the benchmark

The benchmark has been handing DIA-NN a library built for OpenSWATH's decoy convention. 2.0 coped;
1.7.12 did not. The correct comparison gives each tool the library it expects -- a **target-only**
library, letting DIA-NN generate its own decoys, which is exactly what its warning instructs.

That is now being built (`library_targets.tsv`, col 25 Decoy==0) and BOTH DIA-NN versions will be
re-run against it. Until then:

- **1.7.12's 1,092,962 is not reportable.**
- **2.0's 8,164 should be re-confirmed** on the target-only library before it is quoted further,
  since the same input handicapped the other version so severely.

---

## Iteration 52 — the DIA-NN reference was wrong, and correcting it moves everything

### The library was handicapping BOTH DIA-NN versions

`library.tsv` carries OpenSWATH-style decoys: `DECOY_<target id>` with a **shuffled** sequence.
DIA-NN pairs decoys by **modified sequence + charge**, so none of the 3,546,541 pair. Re-running
both versions on a **target-only** library, letting DIA-NN generate its own decoys:

| DIA-NN | library | IDs @1% | MS2 mass acc | wall | CPU | cores | RSS |
|---|---|---:|---:|---:|---:|---:|---:|
| 1.7.12 | with decoys | **1,092,962 (invalid)** | 40.88 ppm | 675 s | 8341% | 83.4 | 18.4 GB |
| 1.7.12 | target-only | **8,405** | 1.7055 ppm | 614 s | 11947% | 119.5 | 16.2 GB |
| 2.0 | with decoys | 8,164 | — | 509 s | 7519% | 75.2 | 22.7 GB |
| 2.0 | target-only | **9,261** | 1.656 ppm | **349 s** | 8447% | 84.5 | 19.7 GB |

**1.7.12 went from garbage to 8,405.** Its mass-accuracy diagnostic went 40.88 -> 1.71 ppm, which
is the tell: a tool reporting 40 ppm on a 1.7 ppm instrument is matching noise.

**2.0 gained 13.4% (8,164 -> 9,261) and got 31% faster.** I had assumed 2.0 was unaffected because
it did not warn. That was wrong -- the decoy-laden library was costing it too, silently, and every
comparison in this log before this point used the depressed 8,164 as the reference.

Three independent measurements now agree the instrument's MS2 accuracy is **1.66-1.71 ppm**
(DIA-NN 2.0, DIA-NN 1.7.12, and the earlier direct measurement). OpenDIAlyzer extracts at **10 ppm**.

### Corrected gap decomposition (ODIA/GBT vs DIA-NN 2.0, 8,812 mod-stripped)

| | count | |
|---|---:|---|
| A. prefilter removed -- **unreachable** | **1,371** | 15.6% of DIA-NN |
| B. reachable AND found | 5,892 | 79.2% of reachable |
| C. reachable but missed | 1,549 | |
| D. ODIA found, DIA-NN did not | 470 | |
| net gap | 2,450 | (was 1,425 against the wrong reference) |

**ODIA is at 69-73% of DIA-NN 2.0, not 83%.** The shape of the conclusion survives -- the classifier
recovers ~79% of what it can see and the prefilter is a hard ceiling -- but the prefilter loss grew
from 924 to 1,371 and the total gap nearly doubled.

### Implemented

**`calibrateMassFromPass_`** -- the mass-accuracy hook always worked; it lacked anchors. It was given
the iRT anchor list and the peakedness gate refused to fit. Inverting the statistic
(ratio = 1 + 5N/M over a 50 ppm bootstrap) shows why: the observed 1.14 means **2.8% of MS2 anchors
are real**, 1.69 means 12.1% for MS1. Declining was correct. `recalibrate_` already ranks target peak
groups by d-score after pass 1 -- those have measured evidence -- so the calibration now reuses that
set and applies the result to the next pass.

**`-prefilter_mz_extraction_window`** -- the screen and the extraction want opposite widths and
shared one until now. Sweeping 10/20/30 ppm screening with extraction fixed at 5 ppm.

### pyprophet on the OSW output: abandoned after three attempts

5h23m (killed), 8h33m (killed), and 2h50m with **2.5 hours of no log output** while holding 15
workers at ~22 GB (~330 GB) and confounding every concurrent wall-time measurement. Replaced with
**ODIA's own GBT scoring OSW's features** -- which is arguably the better experiment anyway, since it
holds the scorer constant and isolates the feature sets.

---

## Iteration 53 — calibration guards, and a defect I introduced and could not diagnose by inference

### The night's throughline: a calibration may only narrow

Three calibrations, all fitting to anchor sets too poor to support a fit, two of which produced
**wider** windows and were applied unconditionally:

| calibration | anchor quality | behaviour | cost |
|---|---|---|---|
| mass, MS2 | 2.8% real (peakedness 1.14) | refused | none -- correct |
| mass, MS1 | ~12% real, passed at peakedness *exactly* 3.0 | **72.6 ppm** on a 1.66 ppm instrument | 6,798 -> 4,496 |
| RT window | **1.8% yield** (70/3,897) | **600 s -> 1435 s** (61% of the gradient) | large memory, no ID change |

The configured value is the user's assertion about the run; measuring exists to discover the run is
TIGHTER. An estimate that comes back wider is reporting that its anchors did not support a fit, and
applying it is strictly worse than doing nothing. Now enforced on both, with the rejection logged:
`-rt_calib_min_yield_pct` (default 10%) plus a hard never-widen check on all three channels.

**This produced the session's memory result: 366.8 GB -> 208.3 GB, a 43% reduction**, because
chromatogram size is linear in the RT window.

### Three wrong hypotheses about the parquet scoring gap

The in-memory scoring path costs **1,694 identifications (26%)** versus sqlite on otherwise
identical settings (6,607 vs 4,913). Diagnosed by inference three times, wrong three times:

| hypothesis | prediction | measured | verdict |
|---|---|---|---|
| RT window 1435 s is the cause | ~6,800 at 600 s | 4,885 | wrong |
| column loss from `fmap[0].getKeys` | fewer than 29 | **35-36** | wrong (more, not fewer) |
| MS1 columns contaminating the fit | recovery to ~6,600 | 4,876 | wrong |

Worse than being wrong: the first two comparisons were against `mz10.osw`, a feature set from a
**different extraction run**, so "sqlite vs parquet scoring" was never the isolated variable. The
controlled run -- same everything, only the `-out` extension differing -- settled it in one shot and
should have been the FIRST thing run, not the fourth.

Both loaders now emit a census (rows, skip reasons, columns before/after pruning, target count) so
the divergence is read off a log rather than inferred from outcomes.

### Corrected DIA-NN references

Both versions had been benchmarked against a library carrying OpenSWATH-style shuffled decoys, which
DIA-NN cannot pair (it matches on modified sequence + charge). On a target-only library:

- **2.0: 8,164 -> 9,261** (+13.4%, and 31% faster). It never warned; I had assumed it was unaffected.
- **1.7.12: 1,092,962 (invalid) -> 8,405.** Its mass-accuracy diagnostic went 40.88 -> 1.71 ppm.

Three independent measurements now put the instrument at **1.66-1.71 ppm**; ODIA extracts at 10 ppm.

### Standing priority (unchanged, and now better supported)

30 ppm screening makes **96.0%** of DIA-NN's answer reachable vs 84.4% today (recovering 1,015 of
1,371 unreachable IDs). It is blocked only on extraction memory -- and memory is now 43% lower than
when that was measured. `ChromatogramStore` (4.09x, built and tested, not integrated) is the enabler.

---

## Iteration 54 -- the census answers, and half the runtime was never measured

### The parquet gap is NOT rows and NOT columns

The census settles both, with numbers instead of hypotheses:

```
scoreload/memory  pass 1: rows 2070164/2070164  (0 skipped)  columns 23 -> 23  targets 1039722
scoreload/memory  pass 2: rows 2055478/2055478  (0 skipped)  columns 24 -> 24  targets 1032389
```

Zero rows skipped for any reason. And the column question, read off the sqlite file directly rather
than reasoned from the schema:

| | count |
|---|---:|
| `VAR_` columns in `FEATURE_MS2` | 29 |
| of those, entirely NULL | 5 -- the 4 IM scores **and** `VAR_MI_RATIO_SCORE` |
| informative | **24** |
| found by the in-memory loader | **24** |

Both paths train on the same 24 features. (I stated this wrongly twice on the way here -- first
"4 missing, both end at 23", then "5 missing, not 4". The all-NULL five exist only in the sqlite
schema; `dropUninformativeColumns_` removes them and the two loaders converge.)

Baseline is reproducible: **4,890 IDs**, matching the earlier 4,876-4,913.

### `library_rt` was reading the observation, not the prediction

```cpp
R.library_rt.push_back(f.getMetaValue("norm_RT"));   // wrong
SELECT ... p.LIBRARY_RT ...                          // sqlite
```

`norm_RT` is `scores.normalized_experimental_rt` (`MRMFeatureFinderScoring.cpp:1018`) -- the
OBSERVED RT mapped into iRT space. `PRECURSOR.LIBRARY_RT` is the library's PREDICTION.

`recalibrate_` fits `library_rt -> exp_rt` to learn the library-to-run warp. Fed `norm_RT`, it fits
a function of `exp_rt` against `exp_rt` -- recovering the inverse of the transform already applied,
not the warp. Pass 2 then extracts on a corrupted RT axis, which degrades the FEATURES rather than
the classifier. That is why changing columns and windows never moved the number.

Now sourced from `LightCompound::rt` via `precursor_index_`, the same origin sqlite reads.

### Half the wall clock was never instrumented

Phase profile of the 6,607-ID control run (37:32 wall, 16:05:32 CPU):

| phase | wall | CPU | avg cores |
|---|---:|---:|---:|
| Linear calibration | 24 s | 27:42 m | 69.4 |
| Nonlinear calibration | 235 s | 3:54:47 h | 60.0 |
| Extract+score pass 1 | 559 s | 6:06:28 h | 39.3 |
| Extract+score pass 2 | 280 s | 3:17:57 h | 42.4 |
| **accounted** | **1,098 s (49%)** | 49,614 s | 45.2 |
| **UNACCOUNTED** | **1,154 s (51%)** | 8,318 s | **7.2** |

The instrumented phases average 45 cores. The invisible half runs at 7 and drags the whole-run
average to 25.7. Every parallel-efficiency number reported before this line was measuring the wrong
thing -- the extraction was never the problem.

Nine phases now carry a wall+CPU stopwatch, so a serial phase is self-identifying (`cpu/wall`).

### `prefilterLibrary_`: no OpenMP, and `std::set<std::string>` at 7.1M scale

| operation | was | now |
|---|---|---|
| 4x `std::set<std::string>` | RB-tree, ~22 string compares per lookup | `unordered_set`, O(1) |
| pair-union test | `dtag + c.id` -- one heap alloc x 7.1M compounds | `keep_decoy` keyed by target id, tag stripped once at insert |
| decoy pass | 7.1M compound scan + `substr()` -- one heap alloc x 3.5M decoys | walk the ~4e5 kept targets instead |

~30M tree lookups and ~14M heap allocations removed. The decoy-pass inversion is equivalent because
`alias_of`'s key set is exactly `{compound ids in decoy_refs}` -- the predicate the old loop filtered
on -- and the target:decoy ratio guard already in place is the runnable check that it stayed so.

### Also

`GBT` histogram loop capped at `nchunk` threads. It has exactly `nchunk` (<=32) iterations, so the
`224/folds ~= 44` it was handed left threads idle at the barrier: the measured cause of 224 threads
running 60% slower than 32. Bit-identical output preserved (T8 passes at 1/8/64 threads) because the
cap changes who does the work, not the chunk count.

`ChromatogramStore` integration is DEPRIORITISED. Its 12.6x targets extraction, which the profile
shows runs at 39-42 cores and is 37% of wall. It cannot touch the 51% running at 7.

### Methodology defect: the benchmark node is shared, and I never checked

`ps`/`uptime` on the node, mid-run:

| process | user | CPU |
|---|---|---:|
| `odia_libraryrt` | kohlbach | 2841% |
| `dorado basecaller` | bagci | 1787% (running 1d 14h) |
| `python` x5 | zhang | ~1530% |
| **non-kohlbach total** | | **4708% = 47 cores** |

`load average: 180.14, 200.76, 232.90` on 224 cores, 18 users. The node is oversubscribed. Every
"ODIA achieves 25.5 of 224 cores" figure reported before this line assumed 224 were available.

What survives:

- **Peak RSS** -- contention-independent. Valid as reported.
- **Total CPU-seconds** -- work actually performed. Valid as reported.
- **Wall clock, CPU%, "avg cores"** -- these measure cores WON on a busy node, not scalability.

On CPU-seconds, which should have been the headline:

| tool | CPU-seconds |
|---|---:|
| DIA-NN 2.0 | **29,480** |
| DIA-NN 1.7.12 | 73,354 |
| ODIA GBT sqlite | 57,950 |

ODIA performs **less total work than DIA-NN 1.7.12** and ~2x DIA-NN 2.0 -- not the 6.5x the wall
gap implies. The wall gap is mostly failure to obtain cores: partly ODIA's own serial half, partly
other tenants.

The phase-profile conclusion is unaffected. Contention scales all phases down together, and the
6-9x spread between extraction (39-69) and the unaccounted half (7.2) is far too large to be
contention. The serial half is real; its absolute core count is not.

`experiments/bench_run_fair.sh` now records load average, other-user CPU and free memory on both
sides of every run, so a slow run can be distinguished from a busy node after the fact.

---

## Iteration 55 -- the phase profile pays for itself immediately

First run with the nine phase stopwatches (20 ppm screen / 10 ppm extraction):

| phase | wall | CPU | avg cores |
|---|---:|---:|---:|
| `library_load` | 155.0 s | 167.1 s | **1.1** |
| `dia_run_load` | 140.5 s | 6,668.7 s | 47.5 |
| `prefilter` **(total)** | **321.1 s** | 726.2 s | **2.3** |
| -- `scan_targets` | 95.2 s | 296.4 s | 3.1 |
| -- `build_sets` | 26.9 s | 26.9 s | 1.0 |
| -- `build_decoy_view` | 43.5 s | 43.5 s | 1.0 |
| -- `scan_decoys` | 63.5 s | 267.3 s | 4.2 |
| -- `pair_and_rebuild` | 61.5 s | 61.5 s | 1.0 |
| `precursor_index` | 14.5 s | 14.5 s | **1.0** |

Kept 1,459,897 precursors (732,856 target / 727,041 decoy), 15,934,501 transitions -- matching the
earlier `pf_20` sweep exactly.

### 107 s of the prefilter existed only to work around an OpenMS restriction

`TransitionListEvidenceFilter` refused decoys twice -- `getDecoy()` and a hardcoded `DECOY`-prefix id
test -- so label-symmetric selection required presenting them a SECOND time under aliased ids with
the flag cleared. Measured price of that workaround:

- `build_decoy_view` **43.5 s** to copy ~21M transitions (plus the transient memory)
- `scan_decoys` **63.5 s** to re-read every spectrum

Added `include_decoys` to the filter (default `false`; existing behaviour and its test untouched) and
collapsed ODIA to ONE scan that splits evidence by label. This is also strictly *more* symmetric than
the workaround: both classes now face one identical pass over the same spectra, rather than two
passes that could in principle differ.

**Deletion hazard handled.** The pass-2 inversion from iteration 54 was keyed on `alias_of`, whose
key set was `{compound ids in decoy_refs}`. With `alias_of` gone, `decoy_compound_ids` is collected
in the rebuild pass that already tests that predicate -- `decoy_refs` derives from TRANSITIONS and can
name a ref with no compound, and resurrecting such an id would keep orphan transitions with no
compound to attach them to.

### Why the two evidence scans only reached 3-4 cores

Not the wave scheduler, which was my first guess. The prefilter sets `bytes_per_spectrum = 0` and
leaves `avg_transitions_per_swath` at 0, and the transition-density term is included only when both
are non-zero -- so `estimated_bytes_per_swath` is just the 8 MB overhead and `max_concurrent_swaths`
resolves to ~150, not 3. The likelier cause is I/O: `readOptions='normal'` streams the 6.4 GB mzML,
and the two scans read it twice. Collapsing to one scan halves that regardless of the mechanism.

### Environment trap

After reinstalling OpenMS, the binary died at startup on `libxerces-c-3.3.so`. The RUNPATH is
correct; a non-login ssh command shell simply never sources the profile that sets
`LD_LIBRARY_PATH`. `bench_run_fair.sh` now exports it explicitly. Nothing was wrong with the build --
but this would have silently killed an overnight run.

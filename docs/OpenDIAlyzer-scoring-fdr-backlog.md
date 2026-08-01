# OpenDIAlyzer scoring / FDR / recalibration — literature backlog (2022–2026)

What to prototype **after** the baseline in-process semi-supervised LDA
(mProphet/PyProphet style), for (1) peak-group scoring, (2) FDR/error estimation,
(3) on-the-fly recalibration + predictor fine-tuning. From a background literature
scan. **Read gains skeptically:** they are largest for immunopeptidomics / PTM /
non-tryptic, and *much smaller for vanilla tryptic bulk DIA* — OpenDIAlyzer's arena.

Current decision: **LDA first** (deterministic, dependency-light, in-process),
boosted trees to the backlog per the ranking below.

## Area 1 — peak-group scoring / rescoring
Same ~20 OpenSWATH `VAR_*` sub-scores as input.
- mProphet LDA (Reiter 2011) — the baseline; linear boundary leaves IDs on the table.
- PyProphet (Rosenberger 2017) — ships LDA / XGBoost / SVM; boosted beats LDA on ID rate.
- Percolator SVM + nested cross-validation (Käll 2007; The 2016) — the CV "brew" is what lets you train on decoys without inflating FDR (couples to Area 2).
- mokapot (Fondrie & Noble 2021) — Percolator with XGBoost: +15% mod PSMs / +19% peptides / +11% proteins vs linear.
- DL rescoring (adds *new* features → classifier): MS2Rescore (Declercq 2022), MS²Rescore 3.0 (Buur 2024), TIMSRescore (Declercq 2025, timsTOF/IM), Oktoberfest/Prosit (Picciani 2023), MSBooster (Yang 2023).
- DIA-NN's own scorer (Demichev 2020): a per-run **ensemble of small feed-forward NNs** on the peak score vector — the non-linear bar to clear.
- XIC neural scorers: Alpha-XIC (Song 2022, +9–16% precursors), SeFilter-DIA (He 2024), WinnowNet (Feng 2025). Best under interference; GPU + least deterministic.

**Prototype-next:** (1) **cross-validated GBT (XGBoost/LightGBM or in-house histogram GBT) on the existing sub-scores, fixed seed** — ~10–20% ID lift, drop-in, deterministic when pinned. (2) DL-predicted MS2/RT/IM agreement features into that GBT (MS2Rescore/MSBooster pattern) — via ONNX Runtime C++ or Koina; smaller lift for tryptic. (3) XIC neural scorer — research backlog.

### Boosting / mokapot arm — BACKLOG ITEM (re-implement in OpenMS 3.6)
The concrete next classifier after the MVP LDA. Two coupled ideas:
- **Gradient-boosted trees as the semi-supervised learner.** Swap LDA's linear
  boundary for XGBoost / LightGBM / HistGradientBoosting on the *same* ~20 `VAR_*`
  sub-scores. Evidence: pyprophet already ships XGBoost as a non-default classifier;
  **mokapot** (Fondrie & Noble 2021) = a Percolator reimplementation whose headline
  is exactly this swap → **+15% modified PSMs / +19% peptides / +11% proteins** at 1%
  FDR over the linear model. Expected ~10–20% more IDs; captures the non-linear
  sub-score interactions LDA cannot.
- **The mokapot/Percolator "brew":** semi-supervised iteration with **nested k-fold
  cross-validation by group** so a peak group is never scored by a model trained on
  its own precursor — this is what keeps a decoy-trained booster from producing the
  optimistic FDR the 2025 entrapment critiques document. The CV harness, not the
  learner, is the load-bearing part; ship it **with** the entrapment harness (Area 2).

**Implementation target: OpenMS 3.6** (the version OpenDIAlyzer links). What 3.6
already provides to build on (`src/openms/include/OpenMS/ML/`):
- `ML/SVM/SimpleSVM` — a ready SVM classifier → the Percolator-*SVM* arm with little new code.
- `ML/CROSSVALIDATION` — reuse for the semi-supervised CV brew (don't hand-roll folds).
- `ML/ROCCURVE`, `ML/REGRESSION`, `ML/GRIDSEARCH` — evaluation, PEP/q monotonization, hyper-search.
- `ML/ONNX` — the path for any DL-rescoring / NN-scorer feature model.
- **Eigen** (contrib) for the linear algebra.
- **No bundled GBT in 3.6** — so the boosting arm must either vendor LightGBM/XGBoost
  into the OpenDIAlyzer build or add an in-house histogram GBT (Eigen-assisted), written
  as an OpenMS-3.6-native component (NOT pushed to the OpenMS repo — it lives in
  OpenDIAlyzer, but uses 3.6 facilities and conventions).
- **Plug point:** the classifier-agnostic interface `src/odia_lda.h`
  (`features, labels, groups -> dscore, qvalue`) is the seam. The MVP LDA is
  deliberately pure-C++/standalone (portable + self-testable); the SVM and GBT arms
  reuse OpenMS 3.6's `ML/` modules behind the same interface, selected by a
  `-classifier lda|svm|gbt` switch. A/B each against LDA + the entrapment harness.

## Area 2 — FDR / error estimation beyond plain target-decoy
2024–25 message is **skeptical**: DIA target-decoy FDR is frequently *optimistic*.
- Entrapment validation is the gold standard (Wen/Freestone/Keich/Noble, Nat Methods 2025): of three field-standard methods, one is **invalid**, one a **lower bound**, one **valid but under-powered**; they give a more powerful paired estimator and find **no DIA tool consistently controls peptide-level FDR — worst on single-cell**. Madej & Lam (2024): entrapment-*query* protocols violate their own assumptions — build the harness correctly (paired/combined).
- DIA-NN library-free FDR is measurably optimistic (J Proteome Res 2025, 10.1021/acs.jproteome.5c00036).
- Competition-based decoys / double competition (Elias & Gygi 2007; Keich/Noble line); PEP (Käll 2008); context FDR run/experiment/global (Rosenberger 2017); picked protein / protein-group FDR (Savitski 2015; The 2022).
- Decoy quality biases the null: `OpenSwathDecoyGenerator` shuffle/reverse/shift (Röst 2014); shuffle+mutate with identity control is a first-class knob — *this is the null our pyprophet-style FDR trusts* (ties to our P2 decoys).

**Prototype-next:** (1) **a correct entrapment harness first** (paired/combined estimator) — not a model, a test rig; the only way to know reported FDR is real, and it usually isn't. Pure C++/bookkeeping. (2) context/group FDR + picked protein-group rollup. (3) PEP + deliberate decoy-quality control. Guardrail: cross-validated classifier scoring (Area 1) is what prevents the optimistic-FDR trap.

### Status — IMPLEMENTED (`src/odia_fdr.h`, self-check `src/odia_fdr_test.cpp`, ctest `odia-fdr`)

Header-only, no new dependencies, reusing the precursor level's own step-down target-decoy
estimator so there is one implementation of the q/PEP mathematics across all three levels.

| Item | State | Where |
|---|---|---|
| Context FDR — peptide + protein q-values | **done** | `rollUp` + `assignQValues`; `-fdr_context global\|none` |
| Picked protein-group competition | **done** | `pickedCompetition`; `-picked_protein true\|false` |
| Entrapment FDP — **combined** estimator | **done** | `entrapmentFdp`; `-entrapment_tag <prefix>` |
| Entrapment — **paired** estimator | *not done* | needs explicit target↔entrapment pairing in the library; see below |
| PEP | done earlier | `lda_detail::assignQValues` |
| Decoy-quality control | open | `OpenSwathDecoyGenerator` knob, not an ODIA change |

Written into the `.osw` as `SCORE_PEPTIDE` / `SCORE_PROTEIN` in the schema pyprophet produces, so
downstream tools read them without knowing which tool wrote the file. Previously ODIA wrote
`SCORE_MS2` and nothing else, which invited readers to treat a precursor q-value as a protein one.

Three notes on what was deliberately *not* claimed:

- **The paired entrapment estimator is not implemented.** It is more powerful than the combined one
  but needs a different formula and an explicit pairing; writing it from a half-remembered
  description would be worse than not having it. The combined estimator is well defined from the
  counts alone and is the conservative choice.
- **Protein level is per-accession, not parsimony.** A shared peptide contributes to every protein
  it maps to. This is the "protein" context of Rosenberger 2017, not protein-group inference; the
  run log says so explicitly rather than letting the table name imply more than was computed.
- **Ties in picked competition go to the decoy.** A tie carries no evidence either way, and
  resolving it toward the target would bias the estimate in exactly the optimistic direction the
  2025 entrapment critiques document.

## Area 3 — on-the-fly recalibration & predictor fine-tuning
- Per-run recalibration (everyone does it): DIA-NN / Spectronaut (Bruderer 2015) recenter mass (dynamic ppm), fit non-linear RT (iRT→observed), recalibrate IM/CCS; pick tolerances by maximizing IDs@1%. Cheap, deterministic — do **before** fine-tuning.
- Predictor fine-tuning at inference (the frontier): DIA-NN fine-tunes RT/IM (helps unseen mods); **alphaDIA transfer learning** (Wallmann, Nat Biotechnol 2025) re-fits PeptDeep per experiment (~120k precursors @60SPD Astral), gains largest for non-standard chemistries, modest for tryptic Orbitrap; AlphaPeptDeep (Zeng 2022); Chronologer (Wilburn 2023); IM2Deep (Declercq 2025). DeepLC transfer-learning (Bouwmeester, Nat Commun 2026) cleanly isolates: a few hundred run-specific peptides beat a fixed global calibration, *especially across modification/setup shifts*.
- **Determinism caveat:** predictor fine-tuning is the least reproducible step (data-order, GPU, early-stopping) — freeze the training set, pin seeds, single-thread, version the weights.

**Prototype-next:** (1) **robust iterative per-run RT + m/z (+IM) recalibration** — what OpenDIAlyzer builds now (monotone/LOESS RT from confident anchors, dynamic mass recentering). (2) on-the-fly fine-tune of the RT (then IM) predictor on confident IDs — ONNX Runtime training / LibTorch / small in-house net; best for mods/labels/ToF. (3) MS2-intensity predictor fine-tuning — backlog.

## Consolidated "try-first-after-LDA"
| Area | First prototype | Why | Main cost |
|---|---|---|---|
| Scoring | CV GBT on existing sub-scores | ~10–20% IDs, reuses features | seed/threading determinism |
| FDR | Correct entrapment harness + context/group FDR | catches optimistic DIA FDR | none (algorithmic) |
| Recalibration | Iterative RT+m/z(+IM); fine-tune predictor 2nd | robustness first | fine-tuning nondeterminism |

Through-line: the GBT and the entrapment harness are **coupled** — cross-validated GBT scoring is exactly what keeps a decoy-trained classifier from producing the optimistic FDR the 2025 critiques document. Do those two together; DL rescoring + predictor fine-tuning next; XIC-neural + MS2 fine-tuning are GPU-dependent research.

## Full citations
Reiter 2011 mProphet 10.1038/nmeth.1584 · Rosenberger 2017 PyProphet 10.1038/nmeth.4398 · Käll 2007 Percolator 10.1038/nmeth1113 · The 2016 Percolator 3.0 10.1007/s13361-016-1460-7 · Fondrie & Noble 2021 mokapot 10.1021/acs.jproteome.0c01010 · Gessulat 2019 Prosit 10.1038/s41592-019-0426-7 · Picciani 2023 Oktoberfest 10.1002/pmic.202300112 · Declercq 2022 MS2Rescore 10.1016/j.mcpro.2022.100266 · Buur 2024 MS²Rescore 3.0 10.1021/acs.jproteome.3c00785 · Declercq 2025 TIMSRescore 10.1021/acs.jproteome.4c00609 · Yang 2023 MSBooster 10.1038/s41467-023-40129-9 · Bouwmeester 2021 DeepLC 10.1038/s41592-021-01301-5 · Demichev 2020 DIA-NN 10.1038/s41592-019-0638-x · Song 2022 Alpha-XIC 10.1093/bioinformatics/btab544 · He 2024 SeFilter-DIA 10.1007/s12539-024-00611-4 · Feng 2025 WinnowNet 10.1038/s41467-025-63977-z · Wen 2025 entrapment 10.1038/s41592-025-02719-x · Madej & Lam 2024 10.1002/pmic.202300398 · DIA-NN library-free FDR 2025 10.1021/acs.jproteome.5c00036 · Elias & Gygi 2007 10.1038/nmeth1019 · Käll 2008 PEP 10.1021/pr700739d · Savitski 2015 picked protein 10.1074/mcp.M114.046995 · The 2022 picked protein-group 10.1016/j.mcpro.2022.100437 · Röst 2014 OpenSWATH 10.1038/nbt.2841 · Wallmann 2025 alphaDIA 10.1038/s41587-025-02791-w · Zeng 2022 AlphaPeptDeep 10.1038/s41467-022-34904-3 · Bouwmeester 2026 DeepLC-TL 10.1038/s41467-026-68981-5 · Wilburn 2023 Chronologer 10.1101/2023.05.30.542978 · Declercq 2025 IM2Deep 10.1021/acs.analchem.5c01142 · Bruderer 2015 Spectronaut 10.1074/mcp.M114.044305

Related: docs/OpenDIAlyzer-pyprophet-integration.md (in-process LDA) · docs/OpenDIAlyzer-plan-synthesized.md.

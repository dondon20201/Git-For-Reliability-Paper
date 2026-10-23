# E4b exploratory mechanism re-test - specification (EXPLORATORY addendum)

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md.

## 0. Status and purpose (READ FIRST)
- EXPLORATORY. The confirmatory record for E4 comes from the previously locked H2 result that caused E4 to become FAIL. E4b does not replay, replace or rescore the previously locked H2 result. Instead, every decision in E4b will be labelled EXPLORATORY. E4b seeks to answer a question that arose in the review of E4: is the mechanism genuinely absent, or is the failure caused by the presence of two removable artifacts (F2/F4 collinearity and outcome noise)? This is then used in the decision for Study-3 (reframe versus fitted-model E5) which is made by the researcher during the review. No designs are created and no test or freeze subjects are read.
- BUILD-ONLY, modeling only. E4b reads COMMITTED outputs (E3b denoised reliabilities + E2 features + E3
  set lists); it reads NO raw data, NO participant choices, recomputes NO reliability. Test sealed.
- E4b does NO optimization and constructs NO candidate designs (E5).

## 1. Inputs (committed; build-only; additions-only)
- `elevation/results/e3b/reliability_subject_split_real.csv` (real Set-1 per-set r_sb_A, r_sb_B).
- `elevation/results/e3b/reliability_highprec_sim.csv` (sim high-precision per-set r_sb_hp).
- `elevation/results/e3b/outcome_reliability_real.csv` (cross-half outcome reliability = attenuation ceiling).
- `elevation/results/e2/problem_information_dictionary.csv` (per-problem features).
- `elevation/results/e4/set_features.csv` (the committed per-set features; reuse to avoid re-aggregation,
  but assert a sample re-derives from E2 as in the E4 review).

## 2. The two artifacts E4b removes (pre-declared)
- ARTIFACT 1 outcome noise: replace E4's noisy per-set r_sb with E3b's denoised outcomes - the sim
  high-precision r_sb_hp (near-noise-free) and the real subject-split (fit on half A, evaluate on the
  INDEPENDENT half B).
- ARTIFACT 2 collinearity: E4 found corr(F2 loss_exposure, F4 expected_information) = 0.90 (real) / 0.96
  (sim), which makes the individual partial signs of F2 and F4 unstable - the direct cause of the G-real
  sign-consistency failure. Decorrelate:
  - PRIMARY treatment: orthogonalize F4 on F2 within the analysis sample -> `F4_resid` = residual of F4
    regressed on F2 = "expected information BEYOND loss exposure". Model uses {F1, F2, F3, F4_resid, F5}.
    This directly answers the scientific question: does information add anything once loss exposure is
    accounted for? (Frozen H2's F4 tested this implicitly but unstably.)
  - ROBUSTNESS treatment: drop F4 entirely -> {F1, F2, F3, F5}. Report alongside.

## 3. Models and evaluation (pre-declared; linear; NOT black-box)
Per target (loss_avoid_adj primary), the pre-declared low-dimensional linear model with the decorrelated
feature set, plus BASE_category (F2 loss exposure) and BASE_count (K) baselines, as in E4.
- SIM high-precision track: within each K in {8,10,12}, fit FULL_decorr and report (a) within-K
  out-of-fold Spearman(predicted, r_sb_hp) via repeated 70/30 set holdout, (b) the FULL-minus-category
  delta, (c) the standardized coefficient of each feature and whether its SIGN is STABLE across K and vs
  real. The clean-outcome question: with outcome noise negligible, do the features predict AND do the
  signs stabilize (does F4_resid have a consistent sign)?
- REAL subject-split track (honest out-of-subject): fit FULL_decorr on the half-A reliabilities and
  predict the INDEPENDENT half-B reliabilities; report Spearman(predicted-from-A, actual-B) within K, and
  the reverse (B->A). Compare to the cross-half OUTCOME RELIABILITY ceiling from E3b (the max attainable).
  This removes the shared-subject inflation that made E4's real anchor an upper bound.
- Report the marginal (univariate) feature-reliability correlations alongside the partial coefficients, so
  the interpretation is unambiguous (the E4 review showed F4 is marginally POSITIVE but partially unstable
  under collinearity).

## 4. Exploratory verdict (NOT a gate; labeled exploratory)
Report a clear exploratory read, explicitly NOT the frozen H2 gate:
- MECHANISM-SUPPORTED-EXPLORATORY if, on the denoised outcomes with decorrelated features: FULL_decorr
  predicts within-K out-of-fold (sim high-precision Spearman >= 0.30 at K in {8,10,12}) AND beats
  BASE_category AND the decorrelated coefficient signs are stable across K and consistent between the sim
  high-precision and real subject-split tracks (in particular F1 and F4_resid). Interpretation: the frozen
  H2 failure was substantially artifact (noise + collinearity); a real feature -> reliability law exists;
  E5 may build designs from the FITTED empirical model (exploratory-derived, honestly labeled), NOT from a
  theoretical single-criterion information objective.
- MECHANISM-NOT-SUPPORTED otherwise. Interpretation: the mechanism does not hold even with a clean outcome
  and decorrelated features; Study 3 reframes to a transparent balanced/heuristic comparison (Option B).
- Either way this is EXPLORATORY and does NOT change the frozen H2 FAIL. The verdict is a recommendation to
  the author, not a pre-specified result.

## 5. Outputs (all new, under elevation/results/e4b/; additions-only)
- `elevation/analysis/e4b_mechanism_reanalysis.R`.
- `elevation/results/e4b/coefficients_decorrelated.csv` (FULL_decorr + drop-F4 standardized coeffs +
  marginal correlations, per track/target/K).
- `elevation/results/e4b/cv_denoised.csv` (within-K OOF Spearman on r_sb_hp; real A->B and B->A; vs the
  outcome-reliability ceiling; FULL-minus-category deltas).
- `elevation/results/e4b/sign_stability.csv` (decorrelated coefficient signs across K + sim-vs-real).
- `elevation/results/e4b/exploratory_verdict.json` (the exploratory read + its basis; explicitly
  `replaces_frozen_h2 = false`).
- `elevation/results/e4b/e4b_config.json` + seeds.
- the run report.

## 6. Pre-declared PASS gates (Part D)
- D1 inputs reused, not recomputed: reads only committed E3b + E2 + E4 CSVs; no raw/choices; feature
  re-aggregation sample matches E2.
- D2 build-only / no leakage: no participant-level data enters (set-level only); test sealed.
- D3 exploratory verdict computed + reported, explicitly labeled NOT the frozen H2 gate
  (`replaces_frozen_h2 = false`); the frozen H2 FAIL is restated as standing.
- D4 decorrelation done as specified (F4_resid = F4 orthogonalized on F2; drop-F4 robustness); marginal +
  partial reported together; models linear/pre-declared.
- D5 reproducibility: seeds committed; CV reproducible in distribution.
- D6 additions-only + frozen integrity (cryptographic): diff <E3b HEAD> -> HEAD touches ONLY {the session log,
  the coordination directory, elevation}; 13 frozen Study-1 subtrees byte-identical; e0/e1/e2/e3/e4/e3b subtrees IDENTICAL;
  all prior spec blobs unchanged; master + tags untouched; branch pushed.
- D7 no forbidden computation: NO optimization, NO candidate designs, NO freeze/test read; the frozen H2
  gate is not re-run or amended; simulation labeled.
- D8 hygiene: EXPLORATORY labeled throughout.

## STOP
After D1-D8: STOP and bring the exploratory verdict to the author. If MECHANISM-SUPPORTED-EXPLORATORY: the next
step is E5 authored on a fitted-model objective (the author’s decision). If MECHANISM-NOT-SUPPORTED: reframe Study 3
(Option B). Do NOT proceed to E5, open freeze (E6), or read test (E7) without the author’s ruling.

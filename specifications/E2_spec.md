# E2 problem-information dictionary - field specification (Study 2)

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md.
This is an analysis specification, not a new pre-specification: the Studies 2-3 contract was
already frozen at E1. It documents the reusable Table-1 problem-information dictionary
(execution spec, Table 1) whose fields feed H2's mechanism model at E4.

## 0. What E2 is (and is NOT)
- E2 builds ONE item-level table: one row per clean problem in the CPC18 clean-problem
  universe, carrying the predictor / structure / entropy / Fisher-information fields that
  later analyses aggregate into set-level mechanism predictors (E4) and design criteria (E5).
- Choice-dependent fields are computed on the 344 BUILD subjects ONLY (E1 split). The 137
  freeze and 205 test subjects' choices are NEVER read. The 205 test subjects stay sealed
  until E7 (E1 contract section 3; hard rule).
- E2 computes NO reliability, NO recovery, NO prediction, NO optimization, NO mechanism model.
  Those are E3 and later. E2 is descriptive only.
- E2 is FULLY DETERMINISTIC: no random number generation, no bootstrap, no seed. Every field
  is a closed-form function of the raw data restricted to build subjects. This makes the
  independent feature audit (the E2 gate) an exact re-derivation, not an expectation-level check.
- Targets, per the frozen E1 contract: loss_avoid_adj is PRIMARY, ev_sens_agn is SECONDARY,
  both metric scale; b_ev / b_risk / b_loss are carried on the RANK scale only (their metric
  ceiling is untrustworthy per E0). E2 provides the design Fisher atoms for all three
  coefficients (they are legitimate design quantities), but makes NO metric reliability claim
  for them; that policy is unchanged.

## 1. Inputs (build-only; additions-only)
- `data/raw/All estimation raw data.csv` - md5 must equal `38a26407a36c690be969f9d9d3fe3e5a`
  (assert before use, as E0/E1 did). 686 subjects, 210 problems, 25 trials per problem.
- `elevation/results/e1/participant_split.csv` - the frozen split. Filter to `split == "build"`
  (344 subjects) BEFORE reading any choice. Assert the build count is 344 and that no
  `split == "test"` (or `held_out_until_E7 == TRUE`) row is ever used for a choice-based field.
- Reuse the FROZEN loaders `R/cpc18_features.R` (`load_cpc18`, `is_clean_row`,
  `clean_problem_features`, `add_scaled_features`) so structure/predictor fields are identical
  to Study 1 and E0. No frozen Study-1 file is modified. All outputs live under `elevation/`.

## 2. Clean-problem universe (the rows)
- Clean = `LotNumA == 1 & LotNumB == 1 & Amb == 0` (frozen `is_clean_row`): two-outcome,
  no-ambiguity gambles. The universe is the set of clean GameIDs that at least one BUILD
  subject saw. Assert this equals the full clean set (expected: 105 problems; the build split
  is cohort-respecting across all 7 forms, so every clean problem is build-visible). Report and
  assert 0 clean problems unseen by build.
- Scaling constants for z_ev / z_sd / z_eloss = SD of d_ev / d_sd / d_eloss over this
  clean-problem universe (the frozen `add_scaled_features` convention). Because build sees the
  full clean set, these equal the Study-1 / E0 constants; record them in the config for audit.

## 3. Field dictionary (columns of `problem_information_dictionary.csv`, one row per clean problem)

### A. Identity and provenance
- `GameID` (int, unique key).
- `in_Set1` ... `in_Set7` (logical) - which of the 7 fixed forms include this clean problem
  (from raw: distinct `Set` values where the problem is clean and present). Supports
  leave-one-form-out (H6) and per-form stratification at E3/E7.
- `n_forms` (int) - number of Sets containing it.
- `n_build_subjects` (int) - distinct BUILD subjects who saw this problem.
- `n_build_trials_exp` (int) - build clean trials at Trial 6-25 that back the experience-regime
  entropy/proportion fields.
- `n_build_trials_desc` (int) - build clean trials at Trial 1-5 (description regime).

### B. Structure (description-based; subject-independent)
- Raw gamble parameters: `Ha, pHa, La, Hb, pHb, Lb`.
- Option moments (frozen `two_outcome_moments`): `ev_a, sd_a, eloss_a, has_loss_a` and
  `ev_b, sd_b, eloss_b, has_loss_b`, where for a two-outcome option (H w.p. p, else L):
  EV = pH + (1-p)L; SD = sqrt(max(0, pH^2 + (1-p)L^2 - EV^2));
  Eloss = p*max(-H,0) + (1-p)*max(-L,0); has_loss = (H < 0 or L < 0).
- `domain` in {gain, loss, mixed}: gain = no negative outcome in either option; loss = no
  positive outcome in either option; mixed = otherwise.
- `loss_informative` (logical) = `xor(has_loss_a, has_loss_b)` (loss in exactly one option;
  the loss_avoid construct's eligibility gate).
- `near_parity` (logical) = `abs(d_ev) <= median(abs(d_ev))` AND `abs(d_sd) >= median(abs(d_sd))`,
  medians taken over the clean-problem universe (reproduces Study-1 `near_parity_sensitivity`).
- `dominance` in {A, B, none}: first-order stochastic dominance of one option's two-outcome
  distribution over the other (deterministic); `none` if neither dominates. Supports the
  "no dominance by extreme items" optimization constraint (execution spec 10.3).

### C. Predictor (the difference features / model regressors)
- `d_ev = ev_b - ev_a`, `d_sd = sd_b - sd_a`, `d_eloss = eloss_b - eloss_a` (frozen
  `clean_problem_features`).
- `z_ev, z_sd, z_eloss` = d_* divided by their clean-universe SDs (frozen `add_scaled_features`);
  these are the exact per-subject Firth regressors:
  P(choose B) = plogis(b0 + b_ev*z_ev + b_risk*z_sd + b_loss*z_eloss).
- `abs_d_ev, abs_d_sd, abs_d_eloss` - magnitudes (stratification convenience).
- `d_ev_noloss` (loss-informative rows only; NA otherwise) = EV(no-loss option) - EV(loss option)
  = `ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a)`. The loss_avoid_adj predictor.

### D. Entropy (observed informativeness; BUILD choices; deterministic)
Experience regime (Trial 6-25) is primary (E1 contract section 4); description regime
(Trial 1-5) is a secondary generalization field.
- `p_choose_B_exp` = build empirical P(choose B) on this problem, pooled over build subjects'
  clean trials at Trial 6-25.
- `choice_entropy_exp` = binary entropy in bits = `-p*log2(p) - (1-p)*log2(1-p)` at
  p = p_choose_B_exp (0 when p in {0,1}). High entropy (p near 0.5) = maximally discriminating.
- `p_choose_noloss_exp` (loss-informative rows only; NA otherwise) = build empirical
  P(choose the no-loss option), experience regime.
- `p_choose_B_desc`, `choice_entropy_desc` - description-regime counterparts (Trial 1-5).

### E. Fisher information (design-expected information at a stable BUILD operating point)
Distinct from entropy: entropy is observed choice spread; Fisher is the design's EXPECTED
information at a population operating point (a design property usable prospectively, the
optimization-relevant quantity). Reference operating point = the build-only pooled logistic
theta_bar_build (section 4). For design row x = (1, z_ev, z_sd, z_eloss):
- `fisher_w_exp` = p_ref*(1 - p_ref), p_ref = plogis(x . theta_bar_build). The per-problem
  Fisher weight.
- `finfo_b_ev = fisher_w_exp * z_ev^2`, `finfo_b_risk = fisher_w_exp * z_sd^2`,
  `finfo_b_loss = fisher_w_exp * z_eloss^2` - the diagonal own-information each problem
  contributes about each coefficient. The full 4x4 per-problem contribution is
  fisher_w_exp * x x^T; any set's information matrix = sum over its problems of that, fully
  reconstructable from `fisher_w_exp` + the z_* columns (documented so E4 needs no extra field).
- `w_noloss_exp` (loss-informative rows only; NA otherwise) = p_nl*(1 - p_nl), p_nl =
  build loss-baseline prediction (section 4). The design-expected discriminating weight for the
  loss residual index.
- `info_loss_avoid_adj` (primary target's item information) = `loss_informative * w_noloss_exp`
  (NA/0 when not loss-informative). Documented as a design proxy for how much a problem informs
  loss_avoid_adj; reliability itself is measured empirically at E3/E7, never asserted from Fisher.
- Target crosswalk (no duplicated columns): the ev_sens_agn (secondary target) item information
  is `finfo_b_ev` (EV-slope information); the loss_avoid_adj (primary) item information is
  `info_loss_avoid_adj`. Stated here rather than stored twice.

## 4. Build reference operating point (build-only; deterministic; exported for audit)
Written to `elevation/results/e2/e2_build_reference.csv` and echoed in the config:
- `theta_bar_build` = coefficients of the pooled logistic `B ~ z_ev + z_sd + z_eloss` fit on
  ALL build clean trials in the experience regime (Trial 6-25), one Bernoulli row per trial.
  A single stable population-average operating point (pools ~344 subjects x clean problems x 20
  trials, so no separation, unlike the per-subject Firth coefficients E0 flagged). It is a
  reference for expected information ONLY; it is not a per-subject estimate and makes no
  reliability claim.
- Build loss baseline = coefficients of `chose_noloss ~ d_ev_noloss` fit on build
  loss-informative clean trials, experience regime (the build-only analogue of Study-1
  `fit_loss_baselines$exp`, used here only to weight the loss information atom; it does NOT
  replace the frozen Study-1 baseline or the frozen reliability engine).
- Report each fit's coefficients, n, and convergence. Both are deterministic (no seed).

## 5. Mechanism-predictor crosswalk (how E2 atoms feed H2 at E4)
H2 (E1 contract section 7) pre-declares five set-level information predictors. Each aggregates
E2 item fields over a candidate set; E2 supplies the atoms, E4 does the aggregation under
grouped cross-validation:
- target-predictor variance -> Var over the set of the target's predictor: z_ev (ev_sens_agn) /
  d_ev_noloss on loss-informative items (loss_avoid_adj).
- loss exposure -> count/fraction `loss_informative`, mean `eloss_a`/`eloss_b`, mean `d_eloss`.
- predictor collinearity -> correlations among z_ev / z_sd / z_eloss across the set.
- expected information -> functionals of sum over the set of fisher_w_exp * x x^T (e.g. the
  target-parameter information, or log-det); `finfo_*` are the diagonal atoms.
- choice entropy -> mean/sum of `choice_entropy_exp` over the set.
E2 is complete for H2 iff every one of these five is reconstructable from the dictionary; the
gate G8 checks that mapping.

## 6. Output files (all new, under `elevation/`; additions-only)
- `elevation/analysis/e2_problem_features.R` - the deterministic build script (written per
  this spec; reuses the frozen loaders).
- `elevation/results/e2/problem_information_dictionary.csv` - THE Table-1 dictionary (one row per
  clean problem; all fields in section 3).
- `elevation/results/e2/e2_build_reference.csv` - theta_bar_build + build loss baseline
  (coefficients, n, convergence).
- `elevation/results/e2/e2_universe_audit.csv` - universe checks: n_clean, n_build_subjects
  (344), n_freeze (137), n_test (205) as read from the split, count of clean problems unseen by
  build (must be 0), per-Set clean counts, scaling constants.
- `elevation/results/e2/e2_feature_config.json` - raw md5, regimes (exp 6-25 / desc 1-5),
  scaling constants, near_parity medians, build N, universe size, theta_bar_build + loss-baseline
  coefficients, input provenance (reads participant_split.csv build subset + raw CSV; asserts NO
  freeze/test choice read), determinism note (no RNG/seed).
- The run report. one appended session-log row.

## 7. Pre-declared review gates (independent feature audit; committed BEFORE results exist)
Verification grades by re-deriving from an independent copy of the data (generator != grader). Deterministic
quantities must match essentially exactly (float epsilon); the fitted reference matches to
solver tolerance.
- G1 Universe and shape: exactly 105 rows, unique GameID; the row GameIDs equal the clean set
  (LotNumA==1 & LotNumB==1 & Amb==0) re-derived from raw; every row has n_build_subjects >= 1;
  0 clean problems unseen by build.
- G2 Build-only / no leakage: build N used = 344 and equals the participant_split build set; NO
  freeze (137) or test (205) subject enters any choice-based column or either reference fit; the
  script filters split=="build" before any choice read; the config records this. Positive check:
  the independently recomputed p_choose_B_exp recomputed on the 344 build subjects matches; a whole-sample or
  test-inclusive recomputation would differ, so the match is evidence of correct restriction.
  The 205 test subjects remain sealed.
- G3 Structure/predictor exact: d_ev, d_sd, d_eloss, z_ev, z_sd, z_eloss, ev_*, sd_*, eloss_*,
  has_loss_*, loss_informative, domain, near_parity, dominance reproduced from raw for all 105
  (floats to < 1e-9; logicals/ints exact). Scaling constants = SD over the 105 clean set; match.
- G4 Entropy/proportion (build-only, deterministic): p_choose_B_exp and choice_entropy_exp
  reproduced on an independent 344-build restriction (Trial 6-25) to max abs diff < 1e-9 across all
  105; binary-bits entropy formula confirmed; p_choose_noloss_exp matches on loss-informative
  rows; description-regime fields likewise.
- G5 Fisher reference (build-only): theta_bar_build reproduced by an independent build-only
  pooled logistic (converged; coefficients match to < 1e-6); fisher_w_exp and finfo_* reproduced
  to < 1e-6; build loss baseline and info_loss_avoid_adj reproduced; neither fit uses freeze/test.
- G6 Additions-only + frozen integrity (cryptographic): branch HEAD advances from
  f6780ca9cfe9b696a98c4fc04658baf1fe6daea4; `git diff --name-only f6780ca HEAD` touches ONLY
  {the session log, the coordination directory, elevation}; all 13 frozen Study-1 top-level paths byte-identical;
  elevation/results/e0 AND elevation/results/e1 subtrees IDENTICAL (E0/E1 untouched);
  elevation/E1_contract.md byte-identical (blob 93de88d2); master (0808abc) + both tags untouched.
- G7 No forbidden computation / test sealed: no reliability, recovery, or prediction is computed;
  no output file contains any test-subject-level choice or any test outcome; the dictionary rows
  are per-problem (no subject-level column), so there is no subject leakage vector in the output.
- G8 Contract consistency + hygiene: targets match the E1 contract (loss_avoid_adj primary +
  ev_sens_agn secondary metric; b_* rank-only, no metric reliability claim); all five H2
  mechanism predictors (section 5) are reconstructable from the dictionary.

## 8. Confirmations recorded here (forward planning questions, from the E1 report)
- Q1 (E2 scope): YES. E2 restricts to the 344 build subjects and the clean-problem universe;
  its predictor / structure / entropy / Fisher-information fields are the Table-1 dictionary that
  feeds H2's mechanism model at E4, with freeze and test participants untouched. Matches the
  pre-declared E2 scope in the execution spec.
- Q2 (E7 evaluation): ACKNOWLEDGED, deferred to E7. Since the split is cohort-respecting and
  every form contributes ~30% to the test split, E7 may report BOTH pooled AND leave-one-form-out
  held-out evaluation (H6 cross-form generalization). E2 provides the per-form membership atoms
  (in_Set1..in_Set7) that make leave-one-form-out possible; the decision itself is an E7 matter,
  not an E2 action.

## STOP
E2 ends at the committed dictionary + reference + audit + config, the report, and the gate
evidence. Do NOT build the mechanism model (E4), candidate designs (E5), or read any freeze/test
outcome. E3 (same-count subset benchmark, the H1 premise gate) is the next analysis, authored by
verification after E2 review.

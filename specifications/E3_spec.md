# E3 same-count subset benchmark - specification (Study 2). The H1 premise gate.

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md. E3 is the premise gate for H1: at a fixed number of
unique clean problems K, does held-out reliability vary materially across item sets? If item-set
variance is trivial, H3 and H4 are moot and the program stops at exploratory (INTEGRITY.md: a failed
premise is a valid outcome, honestly reported).

## 0. What E3 is (and is NOT)
- E3 quantifies reliability heterogeneity across equal-K item sets, and reports recovery and
  prediction distributions alongside it, so E4 can fit the H2 mechanism model. The GATE is on
  reliability heterogeneity (H1) only.
- BUILD-ONLY. Every quantity uses the 344 BUILD subjects of the frozen E1 split. FREEZE (137) is
  reserved for E6 selection and is NOT read at E3. TEST (205) is sealed until E7. This is stricter
  than the E1-report Q3 guess ("build and freeze"): the E1 contract section 3 assigns all
  exploration to BUILD and reserves FREEZE for E6, so E3 reads BUILD only. Held-out evaluation at
  E3 is achieved by within-BUILD resampling (split-half over subjects; grouped folds over
  problems), never by opening freeze or test.
- COMPOSITION FRAME (approved before E3): a candidate design may compose problems FREELY across the
  105-problem clean universe; in_Set membership is a stratification/provenance variable, not a hard
  constraint (single-form administrability is a deployment caveat, noted, not a limit on the
  optimizer). The E2 finding that clean problems PARTITION across the 7 forms (n_forms == 1 for all
  105) means real-data split-half reliability of a FIXED set is only computable WITHIN a form (only
  a form's subjects saw that form's problems). E3 therefore measures reliability two ways:
  real-data WITHIN-FORM (Part A, the primary gate, anchored in observed choices) and SIMULATION
  FREE-COMPOSITION across the whole universe (Part C, the E0-determined validation path, which
  realizes the free-composition frame and supplies recovery). Part C is validated against Part A
  where both are defined (within Set 1).
- E3 does NO optimization, NO mechanism-model fitting (that is E4), and reads NO freeze/test choice.

## 1. Inputs (build-only; additions-only under elevation/)
- `data/raw/All estimation raw data.csv` (assert md5 `38a26407a36c690be969f9d9d3fe3e5a`).
- `elevation/results/e1/participant_split.csv` - filter `split == "build"` (344); assert 344 before
  any choice read; never read freeze/test choices.
- `elevation/results/e2/problem_information_dictionary.csv` - the item universe + features
  (loss_informative, in_Set*, z_*, d_ev_noloss, etc.). E3 selects subsets from this universe.
- Reuse the FROZEN engines unchanged (source them; do NOT redefine): `R/step2_reliability.R`
  (`build_subject_data`, `fit_loss_baselines`, `estimate_half`, `split_correlations`, `sb_correct`,
  `reliability_point`), `R/cpc18_features.R`, and `R/simulate_recovery.R` (`draw_true_params`,
  `analytic_p_b`, `simulate_choices`). Build-only loss baselines are fit on BUILD trials (as E2).

## 2. Targets and estimand (from the frozen E1 contract)
- PRIMARY target: loss_avoid_adj. SECONDARY: ev_sens_agn. Both metric scale (Study-1 definitions,
  frozen `estimate_half`). b_ev/b_risk/b_loss are NOT E3 targets (rank-only per E0/E1).
- Reliability of a K-set S = between-person split-half Spearman, Spearman-Brown corrected, the
  frozen Study-1 engine, experience regime (Trial 6-25, r = 20). For a FIXED set S the frozen
  `estimate_half` is called on S's problems, split into two halves of K/2; the target is estimated
  per half per subject and correlated across subjects (frozen `split_correlations` + `sb_correct`),
  averaged over N_split random halvings. Only the set is fixed; the index computation is identical
  to Study 1, so the grader re-derives with the same definitions.

## 3. Part A (PRIMARY): real-data within-form reliability heterogeneity + H1 gate
For each cell (form F in 1..7, budget K, target T):
- Universe = F's clean problems (from E2: Set1 24, Set2 16, Set3 17, Set4 14, Set5 14, Set6 9,
  Set7 11). Feasible K (K <= |F|, from the frozen grid {6,8,10,12}): Sets 1-5 all of {6,8,10,12};
  Set6 {6,8}; Set7 {6,8,10}.
- Sample M distinct K-subsets uniformly without replacement: M = min(2000, C(|F|,K)); if
  C(|F|,K) <= 2000 enumerate ALL (e.g. Set6 K=8 has C(9,8)=9). Seeds pre-declared (section 6).
- For each set S: r_sb(S) = frozen reliability over N_split = 100 random halvings across F's build
  subjects (all of whom saw all of F's clean problems). loss_avoid_adj: a split whose half lacks a
  loss-informative problem returns NA for that half (frozen behaviour) and is dropped; report the
  identifiable-split fraction (coverage) per set and per cell. ev_sens_agn: always identifiable.
- Heterogeneity + variance decomposition per cell:
  - sigma2_set_raw = Var over the M sets of r_sb(S).
  - Noise floor: on a pre-declared subsample of 50 sets, run the frozen subject-cluster bootstrap
    (`reliability_point`, N_boot = 1000) to get each set's estimation variance; sigma2_noise =
    mean of those. (This captures subject-sampling + split noise for a single set.)
  - sigma2_set_true = max(0, sigma2_set_raw - sigma2_noise); ICC_set = sigma2_set_true /
    sigma2_set_raw. Report sigma_set_true, ICC_set, and the r_sb distribution (min, quartiles, max,
    IQR) and the best-minus-median gap (foreshadows E5 headroom; NO optimization performed).
- PRIMARY GATE CELL: Set 1, K = 12, loss_avoid_adj (richest cell: 24 clean, 63 build subjects, 6
  loss-informative -> most K-sets identifiable). Contingency: if loss_avoid_adj identifiable-set
  coverage at this cell is < 50%, the primary gate target falls back to ev_sens_agn (pre-declared).
- H1 pre-declared PASS: at the primary gate cell, sigma_set_true >= 0.05 (SB units) AND ICC_set >=
  0.50. Report the same statistics for EVERY feasible cell (all forms x feasible K x both targets)
  as the full heterogeneity map. H1 is the premise gate: a FAIL is a valid scientific result and is
  reported honestly (the extension then stays exploratory per the E1 contract).

## 4. Part B: held-out choice-prediction distribution across sets (real data, build folds)
At the primary cells only (Set 1, K in {8, 12}; both targets' sets shared), to show the
which-problems gap in a validity metric (not just reliability), bounded for compute:
- For each of a pre-declared 300-set subsample: fit person-level coefficients on the K-set using
  the frozen per-subject logistic (`fit_subject_firth` from R/simulate_recovery.R) for each Set-1
  build subject, then predict that subject's choices on the Set-1 clean problems NOT in S (held-out
  WITHIN build, grouped by problem so no trained problem leaks), and record mean held-out log loss
  and Brier score. Report the across-set distribution and the best-minus-median gap.
- No freeze/test. Held-out = build subjects' own unseen-in-S problems, not a held-out cohort.

## 5. Part C: simulation free-composition reliability + recovery (E0 validation path)
Realizes the free-composition frame across the whole 105-universe and supplies recovery (which needs
planted truth). Uses the FROZEN generative model (`draw_true_params`, `simulate_choices`) that E0
validated; the parameter distribution is the frozen uniform-range population (E0), a controlled
population with fixed between-person variance - documented as SIMULATION, never an empirical result.
Calibrating the parameter distribution to build is a pre-declared E4/E5 step, not done here.
- For each K in {6,8,10,12}: sample N_sets_sim = 500 K-subsets FREELY from the 105 universe (seeds
  in section 6). Draw N_sim_subj = 300 simulated subjects (frozen `draw_true_params`).
- For each set S: simulate n_reps = 20 choices per problem per subject (`simulate_choices`), compute
  each target index per subject on each half of S (frozen index formulas), and:
  - simulated reliability(S) = split-half SB Spearman across the 300 sim subjects, N_split = 50.
  - recovery(S) = Spearman between the estimated index (from simulated choices) and the TRUE index
    (computed from the noise-free choice probabilities `analytic_p_b` under each subject's planted
    theta on S), across subjects; report correlation + RMSE-of-ranks.
- Report the across-set distributions of simulated reliability and recovery per K, and the joint
  (reliability, recovery) per set (E4 will model these against E2 features).
- ANCHOR CHECK (validates the simulation against real data): restrict Part C to sets drawn only from
  Set 1's 24 clean problems and confirm the simulated reliability heterogeneity agrees with Part A's
  real heterogeneity at Set 1 in magnitude and rank (report Spearman between real and simulated
  r_sb over the shared sets, and the two sigma_set values). This is a reported diagnostic, not a
  hard gate.

## 6. Pre-declared sampling and seeds (committed before results)
- Regime experience, pos = Trial 6-25 (r = 20). Targets loss_avoid_adj (primary), ev_sens_agn.
- Part A: M = min(2000, C(|F|,K)); N_split = 100; noise-floor subsample 50 sets, N_boot = 1000.
- Part B: 300-set subsample at Set 1, K in {8,12}; frozen `fit_subject_firth`.
- Part C: N_sets_sim = 500 per K; N_sim_subj = 300; n_reps = 20; N_split = 50.
- Seeds (deterministic, reproducible): base SEED_E3 = 20260719. Part A cell seed =
  SEED_E3 + 1000*form + K (subset draw) and + 500000 (split draw). Part B seed = SEED_E3 + 900000.
  Part C set-draw seed = SEED_E3 + 700000 + K; param-draw seed = SEED_E3 + 800000 + K. All seeds
  recorded in `e3_config.json` and in a per-cell seeds CSV, so the grader reproduces the exact set
  membership (deterministic given seed) and matches stochastic reliabilities at expectation.

## 7. Outputs (all new, under elevation/; additions-only)
- `elevation/analysis/e3_subset_benchmark.R` - the build script (reuses frozen engines).
- `elevation/results/e3/reliability_by_set_real.csv` - per (form, K, target, set_id): the sampled
  GameIDs, r_sb, coverage, n_subjects. (The per-set training rows for E4.)
- `elevation/results/e3/heterogeneity_summary.csv` - per cell: sigma_set_raw, sigma_noise,
  sigma_set_true, ICC_set, r_sb min/q25/median/q75/max/IQR, best_minus_median, coverage, M.
- `elevation/results/e3/prediction_by_set_real.csv` - Part B per-set held-out log loss + Brier.
- `elevation/results/e3/reliability_recovery_by_set_sim.csv` - Part C per (K, set_id): sim
  reliability, recovery corr + RMSE, sampled GameIDs.
- `elevation/results/e3/sim_real_anchor.csv` - the Set-1 anchor comparison.
- `elevation/results/e3/h1_gate.json` - the primary-cell verdict (sigma_set_true, ICC, thresholds,
  PASS/FAIL) + the full-map summary.
- `elevation/results/e3/e3_config.json` + `elevation/results/e3/e3_seeds.csv` - all constants/seeds.
- the run report.

## 8. Pre-declared PASS gates (Part D)
- D1 engine provenance: `e3_subset_benchmark.R` SOURCES the frozen `estimate_half`,
  `split_correlations`, `sb_correct`, `reliability_point`, `build_subject_data`,
  `fit_loss_baselines`, and the frozen generative functions; it does NOT redefine any of them
  (grep-level check + a numeric spot-check that one fixed set's per-split Spearman matches an
  independent recomputation).
- D2 build-only / no leakage: build N = 344; no freeze/test SubjID enters any frame; the split is
  filtered to build before any choice read; freeze and test choices are never read; 205 test sealed.
- D3 H1 gate computed and reported: sigma_set_true + ICC at the primary cell vs the pre-declared
  thresholds (0.05, 0.50); PASS/FAIL stated plainly; full heterogeneity map over all feasible cells;
  loss coverage reported; contingency applied if primary coverage < 50%.
- D4 distributions present: real reliability per set (Part A), held-out prediction per set (Part B),
  simulation reliability + recovery per set (Part C), and the Set-1 sim-vs-real anchor.
- D5 determinism of design + reproducibility: set membership is reproducible from the recorded
  seeds (deterministic); stochastic reliabilities are reproducible in distribution; seeds committed.
- D6 additions-only + frozen integrity (cryptographic): `git diff --name-only <E2 HEAD> HEAD`
  touches ONLY {the session log, the coordination directory, elevation}; all 13 frozen Study-1 subtrees byte-identical;
  `elevation/results/e0`, `e1`, and `e2` subtrees IDENTICAL (E0/E1/E2 untouched); `E1_contract.md`
  blob `93de88d2`, `E2_spec.md` unchanged; master + tags untouched; branch pushed.
- D7 no forbidden computation: NO optimization, NO mechanism-model fit (E4), NO reading of any
  freeze/test outcome. Simulation is labeled SIMULATION everywhere.
- D8 hygiene: simulation labeled.

## STOP
After D1-D8: STOP. E4 (mechanism model: predict held-out reliability from the E2 information features
under grouped cross-validation) is the next analysis, authored after E3 review. Do NOT fit
the mechanism model, construct candidate designs (E5), open freeze (E6), or read any test outcome
(E7). If the H1 premise FAILS, report it honestly and stop for a planning decision (the extension
stays exploratory per the E1 contract); do not repurpose the analysis to manufacture a pass.

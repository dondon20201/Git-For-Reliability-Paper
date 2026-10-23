# E4 mechanism model - specification (Study 2 core). The H2 gate.

A pre-specification was made before any result it might influence. Integrity is documented in INTEGRITY.md. The pre-specification for E4 H2 is locked in the file elevation/E1_contract.md. For E4 H2 I investigate if the features of previously declared information generalise to out-of-sample reliability. Do the features generalise under grouped cross validation beyond just raw count and problem type? If not then there is not support for the mechanism claim and Study 2's mechanism is labelled exploratory (which is also acceptable).

## 0. What E4 is (and is NOT)
- E4 fits a LOW-DIMENSIONAL, PRE-DECLARED linear model that predicts a problem set's reliability from
  five information features aggregated over the set. It is the main analysis; a regularized model and
  PLS are robustness checks only (roadmap 10.2). NO black-box learner is the primary.
- E4 is a MODELING step on COMMITTED outputs: it reads the E2 problem-information dictionary (per-problem
  features) and the E3 per-set reliabilities (real within-form and simulation free-composition). It reads
  NO raw data and NO participant choices, so it is build-only by construction and carries no leakage
  vector; the 205 test subjects are untouched (E7 sealed). It recomputes NO reliabilities (it consumes
  E3's) and fits NO new choice model.
- E4 does NO optimization and constructs NO candidate designs (that is E5). It only characterizes the
  feature -> reliability relationship and gates H2.
- Estimand note carried from E3 review: E3's pre-declared variance decomposition is CONSERVATIVE (the
  subject-cluster bootstrap subtracts subject noise that is common across sets; the covariance-based
  estimand gave a LARGER true variance). E4 therefore predicts per-set r_sb directly and, as a
  pre-declared robustness, ALSO predicts the covariance-based (subject-split) reliability so the H2
  verdict does not depend on the conservative decomposition.

## 1. Inputs (committed; build-only; additions-only)
- `elevation/results/e2/problem_information_dictionary.csv` - per-problem features (z_ev, z_sd, z_eloss,
  d_ev_noloss, eloss_a, eloss_b, loss_informative, choice_entropy_exp, finfo_b_ev, info_loss_avoid_adj).
- `elevation/results/e3/reliability_by_set_real.csv` - real within-form per-set rows
  (form, K, target, set_id, game_ids ";"-joined, n_loss_informative, n_subjects, r_sb, coverage).
- `elevation/results/e3/reliability_recovery_by_set_sim.csv` - simulation free-composition per-set rows
  (K, target, set_id, game_ids, sim r_sb, recovery), labeled SIMULATION.
- Optional provenance: assert the E2 dictionary is 105 rows and the E3 files parse; read NO raw choices.

## 2. Pre-declared features (per set S, per target T; aggregated from E2 over S's game_ids)
Standardize each feature (z-score) within the analysis sample before modeling. Definitions:
- F1 target_predictor_variance: loss_avoid_adj -> Var(d_ev_noloss) over the loss-informative problems in
  S (NA if < 2 loss-informative); ev_sens_agn -> Var(z_ev) over all problems in S. (How spread the
  target's own predictor is across the set.)
- F2 loss_exposure: n_loss_informative in S. This is the BROAD PROBLEM CATEGORY covariate for the loss
  target and the primary baseline predictor. Also record frac_loss = n_loss/K and mean_eloss =
  mean(eloss_a + eloss_b) over S (reported, not in the primary model).
- F3 predictor_collinearity: max over the three pairs of |Pearson corr| among {z_ev, z_sd, z_eloss}
  across S's problems (mean |corr| reported as robustness). High collinearity = features hard to separate.
- F4 expected_information: loss_avoid_adj -> sum over S of info_loss_avoid_adj; ev_sens_agn -> sum over S
  of finfo_b_ev. (The set's total Fisher information for the target parameter, at the E2 build operating
  point.)
- F5 choice_entropy: mean over S of choice_entropy_exp. (Average observed discriminating power.)
- Covariate K (budget) is included only in cross-K pooled fits; the PRIMARY H2 test is WITHIN-K (K
  constant), so "beyond raw count" holds automatically.

## 3. Pre-declared models (per target; linear, low-dimensional; NOT black-box)
- FULL: r_sb ~ F1 + F2 + F3 + F4 + F5 (within a K; add K when pooling across K).
- BASE_count: r_sb ~ 1 within-K (a constant; the within-K count baseline, since K is fixed) and r_sb ~ K
  pooled (the raw-count baseline).
- BASE_category: r_sb ~ F2 (loss_exposure) within-K [+ K pooled]. "Broad problem category" = loss count.
- Report standardized coefficients + signs for FULL (which features drive reliability). A ridge fit and a
  PLS fit are reported as ROBUSTNESS only (roadmap 10.2), never as the primary.

## 4. Cross-validation (grouped; cluster-aware; the H2 evidence)
Reliability estimates of overlapping sets are correlated through shared items, so ordinary SEs are
invalid; report out-of-fold PREDICTIVE performance with cluster-aware uncertainty (roadmap 10.2), not
p-values. Only identifiable sets (finite r_sb) enter.
- Primary metric: WITHIN-K out-of-fold Spearman(predicted r_sb, actual r_sb). Within-K makes "beyond raw
  count" automatic; comparison to BASE_category tests "beyond broad category".
- CV scheme: repeated random SET-level holdout within each (track, target, K): train FULL / BASE_category
  on a pre-declared 70% of sets, predict the held-out 30%, record out-of-fold Spearman and R^2; repeat
  over pre-declared folds (seed committed). Report the mean and a cluster bootstrap CI where the bootstrap
  resamples PROBLEMS (not sets) and propagates to the sets that contain them (the cluster-aware interval).
- Stricter problem-disjoint robustness (reported, feasibility-flagged): partition the track's problem pool
  into two halves; where enough sampled sets lie entirely within a half, evaluate a model trained on
  one half's sets on the other half's sets (fully item-disjoint). Report the out-of-fold Spearman and the
  usable-set count; if underpowered on a track, say so rather than omit.
- Tracks:
  - Part A REAL anchor (Set 1, empirical): real within-form r_sb, K in {6,8,10,12} (K=12 included for the
    set-level CV; the item-disjoint robustness is limited at high K in 24 problems and is flagged).
  - Part B SIMULATION free-composition (primary for range + the E5-relevant regime): sim r_sb over the
    105-universe, K in {6,8,10,12}. Labeled SIMULATION throughout.

## 5. Pre-declared H2 gate (evaluated build-only; committed before results)
H2 PASSES iff ALL hold:
- G-sim (primary): in the SIMULATION free-composition track, the FULL model's within-K out-of-fold
  Spearman(predicted, actual sim r_sb) >= 0.30 for the PRIMARY target loss_avoid_adj at K in {8,10,12}
  (the well-identified budgets), pooled with a cluster-bootstrap CI whose lower bound > 0.
- G-beyond: the FULL model beats BASE_category out-of-fold on the same metric by a positive margin
  (Delta out-of-fold Spearman > 0, or Delta R^2 > 0) with a cluster bootstrap interval of at least 90%
  supporting the improvement, for loss_avoid_adj at K in {8,10,12} in the simulation track. (Features add
  value beyond loss count alone.)
- G-real (empirical anchor): in the REAL Set-1 track, the FULL model's within-K out-of-fold
  Spearman(predicted, actual real r_sb) > 0 and beats BASE_category, with FULL coefficient SIGNS
  consistent with the simulation track for the features that carry the simulation model.
Secondary target ev_sens_agn is reported analogously (informational, not gating).
- If G-sim or G-beyond FAILS: H2 is NOT supported; the mechanism is reported as EXPLORATORY and E5's
  feature-based optimization has no validated basis. Bring the decision to the author; do not repurpose the
  analysis to pass.
- If G-sim + G-beyond pass but G-real FAILS: the mechanism holds under simulation but is NOT confirmed on
  real data; flag as a limitation (simulation-only mechanism), and E5 proceeds with that caveat explicit.

## 6. Outputs (all new, under elevation/; additions-only)
- `elevation/analysis/e4_mechanism_model.R` - the modeling script (reads committed E2 + E3 outputs only).
- `elevation/results/e4/set_features.csv` - per set (track, target, K, set_id): the five features + K +
  n_loss_informative + r_sb (the modeling table).
- `elevation/results/e4/model_coefficients.csv` - FULL / BASE_category / ridge / PLS standardized
  coefficients per (track, target, K and pooled).
- `elevation/results/e4/cv_performance.csv` - per (track, target, K, model): out-of-fold Spearman + R^2 +
  cluster-bootstrap CI; the FULL-minus-BASE_category deltas.
- `elevation/results/e4/problem_disjoint_robustness.csv` - the item-disjoint check + usable-set counts.
- `elevation/results/e4/h2_gate.json` - G-sim / G-beyond / G-real verdicts vs thresholds; H2 PASS/FAIL.
- `elevation/results/e4/e4_config.json` + seeds.
- the run report.

## 7. Pre-declared PASS gates (Part D)
- D1 inputs reused, not recomputed: E4 reads ONLY committed E2 + E3 CSVs; it recomputes NO reliability and
  reads NO raw participant choices (grep: no raw-CSV read, no step2 engine call for new reliabilities).
  Feature aggregation from E2 is deterministic and re-derivable.
- D2 build-only / no leakage: no participant-level data enters; the 205 test subjects are untouched; E4
  operates on set-level features + set-level reliabilities only.
- D3 H2 gate computed and reported: G-sim, G-beyond, G-real vs the pre-declared thresholds; PASS/FAIL
  stated plainly; the full per-K, per-target, per-track table delivered; failure reported honestly.
- D4 model is the pre-declared low-dimensional linear one; ridge/PLS are labeled robustness; standardized
  coefficients + signs reported; features match section 2 exactly (re-derivable from E2).
- D5 reproducibility: feature aggregation exact from E2; CV folds reproducible from committed seeds;
  out-of-fold metrics reproducible in distribution.
- D6 additions-only + frozen integrity (cryptographic): diff <E3 HEAD> -> HEAD touches ONLY {the session log,
  the coordination directory, elevation}; 13 frozen Study-1 subtrees byte-identical; e0/e1/e2/e3 subtrees IDENTICAL
  (untouched); E1_contract/E2_spec/E3_spec blobs unchanged; master + tags untouched; branch pushed.
- D7 no forbidden computation: NO optimization, NO candidate-design construction (E5), NO freeze/test
  read; simulation track labeled SIMULATION.
- D8 hygiene: simulation labeled.

## STOP
After D1-D8: STOP. E5 (candidate design construction: random, stratified, balanced, D-/A-/c-optimal,
robust average-case, multi-objective, at equal search budget) is the next analysis, authored
after E4 review, and only if H2 is supported (else bring the exploratory-path decision to the author first). Do
NOT construct designs, open freeze (E6), or read any test outcome (E7). If H2 FAILS, report honestly and
STOP for a planning decision.

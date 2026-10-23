# E8 hierarchical re-estimation - specification (Study 2 third act). Does partial pooling lower the wall?

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md.

E8 tests one pre-specified question for the BUILD group. Specifically it tests if the split-half reliability of the two trustworthy target indices increases if the two trustworthy target indices are re-estimated using a hierarchical (partial-pooling) estimator versus the previously frozen per-subject (no-pooling) estimator used in Study 1. This question tests the REDUCIBLE component of measurement noise present in the E0 reliability analysis and not the population true-variance ceiling. This question is a re-analysis of previously collected build data using a different estimator. No designs are created. No mechanism model is fit. No test or freeze subjects are read.

## 0. What E8 is (and is NOT)
- E8 re-estimates loss_avoid_adj (PRIMARY) and ev_sens_agn (SECONDARY) with a multilevel partial-pooling
  model, then re-runs the FROZEN reliability estimand (R/step2_reliability.R) on the pooled estimates and
  compares, split-for-split, against the frozen per-subject estimator. The ONLY thing that differs between
  the two arms is the per-half estimator; the split draws, budgets, regime, Spearman-Brown correction, and
  bootstrap are identical.
- E8 is BUILD-ONLY (344 subjects). The 137 freeze subjects are not read; the 205 test subjects stay sealed
  (the E7 rule: never read before a deliberate confirmatory step). This is a re-analysis, not a selection
  step, so freeze is not opened.
- E8 does NO optimization, constructs NO candidate designs (the retired E5), fits NO information-mechanism
  model (that was E4 / E4b; the frozen H2 FAIL stands untouched), and reads NO problems outside the frozen
  clean-subset definition.
- E8 does not reuse the reserved E5-E7 pre-specification slots. Those named the confirmatory optimization
  capstone, which was retired in its confirmatory form. E8 is a distinct, newly pre-declared
  re-analysis, numbered outside E5-E7.
- Ceilings bind (E0, LOCKED): no estimator can exceed the population true-variance ceiling (loss_avoid_adj
  reliability 0.5446 at 14 problems, growing to 0.70 at ~27; ev_sens_agn 0.2811 at 14). "Raises
  reliability" means moving the ACHIEVED reliability up toward that ceiling by removing estimation noise,
  never past it. An exceedance is a leakage flag, not a result.

## 1. Inputs (committed; build-only; additions-only)
- Raw CPC18 clean-subset trials (md5 38a26407a36c690be969f9d9d3fe3e5a asserted in-script), restricted to
  split == "build" via elevation/results/e1/participant_split.csv (assert 344 / 137 / 205; build filtered
  BEFORE any choice is read; choice columns dropped from the non-build frame; rm(raw)).
- elevation/results/e1/participant_split.csv - the frozen three-way split (seed 20260718).
- Frozen engines SOURCED, never redefined: R/step2_reliability.R (build_subject_data, fit_loss_baselines,
  run_cell, estimate_half, split_correlations, sb_correct, reliability_point), R/step1_measures.R,
  R/cpc18_features.R. The clean-subset and feature construction is the frozen Study-1 pipeline.
- E0 ceiling anchors, read-only, from elevation/results/e0/variance_ceiling_final.csv, used ONLY as the HG1
  guard reference (never refit): loss_avoid_adj ceiling_fulldata 0.5446 (implied single-problem reliability
  r1 = 0.0787); ev_sens_agn 0.2811 (r1 = 0.02717). The k_0.7 / k_0.8 columns (loss 27.31 / 46.82) reproduce
  these r1 values exactly under Spearman-Brown in problem count.

## 2. Targets and the two estimators (the only thing that changes between arms)
Both indices are the frozen Study-1 constructs. For each budget and split a subject's clean problems are
partitioned into two disjoint halves (the frozen scheme). The RAW arm estimates each half per subject
independently (frozen estimate_half). The HIER arm fits ONE multilevel model across ALL eligible subjects
for that half and returns each subject's partial-pooled estimate. The two halves are estimated in SEPARATE
fits sharing no data (anti-leakage, HG3).

- loss_avoid_adj (PRIMARY).
  - RAW (frozen): per-subject mean residual (chose_noloss - pooled-baseline p-hat) over the half's
    loss-informative trials.
  - HIER: multilevel logistic regression on the half's loss-informative trials, aggregated to per
    subject-by-problem binomial counts of chose_noloss out of the r sampled experience trials:
    chose_noloss ~ 1 + d_ev_noloss + (1 | SubjID). The fixed intercept + d_ev_noloss slope are the
    partial-pooling analogue of the frozen pooled EV baseline; the subject RANDOM INTERCEPT (conditional
    mode / posterior mean) is the partial-pooled loss-avoidance index. Per-subject estimate for reliability
    = that random intercept. Fit separately on each half.
- ev_sens_agn (SECONDARY; informational, not gating).
  - RAW (frozen): per-subject Spearman(choice, d_ev) over the half's clean trials.
  - HIER: multilevel logistic regression on the half's clean trials, aggregated to subject-by-problem
    binomial counts of choice B: B ~ 1 + z_ev + (1 + z_ev | SubjID). The subject RANDOM SLOPE on z_ev is
    the partial-pooled EV-sensitivity. LOCKED caveat: the hierarchical operationalization is a logistic
    slope, a monotone re-scaling of the frozen Spearman index; because reliability is rank-based across
    subjects the scale is immaterial, but the two are not the identical statistic, so ev_sens_agn is
    reported as secondary robustness only.
- Estimators (LOCKED):
  - Workhorse across the full budget grid: frequentist partial pooling via lme4::glmer (empirical-Bayes
    conditional modes), family binomial, weights = trials-per-cell, optimizer bobyqa. A fit that fails to
    converge is refit once with a longer iteration budget, then that subject-half is treated as the frozen
    engine treats a missing estimate (NA, excluded from that split's complete pairs, counted). Report the
    convergence-failure rate per budget.
  - Bayesian confirmation on the primary cell only (loss_avoid_adj, reference budget K = 12):
    rstanarm::stan_glmer or brms with weakly-informative priors (fixed effects Normal(0, 5) on the logit
    scale; random-intercept SD half-Normal(0, 2)); at least 20 splits; report Rhat, bulk-ESS, divergence
    count. Purpose: confirm the frequentist partial-pooling gain is not an artifact of the Laplace
    approximation, and deliver the fully Bayesian view the literature (Haines 2020) uses.

## 3. Estimand and the paired comparison (frozen engine; committed before results)
- Estimand (LOCKED, unchanged from Study 1 / E1): between-person split-half reliability, Spearman,
  Spearman-Brown corrected (frozen split_correlations + sb_correct + reliability_point), experience regime
  (Trial 6-25, r = 20), min_pairs = 10, subject-cluster percentile bootstrap CI (bootstrap outside, splits
  inside).
- Budgets B (problems per half): {6, 8, 10, 12}; unstratified within-subject partition (lambda = NA),
  matching E3's free within-form partition; the frozen run_cell eligibility (subject needs >= 2B clean
  problems) selects the eligible BUILD subjects per budget (K = 12 is effectively Set-1 build; smaller B
  pool more forms). Report eligible-N and coverage per budget.
- Pairing (LOCKED): a SINGLE split loop with one committed seed (SEED_E8) draws h1 / h2 / pos once per
  split; BOTH arms are evaluated on those identical indices, so the per-split difference isolates the
  estimator. N_SPLITS_HIER = 100 (fixed; the frozen curve's 200 is halved only to bound glmer runtime, and
  the RAW arm is recomputed here at the same 100 splits so the comparison is exactly paired; the frozen
  full-200 raw numbers are cited only as context).
- Primary reported quantity per (target, budget): rel_raw = SB reliability (raw arm); rel_hier = SB
  reliability (hier arm); Delta = rel_hier - rel_raw, with a subject-cluster bootstrap CI on Delta
  (resample subjects, recompute BOTH arms' per-split Spearman on the resampled subjects, average over
  splits, difference the two arms within each bootstrap draw so the CI is on the paired difference). Also
  gap_closed = (rel_hier - rel_raw) / (rel_ceiling(2B) - rel_raw), where rel_ceiling(k) = k * r1 / (1 +
  (k - 1) * r1) from the E0 anchors (loss r1 = 0.0787, ev r1 = 0.02717); report gap_closed as NA / "at
  ceiling" when rel_raw >= rel_ceiling(2B).

## 4. Pre-declared gates (evaluated build-only; committed before results)
Primary scientific gate:
- G1 (PRIMARY, loss_avoid_adj): at AT LEAST ONE well-identified budget B in {8, 10, 12}, Delta = rel_hier -
  rel_raw >= 0.05 (SB units) with a >= 90% subject-cluster bootstrap interval on Delta excluding 0. K = 6 is
  reported but flagged thin (E3b: loss_avoid_adj at K = 6 has ~34% usable halvings) and does not by itself
  satisfy or fail G1.
- G1-shape (REPORTED, not gated): the borrowing-strength prediction is that Delta is larger at smaller B
  (raw furthest below ceiling where per-subject data is thinnest) and shrinks toward 0 as raw approaches the
  ceiling at B = 12. Report the shape; do not gate on it.
- G2 (SECONDARY, ev_sens_agn): the same Delta and CI, reported at all budgets; informational, not gating.

Honesty guards (LOCKED; hard; a violation blocks the positive reading and is brought to the author, never worked
around):
- HG1 ceiling guard: rel_hier must not exceed rel_ceiling(2B) beyond bootstrap noise at any budget; a
  clear, CI-supported exceedance is a shrinkage-leakage flag (echoing E3b, where E4's real anchor exceeded
  its honest ceiling), reported as such, never as a gain.
- HG2 null-permutation control: with subject identities permuted between the two half-estimate vectors
  (breaking the true correspondence), rel_hier must fall to ~0 (within Monte-Carlo noise of 0) at every
  budget. A materially non-zero permuted reliability means partial pooling manufactures cross-half
  correlation and voids the hier arm. This is the E0 tau = 0 recovery logic applied to the estimator.
- HG3 no cross-half leakage: the two halves are fit in SEPARATE model calls on DISJOINT problems; assert
  in-script that no subject-half fit sees the other half's rows and that hyper-parameters are estimated
  within-half.
- HG4 shrinkage diagnostics: per budget, report the estimated random-effect SD and the mean shrinkage
  factor (ratio of pooled to raw between-subject SD); a reliability gain must coincide with sensible,
  non-degenerate shrinkage (report that the pooled between-subject variance did not collapse to a constant).
- Convergence: glmer convergence-failure rate reported per budget (target < 5%; handled as in section 2).
  brms confirmation: Rhat < 1.01, 0 divergences, bulk-ESS > 400 on the reported parameters.

Bayesian confirmation gate:
- G3 (primary cell, loss_avoid_adj, B = 12): sign(Delta_brms) == sign(Delta_glmer) and |Delta_brms -
  Delta_glmer| < 0.05. Confirms the workhorse.

Interpretation (LOCKED, so neither outcome is a moved goalpost):
- If G1 holds: hierarchical re-estimation recovers reducible measurement noise; the paper gains a
  constructive third act (a partial remedy to Study 1's measurement-noise component), reported with the
  honest caveat that the population true-variance ceiling is untouched.
- If G1 fails (Delta ~ 0 everywhere): the frozen per-subject indices are already near-efficient at CPC18
  budgets; the wall is a true-variance ceiling, not an estimation artifact. This is an equally publishable,
  arguably sharper, result for Study 1 and is reported as such. Either way E8 is a valid scientific outcome
  (INTEGRITY.md). Do NOT repurpose the analysis to pass.

## 5. Outputs (all new, under elevation/; additions-only)
- elevation/analysis/e8_hierarchical_reestimation.R - the analysis script (sources the frozen engines;
  reads only build data + the frozen split + the E0 ceiling anchors).
- elevation/results/e8/reliability_by_budget.csv - per (target, budget): rel_raw, rel_hier, Delta, bootstrap
  CIs, gap_closed, rel_ceiling_2B, eligible_N, coverage, n_splits.
- elevation/results/e8/delta_bootstrap.csv - the bootstrap distribution summaries for Delta per (target,
  budget).
- elevation/results/e8/null_permutation.csv - HG2 permuted rel_hier per (target, budget).
- elevation/results/e8/shrinkage_diagnostics.csv - HG4 random-effect SD + shrinkage factor + glmer
  convergence-failure rate per (target, budget).
- elevation/results/e8/ceiling_guard.csv - HG1: rel_hier vs rel_ceiling(2B) per (target, budget), with the
  exceedance flag.
- elevation/results/e8/brms_confirmation.csv - G3: Delta_glmer vs Delta_brms at the primary cell + Rhat /
  ESS / divergences.
- elevation/results/e8/per_subject_estimates_sample.csv - for a committed sample of splits, the paired
  per-subject (raw, hier) estimates for BOTH halves, so independent verification re-derives the cross-half Spearman from
  committed values.
- elevation/results/e8/g1_gate.json - G1 / G2 / G3 verdicts vs thresholds; HG1-HG4 flags; PASS/FAIL stated
  plainly.
- elevation/results/e8/e8_config.json + elevation/results/e8/e8_seeds.csv.
- elevation/logs/e8_run.log; the run report; the session log.

## 6. Exports for independent re-derivation (generator != grader)
The run EXPORTS enough committed per-unit values that independent verification re-derives the headline without rerunning R: for a
committed subset of splits, the per-subject (raw, hier) estimates on each half
(per_subject_estimates_sample.csv) so independent verification recomputes the cross-half Spearman, SB, and Delta from the
committed numbers; the full per (target, budget) rel_raw / rel_hier / Delta table; the null-permutation and
ceiling-guard tables. All committed CSVs at full round-trippable precision (>= 15 significant digits); note
the readr float-parse caveat (use utils::read.csv for bit checks) that E3b documented.

## 7. Pre-declared PASS gates (Part D)
- D1 engine provenance: R/step2_reliability.R + R/step1_measures.R + R/cpc18_features.R SOURCED, ZERO
  redefinitions (grep the frozen function names); the frozen reliability_point / split_correlations /
  sb_correct produce EVERY reported reliability in BOTH arms (the hier arm changes only the per-half
  per-subject input vector, not the reliability computation); a G_EQ check re-proves that on the RAW arm the
  new script reproduces the frozen estimate_half bit-for-bit on a committed sample of halves before any
  result is produced.
- D2 build-only / no leakage: split == "build" filtered before any choice read; choice columns dropped from
  the non-build frame; rm(raw); split counts asserted 344 / 137 / 205; freeze not read; the 205 test sealed.
  Both arms use only build.
- D3 estimator equivalence on the raw arm + hier definition: the raw arm equals the frozen estimator
  (G_EQ); the hier estimators are exactly the two models in section 2 (no extra covariates; no form fixed
  effect, matching the frozen index construction); the halves share no data (HG3).
- D4 outputs present: every file in section 5, with the paired per-subject sample for re-derivation.
- D5 reproducibility: all seeds committed (e8_seeds.csv); the split draws and both arms reproduce in
  distribution; the committed per-subject sample re-derives the reported Delta exactly.
- D6 additions-only + frozen integrity (CRYPTOGRAPHIC): git diff <pre-039 HEAD> -> HEAD touches ONLY
  {the session log, the coordination directory, elevation}; the 13 frozen Study-1 subtrees object-hash IDENTICAL to freeze
  387f9adc; elevation/results e0 446b142047e9 / e1 0e020fb0e7ea / e2 89c7d2f0d0a1 / e3 36120997f9f8 /
  e4 307fc9c6b49f and e3b all IDENTICAL to the pre-039 tip (E0-E4 + E3b untouched; e8 the only results
  addition); E1_contract / E2_spec / E3_spec / E3b_spec / E4_spec / E4b_spec blobs unchanged; master
  0808abc + both tags untouched; branch pushed.
- D7 no forbidden computation: NO optimization, NO candidate-design construction, NO information-mechanism
  model, NO freeze/test read; the frozen H2 gate (e4) neither re-run nor amended; no reserved E5-E7 slot
  reused.
- D8 hygiene: the pre-declared status labeled honestly (E8 is a
  newly pre-declared re-analysis, not part of the frozen E1 confirmatory arc).

## STOP
After D1-D8: STOP. Report rel_raw / rel_hier / Delta per (target, budget), the honesty-guard flags, and the
plain G1 verdict. Do NOT construct designs, open freeze, or read any test outcome; the 205 test subjects
remain sealed. Whatever the Delta, the frozen Study 1 (historical freeze tag) and the frozen H2 FAIL
stand untouched. Verification proceeds against an independently prepared answer key, re-deriving the paired
reliabilities before grading.

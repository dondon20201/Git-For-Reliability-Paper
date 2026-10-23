# Battery B - Study 2 robustness battery: heterogeneity stability, optimism, leave-one-feature-out (SPEC; committed BEFORE any result)

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md. Companions: E3_spec.md (H1), E4_spec.md (H2),
E3b_spec.md (subject-split outcome), E8_spec.md (robustness-analysis template).

Battery B is a THREE-part robustness battery on the frozen Study 2 result. It re-analyses BUILD data only; it
constructs no designs, opens no freeze, and reads no test subject (the 205 test group stays sealed). It NEVER
overwrites the frozen pre-specification: the E3 H1 PASS (same-count heterogeneity real) and the E4 H2 FAIL
(information mechanism not supported) stand untouched. A robustness pass strengthens them; a robustness surprise is
reported EXPLORATORY and cannot re-write a frozen gate (INTEGRITY.md: a pre-declared gate returning a null is a valid,
successful analysis and must not be repurposed to pass). This spec is committed before results; its gates are not
amended after results.

## 0. What Battery B is (and is NOT)
- B1 STABILITY re-estimates the H1 between-set true-variance decomposition (sigma_set_true, ICC_set) at the E3
  primary cell and across the feasible map, under multiple Monte-Carlo set samples, TWO noise-floor estimands, and
  many subject partitions, plus a null-permutation control. It tests whether H1 is stable to analyst choices E3 fixed
  at single values.
- B2 OPTIMISM measures the winner's curse in the E3 best-minus-median gap: it selects high-reliability equal-length
  sets on one subject half and evaluates them on a disjoint subject half, quantifying how much of the in-sample
  advantage transfers out-of-sample.
- B3 LOFO re-runs the frozen E4 grouped cross-validation dropping each information feature in turn, to test whether
  the H2 null depends on the multi-feature structure, whether any single-feature removal rescues a clean real signal,
  and whether expected information (F4) is redundant with loss exposure (F2).
- Battery B is BUILD-ONLY (344 subjects; the Set-1 cells use the 63 Set-1 build subjects). FREEZE (137) is not read;
  TEST (205) stays sealed (the E7 rule; this is a re-analysis, not a selection step, so freeze is not opened).
- Battery B does NO optimization, constructs NO candidate designs (retired E5), fits NO new choice model, and reads
  NO problems outside the frozen clean-subset definition. B1/B2 recompute reliabilities from raw build choices with
  the FROZEN engine; B3 consumes committed E2/E3/E4 outputs and recomputes NO reliability.
- Ceilings bind (E0, LOCKED): no reliability may exceed the E0 true-variance ceiling beyond bootstrap noise; an
  exceedance is a leakage flag, not a result.

## 1. Inputs (committed; build-only; additions-only under elevation/)
- `data/raw/All estimation raw data.csv` (assert md5 `38a26407a36c690be969f9d9d3fe3e5a` in-script), restricted to
  `split == "build"` via `elevation/results/e1/participant_split.csv` (assert 344/137/205; build filtered BEFORE any
  choice read; choice columns dropped from the non-build frame; rm(raw)). Needed by B1 and B2.
- `elevation/results/e1/participant_split.csv` - the frozen three-way split (seed 20260718).
- `elevation/results/e2/problem_information_dictionary.csv` - the 105-problem item universe + features (F1-F5 source
  columns z_ev, z_sd, z_eloss, d_ev_noloss, loss_informative, finfo_b_ev, info_loss_avoid_adj, choice_entropy_exp).
- `elevation/results/e3/e3_seeds.csv` + `elevation/results/e3/reliability_by_set_real.csv` - the committed Set-1
  set membership (reproduce the exact M-set draw) and the E3 per-set r_sb (B1 reference; B2 in-sample reference).
- `elevation/results/e3b/reliability_subject_split_real.csv` - the committed subject-split per-set r_sb_A / r_sb_B
  (B2 reference / anchor).
- `elevation/results/e4/set_features.csv`, `elevation/results/e4/cv_performance.csv`,
  `elevation/results/e4/coefficient_sign_consistency.csv`, `elevation/results/e4/e4_seeds.csv` - the B3 modeling table,
  the FULL-model anchors, and the committed E4 CV seeds so LOFO is PAIRED to the frozen FULL fit.
- Frozen engines SOURCED, never redefined: `R/step2_reliability.R` (build_subject_data, fit_loss_baselines,
  estimate_half, split_correlations, sb_correct, reliability_point), `R/cpc18_features.R`. Build-only loss baselines
  are the committed E3 constants (exp intercept 0.4301295476692, slope 0.165107793793507).

## 2. The three probes

### B1 - heterogeneity stability
Reuse the E3 Part-A machinery. PRIMARY cell = Set 1, K = 12, loss_avoid_adj (secondary ev_sens_agn reported). For each
of R_SEED independent Monte-Carlo set samples (M = min(2000, C(24,K)) distinct K-subsets, fresh seeds), and each K in
{6,8,10,12}, compute r_sb(S) over N_SPLIT_B1 common halvings on the 63 Set-1 build subjects (frozen estimate_half;
loss halves lacking a loss-informative problem return NA and are dropped, coverage reported). Decompose the across-set
spread under TWO pre-declared estimands:
- (i) bootstrap-floor (E3 primary): sigma2_noise = mean over a 50-set subsample of the frozen reliability_point
  percentile-CI normal-approximation variance; sigma2_set_true = max(0, sigma2_set_raw - sigma2_noise).
- (ii) covariance-based (subject-split): partition the 63 subjects into disjoint halves (~32/31), estimate each set's
  r_sb independently on each half with independent halvings, and take sigma2_set_true_cov = between-set Cov(r_sb_A,
  r_sb_B); ICC_cov = Cov / mean(Var_A, Var_B). Repeat over B_PART_B1 >= 200 random subject partitions; report the
  distribution and a 90% interval of sigma_set_true_cov.
Full-map (REPORTED, not gated): the count of feasible cells (all forms x feasible K x both targets) clearing both
thresholds under the covariance estimand, vs E3's 6/50 under the bootstrap floor. Null control: permute the set<->set
correspondence between the two subject halves; recover sigma_set_true (must collapse to ~0).

### B2 - optimism / winner's curse
Reuse the E3b subject-split machinery over MANY partitions. For K in {8,10,12} (K=12 primary), partition the 63 Set-1
build subjects into a SELECTION half A and EVALUATION half B (disjoint ~32/31), over B_PART_B2 >= 200 pre-declared
partitions. For the same M = min(2000, C(24,K)) E3 sets, compute r_sb^A and r_sb^B (frozen engine, independent
halvings per half; identifiable sets only). Per partition and averaged over both selection directions, report:
- single-best out-of-sample gain (vs median) = r_sb^B(argmax_A r_sb^A) - median_S r_sb^B (winner's-curse target);
- top-decile out-of-sample gain (vs grand mean, null-calibrated) = mean over the top-A-decile of r_sb^B - mean_S r_sb^B
  (the PRIMARY, stable, null-calibrated selectability statistic);
- optimism = in-sample apparent gain - out-of-sample gain; shrinkage ratio = out-of-sample / in-sample.
Reference the naive E3 best-minus-median (full 63) as the in-sample single-sample number. Null control: permute the
evaluation-half set labels; both gains must collapse to ~0.

### B3 - leave-one-feature-out on the frozen E4 mechanism model
Reuse the frozen E4 grouped-CV EXACTLY (200 repeats, 70/30 set-level holdout within (track, target, K), committed E4
seeds so folds match the FULL fit; features z-scored within the training fold; gate = min(oof_spearman_pooled,
oof_spearman_repeatmean)). Tracks real_form1 + sim_free105; target loss_avoid_adj (primary), ev_sens_agn (secondary);
K in {8,10,12} + pooled. For the FULL 5-feature model and each of the five leave-one-out variants:
- recompute the within-K out-of-fold Spearman and the drop-one delta = OOF(FULL) - OOF(drop f);
- report the retained features' standardized coefficient signs, and re-run the frozen H2 three-clause gate (G-sim,
  G-beyond, G-real sign-consistency) for the variant;
- additionally report the joint drop-{F2,F4} model and each single-feature-only model, and corr(F2,F4) on the real
  Set-1 table, to characterize the collinearity.

## 3. Pre-declared constants and seeds (committed before results)
- Regime experience (Trial 6-25, r = 20); targets loss_avoid_adj (primary), ev_sens_agn (secondary); K in {6,8,10,12}
  for B1, {8,10,12} for B2/B3.
- B1: R_SEED = 8; M = min(2000, C(|F|,K)); N_SPLIT_B1 = 100; noise-floor subsample 50 sets, N_boot = 1000;
  B_PART_B1 = 200. Base SEED_B1 = 20260722. Set-draw seed = SEED_B1 + 1000*form + K + 10000*seed_index; split-draw seed
  = that + 500000; partition seed = SEED_B1 + 700000 + partition_index.
- B2: B_PART_B2 = 200; N_SPLIT_B2 = 50 (bounds runtime; MC error reported); decile = top 10% by selection-half r_sb.
  Base SEED_B2 = 20260723; partition seed = SEED_B2 + partition_index; per-half halving seed derived as in E3b.
- B3: reuse the committed E4 CV seeds (`e4/e4_seeds.csv`) and E4 constants (n_cv_repeats 200, train_fraction 0.7);
  no new stochastic choice. All seeds committed to `dispb_config.json` + `dispb_seeds.csv`; the grader reproduces set
  membership deterministically and matches stochastic quantities at expectation.

## 4. Pre-declared gates (evaluated build-only; committed before results; LOCKED interpretation)
B1:
- G-B1a (primary-cell stability): the MINIMUM sigma_set_true over the R_SEED x 2-estimand battery at (Set1, K12,
  loss_avoid_adj) is >= 0.05 AND the minimum ICC_set >= 0.50 (the frozen E3 thresholds; nothing re-tuned).
- G-B1b (partition stability): the 90% interval lower bound of sigma_set_true_cov over B_PART_B1 partitions >= 0.05.
- G-B1c (null control): permuted sigma_set_true (mean over permutations) < 0.05 and far below the real value.
- G-B1d (map consistency, REPORTED not gated): count of cells clearing both thresholds under the covariance estimand.
B2:
- G-B2a (optimism exists): out-of-sample gain < in-sample apparent gain (winner's curse present).
- G-B2b (residual selectability, SCIENTIFIC, interval-vs-zero RULE): the top-decile out-of-sample gain-vs-mean 90%
  interval over partitions. Lower bound > 0 -> report a modest, transferable composition advantage from direct
  empirical selection at build scale (state the magnitude and that it is roughly half the naive gap). Interval covers
  0 -> report that selection does not reliably beat the average out-of-sample at these N and the best-minus-median gap
  is winner's curse. EITHER outcome is honest; the rule is fixed here, not the magnitude.
- G-B2c (null control): permuted out-of-sample gain within Monte-Carlo noise of 0.
B3:
- G-B3a (null robustness, SCIENTIFIC): NO leave-one-out variant that RETAINS F4 (expected information) restores the
  full frozen H2 PASS (G-sim + G-beyond + G-real sign-consistency across K in {8,10,12}) on the real track. A passing
  F4-retaining variant would be reported EXPLORATORY and could NOT overwrite the frozen H2 FAIL (INTEGRITY.md).
- G-B3b (information redundancy): |Delta OOF(drop F4)| <= 0.05 at K in {8,10,12} on the real track; report the joint
  drop-{F2,F4} delta and single-feature-only fits so mutual redundancy is explicit.
- G-B3c (F1 stability): F1 (target_predictor_variance) keeps a NEGATIVE standardized coefficient in FULL and every
  F1-retaining variant on the real track.
- G-B3d (sim sanity control): sim_free105 FULL pooled OOF >= 0.30 and each sim LOFO variant >= 0.20.

Honesty guards (LOCKED; a violation blocks the positive reading and is brought to the author, never worked around):
- HG-ceiling: no recomputed reliability exceeds the E0 true-variance ceiling beyond bootstrap noise.
- HG-frozen: the E4 H2 gate JSON and the E3 h1_gate JSON are neither re-run nor amended; B3 consumes E4 outputs and
  recomputes no reliability; B1/B2 recompute reliabilities only from BUILD raw with the frozen engine.
- HG-null: the two B1/B2 permutation controls must collapse to ~0, else the corresponding statistic is voided.

Interpretation (LOCKED, so no outcome is a moved goalpost): whatever B1/B2/B3 return, the frozen Study 1 (historical
freeze tag), the frozen H1 PASS, and the frozen H2 FAIL stand untouched. B1 defends H1's robustness; B2
calibrates how much of the composition advantage is exploitable out-of-sample; B3 localizes and pressure-tests the H2
null. A pre-declared-gate FAIL is a valid result and is reported as-is.

## 5. Outputs (all new, under elevation/; additions-only)
- `elevation/analysis/dispb_stability.R`, `dispb_optimism.R`, `dispb_lofo.R`, `dispb_verify.R` - sourcing the frozen
  engines; reading only build data + the frozen split + committed E2/E3/E3b/E4 outputs.
- `elevation/results/dispb/b1_stability_by_seed.csv` - per (K, seed, estimand): sigma_set_raw, sigma_noise,
  sigma_set_true, ICC_set, coverage.
- `elevation/results/dispb/b1_partition_distribution.csv` - per partition: sigma_set_true_cov, ICC_cov (+ the 90%
  interval summary) and the null-permutation values.
- `elevation/results/dispb/b1_map_cov.csv` - the full feasible-cell map under the covariance estimand.
- `elevation/results/dispb/b2_selection_gain.csv` - per (K, partition, direction): in/out-of-sample single-best and
  top-decile gains, cross-half rho, and the null-permutation values; plus a per-K summary with the 90% interval.
- `elevation/results/dispb/b3_lofo_cv.csv` - per (track, target, K, variant): oof_spearman pooled + repeatmean + gate,
  drop-one delta, and the joint/single-feature fits.
- `elevation/results/dispb/b3_lofo_signs.csv` - per (track, target, K, variant): retained standardized coefficients +
  signs + the per-variant H2 three-clause verdict; corr(F2,F4).
- `elevation/results/dispb/dispb_gate.json` - G-B1a..d / G-B2a..c / G-B3a..d verdicts vs thresholds; PASS/FAIL plain.
- `elevation/results/dispb/dispb_config.json` + `elevation/results/dispb/dispb_seeds.csv`.
- `elevation/logs/dispb_*.log`; the run report.

## 6. Exports for independent re-derivation (generator != grader)
The run EXPORTS enough committed per-unit values that independent verification re-derives each headline without rerunning R: for a
committed sample of B1 seeds/partitions, the per-set r_sb and the two subject-half r_sb (so independent verification recomputes the
decomposition and the covariance estimand); for B2, the per-set r_sb^A / r_sb^B for a committed sample of partitions
(so independent verification recomputes the selection gain); for B3, the per-(track,K,variant) OOF metrics and the full standardized
coefficient table. All committed CSVs at >= 15 significant digits; note the readr float-parse caveat (use
utils::read.csv for bit checks) that E3b documented. Set membership is reproducible from dispb_seeds.csv.

## 7. Pre-declared PASS gates (Part D)
- D1 engine provenance: `R/step2_reliability.R` + `R/cpc18_features.R` SOURCED, ZERO redefinitions (grep the frozen
  function names); B1/B2 reliabilities come from the frozen estimate_half / split_correlations / sb_correct /
  reliability_point; a G_EQ check re-proves the B1/B2 accelerator (if any) reproduces the frozen estimate_half
  bit-for-bit on a committed sample of halves before any result. B3 recomputes NO reliability and reads NO raw.
- D2 build-only / no leakage: split == "build" filtered before any choice read; choice columns dropped from the
  non-build frame; rm(raw); split counts asserted 344/137/205; freeze not read; the 205 test SEALED.
- D3 gates computed and reported: G-B1a..d, G-B2a..c, G-B3a..d vs the pre-declared thresholds; PASS/FAIL stated
  plainly; the full per-cell / per-partition / per-variant tables delivered; any failure reported honestly.
- D4 outputs present: every file in section 5, with the per-unit samples for re-derivation.
- D5 reproducibility: all seeds committed; set membership reproducible; stochastic quantities reproducible in
  distribution; the committed per-unit samples re-derive the reported headline values.
- D6 additions-only + frozen integrity (CRYPTOGRAPHIC): git diff <pre-043 HEAD e2b1a2f8> -> HEAD touches ONLY
  {the session log, the coordination directory, elevation}; the 13 frozen Study-1 subtrees object-hash IDENTICAL to freeze 387f9adc;
  elevation/results e0 446b1420 / e1 0e020fb0 / e2 89c7d2f0 / e3 36120997 / e3b 84703371 / e4 307fc9c6 / e8 e0fc1688 /
  p2_irt 0a0aec0b / p3_frey 2b305379 / sens all IDENTICAL to the pre-043 tip (nothing prior touched; dispb the only
  results addition); E1_contract / E2_spec / E3_spec / E3b_spec / E4_spec / E8_spec / SENS_spec blobs unchanged;
  master 0808abc + both tags untouched; branch pushed.
- D7 no forbidden computation: NO optimization, NO candidate-design construction, NO new choice-model fit, NO
  freeze/test read; the frozen E3 h1_gate and E4 h2_gate neither re-run nor amended; no reserved E5-E7 slot reused;
  simulation labeled SIMULATION.
- D8 hygiene: the pre-declared status labeled honestly (Battery B is a newly pre-declared
  robustness re-analysis, not part of the frozen E1 confirmatory arc).

## STOP
After D1-D8: STOP. Report B1 (sigma_set_true stability + null), B2 (single-best optimism + top-decile out-of-sample
gain + interval + null), and B3 (LOFO deltas + the F4 sign-locus + redundancy) with plain gate verdicts. Do NOT
construct designs, open freeze, or read any test outcome; the 205 test subjects remain sealed. Whatever the results,
the frozen Study 1, the frozen H1 PASS, and the frozen H2 FAIL stand untouched. independent verification reviews against the
independently prepared answer key, independently re-deriving the headline
B-statistics from raw before grading.

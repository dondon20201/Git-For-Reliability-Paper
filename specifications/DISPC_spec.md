# Battery C - Study 2 defensive diagnostics: is the equal-length heterogeneity the construct or a confound? (SPEC; committed BEFORE any result)

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md. Companions: E3_spec.md (H1), DISPB_spec.md
(stability/optimism/LOFO), E8_spec.md (robustness-analysis template).

C is the final optional robustness study in Study-2. A = SENS, B = Battery B, C = Battery C. Battery C tests for the influence of response style contamination, local item dependence and scoring rule effects on heterogeneity in equal-lengths between E3 H1 using three pre-specified analyses. These analyses are interpretable in both directions. If a confound is present then it will realistically narrow H1. If a confound is not present then it will support the interpretation of the construct in the face of an attempted criticism by a reviewer. Battery C does not create designs, open a freeze or read any test taker (205 sealed). It does not overwrite the pre-specification that was frozen. The E3 H1 PASS and the E4 H2 FAIL remain unchanged. This spec is made before the results it influences and the gates are not changed after results are available.

## 0. What Battery C is (and is NOT)
- Three PRE-DECLARED diagnostics on the E3 equal-length sets at the primary cell (Set 1, K=12, loss_avoid_adj), each
  with a fixed gate and a stopping rule. NOT a feature search: no additional response-style, dependence, or estimator
  variant is added after seeing results (that is the feature-fishing the development-options not-helpful list warns
  against).
- BUILD-ONLY (Set-1 build subjects; the frozen 63). FREEZE (137) not read; TEST (205) sealed. A re-analysis, not a
  selection step.
- Reuses the E3 sets (reproduced from E3's committed seeds) and the FROZEN scoring engine. No new item features, no new
  mechanism model, no optimization, no candidate designs.
- All ICC comparisons use the COVARIANCE (subject-split) estimand as PRIMARY, per Battery B (the bootstrap-floor ICC
  is seed-unstable: B found ICC 0.37-0.80 across 8 seeds; the covariance ICC is the stable one). The bootstrap-floor ICC
  is reported for continuity only.
- Ceilings bind (E0): within-form reliabilities are shared-subject / shared-instrument inflated (the B HG-ceiling
  finding); Battery C compares RELATIVE quantities (ICC before/after, rankings), which that inflation does not affect,
  but the flag is carried.

## 1. Inputs (committed; build-only; additions-only under elevation/)
- `data/raw/All estimation raw data.csv` (assert md5 `38a26407a36c690be969f9d9d3fe3e5a`), restricted to
  `split == "build"` and Set 1 via `elevation/results/e1/participant_split.csv` (assert 344/137/205; build filtered
  BEFORE any choice read; non-build choice columns dropped; rm(raw)). Needed for the response-style scores, the robust
  score, and the within-set residual correlations, which are not in any committed file.
- `elevation/results/e2/problem_information_dictionary.csv` - per-problem features (in_Set1, loss_informative,
  has_loss_b, d_ev_noloss, ev_a, ev_b, sd_a, sd_b) for the 24 Set-1 clean problems.
- `elevation/results/e3/e3_seeds.csv` - to reproduce the exact E3 Set-1 set membership.
- Frozen engines SOURCED, not redefined: `R/step2_reliability.R` (build_subject_data, fit_loss_baselines,
  estimate_half, split_correlations, sb_correct, reliability_point), `R/cpc18_features.R`. The exp loss baselines are
  the committed E3 constants (intercept 0.4301295476692, slope 0.165107793793507).

## 2. The three diagnostics (primary cell Set 1, K=12, loss_avoid_adj; exp Trial 6-25, r=20; ev_sens_agn informational)

### C1 - response-style contamination
Per subject compute two stable response-style traits over ALL 24 Set-1 clean problems (experience trials):
- B-side bias = mean B-choice rate (tendency to choose option B irrespective of EV or loss content);
- safe-choice rate = mean rate of choosing the lower-SD option (0.5 on an SD tie; tendency to choose the safer option
  irrespective of EV).
STRUCTURAL NOTE (pre-declared handling): in Set 1 option B is the higher-SD option for every clean problem, so
safe-choice rate == 1 - B-side bias exactly (corr = -1). Enter both traits; where |corr(B-side, safe)| > 0.95 (as in
Set 1) residualize on the COMBINED single side/risk axis, reported, to avoid a degenerate two-regressor fit. For each
K=12 set and each subject-half, residualize the per-subject loss_avoid_adj on the response-style trait(s) across
subjects, and recompute the between-set covariance sigma_true_cov and ICC_cov on the residualized score.

### C2 - local item dependence
Per set, wpc = within-set mean pairwise Spearman of the per-problem loss residuals (chose-no-loss rate minus the
committed EV baseline) across subjects, over the set's loss-informative problems. Report: Spearman(set r_sb, wpc) across
sets; the low-vs-high-wpc tercile reliability comparison; and a dependence-adjusted ICC (the set reliability compared to
the Spearman-Brown prediction from the mean single-problem reliability, where computable). The discriminating question
is whether wpc exceeds what construct coherence alone predicts; report honestly that r_sb tracking wpc is expected under
BOTH construct coherence and nuisance dependence, so a strong raw relationship is not by itself evidence of a dependence
artifact.

### C3 - alternative-estimator robustness of the set ranking
Recompute set-level reliability under a pre-declared robust variant of the loss score - the MEDIAN (not mean) of the
per-problem loss residuals over a half's loss-informative problems (robust to a single outlier problem) - and correlate
the set ranking (by r_sb) with the primary mean-based ranking across the M sets. A rank-transformed loss variant may be
reported alongside; the median already establishes robustness.

## 3. Pre-declared constants and seeds (committed before results)
Regime experience (Trial 6-25, r=20); target loss_avoid_adj (primary), ev_sens_agn (informational); K=12 primary (K in
{6,8,10} reported). M = min(2000, C(24,K)); N_SPLIT = 100; covariance estimand over B_PART >= 200 subject partitions
(the robust ICC per Battery B); bootstrap-floor ICC reported for continuity. Response-style traits over all 24 Set-1
problems; robust score = median per-problem residual; wpc = mean pairwise Spearman over loss-informative problems. Base
SEED_C = 20260724; set-draw reproduces E3's committed Set-1 seeds; partition and split seeds derived and committed to
`dispc_config.json` + `dispc_seeds.csv`. The grader reproduces set membership deterministically and matches stochastic
quantities at expectation.

## 4. Pre-declared gates (evaluated build-only; committed before results; LOCKED interpretation)
- G-C1 (response style): relative cov-ICC drop = (ICC_raw - ICC_resid) / ICC_raw. **>= 50%** -> response style is a
  material driver; narrow H1 to "composition affects reliability partly through stable response tendencies." **< 15%**
  -> not a response-style artifact; strengthens the construct reading. **15-50%** -> partial; report the qualified
  reading. Report the R^2 of loss on the response-style axis as context. The Set-1 |corr|>0.95 collinearity is handled
  as in section 2 (combined-axis residualization) and reported as a structural caveat.
- G-C2 (local dependence): a dependence-adjusted ICC reduction past the pre-declared threshold (default 50% relative,
  matching C1) OR a low-vs-high-wpc reliability gap that exceeds what construct coherence predicts -> the heterogeneity
  is partly a dependence artifact; narrow. Otherwise -> robust to local dependence. If the outcome is too noisy to
  discriminate (the E3b outcome-noise limit), report INCONCLUSIVE - a valid pre-declared outcome; do NOT expand into a
  search.
- G-C3 (alt-estimator): set-ranking Spearman(primary, robust) **>= 0.30** -> the heterogeneity is not scoring-specific.
  Below -> the "which sets" result is estimator-specific and must be stated as such.

Honesty guards (LOCKED): HG-frozen (the E3 h1_gate and E4 h2_gate are read-only, never re-run or amended; no new
mechanism model; no optimization); HG-ceiling (within-form reliabilities are shared-instrument inflated per B; C uses
relative comparisons, but the flag is carried); HG-build (build-only; freeze not read; 205 sealed). Stopping rule: the
three diagnostics, their constructs, and their thresholds are fixed before results; no extra variant is added post hoc;
a weak or inconclusive result NARROWS the claim rather than spawning more models. Interpretation LOCKED so no outcome is
a moved goalpost: whatever C returns, the frozen Study 1, the H1 PASS, and the H2 FAIL stand untouched; each diagnostic
sharpens H1 whichever way it lands.

## 5. Outputs (all new, under elevation/; additions-only)
- `elevation/analysis/dispc_diagnostics.R`, `dispc_verify.R` - sourcing the frozen engines; reading only build data +
  the frozen split + committed E2/E3 seeds.
- `elevation/results/dispc/response_style_adjusted_icc.csv` - per (K, estimand): sigma_true + ICC before/after
  residualizing each response-style axis; corr(B-side, safe); R^2(loss ~ styles).
- `elevation/results/dispc/local_dependence_diagnostic.csv` - per set: wpc, r_sb; the Spearman(r_sb, wpc), the tercile
  comparison, and the dependence-adjusted ICC.
- `elevation/results/dispc/alt_estimator_set_ranking.csv` - per set: primary r_sb, robust r_sb; the ranking Spearman.
- `elevation/results/dispc/dispc_gate.json` - G-C1 / G-C2 / G-C3 verdicts vs thresholds; PASS / narrow / strengthen /
  inconclusive stated plainly; the Set-1 collinearity flag.
- `elevation/results/dispc/dispc_config.json` + `dispc_seeds.csv`; the D1 G_EQ + spotcheck files.
- `elevation/logs/dispc_*.log`; the run report.

## 6. Exports for independent re-derivation (generator != grader)
The run EXPORTS enough committed per-unit values that independent verification re-derives each headline without rerunning R: for a committed
sample of sets/partitions, the per-subject loss score, the two response-style traits, the residualized score, the
per-problem residuals, and the per-set r_sb under both scorings (so independent verification recomputes the residualized cov-ICC, the wpc
relationship, and the ranking Spearman). All committed CSVs at >= 15 significant digits; use utils::read.csv for bit
checks (the readr float-parse caveat E3b documented). Set membership reproducible from dispc_seeds.csv.

## 7. Pre-declared PASS gates (Part D)
- D1 engine provenance: `R/step2_reliability.R` + `R/cpc18_features.R` SOURCED, ZERO redefinitions (grep the frozen
  function names); the loss_avoid_adj score comes from the frozen estimate_half; a G_EQ check re-proves the from-raw
  loss score reproduces the frozen engine bit-for-bit on a committed sample of halves before any result.
- D2 build-only / no leakage: split == "build" filtered before any choice read; non-build choice columns dropped;
  rm(raw); split counts asserted 344/137/205; freeze not read; the 205 test SEALED.
- D3 gates computed and reported: G-C1 / G-C2 / G-C3 vs the pre-declared thresholds; verdicts (incl. any inconclusive)
  stated plainly; the Set-1 response-style collinearity reported.
- D4 outputs present: every file in section 5, with the per-unit re-derivation sample.
- D5 reproducibility: all seeds committed; set membership reproducible from E3's seeds; stochastic quantities
  reproducible in distribution; the committed per-unit sample re-derives the reported headline values.
- D6 additions-only + frozen integrity (CRYPTOGRAPHIC): git diff <pre-044 HEAD 4585b9ac> -> HEAD touches ONLY
  {the session log, the coordination directory, elevation}; the 13 frozen Study-1 subtrees object-hash IDENTICAL to freeze 387f9adc; all prior
  elevation result subtrees (e0/e1/e2/e3/e3b/e4/e8/p2_irt/p3_frey/sens/dispb) IDENTICAL to the pre-044 tip (nothing
  prior touched; dispc the only results addition); the frozen spec blobs (E1_contract / E2_spec / E3_spec / E3b_spec /
  E4_spec / E8_spec / SENS_spec / DISPB_spec) unchanged; master 0808abc + both tags untouched; branch pushed.
- D7 no forbidden computation: NO optimization, NO candidate-design construction, NO new mechanism-model fit, NO
  freeze/test read; the frozen E3 h1_gate and E4 h2_gate neither re-run nor amended; no reserved E5-E7 slot reused.
- D8 hygiene: the pre-declared status labeled honestly (Battery C is a newly pre-declared
  robustness re-analysis, not part of the frozen E1 confirmatory arc).

## STOP
After D1-D8: STOP. Report C1 (response-style-adjusted cov-ICC + the Set-1 collinearity), C2 (dependence relationship +
whether it discriminates), and C3 (the alt-estimator ranking Spearman) with plain gate verdicts. Do NOT construct
designs, open freeze, or read any test outcome; the 205 test subjects remain sealed. Whatever the results, the frozen
Study 1, the H1 PASS, and the H2 FAIL stand untouched. Verification proceeds against an independently prepared answer key, re-deriving the headline diagnostics from raw before grading.

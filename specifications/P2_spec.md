# P2 information-IRT mechanism model - specification (joint item-information / IRT). Does joint estimation recover the H2 mechanism?

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md. This is the P2 leg of
the third-act package (risk ordering P1/E8 done -> P3 done -> P2/IRT now). "P2" is the IRT path
label from the third-act decision record; it runs THIRD by risk order, after P3.

P2 seeks to explore a JOINT item-response / information modeling strategy in which individual differences on the latent risk trait are estimated along with item discrimination (the ability of an item to discriminate between individuals on the latent risk trait) and then item discrimination is related to the previously frozen E2 values of information. The aim is to see if such a strategy can capture the mechanism proposed in Study-2 that was not captured by the previously frozen two-stage H2 test (E4). Unlike E4 that first estimated per-set reliability and then regressed on features, E3b found that this reliability measure was noisy (out-of-subject reliability = 0.28-0.48 and partition-unstable = 0.25-0.83) and therefore could not be used to either affirm or deny the mechanism. Additionally E4's real anchor exceeded its own outcome ceiling (leakage). P2 avoids the noisy per-set reliability measure by fitting the mechanism at the item level in a single JOINT model. This model is BUILD-only and does not read test or freeze subjects. Additionally the previously frozen Study 1, the previously frozen H2 FAIL or E8 remains unchanged.

## 0. What P2 is (and is NOT)
- P2 fits a one-factor item-response model to the BUILD loss-avoidance choice data: each problem j gets a
  discrimination/loading lambda_j on the latent loss-avoidance trait u_i, estimated JOINTLY with item
  easiness and the frozen EV control; then regresses the estimated lambda_j on the frozen E2 per-problem
  information features (decorrelated). The mechanism claim under test is: features -> lambda_j -> set
  reliability (one-factor test information = sum of lambda_j^2 over a set).
- PRIMARY target = loss_avoid_adj (E0 ceiling 0.5446, the trustworthy-variance index). SECONDARY =
  ev_sens_agn (E0 ceiling 0.2811), the same loading structure on the EV-sensitivity trait, reported as
  robustness only, never the headline.
- P2 is EXPLORATORY and NEW. The frozen pre-specified H2 FAIL (E4) STANDS and is NOT replaced
  (INTEGRITY.md: a re-analysis is new + exploratory and cannot overwrite a frozen result). A positive P2 is an
  exploratory mechanism finding that COMPLEMENTS the frozen null ("the mechanism is visible only under joint
  estimation"); a null or non-identified P2 further supports the honest measurement-limits framing.
- P2 anchors ONLY to the two trustworthy constructs. It does NOT build a latent trait on the Firth
  coefficients b_ev / b_risk / b_loss: E0 LOCKED that their metric scale quasi-separates (within-subject SE
  up to 1e15-1e16; 0/3 simulation recovery; see elevation/results/e0/coefficient_separation_diagnostic.csv),
  so an item-loading model on them is unidentified by construction and is out of scope.
- P2 is BUILD-ONLY (344 subjects). The 137 freeze subjects are not read; the 205 test subjects stay sealed
  (the E7 rule: never read before a deliberate confirmatory step). P2 is a re-analysis, not a selection step,
  so freeze is not opened. It constructs NO candidate designs (the retired E5), opens NO freeze, and reads NO
  problem outside the frozen clean-subset definition.
- Ceilings bind (E0, LOCKED). No estimator exceeds the population true-variance ceiling; E8 already showed
  partial pooling does not raise loss_avoid_adj reliability (the wall is a true-variance ceiling, not
  estimation noise). P2 does NOT re-litigate the wall. It asks the ORTHOGONAL question of whether item
  discrimination is STRUCTURED by information features. An IRT-implied reliability above the E0 ceiling is a
  leakage flag, not a result.

## 1. Inputs (committed; build-only; additions-only; frozen artifacts SOURCED read-only)
- Raw CPC18 clean-subset trials (md5 38a26407a36c690be969f9d9d3fe3e5a asserted in-script), restricted to
  split == "build" via elevation/results/e1/participant_split.csv (assert 344 / 137 / 205; build filtered
  BEFORE any choice is read; choice columns dropped from the non-build frame; rm(raw)). Experience regime,
  frozen (Trial 6-25, r = 20), matching Study 1 / E8.
- Frozen engines SOURCED, never redefined: R/step1_measures.R, R/step2_reliability.R (build_subject_data,
  estimate_half, split_correlations, sb_correct, reliability_point), R/cpc18_features.R. The clean-subset,
  the loss-informative trial definition, and the EV controls are the frozen Study-1 pipeline; P2 adds ONLY
  the new joint estimator on top.
- Frozen E2 problem-information dictionary, read-only:
  elevation/results/e2/problem_information_dictionary.csv (per-problem features F1 target-predictor spread,
  F2 loss exposure, F3 predictor collinearity, F4 expected/Fisher information, F5 choice entropy) plus
  elevation/results/e2/e2_feature_config.json. Features are SOURCED, never recomputed; assert the file md5 in
  the config manifest.
- Frozen E3b denoised outcome, read-only (for the G3 calibration gate and the HG3 outcome ceiling):
  elevation/results/e3b/outcome_reliability_real.csv (cross-half Spearman of per-set reliability at
  K = 6/8/10/12) and elevation/results/e3b/reliability_subject_split_real.csv (the per-set subject-split
  reliabilities). The clean simulation testbed is elevation/results/e3b/reliability_highprec_sim.csv +
  e3b_set_identity.csv + e3b_seeds.csv (so the P2 sim track reuses E3b's committed sets and seeds; no new
  set universe is invented).
- E0 ceiling anchors, read-only, from elevation/results/e0/variance_ceiling_final.csv, used ONLY as the HG1b
  ceiling reference (never refit): loss_avoid_adj ceiling_fulldata 0.5446 (r1 = 0.0787); ev_sens_agn 0.2811
  (r1 = 0.02717).

## 2. The joint model (LOCKED) - the only new estimation
For subject i and problem j, aggregate the frozen loss-informative experience trials to a
subject-by-problem binomial count of chose_noloss out of the r sampled trials (the frozen regime; identical
inputs to E8's HIER arm). The PRIMARY one-factor item model is:

  logit( E[chose_noloss_ij] ) = beta0_j + lambda_j * u_i + gamma * d_ev_noloss_ij,
  u_i ~ Normal(0, 1)   (identification: latent trait variance fixed to 1; sign fixed so the mean loading is
                        positive; equivalently a single reference item loading fixed positive),
  lambda_j            = item discrimination / loading (the quantity the mechanism regression explains),
  beta0_j             = item easiness (free per problem),
  gamma               = shared EV-baseline slope on d_ev_noloss (the frozen EV control, matching E8).

This is a confirmatory one-factor item-factor (2PL-with-covariate) model: u_i is the latent loss-avoidance
trait, lambda_j is how strongly problem j discriminates subjects on it. It EXTENDS E8's random-intercept
model by letting each problem carry its own loading instead of a common one, which is exactly the item-level
structure the H2 mechanism is about.

- SECONDARY (ev_sens_agn), analogous on clean trials: the EV-sensitivity trait v_i ~ Normal(0,1) enters
  through an item-specific amplification of the EV contrast,
  logit(E[choice_B_ij]) = beta0_j + kappa_j * z_ev_ij * v_i + delta * z_ev_ij, with kappa_j the item loading.
  Reported as robustness only (LOCKED caveat, as in E8: this is a monotone re-scaling of the frozen Spearman
  ev index, not the identical statistic).

- Estimator family (LOCKED):
  - PRIMARY engine: marginal-maximum-likelihood item-factor analysis via `mirt` (its EM + the fixed-quadrature
    or QMCEM estimator; `mirt` ships compiled and needs NO external C++ toolchain on the analysis machine). Missing-by-design
    across the 7 fixed forms is handled natively (each subject contributes only the items on their form). The
    r experience trials per cell are represented by the frozen binomial count via response-pattern weighting
    (or an equivalent expansion), asserted equivalent to the binomial likelihood on a committed spot sample.
    Standard errors from the observed information matrix.
  - CONFIRMATION engine (PRIMARY target, reference budget only): a Bayesian one-factor item model
    (brms/Stan). If the C++ toolchain gap that blocked brms in E8 recurs, substitute a second
    frequentist item-factor engine (`TAM` or a `glmmTMB` loading approximation) and say so; report
    loading agreement Spearman(lambda_primary, lambda_confirm) and feature-beta sign/magnitude agreement.
  - Engine risk is folded into G0: if NO available engine fits the loading structure with identified,
    finite standard errors, that IS a G0 identifiability FAIL. STOP and report it as the pre-declared result.

## 3. The mechanism regression and the calibration bridge (committed before results)
- Feature link (the mechanism): regress the estimated PRIMARY loadings lambda_j on the frozen E2 features,
  DECORRELATED. Decorrelation is LOCKED because E4 found corr(F2 loss exposure, F4 expected information) =
  0.90 (real) / 0.96 (sim): orthogonalize in the pre-declared order [F1, F2, F3, F5, then F4 residualized on
  {F1,F2,F3,F5}], and ALSO report a single combined information factor (first PC of {F2,F4}) as the
  transparent alternative. REPORT BOTH the multiple-regression partial coefficients AND the marginal
  (univariate) rank correlations of each raw feature with lambda_j (E4's lesson: F4's marginal was +0.668
  while its partial flipped sign across folds). Evaluate out-of-fold with GROUPED cross-validation, folds =
  whole forms / assignment blocks held out (execution-spec 10.2: subsets share items, so independent SEs are
  invalid; use grouped CV + report predictive performance, not p-values). Baseline = problem count + broad
  category; the mechanism must beat it.
- Calibration bridge (loadings -> set reliability): for a K-problem set S, the one-factor IRT-implied
  reliability is rho_hat(S) = I(S) / (I(S) + 1) with test information I(S) = sum_{j in S} lambda_j^2 at the
  population mean (report the exact mapping used). Compare rho_hat(S) by RANK against the E3b denoised
  subject-split reliability of the SAME sets (reliability_subject_split_real), across K in {8,10,12}. This is
  the bridge from "features structure loadings" to "features predict set reliability", using E3b's committed
  sets so no new outcome is invented.
- Simulation testbed (the clean track, primary per E3b): reuse E3b's committed high-precision simulation
  (reliability_highprec_sim + e3b_set_identity + e3b_seeds; empirical parameter distribution), where the
  across-set spread is mostly true composition signal, to test the feature -> loading map free of the real
  outcome's noise. The real BUILD track is the low-ceiling honesty check, bounded by E3b (below).

## 4. Pre-declared gates (evaluated build-only; committed before results)

G0 IDENTIFIABILITY (HARD; the EARLY STOP gate; declared blind). PASS only if ALL of G0a-G0d hold:
- G0a convergence: the primary `mirt` fit reaches its EM convergence criterion and returns a finite observed
  information matrix with finite loading SEs; if the Bayesian confirmation engine runs, max Rhat < 1.01, 0
  divergent transitions, and min bulk-ESS > 400 on the reported lambda_j and feature betas.
- G0b bounded, non-degenerate loadings: every lambda_j is finite with a finite SE; the fraction of items
  flagged DEGENERATE (near-zero between-subject response variance, or |lambda_j| SE > 5x the median loading
  SE, or |lambda_j| > 6 on the standardized scale) is BELOW 10% of items. At or above 10% degenerate,
  quasi-separation dominates -> G0 FAIL.
- G0c scale + conditioning: the u ~ N(0,1) identification is imposed; the decorrelated feature-regression
  design matrix has condition number < 30 (residual collinearity beyond the F2/F4 decorrelation is
  controlled); the feature betas have finite SEs.
- G0d simulation recovery AND null (LOCKED; the E0 critical check in INTEGRITY.md):
  - Recovery: simulate BUILD-shaped data with a PLANTED feature -> loading map (a known positive beta on a
    designated information feature); refit the full pipeline; the recovered beta's 90% interval covers the
    planted value AND corr(true lambda, estimated lambda) > 0.60.
  - Null-A (no mechanism): simulate loadings drawn INDEPENDENTLY of features; the estimated feature -> loading
    beta 90% interval MUST include 0.
  - Null-B (no trait): simulate var(u) = 0; the loadings must be near-zero / non-identified and the feature
    beta interval MUST include 0. The model must not manufacture a mechanism from noise (E0 tau = 0 logic).
- G0 verdict: if any of G0a-G0d fails, STOP after Phase 1a. Report the identifiability failure as the
  pre-declared RESULT (a clean methodological null consistent with E0's quasi-separation and the fixed-form
  no-anchor design). Do NOT run G1-G3 on real data, do NOT tune the model to force identification.

Scientific gates (evaluated ONLY if G0 passes):
- G1 (MECHANISM, SIMULATION primary; the clean testbed): on the E3b simulation testbed, the decorrelated
  features predict held-out lambda_j out-of-fold (grouped CV) with Spearman >= 0.30, beating the
  count+category baseline (delta > 0), sign-consistent across the primary target and across the K used in
  calibration. This is the pre-specified pass condition for "information structures item discrimination".
- G2 (MECHANISM, REAL honesty check; low-ceiling; NOT a headline): the same feature -> loading relationship
  on real BUILD, reported bounded by the E3b outcome-reliability ceiling. Reported for honesty and
  sign-agreement with G1; a real null does NOT fail P2 by itself (E3b: the real outcome is too noisy at these
  N to confirm OR refute), and a real "success" that exceeds its own outcome ceiling is leakage (HG3), not a
  win.
- G3 (RELIABILITY CALIBRATION): the IRT-implied set reliability rho_hat(S) rank-correlates with the E3b
  denoised subject-split reliability of the same sets, Spearman >= 0.30 at K in {8,10,12} (K = 6 reported but
  flagged thin, per E3b). This is the bridge gate; without it, loadings need not translate into set
  reliability.

Honesty guards (LOCKED; hard; a violation blocks the positive reading and is brought to the author, never worked
around):
- HG1a cross-form anchor guard: the 7 fixed forms share NO common items, so the latent scale is linked only
  through the population-exchangeability assumption. Report a covariate-exchangeability check across forms
  (frozen experience / demographic distributions) AND a within-form sensitivity fit on the largest single
  BUILD form (no cross-form linking needed there); the within-form feature -> loading SIGN must agree with the
  pooled result. Disagreement means the pooled mechanism is an artifact of the linking assumption -> flagged,
  not banked.
- HG1b ceiling guard: rho_hat(S) must not exceed rel_ceiling(K) from the E0 anchors beyond bootstrap noise;
  a CI-supported exceedance is a leakage flag (echoing E3b, where E4's real anchor exceeded its honest
  ceiling), reported as such, never as a gain.
- HG2 null-permutation: permuting the feature rows across items must drop BOTH the feature -> loading beta and
  the CV rank correlation to ~0 (within Monte-Carlo noise). A materially non-zero permuted association voids
  the mechanism (the E0 tau = 0 recovery logic applied to the item side).
- HG3 outcome-ceiling / leakage: the real-track out-of-fold predictive correlation for set reliability must
  NOT exceed sqrt(E3b outcome reliability) at the matched budget; an exceedance is flagged leakage and not
  banked (this is exactly what E3b caught in E4's real anchor).
- HG4 non-degeneracy: report that the latent-trait variance did not collapse and the loading distribution is
  non-degenerate (not all lambda_j equal); mirrors E8's HG4.

Interpretation (LOCKED, so no outcome is a moved goalpost; the frozen H2 FAIL stands throughout):
- G0 FAIL: the joint mechanism is not identifiable on CPC18 BUILD at these budgets (quasi-separation / the
  fixed-form no-anchor design / item degeneracy). A clean pre-declared methodological null that REINFORCES
  the measurement-limits contribution: even a principled joint model cannot extract the mechanism here.
- G0 PASS, G1(sim) + G3(calib) PASS, G2(real) null: the mechanism is identifiable and recovered in
  simulation (information features structure item discrimination) but the real BUILD outcome is too noisy to
  confirm it out-of-sample (E3b ceiling). Reported as a positive STRUCTURAL / simulation result, honestly
  bounded, an EXPLORATORY complement to the frozen H2 null - never a rescue of it.
- G0 PASS, G1 + G2 + G3 PASS: the strongest outcome. Joint estimation recovers the information mechanism the
  two-stage E4 could not see; reported as an EXPLORATORY positive that reframes the H2 story from "no
  mechanism" to "the mechanism requires joint estimation to detect", still labeled exploratory and NOT
  overwriting the frozen pre-specified H2 FAIL.
- A failed gate is a valid scientific outcome (INTEGRITY.md). Do NOT repurpose the analysis to pass.

## 5. Outputs (all new, under elevation/; additions-only)
- elevation/analysis/p2_irt_mechanism.R - the analysis script (sources the frozen engines; reads only build
  data + the frozen split + the frozen E2 / E3b / E0 result files named in section 1).
- elevation/results/p2_irt/loadings_primary.csv - per problem j: lambda_j, SE, easiness beta0_j, degeneracy
  flag, form/block, and the E2 features joined (loss_avoid_adj).
- elevation/results/p2_irt/loadings_secondary.csv - the ev_sens_agn kappa_j analog.
- elevation/results/p2_irt/mechanism_regression.csv - decorrelated partial coefficients + marginal rank
  correlations + grouped-CV out-of-fold Spearman + count+category baseline delta (real AND sim tracks).
- elevation/results/p2_irt/calibration.csv - G3: rho_hat(S) vs E3b subject-split reliability, per K, with the
  rank correlation and the HG1b ceiling check.
- elevation/results/p2_irt/identifiability.csv - G0a-G0d diagnostics: convergence, degeneracy fraction, SE
  summaries, condition number, and the recovery + null-A + null-B simulation outcomes.
- elevation/results/p2_irt/null_permutation.csv - HG2 permuted feature -> loading beta + CV correlation.
- elevation/results/p2_irt/crossform_sensitivity.csv - HG1a covariate-exchangeability + within-form refit
  sign agreement.
- elevation/results/p2_irt/confirmation_engine.csv - loading + feature-beta agreement between the primary and
  confirmation engines (or the documented engine substitution).
- elevation/results/p2_irt/per_item_export.csv + per_subject_trait_sample.csv - committed per-item loadings
  and a committed sample of per-subject trait estimates, so independent verification re-derives the mechanism regression and the
  calibration from committed values without rerunning the fit.
- elevation/results/p2_irt/p2_gate.json - G0 (a-d) / G1 / G2 / G3 verdicts vs thresholds; HG1-HG4 flags;
  PASS/FAIL stated plainly; the LOCKED interpretation branch selected.
- elevation/results/p2_irt/p2_config.json (with input md5s) + elevation/results/p2_irt/p2_seeds.csv.
- elevation/logs/p2_irt.log; the run report; the session log.

## 6. Exports for independent re-derivation (generator != grader)
The run EXPORTS enough committed per-unit values that independent verification re-derives the headline without rerunning the fit:
the full per-item loadings (per_item_export.csv) so independent verification re-runs the decorrelated feature regression and the
grouped-CV rank correlation from the committed lambda_j; a committed sample of per-subject trait estimates so
the trait scale is checkable; the calibration table so rho_hat vs E3b is re-derivable; the null-permutation
and identifiability tables. All committed CSVs at full round-trippable precision (>= 15 significant digits);
note the readr float-parse caveat (use utils::read.csv for bit checks) that E3b documented.

## 7. Housekeeping (separate, labeled; the one allowed top-level change)
If P2 installs new packages (mirt, and any confirmation engine), run `renv::snapshot()` as a SEPARATE
`chore:` commit touching ONLY renv.lock (the P3 precedent). If snapshotting would touch anything other than
renv.lock, STOP and report. The P2 analysis commit itself stays additions-only under elevation/.

## 8. Pre-declared PASS gates (Part D)
- D1 engine provenance: R/step1_measures.R + R/step2_reliability.R + R/cpc18_features.R SOURCED, ZERO
  redefinitions (grep the frozen function names); the E2 features, E3b outcomes, and E0 ceilings are READ from
  the committed frozen files named in section 1 (asserted by md5), never recomputed; the loss-informative
  trial construction matches the frozen pipeline on a committed spot sample before any loading is fit.
- D2 build-only / no leakage: split == "build" filtered before any choice read; choice columns dropped from
  the non-build frame; rm(raw); split counts asserted 344 / 137 / 205; freeze not read; the 205 test sealed.
- D3 model as specified: the primary and secondary item models are exactly section 2 (one factor,
  per-item loading, the frozen EV control, u ~ N(0,1) identification, no extra covariates, no form fixed
  effect beyond the item easiness); the feature regression is the section-3 decorrelated form with grouped CV.
- D4 outputs present: every file in section 5, including the per-item / per-subject exports.
- D5 reproducibility: all seeds committed (p2_seeds.csv); the fit and the simulation tracks reproduce; the
  committed per-item export re-derives the reported mechanism regression and calibration exactly.
- D6 additions-only + frozen integrity (CRYPTOGRAPHIC): the P2 feat commit's diff touches ONLY elevation/
  (plus the session log and coordination directory in the docs commit); the separate chore commit (if any) touches ONLY renv.lock; the
  13 frozen Study-1 subtrees object-hash IDENTICAL to freeze 387f9adc (renv.lock excepted only if snapshotted);
  elevation/results e0 446b142047e9 / e1 0e020fb0e7ea / e2 89c7d2f0d0a1 / e3 36120997f9f8 / e3b 847033713ea5 /
  e4 307fc9c6b49f / e8 e0fc1688 / p3_frey 2b305379 subtrees IDENTICAL to the pre-041 tip (p2_irt the only
  results addition); E1_contract / E2_spec / E3_spec / E3b_spec / E4_spec / E4b_spec / E8_spec / P3_spec blobs
  unchanged; master 0808abc1 + both historical freeze tags (at commits 92008c74 and 83ebfd82) untouched; pushed.
- D7 no forbidden computation: NO optimization, NO candidate-design construction, NO freeze/test read; the
  frozen H2 gate (e4) neither re-run nor amended; no reserved E5-E7 slot reused; P2 labeled EXPLORATORY, a new
  re-analysis, not part of the frozen E1 confirmatory arc.
- D8 hygiene: the frozen-H2-stands and exploratory labels
  present wherever the mechanism result is stated; the temporal-vs-internal and precision-vs-reliability
  distinctions kept explicit.

## STOP
After Phase 1a (fit + G0): if G0 FAILS, write the identifiability result + REPORT, commit, push, and STOP -
do not run G1-G3. If G0 PASSES, complete G1-G3 + HG1-HG4, write all outputs, commit, push, and STOP. Either
way: construct NO designs, open NO freeze, read NO test outcome; the 205 test subjects remain sealed. Whatever
the verdict, the frozen Study 1 (historical freeze tag), the frozen H2 FAIL, E8, and P3 stand untouched,
and the P2 result is EXPLORATORY and cannot replace the frozen H2 FAIL. independent verification reviews against the
independently prepared answer key, independently re-deriving the loadings and the
mechanism regression before grading.

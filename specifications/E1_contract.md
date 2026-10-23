# Studies 2-3 pre-specification contract (E1)

FROZEN before any optimization output. (git commit = timestamp) Binds to all analyses in Study 2-3 for E2-E12. Governed by INTEGRITY.md, the complete execution plan, and the claim-control rules in INTEGRITY.md. Nothing in this section can be changed after commit except by an explicitly designated exploratory addendum, NEW. The main analysis in this section cannot be silently replaced (INTEGRITY.md).

## 0. What E0 established (the inputs to this freeze)
- CPC18 assignment is 7 fixed problem forms keyed on `Set` (largest = Set 1, 125 subjects x 24 clean problems); no clean problem is seen by >= 150 subjects, so there is NO common held-out fixed form. Study-3 validation path = SIMULATION-BACKED (primary), anchored by within-form real-data validation where a form is large enough (Set 1).
- Only 2 of 5 constructs have trustworthy, positive between-person trait variance on the metric scale: loss_avoid_adj (ceiling 0.545 at 14 clean problems; ~27 problems for reliability .7) and ev_sens_agn (ceiling 0.281; ~84 for .7). The 3 Firth coefficients' metric ceilings are NOT TRUSTWORTHY at CPC18 budgets (quasi-separation); their reliability question stays on Study 1's rank scale.

## 1. Targets (LOCKED)
- PRIMARY optimization target: loss_avoid_adj (metric scale; highest trustworthy ceiling).
- SECONDARY optimization target: ev_sens_agn (metric scale).
- b_ev, b_risk, b_loss: carried on the RANK scale only, as secondary robustness (per Study 1). They are NOT primary optimization targets; no metric-scale reliability claim is made for them (E0: metric ceiling untrustworthy).
- Rationale committed: E0 recovery-validated trait variance for the two indices and showed the coefficient metric estimator fails recovery. Scope narrows from five constructs to two by evidence, not by convenience.

## 2. Estimand (LOCKED)
For a target index and a problem SET of size K (K clean problems), the estimand is the between-person split-half reliability (Spearman, Spearman-Brown corrected, the frozen Study-1 engine R/step2_reliability.R) that the set achieves in held-out participants. The scientific question: does an information-optimized set of K problems achieve HIGHER held-out reliability than random and stratified-random sets of the SAME K, and is it NON-INFERIOR to a random set of 2K? Co-primary with reliability: parameter recovery of the target under the design (simulation). Secondary: held-out choice prediction on problems NOT in the set.
- Scope caveat (LOCKED, from E0): reliability is defined within CPC18's item universe as sampled by its 7 forms, not an arbitrary item universe. Ceilings bound achievability: no design can exceed the true-variance ceiling (loss_avoid_adj 0.545, ev_sens_agn 0.281 at 14 problems), so "more reliable" means relative to random at matched length, within that ceiling, NOT reaching conventional .70.

## 3. Participant split (LOCKED; the actual seeded assignment is committed with this contract)
Three-way, COHORT-RESPECTING: within EACH form (Set 1-7), participants are randomly assigned build / freeze / untouched-test in proportions 0.50 / 0.20 / 0.30, with pre-declared seed 20260718 (stratified by Set so every form contributes to every split). Analysis pools each split across forms.
- BUILD: develop the information-feature mechanism model (Study 2) and fit the generative choice model for simulation. All exploration happens here.
- FREEZE (selection): select the final optimized designs / selection rule. No tuning after this closes (E6).
- UNTOUCHED-TEST: reserved for E7 confirmatory validation ONLY. Its reliability, recovery, and prediction outcomes are NEVER computed or inspected before E7. This contract defines the group and locks it; reading it early voids the pre-specification.

## 4. Metrics (LOCKED)
- Reliability: split-half Spearman, SB-corrected, subject-cluster percentile bootstrap CI (frozen step2 engine), experience regime (Trial 6-25, r=20), on the target index.
- Recovery (co-primary): parameter recovery of the target index under the design on simulated participants drawn from the build-group generative model (labeled SIMULATION); report correlation + RMSE vs planted.
- Prediction (secondary): held-out choice log loss and Brier score on problems NOT in the design, vs population-only and vs longer-random baselines.

## 5. Baselines (LOCKED; equal sampling/search budget for every design criterion)
random-K (uniform), stratified-random-K (random within pre-declared feature strata of d_ev / d_sd / d_eloss), and random-2K (double length, for the efficiency claim). Problem-count budgets B (total clean problems): {6, 8, 10, 12}, operationalized within each form's feasibility in E3; the "half as long" comparison pits optimized-B against random-2B where 2B <= the form's clean-problem count (feasible in Set 1, 24 clean).

## 6. Success thresholds (LOCKED; evaluated at E7 in the UNTOUCHED test only)
Merge into the main paper as CONFIRMATORY only if ALL core conditions hold in the untouched test:
- C1 gain: optimized-B median held-out reliability > random-B median, at the pre-declared B.
- C2 support: a bootstrap interval of at least 90% supports (optimized-B reliability - random-B reliability) > 0.
- C3 recovery preserved: optimized-B parameter recovery is not worse than random-B by more than 0.05 (correlation units).
- C4 no validity cost: optimized-B held-out prediction (log loss) is not worse than random-B by more than 1% relative.
The stronger "half as long AND more reliable" headline (H4) ADDITIONALLY requires:
- C5 non-inferiority: optimized-B reliability CI lower bound >= random-2B median reliability - delta, with delta = 0.05 (SB reliability units), pre-declared.
If the core (C1-C4) does NOT hold: keep Paper 1 unchanged and report the extension as EXPLORATORY (honestly labeled, still in the paper). A failed optimization is a valid scientific outcome (INTEGRITY.md).

## 7. Frozen hypotheses (LOCKED)
- H1 same-count heterogeneity: at fixed K, held-out reliability varies materially across item sets (premise gate; if item-set variance is trivial, H3-H4 are moot).
- H2 mechanism: pre-declared information features (target-predictor variance, loss exposure, predictor collinearity, expected information, choice entropy) predict held-out reliability beyond raw count and broad problem category, under grouped cross-validation.
- H3 matched-length superiority: optimized-B > random-B AND > stratified-random-B in the untouched test.
- H4 efficiency (FRAGILE; requires C5): optimized-B non-inferior to random-2B. This is the "half as long" claim; it is the one most likely to fail and is reported honestly either way.
- H5 no validity cost: optimization does not worsen held-out prediction (C4) or recovery (C3).
- H6 robustness: advantages survive alternative seeds, alternative parameter distributions (including a LOW-variance population), a secondary estimator, and leave-one-form-out cross-form generalization.

## 8. What this contract does NOT claim (limitations frozen up front)
- Not "reaches conventional .70": neither target reaches .70 within a single current form; the claim is relative efficiency (more reliable at matched length / non-inferior at half length), within the true-variance ceiling.
- Not external: validation is internal to CPC18 (Study 4 scratched); a prospective external follow-up is only proposed in Discussion, never run.
- Not adaptive / not the invention of information-optimal design: the claim is fixed/frozen information-optimized forms validated for between-person reliability on CPC18, per the novelty framing in the execution spec.
- Precision != reliability: an information criterion maximizes within-person precision; reliability also depends on population heterogeneity, so reliability is claimed as a validated empirical outcome, never as something the criterion guarantees.

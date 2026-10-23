# SENS - mechanism sensitivity and detection curve (SPEC; committed BEFORE any result)

Pre-specification for the sensitivity analysis of the Study 2 mechanism null (frozen H2 FAIL, E4). Purpose: calibrate WHICH feature-to-reliability mechanisms the E4 design could and could not have detected under the empirical set-level outcome noise (E3b) and the expected-information / loss-exposure collinearity (~0.90), so the null is stated as informative or as detectability-limited rather than left ambiguous. Build-only. The sealed 205 test group is NEVER read. This spec is committed before results exist; its gates are not amended after results.

## Scientific question
Given a KNOWN, planted feature-to-reliability mechanism of graded strength, what is the probability the frozen E4 gate would detect it, under the empirical outcome noise, feature collinearity, and seven-form structure of CPC18 BUILD?

## Data and engines
BUILD group only (344 subjects). Reuse, by sourcing (not redefining): the frozen Study-1 scoring engine, `R/cpc18_features.R`, `simulate_recovery.R`, the E3 set-sampling and reliability machinery, the E3b calibrated generator and its noise model, the E2 problem-information dictionary (F1-F5), and the frozen E4 grouped-by-form cross-validation gate. No raw freeze/test choices; freeze not read; 205 sealed.

## Design (Phase 1)
1. Planted-mechanism generator: extend the E3b build-calibrated generator so each problem's discrimination on the latent loss-avoidance trait depends on its information features with coefficient beta (beta = 0 means features do not affect discrimination; larger beta means more-informative problems discriminate individual differences more). Preserve the empirical choice rates, participant heterogeneity, item distribution, seven-form structure, and the observed feature covariance including corr(F2,F4) ~ 0.90.
2. Effect grid (PRE-DECLARED, fixed): beta in {0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40}. For each beta and each of N_sim replicates (pre-declare N_sim; report Monte Carlo error): generate data, sample equal-length sets as in E3, compute each set's split-half reliability carrying the empirical noise, run the FULL frozen E4 gate (grouped-by-form out-of-fold Spearman, sign-consistency, beyond raw-count and category).
3. Per beta, record: the induced feature-to-reliability correlation r_true; the induced set-level reliability change (to locate the "matters" effect); P(correct sign); P(out-of-fold Spearman exceeds the frozen E4 threshold); P(pass the complete gate). Build the detection curve P(pass) vs beta.
4. Attenuation bound: compute observed_ceiling(r_true) = r_true * sqrt(rho_o) using the E3b out-of-subject rho_o at each budget; report the max observable association and the r_true required to reach out-of-fold 0.20 and 0.30 at K=12. Mark the honest real-data association on the same axis.
5. Collinearity scenarios: repeat the detection curve at corr(F2,F4) in {0.5, 0.7, 0.9, 0.95} (re-mix / orthogonalize the two features), to show whether the design can attribute reliability to information specifically as collinearity rises toward the observed value.
6. Feature-definition robustness: repeat the real-data gate verdict and the minimum detectable effect under the alternative per-item F1/F2/F3 definitions ratified in the P2 erratum, to show the null and the detectability floor are not artifacts of one operationalization.
7. Uncertainty: repeated / nested cross-validation for honest uncertainty on the out-of-fold Spearman without treating overlapping item sets as independent.
8. Controls: POSITIVE control = a large planted mechanism (top of the grid or a direct strong reliability effect), which must be detected; NEGATIVE control = a shuffled / irrelevant feature, which must not be detected.

## Pre-declared "matters" threshold
A feature effect matters for task design if it would change set-level reliability by at least 0.05 SB units at K=12. Fixed here. The smallest beta whose induced reliability change reaches this is the design-relevant ("matters") effect.

## Pre-declared gates
- G-power: the design detects the design-relevant effect with power >= 0.80 under empirical noise and collinearity, stable across feature definitions -> the H2 null is INFORMATIVE (design-relevant mechanisms are excluded).
- G-detect-floor: report the minimum detectable effect (smallest true effect detected at power 0.80). If it exceeds the design-relevant effect -> the null is DETECTABILITY-limited; the conclusion narrows to "cannot confirm or exclude a mechanism smaller than the minimum detectable effect."
- G-controls: positive control passes and negative control fails, else the run is void.

## Stopping rule
The effect grid, N_sim, the "matters" threshold, the gate, the collinearity grid, and the feature-definition set are fixed before running and are not expanded post hoc. If constructing the generator requires misspecification choices, report sensitivity to them; if a conclusion depends on a contestable choice, demote the deliverable to the attenuation bound alone (which needs no generator). This analysis estimates DETECTABILITY only; it does not fit a new real-data mechanism and does not optimize anything.

## Exports (elevation/results/sens/)
`sensitivity_detection_by_beta.csv` (per beta: r_true, induced reliability change, P(sign), P(threshold), P(gate), MC error), `attenuation_bound.csv`, `collinearity_scenarios.csv`, `feature_definition_robustness.csv`, `sens_config.json`, `sens_seeds.csv`. Plus the run report and the session-log row.

## Integrity
build_only=TRUE, freeze_read=FALSE, test_read=FALSE (sealed 205 never read). Additions-only under elevation/ (no frozen Study-1 file or prior elevation result subtree touched). Frozen engines sourced, not redefined. Simulation labeled SIMULATION throughout. A pre-declared gate returning a detectability-limited null is a valid, successful analysis and must not be repurposed to pass.

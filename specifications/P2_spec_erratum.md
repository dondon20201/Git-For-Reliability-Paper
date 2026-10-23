# P2_spec.md erratum (post-hoc; P2 review, 2026-07-20)

This note details problems observed in the P2 pre-specification (elevation/P2_spec.md, dc637a2) during running P2 and during the run as well as problems that have been independently verified. Integrity of the pre-specification requires that the checked in pre-spec be not edited. This note details problems that were observed and the actions that were taken. There is no impact on the outcome of P2 (the locked H2 FAIL remains and is enhanced; no mechanism is recovered).

## 1. The spec's EV control `gamma` is not separately identifiable
Section 2 specifies logit(p_ij) = beta0_j + lambda_j * u_i + gamma * d_ev_noloss_ij with FREE per-item easiness
beta0_j. Because d_ev_noloss is a problem-level constant, the gamma * d_ev_noloss term is exactly collinear
with the item easiness (rank(item dummies + d_ev_noloss) == rank(item dummies) == 42). The model therefore
reduces algebraically to the one-factor 2PL with the EV control absorbed into beta0_j. This reduction is the
correct object for a mechanism-on-loadings analysis; it was fitted, and for reference the E8-style
common-intercept variant was also reported, which DOES identify gamma = 0.165108 (matching E2's build_loss_baseline d_ev_noloss
0.1651; the independent baseline 0.1646). Scientific impact: none.

## 2. The G3 calibration threshold sits below its own permutation-null floor
Section 4 sets G3 to pass at Spearman(rho_hat, E3b) >= 0.30, a fixed constant. The added null-floor
diagnostic shows the feature-permutation null's OWN median calibration is 0.316 / 0.390 / 0.445 at K = 8 / 10 /
12 - above 0.30 at every gating budget - so a model with no common trait passes G3 as written. G3 is therefore
recorded PASS-but-VOID (no evidential value). The correct standard would have been a permutation-relative
threshold. Verdict impact: none - G2 is the mechanism gate and it is null; G3 was only the loadings ->
reliability bridge, moot once the mechanism link is absent. (Note: the real calibration does exceed the null's
95th percentile and the partial-given-count is a stable ~0.53, so the loadings carry weak above-chance RANK
information; but the implied reliability LEVEL is uninterpretable - HG1b: 100% of sets above the E0 ceiling,
driven by ~8.5x overdispersion - so it is NOT treated as a mechanism result.)

## 3. F1 / F2 / F3 per-item feature definitions were analysis-side gap-fills (RATIFIED)
E4's set-level formulas for F1 (target-predictor spread), F2 (loss exposure), and F3 (predictor collinearity)
do not decompose to a single item; reasonable per-item analogues were used (documented in p2_config.json).
Verification RATIFIES them: the H2 mechanism null is robust to these definitions because F4 (expected/Fisher
information, the H2 star feature) and F5 (choice entropy) DO decompose exactly and already show no signal
(marginals +0.033 and +0.163 on the primary loadings; +0.016 and +0.152 on the independent re-fit).

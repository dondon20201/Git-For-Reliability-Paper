# Erratum to SENS_spec.md

This erratum describes three issues in the pre-specification for the SENS analysis that resulted in the detection curve deliverable being VOID due to a pre-specified G-controls condition (positive control was detected with probability 0.000 against a required 0.80 bar). These issues were made during the pre-specification prior to result availability. These issues were never corrected after the fact. Instead the run was placed under its pre-specified backup (attenuation bound) in accordance with the pre-specified stop condition. This investigation does not impact the frozen H2 FAIL.

## Defect 1 - the beta grid was never calibrated in reliability (SB) units

The grid beta in {0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40} was fixed before
anyone knew what reliability change a given beta buys. Under the committed
parameterization (each problem's discrimination scaled by 1 + beta * g_j, with
g_j the standardized per-item information) the entire grid is sub-threshold: the
largest induced change at K = 12 is 0.0143 SB units against the pre-declared
0.05 "matters" bar, and the underlying true F4-to-reliability Spearman moves
only from 0.535 to 0.565 across the whole grid. Hence beta_matters = NA and the
"large planted mechanism" the positive control assumes never existed. The deeper
reason the dial barely turns: the baseline association is already 0.535 at
beta = 0, because the frozen choice model by itself implies that more-informative
item sets are more reliable, and set composition adds more on top. The planted
term is a third layer over two that already saturate the achievable gradient.

Fix for any future attempt: parameterize the plant DIRECTLY in SB units,
calibrating beta against the high-precision arm BEFORE fixing the grid, so the
grid spans roughly 0.02 / 0.05 / 0.10 / 0.20 induced change and the positive
control is genuinely large. The concrete lesson is one of sequencing: calibrate,
then freeze the grid.

## Defect 2 - there is no clean null world to calibrate detection against

At beta = 0 the generator equals the frozen model exactly (G_GEN-proven, max abs
diff 6.66e-16). But that world is NOT mechanism-free: it carries the information
structure the frozen choice model itself implies, plus composition effects. With
no additional planted coupling, the frozen gate's two predictive clauses still
fire in 100 of 100 replicates (pooled statistic 0.699, delta over the loss-count
baseline 0.494 at K = 12, bootstrap support 1.000). The pre-declared 0.30
threshold therefore does not discriminate planted-mechanism-present from
planted-mechanism-absent on this design; both sit far above it. This is the same
defect class the P2 erratum recorded for G3 (a threshold sitting below its own
null floor), now visible in E4's predictive clauses. It does not overturn the H2
verdict, which was driven by G-real (the clause with evidential content); it
means the two clauses that passed in E4 (G-sim, G-beyond) had little
discriminating power to begin with.

## Defect 3 - the sign clause is not a usable detector under collinearity

I mapped the frozen G-real sign-consistency clause onto "the expected-information
coefficient is a carrier and points in the planted direction". With
corr(F2, F4) = 0.937 at K = 12, F4 acts as a suppressor and its partial
coefficient is negative essentially always (mean -0.266 at beta = 0.40, positive
in only 4 of 100 replicates). The decisive evidence that this is parameterization
and not truth: on the SAME committed real form-1 reliabilities the F4 coefficient
is -0.157 under the canonical feature definition and +0.085 under the ratified
P2-erratum alternative at K = 12. A clause whose sign flips between two
operationalizations both already ratified cannot certify a mechanism.

Fix: replace the sign clause with a composition-matched or permutation-relative
criterion, which is the same fix the P2 erratum recommended for G3. Note that
this changes the detection target from "the frozen gate as written" to "a gate
with a defensible null floor", which is a different scientific question and must
be named as such if pursued.

## Disposition

Per the committed stopping rule ("if a conclusion depends on a contestable
choice, demote the deliverable to the attenuation bound alone, which needs no
generator"), the deliverable was demoted to the attenuation bound. That bound
stands and is generator-free: the committed real-track association exceeds its
own out-of-subject reliability ceiling at every budget (0.745 vs 0.689 at
K = 12), so the excess is shared-subject and shared-item fit, not transferable
mechanism signal. The H2 null is reported WITHOUT a detectability calibration,
WITH the attenuation bound as the ceiling on what any mechanism could have shown
at these noise levels. No re-run with a retuned grid was performed; moving a
pre-declared bar after seeing results is forbidden by INTEGRITY.md and the
spec's own stopping rule.

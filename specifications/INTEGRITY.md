# Integrity rules

These are the current methodological commitments of this project. These
commitments predate the analysis and were committed to the development repository
before the results they describe. These commitments were not
relaxed at any point during the project.

This file is referenced by the pre-specification files in the `elevation/`
folder. The application of each rule is tracked in the appendix to the manuscript
that details the status of each claim.

These commitments are documented for two reasons. A user who must assign a level
of trust to a result should be able to learn the standard to which a result
was held. Additionally, a rule specified prior to data observation is a
constraint. A rule specified after data observation is only a description.

## 1. All data are real and verified

Only genuine, publicly documented, citable open datasets from their official sources are
used. No invented data, no synthetic data presented as real, no illustrative numbers, and
no placeholder values that later harden into findings.

Every dataset is downloaded from its official source and verified on arrival: md5 checksum
where one is published, row, subject and column counts, and reproduction of at least one
basic summary already known from the source. A checksum or structure that does not match
expectation is a stop condition, reported rather than worked around.

Simulated data appear only in parameter-recovery validation, are labelled unmistakably as
simulation wherever they appear, and are never reported as an empirical result.

## 2. No number is fabricated or recalled

Every quantitative statement in the manuscript comes from code run on the real data. No
number is taken from memory, from a planning document, from a previous run, or from what
the literature approximately says.

This rule is enforced structurally rather than by discipline. Prose in
`manuscript/sections/` may not contain a hand-typed numeral; every value is emitted by an
inline `num()` lookup against `results/step5/manuscript_numbers.csv`, which is itself
extracted programmatically from the committed results files. An unknown key fails the
render, and a stray numeral outside the permitted carriers fails the render. The
manuscript therefore cannot state a number that is not present in a committed results
file. `manuscript/README.md` documents the mechanism.

Where a value could not be computed, the gap is reported as a gap.

## 3. Claims are calibrated to the evidence

Claims state what the data support, with the hedges the data require, and prefer the
scope actually studied over generalisation to the field. Every estimate carries an
uncertainty interval; there are no naked point estimates.

Suggestive results are not upgraded to definitive ones. Consistency with a hypothesis is
not reported as demonstration of it. A single dataset is not reported as the field. Where
a result is fragile, unstable or borderline, that is stated in the same breath as the
result rather than deferred to a limitations section.

A smaller claim that survives scrutiny is preferred to a larger one that does not.

## 4. Every claim carries a status label

Each quantitative claim is labelled VERIFIED (re-derived from the real data in code that
is shown), REPRODUCED (matched against a named external source), EXPLORATORY, SIMULATION,
UNVERIFIED, or FAILED (tested, and did not hold). The full table is the appendix at
`manuscript/sections/07-claim-status.qmd`, rendered from
`results/step5/claim_status_skeleton.csv`.

Did not replicate, did not converge, too unstable to report, and the effect is not present
are complete and acceptable results. An honest null is an outcome, not a shortfall, and is
never presented as though it were a positive finding.

## 5. Analysis is pre-specified, and pre-specified failures stand

Each analytical step is specified in writing and committed before the results it governs
exist. The specifications are in `elevation/`, and
`elevation/results/e11/provenance.md` tabulates, for each one, the commit that introduced
the specification against the commit that introduced its results.

A threshold declared in advance and then missed is reported as a missed threshold. Such
outcomes have occurred in this project and are reported as results. Thresholds are not
revised after the fact to accommodate an outcome.

This is pre-specification, not pre-registration. Commit timestamps come from the author's
own clock and version history can be rewritten, so the record evidences ordering to a
reader who accepts the remote's push history; it is not an independent registry, and the
manuscript does not claim that it is. The specification documents originally used the
word "pre-registration" and were corrected at release; see
`elevation/TERMINOLOGY_ERRATUM.md`.

## 6. Estimators are validated before they are trusted

Every estimator must be shown to recover a known planted value before it is applied to
real data, and that gate precedes all reliability claims. Reliability, precision and
parameter estimates carry bootstrap intervals. Hierarchical and Bayesian results report
Rhat and effective sample size, and an Rhat at or above 1.01 is treated as untrustworthy
rather than borderline.

## 7. Results are generated and graded separately

The type of error this project is primarily concerned with is the type in which
an identical process generates a result, validates that result, and accepts it,
thus allowing errors that are flattering to the hypothesis to persist. As a result,
analysis and validation were separated. Important measures were independently
re-derived using the same original data but in different code and by an entity
that had not seen the original code. These derivations were then reconciled
prior to the acceptance of any result. Deterministic measures were required to
be identical while stochastic measures were required to match in expectation.

Alongside that, results were actively attacked rather than confirmed: re-run under
different splits and seeds, checked for recovery of known values, and examined for what a
reader trying to falsify them would reach for first. What that turned up is reported,
including where it weakens the account.

The repository is not optimised for an impressive paper. It is optimised for claims that
survive adversarial re-analysis. Where honesty and impressiveness conflict, the rule is
honesty, said plainly.

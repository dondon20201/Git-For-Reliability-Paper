# Integrity rules

These are the standing methodological commitments of this project. They were written
before the analysis began, committed to the repository ahead of the results they govern,
and were not relaxed at any point afterwards. The pre-specification documents in
`elevation/` cite this file, and the claim-status appendix of the manuscript is the
audit trail showing where each rule was applied.

They are recorded here for two reasons. A reader deciding how much weight to place on a
result is entitled to know what standard the result was held to. And a rule that is
written down before the data are seen is a constraint; a rule described afterwards is
only a description.

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

The failure this project is most concerned with is the one where the same process
produces a result, checks it, and passes it, so that errors flattering to the hypothesis
survive. Analysis and verification were therefore kept apart: key quantities were
independently re-derived, from the same raw data but in different code, by a party that
had not seen the original implementation, and the two were reconciled before any result
was accepted. Deterministic quantities were required to match exactly; stochastic ones
were required to agree at the expectation level.

Alongside that, results were actively attacked rather than confirmed: re-run under
different splits and seeds, checked for recovery of known values, and examined for what a
reader trying to falsify them would reach for first. What that turned up is reported,
including where it weakens the account.

The repository is not optimised for an impressive paper. It is optimised for claims that
survive adversarial re-analysis. Where honesty and impressiveness conflict, the rule is
honesty, said plainly.

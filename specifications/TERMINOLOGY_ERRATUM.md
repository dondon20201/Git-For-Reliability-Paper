# Erratum: pre-registration corrected to pre-specification

Date of correction: 2026-07-31, at public release of this repository.

## What was wrong

The specification documents in this directory, several analysis scripts, and the
generated provenance record described this project's practice as *pre-registration*, and
described individual specifications as *pre-registered*. That wording appeared 41 times
across 22 files.

This is not correct. None of the specifications in this project were registered in a public
registry. Rather, this project wrote the specifications, hypotheses, and pass/fail criteria,
and checked them into version control before running the code that implemented them.
This is pre-specification. This is weaker than registration and should not be given the name
of registration.

The difference is that these two assertions are grounded in different things. A registration
is timestamped by a third party that is independent of the researcher and has an immutable
history in the registration. A commit is timestamped by the researcher’s own computer and
has a history that can be rewritten. However, to a person who will accept the development
repository’s push history as a witness, the commits provide evidence of their ordering.
They do not provide ordering to a person who will not accept that push history. It is not
the same as a registration.

## What was changed

Every occurrence of *pre-registration*, *pre-registered*, and their unhyphenated variants
was replaced with *pre-specification* and *pre-specified* in the specification documents,
the analysis scripts, and `elevation/results/e11/provenance.md`. The heading of the
timeline table in that file changed from "Pre-registration timeline" to
"Pre-specification timeline".

The manuscript was corrected earlier, on 2026-07-23, and already used the accurate
wording before this repository was released.

**Nothing else changed.** No threshold was altered, no analysis was re-run, no result
moved, and no verdict was revised. Every pre-declared gate stands exactly as it was
declared, including the ones that failed. The substitution was a wording correction to
documents whose substance is unchanged.

## Why the documents were edited rather than left alone

There is a reasonable argument for leaving a pre-specification document exactly as it was
written, on the grounds that editing a frozen record is itself a kind of revision. The
argument does not apply here, because the wording was not a claim that was true when made
and later became false. No registration ever existed. Correcting the word removes an
inaccuracy; leaving it would have preserved one, and would have put these documents in
direct contradiction with the manuscript that cites them.

This erratum exists so that the correction is on the record rather than silent.

## Known remaining occurrences

A small number of descriptive string fields inside committed results files still carry
the old wording, because those strings are emitted by the analysis scripts and the two
must change together. Correcting them requires re-running the scripts that produce them,
which is a separate exercise and is not done here rather than hand-edited, so that every
results file in this repository remains exactly what its script produced. The affected
fields are notes and dispositions; none is a number, and none appears in the manuscript.

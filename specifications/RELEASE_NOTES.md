# Release notes for the specification documents

The files in this folder are copies of the specification documents from the
development repository, where each one was committed before the results it
governs (see `../provenance/provenance.md` for the dates).

Before release, two kinds of change were made to these copies:

1. Internal project-management notes were removed or replaced: run numbers in
   titles and cross-references (for example "dispatch 037", now "E4"), the
   labels "Dispatch B" and "Dispatch C" (now "Battery B" and "Battery C"),
   working-session labels, and formatting rules for the project's own reports
   (dash and quote characters, wording of internal notes).
2. Two sentences in `P2_spec_erratum.md` and one in `P2_spec.md` were reworded
   so that they read correctly after change 1.

No hypothesis, threshold, gate, estimand, sample definition, seed or decision
rule was changed. The earlier wording correction from "pre-registration" to
"pre-specification" is described separately in `TERMINOLOGY_ERRATUM.md`.

Because of these changes, the copies here are not byte-identical to the
versions committed in the development repository. The analysis content is
identical.

File paths inside the documents refer to the layout of the development
repository. The table in `../ANALYSIS_MATERIALS.md` maps each specification to the step of
this repository that carries it out. At integration, historical freeze tag labels
were replaced with references to the same freeze points or commit IDs in five
places. The reported analysis sequence and thresholds were not changed.

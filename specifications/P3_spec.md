# P3 external replication - specification (Frey 2017 Basel-Berlin Risk Study). Does the wall generalize?

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md. This is the P3 leg of the approved third-act package (P1/E8 done -> P3 now ->
P2/IRT later).

P3 asks whether Study 1's reliability wall - single risk-preference measures are unreliable while the
aggregate trait is reliable - REPLICATES in an independent dataset (Frey et al. 2017, N=1507, plus a
6-month retest), and whether E8's null (better estimation does not rescue reliability) holds there too. It
touches NO CPC18 repo analysis, does NOT read the sealed 205 CPC18 test subjects, and does not modify the
frozen Study 1 or Study 2.

## 0. What P3 is (and is NOT)
- P3 is a re-analysis of the PUBLIC Basel-Berlin Risk Study data (OSF rce7g), a local copy of which the author
  placed at `data/raw/Frey Data/` (gitignored via the existing `data/raw/` rule; NEVER committed - assert
  `git check-ignore` succeeds, and if for any reason it does not, STOP and report rather than editing
  .gitignore).
- P3 REPLICATES two things established on CPC18: (A) Study 1's message - single behavioral risk-preference
  measures are unreliable, the aggregate is reliable; (B) E8's method-level finding - sophisticated
  (model-based / hierarchical) estimation does not rescue single-measure reliability, because the limit is a
  population true-variance ceiling, not estimation noise.
- P3 is NOT a literal same-estimand replication. Frey's headline reliability is 6-month TEST-RETEST
  (temporal stability) of single summary measures; Study 1's is SPLIT-HALF (internal). P3 reports BOTH and
  never conflates them; showing the wall holds under temporal stability is a STRENGTHENING, stated as such.
- P3 changes nothing in the CPC18 line. Frozen Study 1 (freeze commit 387f9adc), Study 2
  (H2 FAIL), and E8 stand untouched. The 205 CPC18 test subjects stay sealed.

## 1. Inputs (local public data; gitignored; provenance recorded)
Base path `data/raw/Frey Data/`. Three samples: `main/` (N=1507), `retest_berlin/` (n=109, 6-month retest,
primary), `retest_basel/` (n=64, 3-month, secondary robustness only). Key files:
- `<sample>/measures/measures.Rdata` - objects: `measures` (subjects x 238 columns: the risk measures +
  demographics), `grouping` (180 labels -> `group` in {beh, prop, freq} + task/subtask), `participants`
  (IDs), `cbalpha` (Cronbach alpha for the 53 representative measures Frey used).
- `<sample>/dfd/dfd_perpers.csv` - partid, R (proportion risky over 8 description gambles), R_sh1, R_sh2
  (Frey's pre-computed split-halves). `dfd/dfd_perprob.csv` - 8 rows/subject: partid, gamble_lab, domain
  (gain/loss), gamble_ind, H, R, decision (A/B) - the problem-level choices.
- `<sample>/dfe/dfe_perpers.csv` - partid, samples, switches, R, Rexp. `dfe/dfe_perprob.csv` - problem-level
  DfE choices (compute a DfE split-half from these). `dfe/dfe_samples.csv` - raw sampled outcomes.
- `<sample>/modelpars/modelpars.csv` - per-subject fitted choice-model parameters for the lottery / MPL /
  CCT tasks (ABCT/MPL/CCT families x EU/4P/6P models; parameters include outcome sensitivity OS, probability
  weighting PW, loss aversion LA, choice sensitivity CS). The external analog of CPC18's per-subject params.
- Participant linkage (LOCKED, independently verified): main IDs end "01", retest IDs end "02"; dropping the last
  two digits yields a 6-digit core that links retest -> main (108 of 109 Berlin retest subjects link; the
  one non-linking id is reported, not force-matched).
- Provenance: record the OSF project (rce7g), the file list used, each file's md5, and N per sample into a
  committed manifest. The raw data itself is never committed.

## 2. Targets and estimands (LOCKED)
- Representative measure set = the 53 `cbalpha` measures, classified into traditions by `grouping$group`
  (beh / prop / freq). CPC18 analogs flagged explicitly: DFD, DFDr (Decisions from Description) and DFEre,
  DFEr (Decisions from Experience).
- Reliability estimands:
  - TEST-RETEST (temporal): Spearman correlation of a measure at T1 (main) vs T2 (retest_berlin) across the
    linked ~108 subjects, per measure and aggregated by tradition. (Frey's headline estimand.)
  - SPLIT-HALF (internal, Study-1 style): for DfD, Spearman(R_sh1, R_sh2) across N=1507, Spearman-Brown
    corrected; for DfE, compute two disjoint problem-halves from `dfe_perprob` and do the same.
  - GENERAL FACTOR: a simple, pre-declared general risk factor - the first principal factor (principal-axis
    or first PC) across the standardized representative measures, fit on T1; score T2 with the T1 loadings;
    test-retest = Spearman(T1 score, T2 score). (A transparent aggregate, NOT a re-fit of Frey's bifactor
    model; Frey's published bifactor R = 0.85 is the external benchmark.)

## 3. Analysis (LOCKED)
Part A - replicate Frey's reliability structure (temporal):
- A1: per-measure and per-tradition 6-month test-retest (section 2). Report all 53; aggregate beh/prop/freq.
- A2: general-factor test-retest (section 2). Report loadings + score reliability.
- A3: context - Frey's shipped Cronbach alphas (cbalpha) by tradition (internal consistency), reported
  alongside, not re-derived.

Part B - replicate Study 1 + E8's method on the CPC18-analog tasks:
- B1 (Study-1 wall, external): DfD and DfE single-index reliability, BOTH split-half (internal) and 6-month
  test-retest. Pre-declared prediction: LOW (well below the 0.70 convention; the DfD/DfE test-retest is
  particularly low).
- B2 (parameter wall, external): 6-month test-retest of the model-based risk-preference parameters
  (modelpars: loss aversion, probability weighting, outcome sensitivity for the lottery/MPL/CCT tasks).
  Pre-declared prediction: LOW for single parameters - the CPC18 per-subject-parameter noise reproduces on
  independent data.
- B3 (E8 parallel - does better estimation rescue it?): for the DfD task, compare the reliability of a
  SOPHISTICATED per-subject estimate to the RAW proportion. Sophisticated = a hierarchical partial-pooling
  logistic on the 8 problem-level choices, `decision ~ 1 + domain + (1 | partid)` (partial-pooled per-subject
  risk propensity = random intercept), split-half by fitting each 4-problem half separately and correlating
  BLUPs (SB-corrected); raw = Spearman(R_sh1, R_sh2). Report Delta = rel_hier - rel_raw with a bootstrap CI.
  Pre-declared prediction (from E8): Delta ~ 0 - better estimation does not rescue reliability. Caveat
  LOCKED: DfD is only 8 one-shot choices (4 per half), so B3 is thin and reported as the E8-parallel with
  that limitation stated; B1/B2 are the robust Part-B claims.
- B4 (the contrast): the general factor (A2, reliable ~0.8) vs the single-measure lows (B1/B2) - the
  aggregate-reliable / single-measure-noisy structure that is Study 1's decomposition, now external.

## 4. Pre-declared gates (committed before results)
- G-repro (Part A integrity + external replication): the recomputed per-tradition test-retest reliabilities
  reproduce Frey within tolerance - behavioral tradition mean in [0.36, 0.56] (Frey 0.46); propensity mean
  and frequency mean each strictly GREATER than the behavioral mean and in [0.52, 0.78] (Frey 0.68 / 0.65);
  general-factor test-retest > 0.75 (Frey bifactor 0.85). If these fail, the harmonization/linkage is wrong -
  STOP and diagnose; do not proceed to Part B.
- G-wall (PRIMARY replication): the DfD and DfE single-index reliabilities are LOW - split-half SB < 0.70 AND
  6-month test-retest < 0.35 for the DfD/DfE risk indices - and each is materially below the general factor.
  (External + temporal generalization of Study 1's wall.)
- G-param (Part B2): the single model-based risk-preference parameters (loss aversion, probability weighting,
  outcome sensitivity) have low 6-month test-retest (report each; expected mostly < 0.6), corroborating the
  CPC18 per-subject-parameter result.
- G-nolift (Part B3, the E8 parallel): hierarchical partial-pooling does NOT raise the DfD split-half
  reliability by >= 0.05 with a >= 90% bootstrap interval on Delta excluding 0. Pre-declared prediction: it
  does not (Delta ~ 0), mirroring E8. If it DID lift materially, that is a genuine external DIFFERENCE from
  E8, reported honestly - not suppressed.
- Interpretation (LOCKED, both directions valid): the expected pattern (single measures unreliable
  internally AND temporally; parameters unreliable; hierarchical no rescue; aggregate reliable) is a clean
  external replication of Study 1 + E8. Any DEVIATION (e.g., a behavioral measure that is reliable, or
  hierarchical rescuing DfD) is a real finding reported as-is, never reshaped to fit. A trivial/constant
  measure can be "reliable" without capturing individual differences (Frey's own caveat) - flag any measure
  whose reliability is high only because its between-person variance is degenerate.

## 5. Outputs (all new, additions-only under elevation/)
- `elevation/analysis/p3_frey_replication.R` - the analysis (reads only `data/raw/Frey Data/`; no CPC18 data,
  no repo results except nothing).
- `elevation/results/p3_frey/reliability_by_measure.csv` - per measure: tradition, test-retest r, n, (DfD/DfE
  also split-half).
- `elevation/results/p3_frey/reliability_by_tradition.csv` - beh/prop/freq mean+median test-retest + the
  general-factor reliability; the Frey published benchmark alongside.
- `elevation/results/p3_frey/general_factor.csv` - loadings + T1/T2 score reliability.
- `elevation/results/p3_frey/dfd_dfe_reliability.csv` - split-half (SB) + test-retest for the DfD/DfE indices.
- `elevation/results/p3_frey/dfd_hier_reestimation.csv` - B3: rel_raw, rel_hier, Delta, bootstrap CI, plus
  null-permutation and ceiling sanity checks (E8-style honesty guards).
- `elevation/results/p3_frey/modelpar_reliability.csv` - B2 per-parameter test-retest.
- `elevation/results/p3_frey/provenance_manifest.csv` - OSF source, file md5s, N per sample, linkage rule.
- `elevation/results/p3_frey/per_subject_export.csv` - the harmonized per-subject T1/T2 values for the 53
  measures + DfD R_sh1/R_sh2 + hier BLUPs + modelpar T1/T2, so independent verification re-derives every reliability from
  committed values (generator != grader).
- `elevation/results/p3_frey/p3_gate.json` - all gates vs thresholds; PASS/FAIL plainly.
- `elevation/results/p3_frey/p3_config.json` + seeds; `elevation/logs/p3_frey.log`.
- the run report.

## 6. Housekeeping (separate, labeled; the one allowed top-level change)
Run `renv::snapshot()` to record the packages E8 installed (lme4, rstanarm) plus any P3 dependencies, as a
SEPARATE `chore:` commit touching ONLY `renv.lock`. This is the single permitted top-level change; the P3
analysis commit itself stays additions-only under elevation/. If snapshotting would touch anything other
than renv.lock, STOP and report.

## 7. Pre-declared PASS gates (Part D)
- D1 provenance: the OSF source + file md5s + N per sample recorded in the manifest; the raw data is read from
  `data/raw/Frey Data/` and asserted gitignored; N=1507 main and ~108 linked Berlin retest asserted.
- D2 no CPC18 contamination: P3 reads ONLY the Frey data; it does not read any CPC18 raw data, the CPC18
  participant split, or any elevation/results/e* file; the 205 CPC18 test subjects are untouched by
  construction.
- D3 estimands as specified: test-retest = Spearman T1-vs-T2 on linked subjects; split-half = frozen SB form
  (2r/(1+r)); general factor = the pre-declared simple factor (not Frey's bifactor); hier model exactly the
  section-3 formula.
- D4 outputs present: every file in section 5, including the per-subject export for independent re-derivation.
- D5 reproducibility: seeds committed; reliabilities reproduce from the committed per-subject export; the
  linkage rule reproduces the 108 matches.
- D6 additions-only + frozen integrity (CRYPTOGRAPHIC): the P3 feat commit's diff touches ONLY
  elevation/ (plus the session log and coordination directory in the docs commit); the separate chore commit touches ONLY renv.lock;
  the 13 frozen Study-1 subtrees (R, analysis, data, figures, manuscript, results, renv, ...) byte-identical
  to freeze 387f9adc EXCEPT renv.lock (intended snapshot); elevation/results e0/e1/e2/e3/e4/e3b/e8 subtrees
  IDENTICAL to the pre-040 tip (p3_frey the only results addition); master 0808abc + both tags untouched;
  branch pushed.
- D7 no forbidden computation: no CPC18 test/freeze read; no optimization; the frozen CPC18 analyses are not
  re-run or altered; P3 is clearly labeled an external replication, not part of the CPC18 pre-specification.
- D8 hygiene: the temporal-vs-internal distinction stated wherever reliability is reported.

## STOP
After D1-D8: STOP and report the per-tradition reliabilities (vs Frey), the DfD/DfE single-index lows
(split-half + test-retest), the model-parameter reliabilities, the B3 hierarchical Delta, and the
general-factor reliability, with the plain gate verdicts. Do not start P2/IRT (the next analysis, authored
after P3 review). The CPC18 frozen Study 1, Study 2 (H2 FAIL), and E8 stand untouched; the 205 CPC18 test
subjects remain sealed. Verification proceeds against an independently prepared answer key, re-deriving from
its own staged copy of the Frey data AND cross-checking against Frey's published values.

# E3b denoised reliability outcome - specification (EXPLORATORY addendum)

Specification written and committed before any result it governs exists. Governed by the integrity rules in INTEGRITY.md and by the frozen pre-specification in elevation/E1_contract.md.

## 0. Status and purpose (READ FIRST)
- - EXPLORATORY. The frozen pre-specification (E1) and the frozen result at the gate for H2 (E4
  = FAIL) are the confirmatory record. Nothing in this section re-runs, re-scores, or replaces
  the frozen result at the gate for H2. This addendum is intended for DIAGNOSING the failure of
  H2 by removing two artifacts pointed out in the review at E4: (i) the per-set reliability
  OUTCOME used in E4 is noisy (sampling error in subjects ~0.145 SB units at K=12 in Set 1,
  which is relatively large compared to the across-set true SD ~0.19), which degrades
  prediction and also destabilizes the partial coefficients whose signs were tested in the
  G-real clause; and (ii) the shared real anchor in E4 puts all 63 subjects in both train and
  test sets. E3b creates a lower-noise reliability result that can be evaluated honestly for
  the EXPLORATORY re-test in E4b. E3b does NOT fit any mechanism model (E4b) and does NOT
  undergo any optimization.
- BUILD-ONLY. Part A uses only the 63 Set-1 BUILD subjects (E1 split); Part B uses the frozen build
  generative model. FREEZE (137) is not read (reserved for E6); TEST (205) is sealed (E7). Reading build
  choices for exploration is squarely within BUILD's purpose (E1 contract section 3).
- Reuses the SAME item sets E4 modeled (from E3's committed per-set files), so E4b models a cleaner
  outcome on identical sets. Reuses the FROZEN engines unchanged.

## 1. Inputs (build-only; additions-only)
- `data/raw/All estimation raw data.csv` (assert md5 `38a26407a36c690be969f9d9d3fe3e5a`).
- `elevation/results/e1/participant_split.csv` (filter `split == "build"`; assert 344; Set 1 build = 63).
- `elevation/results/e3/reliability_by_set_real.csv` (the committed real Set-1 sets + game_ids, per K/target).
- `elevation/results/e3/reliability_recovery_by_set_sim.csv` (the committed sim free-composition sets +
  game_ids, per K/target).
- `elevation/results/e2/problem_information_dictionary.csv` (features, for the loss baseline join as E3).
- Reuse FROZEN unchanged (source, do not redefine): `R/step2_reliability.R` (build_subject_data,
  fit_loss_baselines, estimate_half, split_correlations, sb_correct, reliability_point),
  `R/cpc18_features.R`, `R/simulate_recovery.R` (draw_true_params, analytic_p_b, simulate_choices).
  Build-only loss baselines fit on BUILD trials (as E2/E3).

## 2. Part A: real Set-1 subject-split reliabilities (honest out-of-subject outcome)
Goal: two INDEPENDENT estimates of each set's reliability on disjoint subject halves, so E4b can (a) fit
on one half and evaluate on the other (removing shared-subject inflation) and (b) measure how noisy the
outcome is (its cross-half reliability = the ceiling on predictability).
- Partition the 63 Set-1 build subjects into two disjoint halves A (32) and B (31) by a pre-declared seed
  (section 5). The partition is FIXED across all sets and both targets (so halves are a property of
  subjects, not of sets).
- For each committed E3 Set-1 set S (K in {6,8,10,12}; targets loss_avoid_adj primary + ev_sens_agn),
  reusing S's exact game_ids from `reliability_by_set_real.csv`:
  - r_sb_A(S) = frozen `reliability_point` over N_split_half = 200 random halvings of S, across the
    subjects in half A ONLY (experience regime, r = 20, frozen estimate_half + split_correlations +
    sb_correct). r_sb_B(S) = the same across half B ONLY, with an INDEPENDENT halving seed (A and B share
    no randomness). NA where identifiability fails (frozen behaviour); record coverage per half.
- Commit per set: form, K, target, set_id, game_ids, r_sb_A, r_sb_B, n_subj_A, n_subj_B, cov_A, cov_B.
- Report per (K, target): cross-set Spearman(r_sb_A, r_sb_B) over identifiable sets = the OUTCOME
  RELIABILITY. This is the honest denominator: a mechanism model can at best predict half B to the extent
  half A already agrees with half B. Report it plainly; NO disattenuation is applied (project convention).

## 3. Part B: simulation high-precision reliabilities (near-noise-free outcome)
Goal: drive the SIMULATION outcome noise low so E4b can test whether the mechanism and its coefficient
signs are stable once outcome noise is removed (the simulation is where free composition + arbitrary
precision are available; E0-determined path).
- Reuse the committed E3 sim free-composition set game_ids (500 per K, per target) from
  `reliability_recovery_by_set_sim.csv`. Same frozen generative model as E3 (draw_true_params,
  simulate_choices); the frozen uniform-range population is unchanged (calibrating it to build is a
  separate future step, deliberately NOT done here so this isolates the noise+precision effect).
- Recompute per set the simulated split-half SB reliability with N_sim_subj = 2000 (vs E3's 300) and
  N_split = 100 (vs 50), n_reps = 20. Commit per set: K, target, set_id, game_ids, r_sb_hp (+ recovery_hp
  optional). Label SIMULATION.
- Report per (K, target): the noise reduction vs E3 (the estimation SD of the per-set reliability should
  shrink ~sqrt(2000/300) ~ 2.6x); confirm the high-precision point estimates track E3's at expectation
  (no systematic shift, only less noise).

## 4. Outputs (all new, under elevation/results/e3b/; additions-only)
- `elevation/analysis/e3b_denoised_reliability.R` - the build script (reuses frozen engines).
- `elevation/results/e3b/reliability_subject_split_real.csv` - Part A per set (r_sb_A, r_sb_B, ...).
- `elevation/results/e3b/outcome_reliability_real.csv` - per (K,target) cross-half Spearman + n.
- `elevation/results/e3b/reliability_highprec_sim.csv` - Part B per set (r_sb_hp, ...).
- `elevation/results/e3b/highprec_vs_e3_sim.csv` - per (K,target) noise-reduction + agreement vs E3.
- `elevation/results/e3b/e3b_config.json` + `e3b_seeds.csv`.
- the run report.

## 5. Pre-declared constants and seeds
- Regime experience, pos Trial 6-25 (r = 20). Targets loss_avoid_adj (primary), ev_sens_agn.
- Part A: subject partition seed = SEED_E3B = 20260721; N_split_half = 200; half A = 32 / half B = 31
  (largest-remainder); halving seeds = SEED_E3B + 1 (A) / SEED_E3B + 2 (B).
- Part B: N_sim_subj = 2000; N_split = 100; n_reps = 20; param-draw seed = SEED_E3B + 100 + K;
  choice-sim seed = SEED_E3B + 200 + K. All seeds recorded in e3b_seeds.csv.

## 6. Pre-declared PASS gates (Part D)
- D1 engine provenance: frozen estimate_half/split_correlations/sb_correct/reliability_point + generative
  functions SOURCED, not redefined (grep + one fixed-set spot-check vs an independent recompute).
- D2 build-only / no leakage: Part A uses only the 63 Set-1 build subjects (partitioned); no freeze/test
  subject enters; Part B uses the build generative model; 205 test sealed; freeze not read.
- D3 same sets: the Part A set game_ids equal E3's committed real Set-1 sets exactly; Part B game_ids equal
  E3's committed sim sets exactly (set-identity check).
- D4 outputs present: subject-split (r_sb_A, r_sb_B) + outcome reliability (cross-half Spearman); sim
  high-precision r_sb_hp + noise-reduction vs E3.
- D5 reproducibility: seeds committed; deterministic partition; stochastic quantities reproducible in
  distribution.
- D6 additions-only + frozen integrity (cryptographic): diff 771782d4 -> HEAD touches ONLY {the session log,
  the coordination directory, elevation}; 13 frozen Study-1 subtrees byte-identical; elevation/results e0/e1/e2/e3/e4 subtrees
  IDENTICAL (E0-E4 untouched); E1/E2/E3/E4 spec blobs unchanged; master + tags untouched; branch pushed.
- D7 no forbidden computation: NO mechanism model (that is E4b), NO optimization, NO freeze/test read;
  simulation labeled SIMULATION; the frozen H2 gate is NOT re-run or amended.
- D8 hygiene: EXPLORATORY labeled.

## STOP
After D1-D8: STOP. E4b (exploratory mechanism re-test on these denoised outcomes with
decorrelated features) is authored after E3b review. Do NOT fit any mechanism model, construct
designs (E5), open freeze (E6), or read any test outcome (E7). The frozen H2 FAIL stands regardless of
what E3b/E4b find.

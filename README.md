# Git For Reliability Paper

Replication code for the paper "How Many Problems Are Enough, and Which? A
Reliability Ceiling for Individual-Difference Measures of Risky Choice".

The numbered R scripts in `steps/` run in order and write results to `output/`.
`check.R` compares the selected result files stored in `expected/` with the
corresponding outputs.

## What you need

- R 4.6.1 and the package versions recorded in `renv.lock` (installed with
  renv 1.2.3, see below).
- A machine with at least 12 cores. Most long steps start 12 worker processes,
  and step 34 starts 14.
- About 200 MB of disk space for the raw data (135 MB) and the output, plus the
  R package library.
- The raw data in `data/raw/`. See [data/README.md](data/README.md) for where to
  download each dataset, the file names and their md5 checksums. Steps 01, 09
  and 12 download their datasets automatically; the Frey et al. (2017) data
  must be downloaded by hand.
- A git clone of this repository (not a zip download): step 27 checks with
  git that the raw data folder is ignored.

## How to run

From the repository root:

```r
install.packages("renv")   # if renv is not installed
renv::restore()
```

then, in a terminal:

```
Rscript run_all.R
```

`run_all.R` runs each step in a fresh R process and stops at the first failure.
It writes console output to `output/logs/` and runs `check.R` when the selected
range ends at step 39. To split a run, use `Rscript run_all.R 1 23` followed by
`Rscript run_all.R 24 39`.

On Windows without Rtools, `renv::restore()` may fail to build `xfun` 0.59 from
source. In that case run `renv::install("xfun")` and then `renv::restore()`
again; xfun is only used by knitr to print tables.

## Runtime

About 21 hours on 12 cores when nothing else is running (a full run on a busy
machine took 33 hours). The long steps are 29 (about 7.3 hours), 24 (about 6
hours), 30 (about 2 hours), 28 (about 75 minutes), 38 and 31 (about 50 minutes
each), and 05, 22, 34, 15, 16, 11, 07 and 06 (13 to 30 minutes each). All other
steps take a few seconds to a few minutes. Each step's header gives its runtime.

## How to check

```
Rscript check.R
```

compares the CSV and JSON files in `expected/` with the files of the same name
in `output/`. Numeric values must agree within `1e-10 * max(1, abs(expected))`;
infinite values and NaN must match exactly. Numeric values cannot become
missing, missing values cannot become numeric, and TRUE/FALSE flags must
match. CSV column names and row counts must also match. Other text and image
files are not checked. The script prints one line per file and exits with an
error if a checked value differs or an output file is missing. `expected/`
contains the reference files from the original analysis.

Paths written inside output files and log messages (for example
`results/step2/reliability_by_cell_construct.csv`) are relative to `output/`.

## Steps

| Step | Script | What it computes | Paper |
|---|---|---|---|
| 01 | `01_provenance_and_recovery.R` | Data checks, sample description, parameter recovery (simulation) | Methods; sample sizes |
| 02 | `02_firth_estimator_checks.R` | Firth versus plain logistic estimator, negative control (simulation) | Methods (estimator choice) |
| 03 | `03_firth_noise_correlation.R` | Correlation of estimation noise between parameters (simulation) | Methods |
| 04 | `04_bias_measures.R` | Per-participant bias measures and same-trial agreement | Methods; Results |
| 05 | `05_split_half_reliability.R` | Split-half reliability curves, strict problem pool | Results; Figure 1 |
| 06 | `06_cross_half_agreement.R` | Same-trial versus cross-half agreement | Results |
| 07 | `07_fixed_pool_check.R` | Reliability on a fixed problem pool | Results |
| 08 | `08_prior_study_reconciliation.R` | Prior studies on the problem-count axis | Relation to prior findings |
| 09 | `09_reconstruction_validation.R` | Multi-outcome lotteries checked against the official CPC18 code | Methods (problem pools) |
| 10 | `10_expanded_pool_recovery.R` | Expanded problem pool and recovery on it (simulation) | Methods; Results |
| 11 | `11_expanded_pool_reliability.R` | Split-half reliability curves, expanded problem pool | Results; Figure 1 |
| 12 | `12_choices13k_checks.R` | Aggregate agreement with choices13k | Results |
| 13 | `13_reconciliation_figure.R` | Figure 1 and its data | Results; Figure 1 |
| 14 | `14_summary_tables.R` | Trials-needed table and claim-status table | Results; Appendix; Supplement |
| 15 | `15_aggregation_sensitivity.R` | Sensitivity to the aggregation estimator | Supplement |
| 16 | `16_fixed_cohort_reliability.R` | Reliability on a fixed participant cohort | Supplement |
| 17 | `17_problem_assignment_audit.R` | Which problems each participant saw (problem forms) | Methods |
| 18 | `18_variance_ceiling.R` | Variance decomposition and reliability ceiling | Study 2 |
| 19 | `19_variance_ceiling_recovery.R` | Recovery check (simulation) and final ceiling values | Study 2; Abstract |
| 20 | `20_participant_split.R` | Build, freeze and test participant groups | Methods; Study 2 |
| 21 | `21_problem_features.R` | Problem-information dictionary | Methods |
| 22 | `22_problem_subset_reliability.R` | Reliability of equal-length problem subsets (H1) | Study 2 |
| 23 | `23_noise_floor_sensitivity.R` | Sensitivity of the subset noise floor | Study 2 |
| 24 | `24_denoised_subset_reliability.R` | High-precision outcome reliability of the subsets | Study 2 |
| 25 | `25_mechanism_model.R` | Do problem features predict subset reliability (H2) | Study 2 |
| 26 | `26_hierarchical_reestimation.R` | Hierarchical re-estimation, reliability by budget | Ceiling section; Abstract |
| 27 | `27_frey_external_replication.R` | Re-analysis of the Frey et al. (2017) data | Ceiling section; Abstract |
| 28 | `28_item_information_model.R` | Joint item-information model | Ceiling section; Supplement |
| 29 | `29_mechanism_sensitivity.R` | Detectability of the mechanism, attenuation bound (simulation) | Study 2; Supplement |
| 30 | `30_partition_stability.R` | Stability of the subset ranking across partitions | Study 2; Supplement |
| 31 | `31_selection_optimism.R` | Optimism of selecting the best subset | Study 2; Supplement |
| 32 | `32_leave_one_feature_out.R` | Leave-one-feature-out check of the mechanism model | Study 2; Supplement |
| 33 | `33_partition_checks_verify.R` | Verification and assembly of steps 30-32 | Study 2; Supplement |
| 34 | `34_local_dependence_diagnostics.R` | Local dependence and response-style diagnostics | Study 2; Supplement |
| 35 | `35_local_dependence_verify.R` | Verification and assembly of step 34 | Study 2; Supplement |
| 36 | `36_summary_cells.R` | Summary values used in the text | Study 2; Ceiling section; Supplement |
| 37 | `37_hypothesis_ledger.R` | Table of pre-specified hypotheses and outcomes | Hypotheses table |
| 38 | `38_within_form_robustness.R` | Within-form variance decomposition and reliability | Results; Study 2; Supplement |
| 39 | `39_construct_correlation.R` | Correlation between loss avoidance and EV sensitivity | Methods |

"Ceiling section" is the results section "The reliability ceiling reflects
true-score variance, not estimation". Each script starts with a short header
listing what it reads and writes. Recovery and calibration results are labeled
as simulations; those scripts also read the empirical problem features used
to define the simulated designs. Step 19 combines the recovery verdicts with
the empirical ceiling estimates from step 18. The paper's other two figures
plot values produced by steps 08, 13, 18 and 19; their plotting scripts are
not included.

## License

MIT (see `LICENSE`).

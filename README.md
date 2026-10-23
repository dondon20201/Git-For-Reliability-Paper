# Git For Reliability Paper

Code used to replicate the results in the paper "How Many Problems Are Enough, and Which? A Reliability Ceiling for Individual-Difference Measures of Risky Choice".
The R scripts in steps/ are executed sequentially and output results to output/. The script check.R then compares the files containing the expected results from expected/ to the actual output files produced.

## What you need

- R 4.6.1 and the package versions stored in renv.lock (from renv 1.2.3, as described below)
- A 12-core minimum machine. Steps involving long computations initiate 12 worker processes initially, and step 34 initiates 14.
- Approximately 200 MB of disk space for raw data (135 MB) and outputs as well as additional space for R package libraries.
- The raw data in data/raw/. Information on download links, file names, and md5 checksums for each dataset is in data/README.md. Datasets used in steps 01, 09 and 12 are automatically downloaded. The data from Frey et al. (2017) must be downloaded manually.
- A git clone of this repository (not a zip file): step 27 uses git to verify that the raw data directory is ignored.

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

run_all.R runs each step in a new R process, stopping at the first failed step. Console output is saved in output/logs/. If the final step number is 39, then check.R will be run. To run half of a series then continue later, do Rscript run_all.R 1 23 then Rscript run_all.R 24 39.

If you are on Windows and do not have Rtools installed, you may encounter problems with renv::restore() trying to install xfun 0.59 from source. If this happens, do renv::install("xfun") then renv::restore(). xfun is only used by knitr for table output.

## Runtime

Roughly 21 hours on 12 cores in idle conditions (33 hours with a busy system). The long execution times are 29 (~7.3 hours), 24 (~6 hours), 30 (~2 hours), 28 (~75 minutes), 38 and 31 (~50 minutes each), and finally 05, 22, 34, 15, 16, 11, 07 and 06 (13-30 minutes each). Everything else is a few seconds to a few minutes. The time for each step is indicated in the header for that step.

## How to check

```
Rscript check.R
```

The script Rscript check.R compares the CSV and JSON files in expected/ with the corresponding files in output/. For numeric data, the values must be equal within 1e-10 * max(1, abs(expected)). NaN and infinite values must be equal. Missing values must not turn into numeric, numeric must not turn into missing, and TRUE/FALSE must remain the same. The number of rows and the column names in the CSV files must also be equal. It skips PNG images and does not compare text in CSV or JSON files. One line is printed per file that is tested, and if a test value is found to be different or if one of the output files is missing an error is returned. The files in expected/ are the results from the original analysis.

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

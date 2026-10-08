# Step 17: audit of which problems each participant saw (problem forms)
# Reads:   data/raw/All estimation raw data.csv
# Writes:  output/elevation/results/e0/ (participant_problem_incidence.csv, problem and form coverage tables,
#          fixed_form_feasibility.json); log in output/elevation/logs/
# Paper:   Methods (Participant split and held-out groups), through step 20
# Runtime: a few seconds

# E0 Part A: describe participant assignment and overlap in the md5-verified
# CPC18 estimation data. The audit assesses whether the observed forms can
# support a fixed short form, feature-based selection, or a simulation-backed
# design; it does not optimize a design.
#
# Reuses load_cpc18 and is_clean_row from R/cpc18_features.R.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})
source(here("R", "cpc18_features.R"))

res_dir <- here("output", "elevation", "results", "e0")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  df
}
md5_vec <- function(gameids) {
  # md5 of the subject's sorted GameID vector, canonicalised as a comma string.
  digest::digest(paste(sort(as.integer(gameids)), collapse = ","),
                 algo = "md5", serialize = FALSE)
}

# ---- load + verify ---------------------------------------------------------
raw_path <- here("data", "raw", "All estimation raw data.csv")
stopifnot(file.exists(raw_path))
md5_disk <- unname(tools::md5sum(raw_path))
cat("raw md5:", md5_disk, "\n")
stopifnot(identical(md5_disk, "38a26407a36c690be969f9d9d3fe3e5a"))

raw <- load_cpc18(raw_path)
raw$is_clean <- is_clean_row(raw)

n_subjects <- dplyr::n_distinct(raw$SubjID)
n_problems <- dplyr::n_distinct(raw$GameID)
cat(sprintf("N_subjects = %d ; N_problems = %d\n", n_subjects, n_problems))
stopifnot(n_subjects == 686L, n_problems == 210L)

# ---- A.1 participant x problem incidence -----------------------------------
# One (SubjID, GameID) pair per problem the subject saw (25 trials each).
subj_games <- raw |>
  dplyr::distinct(SubjID, GameID, is_clean)

subj_meta <- raw |>
  dplyr::distinct(SubjID, Set, Condition)
stopifnot(nrow(subj_meta) == n_subjects)  # Set/Condition constant within subject

incidence <- subj_games |>
  dplyr::group_by(SubjID) |>
  dplyr::summarise(
    n_problems_total = dplyr::n_distinct(GameID),
    n_problems_clean = dplyr::n_distinct(GameID[is_clean]),
    gameid_set_hash = md5_vec(GameID),
    clean_gameid_set_hash = md5_vec(GameID[is_clean]),
    .groups = "drop"
  ) |>
  dplyr::left_join(subj_meta, by = "SubjID") |>
  dplyr::select(SubjID, Set, Condition, n_problems_total, n_problems_clean,
                gameid_set_hash, clean_gameid_set_hash)
save_csv(incidence, "participant_problem_incidence")

cat(sprintf("problems_total: min=%d median=%g max=%d\n",
            min(incidence$n_problems_total), median(incidence$n_problems_total),
            max(incidence$n_problems_total)))
cat(sprintf("problems_clean: min=%d median=%g max=%d\n",
            min(incidence$n_problems_clean), median(incidence$n_problems_clean),
            max(incidence$n_problems_clean)))

# ---- A.2 cohort structure --------------------------------------------------
# Cluster subjects by identical FULL problem-set and, separately, by identical
# CLEAN problem-set. Cross-tabulate against Set and Condition.
modal <- function(x) {
  tt <- sort(table(x), decreasing = TRUE)
  names(tt)[1]
}
sorted_ids <- function(gameids) paste(sort(unique(as.integer(gameids))), collapse = " ")

# Full-set clusters
full_clusters <- incidence |>
  dplyr::group_by(gameid_set_hash) |>
  dplyr::summarise(n_subjects = dplyr::n(),
                   n_problems_total = dplyr::first(n_problems_total),
                   modal_Set = modal(Set), modal_Condition = modal(Condition),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n_subjects))

# Clean-set clusters (the reported cohort_summary)
clean_ids_by_subj <- subj_games |>
  dplyr::filter(is_clean) |>
  dplyr::group_by(SubjID) |>
  dplyr::summarise(clean_hash = md5_vec(GameID),
                   clean_ids = sorted_ids(GameID),
                   n_clean = dplyr::n_distinct(GameID), .groups = "drop") |>
  dplyr::left_join(subj_meta, by = "SubjID")

cohort_summary <- clean_ids_by_subj |>
  dplyr::group_by(set_hash = clean_hash) |>
  dplyr::summarise(n_subjects = dplyr::n(),
                   n_problems_clean = dplyr::first(n_clean),
                   modal_Set = modal(Set), modal_Condition = modal(Condition),
                   sorted_gameids = dplyr::first(clean_ids),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n_subjects))
save_csv(cohort_summary, "cohort_summary")

# Cross-tabs: how clean-set clusters map onto Set and Condition.
clean_cluster_id <- clean_ids_by_subj |>
  dplyr::mutate(cluster = match(clean_hash, cohort_summary$set_hash))
xtab_set <- clean_cluster_id |>
  dplyr::count(cluster, Set) |>
  tidyr::pivot_wider(names_from = Set, values_from = n, values_fill = 0,
                     names_prefix = "Set_")
xtab_cond <- clean_cluster_id |>
  dplyr::count(cluster, Condition) |>
  tidyr::pivot_wider(names_from = Condition, values_from = n, values_fill = 0,
                     names_prefix = "Cond_")
cohort_crosstab <- cohort_summary |>
  dplyr::mutate(cluster = dplyr::row_number()) |>
  dplyr::select(cluster, set_hash, n_subjects, n_problems_clean,
                modal_Set, modal_Condition) |>
  dplyr::left_join(xtab_set, by = "cluster") |>
  dplyr::left_join(xtab_cond, by = "cluster")
save_csv(cohort_crosstab, "cohort_crosstab")

cat(sprintf("distinct FULL problem-sets: %d\n", nrow(full_clusters)))
cat(sprintf("distinct CLEAN problem-sets: %d\n", nrow(cohort_summary)))
cat(sprintf("distinct Set values: %d ; distinct Condition values: %d\n",
            dplyr::n_distinct(subj_meta$Set), dplyr::n_distinct(subj_meta$Condition)))

# ---- A.3 problem coverage --------------------------------------------------
problem_coverage <- subj_games |>
  dplyr::group_by(GameID) |>
  dplyr::summarise(is_clean = any(is_clean),
                   n_subjects_total = dplyr::n_distinct(SubjID),
                   .groups = "drop") |>
  dplyr::arrange(GameID)
save_csv(problem_coverage, "problem_coverage")
cat(sprintf("clean problems: %d ; coverage/problem: min=%d median=%g max=%d\n",
            sum(problem_coverage$is_clean),
            min(problem_coverage$n_subjects_total),
            median(problem_coverage$n_subjects_total),
            max(problem_coverage$n_subjects_total)))

# ---- A.4 fixed-form feasibility --------------------------------------------
# (1) Largest cohort sharing an IDENTICAL clean-problem subset.
top_clean <- cohort_summary |> dplyr::slice(1)
S_max <- top_clean$n_subjects
C_max <- top_clean$n_problems_clean

# (2) Common core: number of CLEAN problems seen by >= T subjects.
clean_cov <- problem_coverage |> dplyr::filter(is_clean)
Ts <- c(686, 500, 300, 200, 150, 100)
common_core <- vapply(Ts, function(t) sum(clean_cov$n_subjects_total >= t),
                      integer(1))
names(common_core) <- as.character(Ts)

# (3) Gate boolean (literal reading): some IDENTICAL clean subset of size >= 8
#     shared by >= 150 subjects, i.e. a clean-set cohort with both properties.
gate_identical_subset <- any(cohort_summary$n_subjects >= 150 &
                               cohort_summary$n_problems_clean >= 8)
# Secondary common-core reading (weaker; problems need not co-occur in one form).
gate_commoncore_150_ge8 <- common_core[["150"]] >= 8

feas <- list(
  n_subjects = n_subjects,
  n_problems = n_problems,
  n_clean_problems = sum(problem_coverage$is_clean),
  n_distinct_full_sets = nrow(full_clusters),
  n_distinct_clean_sets = nrow(cohort_summary),
  largest_identical_clean_cohort = list(
    S_max = S_max, C_max = C_max, set_hash = top_clean$set_hash),
  common_core_clean_by_T = as.list(common_core),
  gate_identical_subset_ge8_ge150 = unname(gate_identical_subset),
  gate_commoncore_150_ge8 = unname(gate_commoncore_150_ge8)
)
jsonlite::write_json(feas, file.path(res_dir, "fixed_form_feasibility.json"),
                     auto_unbox = TRUE, pretty = TRUE)
cat("S_max =", S_max, " C_max =", C_max, "\n")
cat("common_core (>=T subjects share a clean problem):\n")
print(common_core)
cat("gate identical-subset (>=8 clean shared by >=150):", gate_identical_subset, "\n")
cat("gate common-core (>=8 clean each seen by >=150):", gate_commoncore_150_ge8, "\n")

# ---- Per-cohort feature coverage (SELECTION-RULE path evidence) ------------
# For the modal (largest) clean cohort, report min/max/IQR of the clean
# difference-features spanned WITHIN that cohort's problem set.
features <- clean_problem_features(raw)
top_ids <- as.integer(strsplit(top_clean$sorted_gameids, " ")[[1]])
top_feat <- features |> dplyr::filter(GameID %in% top_ids)
feat_cov <- tibble::tibble(
  feature = c("d_ev", "d_sd", "d_eloss"),
  min = c(min(top_feat$d_ev), min(top_feat$d_sd), min(top_feat$d_eloss)),
  q25 = c(quantile(top_feat$d_ev, .25), quantile(top_feat$d_sd, .25),
          quantile(top_feat$d_eloss, .25)),
  median = c(median(top_feat$d_ev), median(top_feat$d_sd),
             median(top_feat$d_eloss)),
  q75 = c(quantile(top_feat$d_ev, .75), quantile(top_feat$d_sd, .75),
          quantile(top_feat$d_eloss, .75)),
  max = c(max(top_feat$d_ev), max(top_feat$d_sd), max(top_feat$d_eloss)),
  n_problems = nrow(top_feat)
) |>
  dplyr::mutate(IQR = q75 - q25)
save_csv(feat_cov, "modal_cohort_feature_coverage")
# Full clean-pool feature ranges for comparison.
full_feat_cov <- tibble::tibble(
  feature = c("d_ev", "d_sd", "d_eloss"),
  pool_min = c(min(features$d_ev), min(features$d_sd), min(features$d_eloss)),
  pool_max = c(max(features$d_ev), max(features$d_sd), max(features$d_eloss)),
  pool_IQR = c(IQR(features$d_ev), IQR(features$d_sd), IQR(features$d_eloss)),
  n_clean_problems = nrow(features)
)
save_csv(full_feat_cov, "pool_feature_coverage")

cat("\n== modal cohort feature coverage ==\n"); print(feat_cov)
cat("\n== full clean-pool feature coverage ==\n"); print(full_feat_cov)

# ---- Assignment-key design + per-form informative coverage -----------------
# Discovered key: Set (7 forms). Condition (feedback regime) is crossed within
# Sets 1-5 on the same problem set; Sets 6-7 are ByProb-only.
cohort_design <- subj_meta |>
  dplyr::count(Set, Condition) |>
  tidyr::pivot_wider(names_from = Condition, values_from = n, values_fill = 0)
save_csv(cohort_design, "cohort_design_set_by_condition")

per_form <- raw |>
  dplyr::filter(is_clean) |>
  dplyr::distinct(Set, GameID) |>
  dplyr::left_join(features, by = "GameID") |>
  dplyr::group_by(Set) |>
  dplyr::summarise(
    n_clean = dplyr::n(),
    n_loss_informative = sum(loss_informative),
    n_dev_gt1 = sum(abs(d_ev) > 1),
    n_dsd_neg = sum(d_sd < 0),
    n_deloss_nonzero = sum(d_eloss != 0),
    d_ev_min = min(d_ev), d_ev_max = max(d_ev),
    d_sd_min = min(d_sd), d_sd_max = max(d_sd),
    d_eloss_min = min(d_eloss), d_eloss_max = max(d_eloss),
    .groups = "drop")
save_csv(per_form, "per_form_feature_coverage")
cat("\n== per-form informative coverage ==\n"); print(per_form)

cat("\nDONE Part A\n")

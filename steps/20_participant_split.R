# Step 20: participant split into build, freeze and test groups
# Reads:   output/elevation/results/e0/participant_problem_incidence.csv (step 17)
# Writes:  output/elevation/results/e1/ (participant_split.csv, split_summary.csv, baseline_grid.csv, configuration)
# Paper:   Methods (Participant split and held-out groups); Study 2; Limitations
# Runtime: a few seconds

# E1 Part B: three-way participant split and fixed baseline grid.
# This pre-specification step reads only the participant-to-Set incidence table
# from step 17 (SubjID, Set), not trial choices. It performs no optimization,
# mechanism modeling, recovery, or prediction.
#
# Within each form (Set 1-7), assign participants to build / freeze / test
# in target proportions 0.50 / 0.20 / 0.30. Largest-remainder rounding gives
# integer counts and representation of every form in each group.
# Pre-declared seed: 20260718.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

res_dir <- here("output", "elevation", "results", "e1")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  df
}

# ---- inputs (step 17 incidence only) ---------------------------------------
inc_path <- here("output", "elevation", "results", "e0", "participant_problem_incidence.csv")
stopifnot(file.exists(inc_path))
inc <- readr::read_csv(inc_path, show_col_types = FALSE)

# One row per subject; each subject belongs to exactly one Set (form).
subs <- inc |>
  dplyr::distinct(SubjID, Set) |>
  dplyr::arrange(SubjID)
stopifnot(nrow(subs) == 686L)                       # E0: 686 subjects
stopifnot(!anyDuplicated(subs$SubjID))              # SubjID unique
stopifnot(all(!is.na(subs$Set)))                    # every subject has a Set

# ---- pre-declared split parameters -----------------------------------------
SEED       <- 20260718L
split_lvls <- c("build", "freeze", "test")
target_p   <- c(build = 0.50, freeze = 0.20, test = 0.30)
stopifnot(abs(sum(target_p) - 1) < 1e-12)

# Largest-remainder allocation of n units to the 3 splits at target_p.
# Ties in the fractional remainder are broken in split_lvls order (build,
# then freeze, then test) via a stable descending order on the remainder.
alloc_counts <- function(n, p) {
  raw   <- n * p
  base  <- floor(raw)
  rem   <- n - sum(base)
  if (rem > 0) {
    frac <- raw - base
    ord  <- order(-frac)              # stable: ties keep split_lvls order
    base[ord[seq_len(rem)]] <- base[ord[seq_len(rem)]] + 1L
  }
  stats::setNames(as.integer(base), names(p))
}

# ---- seeded within-Set assignment ------------------------------------------
set.seed(SEED)
assign_one_set <- function(ids) {
  n  <- length(ids)
  cn <- alloc_counts(n, target_p)
  stopifnot(sum(cn) == n)
  labels <- rep(split_lvls, times = cn[split_lvls])
  # Randomly pair the fixed label composition to the (sorted) subject ids.
  tibble::tibble(SubjID = ids, split = sample(labels, size = n, replace = FALSE))
}

split_tbl <- subs |>
  dplyr::group_by(Set) |>
  dplyr::group_modify(~ assign_one_set(sort(.x$SubjID))) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    split = factor(split, levels = split_lvls),
    held_out_until_E7 = split == "test"
  ) |>
  dplyr::arrange(SubjID) |>
  dplyr::select(SubjID, Set, split, held_out_until_E7)

# ---- summaries --------------------------------------------------------------
per_set <- split_tbl |>
  dplyr::count(Set, split, name = "n") |>
  tidyr::complete(Set, split = factor(split_lvls, levels = split_lvls),
                  fill = list(n = 0L)) |>
  dplyr::group_by(Set) |>
  dplyr::mutate(set_n = sum(n), prop = n / set_n) |>
  dplyr::ungroup()

pooled <- split_tbl |>
  dplyr::count(split, name = "n") |>
  dplyr::mutate(Set = "pooled", set_n = sum(n), prop = n / set_n) |>
  dplyr::select(Set, split, n, set_n, prop)

split_summary <- dplyr::bind_rows(
  per_set |> dplyr::mutate(Set = as.character(Set)),
  pooled
) |>
  dplyr::arrange(match(Set, c(as.character(sort(unique(subs$Set))), "pooled")),
                 match(split, split_lvls))

# ---- fixed baseline grid (grid only here; no designs are constructed) ------
baseline_grid <- tidyr::expand_grid(
  B = c(6L, 8L, 10L, 12L),
  design_type = c("optimized", "random", "stratified_random", "random_double")
) |>
  dplyr::mutate(
    # random_double is the 2K efficiency comparator; it samples 2*B.
    n_problems = dplyr::if_else(design_type == "random_double", 2L * B, B)
  ) |>
  dplyr::arrange(B, match(design_type,
    c("optimized", "random", "stratified_random", "random_double")))

# ---- ASSERTIONS (pre-declared PASS gates) ----------------------------------
# partition of all 686 subjects, each exactly once
stopifnot(nrow(split_tbl) == 686L)
stopifnot(!anyDuplicated(split_tbl$SubjID))
stopifnot(setequal(split_tbl$SubjID, subs$SubjID))
# every Set appears in all three splits (cohort-respecting)
cov <- split_tbl |> dplyr::count(Set, split) |> dplyr::filter(n > 0)
stopifnot(nrow(cov) == length(unique(subs$Set)) * length(split_lvls))
# pooled test proportion ~ 0.30, in [0.28, 0.32]
test_prop <- mean(split_tbl$split == "test")
stopifnot(test_prop >= 0.28, test_prop <= 0.32)
# test rows flagged for E7 hold-out; build/freeze not flagged
stopifnot(all(split_tbl$held_out_until_E7 == (split_tbl$split == "test")))
# no outcome / choice column computed anywhere
stopifnot(identical(names(split_tbl),
                    c("SubjID", "Set", "split", "held_out_until_E7")))

# ---- write outputs ----------------------------------------------------------
save_csv(split_tbl |> dplyr::mutate(split = as.character(split)),
         "participant_split")
save_csv(split_summary |> dplyr::mutate(split = as.character(split)),
         "split_summary")
save_csv(baseline_grid, "baseline_grid")

cfg <- list(
  seed = SEED,
  target_proportions = as.list(target_p),
  rounding = "largest-remainder within each Set (ties in split_lvls order)",
  stratify_by = "Set (form 1-7); cohort-respecting",
  n_subjects = nrow(split_tbl),
  pooled_counts = as.list(table(split_tbl$split)[split_lvls]),
  pooled_test_proportion = round(test_prop, 6),
  held_out_group = "test (reserved for E7 confirmatory; never inspected before E7)",
  baseline_budgets_B = c(6, 8, 10, 12),
  baseline_design_types = c("optimized", "random", "stratified_random", "random_double"),
  data_input = "elevation/results/e0/participant_problem_incidence.csv (SubjID, Set only)",
  no_outcome_computed = TRUE
)
jsonlite::write_json(cfg, file.path(res_dir, "e1_split_config.json"),
                     pretty = TRUE, auto_unbox = TRUE)

cat("E1 split written. pooled counts:\n")
print(table(split_tbl$split)[split_lvls])
cat(sprintf("pooled test proportion = %.4f\n", test_prop))
cat("PASS: partition of 686; every Set in all 3 splits; no outcome computed.\n")

# Step 07: fixed problem pool check of the reliability curve shape
# Reads:   data/raw/All estimation raw data.csv
# Writes:  output/results/step2/fixed_pool_check.csv, output/results/step2/set_by_clean_count.csv
# Paper:   Results (Reliability generally increases with distinct-problem count)
# Runtime: about 15 minutes with 6 parallel workers

# Sanity-anchor investigation of a non-monotone step 05 reliability curve.
#
# Finding to explain: ev_sens_agn's axis-(i) curve dips at k = 5 in BOTH
# regimes, exactly where the eligible pool changes (678 -> 556 subjects), so
# budget is confounded with subject composition (eligibility tracks Set
# membership: each Set has a fixed problem list, hence a fixed clean count).
# Test: recompute the experience problems-axis on the FIXED pool of subjects
# eligible at every k (>= 14 clean problems), k in 2..7. If the dip disappears
# on the fixed pool, composition — not budget — explains it.
# REAL DATA.

suppressMessages({library(tidyverse); library(here); library(parallel)})
source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))

N_SPLITS <- 200
N_BOOT <- 2000
SEED_BASE <- 187000  # cells 187001..187006; bootstrap seeds 187101..

raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))
features <- add_scaled_features(clean_problem_features(raw))
trials <- raw |>
  filter(is_clean_row(raw)) |>
  select(SubjID, GameID, Trial, B) |>
  left_join(features |>
              select(GameID, z_ev, z_sd, z_eloss, d_ev, loss_informative,
                     has_loss_b, ev_a, ev_b), by = "GameID")
baselines <- fit_loss_baselines(trials)
subj_data <- build_subject_data(trials, baselines)

# Set-composition table for the eligibility pools (Set = leading SubjID digit,
# per the data dictionary; verified against the Set column)
subj_info <- raw |> distinct(SubjID, Set)
clean_counts <- tibble(SubjID = as.numeric(names(subj_data)),
                       n_clean = vapply(subj_data, function(s) nrow(s$feat), integer(1))) |>
  left_join(subj_info, by = "SubjID")
set_table <- clean_counts |>
  count(Set, n_clean) |>
  pivot_wider(names_from = n_clean, values_from = n, values_fill = 0,
              names_prefix = "clean_")
readr::write_csv(set_table, here("output", "results", "step2", "set_by_clean_count.csv"))
print(set_table)

fixed_pool <- names(subj_data)[vapply(subj_data, function(s) nrow(s$feat), integer(1)) >= 14]
cat("fixed pool (>= 14 clean problems):", length(fixed_pool), "subjects\n")
sd_fixed <- subj_data[fixed_pool]

cl <- makeCluster(6)
invisible(clusterCall(cl, function(lp, root) {
  .libPaths(lp); source(file.path(root, "R", "step2_reliability.R")); NULL
}, .libPaths(), here::here()))
clusterExport(cl, c("sd_fixed", "N_SPLITS", "SEED_BASE"))
cells <- parLapplyLB(cl, 2:7, function(k) {
  run_cell(list(cell_id = paste0("fixed_k", k), axis = "i_exp_fixed",
                regime = "exp", k = k, r = 20, lambda = NA),
           sd_fixed, N_SPLITS, seed = SEED_BASE + k)
})
stopCluster(cl)

out <- purrr::map_dfr(cells, function(cell) {
  purrr::map_dfr(c("ev_sens_agn", "loss_avoid_adj", "b_ev"), function(con) {
    E1 <- cell$est[, , 1, con]; E2 <- cell$est[, , 2, con]
    res <- if (!any(is.finite(E1) & is.finite(E2))) {
      tibble(r_raw = NA, ci_lo_raw = NA, ci_hi_raw = NA, r_sb = NA,
             ci_lo_sb = NA, ci_hi_sb = NA, n_subjects = nrow(E1),
             n_splits_used = 0L, n_splits = ncol(E1), finite_frac = 0,
             n_boot = N_BOOT)
    } else {
      reliability_point(E1, E2, N_BOOT, seed = SEED_BASE + 100 + cell$spec$k)
    }
    bind_cols(tibble(k = cell$spec$k, construct = con,
                     pool = "fixed (>=14 clean problems)"), res)
  })
})
readr::write_csv(out, here("output", "results", "step2", "fixed_pool_check.csv"))
print(out |> select(k, construct, r_sb, ci_lo_sb, ci_hi_sb, n_subjects), n = 24)
cat("ADDENDUM_DONE\n")

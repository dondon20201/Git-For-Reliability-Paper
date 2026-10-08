# Step 15: sensitivity of reliability to the aggregation estimator
# Reads:   data/raw/All estimation raw data.csv; output/results/step2/, step4/, step5/ (steps 05, 09-11, 14)
# Writes:  output/results/step2/ (aggregation_sensitivity.csv, aggregation_summary.csv, aggregation_split_r_headline.csv)
# Paper:   Supplement (Aggregation-estimator sensitivity)
# Runtime: about 24 minutes with 12 parallel workers

# Compare aggregation rules for the same split-half correlations.
# Steps 05 and 11 use the arithmetic mean of up to 200 per-split Spearman
# correlations, dropping and counting splits with fewer than 10 complete
# pairs, then apply Spearman-Brown (R/step2_reliability.R::reliability_point).
#
# Rerun run_cell / split_correlations / sb_correct with the saved seeds for
# every reported cell in both pools. Assert that the original point estimates
# are reproduced within 1e-9, then compare median and Fisher-z aggregation
# of those same correlations. Check whether these alternatives change the
# threshold crossings in the step 14 trials-needed table.
#
# Only the three aggregation_* CSVs are written; earlier tables and figures
# are retained.

library(tidyverse)
library(here)
library(parallel)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "fetch_official_cpc18.R"))
source(here("R", "lotshape_reconstruction.R"))
source(here("R", "step2_reliability.R"))

res2_dir <- here("output", "results", "step2")
res4_dir <- here("output", "results", "step4")

N_SPLITS <- 200   # pre-declared

# ---- saved reference files --------------------------------------------------
committed_strict <- read_csv(file.path(res2_dir, "reliability_by_cell_construct.csv"),
                             show_col_types = FALSE)
committed_ext <- read_csv(file.path(res4_dir, "reliability_extended.csv"),
                          show_col_types = FALSE)
seeds_strict <- read_csv(file.path(res2_dir, "seeds_cells.csv"),
                         show_col_types = FALSE)
seeds_ext <- read_csv(file.path(res4_dir, "seeds_cells_extended.csv"),
                      show_col_types = FALSE)
trials_needed_committed <- read_csv(here("output", "results", "step5", "trials_needed_combined.csv"),
                                    show_col_types = FALSE)

# ---- raw data (md5-verified fetch, as in every step) -------------------------
fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))
raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))

# ---- STRICT pool subject data (exact step 05 path) ---------------------------
features <- add_scaled_features(clean_problem_features(raw))
trials_strict <- raw |>
  filter(is_clean_row(raw)) |>
  select(SubjID, GameID, Trial, B) |>
  left_join(features |>
              select(GameID, z_ev, z_sd, z_eloss, d_ev, loss_informative,
                     has_loss_b, ev_a, ev_b), by = "GameID")
baselines_strict <- fit_loss_baselines(trials_strict)
subj_data_strict <- build_subject_data(trials_strict, baselines_strict)
cat(sprintf("strict pool: %d subjects\n", length(subj_data_strict)))

# ---- EXPANDED pool subject data (exact step 11 path, gates included) ---------
gate <- read_csv(file.path(res4_dir, "reconstruction_validation.csv"),
                 show_col_types = FALSE)
stopifnot(all(gate$pass))
verdict_a4 <- read_csv(file.path(res4_dir, "expanded_recovery_verdict.csv"),
                       show_col_types = FALSE)
stopifnot(all(verdict_a4$pass))

param_cols <- c("Ha", "pHa", "La", "LotShapeA", "LotNumA",
                "Hb", "pHb", "Lb", "LotShapeB", "LotNumB", "Amb", "Corr")
problems_exp <- raw |>
  distinct(GameID, across(all_of(param_cols))) |>
  arrange(GameID) |>
  filter(Amb == 0)
feat_exp <- expanded_problem_features(problems_exp,
                                      reconstruct_all_problems(problems_exp))
trials_ext <- raw |>
  semi_join(problems_exp, by = "GameID") |>
  select(SubjID, GameID, Trial, B) |>
  left_join(feat_exp |>
              select(GameID, z_ev, z_sd, z_eloss, d_ev, loss_informative,
                     has_loss_b, ev_a, ev_b), by = "GameID")
baselines_ext <- fit_loss_baselines(trials_ext)
subj_data_ext <- build_subject_data(trials_ext, baselines_ext)
cat(sprintf("expanded pool: %d subjects\n", length(subj_data_ext)))

# ---- cell specs: SAME construction order as steps 05 and 11, with seeds
# bound from their saved seeds CSVs by cell_id (never re-derived) --------------
strict_specs <- c(
  lapply(2:7, function(k) list(axis = "i_exp", regime = "exp", k = k, r = 20, lambda = NA)),
  lapply(c(5, 10), function(r) list(axis = "ii_reps", regime = "exp", k = 5, r = r, lambda = NA)),
  lapply(1:3, function(l) list(axis = "iii_lambda", regime = "exp", k = 5, r = 20, lambda = l)),
  lapply(2:7, function(k) list(axis = "i_desc", regime = "desc", k = k, r = 5, lambda = NA))
)
strict_specs <- lapply(strict_specs, function(s) {
  s$cell_id <- sprintf("%s_k%d_r%d%s", s$axis, s$k, s$r,
                       ifelse(is.na(s$lambda), "", paste0("_l", s$lambda)))
  hit <- seeds_strict$seed[seeds_strict$cell_id == s$cell_id]
  stopifnot(length(hit) == 1L)
  s$seed <- hit
  s$pool <- "strict"
  s
})
stopifnot(identical(vapply(strict_specs, `[[`, character(1), "cell_id"),
                    seeds_strict$cell_id))

realized <- read_csv(file.path(res4_dir, "realized_grid.csv"),
                     show_col_types = FALSE)
ext_specs <- lapply(realized$k, function(k) {
  s <- list(axis = "i_exp_ext", regime = "exp", k = k, r = 20, lambda = NA,
            cell_id = sprintf("ext_exp_k%d_r20", k))
  hit <- seeds_ext$seed[seeds_ext$cell_id == s$cell_id]
  stopifnot(length(hit) == 1L)
  s$seed <- hit
  s$pool <- "expanded"
  s
})
stopifnot(identical(vapply(ext_specs, `[[`, character(1), "cell_id"),
                    seeds_ext$cell_id))

all_specs <- c(strict_specs, ext_specs)

# ---- re-run the original estimation (same engine, same seeds; parallel is
# safe because run_cell seeds internally) --------------------------------------
n_workers <- max(1, min(12, detectCores() - 2))
cat("workers:", n_workers, "\n")
cl <- makeCluster(n_workers)
invisible(clusterCall(cl, function(lp, root) {
  .libPaths(lp)
  source(file.path(root, "R", "step2_reliability.R"))
  NULL
}, .libPaths(), here::here()))
clusterExport(cl, c("subj_data_strict", "subj_data_ext", "N_SPLITS"),
              envir = environment())

t0 <- Sys.time()
cells <- parLapplyLB(cl, all_specs, function(spec) {
  sd_pool <- if (spec$pool == "strict") subj_data_strict else subj_data_ext
  out <- run_cell(spec, sd_pool, N_SPLITS, seed = spec$seed)
  out$pool <- spec$pool
  out
})
stopCluster(cl)
names(cells) <- vapply(cells, function(x) x$spec$cell_id, character(1))
cat("estimation wall time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---- per-cell-construct aggregation (mean reproduces the saved point) --------
# Reported constructs mirror the saved files: strict iii_lambda cells report
# only the two loss constructs; every other cell reports all five.
fisher_mean <- function(rs) {
  tanh(mean(atanh(pmin(pmax(rs, -0.999), 0.999)), na.rm = TRUE))
}

rows <- list()
split_rows <- list()
headline <- tibble(
  pool = c("strict", "strict", "expanded", "expanded"),
  cell_id = c("i_exp_k6_r20", "i_exp_k7_r20", "ext_exp_k12_r20", "ext_exp_k13_r20"),
  construct = c("loss_avoid_adj", "loss_avoid_adj", "loss_avoid_adj", "b_ev")
)

for (cell in cells) {
  spec <- cell$spec
  constructs <- if (spec$axis == "iii_lambda") {
    c("b_loss", "loss_avoid_adj")
  } else step2_constructs
  for (con in constructs) {
    E1 <- cell$est[, , 1, con]
    E2 <- cell$est[, , 2, con]
    rs <- split_correlations(E1, E2, min_pairs = 10)
    n_used <- sum(!is.na(rs))
    if (n_used == 0) {
      agg <- c(raw_mean = NA_real_, med = NA_real_, fis = NA_real_)
    } else {
      agg <- c(raw_mean = mean(rs, na.rm = TRUE),
               med = median(rs, na.rm = TRUE),
               fis = fisher_mean(rs))
    }
    rows[[length(rows) + 1]] <- tibble(
      pool = cell$pool, cell_id = spec$cell_id, axis = spec$axis,
      regime = spec$regime, k = spec$k, r = spec$r, lambda = spec$lambda,
      construct = con, n_splits_used = n_used,
      r_sb_mean = sb_correct(agg[["raw_mean"]]),
      r_sb_median = sb_correct(agg[["med"]]),
      r_sb_fisher = sb_correct(agg[["fis"]])
    )
    if (any(headline$cell_id == spec$cell_id & headline$construct == con)) {
      split_rows[[length(split_rows) + 1]] <- tibble(
        pool = cell$pool, construct = con, k = spec$k,
        split_index = seq_along(rs), r_split = rs
      )
    }
  }
}
sensitivity <- bind_rows(rows) |>
  mutate(d_median = r_sb_median - r_sb_mean,
         d_fisher = r_sb_fisher - r_sb_mean)
split_headline <- bind_rows(split_rows)
stopifnot(nrow(split_headline) == 4 * N_SPLITS)

# ---- REPRODUCTION GATE (pre-declared): sb(raw_mean) must equal the saved
# r_sb to < 1e-9 for EVERY reported cell-construct; NA must match NA. Stop on
# any failure. -----------------------------------------------------------------
committed <- bind_rows(
  committed_strict |> mutate(pool = "strict"),
  committed_ext |> mutate(pool = "expanded")
) |> select(pool, cell_id, construct, r_sb_committed = r_sb)

check <- sensitivity |>
  left_join(committed, by = c("pool", "cell_id", "construct"))
stopifnot(nrow(check) == nrow(sensitivity),
          !any(is.na(check$r_sb_committed) & is.finite(check$r_sb_mean)),
          !any(is.finite(check$r_sb_committed) & !is.finite(check$r_sb_mean)))
both <- is.finite(check$r_sb_committed) & is.finite(check$r_sb_mean)
max_repro_dev <- max(abs(check$r_sb_mean[both] - check$r_sb_committed[both]))
cat(sprintf("reproduction gate: %d cell-constructs (%d finite), max |dev| = %.3g\n",
            nrow(check), sum(both), max_repro_dev))
if (!(max_repro_dev < 1e-9)) {
  stop(sprintf("REPRODUCTION GATE FAILED: max |sb(raw_mean) - committed r_sb| = %.3g >= 1e-9; the re-run is not faithful — STOP",
               max_repro_dev))
}

# ---- CROSSING-INVARIANCE check (pre-declared): first-crossing budget per
# axis x construct x threshold under each aggregator, compared to the
# step 14 trials-needed table; same rule (smallest budget with finite
# r_sb >= threshold, no extrapolation) as the original construction. -----------
thresholds <- c(0.5, 0.7, 0.8)
axis_budget <- list(
  i_exp = "k", ii_reps = "r", iii_lambda = "lambda", i_desc = "k",
  i_exp_ext = "k"
)
first_cross <- function(df, value_col) {
  # df: rows of one axis (ii_reps already augmented), with a budget column
  df |>
    group_by(construct) |>
    group_modify(function(g, key) {
      map_dfr(thresholds, function(th) {
        v <- g[[value_col]]
        hit <- g |> filter(is.finite(v), v >= th) |> arrange(budget) |> slice(1)
        tibble(threshold = th,
               budget = if (nrow(hit) == 0) NA_real_ else hit$budget)
      })
    }) |>
    ungroup()
}
crossings <- imap_dfr(axis_budget, function(bvar, ax) {
  m <- sensitivity |> filter(axis == ax |
                               (ax == "ii_reps" & axis == "i_exp" & k == 5))
  m$budget <- m[[bvar]]
  map_dfr(c(mean = "r_sb_mean", median = "r_sb_median", fisher = "r_sb_fisher"),
          function(col) {
            first_cross(m, col) |> mutate(aggregator = col)
          }, .id = "agg_name") |>
    mutate(axis = ax, .before = 1)
})
cross_wide <- crossings |>
  select(axis, construct, threshold, agg_name, budget) |>
  pivot_wider(names_from = agg_name, values_from = budget, names_prefix = "budget_")

committed_cross <- trials_needed_committed |>
  select(axis, construct, threshold, budget_committed = budget)
cross_cmp <- committed_cross |>
  left_join(cross_wide, by = c("axis", "construct", "threshold"))
stopifnot(nrow(cross_cmp) == nrow(committed_cross))
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
# sanity: the mean aggregator must reproduce every saved crossing
stopifnot(all(same(cross_cmp$budget_committed, cross_cmp$budget_mean)))
changed <- !(same(cross_cmp$budget_committed, cross_cmp$budget_median) &
               same(cross_cmp$budget_committed, cross_cmp$budget_fisher))
n_crossings_changed <- sum(changed)
cat(sprintf("crossing invariance: %d/%d crossings changed under median/Fisher-z\n",
            n_crossings_changed, nrow(cross_cmp)))
if (n_crossings_changed > 0) {
  print(cross_cmp[changed, ], n = 100)
}

# ---- outputs -----------------------------------------------------------------
out_sensitivity <- sensitivity |>
  select(pool, axis, regime, k, r, lambda, construct, n_splits_used,
         r_sb_mean, r_sb_median, r_sb_fisher, d_median, d_fisher)
write_csv(out_sensitivity, file.path(res2_dir, "aggregation_sensitivity.csv"))

fin <- is.finite(out_sensitivity$d_median) & is.finite(out_sensitivity$d_fisher)
max_abs_delta_sb_median <- max(abs(out_sensitivity$d_median[fin]))
max_abs_delta_sb_fisher <- max(abs(out_sensitivity$d_fisher[fin]))
max_abs_delta_sb <- max(max_abs_delta_sb_median, max_abs_delta_sb_fisher)
summary_df <- tibble(
  key = c("max_abs_delta_sb", "max_abs_delta_sb_fisher",
          "max_abs_delta_sb_median", "n_crossings_changed", "n_cells"),
  value = c(max_abs_delta_sb, max_abs_delta_sb_fisher,
            max_abs_delta_sb_median, n_crossings_changed,
            nrow(out_sensitivity))
)
write_csv(summary_df, file.path(res2_dir, "aggregation_summary.csv"))

write_csv(split_headline, file.path(res2_dir, "aggregation_split_r_headline.csv"))

cat("\n---- GATE G1 (pre-declared) ----\n")
cat(sprintf("n_crossings_changed = %d (bound: == 0)\n", n_crossings_changed))
cat(sprintf("max_abs_delta_sb    = %.6g (bound: < 0.02)\n", max_abs_delta_sb))
cat(sprintf("G1 %s\n",
            if (n_crossings_changed == 0 && max_abs_delta_sb < 0.02) "PASSES" else "FAILS"))
cat("\nDone: aggregation_sensitivity.csv, aggregation_summary.csv, aggregation_split_r_headline.csv\n")

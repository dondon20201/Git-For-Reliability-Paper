# Step 01: data checks, sample description and parameter recovery (simulation)
# Reads:   data/raw/All estimation raw data.csv (downloaded and md5-checked if absent),
#          the dictionary PDF; R/fetch_cpc18.R, R/cpc18_features.R, R/simulate_recovery.R
# Writes:  output/results/step0/ (structure_checks.csv, problems_per_subject_distribution.csv,
#          trials_per_subject_summary.csv, recovery_verdict.csv, sanity checks, seeds),
#          output/figures/step0_*.png
# Paper:   Methods (Data; Validation gates); sample sizes in the Abstract and Introduction
# Runtime: about 1 minute, single process

# ---- setup ----
library(tidyverse)
library(here)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "simulate_recovery.R"))

res_dir <- here("output", "results", "step0")
fig_dir <- here("output", "figures")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
save_result <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  df
}

# All RNG seeds used in this step (recorded to output/results/step0/seeds.json).
seeds <- list(
  sim_validation = 181001,
  recovery       = 181002,
  loss_sweep     = 181003,
  bootstrap      = 181004
)
jsonlite::write_json(seeds, file.path(res_dir, "seeds.json"), auto_unbox = TRUE)

# Simulation design constants
N_SIM       <- 500          # simulated subjects per recovery cell
N_BOOT      <- 2000         # bootstrap resamples per CI
K_GRID      <- c(5, 10, 14) # clean problems per simulated subject
REPS_GRID   <- c(5, 25)     # repetitions per problem (5 = description-only)
L_VALUES    <- 1:7          # loss-informative problems in the targeted sweep
K_SWEEP     <- 14           # total problems in the loss sweep
VALID_REPS  <- 10000        # repetitions per problem in simulator validation

# ---- fetch ----
fetch_res <- save_result(fetch_cpc18(here("data", "raw")), "fetch_verification")
knitr::kable(fetch_res)
stopifnot(all(fetch_res$verified))

# ---- load ----
raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))

# ---- structure ----
reps_per_cell <- raw |> count(SubjID, GameID, name = "n_reps")
feedback_map <- raw |> distinct(Trial, Feedback) |> arrange(Trial)
one_tuple_per_game <- raw |>
  distinct(GameID, Ha, pHa, La, LotShapeA, LotNumA,
           Hb, pHb, Lb, LotShapeB, LotNumB, Amb, Corr) |>
  count(GameID) |>
  summarise(max_tuples = max(n)) |>
  pull(max_tuples)

structure_checks <- tribble(
  ~check, ~expected, ~observed,
  "rows", "510750", as.character(nrow(raw)),
  "subjects (distinct SubjID)", "686", as.character(n_distinct(raw$SubjID)),
  "problems (distinct GameID)", "210", as.character(n_distinct(raw$GameID)),
  "columns match dictionary list", "TRUE",
    as.character(identical(names(raw), cpc18_expected_columns)),
  "reps per subject x problem", "all 25",
    paste0("min ", min(reps_per_cell$n_reps), ", max ", max(reps_per_cell$n_reps)),
  "one feature tuple per GameID", "TRUE", as.character(one_tuple_per_game == 1),
  "feedback: trials 1-5 = 0, 6-25 = 1", "TRUE",
    as.character(nrow(feedback_map) == 25 &&
                 all(feedback_map$Feedback == as.integer(feedback_map$Trial > 5))),
  "block = ceiling(Trial/5)", "TRUE",
    as.character(all(raw$block == ceiling(raw$Trial / 5)))
) |>
  mutate(pass = mapply(function(e, o) e == o || grepl("^all 25$", e) &&
                         o == "min 25, max 25", expected, observed))
save_result(structure_checks, "structure_checks")
knitr::kable(structure_checks)
stopifnot(all(structure_checks$pass))

# ---- per-subject-counts ----
per_subject <- raw |>
  group_by(SubjID) |>
  summarise(n_trials = n(), n_problems = n_distinct(GameID), .groups = "drop")

problems_dist <- per_subject |> count(n_problems, name = "n_subjects")
save_result(problems_dist, "problems_per_subject_distribution")
knitr::kable(problems_dist)

trials_summary <- per_subject |>
  summarise(min = min(n_trials), q25 = quantile(n_trials, .25),
            median = median(n_trials), mean = mean(n_trials),
            q75 = quantile(n_trials, .75), max = max(n_trials))
save_result(trials_summary, "trials_per_subject_summary")
knitr::kable(trials_summary)

# ---- clean-subset ----
clean_rows <- raw |> filter(is_clean_row(raw))
features <- add_scaled_features(clean_problem_features(raw))
scales <- attr(features, "scales")

clean_share <- tibble(
  quantity = c("share of rows", "clean problems", "share of problems",
               "loss-informative clean problems (loss in exactly one option)"),
  value = c(nrow(clean_rows) / nrow(raw), nrow(features),
            nrow(features) / n_distinct(raw$GameID),
            sum(features$loss_informative))
)
save_result(clean_share, "clean_subset_share")
knitr::kable(clean_share, digits = 4)

# Corr within the clean subset (Corr != 0 shapes feedback payoffs only; the
# description-based features used here are unaffected).
knitr::kable(clean_rows |> distinct(GameID, Corr) |> count(Corr))

# ---- near-parity ----
thr <- tibble(
  threshold = c("q25 |dEV| (primary EV cutoff)", "median |dEV| (sensitivity cutoff)",
                "median |dSD| (variance cutoff)"),
  value = c(quantile(abs(features$d_ev), .25), median(abs(features$d_ev)),
            median(abs(features$d_sd)))
)
save_result(thr, "near_parity_thresholds")
knitr::kable(thr, digits = 3)

features <- features |>
  mutate(
    near_parity_primary = abs(d_ev) <= quantile(abs(d_ev), .25) &
      abs(d_sd) >= median(abs(d_sd)),
    near_parity_sensitivity = abs(d_ev) <= median(abs(d_ev)) &
      abs(d_sd) >= median(abs(d_sd))
  )

# ---- budgets ----
budgets <- clean_rows |>
  distinct(SubjID, GameID) |>
  left_join(features |>
              select(GameID, loss_informative,
                     near_parity_primary, near_parity_sensitivity),
            by = "GameID") |>
  group_by(SubjID) |>
  summarise(
    clean_problems = n(),
    loss_informative_problems = sum(loss_informative),
    near_parity_problems = sum(near_parity_primary),
    near_parity_problems_sensitivity = sum(near_parity_sensitivity),
    .groups = "drop"
  )
save_result(budgets, "identifiability_budgets_per_subject")

budget_summary <- budgets |>
  pivot_longer(-SubjID) |>
  group_by(name) |>
  summarise(min = min(value), q25 = quantile(value, .25),
            median = median(value), mean = mean(value),
            q75 = quantile(value, .75), max = max(value), .groups = "drop")
save_result(budget_summary, "identifiability_budgets_summary")
knitr::kable(budget_summary, digits = 2)

# ---- sanity-b-rate ----
b_rate <- tibble(
  scope = c("all trials", "clean subset"),
  b_rate = c(mean(raw$B), mean(clean_rows$B)),
  n_trials = c(nrow(raw), nrow(clean_rows))
)
save_result(b_rate, "sanity_b_rate")
knitr::kable(b_rate, digits = 4)

# ---- sanity-monotonicity ----
problem_rates <- clean_rows |>
  group_by(GameID) |>
  summarise(b_rate = mean(B), n_trials = n(), .groups = "drop") |>
  left_join(features |> select(GameID, d_ev), by = "GameID") |>
  mutate(ev_bin = ntile(d_ev, 10))

ev_bins <- problem_rates |>
  group_by(ev_bin) |>
  summarise(
    mean_d_ev = mean(d_ev), b_rate = weighted.mean(b_rate, n_trials),
    n_problems = n(), n_trials = sum(n_trials), .groups = "drop"
  )
save_result(ev_bins, "sanity_ev_monotonicity_bins")
knitr::kable(ev_bins, digits = 3)

mono_rho <- cor(ev_bins$ev_bin, ev_bins$b_rate, method = "spearman")

p_mono <- ggplot(ev_bins, aes(mean_d_ev, b_rate)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_line() + geom_point(aes(size = n_trials), show.legend = FALSE) +
  labs(x = "EV(B) − EV(A), decile-bin mean (clean problems)",
       y = "P(choose B)  [real data]",
       title = "B-choice rate rises with B's EV advantage",
       subtitle = sprintf("Spearman rho over decile bins = %.3f", mono_rho)) +
  theme_minimal()
ggsave(file.path(fig_dir, "step0_ev_monotonicity.png"), p_mono,
       width = 7, height = 4.5, dpi = 150)
p_mono

# ---- sim-validation ----
set.seed(seeds$sim_validation)
val_grid <- tribble(
  ~b0,  ~b_ev, ~b_risk, ~b_loss,
   0.0,   0.2,   -1.5,    -2.0,
   0.0,   0.2,    1.0,     0.5,
   0.0,   2.5,   -1.5,     0.5,
   0.0,   2.5,    1.0,    -2.0,
   0.0,   1.35,  -0.25,   -0.75,
   0.5,   1.0,    0.0,    -1.0
)
val <- validate_simulator(val_grid, features, n_reps = VALID_REPS)
save_result(val, "simulator_validation_full")

val_summary <- val |>
  summarise(
    n_comparisons = n(), max_abs_z = max(abs(z)),
    frac_abs_z_gt3 = mean(abs(z) > 3), mean_z = mean(z),
    max_abs_dev = max(abs(p_sim - p_analytic))
  ) |>
  mutate(pass = frac_abs_z_gt3 < 0.01 & abs(mean_z) < 0.1)
save_result(val_summary, "simulator_validation_summary")
knitr::kable(val_summary, digits = 4)
stopifnot(val_summary$pass)

p_val <- ggplot(val, aes(p_analytic, p_sim)) +
  geom_abline(color = "grey60") +
  geom_point(alpha = 0.3, size = 0.8) +
  labs(x = "Analytic P(choose B)", y = "Simulated frequency",
       title = "Simulator validation [simulation]",
       subtitle = sprintf("%d problem x parameter-set comparisons, %d reps each",
                          nrow(val), VALID_REPS)) +
  theme_minimal()
ggsave(file.path(fig_dir, "step0_simulator_validation.png"), p_val,
       width = 6, height = 5, dpi = 150)
p_val

# ---- recovery-run ----
set.seed(seeds$recovery)
true_params <- draw_true_params(N_SIM)
subj_sets <- split(
  (clean_rows |> distinct(SubjID, GameID))$GameID,
  (clean_rows |> distinct(SubjID, GameID))$SubjID
)
knitr::kable(
  tibble(k = K_GRID,
         eligible_donors = sapply(K_GRID, function(k)
           sum(vapply(subj_sets, length, integer(1)) >= k))),
  caption = "Real subjects eligible as problem-set donors per budget"
)

cells <- expand_grid(k = K_GRID, reps = REPS_GRID)
recovery_raw <- purrr::pmap_dfr(cells, function(k, reps) {
  run_recovery_cell(true_params, subj_sets, features, k, reps)
})
saveRDS(recovery_raw, file.path(res_dir, "recovery_fits_raw.rds"))

flag_summary <- recovery_raw |>
  group_by(k_problems, n_reps) |>
  summarise(
    n = n(), frac_separation = mean(separation),
    frac_degenerate = mean(degenerate), frac_extreme = mean(extreme),
    frac_nonconverged = mean(!converged), .groups = "drop"
  )
save_result(flag_summary, "recovery_fit_flags")
knitr::kable(flag_summary, digits = 3)

# ---- recovery-metrics ----
set.seed(seeds$bootstrap)
param_map <- tribble(
  ~param, ~true_col, ~est_col,
  "b_ev",   "b_ev",   "est_b_ev",
  "b_risk", "b_risk", "est_b_risk",
  "b_loss", "b_loss", "est_b_loss"
)
recovery_metrics <- recovery_raw |>
  group_by(k_problems, n_reps) |>
  group_modify(function(cell, key) {
    purrr::pmap_dfr(param_map, function(param, true_col, est_col) {
      bind_rows(
        boot_cor(cell[[true_col]], cell[[est_col]], "spearman", N_BOOT),
        boot_cor(cell[[true_col]], cell[[est_col]], "pearson", N_BOOT)
      ) |> mutate(param = param, .before = 1)
    })
  }) |>
  ungroup()
save_result(recovery_metrics, "recovery_by_cell")
knitr::kable(recovery_metrics |> filter(method == "spearman"), digits = 3,
             caption = "Recovery (Spearman, primary) [simulation]")
knitr::kable(recovery_metrics |> filter(method == "pearson"), digits = 3,
             caption = "Recovery (Pearson, secondary) [simulation]")

p_rec <- recovery_metrics |>
  filter(method == "spearman") |>
  mutate(budget = paste0(k_problems, " problems x ", n_reps, " reps"),
         trials = k_problems * n_reps) |>
  ggplot(aes(trials, r, color = param)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = 0.08)) +
  geom_line(position = position_dodge(width = 0.08)) +
  scale_x_log10(breaks = sort(unique(as.numeric(outer(K_GRID, REPS_GRID))))) +
  labs(x = "Total trials (problems x reps, log scale)",
       y = "Spearman r (true vs recovered)",
       title = "Parameter recovery vs realistic trial budget [simulation]",
       subtitle = sprintf("Points: %d simulated subjects; bars: 95%% bootstrap CI", N_SIM)) +
  theme_minimal()
ggsave(file.path(fig_dir, "step0_recovery_curves.png"), p_rec,
       width = 8, height = 5, dpi = 150)
p_rec

# ---- separability ----
largest <- recovery_raw |> filter(k_problems == max(K_GRID), n_reps == max(REPS_GRID))
sep_matrix <- sapply(recovery_param_names, function(tc)
  sapply(paste0("est_", recovery_param_names), function(ec)
    cor(largest[[tc]], largest[[ec]], method = "spearman",
        use = "pairwise.complete.obs")))
sep_df <- as_tibble(sep_matrix, rownames = "recovered") |>
  rename_with(~ paste0("true_", .x), -recovered)
save_result(sep_df, "separability_matrix")
knitr::kable(sep_df, digits = 3,
             caption = "Spearman true x recovered, 14 problems x 25 reps [simulation]")

# ---- loss-sweep ----
set.seed(seeds$loss_sweep)
sweep_raw <- run_loss_sweep(true_params, features, L_VALUES, K_SWEEP, 25)
saveRDS(sweep_raw, file.path(res_dir, "loss_sweep_fits_raw.rds"))

set.seed(seeds$bootstrap)
sweep_metrics <- sweep_raw |>
  group_by(l_loss) |>
  group_modify(~ bind_rows(
    boot_cor(.x$b_loss, .x$est_b_loss, "spearman", N_BOOT),
    boot_cor(.x$b_loss, .x$est_b_loss, "pearson", N_BOOT)
  )) |>
  ungroup()
save_result(sweep_metrics, "loss_recovery_by_informative_problems")
knitr::kable(sweep_metrics |> filter(method == "spearman"), digits = 3,
             caption = "b_loss recovery by loss-informative problem count [simulation]")

# Context: realized loss-informative counts in the donor-based 14 x 25 cell
knitr::kable(largest |> count(n_loss_informative),
             caption = "Realized loss-informative counts, donor-based sets (k = 14)")

p_sweep <- sweep_metrics |>
  filter(method == "spearman") |>
  ggplot(aes(l_loss, r)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi)) + geom_line() +
  scale_x_continuous(breaks = L_VALUES) +
  labs(x = "Loss-informative problems in set (of 14 total, 25 reps)",
       y = "Spearman r (true vs recovered b_loss)",
       title = "Loss-aversion recovery vs loss-informative problems [simulation]") +
  theme_minimal()
ggsave(file.path(fig_dir, "step0_loss_sweep.png"), p_sweep,
       width = 7, height = 4.5, dpi = 150)
p_sweep

# ---- verdict ----
crit <- 0.7
verdict <- recovery_metrics |>
  filter(method == "spearman", k_problems == max(K_GRID), n_reps == max(REPS_GRID)) |>
  transmute(param, spearman_r = r, ci_lo, ci_hi,
            recoverable_at_14x25 = r >= crit,
            criterion = paste0("Spearman r >= ", crit, " (pre-declared)"))
save_result(verdict, "recovery_verdict")
knitr::kable(verdict, digits = 3)

# ---- repro ----
knitr::kable(enframe(unlist(seeds), name = "seed_name", value = "seed"))

# ---- appendix ----
for (f in c("R/fetch_cpc18.R", "R/cpc18_features.R", "R/simulate_recovery.R")) {
  cat("\n### ", f, "\n\n```r\n", paste(readLines(here(f)), collapse = "\n"),
      "\n```\n", sep = "")
}

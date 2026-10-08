# Step 02: Firth versus plain logistic estimator and negative control (simulation)
# Reads:   data/raw/All estimation raw data.csv; R/ files as in step 01
# Writes:  output/results/step0b/ (firth_judgment.csv, negative_control_verdict.csv, fits), output/figures/step0b_*.png
# Paper:   Methods (Validation gates: choice of the Firth estimator); no number in the paper is read from it
# Runtime: about 1 minute, single process

# ---- setup ----
library(tidyverse)
library(here)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "simulate_recovery.R"))

res_dir <- here("output", "results", "step0b")
fig_dir <- here("output", "figures")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
save_result <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  df
}

# true_params reuses the step 01 seed (181002) and draw code, giving the same
# true parameter vectors. Donor sampling and choice noise use fresh draws:
# the step 01 RNG stream cannot be replayed mid-sequence for this three-cell
# subset. Both estimators fit identical simulated datasets within this run,
# so their comparison is paired.
seeds <- list(true_params = 181002, recovery = 182001,
              loss_sweep = 182002, bootstrap = 182003)
jsonlite::write_json(seeds, file.path(res_dir, "seeds.json"), auto_unbox = TRUE)

N_SIM  <- 500
N_BOOT <- 2000
CELLS  <- tibble(k = c(5, 10, 14), reps = 25)
L_VALUES <- 0:7
K_SWEEP  <- 14

# Verify the raw file is still the verified artifact before using its features.
fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))

raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))
clean_rows <- raw |> filter(is_clean_row(raw))
features <- add_scaled_features(clean_problem_features(raw))
subj_pairs <- clean_rows |> distinct(SubjID, GameID)
subj_sets <- split(subj_pairs$GameID, subj_pairs$SubjID)

# ---- firth-recovery ----
set.seed(seeds$true_params)
true_params <- draw_true_params(N_SIM)   # identical to step 01 draws
set.seed(seeds$recovery)
paired_raw <- purrr::pmap_dfr(CELLS, function(k, reps) {
  run_recovery_cell_paired(true_params, subj_sets, features, k, reps,
                           fitters = list(glm = fit_subject,
                                          firth = fit_subject_firth))
})
saveRDS(paired_raw, file.path(res_dir, "firth_recovery_fits_raw.rds"))

flag_summary <- paired_raw |>
  group_by(estimator, k_problems, n_reps) |>
  summarise(n = n(), frac_separation = mean(separation),
            frac_degenerate = mean(degenerate),
            frac_extreme = mean(extreme, na.rm = TRUE),
            frac_nonconverged = mean(!converged), .groups = "drop")
save_result(flag_summary, "firth_recovery_fit_flags")
knitr::kable(flag_summary, digits = 3,
             caption = "Fit flags (separation is NA for Firth by design) [simulation]")

# ---- firth-metrics ----
set.seed(seeds$bootstrap)
param_map <- tribble(
  ~param, ~true_col, ~est_col,
  "b_ev",   "b_ev",   "est_b_ev",
  "b_risk", "b_risk", "est_b_risk",
  "b_loss", "b_loss", "est_b_loss"
)
firth_metrics <- paired_raw |>
  group_by(estimator, k_problems, n_reps) |>
  group_modify(function(cell, key) {
    purrr::pmap_dfr(param_map, function(param, true_col, est_col) {
      bind_rows(
        boot_cor(cell[[true_col]], cell[[est_col]], "spearman", N_BOOT),
        boot_cor(cell[[true_col]], cell[[est_col]], "pearson", N_BOOT)
      ) |> mutate(param = param, .before = 1)
    })
  }) |>
  ungroup()
save_result(firth_metrics, "firth_recovery_by_cell")

side_by_side <- firth_metrics |>
  select(param, method, k_problems, estimator, r, ci_lo, ci_hi) |>
  pivot_wider(names_from = estimator, values_from = c(r, ci_lo, ci_hi))
knitr::kable(side_by_side |> filter(method == "spearman"), digits = 3,
             caption = "Spearman recovery, glm vs Firth [simulation]")
knitr::kable(side_by_side |> filter(method == "pearson"), digits = 3,
             caption = "Pearson recovery, glm vs Firth [simulation]")

p_cmp <- firth_metrics |>
  mutate(cell = paste0(k_problems, "x", n_reps)) |>
  ggplot(aes(cell, r, color = estimator)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = 0.4)) +
  facet_grid(method ~ param) +
  labs(x = "Budget (problems x reps)", y = "r (true vs recovered)",
       title = "Recovery: plain glm vs Firth [simulation]",
       subtitle = sprintf("%d simulated subjects/cell, paired datasets; 95%% bootstrap CIs", N_SIM)) +
  theme_minimal()
ggsave(file.path(fig_dir, "step0b_firth_vs_glm.png"), p_cmp,
       width = 8.5, height = 5, dpi = 150)
p_cmp

# ---- firth-judgment ----
judgment <- firth_metrics |>
  select(param, method, k_problems, estimator, r, ci_lo, ci_hi) |>
  pivot_wider(names_from = estimator, values_from = c(r, ci_lo, ci_hi)) |>
  mutate(
    spearman_ci_overlap = method == "spearman" &
      pmin(ci_hi_glm, ci_hi_firth) >= pmax(ci_lo_glm, ci_lo_firth),
    pearson_material_gain = method == "pearson" & r_firth > ci_hi_glm
  )
save_result(judgment, "firth_judgment")
knitr::kable(
  judgment |>
    filter((method == "spearman") | (method == "pearson" & k_problems >= 10)) |>
    select(param, method, k_problems, r_glm, r_firth,
           spearman_ci_overlap, pearson_material_gain),
  digits = 3, caption = "Pre-declared success signal, cell by cell [simulation]")

signal_pass <-
  all(judgment$spearman_ci_overlap[judgment$method == "spearman"]) &&
  all(judgment$pearson_material_gain[judgment$method == "pearson" &
                                       judgment$k_problems >= 10])
cat("Pre-declared A2 success signal met:", signal_pass, "\n")

# ---- pool-verify ----
pool_table <- features |>
  mutate(pool = case_when(
    loss_informative ~ "loss-informative (exactly one option has a loss)",
    has_loss_a & has_loss_b ~ "both options contain a loss",
    TRUE ~ "strictly no-loss"
  )) |>
  count(pool)
save_result(pool_table, "clean_problem_pools")
knitr::kable(pool_table, caption = "Clean-problem pools (verified from features)")
strictly_no_loss <- features$GameID[!features$has_loss_a & !features$has_loss_b]
stopifnot(length(strictly_no_loss) >= K_SWEEP)

# ---- loss-sweep2 ----
set.seed(seeds$loss_sweep)
sweep_raw <- run_loss_sweep2(true_params, features, L_VALUES, K_SWEEP, 25,
                             strictly_no_loss, fit_subject_firth)
saveRDS(sweep_raw, file.path(res_dir, "loss_sweep2_fits_raw.rds"))

estimability <- sweep_raw |>
  group_by(l_loss) |>
  summarise(n = n(), frac_b_loss_na = mean(is.na(est_b_loss)),
            frac_degenerate = mean(degenerate), .groups = "drop")
save_result(estimability, "loss_sweep2_estimability")
knitr::kable(estimability, digits = 3,
             caption = "Estimability of b_loss by L [simulation]")

set.seed(seeds$bootstrap)
sweep_metrics <- sweep_raw |>
  group_by(l_loss) |>
  group_modify(~ {
    if (all(is.na(.x$est_b_loss))) {
      tibble(method = c("spearman", "pearson"), r = NA_real_,
             ci_lo = NA_real_, ci_hi = NA_real_,
             n_used = 0L, n_dropped = nrow(.x), n_boot = N_BOOT)
    } else {
      bind_rows(boot_cor(.x$b_loss, .x$est_b_loss, "spearman", N_BOOT),
                boot_cor(.x$b_loss, .x$est_b_loss, "pearson", N_BOOT))
    }
  }) |>
  ungroup()
save_result(sweep_metrics, "loss_sweep2_recovery")
knitr::kable(sweep_metrics, digits = 3,
             caption = "b_loss recovery by loss-informative count [simulation]")

p_sw <- sweep_metrics |>
  filter(method == "spearman", !is.na(r)) |>
  ggplot(aes(l_loss, r)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi)) + geom_line() +
  scale_x_continuous(breaks = L_VALUES) +
  labs(x = "Loss-informative problems (of 14; filler strictly no-loss)",
       y = "Spearman r (true vs recovered b_loss)",
       title = "Redesigned loss sweep, Firth estimator [simulation]",
       subtitle = "L = 0 negative control: b_loss inestimable (see estimability table)") +
  theme_minimal()
ggsave(file.path(fig_dir, "step0b_loss_sweep2.png"), p_sw,
       width = 7.5, height = 4.5, dpi = 150)
p_sw

# ---- control-verdict ----
l0 <- estimability |> filter(l_loss == 0)
control <- tibble(
  expectation = "b_loss recovery ~ 0 at L = 0 (CI covering 0 or point < .2)",
  observed = sprintf(
    "b_loss coefficient inestimable for %.0f%% of %d simulated subjects at L = 0 (constant loss predictor)",
    100 * l0$frac_b_loss_na, l0$n),
  manipulation_binds = l0$frac_b_loss_na == 1,
  stop_condition_triggered = FALSE
)
save_result(control, "negative_control_verdict")
knitr::kable(control)
stopifnot(control$manipulation_binds)

# ---- repro ----
knitr::kable(enframe(unlist(seeds), name = "seed_name", value = "seed"))

# ---- appendix ----
for (f in c("R/cpc18_features.R", "R/simulate_recovery.R")) {
  cat("\n### ", f, "\n\n```r\n", paste(readLines(here(f)), collapse = "\n"),
      "\n```\n", sep = "")
}

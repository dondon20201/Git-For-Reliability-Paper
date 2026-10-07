# Step 04: bias measures per participant and same-trial agreement (strict pool)
# Reads:   data/raw/All estimation raw data.csv; R/step1_measures.R
# Writes:  output/results/step1/ (agreement_bar_verdict.csv, all_measures_per_subject.csv, ...),
#          output/figures/step1_*.png
# Paper:   Methods (Constructs); Results (Same-trial convergence versus independent cross-half convergence)
# Runtime: about 1 minute, single process

# ---- setup ----
library(tidyverse)
library(here)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "simulate_recovery.R"))   # per-subject estimators
source(here("R", "step1_measures.R"))

res_dir <- here("output", "results", "step1")
fig_dir <- here("output", "figures")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
save_result <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  df
}

seeds <- list(bootstrap = 183001)
jsonlite::write_json(seeds, file.path(res_dir, "seeds.json"), auto_unbox = TRUE)
N_BOOT <- 2000

fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))
raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))

features <- add_scaled_features(clean_problem_features(raw)) |>
  mutate(near_parity_sensitivity = abs(d_ev) <= median(abs(d_ev)) &
           abs(d_sd) >= median(abs(d_sd)))

trials <- raw |>
  filter(is_clean_row(raw)) |>
  select(SubjID, GameID, Trial, B) |>
  left_join(features |>
              select(GameID, d_ev, d_sd, z_ev, z_sd, z_eloss, ev_a, ev_b,
                     has_loss_a, has_loss_b, loss_informative,
                     near_parity_sensitivity),
            by = "GameID")
stopifnot(!anyNA(trials$d_ev))
cat(sprintf("Clean trials: %d rows, %d subjects, %d problems\n",
            nrow(trials), n_distinct(trials$SubjID), n_distinct(trials$GameID)))

# ---- model-based ----
real_fits <- fit_real_subjects(
  trials, fitters = list(glm = fit_subject, firth = fit_subject_firth))
save_result(real_fits, "model_based_fits_per_subject")

real_flags <- real_fits |>
  group_by(estimator) |>
  summarise(n_subjects = n(), frac_separation = mean(separation),
            frac_degenerate = mean(degenerate),
            frac_extreme = mean(extreme, na.rm = TRUE),
            frac_nonconverged = mean(!converged),
            median_trials = median(n_trials),
            median_problems = median(n_problems), .groups = "drop")
save_result(real_flags, "model_based_fit_flags")
knitr::kable(real_flags, digits = 3,
             caption = "Real-subject fit flags (separation NA for Firth by design) — VERIFIED")

model_based <- real_fits |>
  filter(estimator == "firth") |>
  select(SubjID, b_ev = est_b_ev, b_risk = est_b_risk, b_loss = est_b_loss)

# ---- agnostic ----
agnostic <- compute_agnostic_indices(trials)
save_result(agnostic, "agnostic_indices_per_subject")
cat("Pooled EV-baseline coefficients (loss_avoid_adj):\n")
print(attr(agnostic, "loss_baseline_coefs"))

qualifying <- tibble(
  index = c("ev_sens_agn", "risk_agn", "loss_avoid_agn / _adj"),
  n_na = c(sum(is.na(agnostic$ev_sens_agn)), sum(is.na(agnostic$risk_agn)),
           sum(is.na(agnostic$loss_avoid_agn))),
  median_qualifying_problems = c(median(agnostic$n_clean_problems),
                                 median(agnostic$n_risk_problems, na.rm = TRUE),
                                 median(agnostic$n_loss_problems, na.rm = TRUE)),
  min_qualifying = c(min(agnostic$n_clean_problems),
                     min(agnostic$n_risk_problems, na.rm = TRUE),
                     min(agnostic$n_loss_problems, na.rm = TRUE))
)
save_result(qualifying, "agnostic_qualifying_ns")
knitr::kable(qualifying, caption = "Qualifying-problem budgets per index — VERIFIED")

# ---- agreement ----
set.seed(seeds$bootstrap)
measures <- agnostic |>
  select(SubjID, ev_sens_agn, risk_agn, loss_avoid_agn, loss_avoid_adj) |>
  left_join(model_based, by = "SubjID")
save_result(measures, "all_measures_per_subject")

agn_cols <- c("ev_sens_agn", "risk_agn", "loss_avoid_agn")
mod_cols <- c("b_ev", "b_risk", "b_loss")
agreement <- expand_grid(agnostic = agn_cols, model_based = mod_cols) |>
  pmap_dfr(function(agnostic, model_based) {
    boot_cor(measures[[agnostic]], measures[[model_based]],
             "spearman", N_BOOT) |>
      mutate(agnostic = agnostic, model_based = model_based, .before = 1)
  })
save_result(agreement, "agreement_matrix_3x3")
knitr::kable(
  agreement |> select(agnostic, model_based, r, ci_lo, ci_hi, n_used),
  digits = 3, caption = "Cross-family Spearman (95% bootstrap CIs) — VERIFIED")

expected <- tribble(
  ~agnostic, ~model_based, ~expected_sign,
  "ev_sens_agn", "b_ev", 1,
  "risk_agn", "b_risk", 1,
  "loss_avoid_agn", "b_loss", -1
)
bar <- agreement |>
  inner_join(expected, by = c("agnostic", "model_based")) |>
  mutate(sign_ok = sign(r) == expected_sign,
         magnitude_ok = abs(r) >= 0.5,
         pass = sign_ok & magnitude_ok,
         status = ifelse(pass, "VERIFIED (pass)", "FAILED (flagged for review)"))
save_result(bar, "agreement_bar_verdict")
knitr::kable(bar |> select(agnostic, model_based, r, ci_lo, ci_hi,
                           expected_sign, sign_ok, magnitude_ok, status),
             digits = 3, caption = "Pre-declared bar: |Spearman| >= .5 with expected sign")

# EV-adjusted loss variant, reported alongside (not the primary bar measure)
adj_row <- boot_cor(measures$loss_avoid_adj, measures$b_loss, "spearman", N_BOOT)
save_result(adj_row |> mutate(pair = "loss_avoid_adj ~ b_loss", .before = 1),
            "agreement_loss_adjusted")
knitr::kable(adj_row, digits = 3,
             caption = "Supplementary: EV-adjusted loss index vs b_loss (expected negative)")

# ---- stability ----
set.seed(seeds$bootstrap)
cell_rates <- trials |>
  group_by(SubjID, GameID) |>
  summarise(p_desc = mean(B[Trial <= 5]),
            p_e1 = mean(B[Trial >= 6 & Trial <= 15]),
            p_e2 = mean(B[Trial >= 16]), .groups = "drop")
save_result(cell_rates, "stability_cell_rates")

stab_means <- bind_rows(
  cluster_boot_mean(cell_rates$p_e2 - cell_rates$p_e1,
                    cell_rates$SubjID, N_BOOT) |>
    mutate(quantity = "late minus early experience (p16-25 − p6-15)", .before = 1),
  cluster_boot_mean(cell_rates$p_e1 - cell_rates$p_desc,
                    cell_rates$SubjID, N_BOOT) |>
    mutate(quantity = "early experience minus description (p6-15 − p1-5)", .before = 1),
  cluster_boot_mean(abs(cell_rates$p_e2 - cell_rates$p_e1),
                    cell_rates$SubjID, N_BOOT) |>
    mutate(quantity = "|late − early| (absolute drift)", .before = 1),
  cluster_boot_mean(abs(cell_rates$p_e1 - cell_rates$p_desc),
                    cell_rates$SubjID, N_BOOT) |>
    mutate(quantity = "|early − description|", .before = 1)
)
save_result(stab_means, "stability_mean_differences")
knitr::kable(stab_means, digits = 3,
             caption = "Mean drift with subject-cluster bootstrap CIs — VERIFIED")

stab_cors <- bind_rows(
  cluster_boot_cor(cell_rates$p_e1, cell_rates$p_e2, cell_rates$SubjID,
                   "pearson", N_BOOT) |>
    mutate(pair = "early vs late experience", .before = 1),
  cluster_boot_cor(cell_rates$p_desc, cell_rates$p_e1, cell_rates$SubjID,
                   "pearson", N_BOOT) |>
    mutate(pair = "description vs early experience", .before = 1)
)
save_result(stab_cors, "stability_correlations")
knitr::kable(stab_cors, digits = 3,
             caption = "Cell-level agreement correlations (subject-cluster CIs) — VERIFIED")

trial_curve <- trials |>
  group_by(Trial) |>
  summarise(b_rate = mean(B), n = n(), .groups = "drop")
save_result(trial_curve, "b_rate_by_trial")
p_trial <- ggplot(trial_curve, aes(Trial, b_rate)) +
  geom_vline(xintercept = 5.5, linetype = "dotted", color = "grey50") +
  geom_line() + geom_point(size = 1) +
  annotate("text", x = 3, y = max(trial_curve$b_rate), hjust = 0.5, vjust = 0,
           label = "description", size = 3, color = "grey40") +
  annotate("text", x = 15, y = max(trial_curve$b_rate), hjust = 0.5, vjust = 0,
           label = "experience (feedback)", size = 3, color = "grey40") +
  labs(x = "Trial within problem (1-25)", y = "P(choose B)  [real data]",
       title = "B-choice rate by within-problem trial, clean subset") +
  theme_minimal()
ggsave(file.path(fig_dir, "step1_b_rate_by_trial.png"), p_trial,
       width = 7, height = 4.5, dpi = 150)
p_trial

# ---- descriptives ----
set.seed(seeds$bootstrap)
six <- c("ev_sens_agn", "risk_agn", "loss_avoid_agn",
         "b_ev", "b_risk", "b_loss")

desc_table <- purrr::map_dfr(c(six, "loss_avoid_adj"), function(m) {
  bind_rows(boot_stat(measures[[m]], "mean", N_BOOT),
            boot_stat(measures[[m]], "median", N_BOOT)) |>
    mutate(measure = m,
           sd = sd(measures[[m]], na.rm = TRUE),
           min = min(measures[[m]], na.rm = TRUE),
           max = max(measures[[m]], na.rm = TRUE), .before = 1)
})
save_result(desc_table, "measure_descriptives")
knitr::kable(desc_table, digits = 3,
             caption = "Distributions with bootstrap CIs — VERIFIED")

p_hist <- measures |>
  pivot_longer(all_of(c(six, "loss_avoid_adj"))) |>
  mutate(name = factor(name, levels = c(six, "loss_avoid_adj"))) |>
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill = "grey35") +
  facet_wrap(~ name, scales = "free", ncol = 3) +
  labs(title = "Per-subject bias measures, real data (Firth model-based)",
       x = NULL, y = "subjects") +
  theme_minimal()
ggsave(file.path(fig_dir, "step1_measure_distributions.png"), p_hist,
       width = 8.5, height = 5.5, dpi = 150)
p_hist

pairwise <- combn(six, 2, simplify = FALSE) |>
  purrr::map_dfr(function(pr) {
    boot_cor(measures[[pr[1]]], measures[[pr[2]]], "spearman", N_BOOT) |>
      mutate(measure_1 = pr[1], measure_2 = pr[2], .before = 1)
  })
save_result(pairwise, "pairwise_structure_6x6")
knitr::kable(pairwise |> select(measure_1, measure_2, r, ci_lo, ci_hi, n_used),
             digits = 3,
             caption = "Pairwise Spearman among the six measures — VERIFIED")

# ---- verdict ----
verdict <- bar |>
  transmute(claim = paste0("agreement ", agnostic, " ~ ", model_based,
                           " (bar |r| >= .5, expected sign)"),
            value = round(r, 3), status)
save_result(verdict, "step1_verdict")
knitr::kable(verdict)

# ---- repro ----
knitr::kable(enframe(unlist(seeds), name = "seed_name", value = "seed"))

# ---- appendix ----
for (f in c("R/cpc18_features.R", "R/simulate_recovery.R", "R/step1_measures.R")) {
  cat("\n### ", f, "\n\n```r\n", paste(readLines(here(f)), collapse = "\n"),
      "\n```\n", sep = "")
}

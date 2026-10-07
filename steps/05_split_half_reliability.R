# Step 05: split-half reliability curves (strict pool)
# Reads:   data/raw/All estimation raw data.csv; R/step2_reliability.R
# Writes:  output/results/step2/ (reliability_by_cell_construct.csv, trials_needed_table.csv,
#          seeds_cells.csv, seeds_bootstrap.csv, sanity checks), output/figures/step2_*.png
# Paper:   Results (Reliability across observed distinct-problem budgets; Repetitions and loss exposure;
#          Description versus experience); Figure 1
# Runtime: about 29 minutes with 12 parallel workers

# ---- setup ----
library(tidyverse)
library(here)
library(parallel)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))

res_dir <- here("output", "results", "step2")
fig_dir <- here("output", "figures")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
save_result <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  df
}

N_SPLITS <- 200
N_BOOT   <- 2000
SEED_CELL_BASE <- 185000   # cell i uses seed SEED_CELL_BASE + i
SEED_BOOT_BASE <- 186000   # metrics task j uses seed SEED_BOOT_BASE + j

fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))
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
cat(sprintf("%d subjects; loss_avoid_adj baselines (intercept, slope): exp = (%.3f, %.3f), desc = (%.3f, %.3f)\n",
            length(subj_data), baselines$exp[1], baselines$exp[2],
            baselines$desc[1], baselines$desc[2]))

# ---- cells ----
cell_specs <- c(
  lapply(2:7, function(k) list(axis = "i_exp", regime = "exp", k = k, r = 20, lambda = NA)),
  lapply(c(5, 10), function(r) list(axis = "ii_reps", regime = "exp", k = 5, r = r, lambda = NA)),
  lapply(1:3, function(l) list(axis = "iii_lambda", regime = "exp", k = 5, r = 20, lambda = l)),
  lapply(2:7, function(k) list(axis = "i_desc", regime = "desc", k = k, r = 5, lambda = NA))
)
cell_specs <- lapply(seq_along(cell_specs), function(i) {
  s <- cell_specs[[i]]
  s$cell_id <- sprintf("%s_k%d_r%d%s", s$axis, s$k, s$r,
                       ifelse(is.na(s$lambda), "", paste0("_l", s$lambda)))
  s$seed <- SEED_CELL_BASE + i
  s
})
seeds_df <- purrr::map_dfr(cell_specs, ~ tibble(cell_id = .x$cell_id, seed = .x$seed))
save_result(seeds_df, "seeds_cells")
knitr::kable(purrr::map_dfr(cell_specs, ~ as_tibble(.x[c("cell_id", "axis", "regime", "k", "r", "lambda", "seed")])))

# ---- estimate ----
n_workers <- max(1, min(12, detectCores() - 2))
cat("workers:", n_workers, "\n")
cl <- makeCluster(n_workers)
invisible(clusterCall(cl, function(lp, root) {
  .libPaths(lp)
  source(file.path(root, "R", "step2_reliability.R"))
  NULL
}, .libPaths(), here::here()))
clusterExport(cl, c("subj_data", "N_SPLITS"), envir = environment())

t0 <- Sys.time()
cells <- parLapplyLB(cl, cell_specs, function(spec) {
  run_cell(spec, subj_data, N_SPLITS, seed = spec$seed)
})
names(cells) <- vapply(cells, function(x) x$spec$cell_id, character(1))
cat("estimation wall time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---- metrics ----
tasks <- purrr::list_flatten(lapply(cells, function(cell) {
  constructs <- if (cell$spec$axis == "iii_lambda") {
    c("b_loss", "loss_avoid_adj")
  } else step2_constructs
  lapply(constructs, function(con) {
    list(cell_id = cell$spec$cell_id, spec = cell$spec, construct = con,
         E1 = cell$est[, , 1, con], E2 = cell$est[, , 2, con])
  })
}))
tasks <- lapply(seq_along(tasks), function(j) {
  t <- tasks[[j]]; t$seed <- SEED_BOOT_BASE + j; t
})
save_result(purrr::map_dfr(tasks, ~ tibble(cell_id = .x$cell_id,
                                           construct = .x$construct,
                                           seed = .x$seed)),
            "seeds_bootstrap")

clusterExport(cl, "N_BOOT", envir = environment())
t0 <- Sys.time()
metrics_list <- parLapplyLB(cl, tasks, function(t) {
  if (!any(is.finite(t$E1) & is.finite(t$E2))) {
    out <- tibble::tibble(
      r_raw = NA_real_, ci_lo_raw = NA_real_, ci_hi_raw = NA_real_,
      r_sb = NA_real_, ci_lo_sb = NA_real_, ci_hi_sb = NA_real_,
      n_subjects = nrow(t$E1), n_splits_used = 0L, n_splits = ncol(t$E1),
      finite_frac = 0, n_boot = N_BOOT)
  } else {
    out <- reliability_point(t$E1, t$E2, N_BOOT, seed = t$seed)
  }
  dplyr::bind_cols(
    tibble::tibble(cell_id = t$cell_id, axis = t$spec$axis,
                   regime = t$spec$regime, k = t$spec$k, r = t$spec$r,
                   lambda = t$spec$lambda, construct = t$construct), out)
})
stopCluster(cl)
metrics <- bind_rows(metrics_list) |>
  mutate(status = case_when(
    finite_frac == 0 ~ "unidentified (k < 4 problems/half for 4-parameter fit)",
    n_splits_used < n_splits ~ sprintf("VERIFIED (%d/%d splits usable)", n_splits_used, n_splits),
    TRUE ~ "VERIFIED"))
cat("metrics wall time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
save_result(metrics, "reliability_by_cell_construct")
knitr::kable(metrics |> select(-cell_id, -n_boot), digits = 3)

# ---- curves-fig ----
axis_i <- metrics |>
  filter(axis %in% c("i_exp", "i_desc"), is.finite(r_sb)) |>
  mutate(regime_lab = ifelse(regime == "exp",
                             "experience (r = 20)", "description (r = 5)"),
         total_problems = 2 * k)

p_i <- ggplot(axis_i, aes(total_problems, r_sb, color = construct)) +
  geom_hline(yintercept = c(.5, .7, .8), linetype = "dotted", color = "grey70") +
  geom_pointrange(aes(ymin = ci_lo_sb, ymax = ci_hi_sb),
                  position = position_dodge(width = 0.5), size = 0.3) +
  geom_line(position = position_dodge(width = 0.5)) +
  facet_wrap(~ regime_lab) +
  scale_x_continuous(breaks = seq(4, 14, 2)) +
  labs(x = "Total problems (2k; halves of k)",
       y = "Split-half reliability (Spearman, SB-corrected)",
       title = "Rank reliability vs problem budget [real data]",
       subtitle = "SB correction labeled approximate for rank correlations; 95% subject-cluster bootstrap CIs; model-based constructs unidentified at k < 4") +
  theme_minimal()
ggsave(file.path(fig_dir, "step2_reliability_problems.png"), p_i,
       width = 9, height = 5.5, dpi = 150)
p_i

# ---- curves-fig2 ----
p_ii <- metrics |>
  filter((axis == "ii_reps" | (axis == "i_exp" & k == 5)), is.finite(r_sb)) |>
  ggplot(aes(r, r_sb, color = construct)) +
  geom_hline(yintercept = c(.5, .7), linetype = "dotted", color = "grey70") +
  geom_pointrange(aes(ymin = ci_lo_sb, ymax = ci_hi_sb),
                  position = position_dodge(width = 0.8), size = 0.3) +
  geom_line(position = position_dodge(width = 0.8)) +
  scale_x_continuous(breaks = c(5, 10, 20)) +
  labs(x = "Reps per problem (experience trials; k = 5 problems/half)",
       y = "SB-corrected Spearman reliability",
       title = "Rank reliability vs within-problem reps [real data]") +
  theme_minimal()

p_iii <- metrics |>
  filter(axis == "iii_lambda", is.finite(r_sb)) |>
  ggplot(aes(lambda, r_sb, color = construct)) +
  geom_hline(yintercept = c(.5, .7), linetype = "dotted", color = "grey70") +
  geom_pointrange(aes(ymin = ci_lo_sb, ymax = ci_hi_sb),
                  position = position_dodge(width = 0.15), size = 0.3) +
  geom_line(position = position_dodge(width = 0.15)) +
  scale_x_continuous(breaks = 1:3) +
  labs(x = "Loss-informative problems per half (stratified; k = 5, r = 20)",
       y = "SB-corrected Spearman reliability",
       title = "Loss constructs vs loss-informative budget [real data]") +
  theme_minimal()
ggsave(file.path(fig_dir, "step2_reliability_reps.png"), p_ii,
       width = 7, height = 4.5, dpi = 150)
ggsave(file.path(fig_dir, "step2_reliability_lambda.png"), p_iii,
       width = 7, height = 4.5, dpi = 150)
p_ii; p_iii

# ---- trials-needed ----
thresholds <- c(0.5, 0.7, 0.8)
axis_order <- list(
  i_exp = list(var = "k", label = "problems (experience, r=20)"),
  ii_reps = list(var = "r", label = "reps/problem (k=5, experience)"),
  iii_lambda = list(var = "lambda", label = "loss-informative problems/half (k=5)"),
  i_desc = list(var = "k", label = "problems (description, r=5)")
)
trials_needed <- purrr::imap_dfr(axis_order, function(ax, ax_name) {
  m <- metrics |> filter(axis == ax_name |
                           (ax_name == "ii_reps" & axis == "i_exp" & k == 5))
  if (ax_name == "ii_reps") m <- m |> mutate(budget = r)
  else if (ax_name == "iii_lambda") m <- m |> mutate(budget = lambda)
  else m <- m |> mutate(budget = k)
  m |>
    group_by(construct) |>
    group_modify(function(g, key) {
      purrr::map_dfr(thresholds, function(th) {
        hit <- g |> filter(is.finite(r_sb), r_sb >= th) |> arrange(budget) |> slice(1)
        if (nrow(hit) == 0) {
          tibble(threshold = th, budget = NA_real_, total_problems = NA_real_,
                 total_trials = NA_real_, r_sb = NA_real_, ci_lo_sb = NA_real_,
                 ci_hi_sb = NA_real_, n_subjects = NA_integer_,
                 entry = "not reached within CPC18 budgets")
        } else {
          tibble(threshold = th, budget = hit$budget,
                 total_problems = 2 * hit$k, total_trials = 2 * hit$k * hit$r,
                 r_sb = hit$r_sb, ci_lo_sb = hit$ci_lo_sb,
                 ci_hi_sb = hit$ci_hi_sb, n_subjects = hit$n_subjects,
                 entry = "reached")
        }
      })
    }) |>
    ungroup() |>
    mutate(axis = ax_name, axis_label = ax$label, .before = 1)
})
save_result(trials_needed, "trials_needed_table")
knitr::kable(trials_needed |>
               select(axis_label, construct, threshold, budget, total_problems,
                      total_trials, r_sb, ci_lo_sb, ci_hi_sb, entry),
             digits = 3)

# ---- sanity ----
anchor <- metrics |>
  filter(axis == "i_exp", construct %in% c("b_ev", "b_loss"), is.finite(r_sb)) |>
  select(k, construct, r_sb) |>
  pivot_wider(names_from = construct, values_from = r_sb) |>
  mutate(b_ev_dominates = b_ev >= b_loss)

mono <- metrics |>
  filter(axis %in% c("i_exp", "i_desc"), is.finite(r_sb)) |>
  group_by(axis, construct) |>
  arrange(k, .by_group = TRUE) |>
  summarise(n_points = n(),
            n_decreasing_steps = sum(diff(r_sb) < 0),
            max_decrease = ifelse(n() > 1, max(0, -min(diff(r_sb))), 0),
            spearman_vs_budget = ifelse(n() > 1, cor(k, r_sb, method = "spearman"), NA),
            .groups = "drop")
save_result(anchor, "sanity_dominance")
save_result(mono, "sanity_monotonicity")
knitr::kable(anchor, digits = 3, caption = "b_ev vs b_loss dominance by budget")
knitr::kable(mono, digits = 3, caption = "Monotonicity diagnostics per curve")

# ---- repro ----
jsonlite::write_json(
  list(seed_cell_base = SEED_CELL_BASE, seed_boot_base = SEED_BOOT_BASE,
       n_splits = N_SPLITS, n_boot = N_BOOT),
  file.path(res_dir, "seeds.json"), auto_unbox = TRUE)

# ---- appendix ----
for (f in c("R/cpc18_features.R", "R/step2_reliability.R")) {
  cat("\n### ", f, "\n\n```r\n", paste(readLines(here(f)), collapse = "\n"),
      "\n```\n", sep = "")
}

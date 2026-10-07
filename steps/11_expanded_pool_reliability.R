# Step 11: split-half reliability curves (expanded pool)
# Reads:   data/raw/All estimation raw data.csv; output/results/step2/ (step 05); output/results/step4/ (steps 09-10)
# Writes:  output/results/step4/ (reliability_extended.csv, extended_interpretation.csv, seeds_*_extended.csv,
#          overlap_clean_vs_expanded.csv), output/figures/step4_reliability_extended.png
# Paper:   Results (Reliability differs across problem pools); Figure 1
# Runtime: about 17 minutes with 12 parallel workers

# Extended reliability curves on the expanded pool (all Amb == 0 problems,
# features from the validated reconstruction) and the clean-vs-expanded
# overlap comparison.
#
# REAL DATA. Same machinery and pre-declared design constants as step 05
# (R/step2_reliability.R): rank-based split-half (problems split
# whole), raw + Spearman-Brown, N_SPLITS = 200, subject-cluster bootstrap
# N_BOOT = 2000 (bootstrap outside, splits inside), seeds recorded, k < 4
# model-based cells reported as unidentified. Measure definitions carry over
# with ONLY the features recomputed from reconstructed distributions:
#   - ev_sens_agn uses reconstructed EVs (d_ev);
#   - loss_avoid_adj's "contains a loss" = any outcome < 0 in the option's
#     reconstructed support (has_loss), with the pooled EV baselines refit on
#     the expanded pool;
#   - model-based b_ev/b_risk/b_loss use the SAME 4-parameter Firth formula
#     with z_ev/z_sd/z_eloss recomputed (scaled by expanded-pool SDs;
#     monotone-linear, cannot affect rank statistics).
# Realized grid (pre-declared rule, computed in step 10): experience regime,
# r = 20, k = 2..13. The step 10 recovery re-gate PASSED for all three model
# parameters, so all five constructs are interpretable.
#
# Overlap comparison, pre-declared stance: clean-vs-expanded differences at
# the same k are RESULTS (problem-mix dependence), never errors to adjust away.

library(tidyverse)
library(here)
library(parallel)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "fetch_official_cpc18.R"))
source(here("R", "lotshape_reconstruction.R"))
source(here("R", "step2_reliability.R"))

res_dir <- here("output", "results", "step4")
fig_dir <- here("output", "figures")

N_SPLITS <- 200
N_BOOT   <- 2000
SEED_CELL_BASE <- 195000   # cell i uses seed SEED_CELL_BASE + i
SEED_BOOT_BASE <- 196000   # metrics task j uses seed SEED_BOOT_BASE + j

# ---- gates: step 09 validation + step 10 recovery re-gate must exist --------
gate <- read_csv(file.path(res_dir, "reconstruction_validation.csv"),
                 show_col_types = FALSE)
stopifnot(all(gate$pass))
verdict_a4 <- read_csv(file.path(res_dir, "expanded_recovery_verdict.csv"),
                       show_col_types = FALSE)
recovery_failed <- verdict_a4$param[!verdict_a4$pass]

# ---- expanded-pool subject data ---------------------------------------------
fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))
raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))

param_cols <- c("Ha", "pHa", "La", "LotShapeA", "LotNumA",
                "Hb", "pHb", "Lb", "LotShapeB", "LotNumB", "Amb", "Corr")
problems_exp <- raw |>
  distinct(GameID, across(all_of(param_cols))) |>
  arrange(GameID) |>
  filter(Amb == 0)
feat_exp <- expanded_problem_features(problems_exp,
                                      reconstruct_all_problems(problems_exp))

trials <- raw |>
  semi_join(problems_exp, by = "GameID") |>
  select(SubjID, GameID, Trial, B) |>
  left_join(feat_exp |>
              select(GameID, z_ev, z_sd, z_eloss, d_ev, loss_informative,
                     has_loss_b, ev_a, ev_b), by = "GameID")

baselines <- fit_loss_baselines(trials)
subj_data <- build_subject_data(trials, baselines)
cat(sprintf("%d subjects; expanded loss_avoid_adj baselines: exp = (%.3f, %.3f), desc = (%.3f, %.3f)\n",
            length(subj_data), baselines$exp[1], baselines$exp[2],
            baselines$desc[1], baselines$desc[2]))

# ---- cells (realized grid from step 10) -------------------------------------
realized <- read_csv(file.path(res_dir, "realized_grid.csv"),
                     show_col_types = FALSE)
cell_specs <- lapply(seq_len(nrow(realized)), function(i) {
  k <- realized$k[i]
  list(axis = "i_exp_ext", regime = "exp", k = k, r = 20, lambda = NA,
       cell_id = sprintf("ext_exp_k%d_r20", k),
       seed = SEED_CELL_BASE + i)
})
write_csv(purrr::map_dfr(cell_specs, ~ tibble(cell_id = .x$cell_id, seed = .x$seed)),
          file.path(res_dir, "seeds_cells_extended.csv"))

# ---- split-half estimation (parallel; same pattern as step 05) --------------
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
cat("estimation wall time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---- reliability metrics (bootstrap outside, splits inside) -----------------
tasks <- purrr::list_flatten(lapply(cells, function(cell) {
  lapply(step2_constructs, function(con) {
    list(cell_id = cell$spec$cell_id, spec = cell$spec, construct = con,
         E1 = cell$est[, , 1, con], E2 = cell$est[, , 2, con])
  })
}))
tasks <- lapply(seq_along(tasks), function(j) {
  t <- tasks[[j]]; t$seed <- SEED_BOOT_BASE + j; t
})
write_csv(purrr::map_dfr(tasks, ~ tibble(cell_id = .x$cell_id,
                                         construct = .x$construct,
                                         seed = .x$seed)),
          file.path(res_dir, "seeds_bootstrap_extended.csv"))

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
    construct %in% recovery_failed ~ "recovery-failed on expanded design — not interpreted",
    n_splits_used < n_splits ~ sprintf("VERIFIED (%d/%d splits usable)", n_splits_used, n_splits),
    TRUE ~ "VERIFIED"))
cat("metrics wall time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
write_csv(metrics, file.path(res_dir, "reliability_extended.csv"))
print(metrics |> select(cell_id, construct, r_raw, r_sb, ci_lo_sb, ci_hi_sb,
                        n_subjects, status), n = 60)

# ---- overlap comparison at the same k ---------------------------------------
step2_metrics <- read_csv(here("output", "results", "step2", "reliability_by_cell_construct.csv"),
                          show_col_types = FALSE)
overlap <- step2_metrics |>
  filter(axis == "i_exp") |>
  select(construct, k, r,
         r_raw_clean = r_raw, r_sb_clean = r_sb,
         ci_lo_sb_clean = ci_lo_sb, ci_hi_sb_clean = ci_hi_sb,
         n_subjects_clean = n_subjects, status_clean = status) |>
  inner_join(
    metrics |>
      select(construct, k, r,
             r_raw_expanded = r_raw, r_sb_expanded = r_sb,
             ci_lo_sb_expanded = ci_lo_sb, ci_hi_sb_expanded = ci_hi_sb,
             n_subjects_expanded = n_subjects, status_expanded = status),
    by = c("construct", "k", "r")) |>
  mutate(diff_sb_expanded_minus_clean = r_sb_expanded - r_sb_clean,
         stance = "pre-declared: same-k differences are RESULTS (problem-mix dependence), not errors") |>
  arrange(construct, k)
write_csv(overlap, file.path(res_dir, "overlap_clean_vs_expanded.csv"))
print(overlap |> select(construct, k, r_sb_clean, r_sb_expanded,
                        diff_sb_expanded_minus_clean), n = 40)

# ---- pre-declared interpretation readouts (no extrapolation) ----------------
first_cross <- function(m, con, thr) {
  ks <- m |> filter(construct == con, is.finite(r_sb), r_sb >= thr) |> pull(k)
  if (length(ks) == 0) "not reached within realized grid" else as.character(min(ks))
}
interp <- bind_rows(
  tibble(question = "loss_avoid_adj crosses SB .5 at k (problems/half)",
         answer = first_cross(metrics, "loss_avoid_adj", 0.5)),
  tibble(question = "loss_avoid_adj crosses SB .7 at k (problems/half)",
         answer = first_cross(metrics, "loss_avoid_adj", 0.7)),
  map_dfr(step2_constructs, function(con) {
    tibble(question = sprintf("%s reaches SB .7 within realized grid (k<=13)", con),
           answer = ifelse(first_cross(metrics, con, 0.7) == "not reached within realized grid",
                           "no", sprintf("yes (k = %s)", first_cross(metrics, con, 0.7))))
  })
) |>
  mutate(rule = "pre-declared; never extrapolate beyond the realized grid")
write_csv(interp, file.path(res_dir, "extended_interpretation.csv"))
print(interp, width = Inf)

# ---- figure: clean (step 05) vs expanded (this step) curves -----------------
plot_df <- bind_rows(
  step2_metrics |> filter(axis == "i_exp", is.finite(r_sb)) |>
    mutate(pool = "clean two-outcome (Step 2, 105-problem pool)"),
  metrics |> filter(is.finite(r_sb)) |>
    mutate(pool = "expanded Amb==0 (Step 4, 182-problem pool)")
) |> mutate(total_problems = 2 * k)
p <- ggplot(plot_df, aes(total_problems, r_sb, color = construct,
                         linetype = pool, shape = pool)) +
  geom_hline(yintercept = c(.5, .7, .8), linetype = "dotted", color = "grey70") +
  geom_pointrange(aes(ymin = ci_lo_sb, ymax = ci_hi_sb),
                  position = position_dodge(width = 0.5), size = 0.25) +
  geom_line(position = position_dodge(width = 0.5)) +
  facet_wrap(~ construct) +
  scale_x_continuous(breaks = seq(4, 26, 4)) +
  guides(color = "none") +
  labs(x = "Total problems (2k; halves of k)",
       y = "Split-half reliability (Spearman, SB-corrected)",
       title = "Extended reliability curves: clean vs expanded problem pool [real data]",
       subtitle = "Experience regime, r = 20. SB approximate for rank correlations; 95% subject-cluster bootstrap CIs;\nmodel-based constructs unidentified at k < 4; same-k differences are results (problem-mix dependence)") +
  theme_minimal()
ggsave(file.path(fig_dir, "step4_reliability_extended.png"), p,
       width = 10, height = 7, dpi = 150)
cat("\nDone: reliability_extended.csv, overlap_clean_vs_expanded.csv, extended_interpretation.csv, figure saved.\n")

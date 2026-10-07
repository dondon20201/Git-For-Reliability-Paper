# Step 06: same-trial versus cross-half agreement between measures
# Reads:   data/raw/All estimation raw data.csv; output/results/step2/reliability_by_cell_construct.csv (step 05)
# Writes:  output/results/step1_addendum/ (cross_half_agreement.csv, ceiling_check.csv, seeds.json)
# Paper:   Results (Same-trial convergence versus independent cross-half convergence)
# Runtime: about 13 minutes, single process

# ---- setup ----
library(tidyverse)
library(here)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))
source(here("R", "cross_half_agreement.R"))

res_dir <- here("output", "results", "step1_addendum")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

N_SPLITS <- 200
N_BOOT   <- 2000
SEED_CELL     <- 187001   # split-half estimation cell
SEED_BOOT_EV   <- 188001  # bootstrap, ev pair
SEED_BOOT_LOSS <- 188002  # bootstrap, loss pair

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

# ---- estimate ----
spec <- list(cell_id = "addendum_i_exp_k7_r20", axis = "i_exp",
             regime = "exp", k = 7, r = 20, lambda = NA)
t0 <- Sys.time()
cell <- run_cell(spec, subj_data, N_SPLITS, seed = SEED_CELL)
cat("n eligible subjects:", cell$n_eligible, "\n")
stopifnot(cell$n_eligible == 446)  # the step 05 k=7 pool, as pre-declared
cat("estimation wall time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---- agreement ----
pairs <- list(
  list(pair = "ev_sens_agn~b_ev",     agn = "ev_sens_agn",
       mod = "b_ev",   seed = SEED_BOOT_EV),
  list(pair = "loss_avoid_adj~b_loss", agn = "loss_avoid_adj",
       mod = "b_loss", seed = SEED_BOOT_LOSS)
)

t0 <- Sys.time()
agreement <- purrr::map_dfr(pairs, function(p) {
  cross_half_point(
    A1 = cell$est[, , 1, p$agn], A2 = cell$est[, , 2, p$agn],
    M1 = cell$est[, , 1, p$mod], M2 = cell$est[, , 2, p$mod],
    n_boot = N_BOOT, seed = p$seed
  ) |>
    mutate(pair = p$pair, agnostic = p$agn, model_based = p$mod,
           cell_id = spec$cell_id, regime = spec$regime,
           k = spec$k, r = spec$r, n_splits = N_SPLITS,
           seed_cell = SEED_CELL, seed_boot = p$seed,
           status = "VERIFIED", .before = 1)
})
cat("bootstrap wall time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
readr::write_csv(agreement, file.path(res_dir, "cross_half_agreement.csv"))
knitr::kable(agreement |>
               select(pair, measure, rho, ci_lo, ci_hi, rho_dir1, rho_dir2,
                      n_subjects, n_split_dirs_used, n_split_dirs),
             digits = 3,
             caption = "Cross-half vs same-half Spearman agreement (raw; 95% subject-cluster bootstrap CIs; no SB correction, no disattenuation)")

# ---- ceiling ----
step2 <- readr::read_csv(here("output", "results", "step2",
                              "reliability_by_cell_construct.csv"),
                         show_col_types = FALSE) |>
  filter(axis == "i_exp", k == 7)

ceiling_tbl <- purrr::map_dfr(pairs, function(p) {
  r1 <- step2 |> filter(construct == p$agn)
  r2 <- step2 |> filter(construct == p$mod)
  cross <- agreement |> filter(pair == p$pair, measure == "cross_half")
  tibble(
    pair = p$pair,
    rel_agn_sb = r1$r_sb, rel_mod_sb = r2$r_sb,
    ceiling_sb = ifelse(r1$r_sb > 0 & r2$r_sb > 0,
                        sqrt(r1$r_sb * r2$r_sb), NA_real_),
    rel_agn_raw = r1$r_raw, rel_mod_raw = r2$r_raw,
    ceiling_raw = ifelse(r1$r_raw > 0 & r2$r_raw > 0,
                         sqrt(r1$r_raw * r2$r_raw), NA_real_),
    observed_cross_half = cross$rho,
    within_sb_ceiling = cross$rho <= ceiling_sb
  )
})
readr::write_csv(ceiling_tbl, file.path(res_dir, "ceiling_check.csv"))
knitr::kable(ceiling_tbl, digits = 3,
             caption = "Observed cross-half agreement vs reliability-implied ceilings (Step 2, k = 7 experience)")

# ---- interp ----
g <- function(pr, ms) agreement |> filter(pair == pr, measure == ms)
ev_c <- g("ev_sens_agn~b_ev", "cross_half")
ev_s <- g("ev_sens_agn~b_ev", "same_half")
lo_c <- g("loss_avoid_adj~b_loss", "cross_half")
lo_s <- g("loss_avoid_adj~b_loss", "same_half")
cat(sprintf(
"- **EV pair:** same-half agreement %.3f [%.3f, %.3f] collapses to **%.3f [%.3f, %.3f] cross-half** — the gap is the shared-sampling-noise component of Step 1's same-trial agreement (.787 on the full budget).
- **Loss pair:** same-half %.3f [%.3f, %.3f] vs cross-half %.3f [%.3f, %.3f].
- Values are raw Spearman between 7-problem half measures; no SB correction or disattenuation is applied or implied.\n",
  ev_s$rho, ev_s$ci_lo, ev_s$ci_hi, ev_c$rho, ev_c$ci_lo, ev_c$ci_hi,
  lo_s$rho, lo_s$ci_lo, lo_s$ci_hi, lo_c$rho, lo_c$ci_lo, lo_c$ci_hi))

# ---- repro ----
jsonlite::write_json(
  list(seed_cell = SEED_CELL, seed_boot_ev = SEED_BOOT_EV,
       seed_boot_loss = SEED_BOOT_LOSS, n_splits = N_SPLITS, n_boot = N_BOOT,
       cell_id = spec$cell_id),
  file.path(res_dir, "seeds.json"), auto_unbox = TRUE)

# ---- appendix ----
for (f in c("R/cross_half_agreement.R", "R/step2_reliability.R")) {
  cat("\n### ", f, "\n\n```r\n", paste(readLines(here(f)), collapse = "\n"),
      "\n```\n", sep = "")
}

# Step 10: expanded problem pool and parameter recovery on it (simulation)
# Reads:   data/raw/All estimation raw data.csv; output/results/step4/reconstruction_validation.csv (step 09)
# Writes:  output/results/step4/ (expanded_pool_summary.csv, expanded_recovery_verdict.csv, realized_grid.csv,
#          eligibility, recovery fits, seeds.json)
# Paper:   Methods (Problem pools; Validation gates); Results (Reliability differs across problem pools)
# Runtime: a few seconds

# Expanded problem pool (all Amb == 0, LotNum unrestricted) and the RECOVERY
# RE-GATE on the expanded design.
#
# Runs ONLY after the step 09 validation gate passed (asserted below by
# re-checking output/results/step4/reconstruction_validation.csv).
#
# Expanded pool: eligibility table per k (step 05 convention: a subject is
#     eligible at k if they have >= 2k expanded problems), Corr != 0 problem
#     count (kept; outcome correlation affects joint payoff draws, never the
#     marginal distributions or EVs used here), realized experience-regime grid
#     (r = 20; k from 2 up to the largest k with n_subjects >= 100).
# Recovery re-gate: SIMULATION (labeled as such everywhere): 500 simulated
#     subjects with known planted parameters on REAL donor subjects' FULL
#     expanded problem sets (step 01 machinery), 25 reps/problem, refit with
#     the step 04 Firth estimator on features computed from the validated
#     reconstruction.
#     PRE-DECLARED BAR: Spearman(true, recovered) >= .7 per parameter at the
#     full expanded budget. Parameters that fail keep their extended curves
#     but are labeled "recovery-failed on expanded design — not interpreted".

library(tidyverse)
library(here)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "fetch_official_cpc18.R"))
source(here("R", "lotshape_reconstruction.R"))
source(here("R", "simulate_recovery.R"))

res_dir <- here("output", "results", "step4")
SEEDS <- list(recovery = 191001)
jsonlite::write_json(SEEDS, file.path(res_dir, "seeds.json"), auto_unbox = TRUE)

# ---- gate check: the step 09 validation must have passed --------------------
gate <- read_csv(file.path(res_dir, "reconstruction_validation.csv"),
                 show_col_types = FALSE)
stopifnot(all(gate$pass))

# ---- verified inputs ---------------------------------------------------------
fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))
raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))

param_cols <- c("Ha", "pHa", "La", "LotShapeA", "LotNumA",
                "Hb", "pHb", "Lb", "LotShapeB", "LotNumB", "Amb", "Corr")
problems_all <- raw |>
  distinct(GameID, across(all_of(param_cols))) |>
  arrange(GameID)
stopifnot(nrow(problems_all) == 210)

# ---- expanded pool -----------------------------------------------------------
problems_exp <- problems_all |> filter(Amb == 0)
dists <- reconstruct_all_problems(problems_exp)
feat_exp <- expanded_problem_features(problems_exp, dists)
scales <- attr(feat_exp, "scales")

pool_summary <- tibble(
  n_problems_total = nrow(problems_all),
  n_problems_expanded = nrow(problems_exp),
  n_problems_ambiguous_excluded = nrow(problems_all) - nrow(problems_exp),
  n_problems_clean_step2 = sum(problems_exp$LotNumA == 1 & problems_exp$LotNumB == 1),
  n_problems_multi_outcome = sum(problems_exp$LotNumA > 1 | problems_exp$LotNumB > 1),
  n_problems_corr_nonzero = sum(problems_exp$Corr != 0),
  corr_note = paste("Corr != 0 problems KEPT: outcome correlation affects the",
                    "joint draw of (Apay, Bpay) under feedback, not the marginal",
                    "distribution or EV of either option; all features here are",
                    "marginal-distribution functionals."),
  n_loss_informative = sum(feat_exp$loss_informative),
  scale_d_ev = scales[["d_ev"]],
  scale_d_sd = scales[["d_sd"]],
  scale_d_eloss = scales[["d_eloss"]]
)
write_csv(pool_summary, file.path(res_dir, "expanded_pool_summary.csv"))
print(pool_summary |> select(-corr_note), width = Inf)

# per-subject expanded problem sets (subjects face 30 games; Amb==0 subset)
subj_sets <- raw |>
  semi_join(problems_exp, by = "GameID") |>
  distinct(SubjID, GameID) |>
  group_by(SubjID) |>
  summarise(games = list(GameID), n_expanded = n(), .groups = "drop")

eligibility <- tibble(k = 2:15) |>
  mutate(problems_needed = 2 * k,
         n_subjects = map_int(k, ~ sum(subj_sets$n_expanded >= 2 * .x)))
write_csv(eligibility, file.path(res_dir, "expanded_pool_eligibility.csv"))
print(eligibility, n = Inf)

k_max <- max(eligibility$k[eligibility$n_subjects >= 100])
realized_grid <- tibble(
  regime = "exp", r = 20, k = 2:k_max,
  n_subjects = eligibility$n_subjects[match(2:k_max, eligibility$k)],
  rule = "k from 2 up to the maximum with n_subjects >= 100 (pre-declared)"
)
write_csv(realized_grid, file.path(res_dir, "realized_grid.csv"))
cat(sprintf("\nRealized grid: experience regime, r = 20, k = 2..%d\n", k_max))
cat(sprintf("Per-subject expanded problems: median %g, min %d, max %d\n",
            median(subj_sets$n_expanded), min(subj_sets$n_expanded),
            max(subj_sets$n_expanded)))

# ---- recovery re-gate (SIMULATION; labeled as such) --------------------------
# Each of 500 simulated subjects inherits the FULL expanded problem set of a
# randomly drawn real donor subject (full expanded budget), 25 reps/problem.
set.seed(SEEDS$recovery)
true_params <- draw_true_params(500)
donor_sets <- setNames(subj_sets$games, subj_sets$SubjID)

fits <- map_dfr(seq_len(nrow(true_params)), function(i) {
  tp <- true_params[i, ]
  donor <- sample(names(donor_sets), 1)
  f <- feat_exp[match(donor_sets[[donor]], feat_exp$GameID), ]
  bind_cols(
    tp, fit_subject_firth(simulate_choices(tp, f, n_reps = 25)),
    tibble(k_problems = nrow(f), n_reps = 25, donor_subj = donor,
           n_loss_informative = sum(f$loss_informative))
  )
})
saveRDS(fits, file.path(res_dir, "expanded_recovery_fits_raw.rds"))

flags <- fits |>
  summarise(n_sim = n(), n_converged = sum(converged),
            n_degenerate = sum(degenerate), n_extreme = sum(extreme),
            median_k = median(k_problems), min_k = min(k_problems),
            max_k = max(k_problems)) |>
  mutate(label = "SIMULATION (recovery re-gate, expanded design, Firth estimator)")
write_csv(flags, file.path(res_dir, "expanded_recovery_fit_flags.csv"))
print(flags, width = Inf)

set.seed(SEEDS$recovery + 1)
recov <- map_dfr(c(b_ev = "b_ev", b_risk = "b_risk", b_loss = "b_loss"),
                 function(p) {
  boot_cor(fits[[p]], fits[[paste0("est_", p)]], method = "spearman",
           n_boot = 2000)
}, .id = "param") |>
  mutate(design = "expanded pool, full budget (donor's complete Amb==0 set), 25 reps",
         estimator = "Firth (brglm2, AS_mean)",
         label = "SIMULATION")
write_csv(recov, file.path(res_dir, "expanded_recovery_by_param.csv"))
print(recov, width = Inf)

verdict <- recov |>
  transmute(param, spearman_r = r, ci_lo, ci_hi,
            pass = r >= 0.7,
            criterion = "Spearman(true, recovered) >= 0.7 at full expanded budget (pre-declared)",
            consequence = ifelse(r >= 0.7, "interpreted",
                                 "recovery-failed on expanded design — not interpreted"))
write_csv(verdict, file.path(res_dir, "expanded_recovery_verdict.csv"))
print(verdict, width = Inf)
cat("\nRecovery re-gate finished. Failed parameters (if any) keep their curves",
    "\nbut are labeled 'recovery-failed on expanded design — not interpreted'.\n")

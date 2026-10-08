# Step 19: recovery check of the variance decomposition (simulation) and final ceiling values
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e0/ (step 18); steps/18_variance_ceiling.R (functions)
# Writes:  output/elevation/results/e0/ (variance_ceiling_final.csv, recovery_*.csv, recovery_config.json)
# Paper:   Study 2 (Locating the ceiling); Abstract; Discussion
# Runtime: about 2 minutes with 12 parallel workers

# E0 Part B.4: simulation recovery check for the variance-ceiling estimator.
# Simulated subjects use the real clean-problem features and several planted
# true-variance levels, including zero. The recovery tables are SIMULATION
# results. The final ceiling table combines these verdicts with the empirical
# estimates from step 18.
#
# The pre-declared checks require sigma2_true's bootstrap CI to cover each
# planted level and to include 0 for the zero-variance population. The latter
# checks that trial/item noise does not produce apparent trait variance.
# If recovery fails, the construct's empirical ceiling is marked NOT TRUSTWORTHY.
#
# The B.1-B.2 decomposition helpers are sourced from step 18 (E0_LIB_ONLY),
# using the same code as the empirical analysis. The construct estimator is
# estimate_half() from R/step2_reliability.R.

E0_LIB_ONLY <- TRUE
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})
source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))
source(here("R", "simulate_recovery.R"))        # analytic_p_b
source(here("steps", "18_variance_ceiling.R"))  # helpers

res_dir <- here("output", "elevation", "results", "e0")
log_dir <- here("output", "elevation", "logs")
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv"))); df
}

# ---- pre-declared SIMULATION constants -------------------------------------
N_SIM      <- 400L          # simulated subjects per level
K_FIXED    <- 14L           # clean problems/subject (= real median k_ref)
R_REPS     <- 25L           # trials/problem simulated; experience 6:25 used
B_INNER    <- 250L          # within-subject problem bootstrap
N_OUT      <- 1000L         # subject-cluster bootstrap for CIs
M_BENCH    <- 120L          # realizations per subject for the MC benchmark
TAUS       <- c(0, 0.35, 0.70, 1.05)   # between-subject SD scale (0 = zero-var)
CENTER     <- c(b0 = 0, b_ev = 1.5, b_risk = -0.2, b_loss = 0.2)
BASE_SD    <- c(b0 = 0.3, b_ev = 1.0, b_risk = 1.0, b_loss = 1.0)
SEED_BASE  <- 900000L

# ---- real features + fixed loss baseline (normalization reused in sim) ------
raw_path <- here("data", "raw", "All estimation raw data.csv")
stopifnot(identical(unname(tools::md5sum(raw_path)),
                    "38a26407a36c690be969f9d9d3fe3e5a"))
raw <- load_cpc18(raw_path)
features <- add_scaled_features(clean_problem_features(raw))
trials <- raw |>
  dplyr::filter(is_clean_row(raw)) |>
  dplyr::select(SubjID, GameID, Trial, B) |>
  dplyr::left_join(features |>
                     dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev,
                                   loss_informative, has_loss_b, ev_a, ev_b),
                   by = "GameID")
baselines <- fit_loss_baselines(trials)   # fixed exp/desc normalization curves

# One clean feature table (the pool a simulated subject draws problems from),
# with the loss baseline p-hats attached (as build_subject_data would).
pool_feat <- features |>
  dplyr::mutate(
    d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a),
    p_base_exp = ifelse(loss_informative,
                        stats::plogis(baselines$exp[1] + baselines$exp[2] * d_ev_noloss), NA),
    p_base_desc = ifelse(loss_informative,
                         stats::plogis(baselines$desc[1] + baselines$desc[2] * d_ev_noloss), NA))
pool_ids <- pool_feat$GameID
stopifnot(length(pool_ids) == 105L)

# ---- build one simulated subject's `sd` object -----------------------------
# theta: named c(b0,b_ev,b_risk,b_loss); games: chosen GameIDs; returns
# list(feat, B_mat[k x R_REPS]) consumable by estimate_half().
build_sim_subject <- function(theta, games, pool_feat, r_reps) {
  f <- pool_feat[match(games, pool_feat$GameID), ]
  p <- as.numeric(stats::plogis(theta[["b0"]] + theta[["b_ev"]] * f$z_ev +
                                  theta[["b_risk"]] * f$z_sd +
                                  theta[["b_loss"]] * f$z_eloss))
  k <- length(games)
  B_mat <- matrix(stats::rbinom(k * r_reps, 1L, rep(p, times = r_reps)),
                  nrow = k, ncol = r_reps)   # column-major: col j ~ Bernoulli(p_i)
  list(feat = f, B_mat = B_mat)
}

# ---- per-subject worker: observed point+var AND MC benchmark ---------------
sim_subject_worker <- function(theta, seed, pool_feat, pool_ids, k_fixed,
                               r_reps, b_inner, m_bench) {
  set.seed(seed)
  # observed realization -> point + within-subject problem-bootstrap var
  games_obs <- sample(pool_ids, k_fixed)
  sd_obs <- build_sim_subject(theta, games_obs, pool_feat, r_reps)
  pv <- compute_subject_point_var(sd_obs, b_inner, seed + 1L)
  # MC benchmark: fresh problem draws + fresh reps at the same design -> E[point]
  bench <- matrix(NA_real_, m_bench, length(E0_CONSTRUCTS),
                  dimnames = list(NULL, E0_CONSTRUCTS))
  for (m in seq_len(m_bench)) {
    g <- sample(pool_ids, k_fixed)
    sdm <- build_sim_subject(theta, g, pool_feat, r_reps)
    bench[m, ] <- estimate_half(sdm, seq_len(k_fixed), 6:25, "exp")[E0_CONSTRUCTS]
  }
  mu_hat <- apply(bench, 2, function(x) mean(x[is.finite(x)]))
  mc_var <- apply(bench, 2, function(x) {
    xf <- x[is.finite(x)]; if (length(xf) < 2) NA_real_ else stats::var(xf) })
  m_valid <- apply(bench, 2, function(x) sum(is.finite(x)))
  list(point = pv$point[E0_CONSTRUCTS], samp_var = pv$samp_var[E0_CONSTRUCTS],
       mu_hat = mu_hat, mc_var = mc_var, m_valid = m_valid)
}

# ---- run one planted level -------------------------------------------------
run_level <- function(tau, level_seed) {
  set.seed(level_seed)
  # draw per-subject generative parameters (tau scales between-subject SD)
  theta_df <- tibble::tibble(
    b0     = stats::rnorm(N_SIM, CENTER[["b0"]],     tau * BASE_SD[["b0"]]),
    b_ev   = stats::rnorm(N_SIM, CENTER[["b_ev"]],   tau * BASE_SD[["b_ev"]]),
    b_risk = stats::rnorm(N_SIM, CENTER[["b_risk"]], tau * BASE_SD[["b_risk"]]),
    b_loss = stats::rnorm(N_SIM, CENTER[["b_loss"]], tau * BASE_SD[["b_loss"]]))
  thetas_local <- lapply(seq_len(N_SIM), function(i) {
    c(b0 = theta_df$b0[i], b_ev = theta_df$b_ev[i],
      b_risk = theta_df$b_risk[i], b_loss = theta_df$b_loss[i]) })
  seeds_local <- level_seed + 10L + 3L * seq_len(N_SIM)

  n_workers <- max(1, min(12, detectCores() - 2))
  cl <- makeCluster(n_workers)
  invisible(clusterCall(cl, function(root, lp) {
    .libPaths(lp)
    setwd(root)                    # so here() inside sourced scripts resolves
    suppressPackageStartupMessages(library(here))
    E0_LIB_ONLY <<- TRUE
    source(file.path(root, "R", "step2_reliability.R"))
    source(file.path(root, "steps", "18_variance_ceiling.R"))
    NULL
  }, here::here(), .libPaths()))
  clusterExport(cl, c("build_sim_subject", "sim_subject_worker", "pool_feat",
                      "pool_ids", "K_FIXED", "R_REPS", "B_INNER", "M_BENCH",
                      "E0_CONSTRUCTS"), envir = globalenv())
  clusterExport(cl, c("thetas_local", "seeds_local"), envir = environment())
  res <- parLapplyLB(cl, seq_len(N_SIM), function(i) {
    sim_subject_worker(thetas_local[[i]], seeds_local[i], pool_feat, pool_ids,
                       K_FIXED, R_REPS, B_INNER, M_BENCH)
  })
  stopCluster(cl)

  points_mat <- t(vapply(res, function(x) x$point, numeric(length(E0_CONSTRUCTS))))
  vars_mat   <- t(vapply(res, function(x) x$samp_var, numeric(length(E0_CONSTRUCTS))))
  mu_mat     <- t(vapply(res, function(x) x$mu_hat, numeric(length(E0_CONSTRUCTS))))
  mcv_mat    <- t(vapply(res, function(x) x$mc_var, numeric(length(E0_CONSTRUCTS))))
  colnames(points_mat) <- colnames(vars_mat) <- colnames(mu_mat) <-
    colnames(mcv_mat) <- E0_CONSTRUCTS

  dec <- cluster_bootstrap_ceiling(points_mat, vars_mat, k_ref = K_FIXED,
                                   N_out = N_OUT, seed = level_seed + 5L)

  # planted true variance benchmarks:
  #  - MC: Var_s(mu_s) with mu_s = E[point at this design], bias-corrected for
  #    finite-M MC noise: Var(mu_hat) - mean(mc_var / M_bench).
  #  - generative: Var_s of the drawn coefficient (exact by construction).
  gen_var <- c(ev_sens_agn = NA, loss_avoid_adj = NA,
               b_ev = stats::var(theta_df$b_ev),
               b_risk = stats::var(theta_df$b_risk),
               b_loss = stats::var(theta_df$b_loss))
  planted_mc <- vapply(E0_CONSTRUCTS, function(cn) {
    mu <- mu_mat[, cn]; mcv <- mcv_mat[, cn]
    ok <- is.finite(mu)
    max(0, stats::var(mu[ok]) - mean(mcv[ok] / M_BENCH, na.rm = TRUE))
  }, numeric(1))

  out <- dec |>
    dplyr::transmute(
      construct, sim_tau = tau,
      planted_mc = planted_mc[construct],
      planted_generative = gen_var[construct],
      est_sigma2_obs = sigma2_obs, est_sigma2_err = sigma2_err,
      est_sigma2_true = sigma2_true,
      est_true_lo = sigma2_true_lo, est_true_hi = sigma2_true_hi,
      est_diff_lo = sigma2_diff_lo,
      ceiling_fulldata,
      optimization_target,
      covered_mc = (sigma2_true_lo <= planted_mc[construct] &
                      planted_mc[construct] <= sigma2_true_hi),
      ci_includes_0 = (sigma2_diff_lo <= 0))
  out
}

# ---- execute all levels ----------------------------------------------------
cat("=== E0 recovery gate (SIMULATION) ===\n")
cat(sprintf("N_sim=%d k=%d r=%d B_inner=%d N_out=%d M_bench=%d taus=%s\n",
            N_SIM, K_FIXED, R_REPS, B_INNER, N_OUT, M_BENCH,
            paste(TAUS, collapse = ",")))
t0 <- Sys.time()
recovery <- purrr::imap_dfr(TAUS, function(tau, i) {
  cat(sprintf("-- level %d/%d: tau=%.2f --\n", i, length(TAUS), tau))
  run_level(tau, SEED_BASE + i * 1000L)
})
cat("recovery wall:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
recovery <- recovery |> dplyr::mutate(data = "SIMULATION", .before = 1)
save_csv(recovery, "recovery_variance_ceiling")

# ---- verdict per construct -------------------------------------------------
verdict <- recovery |>
  dplyr::group_by(construct) |>
  dplyr::summarise(
    zero_var_ci_includes_0 = ci_includes_0[sim_tau == 0],
    n_nonzero_levels = sum(sim_tau > 0),
    n_nonzero_covered = sum(covered_mc[sim_tau > 0]),
    .groups = "drop") |>
  dplyr::mutate(
    verdict = dplyr::case_when(
      zero_var_ci_includes_0 & n_nonzero_covered == n_nonzero_levels ~ "PASS",
      TRUE ~ "FAIL"),
    note = dplyr::case_when(
      verdict == "PASS" ~ "recovers planted variance incl. zero-variance population",
      TRUE ~ "estimator does not recover planted metric variance -> real-data ceiling NOT TRUSTWORTHY"))
save_csv(verdict, "recovery_verdict")

# ---- authoritative per-construct determination (B.3 real ceiling + B.4) -----
# Integrates the real-data ceiling with the recovery verdict into a single
# trustworthiness/optimization-target flag. Three outcomes:
#   optimization target     recovery PASS AND real sigma2_true CI excludes 0
#   honest null             recovery PASS AND real sigma2_true CI includes 0
#   NOT TRUSTWORTHY         recovery FAIL (estimator collapses -> ceiling invalid)
vc_real <- readr::read_csv(file.path(res_dir, "variance_ceiling.csv"),
                           show_col_types = FALSE)
final <- vc_real |>
  dplyr::select(construct, n_subjects, k_ref, sigma2_obs, sigma2_err,
                sigma2_true, sigma2_true_lo, sigma2_true_hi, ceiling_fulldata,
                k_0.7, k_0.8, optimization_target) |>
  dplyr::left_join(verdict |> dplyr::select(construct, recovery = verdict),
                   by = "construct") |>
  dplyr::mutate(determination = dplyr::case_when(
    recovery == "FAIL" ~ "NOT TRUSTWORTHY (recovery FAIL: quasi-separation collapses metric decomposition)",
    optimization_target ~ "OPTIMIZATION TARGET (recovery PASS; real sigma2_true CI excludes 0)",
    TRUE ~ "HONEST NULL (recovery PASS; real sigma2_true CI includes 0)"))
save_csv(final, "variance_ceiling_final")
cat("\n== authoritative determination (real ceiling x recovery) ==\n")
print(final |> dplyr::select(construct, sigma2_true, sigma2_true_lo,
                             sigma2_true_hi, ceiling_fulldata, recovery,
                             determination), width = 220)

cat("\n== recovery (SIMULATION) ==\n")
print(recovery |> dplyr::select(construct, sim_tau, planted_mc, planted_generative,
                                est_sigma2_true, est_true_lo, est_true_hi,
                                covered_mc, ci_includes_0), n = 40, width = 200)
cat("\n== verdict (SIMULATION) ==\n"); print(verdict, width = 200)
jsonlite::write_json(
  list(note = "SIMULATION", n_sim = N_SIM, k_fixed = K_FIXED, r_reps = R_REPS,
       b_inner = B_INNER, n_out = N_OUT, m_bench = M_BENCH, taus = TAUS,
       center = as.list(CENTER), base_sd = as.list(BASE_SD),
       seed_base = SEED_BASE),
  file.path(res_dir, "recovery_config.json"), auto_unbox = TRUE, pretty = TRUE)
cat("\nDONE Part B.4 (SIMULATION)\n")

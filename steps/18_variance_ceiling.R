# Step 18: variance decomposition and reliability ceiling
# Reads:   data/raw/All estimation raw data.csv; R/step2_reliability.R
# Writes:  output/elevation/results/e0/ (variance_ceiling.csv, per_subject_estimates.csv, seeds_inner.csv,
#          coverage and configuration); log in output/elevation/logs/
# Paper:   Methods (Locating the ceiling: variance components); Study 2 (Locating the ceiling)
# Runtime: about 1 minute with 12 parallel workers

# E0 Part B: ceiling on per-construct true between-person variance,
# using md5-checked CPC18 data, via method-of-moments variance
# decomposition.
# For each of the 5 constructs {ev_sens_agn, loss_avoid_adj,
# b_ev, b_risk, b_loss}, in the experience regime (Trials 6-25,
# r = 20), on clean trials:
# B.1 per-subject estimate using all data (estimate_half) plus within-
# subject sampling variance from a nonparametric bootstrap over the
# clean problems experienced by the subject (resample problems with
# replacement, recompute)
# B.2 sigma2_obs = between-person variance of points; sigma2_err =
# average within-person sampling variance; sigma2_true = max(0,
# sigma2_obs - sigma2_err); ceiling = sigma2_true/sigma2_obs.
# Reliability-at-k model (1/k error): rel(k) = sigma2_true / (sigma2_true +
# sigma2_err * (k_ref/k)), where k_ref is the median number of clean
# problems experienced by the subject and .7 and .8 are target reliability levels for
# which k_0.7 and k_0.8 are estimated problem counts.
# B.3 Subject-cluster bootstrap CIs on sigma2_true, ceiling, k_0.7,
# k_0.8. For each construct, there is a flag indicating whether or not
# to optimize for a construct with a CI on sigma2_true that does not
# include 0.
# The estimator for the constructs is the same as for the main
# reliability analysis: it is estimate_half() from R/step2_reliability.R.
# The decomposition functions used are also exported from this file.
# These functions are loaded in by step 19 and are used by the recovery
# gate to run the same B.1-B.2 analysis on generated data. By
# predefining E0_LIB_ONLY to be TRUE before the file is sourced, the
# actual data analysis code can be avoided.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})
source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))   # estimate_half, build_subject_data, ...

E0_CONSTRUCTS <- c("ev_sens_agn", "loss_avoid_adj", "b_ev", "b_risk", "b_loss")

# ---- B.1 per-subject point + within-subject sampling variance ---------------
# sd = list(feat, B_mat) as built by build_subject_data(). Experience regime,
# r = 20 (pos = 6:25). Returns 5-construct point, sampling variance, valid-boot
# counts, and design sizes. Resampling PROBLEMS with replacement (whole 25-rep
# rows carried) mirrors finite problem sampling at the fixed r = 20 regime.
compute_subject_point_var <- function(sd, B_inner, seed) {
  set.seed(seed)
  k <- nrow(sd$feat)
  pos <- 6:25
  point <- estimate_half(sd, seq_len(k), pos, "exp")
  boots <- matrix(NA_real_, B_inner, length(point),
                  dimnames = list(NULL, names(point)))
  for (b in seq_len(B_inner)) {
    rows <- sample.int(k, k, replace = TRUE)
    boots[b, ] <- estimate_half(sd, rows, pos, "exp")
  }
  samp_var <- apply(boots, 2, function(x) {
    xf <- x[is.finite(x)]
    if (length(xf) < 2) NA_real_ else stats::var(xf)
  })
  n_valid <- apply(boots, 2, function(x) sum(is.finite(x)))
  list(point = point, samp_var = samp_var, n_valid = n_valid,
       k = k, n_loss = sum(sd$feat$loss_informative))
}

# ---- B.2 decomposition for one construct ------------------------------------
# points, vars: per-subject vectors (may contain NA). k_ref: reference budget.
# Returns MoM decomposition + k needed at each target reliability.
decompose_one <- function(points, vars, k_ref, targets = c(0.7, 0.8)) {
  ok <- is.finite(points) & is.finite(vars)
  p <- points[ok]; v <- vars[ok]
  n <- length(p)
  if (n < 3) {
    return(c(n = n, s2_obs = NA, s2_err = NA, s2_diff = NA, s2_true = NA,
             ceiling = NA, setNames(rep(NA, length(targets)),
                                    paste0("k_", targets))))
  }
  s2_obs <- stats::var(p)
  s2_err <- mean(v)
  s2_diff <- s2_obs - s2_err
  s2_true <- max(0, s2_diff)
  ceiling <- if (s2_obs > 0) s2_true / s2_obs else NA_real_
  # rel(k) = s2_true / (s2_true + s2_err*(k_ref/k)) = target  ->  solve for k.
  # Asymptote (k -> Inf) is 1 when s2_true > 0; "not reached" (Inf) iff
  # s2_true == 0 (rel is 0 at every k).
  k_need <- vapply(targets, function(tg) {
    if (s2_true <= 0) return(Inf)
    s2_err * k_ref / (s2_true * (1 / tg - 1))
  }, numeric(1))
  c(n = n, s2_obs = s2_obs, s2_err = s2_err, s2_diff = s2_diff,
    s2_true = s2_true, ceiling = ceiling,
    setNames(k_need, paste0("k_", targets)))
}

# ---- B.3 subject-cluster bootstrap over the (point, var) matrices -----------
# points_mat / vars_mat: [n_subj x n_construct]. Resample subjects (rows) with
# replacement; recompute every construct's decomposition per draw. Returns point
# decomposition + percentile CIs + optimization-target flag (untruncated
# sigma2_diff CI lower bound > 0, i.e. sigma2_true CI excludes 0).
cluster_bootstrap_ceiling <- function(points_mat, vars_mat, k_ref, N_out, seed,
                                       targets = c(0.7, 0.8)) {
  constructs <- colnames(points_mat)
  n_subj <- nrow(points_mat)
  point_dec <- lapply(constructs, function(cn)
    decompose_one(points_mat[, cn], vars_mat[, cn], k_ref, targets))
  names(point_dec) <- constructs

  set.seed(seed)
  keys <- c("s2_obs", "s2_err", "s2_diff", "s2_true", "ceiling",
            paste0("k_", targets))
  boot <- array(NA_real_, c(N_out, length(constructs), length(keys)),
                dimnames = list(NULL, constructs, keys))
  for (bidx in seq_len(N_out)) {
    idx <- sample.int(n_subj, replace = TRUE)
    for (cn in constructs) {
      d <- decompose_one(points_mat[idx, cn], vars_mat[idx, cn], k_ref, targets)
      boot[bidx, cn, ] <- d[keys]
    }
  }
  qci <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  rows <- lapply(constructs, function(cn) {
    d <- point_dec[[cn]]
    ci_true <- qci(boot[, cn, "s2_true"])
    ci_ceil <- qci(boot[, cn, "ceiling"])
    ci_k7 <- qci(boot[, cn, "k_0.7"])
    ci_k8 <- qci(boot[, cn, "k_0.8"])
    ci_diff <- qci(boot[, cn, "s2_diff"])
    tibble::tibble(
      construct = cn, n_subjects = d[["n"]], k_ref = k_ref,
      sigma2_obs = d[["s2_obs"]], sigma2_err = d[["s2_err"]],
      sigma2_diff = d[["s2_diff"]], sigma2_true = d[["s2_true"]],
      sigma2_true_lo = ci_true[1], sigma2_true_hi = ci_true[2],
      sigma2_diff_lo = ci_diff[1], sigma2_diff_hi = ci_diff[2],
      ceiling_fulldata = d[["ceiling"]],
      ceiling_lo = ci_ceil[1], ceiling_hi = ci_ceil[2],
      k_0.7 = d[["k_0.7"]], k_0.7_lo = ci_k7[1], k_0.7_hi = ci_k7[2],
      k_0.8 = d[["k_0.8"]], k_0.8_lo = ci_k8[1], k_0.8_hi = ci_k8[2],
      optimization_target = ci_diff[1] > 0
    )
  })
  dplyr::bind_rows(rows)
}

# ============================================================================
# Real-data pipeline (skipped when E0_LIB_ONLY is TRUE, e.g. sourced by recovery)
# ============================================================================
if (!isTRUE(get0("E0_LIB_ONLY", ifnotfound = FALSE))) {

  res_dir <- here("output", "elevation", "results", "e0")
  log_dir <- here("output", "elevation", "logs")
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  save_csv <- function(df, name) {
    readr::write_csv(df, file.path(res_dir, paste0(name, ".csv"))); df
  }

  B_INNER <- 300L
  N_OUT <- 2000L
  SEED_INNER_BASE <- 733000L
  SEED_OUTER <- 734001L

  raw_path <- here("data", "raw", "All estimation raw data.csv")
  stopifnot(identical(unname(tools::md5sum(raw_path)),
                      "38a26407a36c690be969f9d9d3fe3e5a"))
  raw <- load_cpc18(raw_path)
  stopifnot(dplyr::n_distinct(raw$SubjID) == 686L,
            dplyr::n_distinct(raw$GameID) == 210L)

  features <- add_scaled_features(clean_problem_features(raw))
  trials <- raw |>
    dplyr::filter(is_clean_row(raw)) |>
    dplyr::select(SubjID, GameID, Trial, B) |>
    dplyr::left_join(features |>
                       dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev,
                                     loss_informative, has_loss_b, ev_a, ev_b),
                     by = "GameID")
  baselines <- fit_loss_baselines(trials)
  subj_data <- build_subject_data(trials, baselines)

  k_ref <- stats::median(vapply(subj_data, function(sd) nrow(sd$feat), integer(1)))
  cat(sprintf("N subjects = %d ; k_ref (median clean problems) = %g\n",
              length(subj_data), k_ref))

  # ---- B.1 parallel per-subject point + sampling variance -------------------
  ids <- names(subj_data)
  seeds <- SEED_INNER_BASE + seq_along(ids)
  readr::write_csv(tibble::tibble(SubjID = ids, inner_seed = seeds),
                   file.path(res_dir, "seeds_inner.csv"))
  n_workers <- max(1, min(12, detectCores() - 2))
  cat("workers:", n_workers, "  B_inner:", B_INNER, "\n")
  cl <- makeCluster(n_workers)
  invisible(clusterCall(cl, function(root) {
    source(file.path(root, "R", "step2_reliability.R")); NULL
  }, here::here()))
  clusterExport(cl, c("compute_subject_point_var", "B_INNER", "subj_data",
                      "seeds"), envir = environment())

  t0 <- Sys.time()
  per_subj <- parLapplyLB(cl, seq_along(ids), function(i) {
    compute_subject_point_var(subj_data[[i]], B_INNER, seeds[i])
  })
  stopCluster(cl)
  names(per_subj) <- ids
  cat("B.1 wall:",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

  # Assemble per-subject point / sampling-variance matrices [subject x construct]
  points_mat <- t(vapply(per_subj, function(x) x$point[E0_CONSTRUCTS],
                         numeric(length(E0_CONSTRUCTS))))
  vars_mat   <- t(vapply(per_subj, function(x) x$samp_var[E0_CONSTRUCTS],
                         numeric(length(E0_CONSTRUCTS))))
  nvalid_mat <- t(vapply(per_subj, function(x) x$n_valid[E0_CONSTRUCTS],
                         numeric(length(E0_CONSTRUCTS))))
  colnames(points_mat) <- colnames(vars_mat) <- colnames(nvalid_mat) <- E0_CONSTRUCTS

  # per_subject_estimates.csv (long): SubjID, construct, point, se, n_boot_valid
  per_subject_estimates <- purrr::map_dfr(E0_CONSTRUCTS, function(cn) {
    tibble::tibble(SubjID = ids, construct = cn,
                   point = points_mat[, cn],
                   se = sqrt(vars_mat[, cn]),
                   samp_var = vars_mat[, cn],
                   n_boot_valid = nvalid_mat[, cn])
  })
  save_csv(per_subject_estimates, "per_subject_estimates")

  # Per-construct coverage note (how many subjects contribute a finite estimate)
  coverage <- purrr::map_dfr(E0_CONSTRUCTS, function(cn) {
    ok <- is.finite(points_mat[, cn]) & is.finite(vars_mat[, cn])
    tibble::tibble(construct = cn, n_subjects_used = sum(ok),
                   n_subjects_total = length(ids),
                   n_point_na = sum(!is.finite(points_mat[, cn])),
                   n_var_na = sum(!is.finite(vars_mat[, cn])))
  })
  save_csv(coverage, "variance_ceiling_coverage")
  cat("\n== coverage ==\n"); print(coverage)

  # Quasi-separation diagnostic: the per-subject Firth coefficients (b_ev/b_risk/
  # b_loss) have a well-behaved core but a heavy tail of quasi-separated fits;
  # the within-subject problem-bootstrap amplifies this (SE blow-up), which is
  # what dominates sigma2_err for the coefficient constructs. Documented so the
  # metric-scale ceiling breakdown is transparent (consistent with the
  # rank-only reliability statistics). |coef| > 15 mirrors the fitter's
  # `extreme` flag.
  sep_diag <- purrr::map_dfr(E0_CONSTRUCTS, function(cn) {
    x <- points_mat[, cn]; se <- sqrt(vars_mat[, cn])
    tibble::tibble(
      construct = cn, point_median = stats::median(x), point_IQR = stats::IQR(x),
      point_MAD = stats::mad(x), point_min = min(x), point_max = max(x),
      n_extreme_abs_gt15 = sum(abs(x) > 15),
      pct_extreme = 100 * mean(abs(x) > 15),
      se_median = stats::median(se), se_max = max(se))
  })
  save_csv(sep_diag, "coefficient_separation_diagnostic")
  cat("\n== separation diagnostic ==\n"); print(sep_diag, width = 200)

  # ---- B.2 + B.3 decomposition with subject-cluster bootstrap CI ------------
  t0 <- Sys.time()
  variance_ceiling <- cluster_bootstrap_ceiling(points_mat, vars_mat, k_ref,
                                                N_OUT, SEED_OUTER)
  cat("B.3 bootstrap wall:",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  variance_ceiling <- variance_ceiling |>
    dplyr::mutate(
      target_flag = ifelse(optimization_target,
                           "optimization target (sigma2_true CI excludes 0)",
                           "no detectable between-person trait variance (CI includes 0)"),
      .after = optimization_target)
  save_csv(variance_ceiling, "variance_ceiling")

  jsonlite::write_json(
    list(B_inner = B_INNER, n_out = N_OUT, seed_inner_base = SEED_INNER_BASE,
         seed_outer = SEED_OUTER, k_ref = k_ref, regime = "experience (6-25)",
         r = 20, constructs = E0_CONSTRUCTS),
    file.path(res_dir, "variance_ceiling_config.json"),
    auto_unbox = TRUE, pretty = TRUE)

  cat("\n== variance ceiling (per construct) ==\n")
  print(variance_ceiling |>
          dplyr::select(construct, n_subjects, sigma2_obs, sigma2_err,
                        sigma2_true, sigma2_true_lo, sigma2_true_hi,
                        ceiling_fulldata, k_0.7, k_0.8, optimization_target),
        width = 200)
  cat("\nDONE Part B\n")
}


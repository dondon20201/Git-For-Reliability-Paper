# SIMULATION-ONLY machinery for the parameter-recovery gate (steps 01-02).
# Everything produced by these functions is simulated data. It must be labeled
# as simulation everywhere it appears and is never an empirical result.
#
# Generative model (per simulated subject; static across repetitions):
#   P(choose B | problem) =
#     plogis(b0 + b_ev * z_ev + b_risk * z_sd + b_loss * z_eloss)
# where z_* are the scaled description-based difference features of REAL clean
# CPC18 problems (R/cpc18_features.R). The estimator refits the same logistic
# model per subject with glm().

recovery_param_names <- c("b0", "b_ev", "b_risk", "b_loss")

# Analytic P(choose B) for one parameter vector on a feature table.
analytic_p_b <- function(pars, feat) {
  plogis(pars[["b0"]] + pars[["b_ev"]] * feat$z_ev +
           pars[["b_risk"]] * feat$z_sd + pars[["b_loss"]] * feat$z_eloss)
}

# Simulate n_reps Bernoulli choices per problem row of `feat`.
simulate_choices <- function(pars, feat, n_reps) {
  p <- rep(analytic_p_b(pars, feat), each = n_reps)
  out <- feat[rep(seq_len(nrow(feat)), each = n_reps),
              c("GameID", "z_ev", "z_sd", "z_eloss")]
  out$choice <- stats::rbinom(length(p), 1L, p)
  out
}

# Independent uniform draws of true parameters across plausible ranges
# (in logits per problem-set SD).
draw_true_params <- function(n) {
  tibble::tibble(
    sim_id = seq_len(n),
    b0     = stats::runif(n, -0.5, 0.5),
    b_ev   = stats::runif(n, 0.2, 2.5),
    b_risk = stats::runif(n, -1.5, 1.0),
    b_loss = stats::runif(n, -2.0, 0.5)
  )
}

# Simulator self-check: for each parameter vector, simulate many repetitions of
# every problem THROUGH simulate_choices() (the exact production path) and
# compare simulated choice frequencies to the analytic probabilities.
# z is the binomial z-score of the deviation.
validate_simulator <- function(par_grid, feat, n_reps = 10000) {
  purrr::map_dfr(seq_len(nrow(par_grid)), function(i) {
    pars <- par_grid[i, ]
    ch <- simulate_choices(pars, feat, n_reps)
    p_sim <- as.numeric(
      tapply(ch$choice, factor(ch$GameID, levels = feat$GameID), mean)
    )
    p_analytic <- analytic_p_b(pars, feat)
    tibble::tibble(
      par_set = i, GameID = feat$GameID,
      p_analytic = p_analytic, p_sim = p_sim,
      z = (p_sim - p_analytic) /
        sqrt(pmax(p_analytic * (1 - p_analytic), 1e-12) / n_reps)
    )
  })
}

# Per-subject logistic estimator with honesty flags:
#   separation — glm reported fitted probabilities numerically 0 or 1
#   degenerate — some coefficient inestimable (NA; e.g. a constant predictor)
#   extreme    — |coefficient| > max_abs_coef (per-SD logits; MLE blow-up)
fit_subject <- function(choices, max_abs_coef = 15) {
  separation <- FALSE
  fit <- withCallingHandlers(
    stats::glm(choice ~ z_ev + z_sd + z_eloss,
               family = stats::binomial(), data = choices),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
        separation <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  co <- stats::coef(fit)
  stopifnot(identical(names(co), c("(Intercept)", "z_ev", "z_sd", "z_eloss")))
  tibble::tibble(
    est_b0 = co[["(Intercept)"]], est_b_ev = co[["z_ev"]],
    est_b_risk = co[["z_sd"]], est_b_loss = co[["z_eloss"]],
    converged = fit$converged, separation = separation,
    degenerate = anyNA(co),
    extreme = any(abs(co[!is.na(co)]) > max_abs_coef)
  )
}

# One recovery cell: each simulated subject gets the clean-problem set of a
# randomly drawn REAL subject (donor), subsampled to k_problems, so problem
# sets mirror what real subjects actually encounter. Only real subjects with
# at least k_problems clean problems are eligible donors.
run_recovery_cell <- function(true_params, subj_sets, feat, k_problems, n_reps) {
  eligible <- subj_sets[vapply(subj_sets, length, integer(1)) >= k_problems]
  stopifnot(length(eligible) > 0)
  purrr::map_dfr(seq_len(nrow(true_params)), function(i) {
    tp <- true_params[i, ]
    donor <- sample(names(eligible), 1)
    games <- sample(eligible[[donor]], k_problems)
    f <- feat[match(games, feat$GameID), ]
    dplyr::bind_cols(
      tp, fit_subject(simulate_choices(tp, f, n_reps)),
      tibble::tibble(
        k_problems = k_problems, n_reps = n_reps, donor_subj = donor,
        n_loss_informative = sum(f$loss_informative)
      )
    )
  })
}

# Loss-aversion sweep: CONSTRUCTED problem sets with exactly l_loss
# loss-informative clean problems and (k_total - l_loss) other clean problems,
# sampled uniformly from the respective pools. Labeled as constructed (not
# donor-based) wherever reported.
run_loss_sweep <- function(true_params, feat, l_values, k_total, n_reps) {
  pool_loss <- feat$GameID[feat$loss_informative]
  pool_other <- feat$GameID[!feat$loss_informative]
  purrr::map_dfr(l_values, function(l_loss) {
    stopifnot(l_loss <= length(pool_loss), k_total - l_loss <= length(pool_other))
    purrr::map_dfr(seq_len(nrow(true_params)), function(i) {
      tp <- true_params[i, ]
      games <- c(sample(pool_loss, l_loss), sample(pool_other, k_total - l_loss))
      f <- feat[match(games, feat$GameID), ]
      dplyr::bind_cols(
        tp, fit_subject(simulate_choices(tp, f, n_reps)),
        tibble::tibble(l_loss = l_loss, k_total = k_total, n_reps = n_reps)
      )
    })
  })
}

# --- Additions for the step 02 gate conditions -------------------------------
# The functions above are kept unchanged (used by step 01); the paired-estimator
# and redesigned-sweep machinery below extends them without altering behavior.

# Firth-penalized per-subject logistic estimator (brglm2, classic mean-bias
# reduction type = "AS_mean" = Firth 1993 / Jeffreys penalty). Same formula and
# flag semantics as fit_subject(); Firth has no separation failure mode, so
# `separation` is recorded as NA. Rank-deficient designs (e.g. an all-zero
# z_eloss column when a set has no loss-relevant variation) are caught and
# returned as degenerate fits rather than errors.
fit_subject_firth <- function(choices, max_abs_coef = 15) {
  fit <- tryCatch(
    suppressWarnings(
      stats::glm(choice ~ z_ev + z_sd + z_eloss,
                 family = stats::binomial(), data = choices,
                 method = brglm2::brglmFit, type = "AS_mean")
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble::tibble(
      est_b0 = NA_real_, est_b_ev = NA_real_,
      est_b_risk = NA_real_, est_b_loss = NA_real_,
      converged = FALSE, separation = NA, degenerate = TRUE, extreme = NA
    ))
  }
  co <- stats::coef(fit)
  stopifnot(identical(names(co), c("(Intercept)", "z_ev", "z_sd", "z_eloss")))
  tibble::tibble(
    est_b0 = co[["(Intercept)"]], est_b_ev = co[["z_ev"]],
    est_b_risk = co[["z_sd"]], est_b_loss = co[["z_eloss"]],
    converged = isTRUE(fit$converged), separation = NA,
    degenerate = anyNA(co),
    extreme = any(abs(co[!is.na(co)]) > max_abs_coef)
  )
}

# Paired recovery cell: identical simulated data per subject, fitted by every
# estimator in `fitters`, so estimator comparisons are within-dataset.
run_recovery_cell_paired <- function(true_params, subj_sets, feat, k_problems,
                                     n_reps, fitters) {
  eligible <- subj_sets[vapply(subj_sets, length, integer(1)) >= k_problems]
  stopifnot(length(eligible) > 0, length(names(fitters)) == length(fitters))
  purrr::map_dfr(seq_len(nrow(true_params)), function(i) {
    tp <- true_params[i, ]
    donor <- sample(names(eligible), 1)
    games <- sample(eligible[[donor]], k_problems)
    f <- feat[match(games, feat$GameID), ]
    ch <- simulate_choices(tp, f, n_reps)
    purrr::map_dfr(names(fitters), function(nm) {
      dplyr::bind_cols(
        tp, fitters[[nm]](ch),
        tibble::tibble(
          estimator = nm, k_problems = k_problems, n_reps = n_reps,
          donor_subj = donor, n_loss_informative = sum(f$loss_informative)
        )
      )
    })
  })
}

# Redesigned loss sweep: exactly l_loss loss-informative
# problems, filler drawn ONLY from an explicit pool (strictly no-loss clean
# problems), supporting the l_loss = 0 negative control. Estimator injectable.
run_loss_sweep2 <- function(true_params, feat, l_values, k_total, n_reps,
                            filler_ids, fit_fun) {
  pool_loss <- feat$GameID[feat$loss_informative]
  purrr::map_dfr(l_values, function(l_loss) {
    stopifnot(l_loss <= length(pool_loss),
              k_total - l_loss <= length(filler_ids))
    purrr::map_dfr(seq_len(nrow(true_params)), function(i) {
      tp <- true_params[i, ]
      games <- c(if (l_loss > 0) sample(pool_loss, l_loss),
                 sample(filler_ids, k_total - l_loss))
      f <- feat[match(games, feat$GameID), ]
      dplyr::bind_cols(
        tp, fit_fun(simulate_choices(tp, f, n_reps)),
        tibble::tibble(l_loss = l_loss, k_total = k_total, n_reps = n_reps)
      )
    })
  })
}

# Correlation with a percentile bootstrap CI (resampling simulated subjects).
# Non-finite estimate pairs are dropped and counted in n_dropped.
boot_cor <- function(x, y, method = c("spearman", "pearson"), n_boot = 2000) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  xs <- x[ok]; ys <- y[ok]
  n_used <- length(xs)
  point <- stats::cor(xs, ys, method = method)
  boots <- replicate(n_boot, {
    idx <- sample.int(n_used, replace = TRUE)
    suppressWarnings(stats::cor(xs[idx], ys[idx], method = method))
  })
  ci <- stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  tibble::tibble(
    method = method, r = point, ci_lo = ci[1], ci_hi = ci[2],
    n_used = n_used, n_dropped = sum(!ok), n_boot = n_boot
  )
}

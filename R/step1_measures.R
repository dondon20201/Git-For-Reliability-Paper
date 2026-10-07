# Per-subject bias measures on REAL clean-subset trials (used by step 04).
# Model-agnostic index definitions are fixed in advance; exact formulas are
# documented here.

# Expects `trials` = clean-subset trial rows joined with problem features:
# SubjID, GameID, B, Trial, d_ev, d_sd, loss_informative,
# near_parity_sensitivity, has_loss_a, has_loss_b, ev_a, ev_b.
#
# Indices (one row per subject appearing in `trials`):
# - ev_sens_agn: Spearman correlation between B (chose option B) and d_ev
#   across all the subject's clean trials. NA if B is constant.
# - risk_agn: among near-parity problems (|dEV| <= median |dEV| AND
#   |dSD| >= median |dSD|, problem-level), proportion of trials choosing the
#   higher-SD option. NA if the subject has no qualifying problem.
# - loss_avoid_agn: among loss-informative problems (loss in exactly one
#   option), proportion of trials choosing the no-loss option.
# - loss_avoid_adj: EV-adjusted variant — per-trial residual y - p_hat, where
#   y = chose no-loss option and p_hat comes from ONE pooled (all subjects)
#   logistic regression y ~ dEV_noloss (dEV_noloss = EV(no-loss option) -
#   EV(loss option)); the index is the subject's mean residual. The pooled
#   curve is a fixed normalization shared by all subjects, not a per-subject
#   model fit.
compute_agnostic_indices <- function(trials) {
  ev_sens <- trials |>
    dplyr::group_by(SubjID) |>
    dplyr::summarise(
      ev_sens_agn = suppressWarnings(
        stats::cor(B, d_ev, method = "spearman")),
      n_clean_problems = dplyr::n_distinct(GameID), .groups = "drop"
    )

  risk <- trials |>
    dplyr::filter(near_parity_sensitivity) |>
    dplyr::group_by(SubjID) |>
    dplyr::summarise(
      risk_agn = mean(ifelse(d_sd > 0, B, 1 - B)),
      n_risk_problems = dplyr::n_distinct(GameID), .groups = "drop"
    )

  loss_trials <- trials |>
    dplyr::filter(loss_informative) |>
    dplyr::mutate(
      chose_noloss = ifelse(has_loss_b, 1 - B, B),
      d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a)
    )
  baseline <- stats::glm(chose_noloss ~ d_ev_noloss,
                         family = stats::binomial(), data = loss_trials)
  loss_trials$p_baseline <- stats::fitted(baseline)
  loss <- loss_trials |>
    dplyr::group_by(SubjID) |>
    dplyr::summarise(
      loss_avoid_agn = mean(chose_noloss),
      loss_avoid_adj = mean(chose_noloss - p_baseline),
      n_loss_problems = dplyr::n_distinct(GameID), .groups = "drop"
    )

  out <- ev_sens |>
    dplyr::left_join(risk, by = "SubjID") |>
    dplyr::left_join(loss, by = "SubjID")
  attr(out, "loss_baseline_coefs") <- stats::coef(baseline)
  out
}

# Per-subject model-based fits on real trials (both estimators), reusing the
# recovery estimators so the family is identical to what recovery validated.
fit_real_subjects <- function(trials, fitters) {
  by_subj <- split(
    trials |> dplyr::transmute(SubjID, GameID, z_ev, z_sd, z_eloss, choice = B),
    trials$SubjID
  )
  purrr::map_dfr(names(by_subj), function(sid) {
    ch <- by_subj[[sid]]
    purrr::map_dfr(names(fitters), function(nm) {
      dplyr::bind_cols(
        tibble::tibble(SubjID = ch$SubjID[1], estimator = nm,
                       n_trials = nrow(ch),
                       n_problems = dplyr::n_distinct(ch$GameID)),
        fitters[[nm]](ch)
      )
    })
  })
}

# Mean of a per-cell statistic with a subject-cluster percentile bootstrap CI
# (resample subjects with replacement; each draw keeps all of a subject's
# cells, which respects the nesting of cells within subjects).
cluster_boot_mean <- function(values, subj, n_boot = 2000) {
  bysub <- split(values, subj)
  sums <- vapply(bysub, sum, numeric(1))
  ns <- vapply(bysub, length, numeric(1))
  point <- sum(sums) / sum(ns)
  boots <- replicate(n_boot, {
    idx <- sample.int(length(sums), replace = TRUE)
    sum(sums[idx]) / sum(ns[idx])
  })
  ci <- stats::quantile(boots, c(0.025, 0.975), names = FALSE)
  tibble::tibble(mean = point, ci_lo = ci[1], ci_hi = ci[2],
                 n_cells = length(values), n_subjects = length(sums),
                 n_boot = n_boot)
}

# Correlation between two per-cell vectors with a subject-cluster bootstrap CI.
cluster_boot_cor <- function(x, y, subj, method = "pearson", n_boot = 2000) {
  bysub_x <- split(x, subj); bysub_y <- split(y, subj)
  point <- stats::cor(x, y, method = method)
  boots <- replicate(n_boot, {
    idx <- sample.int(length(bysub_x), replace = TRUE)
    suppressWarnings(stats::cor(unlist(bysub_x[idx]), unlist(bysub_y[idx]),
                                method = method))
  })
  ci <- stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  tibble::tibble(method = method, r = point, ci_lo = ci[1], ci_hi = ci[2],
                 n_cells = length(x), n_subjects = length(bysub_x),
                 n_boot = n_boot)
}

# Mean or median of a per-subject vector with a percentile bootstrap CI.
boot_stat <- function(x, stat = c("mean", "median"), n_boot = 2000) {
  stat <- match.arg(stat)
  f <- if (stat == "mean") mean else stats::median
  xs <- x[is.finite(x)]
  point <- f(xs)
  boots <- replicate(n_boot, f(sample(xs, replace = TRUE)))
  ci <- stats::quantile(boots, c(0.025, 0.975), names = FALSE)
  tibble::tibble(stat = stat, value = point, ci_lo = ci[1], ci_hi = ci[2],
                 n_used = length(xs), n_na = sum(!is.finite(x)), n_boot = n_boot)
}

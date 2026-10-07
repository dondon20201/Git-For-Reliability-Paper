# Rank-based split-half reliability versus trial budget (used from step 05 on).
# REAL DATA ONLY. Design constants (N_SPLITS = 200, N_BOOT = 2000, budget
# grids) are pre-declared and fixed.
#
# Scheme:
# - Within each subject, clean problems are randomly partitioned into two
#   halves (problems whole: all sampled trials of a problem stay in one half).
#   Subjects with more than 2k problems are freshly subsampled to 2k each split.
# - Constructs are estimated per half; reliability at a budget point for a
#   construct = Spearman across subjects between half-1 and half-2 estimates,
#   computed per split on complete pairs, averaged over N_SPLITS splits.
# - Trial positions: within a split, one common random subset of r positions
#   from the regime's pool (experience 6-25 or description 1-5) is used for
#   all subjects and problems (a fresh draw each split; positions are ~
#   exchangeable within regime per the stability check in step 04).
# - Bootstrap OUTSIDE, splits INSIDE: each of N_BOOT subject resamples
#   recomputes every per-split Spearman on the resampled subjects and then
#   averages over splits; percentile CI over the N_BOOT averages.
# - Spearman-Brown 2r/(1+r) is reported alongside raw half-half reliability,
#   labeled as an approximation when applied to rank correlations.
#
# Identification note (structural, stated in advance): the model-based
# constructs (b_ev, b_risk, b_loss) come from a 4-parameter Firth fit, so a
# half with k < 4 problems has a rank-deficient design (only k distinct
# feature rows). Those fits return NA and the affected grid points are
# reported as unidentified rather than as reliabilities.

step2_constructs <- c("b_ev", "b_risk", "b_loss", "ev_sens_agn", "loss_avoid_adj")

# ---- per-subject data ------------------------------------------------------

# trials: clean-subset rows with SubjID, GameID, Trial, B + joined features
# (z_ev, z_sd, z_eloss, d_ev, loss_informative, has_loss_b, ev_a, ev_b).
# Returns a named list (by SubjID): feat (one row per problem, with baseline
# p-hats) and B_mat (problems x 25 choice matrix, columns = Trial 1..25).
build_subject_data <- function(trials, baselines) {
  trials <- trials |> dplyr::arrange(SubjID, GameID, Trial)
  n_per_cell <- trials |> dplyr::count(SubjID, GameID)
  stopifnot(all(n_per_cell$n == 25))

  feat_all <- trials |>
    dplyr::distinct(SubjID, GameID, z_ev, z_sd, z_eloss, d_ev,
                    loss_informative, has_loss_b, ev_a, ev_b) |>
    dplyr::mutate(
      d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a),
      p_base_exp = ifelse(
        loss_informative,
        stats::plogis(baselines$exp[1] + baselines$exp[2] * d_ev_noloss), NA),
      p_base_desc = ifelse(
        loss_informative,
        stats::plogis(baselines$desc[1] + baselines$desc[2] * d_ev_noloss), NA)
    )

  ids <- unique(trials$SubjID)
  out <- lapply(ids, function(sid) {
    tr <- trials[trials$SubjID == sid, ]
    feat <- feat_all[feat_all$SubjID == sid, ]
    B_mat <- matrix(tr$B, ncol = 25, byrow = TRUE)  # rows follow feat order
    stopifnot(nrow(B_mat) == nrow(feat))
    list(feat = feat, B_mat = B_mat)
  })
  names(out) <- as.character(ids)
  out
}

# Pooled EV baselines for loss_avoid_adj, one per regime (fixed normalization
# curves shared by all subjects; regime = trial window).
fit_loss_baselines <- function(trials) {
  lt <- trials |>
    dplyr::filter(loss_informative) |>
    dplyr::mutate(chose_noloss = ifelse(has_loss_b, 1 - B, B),
                  d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a))
  fit_window <- function(lo, hi) {
    stats::coef(stats::glm(chose_noloss ~ d_ev_noloss, family = stats::binomial(),
                           data = lt[lt$Trial >= lo & lt$Trial <= hi, ]))
  }
  list(exp = fit_window(6, 25), desc = fit_window(1, 5))
}

# ---- per-half estimation ---------------------------------------------------

# Firth fit on aggregated per-problem counts (glm.fit convention: y =
# proportion, weights = n trials). Returns 4 coefficients or NAs.
firth_half_fit <- function(X, y_prop, w) {
  fit <- tryCatch(
    suppressWarnings(
      brglm2::brglmFit(x = X, y = y_prop, weights = w,
                       family = stats::binomial(),
                       control = list(type = "AS_mean"), intercept = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || anyNA(fit$coefficients)) return(rep(NA_real_, 4))
  as.numeric(fit$coefficients)
}

# Estimate all five constructs on one half.
# sd_feat / rows: subject's feat and the half's row indices; pos: trial columns.
estimate_half <- function(sd, rows, pos, regime) {
  f <- sd$feat[rows, ]
  r <- length(pos)
  nB <- if (r == 1) sd$B_mat[rows, pos] else rowSums(sd$B_mat[rows, pos, drop = FALSE])

  # model-based (NA when k < 4: rank-deficient by construction)
  co <- if (length(rows) >= 4) {
    firth_half_fit(cbind(1, f$z_ev, f$z_sd, f$z_eloss), nB / r, rep(r, length(rows)))
  } else rep(NA_real_, 4)

  # ev_sens_agn: Spearman between per-trial 0/1 choice and d_ev
  b_vec <- as.vector(t(sd$B_mat[rows, pos, drop = FALSE]))
  ev_vec <- rep(f$d_ev, each = r)
  ev_sens <- suppressWarnings(stats::cor(b_vec, ev_vec, method = "spearman"))

  # loss_avoid_adj: mean(chose_noloss - baseline p-hat) over the half's
  # loss-informative trials
  li <- which(f$loss_informative)
  loss_adj <- if (length(li) == 0) NA_real_ else {
    counts <- ifelse(f$has_loss_b[li], r - nB[li], nB[li])
    p_hat <- if (regime == "exp") f$p_base_exp[li] else f$p_base_desc[li]
    (sum(counts) - r * sum(p_hat)) / (r * length(li))
  }

  c(b_ev = co[2], b_risk = co[3], b_loss = co[4],
    ev_sens_agn = ev_sens, loss_avoid_adj = loss_adj)
}

# ---- cell runner -----------------------------------------------------------

# spec: list(cell_id, axis, regime ("exp"|"desc"), k (problems/half),
#            r (reps/problem), lambda (loss-informative/half; NA = unstratified)
# Returns EST arrays [subject x split x half x construct] + eligibility info.
run_cell <- function(spec, subj_data, n_splits, seed) {
  set.seed(seed)
  pool <- if (spec$regime == "exp") 6:25 else 1:5
  stopifnot(spec$r <= length(pool))

  eligible <- names(subj_data)[vapply(subj_data, function(sd) {
    if (is.na(spec$lambda)) {
      nrow(sd$feat) >= 2 * spec$k
    } else {
      sum(sd$feat$loss_informative) >= 2 * spec$lambda &&
        sum(!sd$feat$loss_informative) >= 2 * (spec$k - spec$lambda)
    }
  }, logical(1))]
  n <- length(eligible)
  est <- array(NA_real_,
               dim = c(n, n_splits, 2, length(step2_constructs)),
               dimnames = list(eligible, NULL, NULL, step2_constructs))

  for (s in seq_len(n_splits)) {
    pos <- if (spec$r < length(pool)) sort(sample(pool, spec$r)) else pool
    for (i in seq_len(n)) {
      sd <- subj_data[[eligible[i]]]
      if (is.na(spec$lambda)) {
        take <- sample(seq_len(nrow(sd$feat)), 2 * spec$k)
        h1 <- take[seq_len(spec$k)]
        h2 <- take[(spec$k + 1):(2 * spec$k)]
      } else {
        li <- which(sd$feat$loss_informative)
        ot <- which(!sd$feat$loss_informative)
        li_take <- sample(li, 2 * spec$lambda)
        ot_take <- if (spec$k > spec$lambda) sample(ot, 2 * (spec$k - spec$lambda)) else integer(0)
        h1 <- c(li_take[seq_len(spec$lambda)],
                ot_take[seq_len(spec$k - spec$lambda)])
        h2 <- c(li_take[(spec$lambda + 1):(2 * spec$lambda)],
                if (spec$k > spec$lambda) ot_take[(spec$k - spec$lambda + 1):(2 * (spec$k - spec$lambda))])
      }
      est[i, s, 1, ] <- estimate_half(sd, h1, pos, spec$regime)
      est[i, s, 2, ] <- estimate_half(sd, h2, pos, spec$regime)
    }
  }
  list(spec = spec, est = est, n_eligible = n)
}

# ---- reliability metrics ---------------------------------------------------

# Per-split Spearman on complete pairs; splits with fewer than min_pairs
# complete pairs return NA and are excluded from the average (counted).
split_correlations <- function(E1, E2, min_pairs = 10) {
  vapply(seq_len(ncol(E1)), function(s) {
    x <- E1[, s]; y <- E2[, s]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_pairs) return(NA_real_)
    stats::cor(x[ok], y[ok], method = "spearman")
  }, numeric(1))
}

sb_correct <- function(r) 2 * r / (1 + r)

# Reliability for one cell-construct with subject-cluster bootstrap
# (bootstrap outside, splits inside).
reliability_point <- function(E1, E2, n_boot, seed, min_pairs = 10) {
  set.seed(seed)
  rs <- split_correlations(E1, E2, min_pairs)
  point <- mean(rs, na.rm = TRUE)
  n_subj <- nrow(E1)
  boots <- replicate(n_boot, {
    idx <- sample.int(n_subj, replace = TRUE)
    mean(split_correlations(E1[idx, , drop = FALSE],
                            E2[idx, , drop = FALSE], min_pairs), na.rm = TRUE)
  })
  ci <- stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  pair_frac <- mean(is.finite(E1) & is.finite(E2))
  tibble::tibble(
    r_raw = point, ci_lo_raw = ci[1], ci_hi_raw = ci[2],
    r_sb = sb_correct(point),
    ci_lo_sb = sb_correct(ci[1]), ci_hi_sb = sb_correct(ci[2]),
    n_subjects = n_subj,
    n_splits_used = sum(!is.na(rs)), n_splits = length(rs),
    finite_frac = pair_frac, n_boot = n_boot
  )
}

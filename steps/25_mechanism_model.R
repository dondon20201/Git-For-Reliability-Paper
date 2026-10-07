# Step 25: mechanism model: do problem features predict subset reliability (H2)
# Reads:   output/elevation/results/e2/, e3/ (steps 21-23)
# Writes:  output/elevation/results/e4/ (h2_gate.json, set_features.csv, cv_performance*.csv, robustness tables)
# Paper:   Study 2 (The information mechanism)
# Runtime: about 4 minutes, single process

# E4: mechanism model (Study 2 core) - the H2 gate. Predicts each problem
# set's split-half reliability from pre-declared set-level problem features.
#
# A MODELING step on set-level outputs of earlier steps only: the E2
# problem-information dictionary, the two E3 per-set reliability files (real
# within-form and simulation free-composition), and E3's per-set split-noise
# diagnostic used by the outcome-noise robustness. It reads NO raw data, NO
# participant choices, and NO participant split file. It recomputes NO
# reliability and fits NO choice model. Build-only is therefore STRUCTURAL, not
# filter-dependent: no participant-level record is ever loaded, so the 205
# held-out test subjects and the 137 freeze subjects cannot enter by any path.
# The finest granularity present in any input is the SET.
#
# NO optimization, NO candidate-design construction, NO freeze/test read. The
# primary model is the PRE-DECLARED low-dimensional linear one; ridge and PLS
# are ROBUSTNESS only and never gate.
#
# SIMULATION LABELLING. The sim_free105 and sim_anchor_set1 tracks consume E3
# SIMULATION reliabilities and are labeled SIMULATION in every output.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

res_dir <- here("output", "elevation", "results", "e4")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  invisible(df)
}

t_start <- Sys.time()

# ---- frozen constants (fixed before results exist) --------------------------
SEED_E4     <- 20260720L   # E1 20260718, E3 20260719, E4 20260720
N_REP       <- 200L        # repeated set-level holdout repeats per cell
TRAIN_FRAC  <- 0.70
N_BOOT      <- 2000L       # cluster (problem-level) bootstrap replicates
RIDGE_GRID  <- 10^seq(-4, 4, length.out = 41L)
PLS_NCOMP   <- 1:5
INNER_FOLDS <- 5L
TARGETS     <- c("loss_avoid_adj", "ev_sens_agn")
K_GRID      <- c(6L, 8L, 10L, 12L)
K_GATE      <- c(8L, 10L, 12L)             # the well-identified budgets
FEATURES    <- c("F1_target_predictor_variance", "F2_loss_exposure",
                 "F3_predictor_collinearity", "F4_expected_information",
                 "F5_choice_entropy")
MODELS      <- c("FULL", "BASE_count", "BASE_category", "ridge", "PLS")
GATE_SPEARMAN_MIN <- 0.30                  # G-sim threshold
GATE_BOOT_SUPPORT <- 0.90                  # G-beyond threshold
MIN_CELL_N  <- 40L                         # minimum sets for a CV cell to be evaluable

seeds_log <- tibble::tibble(purpose = character(), scope = character(),
                            seed = integer(), derivation = character())
add_seed <- function(purpose, scope, seed, derivation) {
  seeds_log <<- dplyr::bind_rows(seeds_log, tibble::tibble(
    purpose = purpose, scope = scope, seed = as.integer(seed),
    derivation = derivation))
  invisible(seed)
}

# ---- inputs: E2 + E3 set-level outputs ONLY --------------------------------
p_e2   <- here("output", "elevation", "results", "e2", "problem_information_dictionary.csv")
p_real <- here("output", "elevation", "results", "e3", "reliability_by_set_real.csv")
p_sim  <- here("output", "elevation", "results", "e3", "reliability_recovery_by_set_sim.csv")
p_spl  <- here("output", "elevation", "results", "e3", "e3_split_noise_by_set.csv")

D    <- readr::read_csv(p_e2,   show_col_types = FALSE)
REAL <- readr::read_csv(p_real, show_col_types = FALSE)
SIM  <- readr::read_csv(p_sim,  show_col_types = FALSE)
SPL  <- readr::read_csv(p_spl,  show_col_types = FALSE)

input_md5 <- vapply(c(e2_dictionary = p_e2, e3_real = p_real, e3_sim = p_sim,
                      e3_split_noise = p_spl),
                    function(f) unname(tools::md5sum(f)), character(1))

# provenance assertions
stopifnot(nrow(D) == 105L)
stopifnot(!anyDuplicated(D$GameID))
stopifnot(all(c("z_ev", "z_sd", "z_eloss", "d_ev_noloss", "eloss_a", "eloss_b",
                "loss_informative", "choice_entropy_exp", "finfo_b_ev",
                "info_loss_avoid_adj") %in% names(D)))
stopifnot(all(TARGETS %in% REAL$target), all(TARGETS %in% SIM$target))

D <- D |> dplyr::arrange(GameID)
gid_pos <- match(seq_len(max(D$GameID)), D$GameID)   # GameID -> row index of D

# ---- set inventory ---------------------------------------------------------
# track labels: the two GATING tracks are real_form1 (empirical anchor) and
# sim_free105 (primary). All other
# tracks are reported for transparency and are marked gating = FALSE. No gate
# threshold anywhere reads a non-gating track.
sets_real <- REAL |>
  dplyr::transmute(track = paste0("real_form", form),
                   kind = "REAL", form = as.integer(form), K = as.integer(K),
                   target, set_id = as.integer(set_id), game_ids,
                   n_loss_informative = as.integer(n_loss_informative),
                   n_subjects = as.integer(n_subjects), r_sb)
sets_sim <- SIM |>
  dplyr::transmute(track = dplyr::if_else(composition == "free_105",
                                          "sim_free105", "sim_anchor_set1"),
                   kind = "SIMULATION", form = NA_integer_, K = as.integer(K),
                   target, set_id = as.integer(set_id), game_ids,
                   n_loss_informative = as.integer(n_loss_informative),
                   n_subjects = NA_integer_, r_sb = sim_r_sb)
SETS <- dplyr::bind_rows(sets_real, sets_sim) |>
  dplyr::mutate(gating = track %in% c("real_form1", "sim_free105"))

GATE_TRACK_SIM  <- "sim_free105"
GATE_TRACK_REAL <- "real_form1"

# ---- pre-declared feature aggregation --------------------------------------
# Computed once per unique game_ids string (target-independent parts) and then
# specialized per target. Deterministic; no RNG.
uniq_ids <- unique(SETS$game_ids)
idx_list <- lapply(strsplit(uniq_ids, ";", fixed = TRUE),
                   function(s) gid_pos[as.integer(s)])
stopifnot(!any(vapply(idx_list, function(i) anyNA(i), logical(1))))

z_ev <- D$z_ev; z_sd <- D$z_sd; z_eloss <- D$z_eloss
li_v <- D$loss_informative
dnl  <- D$d_ev_noloss
ent  <- D$choice_entropy_exp
fev  <- D$finfo_b_ev
ila  <- D$info_loss_avoid_adj
elss <- D$eloss_a + D$eloss_b

abs_cor <- function(x, y) {
  if (length(x) < 2L) return(NA_real_)
  sx <- stats::sd(x); sy <- stats::sd(y)
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
  abs(stats::cor(x, y))
}

n_u <- length(idx_list)
base_tbl <- tibble::tibble(
  game_ids = uniq_ids,
  K_obs = vapply(idx_list, length, integer(1)),
  F2_loss_exposure = vapply(idx_list, function(i) sum(li_v[i]), numeric(1)),
  F3_predictor_collinearity = vapply(idx_list, function(i) {
    cc <- c(abs_cor(z_ev[i], z_sd[i]), abs_cor(z_ev[i], z_eloss[i]),
            abs_cor(z_sd[i], z_eloss[i]))
    if (all(is.na(cc))) NA_real_ else max(cc, na.rm = TRUE)
  }, numeric(1)),
  collinearity_mean = vapply(idx_list, function(i) {
    cc <- c(abs_cor(z_ev[i], z_sd[i]), abs_cor(z_ev[i], z_eloss[i]),
            abs_cor(z_sd[i], z_eloss[i]))
    if (all(is.na(cc))) NA_real_ else mean(cc, na.rm = TRUE)
  }, numeric(1)),
  F5_choice_entropy = vapply(idx_list, function(i) mean(ent[i]), numeric(1)),
  mean_eloss = vapply(idx_list, function(i) mean(elss[i]), numeric(1)),
  var_z_ev = vapply(idx_list, function(i) stats::var(z_ev[i]), numeric(1)),
  var_d_ev_noloss_li = vapply(idx_list, function(i) {
    v <- dnl[i][li_v[i]]; v <- v[is.finite(v)]
    if (length(v) >= 2L) stats::var(v) else NA_real_
  }, numeric(1)),
  sum_info_loss = vapply(idx_list, function(i) sum(ila[i]), numeric(1)),
  sum_finfo_b_ev = vapply(idx_list, function(i) sum(fev[i]), numeric(1))
)

FEAT <- SETS |>
  dplyr::left_join(base_tbl, by = "game_ids") |>
  dplyr::mutate(
    F1_target_predictor_variance = dplyr::if_else(
      target == "loss_avoid_adj", var_d_ev_noloss_li, var_z_ev),
    F4_expected_information = dplyr::if_else(
      target == "loss_avoid_adj", sum_info_loss, sum_finfo_b_ev),
    frac_loss = F2_loss_exposure / K_obs)

stopifnot(all(FEAT$K_obs == FEAT$K))
stopifnot(all(FEAT$F2_loss_exposure == FEAT$n_loss_informative))

save_csv(
  FEAT |> dplyr::select(track, kind, gating, form, K, target, set_id,
                        dplyr::all_of(FEATURES), collinearity_mean, frac_loss,
                        mean_eloss, n_loss_informative, n_subjects, r_sb,
                        game_ids),
  "set_features")

# ---- modeling primitives ---------------------------------------------------
# Ridge: closed form on standardized predictors, lambda by deterministic inner
# CV. PLS: SIMPLS, ncomp by the same inner CV. Both hand-implemented in base R
# so no package is added to the lockfile and both are exactly reproducible.
ridge_fit <- function(X, y, lambda) {
  p <- ncol(X)
  b <- solve(crossprod(X) + lambda * diag(p), crossprod(X, y))
  as.numeric(b)
}

simpls_fit <- function(X, y, ncomp) {
  n <- nrow(X); p <- ncol(X)
  ncomp <- min(ncomp, p, n - 1L)
  S <- crossprod(X, y)
  R <- matrix(0, p, ncomp); V <- matrix(0, p, ncomp); Q <- numeric(ncomp)
  for (a in seq_len(ncomp)) {
    if (sqrt(sum(S^2)) < 1e-12) { R <- R[, seq_len(a - 1L), drop = FALSE]
                                  Q <- Q[seq_len(a - 1L)]; break }
    r <- S
    tt <- X %*% r
    nt <- sqrt(sum(tt^2))
    if (nt < 1e-12) { R <- R[, seq_len(a - 1L), drop = FALSE]
                      Q <- Q[seq_len(a - 1L)]; break }
    tt <- tt / nt; r <- r / nt
    pl <- crossprod(X, tt)
    Q[a] <- sum(y * tt)
    v <- pl
    if (a > 1L) {
      Vp <- V[, seq_len(a - 1L), drop = FALSE]
      v <- v - Vp %*% crossprod(Vp, pl)
    }
    nv <- sqrt(sum(v^2))
    if (nv < 1e-12) { R <- R[, seq_len(a - 1L), drop = FALSE]
                      Q <- Q[seq_len(a - 1L)]; break }
    v <- v / nv
    S <- S - v %*% crossprod(v, S)
    R[, a] <- r; V[, a] <- v
  }
  if (ncol(R) == 0L) return(rep(0, p))
  as.numeric(R %*% Q)
}

inner_folds <- function(n, k, seed) {
  set.seed(seed)
  ((sample.int(n) - 1L) %% k) + 1L
}

# Rank-tolerant least squares: an aliased (collinear or constant) column gets a
# zero coefficient rather than aborting the fit, so a degenerate training fold
# degrades gracefully instead of silently voiding the cell.
safe_lsq <- function(X, y) {
  f <- stats::lm.fit(X, y)
  b <- f$coefficients; b[is.na(b)] <- 0
  as.numeric(b)
}

# Fit every model on a training sample and predict a test sample. Feature
# standardization uses TRAINING moments only (no leakage across the split).
fit_predict <- function(Xtr, ytr, Xte, seed) {
  mu <- colMeans(Xtr)
  sg <- apply(Xtr, 2, stats::sd)
  sg[!is.finite(sg) | sg == 0] <- 1
  Ztr <- sweep(sweep(Xtr, 2, mu, "-"), 2, sg, "/")
  Zte <- sweep(sweep(Xte, 2, mu, "-"), 2, sg, "/")
  ybar <- mean(ytr); yc <- ytr - ybar
  out <- list()

  # FULL: all five pre-declared features (PRIMARY)
  out$FULL <- as.numeric(cbind(1, Zte) %*% safe_lsq(cbind(1, Ztr), ytr))

  # BASE_count: within-K this is the intercept only (K is constant), so its
  # predictions are constant and its Spearman is undefined by construction.
  out$BASE_count <- rep(ybar, nrow(Zte))

  # BASE_category: loss exposure only (the broad problem category)
  j2 <- match("F2_loss_exposure", colnames(Ztr))
  out$BASE_category <- as.numeric(
    cbind(1, Zte[, j2]) %*% safe_lsq(cbind(1, Ztr[, j2]), ytr))

  # Inner CV folds, prepared once and reused by both robustness learners.
  fid <- inner_folds(nrow(Ztr), min(INNER_FOLDS, nrow(Ztr)), seed)
  p <- ncol(Ztr); Ip <- diag(p)
  prep <- lapply(unique(fid), function(f) {
    i1 <- fid != f; i2 <- !i1
    if (sum(i1) < 3L || sum(i2) < 1L) return(NULL)
    mu1 <- colMeans(Ztr[i1, , drop = FALSE])
    Z1 <- sweep(Ztr[i1, , drop = FALSE], 2, mu1, "-")
    Z2 <- sweep(Ztr[i2, , drop = FALSE], 2, mu1, "-")
    yb <- mean(ytr[i1]); y1c <- ytr[i1] - yb
    list(XtX = crossprod(Z1), Xty = crossprod(Z1, y1c), Z1 = Z1, y1c = y1c,
         Z2 = Z2, y2 = ytr[i2], yb = yb, n2 = sum(i2))
  })
  prep <- Filter(Negate(is.null), prep)

  # ridge (ROBUSTNESS only)
  mse <- vapply(RIDGE_GRID, function(lam) {
    e <- 0; m <- 0
    for (P in prep) {
      b <- solve(P$XtX + lam * Ip, P$Xty)
      e <- e + sum((P$y2 - (P$yb + as.numeric(P$Z2 %*% b)))^2); m <- m + P$n2
    }
    if (m == 0L) NA_real_ else e / m
  }, numeric(1))
  lam <- if (all(is.na(mse))) RIDGE_GRID[1] else RIDGE_GRID[which.min(mse)]
  Zc <- sweep(Ztr, 2, colMeans(Ztr), "-")
  br <- ridge_fit(Zc, yc, lam)
  out$ridge <- ybar + as.numeric(sweep(Zte, 2, colMeans(Ztr), "-") %*% br)

  # PLS (ROBUSTNESS only)
  msp <- vapply(PLS_NCOMP, function(nc) {
    e <- 0; m <- 0
    for (P in prep) {
      b <- simpls_fit(P$Z1, P$y1c, nc)
      e <- e + sum((P$y2 - (P$yb + as.numeric(P$Z2 %*% b)))^2); m <- m + P$n2
    }
    if (m == 0L) NA_real_ else e / m
  }, numeric(1))
  nc <- if (all(is.na(msp))) 1L else PLS_NCOMP[which.min(msp)]
  bp <- simpls_fit(Zc, yc, nc)
  out$PLS <- ybar + as.numeric(sweep(Zte, 2, colMeans(Ztr), "-") %*% bp)

  list(pred = out, ridge_lambda = lam, pls_ncomp = nc)
}

sp <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
}

# ---- the analysis sample for a cell ----------------------------------------
# FULL and BASE_category are fit on the IDENTICAL sample (every set with a
# finite r_sb and all five features finite) so the FULL-minus-BASE comparison
# is not confounded by differential exclusions.
cell_data <- function(track, target, K) {
  d <- FEAT |> dplyr::filter(track == !!track, target == !!target, K == !!K)
  keep <- is.finite(d$r_sb) &
    Reduce(`&`, lapply(FEATURES, function(f) is.finite(d[[f]])))
  list(all_n = nrow(d), d = d[keep, , drop = FALSE],
       n_dropped_r_sb = sum(!is.finite(d$r_sb)),
       n_dropped_feat = sum(is.finite(d$r_sb) &
                              !Reduce(`&`, lapply(FEATURES, function(f) is.finite(d[[f]])))))
}

# ---- CV: repeated set-level 70/30 holdout within K -------------------------
cv_cell <- function(d, seed) {
  n <- nrow(d)
  X <- as.matrix(d[, FEATURES]); y <- d$r_sb
  n_tr <- max(5L, floor(TRAIN_FRAC * n))
  if (n - n_tr < 3L) return(NULL)
  pred_sum <- matrix(0, n, length(MODELS), dimnames = list(NULL, MODELS))
  pred_cnt <- matrix(0L, n, length(MODELS), dimnames = list(NULL, MODELS))
  rep_sp <- matrix(NA_real_, N_REP, length(MODELS), dimnames = list(NULL, MODELS))
  rep_r2 <- matrix(NA_real_, N_REP, length(MODELS), dimnames = list(NULL, MODELS))
  lams <- numeric(N_REP); ncs <- integer(N_REP)
  for (r in seq_len(N_REP)) {
    set.seed(seed + r)
    tr <- sample.int(n, n_tr); te <- setdiff(seq_len(n), tr)
    fp <- fit_predict(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE],
                      seed = seed + 100000L + r)
    lams[r] <- fp$ridge_lambda; ncs[r] <- fp$pls_ncomp
    sst <- sum((y[te] - mean(y[tr]))^2)
    for (m in MODELS) {
      p <- fp$pred[[m]]
      pred_sum[te, m] <- pred_sum[te, m] + ifelse(is.finite(p), p, 0)
      pred_cnt[te, m] <- pred_cnt[te, m] + as.integer(is.finite(p))
      rep_sp[r, m] <- sp(p, y[te])
      rep_r2[r, m] <- if (sst > 0 && all(is.finite(p))) 1 - sum((y[te] - p)^2) / sst else NA_real_
    }
  }
  oof <- pred_sum / ifelse(pred_cnt == 0L, NA, pred_cnt)
  list(y = y, oof = oof, rep_sp = rep_sp, rep_r2 = rep_r2,
       ridge_lambda_median = stats::median(lams), pls_ncomp_median = stats::median(ncs),
       n = n, n_train = n_tr, n_test = n - n_tr)
}

# ---- cluster bootstrap over PROBLEMS ---------------------------------------
# One replicate: resample the track's problem pool with replacement; give each
# set the mean bootstrap multiplicity of the problems it contains; then draw a
# set-level sample of the same size with probability proportional to that
# weight. This propagates problem-level resampling to the overlapping sets that
# contain them, which ordinary SEs cannot do. Within a replicate the SAME
# problem draw is shared across K and across models, so the FULL-minus-BASE
# deltas and the pooled-across-K statistic are correctly paired.
set_problem_matrix <- function(d, pool) {
  ids <- strsplit(d$game_ids, ";", fixed = TRUE)
  M <- matrix(0, nrow(d), length(pool))
  lk <- stats::setNames(seq_along(pool), as.character(pool))
  for (i in seq_along(ids)) M[i, lk[ids[[i]]]] <- 1
  M / rowSums(M)
}

bootstrap_cell_group <- function(cells, pool, seed) {
  # cells: named list of K -> list(y, oof, W) ; returns per-K and pooled draws
  np <- length(pool)
  ks <- names(cells)
  out_sp <- array(NA_real_, c(N_BOOT, length(ks), length(MODELS)),
                  dimnames = list(NULL, ks, MODELS))
  set.seed(seed)
  for (b in seq_len(N_BOOT)) {
    cnt <- as.numeric(stats::rmultinom(1, np, rep(1 / np, np)))
    for (ki in seq_along(ks)) {
      cc <- cells[[ks[ki]]]
      w <- as.numeric(cc$W %*% cnt)
      if (sum(w) <= 0) next
      idx <- sample.int(length(w), length(w), replace = TRUE, prob = w)
      yb <- cc$y[idx]
      for (m in MODELS) out_sp[b, ki, m] <- sp(cc$oof[idx, m], yb)
    }
  }
  out_sp
}

ci_lo <- function(v) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, 0.025, na.rm = TRUE))
ci_hi <- function(v) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, 0.975, na.rm = TRUE))

# ---- main loop over tracks, targets, K -------------------------------------
tracks <- sort(unique(FEAT$track))
cv_rows <- list(); coef_rows <- list(); boot_store <- list(); cell_store <- list()

for (tk in tracks) {
  pool <- sort(unique(as.integer(unlist(strsplit(
    FEAT$game_ids[FEAT$track == tk], ";", fixed = TRUE)))))
  for (tg in TARGETS) {
    cells <- list()
    for (K in K_GRID) {
      cd <- cell_data(tk, tg, K)
      d <- cd$d
      if (nrow(d) < MIN_CELL_N) {
        cv_rows[[length(cv_rows) + 1L]] <- tibble::tibble(
          track = tk, target = tg, K = K, model = MODELS,
          kind = unique(FEAT$kind[FEAT$track == tk]),
          gating = tk %in% c(GATE_TRACK_SIM, GATE_TRACK_REAL),
          n_sets_total = cd$all_n, n_sets_used = nrow(d),
          n_dropped_r_sb_na = cd$n_dropped_r_sb, n_dropped_feature_na = cd$n_dropped_feat,
          n_train = NA_integer_, n_test = NA_integer_,
          oof_spearman_repeatmean = NA_real_, oof_spearman_pooled = NA_real_,
          oof_r2_repeatmean = NA_real_,
          boot_ci_lo = NA_real_, boot_ci_hi = NA_real_,
          delta_vs_base_category = NA_real_, delta_boot_ci_lo = NA_real_,
          delta_boot_ci_hi = NA_real_, delta_boot_support = NA_real_,
          delta_r2_vs_base_category = NA_real_,
          delta_repeatmean_vs_base_category = NA_real_,
          evaluable = FALSE,
          note = paste0("cell not evaluable: n_sets_used ", nrow(d),
                        " < MIN_CELL_N ", MIN_CELL_N))
        next
      }
      seed_k <- SEED_E4 + 1000L * match(tk, tracks) + 100L * match(tg, TARGETS) + K
      add_seed("cv_holdout", paste(tk, tg, K, sep = "|"), seed_k,
               "SEED_E4 + 1000*track_index + 100*target_index + K; repeat r uses seed + r")
      cv <- cv_cell(d, seed_k)
      if (is.null(cv)) next
      W <- set_problem_matrix(d, pool)
      cells[[as.character(K)]] <- list(y = cv$y, oof = cv$oof, W = W)
      cell_store[[paste(tk, tg, K)]] <- list(d = d, cv = cv)

      # full-sample standardized coefficients
      Xs <- scale(as.matrix(d[, FEATURES]))
      Xs[!is.finite(Xs)] <- 0
      ysd <- stats::sd(d$r_sb)
      bfull <- safe_lsq(cbind(1, Xs), d$r_sb)
      j2 <- match("F2_loss_exposure", FEATURES)
      bbase <- safe_lsq(cbind(1, Xs[, j2]), d$r_sb)
      lam_m <- cv$ridge_lambda_median; nc_m <- cv$pls_ncomp_median
      brid <- ridge_fit(Xs, d$r_sb - mean(d$r_sb), lam_m)
      bpls <- simpls_fit(Xs, d$r_sb - mean(d$r_sb), nc_m)
      coef_rows[[length(coef_rows) + 1L]] <- dplyr::bind_rows(
        tibble::tibble(model = "FULL", term = c("(Intercept)", FEATURES),
                       beta = bfull, beta_over_sd_y = bfull / ysd),
        tibble::tibble(model = "BASE_category", term = c("(Intercept)", "F2_loss_exposure"),
                       beta = bbase, beta_over_sd_y = bbase / ysd),
        tibble::tibble(model = "ridge", term = FEATURES, beta = brid,
                       beta_over_sd_y = brid / ysd),
        tibble::tibble(model = "PLS", term = FEATURES, beta = bpls,
                       beta_over_sd_y = bpls / ysd)) |>
        dplyr::mutate(track = tk, target = tg, K = K, scope = "within_K",
                      kind = unique(FEAT$kind[FEAT$track == tk]),
                      gating = tk %in% c(GATE_TRACK_SIM, GATE_TRACK_REAL),
                      n_sets_used = nrow(d), sd_r_sb = ysd,
                      ridge_lambda = lam_m, pls_ncomp = nc_m,
                      role = dplyr::if_else(model %in% c("ridge", "PLS"),
                                            "ROBUSTNESS", "PRIMARY"),
                      .before = 1)
    }
    if (!length(cells)) next
    seed_b <- SEED_E4 + 500000L + 1000L * match(tk, tracks) + 100L * match(tg, TARGETS)
    add_seed("cluster_bootstrap", paste(tk, tg, sep = "|"), seed_b,
             "SEED_E4 + 500000 + 1000*track_index + 100*target_index; one problem draw shared across K and models per replicate")
    B <- bootstrap_cell_group(cells, pool, seed_b)
    boot_store[[paste(tk, tg)]] <- B

    for (K in as.integer(names(cells))) {
      ki <- match(as.character(K), dimnames(B)[[2]])
      cs <- cell_store[[paste(tk, tg, K)]]
      cv <- cs$cv
      dfull <- B[, ki, "FULL"] - B[, ki, "BASE_category"]
      pooled_full <- vapply(MODELS, function(m) sp(cv$oof[, m], cv$y), numeric(1))
      cv_rows[[length(cv_rows) + 1L]] <- tibble::tibble(
        track = tk, target = tg, K = K, model = MODELS,
        kind = unique(FEAT$kind[FEAT$track == tk]),
        gating = tk %in% c(GATE_TRACK_SIM, GATE_TRACK_REAL),
        n_sets_total = cell_data(tk, tg, K)$all_n, n_sets_used = cv$n,
        n_dropped_r_sb_na = cell_data(tk, tg, K)$n_dropped_r_sb,
        n_dropped_feature_na = cell_data(tk, tg, K)$n_dropped_feat,
        n_train = cv$n_train, n_test = cv$n_test,
        oof_spearman_repeatmean = colMeans(cv$rep_sp, na.rm = TRUE)[MODELS],
        oof_spearman_pooled = pooled_full[MODELS],
        oof_r2_repeatmean = colMeans(cv$rep_r2, na.rm = TRUE)[MODELS],
        boot_ci_lo = vapply(MODELS, function(m) ci_lo(B[, ki, m]), numeric(1)),
        boot_ci_hi = vapply(MODELS, function(m) ci_hi(B[, ki, m]), numeric(1)),
        delta_vs_base_category = pooled_full[MODELS] - pooled_full["BASE_category"],
        delta_boot_ci_lo = vapply(MODELS, function(m)
          ci_lo(B[, ki, m] - B[, ki, "BASE_category"]), numeric(1)),
        delta_boot_ci_hi = vapply(MODELS, function(m)
          ci_hi(B[, ki, m] - B[, ki, "BASE_category"]), numeric(1)),
        delta_boot_support = vapply(MODELS, function(m)
          mean(B[, ki, m] - B[, ki, "BASE_category"] > 0, na.rm = TRUE), numeric(1)),
        delta_r2_vs_base_category = colMeans(cv$rep_r2, na.rm = TRUE)[MODELS] -
          colMeans(cv$rep_r2, na.rm = TRUE)["BASE_category"],
        # The pooled per-set out-of-fold prediction averages 200 leave-out
        # training means, which induces a spurious negative correlation with y
        # (visible as a large non-zero "Spearman" for the constant BASE_count
        # predictor). That artifact depresses a weak model more than a strong
        # one, so it INFLATES the pooled delta. The repeat-mean delta is free of
        # it and is uniformly the smaller, more conservative figure; both are
        # saved and G-beyond is read off both.
        delta_repeatmean_vs_base_category =
          colMeans(cv$rep_sp, na.rm = TRUE)[MODELS] -
          colMeans(cv$rep_sp, na.rm = TRUE)["BASE_category"],
        evaluable = TRUE, note = NA_character_)
    }
  }
}

CV <- dplyr::bind_rows(cv_rows) |>
  dplyr::arrange(dplyr::desc(gating), track, target, K, model)
save_csv(CV, "cv_performance")
COEF <- dplyr::bind_rows(coef_rows) |>
  dplyr::arrange(dplyr::desc(gating), track, target, K, model, term)
save_csv(COEF, "model_coefficients")

# ---- pooled-across-K gate statistics ---------------------------------------
pooled_stat <- function(tk, tg, models = MODELS) {
  B <- boot_store[[paste(tk, tg)]]
  if (is.null(B)) return(NULL)
  ks <- intersect(as.character(K_GATE), dimnames(B)[[2]])
  if (!length(ks)) return(NULL)
  obs <- vapply(models, function(m) mean(vapply(ks, function(k)
    CV$oof_spearman_pooled[CV$track == tk & CV$target == tg &
                             CV$K == as.integer(k) & CV$model == m], numeric(1))),
    numeric(1))
  obs_rm <- vapply(models, function(m) mean(vapply(ks, function(k)
    CV$oof_spearman_repeatmean[CV$track == tk & CV$target == tg &
                                 CV$K == as.integer(k) & CV$model == m], numeric(1))),
    numeric(1))
  draws <- vapply(models, function(m) rowMeans(B[, ks, m, drop = FALSE]),
                  numeric(N_BOOT))
  list(ks = ks, obs = obs, obs_repeatmean = obs_rm, draws = draws)
}

pool_rows <- list()
for (tk in tracks) for (tg in TARGETS) {
  ps <- pooled_stat(tk, tg)
  if (is.null(ps)) next
  dd <- ps$draws[, "FULL"] - ps$draws[, "BASE_category"]
  pool_rows[[length(pool_rows) + 1L]] <- tibble::tibble(
    track = tk, kind = unique(FEAT$kind[FEAT$track == tk]),
    target = tg, K_pooled = paste(ps$ks, collapse = ";"),
    model = MODELS,
    gating = tk %in% c(GATE_TRACK_SIM, GATE_TRACK_REAL),
    pooled_oof_spearman = ps$obs[MODELS],
    pooled_oof_spearman_repeatmean = ps$obs_repeatmean[MODELS],
    pooled_boot_ci_lo = vapply(MODELS, function(m) ci_lo(ps$draws[, m]), numeric(1)),
    pooled_boot_ci_hi = vapply(MODELS, function(m) ci_hi(ps$draws[, m]), numeric(1)),
    pooled_delta_vs_base = ps$obs[MODELS] - ps$obs["BASE_category"],
    pooled_delta_ci_lo = vapply(MODELS, function(m)
      ci_lo(ps$draws[, m] - ps$draws[, "BASE_category"]), numeric(1)),
    pooled_delta_ci_hi = vapply(MODELS, function(m)
      ci_hi(ps$draws[, m] - ps$draws[, "BASE_category"]), numeric(1)),
    pooled_delta_support = vapply(MODELS, function(m)
      mean(ps$draws[, m] - ps$draws[, "BASE_category"] > 0, na.rm = TRUE), numeric(1)))
}
POOL <- dplyr::bind_rows(pool_rows) |> dplyr::arrange(dplyr::desc(gating), track, target, model)
save_csv(POOL, "cv_performance_pooled_K")

# ---- robustness (not pre-declared): delete-one-problem jackknife -----------
# The pre-declared bootstrap above resamples problems and then redraws n sets
# with probability proportional to the mean bootstrap multiplicity of a set's
# problems. Because that weight has mean 1 and standard deviation of order
# 1/sqrt(K), the dominant variance term is the ordinary n-out-of-n set resample,
# so the interval it produces is close to a NON-cluster-aware one and is far too
# narrow for a design where 500 sets are built from 105 problems. That is a
# defect of the pre-declared bootstrap scheme.
#
# The pre-declared interval is NOT replaced or deleted; it is reported as run.
# This adds an honest second estimate of the same quantity: delete each problem
# in turn, recompute the statistic on the sets that do NOT contain it, and form
# the standard delete-one jackknife standard error over problems. It is
# deterministic (no RNG) and it moves the G-beyond evidence in the CONSERVATIVE
# direction, which is the direction a defect of this kind must be resolved in.
jack_rows <- list()
for (tk in c(GATE_TRACK_SIM, GATE_TRACK_REAL)) for (tg in TARGETS) {
  for (K in K_GRID) {
    cs <- cell_store[[paste(tk, tg, K)]]
    if (is.null(cs)) next
    d <- cs$d; cv <- cs$cv
    ids <- strsplit(d$game_ids, ";", fixed = TRUE)
    pool <- sort(unique(as.integer(unlist(ids))))
    contains <- lapply(pool, function(p) vapply(ids, function(s)
      as.character(p) %in% s, logical(1)))
    for (m in c("FULL", "BASE_category")) {
      th <- sp(cv$oof[, m], cv$y)
      ps_i <- vapply(seq_along(pool), function(j) {
        keep <- !contains[[j]]
        if (sum(keep) < 10L) return(NA_real_)
        sp(cv$oof[keep, m], cv$y[keep])
      }, numeric(1))
      dl_i <- vapply(seq_along(pool), function(j) {
        keep <- !contains[[j]]
        if (sum(keep) < 10L) return(NA_real_)
        sp(cv$oof[keep, "FULL"], cv$y[keep]) - sp(cv$oof[keep, "BASE_category"], cv$y[keep])
      }, numeric(1))
      ok <- is.finite(ps_i); nj <- sum(ok)
      se <- if (nj < 3L) NA_real_ else
        sqrt((nj - 1) / nj * sum((ps_i[ok] - mean(ps_i[ok]))^2))
      okd <- is.finite(dl_i); njd <- sum(okd)
      d_obs <- sp(cv$oof[, "FULL"], cv$y) - sp(cv$oof[, "BASE_category"], cv$y)
      se_d <- if (njd < 3L) NA_real_ else
        sqrt((njd - 1) / njd * sum((dl_i[okd] - mean(dl_i[okd]))^2))
      jack_rows[[length(jack_rows) + 1L]] <- tibble::tibble(
        track = tk, kind = unique(FEAT$kind[FEAT$track == tk]), target = tg,
        K = K, model = m, n_problems_deleted = nj, n_sets = nrow(d),
        statistic = th, jackknife_se = se,
        jack_ci_lo = th - 1.959963985 * se, jack_ci_hi = th + 1.959963985 * se,
        delta_vs_base = d_obs, delta_jackknife_se = se_d,
        delta_implied_support = if (is.na(se_d) || se_d == 0) NA_real_ else
          stats::pnorm(d_obs / se_d),
        bootstrap_se_predeclared = stats::sd(
          boot_store[[paste(tk, tg)]][, as.character(K), m], na.rm = TRUE))
    }
  }
}
JACK <- dplyr::bind_rows(jack_rows)
save_csv(JACK, "cluster_uncertainty_jackknife")

# ---- pooled cross-K fits ---------------------------------------------------
# BASE_count is the raw-count baseline "r_sb ~ K pooled"; FULL and
# BASE_category add K when pooling, and coefficients are reported per track,
# target, K and pooled. The within-K analysis above
# is the PRIMARY H2 test and is unaffected by this; these fits quantify how much
# the raw count alone explains, which the within-K design cannot show.
pooled_rows <- list(); pooled_cv_rows <- list()
for (tk in c(GATE_TRACK_SIM, GATE_TRACK_REAL)) for (tg in TARGETS) {
  dp <- dplyr::bind_rows(lapply(K_GRID, function(K) cell_data(tk, tg, K)$d))
  if (nrow(dp) < MIN_CELL_N) next
  Xp <- as.matrix(dp[, FEATURES]); Xp <- cbind(Xp, K = as.numeric(dp$K))
  Xs <- scale(Xp); Xs[!is.finite(Xs)] <- 0
  jK <- ncol(Xs); j2 <- match("F2_loss_exposure", colnames(Xs))
  bfull <- safe_lsq(cbind(1, Xs), dp$r_sb)
  bcnt  <- safe_lsq(cbind(1, Xs[, jK]), dp$r_sb)
  bcat  <- safe_lsq(cbind(1, Xs[, c(j2, jK)]), dp$r_sb)
  pooled_rows[[length(pooled_rows) + 1L]] <- dplyr::bind_rows(
    tibble::tibble(model = "FULL", term = c("(Intercept)", colnames(Xs)), beta = bfull),
    tibble::tibble(model = "BASE_count", term = c("(Intercept)", "K"), beta = bcnt),
    tibble::tibble(model = "BASE_category",
                   term = c("(Intercept)", "F2_loss_exposure", "K"), beta = bcat)) |>
    dplyr::mutate(track = tk, kind = unique(FEAT$kind[FEAT$track == tk]),
                  target = tg, scope = "pooled_K", n_sets_used = nrow(dp),
                  sd_r_sb = stats::sd(dp$r_sb), .before = 1)
  seed_p <- SEED_E4 + 800000L + 1000L * match(tk, tracks) + 100L * match(tg, TARGETS)
  add_seed("pooled_K_cv", paste(tk, tg, sep = "|"), seed_p,
           "SEED_E4 + 800000 + 1000*track_index + 100*target_index; repeat r uses seed + 1000*r")
  n <- nrow(dp); n_tr <- floor(TRAIN_FRAC * n)
  scr <- matrix(NA_real_, N_REP, 3, dimnames = list(NULL, c("FULL","BASE_count","BASE_category")))
  for (r in seq_len(N_REP)) {
    set.seed(seed_p + 1000L * r)
    tr <- sample.int(n, n_tr); te <- setdiff(seq_len(n), tr)
    mu <- colMeans(Xp[tr, , drop = FALSE]); sg <- apply(Xp[tr, , drop = FALSE], 2, stats::sd)
    sg[!is.finite(sg) | sg == 0] <- 1
    Ztr <- sweep(sweep(Xp[tr, , drop = FALSE], 2, mu, "-"), 2, sg, "/")
    Zte <- sweep(sweep(Xp[te, , drop = FALSE], 2, mu, "-"), 2, sg, "/")
    scr[r, "FULL"] <- sp(as.numeric(cbind(1, Zte) %*% safe_lsq(cbind(1, Ztr), dp$r_sb[tr])), dp$r_sb[te])
    scr[r, "BASE_count"] <- sp(as.numeric(cbind(1, Zte[, jK]) %*%
      safe_lsq(cbind(1, Ztr[, jK]), dp$r_sb[tr])), dp$r_sb[te])
    scr[r, "BASE_category"] <- sp(as.numeric(cbind(1, Zte[, c(j2, jK)]) %*%
      safe_lsq(cbind(1, Ztr[, c(j2, jK)]), dp$r_sb[tr])), dp$r_sb[te])
  }
  pooled_cv_rows[[length(pooled_cv_rows) + 1L]] <- tibble::tibble(
    track = tk, kind = unique(FEAT$kind[FEAT$track == tk]), target = tg,
    scope = "pooled_K", n_sets_used = n, model = colnames(scr),
    oof_spearman_repeatmean = colMeans(scr, na.rm = TRUE),
    delta_vs_base_count = colMeans(scr, na.rm = TRUE) - colMeans(scr, na.rm = TRUE)["BASE_count"],
    delta_vs_base_category = colMeans(scr, na.rm = TRUE) - colMeans(scr, na.rm = TRUE)["BASE_category"])
}
POOLED_FIT <- dplyr::bind_rows(pooled_rows)
save_csv(POOLED_FIT, "model_coefficients_pooled_K")
POOLED_FIT_CV <- dplyr::bind_rows(pooled_cv_rows)
save_csv(POOLED_FIT_CV, "cv_performance_pooled_fits")

# ---- item-disjoint robustness (feasibility-flagged) ------------------------
# Literal pre-declared version: halve the track's problem pool, train only on
# sets lying ENTIRELY in one half and test only on sets lying ENTIRELY in the
# other. A relaxed variant (not pre-declared; >= 80% of items in a half) is
# reported alongside and clearly labeled, because the literal version is
# combinatorially starved at these K.
disjoint_rows <- list()
for (tk in c(GATE_TRACK_SIM, GATE_TRACK_REAL)) {
  pool <- sort(unique(as.integer(unlist(strsplit(
    FEAT$game_ids[FEAT$track == tk], ";", fixed = TRUE)))))
  seed_d <- SEED_E4 + 900000L + match(tk, tracks)
  add_seed("item_disjoint_halving", tk, seed_d, "SEED_E4 + 900000 + track_index")
  set.seed(seed_d)
  perm <- sample(pool)
  halfA <- perm[seq_len(floor(length(perm) / 2))]
  for (tg in TARGETS) for (K in K_GRID) {
    d <- cell_data(tk, tg, K)$d
    if (!nrow(d)) next
    ids <- strsplit(d$game_ids, ";", fixed = TRUE)
    fracA <- vapply(ids, function(s) mean(as.integer(s) %in% halfA), numeric(1))
    for (variant in c("strict_disjoint", "relaxed_80pct")) {
      thr <- if (variant == "strict_disjoint") 1.0 else 0.8
      iA <- which(fracA >= thr); iB <- which((1 - fracA) >= thr)
      res <- lapply(list(c("A", "B"), c("B", "A")), function(dir) {
        tri <- if (dir[1] == "A") iA else iB
        tei <- if (dir[2] == "A") iA else iB
        if (length(tri) < 10L || length(tei) < 10L)
          return(tibble::tibble(direction = paste0(dir[1], "->", dir[2]),
                                n_train = length(tri), n_test = length(tei),
                                spearman_full = NA_real_, spearman_base = NA_real_,
                                mean_item_overlap = NA_real_, feasible = FALSE))
        X <- as.matrix(d[, FEATURES])
        fp <- fit_predict(X[tri, , drop = FALSE], d$r_sb[tri],
                          X[tei, , drop = FALSE], seed = seed_d + K)
        ov <- mean(vapply(tei, function(i) mean(vapply(tri, function(j)
          length(intersect(ids[[i]], ids[[j]])) / K, numeric(1))), numeric(1)))
        tibble::tibble(direction = paste0(dir[1], "->", dir[2]),
                       n_train = length(tri), n_test = length(tei),
                       spearman_full = sp(fp$pred$FULL, d$r_sb[tei]),
                       spearman_base = sp(fp$pred$BASE_category, d$r_sb[tei]),
                       mean_item_overlap = ov, feasible = TRUE)
      })
      disjoint_rows[[length(disjoint_rows) + 1L]] <- dplyr::bind_rows(res) |>
        dplyr::mutate(track = tk, kind = unique(FEAT$kind[FEAT$track == tk]),
                      target = tg, K = K, variant = variant,
                      pool_size = length(pool), half_size = length(halfA),
                      n_sets_in_cell = nrow(d),
                      variant_source = dplyr::if_else(
                        variant == "strict_disjoint", "E4_spec s4 (pre-declared)",
                        "Code-added relaxation, reported because the pre-declared version is starved"),
                      .before = 1)
    }
  }
}
DISJ <- dplyr::bind_rows(disjoint_rows)
save_csv(DISJ, "problem_disjoint_robustness")

# ---- outcome-noise robustness (deviation from the pre-declared outcome) ----
# The pre-declared robustness outcome is the covariance-based (subject-split)
# reliability. E3 never computed a per-set reliability on disjoint
# SUBJECT halves, and producing one would require reading raw build choices and
# recomputing reliabilities, which this step's design forbids. The feasible
# substitute from existing E3 outputs is E3's split-draw diagnostic: mean_r_sb averaged over 8
# independent problem-split seeds for 50 form-1 sets at K in {8, 12}, a
# LOWER-NOISE version of the same outcome. If the mechanism still predicts it,
# the H2 evidence is not an artifact of split-draw noise in a single r_sb draw.
# NO disattenuation is applied anywhere (project convention).
noise_rows <- list()
if (nrow(SPL)) {
  spl <- SPL |> dplyr::mutate(K = as.integer(K), set_id = as.integer(set_id))
  for (tg in TARGETS) for (K in sort(unique(spl$K))) {
    s <- spl |> dplyr::filter(target == tg, K == !!K)
    d <- FEAT |> dplyr::filter(track == GATE_TRACK_REAL, target == tg, K == !!K) |>
      dplyr::inner_join(s |> dplyr::select(set_id, mean_r_sb, var_across_split_seeds,
                                           n_finite), by = "set_id")
    d <- d[is.finite(d$mean_r_sb) &
             Reduce(`&`, lapply(FEATURES, function(f) is.finite(d[[f]]))), , drop = FALSE]
    if (nrow(d) < 20L) {
      noise_rows[[length(noise_rows) + 1L]] <- tibble::tibble(
        track = GATE_TRACK_REAL, target = tg, K = K, n_sets = nrow(d),
        outcome = "mean_r_sb_over_8_split_seeds",
        oof_spearman_full = NA_real_, oof_spearman_base = NA_real_,
        oof_spearman_full_on_single_draw = NA_real_, underpowered = TRUE)
      next
    }
    seed_n <- SEED_E4 + 700000L + 100L * match(tg, TARGETS) + K
    add_seed("outcome_noise_robustness", paste(GATE_TRACK_REAL, tg, K, sep = "|"),
             seed_n, "SEED_E4 + 700000 + 100*target_index + K")
    d2 <- d; d2$r_sb_single <- d2$r_sb; d2$r_sb <- d2$mean_r_sb
    cvA <- cv_cell(d2, seed_n)
    d3 <- d; d3$r_sb <- d3$r_sb_single <- d$r_sb
    cvB <- cv_cell(d3, seed_n)
    noise_rows[[length(noise_rows) + 1L]] <- tibble::tibble(
      track = GATE_TRACK_REAL, target = tg, K = K, n_sets = nrow(d),
      outcome = "mean_r_sb_over_8_split_seeds",
      oof_spearman_full = if (is.null(cvA)) NA_real_ else sp(cvA$oof[, "FULL"], cvA$y),
      oof_spearman_base = if (is.null(cvA)) NA_real_ else sp(cvA$oof[, "BASE_category"], cvA$y),
      oof_spearman_full_on_single_draw = if (is.null(cvB)) NA_real_ else
        sp(cvB$oof[, "FULL"], cvB$y),
      underpowered = nrow(d) < 60L)
  }
}
NOISE <- dplyr::bind_rows(noise_rows) |>
  dplyr::mutate(
    literal_spec_outcome = "covariance-based (subject-split) reliability",
    literal_spec_feasible = FALSE,
    reason = "E3 committed no per-set reliability on disjoint SUBJECT halves; producing one needs raw build choices + reliability recomputation, which E4_spec D1/D2 forbid",
    disattenuation_applied = FALSE)
save_csv(NOISE, "outcome_noise_robustness")

# ---- H2 gate ---------------------------------------------------------------
getcv <- function(tk, tg, K, m, col) {
  v <- CV[[col]][CV$track == tk & CV$target == tg & CV$K == K & CV$model == m]
  if (!length(v)) NA_real_ else v[1]
}
# Pre-declared conservative rule fixed in code before results: where the two
# point estimates of the same quantity disagree, the gate uses the SMALLER.
gate_stat <- function(tk, tg, K, m) {
  a <- getcv(tk, tg, K, m, "oof_spearman_pooled")
  b <- getcv(tk, tg, K, m, "oof_spearman_repeatmean")
  min(a, b, na.rm = FALSE)
}

# G-sim
gsim_per_K <- lapply(K_GATE, function(K) list(
  K = K,
  spearman_pooled = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "oof_spearman_pooled"),
  spearman_repeatmean = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "oof_spearman_repeatmean"),
  spearman_gate = gate_stat(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL"),
  ci_lo = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "boot_ci_lo"),
  ci_hi = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "boot_ci_hi"),
  meets_threshold = isTRUE(gate_stat(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL") >= GATE_SPEARMAN_MIN),
  ci_lower_positive = isTRUE(getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "boot_ci_lo") > 0)))
ps_sim <- pooled_stat(GATE_TRACK_SIM, "loss_avoid_adj")
gsim_pooled <- list(
  K_pooled = ps_sim$ks,
  spearman = unname(ps_sim$obs["FULL"]),
  spearman_repeatmean = unname(ps_sim$obs_repeatmean["FULL"]),
  spearman_gate = min(unname(ps_sim$obs["FULL"]), unname(ps_sim$obs_repeatmean["FULL"])),
  ci_lo = ci_lo(ps_sim$draws[, "FULL"]), ci_hi = ci_hi(ps_sim$draws[, "FULL"]))
# The pre-declared "CI lower bound > 0" clause applies wherever the 0.30
# threshold is applied, so it is conjoined at EVERY gate K as well as on the pooled
# statistic. Requiring it per K makes G-sim strictly harder, never easier.
gsim_pass <- all(vapply(gsim_per_K, function(x) x$meets_threshold, logical(1))) &&
  all(vapply(gsim_per_K, function(x) x$ci_lower_positive, logical(1))) &&
  isTRUE(gsim_pooled$spearman_gate >= GATE_SPEARMAN_MIN) && isTRUE(gsim_pooled$ci_lo > 0)

# G-beyond
gb_per_K <- lapply(K_GATE, function(K) list(
  K = K,
  delta_spearman = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_vs_base_category"),
  delta_r2 = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_r2_vs_base_category"),
  boot_support = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_boot_support"),
  delta_ci_lo = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_boot_ci_lo"),
  delta_ci_hi = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_boot_ci_hi"),
  positive_margin = isTRUE(getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_vs_base_category") > 0) ||
    isTRUE(getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_r2_vs_base_category") > 0),
  support_ok = isTRUE(getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_boot_support") >= GATE_BOOT_SUPPORT)))
gb_pooled_support <- POOL$pooled_delta_support[POOL$track == GATE_TRACK_SIM &
                                                 POOL$target == "loss_avoid_adj" & POOL$model == "FULL"]
gb_pooled_delta <- POOL$pooled_delta_vs_base[POOL$track == GATE_TRACK_SIM &
                                               POOL$target == "loss_avoid_adj" & POOL$model == "FULL"]
gbeyond_pass <- all(vapply(gb_per_K, function(x) x$positive_margin && x$support_ok, logical(1)))

# G-real: sign consistency between the real Set-1 FULL model and the sim FULL
# model. The pre-declared gate names "the features that carry the simulation
# model" but fixes no carrier rule, so both a strict all-five rule and a carrier rule
# (|standardized beta| >= 0.05 in the sim model at the same K) are computed and
# reported; the gate uses the carrier rule and the strict rule is disclosed.
CARRIER_MIN_BETA <- 0.05
sign_tbl <- list()
for (K in K_GATE) {
  bs <- COEF |> dplyr::filter(track == GATE_TRACK_SIM, target == "loss_avoid_adj",
                              K == !!K, model == "FULL", term %in% FEATURES)
  br <- COEF |> dplyr::filter(track == GATE_TRACK_REAL, target == "loss_avoid_adj",
                              K == !!K, model == "FULL", term %in% FEATURES)
  if (!nrow(bs) || !nrow(br)) next
  m <- dplyr::inner_join(bs |> dplyr::select(term, beta_sim = beta),
                         br |> dplyr::select(term, beta_real = beta), by = "term") |>
    dplyr::mutate(K = K, carrier = abs(beta_sim) >= CARRIER_MIN_BETA,
                  sign_agree = sign(beta_sim) == sign(beta_real),
                  target = "loss_avoid_adj",
                  track_sim = GATE_TRACK_SIM, kind_sim = "SIMULATION",
                  track_real = GATE_TRACK_REAL, kind_real = "REAL")
  sign_tbl[[length(sign_tbl) + 1L]] <- m
}
SIGN <- dplyr::bind_rows(sign_tbl)
save_csv(SIGN, "coefficient_sign_consistency")

# Carrier-threshold sensitivity. The pre-declared gate names "the features that
# carry the simulation model" but fixes NO numeric rule, so the verdict's dependence on
# the rule chosen in code is swept here and disclosed rather than left implicit.
# The rule is NOT changed on the basis of this sweep.
carrier_sweep <- NULL
if (nrow(SIGN)) {
  thr <- sort(unique(c(0, round(abs(SIGN$beta_sim), 10), CARRIER_MIN_BETA,
                       seq(0, 0.3, by = 0.005))))
  carrier_sweep <- dplyr::bind_rows(lapply(thr, function(t) {
    per_k <- vapply(K_GATE, function(K) {
      s <- SIGN[SIGN$K == K, ]
      cr <- abs(s$beta_sim) >= t
      sum(cr) > 0 && all(s$sign_agree[cr])
    }, logical(1))
    n_car <- vapply(K_GATE, function(K)
      sum(abs(SIGN$beta_sim[SIGN$K == K]) >= t), integer(1))
    tibble::tibble(carrier_threshold = t, g_real_sign_clause_pass = all(per_k),
                   min_n_carrier_over_gate_K = min(n_car),
                   n_gate_K_passing = sum(per_k))
  }))
  save_csv(carrier_sweep, "carrier_rule_sensitivity")
}

greal_per_K <- lapply(K_GATE, function(K) {
  s <- SIGN |> dplyr::filter(K == !!K)
  list(K = K,
       spearman = getcv(GATE_TRACK_REAL, "loss_avoid_adj", K, "FULL", "oof_spearman_pooled"),
       spearman_repeatmean = getcv(GATE_TRACK_REAL, "loss_avoid_adj", K, "FULL", "oof_spearman_repeatmean"),
       spearman_gate = gate_stat(GATE_TRACK_REAL, "loss_avoid_adj", K, "FULL"),
       ci_lo = getcv(GATE_TRACK_REAL, "loss_avoid_adj", K, "FULL", "boot_ci_lo"),
       delta_vs_base = getcv(GATE_TRACK_REAL, "loss_avoid_adj", K, "FULL", "delta_vs_base_category"),
       n_carrier = sum(s$carrier), n_carrier_sign_agree = sum(s$carrier & s$sign_agree),
       n_all_sign_agree = sum(s$sign_agree),
       positive = isTRUE(gate_stat(GATE_TRACK_REAL, "loss_avoid_adj", K, "FULL") > 0),
       beats_base = isTRUE(getcv(GATE_TRACK_REAL, "loss_avoid_adj", K, "FULL", "delta_vs_base_category") > 0),
       carrier_signs_consistent = nrow(s) > 0 && sum(s$carrier) > 0 &&
         all(s$sign_agree[s$carrier]))
})
greal_pass <- all(vapply(greal_per_K, function(x)
  x$positive && x$beats_base && x$carrier_signs_consistent, logical(1)))

h2_pass <- gsim_pass && gbeyond_pass && greal_pass

gate <- list(
  spec = "elevation/E4_spec.md",
  primary_target = "loss_avoid_adj", secondary_target_informational = "ev_sens_agn",
  gate_tracks = list(sim_primary = GATE_TRACK_SIM, real_anchor = GATE_TRACK_REAL),
  K_gate = K_GATE,
  thresholds = list(gsim_spearman_min = GATE_SPEARMAN_MIN,
                    gsim_ci_lower_must_exceed = 0,
                    gbeyond_boot_support_min = GATE_BOOT_SUPPORT,
                    greal_spearman_must_exceed = 0),
  conservative_point_estimate_rule = "where oof_spearman_pooled and oof_spearman_repeatmean disagree, the gate uses the SMALLER (fixed in code before results)",
  g_sim = list(per_K = gsim_per_K, pooled = gsim_pooled, pass = gsim_pass),
  g_beyond = list(per_K = gb_per_K, pooled_delta = unname(gb_pooled_delta),
                  pooled_support = unname(gb_pooled_support), pass = gbeyond_pass),
  g_real = list(per_K = greal_per_K, carrier_rule = paste0("|beta_sim| >= ", CARRIER_MIN_BETA),
                carrier_rule_source = "fixed in Code before results; E4_spec s5 names carrier features but fixes NO numeric rule",
                strict_all_five_rule_pass = if (!nrow(SIGN)) NA else
                  all(vapply(K_GATE, function(K) all(SIGN$sign_agree[SIGN$K == K]), logical(1))),
                carrier_threshold_flip_point = if (is.null(carrier_sweep)) NA_real_ else {
                  p <- carrier_sweep$carrier_threshold[carrier_sweep$g_real_sign_clause_pass]
                  if (!length(p)) NA_real_ else min(p) },
                carrier_sweep_note = "see carrier_rule_sensitivity.csv; the sign clause fails for every per-K carrier threshold at or below the flip point and passes above it only because the carrier set empties, so the verdict is knife-edge and is disclosed rather than adjusted",
                pass = greal_pass),
  # Robustness (not pre-declared) recorded alongside, NOT substituted for, the
  # pre-declared statistics. The pre-declared bootstrap under-propagates
  # problem-level clustering (see the jackknife section of the script); the
  # jackknife-implied G-beyond support is the conservative reading.
  robustness = list(
    g_beyond_delta_on_conservative_estimator = lapply(K_GATE, function(K) list(
      K = K,
      delta_pooled_estimator = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_vs_base_category"),
      delta_repeatmean_estimator = getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_repeatmean_vs_base_category"),
      still_positive = isTRUE(getcv(GATE_TRACK_SIM, "loss_avoid_adj", K, "FULL", "delta_repeatmean_vs_base_category") > 0))),
    g_beyond_support_jackknife = lapply(K_GATE, function(K) {
      j <- JACK[JACK$track == GATE_TRACK_SIM & JACK$target == "loss_avoid_adj" &
                  JACK$K == K & JACK$model == "FULL", ]
      list(K = K,
           implied_support = if (!nrow(j)) NA_real_ else j$delta_implied_support[1],
           meets_threshold = nrow(j) > 0 && isTRUE(j$delta_implied_support[1] >= GATE_BOOT_SUPPORT),
           jackknife_se = if (!nrow(j)) NA_real_ else j$delta_jackknife_se[1],
           predeclared_bootstrap_se = if (!nrow(j)) NA_real_ else j$bootstrap_se_predeclared[1]) }),
    g_beyond_survives_jackknife = {
      v <- vapply(K_GATE, function(K) {
        j <- JACK[JACK$track == GATE_TRACK_SIM & JACK$target == "loss_avoid_adj" &
                    JACK$K == K & JACK$model == "FULL", ]
        nrow(j) > 0 && isTRUE(j$delta_implied_support[1] >= GATE_BOOT_SUPPORT) }, logical(1))
      all(v) },
    note = "the pre-declared G-beyond verdict above is computed from the pre-declared bootstrap and is NOT amended by these figures; they are reported so the report can state plainly whether G-beyond survives a defensible cluster-uncertainty estimate"),
  h2_pass = h2_pass,
  h2_verdict = if (h2_pass) "PASS" else "FAIL",
  # H2 is PASS iff all three clauses hold; the pre-declared rule attaches
  # two DIFFERENT consequences to the two failure modes. Both are recorded here
  # verbatim in substance, together with which one this run's pattern selects,
  # so that the consequence is read off the pre-declared rule rather than chosen
  # after the fact. Neither the verdict nor any threshold depends on this field.
  consequence_branches = list(
    if_g_sim_or_g_beyond_fails = "H2 is NOT supported; the mechanism is reported EXPLORATORY and E5's feature-based optimization has no validated basis; bring the decision to the author; do not repurpose the analysis to pass",
    if_only_g_real_fails = "the mechanism holds under simulation but is NOT confirmed on real data; flag as a limitation (simulation-only mechanism), and E5 proceeds with that caveat explicit"),
  branch_selected = if (h2_pass) "none (H2 passed)" else if (!gsim_pass || !gbeyond_pass)
    "if_g_sim_or_g_beyond_fails" else "if_only_g_real_fails",
  build_only = TRUE, raw_data_read = FALSE, choices_read = FALSE,
  reliability_recomputed = FALSE, freeze_read = FALSE, test_read = FALSE,
  optimization_performed = FALSE, candidate_designs_constructed = FALSE)
jsonlite::write_json(gate, file.path(res_dir, "h2_gate.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 12, na = "null", null = "null")

# ---- config + seeds --------------------------------------------------------
save_csv(seeds_log, "e4_seeds")
cfg <- list(
  spec = "elevation/E4_spec.md",
  seed_base = SEED_E4, n_cv_repeats = N_REP, train_fraction = TRAIN_FRAC,
  n_bootstrap = N_BOOT, inner_folds = INNER_FOLDS,
  ridge_lambda_grid = range(RIDGE_GRID), ridge_grid_length = length(RIDGE_GRID),
  pls_ncomp_grid = range(PLS_NCOMP), min_cell_n = MIN_CELL_N,
  features = FEATURES, models = MODELS,
  primary_model = "FULL", robustness_models = c("ridge", "PLS"),
  targets = TARGETS, K_grid = K_GRID, K_gate = K_GATE,
  tracks_gating = c(GATE_TRACK_SIM, GATE_TRACK_REAL),
  tracks_all = tracks,
  input_files = as.list(input_md5),
  raw_data_read = FALSE, participant_choices_read = FALSE,
  e1_split_read = FALSE, reliability_recomputed = FALSE,
  build_only_mechanism = "structural: no participant-level record is loaded by this script",
  standardization = "features z-scored within the analysis sample for reported coefficients; within the TRAINING fold only for every cross-validated prediction",
  bootstrap_scheme = "problem-level multinomial resample of the track pool; set weight = mean bootstrap multiplicity of its problems; set-level resample of the same size drawn with probability proportional to that weight; one problem draw shared across K and models per replicate",
  base_count_within_K_note = "within K the count baseline is an intercept, so its predictions are constant and its Spearman is undefined by construction; oof_spearman_repeatmean is correctly NA, but oof_spearman_pooled is NOT zero for it because averaging leave-out training means induces a spurious negative correlation with the outcome. Read the BASE_count pooled column as a diagnostic of that artifact, never as a baseline.",
  inner_cv_seed_derivation = "the ridge/PLS inner folds use seed = (cv_holdout seed for the cell) + 100000 + repeat index; recorded here because e4_seeds.csv logs only the outer holdout seed",
  f3_na_convention = "F3 and collinearity_mean take the max/mean over the pairs whose correlation is DEFINED; a pair is undefined when one of z_ev, z_sd, z_eloss is constant across the set (272 of 77470 rows). F3 is provably invariant to this choice since |corr| >= 0; the reported-only collinearity_mean is not, and E4_spec s2 does not fix the convention.",
  disattenuation_applied = FALSE,
  known_defects_reported_not_silently_fixed = list(
    cluster_bootstrap_too_narrow = "the pre-declared bootstrap's set weight (mean bootstrap multiplicity of a set's problems) has SD of order 1/sqrt(K), so the n-out-of-n set resample dominates and the interval is close to a non-cluster-aware one. A delete-one-problem jackknife is reported alongside in cluster_uncertainty_jackknife.csv and is the conservative reading of G-beyond. The pre-declared statistic was NOT deleted or replaced.",
    cv_repeat_seed_stride = "the per-repeat seed stride is 1 while the per-K stride is K, so two cells whose K differ by 1 or 2 and whose n coincide reuse 198 of 200 holdout partitions. This affects the SECONDARY target only (all four gating cells have distinct n at every K, so no gate number is affected). Reported rather than corrected, because re-striding the seed after seeing results would re-roll every fold post hoc.",
    real_track_item_overlap = "the real Set-1 track draws up to 2000 sets from only 24 problems, so train and test sets overlap heavily in items; problem_disjoint_robustness.csv shows the real-track out-of-fold Spearman changing sign with the direction of the item split. The real-track figures are therefore an upper bound on transferable predictive accuracy."),
  deviations = list(
    covariance_subject_split = "E4_spec s0's covariance-based (subject-split) reliability robustness is NOT derivable from committed E3 outputs and cannot be produced without reading raw choices and recomputing reliabilities, which D1/D2 forbid. The substitute in outcome_noise_robustness.csv averages PROBLEM-split draw noise, which is a DIFFERENT noise source from the subject-cluster component the E3 review raised; it therefore does NOT discharge that recommendation, and the gap goes back to planning. Its cells also run at n = 34 to 50, below MIN_CELL_N, and are flagged underpowered.",
    carrier_rule = "E4_spec s5 G-real names carrier features but fixes no rule; |beta_sim| >= 0.05 fixed in code, strict all-five agreement also reported",
    disjoint_relaxation = "E4_spec s4's strict item-disjoint check is combinatorially starved at these K; a labeled 80% relaxation is reported alongside it, never in place of it"),
  runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
jsonlite::write_json(cfg, file.path(res_dir, "e4_config.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 12, na = "null")

cat("\nE4 complete. H2 verdict:", gate$h2_verdict,
    "| G-sim", gsim_pass, "| G-beyond", gbeyond_pass, "| G-real", greal_pass,
    "\nruntime (min):", round(cfg$runtime_minutes, 2), "\n")

# Step 32: leave-one-feature-out check of the mechanism model
# Reads:   output/elevation/results/e2/, e4/ (steps 21, 25)
# Writes:  output/elevation/results/dispb/ (b3_*.csv, dispb_part_b3.json, seeds)
# Paper:   Study 2 (Ruling out alternatives); Supplement (Study 2 robustness detail), through steps 33 and 36
# Runtime: about 2 minutes, single process

# B3: leave-one-feature-out checks of the frozen E4 model, with pre-declared
# gates G-B3a..d. Given the five-feature H2 FAIL, test whether dropping a
# feature changes the empirical signal, whether expected information (F4) is
# redundant with loss exposure (F2), and whether the negative F1 coefficient
# persists in every variant that retains F1.
#
# Inputs are four committed set-level CSVs and E4's saved seeds. No raw data,
# participant choices, or split file is loaded, and no reliability is
# recomputed. The 137 freeze and 205 test subjects therefore cannot enter.
# This step does not optimize designs or fit a new choice model. The original
# H2 FAIL is retained; a passing variant would be reported as EXPLORATORY.
#
# E4's modeling functions are defined inside 25_mechanism_model.R rather than
# the shared Study-1 engine. Sourcing that script would rerun the H2 gate,
# so its functions are restated here. Before any LOFO result is reported, the
# full model must reproduce E4's committed cv_performance.csv statistics and
# bootstrap interval to < 1e-10 using the saved seeds and analysis samples.
# Each variant uses the same folds and bootstrap resamples as the full model,
# making the drop-one differences paired. The shared Study-1 engine is sourced
# without redefinition.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

res_dir <- here("output", "elevation", "results", "dispb")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  invisible(df)
}
t_start <- Sys.time()
say <- function(...) cat(sprintf(...), sep = "")

# =============================================================================
# 0. Frozen E4 constants, restated (the committed E4 CV seeds and E4 constants
#    are reused); NOTHING is re-tuned.
# =============================================================================
SEED_E4     <- 20260720L
N_REP       <- 200L
TRAIN_FRAC  <- 0.70
N_BOOT      <- 2000L
MIN_CELL_N  <- 40L
TARGETS     <- c("loss_avoid_adj", "ev_sens_agn")
K_GRID      <- c(6L, 8L, 10L, 12L)
K_GATE      <- c(8L, 10L, 12L)
FEATURES    <- c("F1_target_predictor_variance", "F2_loss_exposure",
                 "F3_predictor_collinearity", "F4_expected_information",
                 "F5_choice_entropy")
GATE_SPEARMAN_MIN <- 0.30
GATE_BOOT_SUPPORT <- 0.90
CARRIER_MIN_BETA  <- 0.05
GATE_TRACK_SIM  <- "sim_free105"
GATE_TRACK_REAL <- "real_form1"
GATE_TARGET     <- "loss_avoid_adj"

# B3 gate thresholds (fixed before results)
B3B_DELTA_MAX  <- 0.05    # |Delta OOF(drop F4)| on the real track
B3D_FULL_MIN   <- 0.30    # sim FULL pooled OOF
B3D_VARIANT_MIN <- 0.20   # each sim LOFO variant pooled OOF

# Metric convention, fixed in code before results and reported: the frozen E4
# gate reads a cell through the CONSERVATIVE point estimate
# min(oof_spearman_pooled, oof_spearman_repeatmean). Every B3 gate clause uses
# that same statistic, and the pooled and repeat-mean deltas are ALSO reported
# so any dependence on the choice is visible rather than implicit.

# ---- the leave-one-out variant grid ----------------------------------------
VARIANTS <- c(
  list(FULL = FEATURES),
  stats::setNames(lapply(FEATURES, function(f) setdiff(FEATURES, f)),
                  paste0("drop_", sub("_.*$", "", FEATURES))),
  list(drop_F2_F4 = setdiff(FEATURES, c("F2_loss_exposure",
                                        "F4_expected_information"))),
  stats::setNames(lapply(FEATURES, function(f) f),
                  paste0("only_", sub("_.*$", "", FEATURES))),
  list(BASE_category = "F2_loss_exposure"))
VNAMES <- names(VARIANTS)
stopifnot(!anyDuplicated(VNAMES))
say("variants: %s\n", paste(VNAMES, collapse = ", "))

# =============================================================================
# 1. Inputs: committed E2 / E4 outputs ONLY (no raw, no choices, no split)
# =============================================================================
# utils::read.csv is used for the numeric tables so the doubles round-trip
# bit-for-bit (the readr float-parse caveat documented in step 24); the FULL-model
# reproduction check below depends on exact bits.
p_feat  <- here("output", "elevation", "results", "e4", "set_features.csv")
p_cv    <- here("output", "elevation", "results", "e4", "cv_performance.csv")
p_sign  <- here("output", "elevation", "results", "e4", "coefficient_sign_consistency.csv")
p_seeds <- here("output", "elevation", "results", "e4", "e4_seeds.csv")
p_e2    <- here("output", "elevation", "results", "e2", "problem_information_dictionary.csv")
input_md5 <- vapply(c(e4_set_features = p_feat, e4_cv_performance = p_cv,
                      e4_sign_consistency = p_sign, e4_seeds = p_seeds,
                      e2_dictionary = p_e2),
                    function(f) unname(tools::md5sum(f)), character(1))
print(input_md5)

FEAT  <- utils::read.csv(p_feat, stringsAsFactors = FALSE)
E4CV  <- utils::read.csv(p_cv, stringsAsFactors = FALSE)
E4SGN <- utils::read.csv(p_sign, stringsAsFactors = FALSE)
E4SD  <- utils::read.csv(p_seeds, stringsAsFactors = FALSE)
D     <- utils::read.csv(p_e2, stringsAsFactors = FALSE)
stopifnot(nrow(D) == 105L)
stopifnot(all(FEATURES %in% names(FEAT)))
stopifnot(all(c("track", "target", "K", "set_id", "r_sb", "game_ids") %in%
                names(FEAT)))
tracks <- sort(unique(FEAT$track))
stopifnot(GATE_TRACK_SIM %in% tracks, GATE_TRACK_REAL %in% tracks)
say("tracks: %s\n", paste(tracks, collapse = ", "))

# The committed E4 seeds, used verbatim so LOFO is PAIRED to the frozen FULL
# fit; the derivation is re-asserted so a mis-keyed lookup cannot pass silently.
get_seed <- function(purpose, scope) {
  v <- E4SD$seed[E4SD$purpose == purpose & E4SD$scope == scope]
  stopifnot(length(v) == 1L)
  as.integer(v)
}
for (tk in c(GATE_TRACK_SIM, GATE_TRACK_REAL)) for (tg in TARGETS) {
  for (K in K_GRID) {
    s <- get_seed("cv_holdout", paste(tk, tg, K, sep = "|"))
    stopifnot(s == SEED_E4 + 1000L * match(tk, tracks) +
                100L * match(tg, TARGETS) + K)
  }
  s <- get_seed("cluster_bootstrap", paste(tk, tg, sep = "|"))
  stopifnot(s == SEED_E4 + 500000L + 1000L * match(tk, tracks) +
              100L * match(tg, TARGETS))
}
say("committed E4 seeds: all cv_holdout and cluster_bootstrap keys reproduce their derivation\n")

# =============================================================================
# 2. E4 modeling primitives, restated verbatim (fidelity PROVEN in section 4)
# =============================================================================
safe_lsq <- function(X, y) {
  f <- stats::lm.fit(X, y)
  b <- f$coefficients; b[is.na(b)] <- 0
  as.numeric(b)
}

sp <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
}

# E4's analysis sample for a cell: every set with a finite r_sb AND all FIVE
# pre-declared features finite. It is held FIXED across variants, so every
# variant is fit on the identical rows and the folds match the frozen run.
cell_data <- function(track, target, K) {
  d <- FEAT[FEAT$track == track & FEAT$target == target & FEAT$K == K, ,
            drop = FALSE]
  keep <- is.finite(d$r_sb) &
    Reduce(`&`, lapply(FEATURES, function(f) is.finite(d[[f]])))
  list(all_n = nrow(d), d = d[keep, , drop = FALSE])
}

set_problem_matrix <- function(d, pool) {
  ids <- strsplit(d$game_ids, ";", fixed = TRUE)
  M <- matrix(0, nrow(d), length(pool))
  lk <- stats::setNames(seq_along(pool), as.character(pool))
  for (i in seq_along(ids)) M[i, lk[ids[[i]]]] <- 1
  M / rowSums(M)
}

# One cell: E4's repeated set-level 70/30 holdout within K, with EVERY variant
# scored on the SAME train/test partition in every repeat.
cv_cell_variants <- function(d, seed) {
  n <- nrow(d)
  X <- as.matrix(d[, FEATURES]); y <- d$r_sb
  n_tr <- max(5L, floor(TRAIN_FRAC * n))
  if (n - n_tr < 3L) return(NULL)
  nv <- length(VNAMES)
  pred_sum <- matrix(0, n, nv, dimnames = list(NULL, VNAMES))
  pred_cnt <- matrix(0L, n, nv, dimnames = list(NULL, VNAMES))
  rep_sp <- matrix(NA_real_, N_REP, nv, dimnames = list(NULL, VNAMES))
  rep_r2 <- matrix(NA_real_, N_REP, nv, dimnames = list(NULL, VNAMES))
  for (r in seq_len(N_REP)) {
    set.seed(seed + r)                       # E4's outer holdout stream
    tr <- sample.int(n, n_tr); te <- setdiff(seq_len(n), tr)
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    mu <- colMeans(Xtr)
    sg <- apply(Xtr, 2, stats::sd)
    sg[!is.finite(sg) | sg == 0] <- 1
    Ztr <- sweep(sweep(Xtr, 2, mu, "-"), 2, sg, "/")
    Zte <- sweep(sweep(Xte, 2, mu, "-"), 2, sg, "/")
    sst <- sum((y[te] - mean(y[tr]))^2)
    for (v in VNAMES) {
      j <- match(VARIANTS[[v]], FEATURES)
      p <- as.numeric(cbind(1, Zte[, j, drop = FALSE]) %*%
                        safe_lsq(cbind(1, Ztr[, j, drop = FALSE]), y[tr]))
      pred_sum[te, v] <- pred_sum[te, v] + ifelse(is.finite(p), p, 0)
      pred_cnt[te, v] <- pred_cnt[te, v] + as.integer(is.finite(p))
      rep_sp[r, v] <- sp(p, y[te])
      rep_r2[r, v] <- if (sst > 0 && all(is.finite(p)))
        1 - sum((y[te] - p)^2) / sst else NA_real_
    }
  }
  oof <- pred_sum / ifelse(pred_cnt == 0L, NA, pred_cnt)
  list(y = y, oof = oof, rep_sp = rep_sp, rep_r2 = rep_r2,
       n = n, n_train = n_tr, n_test = n - n_tr)
}

# E4's cluster bootstrap over PROBLEMS, with the SAME problem draw shared
# across K and across variants within a replicate, so every delta is paired.
bootstrap_cell_group <- function(cells, pool, seed) {
  np <- length(pool)
  ks <- names(cells)
  out_sp <- array(NA_real_, c(N_BOOT, length(ks), length(VNAMES)),
                  dimnames = list(NULL, ks, VNAMES))
  set.seed(seed)
  for (b in seq_len(N_BOOT)) {
    cnt <- as.numeric(stats::rmultinom(1, np, rep(1 / np, np)))
    for (ki in seq_along(ks)) {
      cc <- cells[[ks[ki]]]
      w <- as.numeric(cc$W %*% cnt)
      if (sum(w) <= 0) next
      idx <- sample.int(length(w), length(w), replace = TRUE, prob = w)
      yb <- cc$y[idx]
      for (v in VNAMES) out_sp[b, ki, v] <- sp(cc$oof[idx, v], yb)
    }
  }
  out_sp
}

ci_lo <- function(v) if (all(is.na(v))) NA_real_ else
  unname(stats::quantile(v, 0.025, na.rm = TRUE))
ci_hi <- function(v) if (all(is.na(v))) NA_real_ else
  unname(stats::quantile(v, 0.975, na.rm = TRUE))

# =============================================================================
# 3. Run every gating track x target x K
# =============================================================================
say("\n=== B3: LOFO cross-validation on the frozen E4 folds ===\n")
cv_rows <- list(); coef_rows <- list(); pooled_rows <- list()
boot_store <- list(); cell_store <- list()

for (tk in c(GATE_TRACK_SIM, GATE_TRACK_REAL)) {
  pool <- sort(unique(as.integer(unlist(strsplit(
    FEAT$game_ids[FEAT$track == tk], ";", fixed = TRUE)))))
  for (tg in TARGETS) {
    cells <- list()
    for (K in K_GRID) {
      cd <- cell_data(tk, tg, K)
      d <- cd$d
      if (nrow(d) < MIN_CELL_N) next
      seed_k <- get_seed("cv_holdout", paste(tk, tg, K, sep = "|"))
      cv <- cv_cell_variants(d, seed_k)
      if (is.null(cv)) next
      cells[[as.character(K)]] <- list(y = cv$y, oof = cv$oof,
                                       W = set_problem_matrix(d, pool))
      cell_store[[paste(tk, tg, K)]] <- list(d = d, cv = cv)

      # full-sample standardized coefficients per variant (E4's construction)
      Xs <- scale(as.matrix(d[, FEATURES]))
      Xs[!is.finite(Xs)] <- 0
      ysd <- stats::sd(d$r_sb)
      for (v in VNAMES) {
        j <- match(VARIANTS[[v]], FEATURES)
        b <- safe_lsq(cbind(1, Xs[, j, drop = FALSE]), d$r_sb)
        coef_rows[[length(coef_rows) + 1L]] <- tibble::tibble(
          track = tk, target = tg, K = K, variant = v,
          kind = if (grepl("^sim", tk)) "SIMULATION" else "REAL",
          n_features = length(j), n_sets_used = nrow(d), sd_r_sb = ysd,
          term = c("(Intercept)", FEATURES[j]), beta = b,
          beta_over_sd_y = b / ysd)
      }
    }
    if (!length(cells)) next
    seed_b <- get_seed("cluster_bootstrap", paste(tk, tg, sep = "|"))
    B <- bootstrap_cell_group(cells, pool, seed_b)
    boot_store[[paste(tk, tg)]] <- B

    for (K in as.integer(names(cells))) {
      ki <- match(as.character(K), dimnames(B)[[2]])
      cv <- cell_store[[paste(tk, tg, K)]]$cv
      pooled <- vapply(VNAMES, function(v) sp(cv$oof[, v], cv$y), numeric(1))
      repmean <- colMeans(cv$rep_sp, na.rm = TRUE)[VNAMES]
      r2mean <- colMeans(cv$rep_r2, na.rm = TRUE)[VNAMES]
      cv_rows[[length(cv_rows) + 1L]] <- tibble::tibble(
        track = tk, target = tg, K = K, variant = VNAMES,
        kind = if (grepl("^sim", tk)) "SIMULATION" else "REAL",
        n_features = vapply(VARIANTS[VNAMES], length, integer(1)),
        n_sets_used = cv$n, n_train = cv$n_train, n_test = cv$n_test,
        oof_spearman_pooled = pooled,
        oof_spearman_repeatmean = repmean,
        oof_spearman_gate = pmin(pooled, repmean),
        oof_r2_repeatmean = r2mean,
        boot_ci_lo = vapply(VNAMES, function(v) ci_lo(B[, ki, v]), numeric(1)),
        boot_ci_hi = vapply(VNAMES, function(v) ci_hi(B[, ki, v]), numeric(1)),
        delta_vs_base_category = pooled - pooled[["BASE_category"]],
        delta_repeatmean_vs_base_category = repmean - repmean[["BASE_category"]],
        delta_r2_vs_base_category = r2mean - r2mean[["BASE_category"]],
        delta_boot_support = vapply(VNAMES, function(v)
          mean(B[, ki, v] - B[, ki, "BASE_category"] > 0, na.rm = TRUE),
          numeric(1)),
        drop_one_delta_pooled = pooled[["FULL"]] - pooled,
        drop_one_delta_repeatmean = repmean[["FULL"]] - repmean,
        drop_one_delta_gate = pmin(pooled, repmean)[["FULL"]] -
          pmin(pooled, repmean))
    }

    # pooled-across-K_GATE statistics, E4's pooled_stat construction
    ks <- intersect(as.character(K_GATE), dimnames(B)[[2]])
    if (length(ks)) {
      obs <- vapply(VNAMES, function(v) mean(vapply(ks, function(k)
        sp(cell_store[[paste(tk, tg, as.integer(k))]]$cv$oof[, v],
           cell_store[[paste(tk, tg, as.integer(k))]]$cv$y), numeric(1))),
        numeric(1))
      obs_rm <- vapply(VNAMES, function(v) mean(vapply(ks, function(k)
        colMeans(cell_store[[paste(tk, tg, as.integer(k))]]$cv$rep_sp,
                 na.rm = TRUE)[[v]], numeric(1))), numeric(1))
      draws <- vapply(VNAMES, function(v) rowMeans(B[, ks, v, drop = FALSE]),
                      numeric(N_BOOT))
      pooled_rows[[length(pooled_rows) + 1L]] <- tibble::tibble(
        track = tk, target = tg, K = NA_integer_, variant = VNAMES,
        kind = if (grepl("^sim", tk)) "SIMULATION" else "REAL",
        K_pooled = paste(ks, collapse = ";"),
        n_features = vapply(VARIANTS[VNAMES], length, integer(1)),
        oof_spearman_pooled = obs, oof_spearman_repeatmean = obs_rm,
        oof_spearman_gate = pmin(obs, obs_rm),
        boot_ci_lo = vapply(VNAMES, function(v) ci_lo(draws[, v]), numeric(1)),
        boot_ci_hi = vapply(VNAMES, function(v) ci_hi(draws[, v]), numeric(1)),
        delta_vs_base_category = obs - obs[["BASE_category"]],
        delta_boot_support = vapply(VNAMES, function(v)
          mean(draws[, v] - draws[, "BASE_category"] > 0, na.rm = TRUE),
          numeric(1)),
        drop_one_delta_pooled = obs[["FULL"]] - obs,
        drop_one_delta_gate = pmin(obs, obs_rm)[["FULL"]] - pmin(obs, obs_rm))
    }
    say("  %s | %s: %d evaluable K cells\n", tk, tg, length(cells))
  }
}
CV <- dplyr::bind_rows(cv_rows)
POOL <- dplyr::bind_rows(pooled_rows)
COEF <- dplyr::bind_rows(coef_rows)

# =============================================================================
# 4. FIDELITY CHECK: the FULL model must reproduce E4's committed numbers
# =============================================================================
say("\n=== fidelity: FULL and BASE_category vs committed E4 cv_performance ===\n")
chk <- CV |>
  dplyr::filter(variant %in% c("FULL", "BASE_category")) |>
  dplyr::select(track, target, K, variant, oof_spearman_pooled,
                oof_spearman_repeatmean, boot_ci_lo, boot_ci_hi) |>
  dplyr::inner_join(
    E4CV |> dplyr::filter(model %in% c("FULL", "BASE_category")) |>
      dplyr::select(track, target, K, variant = model,
                    e4_pooled = oof_spearman_pooled,
                    e4_repmean = oof_spearman_repeatmean,
                    e4_ci_lo = boot_ci_lo, e4_ci_hi = boot_ci_hi),
    by = c("track", "target", "K", "variant"))
d_pooled <- max(abs(chk$oof_spearman_pooled - chk$e4_pooled), na.rm = TRUE)
d_repmean <- max(abs(chk$oof_spearman_repeatmean - chk$e4_repmean), na.rm = TRUE)
d_ci <- max(abs(c(chk$boot_ci_lo - chk$e4_ci_lo,
                  chk$boot_ci_hi - chk$e4_ci_hi)), na.rm = TRUE)
say("rows compared %d; max |d| pooled = %.3g, repeatmean = %.3g, bootstrap CI = %.3g\n",
    nrow(chk), d_pooled, d_repmean, d_ci)
stopifnot(nrow(chk) >= 16L, d_pooled < 1e-10, d_repmean < 1e-10, d_ci < 1e-10)

# only_F2 is, by construction, the same fit as BASE_category.
id_chk <- CV |> dplyr::filter(variant %in% c("only_F2", "BASE_category")) |>
  dplyr::select(track, target, K, variant, oof_spearman_pooled) |>
  tidyr::pivot_wider(names_from = variant, values_from = oof_spearman_pooled)
stopifnot(max(abs(id_chk$only_F2 - id_chk$BASE_category), na.rm = TRUE) < 1e-12)
say("only_F2 == BASE_category confirmed on all %d cells\n", nrow(id_chk))

# The committed E4 sign table must also reproduce from the FULL coefficients.
sgn_chk <- COEF |>
  dplyr::filter(variant == "FULL", target == GATE_TARGET, K %in% K_GATE,
                term %in% FEATURES) |>
  dplyr::select(track, K, term, beta) |>
  tidyr::pivot_wider(names_from = track, values_from = beta) |>
  dplyr::inner_join(E4SGN |> dplyr::select(term, K, beta_sim, beta_real),
                    by = c("term", "K"))
d_sgn <- max(abs(c(sgn_chk[[GATE_TRACK_SIM]] - sgn_chk$beta_sim,
                   sgn_chk[[GATE_TRACK_REAL]] - sgn_chk$beta_real)),
             na.rm = TRUE)
say("committed E4 coefficient_sign_consistency reproduced, max |d| = %.3g over %d terms\n",
    d_sgn, nrow(sgn_chk))
stopifnot(d_sgn < 1e-10)

# =============================================================================
# 5. Per-variant re-run of the frozen H2 three-clause gate
# =============================================================================
getcv <- function(tk, tg, K, v, col) {
  x <- CV[[col]][CV$track == tk & CV$target == tg & CV$K == K & CV$variant == v]
  if (!length(x)) NA_real_ else x[1]
}
getpool <- function(tk, tg, v, col) {
  x <- POOL[[col]][POOL$track == tk & POOL$target == tg & POOL$variant == v]
  if (!length(x)) NA_real_ else x[1]
}
beta_of <- function(tk, K, v, term) {
  x <- COEF$beta[COEF$track == tk & COEF$target == GATE_TARGET & COEF$K == K &
                   COEF$variant == v & COEF$term == term]
  if (!length(x)) NA_real_ else x[1]
}

variant_gate <- function(v) {
  # G-sim: simulation track clears 0.30 with a positive CI lower bound at every
  # gate K and pooled (E4's rule, unchanged).
  gsim_k <- lapply(K_GATE, function(K) list(
    K = K,
    spearman_gate = getcv(GATE_TRACK_SIM, GATE_TARGET, K, v, "oof_spearman_gate"),
    ci_lo = getcv(GATE_TRACK_SIM, GATE_TARGET, K, v, "boot_ci_lo"),
    meets = isTRUE(getcv(GATE_TRACK_SIM, GATE_TARGET, K, v,
                         "oof_spearman_gate") >= GATE_SPEARMAN_MIN),
    ci_positive = isTRUE(getcv(GATE_TRACK_SIM, GATE_TARGET, K, v,
                               "boot_ci_lo") > 0)))
  gsim_pooled_stat <- getpool(GATE_TRACK_SIM, GATE_TARGET, v, "oof_spearman_gate")
  gsim_pooled_ci <- getpool(GATE_TRACK_SIM, GATE_TARGET, v, "boot_ci_lo")
  gsim <- all(vapply(gsim_k, function(x) x$meets && x$ci_positive, logical(1))) &&
    isTRUE(gsim_pooled_stat >= GATE_SPEARMAN_MIN) && isTRUE(gsim_pooled_ci > 0)

  # G-beyond: positive margin over the fixed F2-only baseline at every gate K,
  # with bootstrap support at or above 0.90.
  gb_k <- lapply(K_GATE, function(K) list(
    K = K,
    delta = getcv(GATE_TRACK_SIM, GATE_TARGET, K, v, "delta_vs_base_category"),
    delta_r2 = getcv(GATE_TRACK_SIM, GATE_TARGET, K, v,
                     "delta_r2_vs_base_category"),
    support = getcv(GATE_TRACK_SIM, GATE_TARGET, K, v, "delta_boot_support"),
    positive_margin = isTRUE(getcv(GATE_TRACK_SIM, GATE_TARGET, K, v,
                                   "delta_vs_base_category") > 0) ||
      isTRUE(getcv(GATE_TRACK_SIM, GATE_TARGET, K, v,
                   "delta_r2_vs_base_category") > 0),
    support_ok = isTRUE(getcv(GATE_TRACK_SIM, GATE_TARGET, K, v,
                              "delta_boot_support") >= GATE_BOOT_SUPPORT)))
  gbeyond <- all(vapply(gb_k, function(x) x$positive_margin && x$support_ok,
                        logical(1)))

  # G-real: the real anchor is positive, beats the baseline, and the signs of
  # the SIMULATION model's carrier features agree with the real model's, using
  # E4's committed carrier rule |beta_sim| >= 0.05 over the variant's RETAINED
  # features (the rule is not re-tuned).
  gr_k <- lapply(K_GATE, function(K) {
    terms <- VARIANTS[[v]]
    bs <- vapply(terms, function(tm) beta_of(GATE_TRACK_SIM, K, v, tm), numeric(1))
    br <- vapply(terms, function(tm) beta_of(GATE_TRACK_REAL, K, v, tm), numeric(1))
    carrier <- abs(bs) >= CARRIER_MIN_BETA
    agree <- sign(bs) == sign(br)
    list(K = K,
         spearman_gate = getcv(GATE_TRACK_REAL, GATE_TARGET, K, v,
                               "oof_spearman_gate"),
         delta_vs_base = getcv(GATE_TRACK_REAL, GATE_TARGET, K, v,
                               "delta_vs_base_category"),
         n_carrier = sum(carrier),
         n_carrier_sign_agree = sum(carrier & agree),
         n_all_sign_agree = sum(agree), n_terms = length(terms),
         positive = isTRUE(getcv(GATE_TRACK_REAL, GATE_TARGET, K, v,
                                 "oof_spearman_gate") > 0),
         beats_base = isTRUE(getcv(GATE_TRACK_REAL, GATE_TARGET, K, v,
                                   "delta_vs_base_category") > 0),
         carrier_signs_consistent = sum(carrier) > 0 && all(agree[carrier]))
  })
  greal <- all(vapply(gr_k, function(x)
    x$positive && x$beats_base && x$carrier_signs_consistent, logical(1)))

  list(variant = v, retains_F4 = "F4_expected_information" %in% VARIANTS[[v]],
       retains_F1 = "F1_target_predictor_variance" %in% VARIANTS[[v]],
       g_sim = list(per_K = gsim_k, pooled_gate = gsim_pooled_stat,
                    pooled_ci_lo = gsim_pooled_ci, pass = gsim),
       g_beyond = list(per_K = gb_k, pass = gbeyond),
       g_real = list(per_K = gr_k, carrier_rule = CARRIER_MIN_BETA,
                     pass = greal),
       h2_pass = gsim && gbeyond && greal)
}
VG <- lapply(VNAMES, variant_gate)
names(VG) <- VNAMES

# The FULL variant must reproduce the frozen E4 H2 FAIL (a fidelity check on
# the gate logic itself; the frozen h2_gate.json is NOT re-run or amended).
h2_frozen <- jsonlite::read_json(
  here("output", "elevation", "results", "e4", "h2_gate.json"), simplifyVector = TRUE)
full_match <- identical(VG$FULL$h2_pass, h2_frozen$h2_pass) &&
  identical(VG$FULL$g_sim$pass, h2_frozen$g_sim$pass) &&
  identical(VG$FULL$g_beyond$pass, h2_frozen$g_beyond$pass) &&
  identical(VG$FULL$g_real$pass, h2_frozen$g_real$pass)
say("FULL variant reproduces the frozen H2 clause pattern (G-sim %s, G-beyond %s, G-real %s, H2 %s): %s\n",
    VG$FULL$g_sim$pass, VG$FULL$g_beyond$pass, VG$FULL$g_real$pass,
    VG$FULL$h2_pass, full_match)
stopifnot(full_match)

# =============================================================================
# 6. Outputs
# =============================================================================
b3_lofo_cv <- dplyr::bind_rows(
  CV |> dplyr::mutate(scope = "within_K", K_pooled = NA_character_),
  POOL |> dplyr::mutate(scope = "pooled_K", n_sets_used = NA_integer_,
                        n_train = NA_integer_, n_test = NA_integer_,
                        oof_r2_repeatmean = NA_real_,
                        delta_repeatmean_vs_base_category = NA_real_,
                        delta_r2_vs_base_category = NA_real_,
                        drop_one_delta_repeatmean = NA_real_)) |>
  dplyr::mutate(
    h2_gate_pass = vapply(variant, function(v) VG[[v]]$h2_pass, logical(1)),
    g_sim_pass = vapply(variant, function(v) VG[[v]]$g_sim$pass, logical(1)),
    g_beyond_pass = vapply(variant, function(v) VG[[v]]$g_beyond$pass, logical(1)),
    g_real_pass = vapply(variant, function(v) VG[[v]]$g_real$pass, logical(1))) |>
  dplyr::select(track, kind, target, scope, K, K_pooled, variant, n_features,
                n_sets_used, n_train, n_test, oof_spearman_pooled,
                oof_spearman_repeatmean, oof_spearman_gate, oof_r2_repeatmean,
                boot_ci_lo, boot_ci_hi, delta_vs_base_category,
                delta_repeatmean_vs_base_category, delta_r2_vs_base_category,
                delta_boot_support, drop_one_delta_pooled,
                drop_one_delta_repeatmean, drop_one_delta_gate,
                g_sim_pass, g_beyond_pass, g_real_pass, h2_gate_pass) |>
  dplyr::arrange(track, target, scope, K, variant)
save_csv(b3_lofo_cv, "b3_lofo_cv")

# corr(F2, F4) on the real Set-1 table (the collinearity B3 is characterising)
f24 <- purrr::map_dfr(c(GATE_TRACK_REAL, GATE_TRACK_SIM), function(tk)
  purrr::map_dfr(TARGETS, function(tg) purrr::map_dfr(K_GRID, function(K) {
    d <- cell_data(tk, tg, K)$d
    if (nrow(d) < MIN_CELL_N) return(tibble::tibble())
    tibble::tibble(track = tk, target = tg, K = K, n_sets = nrow(d),
                   pearson_F2_F4 = stats::cor(d$F2_loss_exposure,
                                              d$F4_expected_information),
                   spearman_F2_F4 = stats::cor(d$F2_loss_exposure,
                                               d$F4_expected_information,
                                               method = "spearman"))
  })))
save_csv(f24, "b3_feature_collinearity")
say("\ncorr(F2, F4), real track, primary target:\n")
print(as.data.frame(f24 |> dplyr::filter(track == GATE_TRACK_REAL,
                                         target == GATE_TARGET)))

b3_lofo_signs <- COEF |>
  dplyr::filter(term %in% FEATURES) |>
  dplyr::mutate(sign_beta = sign(beta)) |>
  dplyr::left_join(
    COEF |> dplyr::filter(track == GATE_TRACK_SIM, term %in% FEATURES) |>
      dplyr::select(target, K, variant, term, beta_sim = beta),
    by = c("target", "K", "variant", "term")) |>
  dplyr::mutate(
    carrier = abs(beta_sim) >= CARRIER_MIN_BETA,
    sign_agrees_with_sim = sign(beta) == sign(beta_sim),
    h2_gate_pass = vapply(variant, function(v) VG[[v]]$h2_pass, logical(1)),
    g_sim_pass = vapply(variant, function(v) VG[[v]]$g_sim$pass, logical(1)),
    g_beyond_pass = vapply(variant, function(v) VG[[v]]$g_beyond$pass, logical(1)),
    g_real_pass = vapply(variant, function(v) VG[[v]]$g_real$pass, logical(1))) |>
  dplyr::left_join(f24 |> dplyr::select(track, target, K, pearson_F2_F4,
                                        spearman_F2_F4),
                   by = c("track", "target", "K")) |>
  dplyr::select(track, kind, target, K, variant, n_features, term, beta,
                beta_over_sd_y, sign_beta, beta_sim, carrier,
                sign_agrees_with_sim, pearson_F2_F4, spearman_F2_F4,
                g_sim_pass, g_beyond_pass, g_real_pass, h2_gate_pass) |>
  dplyr::arrange(track, target, K, variant, term)
save_csv(b3_lofo_signs, "b3_lofo_signs")

# =============================================================================
# 7. Gates G-B3a..d
# =============================================================================
f4_variants <- VNAMES[vapply(VNAMES, function(v)
  VG[[v]]$retains_F4 && grepl("^drop_", v), logical(1))]
f4_pass <- vapply(f4_variants, function(v) VG[[v]]$h2_pass, logical(1))
g_b3a <- !any(f4_pass)

dropF4 <- CV |>
  dplyr::filter(track == GATE_TRACK_REAL, target == GATE_TARGET,
                K %in% K_GATE, variant == "drop_F4")
g_b3b_gate <- max(abs(dropF4$drop_one_delta_gate), na.rm = TRUE)
g_b3b_pooled <- max(abs(dropF4$drop_one_delta_pooled), na.rm = TRUE)
g_b3b_repmean <- max(abs(dropF4$drop_one_delta_repeatmean), na.rm = TRUE)
g_b3b <- isTRUE(g_b3b_gate <= B3B_DELTA_MAX)
g_b3b_all_metrics <- isTRUE(max(g_b3b_gate, g_b3b_pooled, g_b3b_repmean) <=
                              B3B_DELTA_MAX)

f1_rows <- COEF |>
  dplyr::filter(track == GATE_TRACK_REAL, target == GATE_TARGET,
                K %in% K_GATE, term == "F1_target_predictor_variance")
g_b3c <- nrow(f1_rows) > 0 && all(f1_rows$beta < 0)

sim_pooled <- POOL |>
  dplyr::filter(track == GATE_TRACK_SIM, target == GATE_TARGET)
sim_full <- sim_pooled$oof_spearman_gate[sim_pooled$variant == "FULL"]
sim_lofo <- sim_pooled |> dplyr::filter(grepl("^drop_F[1-5]$", variant))
g_b3d <- isTRUE(sim_full >= B3D_FULL_MIN) &&
  all(sim_lofo$oof_spearman_gate >= B3D_VARIANT_MIN)

say("\n=== B3 GATES ===\n")
say("G-B3a null robustness: F4-retaining LOFO variants passing the full frozen H2 gate = %d of %d (%s) -> %s\n",
    sum(f4_pass), length(f4_pass), paste(f4_variants, collapse = ","),
    if (g_b3a) "PASS" else "FAIL")
say("G-B3b information redundancy: max |Delta OOF(drop F4)| on the real track over K %s = %.5f (gate metric; pooled %.5f, repeatmean %.5f) vs %.2f -> %s\n",
    paste(K_GATE, collapse = ","), g_b3b_gate, g_b3b_pooled, g_b3b_repmean,
    B3B_DELTA_MAX, if (g_b3b) "PASS" else "FAIL")
say("G-B3c F1 stability: %d of %d F1-retaining real-track coefficients negative -> %s\n",
    sum(f1_rows$beta < 0), nrow(f1_rows), if (g_b3c) "PASS" else "FAIL")
say("G-B3d sim sanity: FULL pooled %.4f (>= %.2f); LOFO min %.4f (>= %.2f) -> %s\n",
    sim_full, B3D_FULL_MIN, min(sim_lofo$oof_spearman_gate), B3D_VARIANT_MIN,
    if (g_b3d) "PASS" else "FAIL")

real_deltas <- CV |>
  dplyr::filter(track == GATE_TRACK_REAL, target == GATE_TARGET,
                K %in% K_GATE, grepl("^drop_F[1-5]$", variant)) |>
  dplyr::select(K, variant, drop_one_delta_gate, drop_one_delta_pooled)
say("\nreal-track drop-one deltas (gate metric), primary target:\n")
print(as.data.frame(tidyr::pivot_wider(
  real_deltas |> dplyr::select(K, variant, drop_one_delta_gate),
  names_from = variant, values_from = drop_one_delta_gate)))

part <- list(
  gate = list(
    probe = "B3 leave-one-feature-out on the frozen E4 mechanism model",
    frozen_h2 = list(verdict = h2_frozen$h2_verdict,
                     rerun_or_amended = FALSE,
                     full_variant_reproduces_clause_pattern = full_match),
    metric_convention = paste("every clause uses the frozen E4 conservative",
                              "statistic min(oof_spearman_pooled,",
                              "oof_spearman_repeatmean); pooled and repeat-mean",
                              "deltas are also reported"),
    `G-B3a` = list(
      criterion = paste("no F4-retaining leave-one-out variant restores the",
                        "full frozen H2 PASS on the gating tracks"),
      f4_retaining_variants = f4_variants,
      n_passing = sum(f4_pass),
      per_variant = lapply(f4_variants, function(v) list(
        variant = v, g_sim = VG[[v]]$g_sim$pass, g_beyond = VG[[v]]$g_beyond$pass,
        g_real = VG[[v]]$g_real$pass, h2_pass = VG[[v]]$h2_pass)),
      verdict = if (g_b3a) "PASS" else "FAIL", pass = g_b3a),
    `G-B3b` = list(
      threshold = B3B_DELTA_MAX,
      max_abs_delta_gate_metric = g_b3b_gate,
      max_abs_delta_pooled = g_b3b_pooled,
      max_abs_delta_repeatmean = g_b3b_repmean,
      per_K = lapply(seq_len(nrow(dropF4)), function(i)
        as.list(dropF4[i, c("K", "drop_one_delta_gate",
                            "drop_one_delta_pooled",
                            "drop_one_delta_repeatmean")])),
      holds_on_all_three_metrics = g_b3b_all_metrics,
      verdict = if (g_b3b) "PASS" else "FAIL", pass = g_b3b),
    `G-B3c` = list(
      criterion = "F1 standardized coefficient negative in FULL and every F1-retaining variant on the real track",
      n_coefficients = nrow(f1_rows), n_negative = sum(f1_rows$beta < 0),
      max_beta = max(f1_rows$beta), min_beta = min(f1_rows$beta),
      verdict = if (g_b3c) "PASS" else "FAIL", pass = g_b3c),
    `G-B3d` = list(
      full_min = B3D_FULL_MIN, variant_min = B3D_VARIANT_MIN,
      sim_full_pooled = sim_full,
      sim_lofo_pooled = stats::setNames(as.list(sim_lofo$oof_spearman_gate),
                                        sim_lofo$variant),
      verdict = if (g_b3d) "PASS" else "FAIL", pass = g_b3d),
    feature_collinearity = list(
      real_track_primary_target = lapply(seq_len(nrow(f24)), function(i)
        as.list(f24[i, ]))[
          which(f24$track == GATE_TRACK_REAL & f24$target == GATE_TARGET)])),
  config = list(
    probe = "B3", seed_source = "elevation/results/e4/e4_seeds.csv (verbatim)",
    seed_base_e4 = SEED_E4, n_cv_repeats = N_REP, train_fraction = TRAIN_FRAC,
    n_bootstrap = N_BOOT, min_cell_n = MIN_CELL_N,
    features = FEATURES, variants = VNAMES,
    variant_features = lapply(VARIANTS, identity),
    targets = TARGETS, K_grid = K_GRID, K_gate = K_GATE,
    tracks = c(GATE_TRACK_SIM, GATE_TRACK_REAL),
    carrier_rule = sprintf("|beta_sim| >= %.2f over the variant's RETAINED features (E4's committed rule, not re-tuned)",
                           CARRIER_MIN_BETA),
    baseline_note = paste("BASE_category stays the fixed F2-only baseline for",
                          "every variant, including drop_F2, so the margin is",
                          "measured against the SAME baseline E4 used"),
    analysis_sample_note = paste("the cell sample is E4's (finite r_sb and all",
                                 "FIVE features finite) and is held FIXED",
                                 "across variants, so n, the holdout folds and",
                                 "the bootstrap resamples are identical to the",
                                 "frozen FULL fit and every delta is paired"),
    fidelity = list(rows_compared = nrow(chk),
                    max_abs_diff_oof_pooled = d_pooled,
                    max_abs_diff_oof_repeatmean = d_repmean,
                    max_abs_diff_bootstrap_ci = d_ci,
                    max_abs_diff_committed_sign_table = d_sgn),
    input_files = as.list(input_md5),
    raw_data_read = FALSE, participant_choices_read = FALSE,
    e1_split_read = FALSE, reliability_recomputed = FALSE,
    build_only_mechanism = "structural: no participant-level record is loaded by this script",
    freeze_read = FALSE, test_read = FALSE,
    optimization_performed = FALSE, candidate_design_constructed = FALSE,
    new_choice_model_fitted = FALSE, e4_h2_gate_rerun_or_amended = FALSE,
    runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
jsonlite::write_json(part, file.path(res_dir, "dispb_part_b3.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

seeds_b3 <- E4SD |>
  dplyr::filter(purpose %in% c("cv_holdout", "cluster_bootstrap")) |>
  dplyr::filter(grepl(paste0("^(", GATE_TRACK_SIM, "|", GATE_TRACK_REAL, ")\\|"),
                      scope)) |>
  dplyr::transmute(probe = "B3", scope = scope, purpose = purpose,
                   seed = as.integer(seed),
                   derivation = paste("committed E4 seed, reused VERBATIM:",
                                      derivation))
readr::write_csv(seeds_b3, file.path(res_dir, "dispb_seeds_b3.csv"))
say("committed E4 seeds reused: %d rows\n", nrow(seeds_b3))

say("\nB3 complete in %.1f min. Outputs in elevation/results/dispb/\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))

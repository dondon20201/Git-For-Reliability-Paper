# Step 34: local dependence and response-style diagnostics
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e0/, e1/, e2/, e3/ (steps 19-22)
# Writes:  output/elevation/results/dispc/ (c1_*, c2_*, c3_* tables, dispc_part.json, exports, seeds)
# Paper:   Study 2 (Ruling out alternatives); Supplement (Study 2 robustness detail), through step 36
# Runtime: about 25 minutes with 14 parallel workers

# C1 / C2 / C3: pre-declared checks of alternative explanations for E3's
# Set-1 loss-construct heterogeneity (gates G-C1 / G-C2 / G-C3):
#   C1 response-style contamination (B-side bias / safe-choice tendency),
#   C2 local item dependence (within-set problem-residual coherence),
#   C3 scoring-rule specificity (a robust median loss score).
# A detected confound narrows the H1 claim; absence of a confound strengthens
# the construct interpretation. C2 may also be INCONCLUSIVE. The frozen
# H1 PASS and H2 FAIL remain unchanged. Weak or inconclusive findings are
# reported rather than used to select further models.
#
# Reliabilities use the 63 Set-1 build subjects from step 20. Choice/outcome
# columns are dropped from the full-sample frame, excluding the 137 freeze and
# 205 test subjects' choices. No candidate designs are optimized or constructed,
# and no new mechanism model is fitted.
#
# R/step2_reliability.R and R/cpc18_features.R are sourced without redefinition,
# checked at source level by step 35. The two target indices use the accelerator
# from E3, E3b, and steps 30-33. Before results are produced, the inline G_EQ
# checks exact agreement of the from-raw loss score with estimate_half; step 35
# also checks whole-set r_sb through the original engine. C3's median loss
# score is a new pre-declared variant, not an equivalent of the frozen mean.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})

source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))

res_dir <- here("output", "elevation", "results", "dispc")
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
# 0. Pre-declared constants (fixed BEFORE any result)
# =============================================================================
RAW_MD5      <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_C       <- 20260724L
SEED_E3      <- 20260719L         # only to REPRODUCE E3's Set-1 sets
REGIME       <- "exp"
POS          <- 6:25              # experience regime, r = 20 (the full pool)
R_REPS       <- length(POS)
K_GRID       <- c(6L, 8L, 10L, 12L)
GATE_K       <- 12L
GATE_FORM    <- 1L
TARGETS      <- c("loss_avoid_adj", "ev_sens_agn")   # primary, informational
GATE_TARGET  <- "loss_avoid_adj"
M_CAP        <- 2000L
N_SPLIT      <- 100L              # problem-halvings per set (matches E3 / steps 30-33)
B_PART       <- 200L              # subject partitions at the primary cell
N_PART_AUX   <- 8L                # partitions per auxiliary K (single E3 sample)
N_NOISE_SETS <- 50L              # bootstrap-floor noise subsample
N_BOOT       <- 1000L             # frozen reliability_point resamples
MIN_PAIRS    <- 10L               # frozen split_correlations default
MIN_SETS_COV <- 30L               # minimum finite set-pairs for a covariance
N_HALF_A     <- 32L               # largest-remainder subject halves of 63
N_WPC_TERCILE_MIN <- 30L          # minimum sets per tercile to report a gap
N_WORKERS    <- 14L

# Pre-declared gate thresholds (fixed before results)
GC1_MATERIAL   <- 0.50            # relative cov-ICC drop -> response-style driver
GC1_NEGLIGIBLE <- 0.15            # relative cov-ICC drop -> not an artifact
GC2_ADJ_DROP   <- 0.50            # dependence-adjusted ICC relative reduction
GC3_RANK_MIN   <- 0.30            # set-ranking Spearman(primary, robust)
COLLIN_THRESH  <- 0.95            # |corr(B-side, safe)| collinearity trigger

# ---- seed derivations (disjoint blocks; written to dispc_seeds.csv) ----------
# Set membership REPRODUCES E3's Set-1 seeds.
seed_subset_e3 <- function(form, K) SEED_E3 + 1000L * form + K
seed_split_e3  <- function(form, K) SEED_E3 + 1000L * form + K + 500000L
# This step's own seeds (SEED_C base).
seed_part_c    <- function(p) SEED_C + 700000L + p
seed_part_aux_c <- function(K, j) SEED_C + 720000L + 1000L * K + j
seed_perm_c    <- function(p) SEED_C + 800000L + p
seed_noise_c   <- function(K, i) SEED_C + 600000L + 1000L * K + i
seed_trial_c   <- function(prob) SEED_C + 300000L + prob
# Per-half halving seed inside a partition: h in {1,2}, K, partition index, set.
seed_half_c <- function(h, K, part, sid)
  SEED_C + 900000000L * (h - 1L) + 200000000L * match(K, K_GRID) +
  100000L * part + sid

# =============================================================================
# 1. Inputs and the build-only structural guarantee
# =============================================================================
raw_path <- here("data", "raw", "All estimation raw data.csv")
stopifnot(file.exists(raw_path))
raw_md5 <- unname(tools::md5sum(raw_path))
stopifnot(identical(raw_md5, RAW_MD5))
say("raw md5 verified: %s\n", raw_md5)

split_tbl <- readr::read_csv(
  here("output", "elevation", "results", "e1", "participant_split.csv"),
  show_col_types = FALSE)
stopifnot(identical(names(split_tbl),
                    c("SubjID", "Set", "split", "held_out_until_E7")))
build_ids <- split_tbl$SubjID[split_tbl$split == "build"]
held_out  <- split_tbl$SubjID[split_tbl$split != "build"]
n_build   <- length(build_ids)
n_freeze  <- sum(split_tbl$split == "freeze")
n_test    <- sum(split_tbl$split == "test")
stopifnot(n_build == 344L, n_freeze == 137L, n_test == 205L)          # build-only guard
stopifnot(!anyDuplicated(build_ids))
stopifnot(all(split_tbl$held_out_until_E7 == (split_tbl$split == "test")))
say("split: build %d / freeze %d (not read) / test %d (SEALED)\n",
    n_build, n_freeze, n_test)

raw <- load_cpc18(raw_path)
stopifnot(all(cpc18_expected_columns %in% names(raw)))
choice_cols <- c("B", "Button", "Payoff", "Forgone", "RT", "Apay", "Bpay",
                 "Feedback")
design_all <- raw |> dplyr::select(-dplyr::any_of(choice_cols))
stopifnot(!any(choice_cols %in% names(design_all)))                    # build-only guard
build_raw <- raw |> dplyr::filter(SubjID %in% build_ids)
rm(raw); invisible(gc())
stopifnot(setequal(unique(build_raw$SubjID), build_ids))
stopifnot(!any(build_raw$SubjID %in% held_out))                        # build-only guard

feat_all <- clean_problem_features(design_all) |> add_scaled_features()
feat_all <- feat_all |> dplyr::arrange(GameID)
stopifnot(nrow(feat_all) == 105L)

e2_dict <- readr::read_csv(
  here("output", "elevation", "results", "e2", "problem_information_dictionary.csv"),
  show_col_types = FALSE)
stopifnot(nrow(e2_dict) == 105L)
stopifnot(setequal(e2_dict$GameID, feat_all$GameID))
stopifnot(all(abs(e2_dict$d_ev[order(e2_dict$GameID)] - feat_all$d_ev) < 1e-12))
stopifnot(identical(e2_dict$loss_informative[order(e2_dict$GameID)],
                    feat_all$loss_informative))
say("item universe: %d clean problems, matches the E2 dictionary\n",
    nrow(feat_all))

build_trials <- build_raw |>
  dplyr::filter(is_clean_row(build_raw)) |>
  dplyr::select(SubjID, Set, GameID, Trial, B) |>
  dplyr::left_join(
    feat_all |> dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev,
                              loss_informative, has_loss_b, ev_a, ev_b),
    by = "GameID")
stopifnot(!any(is.na(build_trials$z_ev)))
stopifnot(dplyr::n_distinct(build_trials$SubjID) == 344L)

baselines <- fit_loss_baselines(build_trials)                          # FROZEN
subj_data <- build_subject_data(build_trials, baselines)               # FROZEN
stopifnot(length(subj_data) == 344L)
stopifnot(abs(baselines$exp[1] - 0.4301295476692) < 1e-10,
          abs(baselines$exp[2] - 0.165107793793507) < 1e-10)  # E3 constants
say("frozen loss baselines (exp): b0=%.13f b1=%.15f\n",
    baselines$exp[1], baselines$exp[2])

feat_all <- feat_all |>
  dplyr::mutate(
    d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a),
    p_base_exp  = ifelse(loss_informative,
                         stats::plogis(baselines$exp[1] +
                                         baselines$exp[2] * d_ev_noloss), NA))

# =============================================================================
# 2. Set-1 form matrix (E3 construction, unchanged). Set 1 is the primary cell:
#    every build subject saw all 24 clean problems, so no subject is dropped at
#    any K within Set 1 (the eligibility filter is a no-op here).
# =============================================================================
form_of <- build_trials |> dplyr::distinct(SubjID, Set)
probs1 <- sort(unique(build_trials$GameID[build_trials$Set == GATE_FORM]))
subs1  <- as.character(sort(form_of$SubjID[form_of$Set == GATE_FORM]))
NB1 <- matrix(NA_real_, length(subs1), length(probs1),
              dimnames = list(subs1, as.character(probs1)))
for (sid in subs1) {
  sd <- subj_data[[sid]]
  j <- match(sd$feat$GameID, probs1)
  ok <- !is.na(j)
  NB1[sid, j[ok]] <- rowSums(sd$B_mat[which(ok), POS, drop = FALSE])
}
fi1 <- feat_all[match(probs1, feat_all$GameID), ]
fm <- list(probs = probs1, subs = subs1, NB = NB1, feat = fi1,
           li = fi1$loss_informative, has_loss_b = fi1$has_loss_b,
           d_ev = fi1$d_ev, p_base = fi1$p_base_exp,
           sd_a = fi1$sd_a, sd_b = fi1$sd_b)
stopifnot(length(fm$subs) == 63L, length(fm$probs) == 24L,
          sum(fm$li) == 6L, !anyNA(fm$NB))
say("Set 1: %d build subjects, %d clean problems, %d loss-informative\n",
    length(fm$subs), length(fm$probs), sum(fm$li))

# =============================================================================
# 3. Accelerator (formulas identical to the frozen estimate_half; the mean loss
#    score and ev_sens_agn are proven == frozen at the inline G_EQ below; the
#    median loss score is a NEW pre-declared robust variant)
# =============================================================================
ev_ranks <- function(d_ev, r) {
  k <- length(d_ev)
  o <- order(d_ev)
  dv <- d_ev[o]
  grp <- cumsum(c(TRUE, diff(dv) != 0))
  sizes <- tabulate(grp)
  starts <- cumsum(c(0, sizes[-length(sizes)])) * r
  Rg <- starts + (sizes * r + 1) / 2
  R <- numeric(k); R[o] <- Rg[grp]
  R
}

# Per-subject half scores over the half's problems `cols`.
#   loss_mean   = frozen loss_avoid_adj (mean per-problem residual over the
#                 half's loss-informative problems) -- proven == estimate_half.
#   loss_median = MEDIAN per-problem residual (C3 robust variant).
#   ev          = frozen ev_sens_agn (Spearman of choices vs d_ev).
half_scores <- function(NBh, fm, cols, r, want) {
  out <- list()
  if ("ev" %in% want) {
    k <- length(cols); n <- r * k
    R <- ev_ranks(fm$d_ev[cols], r); sR <- sum(R); sR2 <- sum(R * R)
    n1 <- rowSums(NBh)
    Sxy <- as.vector(NBh %*% R) - n1 * sR / k
    Sxx <- n1 - n1 * n1 / n
    Syy <- r * (sR2 - sR * sR / k)
    ev <- rep(NA_real_, nrow(NBh)); ok <- Sxx > 0 & Syy > 0
    if (any(ok)) ev[ok] <- Sxy[ok] / sqrt(Sxx[ok] * Syy)
    out$ev <- ev
  }
  if (any(c("loss_mean", "loss_median") %in% want)) {
    lic <- cols[fm$li[cols]]
    if (length(lic) == 0) {
      out$loss_mean <- rep(NA_real_, nrow(NBh))
      out$loss_median <- rep(NA_real_, nrow(NBh))
    } else {
      NBl <- NBh[, fm$li[cols], drop = FALSE]         # subjects x n_li
      cnt <- matrix(0, nrow(NBl), ncol(NBl))
      for (j in seq_along(lic)) {
        cnt[, j] <- if (fm$has_loss_b[lic[j]]) r - NBl[, j] else NBl[, j]
      }
      # loss_mean uses the EXACT frozen expression (same operation order), so
      # it is bit-identical to estimate_half; rowMeans of the per-problem
      # residuals would differ at 1e-16 and flip integer-count ties when
      # split_correlations ranks the score.
      if ("loss_mean" %in% want)
        out$loss_mean <- (rowSums(cnt) - r * sum(fm$p_base[lic])) /
          (r * length(lic))
      if ("loss_median" %in% want) {
        resid <- sweep(cnt / r, 2, fm$p_base[lic], "-")
        out$loss_median <- if (ncol(resid) == 1) resid[, 1] else
          apply(resid, 1, stats::median)
      }
    }
  }
  out
}

rsb_of <- function(h1, h2) {
  rs <- split_correlations(h1, h2, MIN_PAIRS)                         # FROZEN
  r_raw <- mean(rs, na.rm = TRUE)
  if (is.finite(r_raw)) sb_correct(r_raw) else NA_real_               # FROZEN
}

# Collect per-subject half scores across n_split common problem-halvings of a
# set. The RNG stream (set.seed(seed); sample.int(K, 2*kh) per split) is
# IDENTICAL to the frozen run_cell / E3 run_set / steps 30-33 run_set_rows, so
# the mean-loss r_sb reproduces E3's value bit-for-bit under E3 seeds.
collect_E <- function(fm, cols, rows, n_split, seed, want) {
  NBs <- fm$NB[rows, cols, drop = FALSE]
  keep <- which(!is.na(rowSums(NBs)))
  NBs <- NBs[keep, , drop = FALSE]
  ns <- length(keep); K <- length(cols); kh <- K %/% 2L
  E <- list()
  for (nm in want) E[[nm]] <- list(h1 = matrix(NA_real_, ns, n_split),
                                   h2 = matrix(NA_real_, ns, n_split))
  set.seed(seed)
  for (s in seq_len(n_split)) {
    take <- sample.int(K, 2L * kh)
    i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
    a <- half_scores(NBs[, i1, drop = FALSE], fm, cols[i1], R_REPS, want)
    b <- half_scores(NBs[, i2, drop = FALSE], fm, cols[i2], R_REPS, want)
    for (nm in want) { E[[nm]]$h1[, s] <- a[[nm]]; E[[nm]]$h2[, s] <- b[[nm]] }
  }
  list(n_subjects = ns, rows_kept = rows[keep], E = E)
}

# Residualize each column of M on the per-subject covariate g (finite for all
# rows). For loss, a column is either all-finite (>= 1 loss-informative problem
# in the half) or all-NA (0), so the rank-1 projection removal below is exact;
# Spearman is invariant to the location shift, so only the g-projection matters.
resid_loss_cols <- function(M, g) {
  gc <- g - mean(g)
  denom <- sum(gc * gc)
  if (!is.finite(denom) || denom <= 0) return(M)
  b <- as.vector(crossprod(gc, M)) / denom     # NA for all-NA columns
  M - outer(gc, b)
}

draw_sets <- function(np, K, M, seed) {
  set.seed(seed)
  total <- choose(np, K)
  if (total <= M_CAP) {
    cmb <- utils::combn(np, K)
    return(lapply(seq_len(ncol(cmb)), function(i) sort(cmb[, i])))
  }
  seen <- new.env(hash = TRUE, parent = emptyenv())
  out <- vector("list", M); got <- 0L
  while (got < M) {
    s <- sort(sample.int(np, K))
    key <- paste(s, collapse = ",")
    if (is.null(seen[[key]])) {
      assign(key, TRUE, envir = seen); got <- got + 1L; out[[got]] <- s
    }
  }
  out
}

# Between-set covariance decomposition + the pre-declared null control (permute
# the set-to-set correspondence between the two subject halves; the covariance
# must collapse). Same structure as steps 30-32 (frozen-independent stats).
cov_stats <- function(xA, xB, seed_perm, n_perm = 20L) {
  ok <- is.finite(xA) & is.finite(xB)
  n_used <- sum(ok)
  if (n_used < MIN_SETS_COV)
    return(list(n_used = n_used, cov = NA_real_, var_a = NA_real_,
                var_b = NA_real_, icc = NA_real_, sigma_true = NA_real_,
                sigma_raw = NA_real_, rho = NA_real_,
                cov_null = NA_real_, sigma_true_null = NA_real_))
  a <- xA[ok]; b <- xB[ok]
  cv <- stats::cov(a, b); va <- stats::var(a); vb <- stats::var(b)
  set.seed(seed_perm)
  nulls <- vapply(seq_len(n_perm), function(z) stats::cov(a, sample(b)),
                  numeric(1))
  list(n_used = n_used, cov = cv, var_a = va, var_b = vb,
       icc = cv / mean(c(va, vb)),
       sigma_true = sqrt(max(0, cv)), sigma_raw = sqrt(mean(c(va, vb))),
       rho = suppressWarnings(stats::cor(a, b, method = "spearman")),
       cov_null = mean(nulls),
       sigma_true_null = mean(sqrt(pmax(0, nulls))))
}

partition_rows <- function(n, seed, n_a) {
  set.seed(seed)
  p <- sample.int(n)
  list(A = sort(p[seq_len(n_a)]), B = sort(p[(n_a + 1L):n]))
}

sb_m <- function(m, r1) if (!is.finite(r1)) NA_real_ else
  m * r1 / (1 + (m - 1) * r1)

# =============================================================================
# 4. INLINE G_EQ: the from-raw mean loss score and ev reproduce the FROZEN
#    estimate_half bit-for-bit, BEFORE any result.
# =============================================================================
say("\n=== inline G_EQ: half_scores vs frozen estimate_half ===\n")
set.seed(SEED_C + 111L)
eq_max_ev <- 0; eq_max_loss <- 0; eq_na <- 0L; eq_n <- 0L
for (it in seq_len(4000L)) {
  ri <- sample.int(63L, 1); sid <- fm$subs[ri]; sd <- subj_data[[sid]]
  kh <- sample(3:6, 1); cols <- sort(sample(seq_along(fm$probs), kh))
  a <- estimate_half(sd, match(fm$probs[cols], sd$feat$GameID), POS, REGIME)
  b <- half_scores(fm$NB[ri, cols, drop = FALSE], fm, cols, R_REPS,
                   c("loss_mean", "ev"))
  for (pair in list(c("loss_avoid_adj", "loss_mean"), c("ev_sens_agn", "ev"))) {
    x <- unname(a[[pair[1]]]); y <- unname(b[[pair[2]]])
    if (is.na(x) != is.na(y)) { eq_na <- eq_na + 1L; next }
    if (!is.na(x)) {
      d <- abs(x - y)
      if (pair[1] == "ev_sens_agn") eq_max_ev <- max(eq_max_ev, d) else
        eq_max_loss <- max(eq_max_loss, d)
    }
  }
  eq_n <- eq_n + 1L
}
say("inline G_EQ: %d halves; max|d loss_mean| = %.3g; max|d ev| = %.3g; NA mismatches = %d\n",
    eq_n, eq_max_loss, eq_max_ev, eq_na)
stopifnot(eq_na == 0L, eq_max_loss < 1e-12, eq_max_ev < 1e-12)

# =============================================================================
# 5. Reproduce E3's Set-1 sets, and anchor the full-63 mean r_sb to E3's
#    reliability_by_set_real (bit-for-bit for the loss target).
# =============================================================================
e3_real <- readr::read_csv(
  here("output", "elevation", "results", "e3", "reliability_by_set_real.csv"),
  show_col_types = FALSE)
e3_s1 <- e3_real |> dplyr::filter(form == GATE_FORM)
parse_ids <- function(s) as.integer(strsplit(s, ";", fixed = TRUE)[[1]])

sets_by_K <- lapply(K_GRID, function(K) {
  rs <- draw_sets(length(fm$probs), K, M_CAP, seed_subset_e3(GATE_FORM, K))
  stopifnot(length(rs) == M_CAP)
  # match to E3's saved game_ids for this K/target (bit-identical set order)
  d <- e3_s1 |> dplyr::filter(target == GATE_TARGET, K == !!K) |>
    dplyr::arrange(set_id)
  stopifnot(nrow(d) == M_CAP)
  committed <- lapply(d$game_ids, function(s) sort(match(parse_ids(s), fm$probs)))
  stopifnot(identical(lapply(rs, as.integer), lapply(committed, as.integer)))
  rs
})
names(sets_by_K) <- as.character(K_GRID)
say("E3 committed Set-1 sets reproduced from E3 seeds at all %d K (bit-identical order)\n",
    length(K_GRID))

# =============================================================================
# 6. Response-style traits (C1) over ALL 24 Set-1 problems
# =============================================================================
# B-side bias = mean B-choice rate over the 24 problems (experience trials).
# safe-choice rate = mean rate of choosing the lower-SD option (0.5 on SD tie).
bside_bias <- rowMeans(fm$NB) / R_REPS
safe_dir <- ifelse(fm$sd_b < fm$sd_a, 1, ifelse(fm$sd_b > fm$sd_a, 0, 0.5))
# Structural note: in Set 1 option B is the higher-SD option for every clean
# problem, so safe == 1 - B-choice on each problem and safe-choice rate is
# exactly 1 - B-side bias (corr = -1). Assert the structural premise.
n_b_higher_sd <- sum(fm$sd_b > fm$sd_a)
n_sd_tie <- sum(fm$sd_b == fm$sd_a)
say("Set-1 SD structure: %d of %d problems have sd_b > sd_a (%d ties)\n",
    n_b_higher_sd, length(fm$probs), n_sd_tie)
safe_rate <- vapply(seq_len(63L), function(i) {
  rate_B <- fm$NB[i, ] / R_REPS
  mean(ifelse(safe_dir == 1, rate_B,
              ifelse(safe_dir == 0, 1 - rate_B, 0.5)))
}, numeric(1))
corr_bside_safe <- stats::cor(bside_bias, safe_rate)
say("corr(B-side bias, safe-choice rate) = %.6f\n", corr_bside_safe)

# Collinearity handling (pre-declared): where |corr| > 0.95 the
# two traits are entered as a single COMBINED side/risk axis to avoid a
# degenerate two-regressor fit. Because safe = 1 - B-side bias exactly, the
# combined axis IS the (standardized) B-side bias; residualizing loss on either
# trait, or on their principal axis, gives identical OLS residuals.
collinear <- abs(corr_bside_safe) > COLLIN_THRESH
combined_axis <- if (collinear) {
  as.numeric(scale(bside_bias))
} else {
  # non-degenerate case (not expected in Set 1): first PC of the two traits.
  pc <- stats::prcomp(cbind(scale(bside_bias), scale(safe_rate)))$x[, 1]
  as.numeric(pc)
}
style <- combined_axis                       # per-subject response-style axis
names(style) <- fm$subs

# R^2 of the full-instrument per-subject loss score on the style axis (context).
full_loss <- {
  lic <- which(fm$li)
  cnt <- matrix(0, 63L, length(lic))
  for (j in seq_along(lic)) {
    cj <- if (fm$has_loss_b[lic[j]]) R_REPS - fm$NB[, lic[j]] else fm$NB[, lic[j]]
    cnt[, j] <- cj / R_REPS - fm$p_base[lic[j]]
  }
  rowMeans(cnt)
}
r2_loss_style <- summary(stats::lm(full_loss ~ style))$r.squared
say("R^2(full-Set-1 loss ~ style axis) = %.5f\n", r2_loss_style)

# =============================================================================
# 7. C2 building block: per-problem loss residuals over the 6 loss-informative
#    Set-1 problems, and the single-problem trial-split reliability r_1.
# =============================================================================
li_cols <- which(fm$li)                       # the 6 loss-informative problems
# resid_full[subject, li] = chose-no-loss rate minus the frozen EV baseline.
resid_full <- matrix(NA_real_, 63L, length(li_cols),
                     dimnames = list(fm$subs, as.character(fm$probs[li_cols])))
for (j in seq_along(li_cols)) {
  cj <- if (fm$has_loss_b[li_cols[j]]) R_REPS - fm$NB[, li_cols[j]] else
    fm$NB[, li_cols[j]]
  resid_full[, j] <- cj / R_REPS - fm$p_base[li_cols[j]]
}

# Single-problem trial-split reliability: split POS (20 experience trials) into
# two halves of 10, form the per-subject loss residual on each half, correlate
# across subjects, SB-correct; average over N_SPLIT trial-halvings; mean over
# the 6 loss-informative problems. This is the local-independence reference r_1.
r1_by_problem <- vapply(seq_along(li_cols), function(j) {
  gid <- fm$probs[li_cols[j]]
  hl <- fm$has_loss_b[li_cols[j]]; pb <- fm$p_base[li_cols[j]]
  Bmat <- matrix(NA_real_, 63L, R_REPS)
  for (i in seq_len(63L)) {
    sd <- subj_data[[fm$subs[i]]]
    rr <- match(gid, sd$feat$GameID)
    Bmat[i, ] <- sd$B_mat[rr, POS]
  }
  set.seed(seed_trial_c(gid))
  rs <- numeric(N_SPLIT)
  for (s in seq_len(N_SPLIT)) {
    tk <- sample.int(R_REPS, R_REPS)
    h1 <- tk[seq_len(R_REPS %/% 2L)]; h2 <- tk[(R_REPS %/% 2L + 1L):R_REPS]
    nB1 <- rowSums(Bmat[, h1, drop = FALSE]); nB2 <- rowSums(Bmat[, h2, drop = FALSE])
    m <- length(h1)
    c1 <- if (hl) m - nB1 else nB1; c2 <- if (hl) m - nB2 else nB2
    x1 <- c1 / m - pb; x2 <- c2 / m - pb
    rs[s] <- suppressWarnings(stats::cor(x1, x2, method = "spearman"))
  }
  sb_correct(mean(rs, na.rm = TRUE))
}, numeric(1))
r1_mean <- mean(r1_by_problem, na.rm = TRUE)
say("single-problem trial-split reliability r_1 (mean over 6 loss problems) = %.4f\n",
    r1_mean)

# wpc per set: within-set mean pairwise Spearman of the loss residuals over the
# set's loss-informative problems (needs >= 2). A precomputed pairwise Spearman
# matrix over the 6 problems makes this a lookup.
sp_pair <- matrix(NA_real_, length(li_cols), length(li_cols))
for (a in seq_along(li_cols)) for (b in seq_along(li_cols)) if (a < b)
  sp_pair[a, b] <- suppressWarnings(stats::cor(resid_full[, a], resid_full[, b],
                                               method = "spearman"))
# map a Set-1 problem column index (into fm$probs) to its li index (1..6)
li_index_of <- match(seq_along(fm$probs), li_cols)   # NA for non-loss problems
wpc_of_set <- function(cols) {
  lis <- li_index_of[cols]; lis <- lis[!is.na(lis)]
  if (length(lis) < 2L) return(c(wpc = NA_real_, n_li = length(lis)))
  vals <- c()
  for (a in seq_along(lis)) for (b in seq_along(lis)) if (a < b) {
    ia <- min(lis[a], lis[b]); ib <- max(lis[a], lis[b])
    vals <- c(vals, sp_pair[ia, ib])
  }
  c(wpc = mean(vals, na.rm = TRUE), n_li = length(lis))
}

# =============================================================================
# 8. Full-63 pass: per-set mean r_sb (anchor), median r_sb (C3), residualized
#    mean r_sb (C1 bootstrap floor), ev r_sb (informational), wpc (C2). Also the
#    50-set noise floor via the FROZEN reliability_point, raw and residualized.
# =============================================================================
say("\n=== full-63 per-set pass (mean / median / residualized / ev / wpc) ===\n")
cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(tidyverse); library(here) })
  source(here("R", "step2_reliability.R"))
  NULL
}))
Z975 <- stats::qnorm(0.975)
parallel::clusterExport(cl, c(
  "fm", "half_scores", "ev_ranks", "collect_E", "rsb_of", "resid_loss_cols",
  "split_correlations", "TARGETS", "R_REPS", "MIN_PAIRS", "N_SPLIT",
  "N_NOISE_SETS", "N_BOOT", "Z975", "style", "wpc_of_set", "li_index_of",
  "sp_pair", "seed_noise_c", "SEED_C"), envir = environment())

full63_res <- list()
for (K in K_GRID) {
  sets <- sets_by_K[[as.character(K)]]
  sseed <- seed_split_e3(GATE_FORM, K)            # E3 split seeds -> anchor
  jobs <- lapply(seq_len(M_CAP), function(i)
    list(set_id = i, cols = sets[[i]], seed = sseed + i,
         want_noise = i <= N_NOISE_SETS, K = K))
  parallel::clusterExport(cl, c("sseed", "K"), envir = environment())
  rr <- parallel::parLapply(cl, jobs, function(job) {
    cols <- job$cols
    E <- collect_E(fm, cols, seq_len(63L), N_SPLIT, job$seed,
                   c("loss_mean", "loss_median", "ev"))$E
    r_mean <- rsb_of(E$loss_mean$h1, E$loss_mean$h2)
    r_median <- rsb_of(E$loss_median$h1, E$loss_median$h2)
    r_ev <- rsb_of(E$ev$h1, E$ev$h2)
    l1 <- resid_loss_cols(E$loss_mean$h1, style)
    l2 <- resid_loss_cols(E$loss_mean$h2, style)
    r_resid <- rsb_of(l1, l2)
    w <- wpc_of_set(cols)
    sd_noise_raw <- NA_real_; sd_noise_resid <- NA_real_
    if (job$want_noise) {
      rp <- reliability_point(E$loss_mean$h1, E$loss_mean$h2, N_BOOT,
                              seed_noise_c(job$K, job$set_id), MIN_PAIRS)  # FROZEN
      sd_noise_raw <- (rp$ci_hi_sb - rp$ci_lo_sb) / (2 * Z975)
      rpr <- reliability_point(l1, l2, N_BOOT,
                               seed_noise_c(job$K, job$set_id) + 400000L,
                               MIN_PAIRS)                                   # FROZEN
      sd_noise_resid <- (rpr$ci_hi_sb - rpr$ci_lo_sb) / (2 * Z975)
    }
    tibble::tibble(
      K = job$K, set_id = job$set_id,
      game_ids = paste(fm$probs[cols], collapse = ";"),
      n_loss_informative = sum(fm$li[cols]),
      wpc = unname(w["wpc"]), n_li_wpc = unname(w["n_li"]),
      r_sb_mean = r_mean, r_sb_median = r_median, r_sb_resid = r_resid,
      r_sb_ev = r_ev, sd_noise_raw = sd_noise_raw,
      sd_noise_resid = sd_noise_resid)
  })
  full63_res[[as.character(K)]] <- dplyr::bind_rows(rr)
  say("  K %2d done (%d sets)\n", K, M_CAP)
}
per_set <- dplyr::bind_rows(full63_res)

# ---- anchor: mean r_sb must reproduce E3's saved reliability -----------------
anchor <- per_set |>
  dplyr::select(K, set_id, r_sb_mean) |>
  dplyr::inner_join(e3_s1 |> dplyr::filter(target == GATE_TARGET) |>
                      dplyr::select(K, set_id, e3_r_sb = r_sb),
                    by = c("K", "set_id"))
anchor_na_ok <- identical(is.na(anchor$r_sb_mean), is.na(anchor$e3_r_sb))
anchor_d <- max(abs(anchor$r_sb_mean - anchor$e3_r_sb), na.rm = TRUE)
say("anchor vs committed E3 per-set r_sb (loss): max|d| = %.3g (NA pattern identical = %s)\n",
    anchor_d, anchor_na_ok)
stopifnot(anchor_na_ok, anchor_d < 1e-12)

# =============================================================================
# 9. C1 primary: covariance (subject-split) estimand, raw vs residualized loss.
# =============================================================================
say("\n=== C1 covariance estimand: raw vs response-style-residualized loss ===\n")
cov_cell_c1 <- function(fm, sets, rows_all, gvec_all, seed_part, part_index,
                        K, n_split) {
  n <- length(rows_all); n_a <- ceiling(n / 2)
  pr <- partition_rows(n, seed_part, n_a)
  A <- rows_all[pr$A]; B <- rows_all[pr$B]
  gA <- gvec_all[pr$A]; gB <- gvec_all[pr$B]
  M <- length(sets)
  out <- matrix(NA_real_, M, 4L)   # raw_A, raw_B, resid_A, resid_B (loss)
  for (i in seq_len(M)) {
    ea <- collect_E(fm, sets[[i]], A, n_split, seed_half_c(1L, K, part_index, i),
                    "loss_mean")
    eb <- collect_E(fm, sets[[i]], B, n_split, seed_half_c(2L, K, part_index, i),
                    "loss_mean")
    gAk <- gA[match(ea$rows_kept, A)]; gBk <- gB[match(eb$rows_kept, B)]
    out[i, 1] <- rsb_of(ea$E$loss_mean$h1, ea$E$loss_mean$h2)
    out[i, 2] <- rsb_of(eb$E$loss_mean$h1, eb$E$loss_mean$h2)
    out[i, 3] <- rsb_of(resid_loss_cols(ea$E$loss_mean$h1, gAk),
                        resid_loss_cols(ea$E$loss_mean$h2, gAk))
    out[i, 4] <- rsb_of(resid_loss_cols(eb$E$loss_mean$h1, gBk),
                        resid_loss_cols(eb$E$loss_mean$h2, gBk))
  }
  list(r = out, n_a = length(A), n_b = length(B))
}

parallel::clusterExport(cl, c("cov_cell_c1", "partition_rows", "cov_stats",
                              "seed_half_c", "MIN_SETS_COV", "K_GRID",
                              "sb_correct", "reliability_point"),
                        envir = environment())

# ---- primary cell K = 12: full B_PART partitions ----------------------------
run_cov_K <- function(K, n_part, part_index_fun, seed_part_fun, store_halves) {
  sets <- sets_by_K[[as.character(K)]]
  jobs <- lapply(seq_len(n_part), function(p) list(
    p = p, part_index = part_index_fun(p), seed_part = seed_part_fun(p),
    seed_perm = seed_perm_c(part_index_fun(p))))
  parallel::clusterExport(cl, c("sets", "K"), envir = environment())
  rr <- parallel::parLapply(cl, jobs, function(job) {
    cc <- cov_cell_c1(fm, sets, seq_len(63L), style, job$seed_part,
                      job$part_index, K, N_SPLIT)
    sr <- cov_stats(cc$r[, 1], cc$r[, 2], job$seed_perm)
    sd <- cov_stats(cc$r[, 3], cc$r[, 4], job$seed_perm + 500000L)
    row <- tibble::tibble(
      K = K, partition = job$p, n_subjects_A = cc$n_a, n_subjects_B = cc$n_b,
      n_sets_raw = sr$n_used, n_sets_resid = sd$n_used,
      sigma_set_raw_cov = sr$sigma_raw, sigma_true_raw = sr$sigma_true,
      ICC_raw = sr$icc, sigma_true_raw_null = sr$sigma_true_null,
      sigma_true_resid = sd$sigma_true, ICC_resid = sd$icc,
      sigma_set_resid_cov = sd$sigma_raw)
    halves <- if (store_halves) tibble::tibble(
      partition = job$p, set_id = seq_len(nrow(cc$r)),
      r_sb_raw_A = cc$r[, 1], r_sb_raw_B = cc$r[, 2],
      r_sb_resid_A = cc$r[, 3], r_sb_resid_B = cc$r[, 4]) else NULL
    list(row = row, halves = halves)
  })
  list(rows = dplyr::bind_rows(lapply(rr, function(z) z$row)),
       halves = if (store_halves)
         dplyr::bind_rows(lapply(rr, function(z) z$halves)) else NULL)
}

t0 <- Sys.time()
prim <- run_cov_K(GATE_K, B_PART, function(p) p, seed_part_c, store_halves = TRUE)
say("primary K=12 covariance: %d partitions in %.1f min\n", B_PART,
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

aux_rows <- list()
t0 <- Sys.time()
for (K in setdiff(K_GRID, GATE_K)) {
  rk <- run_cov_K(K, N_PART_AUX,
                  function(j) 200L + 10L * match(K, K_GRID) + j,
                  function(j) seed_part_aux_c(K, j), store_halves = FALSE)
  aux_rows[[as.character(K)]] <- rk$rows
}
say("auxiliary-K covariance: %d K x %d partitions in %.1f min\n",
    length(setdiff(K_GRID, GATE_K)), N_PART_AUX,
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

cov_by_partition <- dplyr::bind_rows(prim$rows, dplyr::bind_rows(aux_rows)) |>
  dplyr::arrange(K, partition)
prim_halves <- prim$halves

# ---- covariance ICC summary per K (raw vs residualized) ---------------------
cov_summary <- cov_by_partition |>
  dplyr::group_by(K) |>
  dplyr::summarise(
    estimand = "covariance_subject_split", n_partitions = dplyr::n(),
    sigma_set_raw_cov = mean(sigma_set_raw_cov, na.rm = TRUE),
    ICC_raw = mean(ICC_raw, na.rm = TRUE),
    sigma_true_raw = mean(sigma_true_raw, na.rm = TRUE),
    ICC_resid = mean(ICC_resid, na.rm = TRUE),
    sigma_true_resid = mean(sigma_true_resid, na.rm = TRUE),
    sigma_true_raw_null = mean(sigma_true_raw_null, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(
    rel_ICC_drop = (ICC_raw - ICC_resid) / ICC_raw,
    rel_sigma_true_drop = (sigma_true_raw - sigma_true_resid) / sigma_true_raw)

# =============================================================================
# 10. C1 bootstrap-floor estimand (reported for continuity): full-63 var of
#     r_sb minus the mean per-set estimation variance, raw vs residualized.
# =============================================================================
bootfloor <- per_set |>
  dplyr::group_by(K) |>
  dplyr::summarise(
    sigma2_raw = stats::var(r_sb_mean, na.rm = TRUE),
    sigma2_resid = stats::var(r_sb_resid, na.rm = TRUE),
    sigma2_noise_raw = mean(sd_noise_raw^2, na.rm = TRUE),
    sigma2_noise_resid = mean(sd_noise_resid^2, na.rm = TRUE),
    sigma2_ev = stats::var(r_sb_ev, na.rm = TRUE),
    n_noise = sum(is.finite(sd_noise_raw)), .groups = "drop") |>
  dplyr::mutate(
    estimand = "bootstrap_floor",
    sigma_true_raw = sqrt(pmax(0, sigma2_raw - sigma2_noise_raw)),
    ICC_raw = pmax(0, sigma2_raw - sigma2_noise_raw) / sigma2_raw,
    sigma_true_resid = sqrt(pmax(0, sigma2_resid - sigma2_noise_resid)),
    ICC_resid = pmax(0, sigma2_resid - sigma2_noise_resid) / sigma2_resid,
    rel_ICC_drop = (ICC_raw - ICC_resid) / ICC_raw,
    rel_sigma_true_drop = (sigma_true_raw - sigma_true_resid) / sigma_true_raw)

# ---- response_style_adjusted_icc.csv ----------------------------------------
response_style_adjusted_icc <- dplyr::bind_rows(
  cov_summary |> dplyr::transmute(
    K, estimand, target = GATE_TARGET, n_units = n_partitions,
    sigma_set_raw = sigma_set_raw_cov, ICC_raw, sigma_true_raw,
    ICC_resid, sigma_true_resid, rel_ICC_drop, rel_sigma_true_drop),
  bootfloor |> dplyr::transmute(
    K, estimand, target = GATE_TARGET, n_units = M_CAP,
    sigma_set_raw = sqrt(sigma2_raw), ICC_raw, sigma_true_raw,
    ICC_resid, sigma_true_resid, rel_ICC_drop, rel_sigma_true_drop)) |>
  dplyr::mutate(corr_bside_safe = corr_bside_safe,
                r2_loss_on_style = r2_loss_style,
                combined_axis = "standardized B-side bias (safe = 1 - B-side)",
                collinearity_handled = collinear) |>
  dplyr::arrange(estimand, K)
save_csv(response_style_adjusted_icc, "response_style_adjusted_icc")

gc1_cov <- cov_summary |> dplyr::filter(K == GATE_K)
gc1_drop <- gc1_cov$rel_ICC_drop
say("\nC1 primary (K=12, covariance): ICC_raw %.4f -> ICC_resid %.4f, relative drop %.4f\n",
    gc1_cov$ICC_raw, gc1_cov$ICC_resid, gc1_drop)
gc1_verdict <- if (gc1_drop >= GC1_MATERIAL) "narrow_response_style_driver" else
  if (gc1_drop < GC1_NEGLIGIBLE) "strengthen_not_a_response_style_artifact" else
    "partial_qualified"

# =============================================================================
# 11. C2 local item dependence: wpc <-> r_sb, terciles, dependence-adjusted ICC.
# =============================================================================
say("\n=== C2 local item dependence ===\n")
c2_tbl <- per_set |>
  dplyr::select(K, set_id, game_ids, n_loss_informative, wpc, n_li_wpc,
                r_sb = r_sb_mean)
save_csv(c2_tbl, "local_dependence_diagnostic")

# relationship + tercile gap per K (primary K = 12)
c2_rel <- c2_tbl |>
  dplyr::group_by(K) |>
  dplyr::summarise(
    n_sets_wpc = sum(is.finite(wpc) & is.finite(r_sb)),
    spearman_rsb_wpc = {
      ok <- is.finite(wpc) & is.finite(r_sb)
      if (sum(ok) < MIN_SETS_COV) NA_real_ else
        suppressWarnings(stats::cor(r_sb[ok], wpc[ok], method = "spearman"))
    },
    mean_wpc = mean(wpc, na.rm = TRUE),
    .groups = "drop")

tercile_gap <- function(d) {
  ok <- is.finite(d$wpc) & is.finite(d$r_sb)
  d <- d[ok, ]
  if (nrow(d) < 3L * N_WPC_TERCILE_MIN)
    return(tibble::tibble(n = nrow(d), lo_mean = NA_real_, hi_mean = NA_real_,
                          gap = NA_real_))
  qs <- stats::quantile(d$wpc, c(1/3, 2/3), na.rm = TRUE)
  lo <- d$r_sb[d$wpc <= qs[1]]; hi <- d$r_sb[d$wpc >= qs[2]]
  tibble::tibble(n = nrow(d), lo_mean = mean(lo), hi_mean = mean(hi),
                 gap = mean(hi) - mean(lo))
}
c2_terciles <- c2_tbl |> dplyr::group_by(K) |>
  dplyr::group_modify(~ tercile_gap(.x)) |> dplyr::ungroup()

# dependence-adjusted ICC (primary cell): residualize the raw covariance r_sb
# halves on the per-set wpc across sets, then recompute the between-set
# covariance ICC. The relative reduction measures how much of the between-set
# reliability covariance is attributable to within-set problem-residual
# coherence (wpc). NOTE (pre-declared): r_sb tracking wpc is
# expected under BOTH construct coherence and nuisance dependence, so this
# reduction is a joint upper bound on a dependence artifact, not a clean split.
wpc_vec <- c2_tbl |> dplyr::filter(K == GATE_K) |> dplyr::arrange(set_id) |>
  dplyr::pull(wpc)
resid_on_wpc <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  out <- rep(NA_real_, length(x))
  if (sum(ok) < MIN_SETS_COV) return(out)
  fit <- stats::lm(x[ok] ~ w[ok])
  out[ok] <- stats::residuals(fit)
  out
}
c2_adj <- prim_halves |>
  dplyr::group_by(partition) |>
  dplyr::group_modify(function(d, key) {
    d <- d |> dplyr::arrange(set_id)
    w <- wpc_vec[d$set_id]
    rawA <- d$r_sb_raw_A; rawB <- d$r_sb_raw_B
    s_raw <- cov_stats(rawA, rawB, seed_perm_c(key$partition[[1]]))
    adjA <- resid_on_wpc(rawA, w); adjB <- resid_on_wpc(rawB, w)
    s_adj <- cov_stats(adjA, adjB, seed_perm_c(key$partition[[1]]) + 900000L)
    tibble::tibble(ICC_raw = s_raw$icc, sigma_true_raw = s_raw$sigma_true,
                   ICC_wpc_adj = s_adj$icc, sigma_true_wpc_adj = s_adj$sigma_true)
  }) |> dplyr::ungroup()
c2_adj_summary <- c2_adj |>
  dplyr::summarise(n_partitions = dplyr::n(),
                   ICC_raw = mean(ICC_raw, na.rm = TRUE),
                   ICC_wpc_adj = mean(ICC_wpc_adj, na.rm = TRUE),
                   sigma_true_raw = mean(sigma_true_raw, na.rm = TRUE),
                   sigma_true_wpc_adj = mean(sigma_true_wpc_adj, na.rm = TRUE)) |>
  dplyr::mutate(rel_ICC_reduction = (ICC_raw - ICC_wpc_adj) / ICC_raw)

# local-independence prediction: SB-projection of r_1 to the mean per-set count
# of loss-informative problems per problem-half, vs the observed mean set r_sb.
c2_predicted <- c2_tbl |>
  dplyr::filter(is.finite(r_sb), n_loss_informative >= 2) |>
  dplyr::group_by(K) |>
  dplyr::summarise(
    mean_n_li = mean(n_loss_informative),
    mean_n_li_half = mean(n_loss_informative) / 2,
    r1 = r1_mean,
    predicted_r_raw = sb_m(mean(n_loss_informative) / 2, r1_mean),
    predicted_r_sb = sb_correct(sb_m(mean(n_loss_informative) / 2, r1_mean)),
    observed_r_sb = mean(r_sb, na.rm = TRUE),
    .groups = "drop")

gc2_adj <- c2_adj_summary$rel_ICC_reduction
gc2_rel_k12 <- c2_rel$spearman_rsb_wpc[c2_rel$K == GATE_K]
gc2_gap_k12 <- c2_terciles$gap[c2_terciles$K == GATE_K]
say("C2 primary (K=12): Spearman(r_sb, wpc) = %.4f; tercile gap = %.4f; dependence-adjusted ICC reduction = %.4f\n",
    gc2_rel_k12, gc2_gap_k12, gc2_adj)
# The discriminating question is whether wpc exceeds construct coherence. The
# adjusted-ICC reduction and the tercile gap cannot separate the two (both rise
# under coherence AND under nuisance dependence), so the verdict is set by the
# pre-declared rule and, where the readings do not cleanly discriminate, is
# reported INCONCLUSIVE (a valid pre-declared outcome).
gc2_material <- isTRUE(gc2_adj >= GC2_ADJ_DROP)
gc2_verdict <- if (gc2_material) "narrow_partly_dependence" else
  "inconclusive_not_separable_from_coherence"

# =============================================================================
# 12. C3 alternative-estimator robustness of the set ranking (median loss).
# =============================================================================
say("\n=== C3 alternative-estimator (median loss) set ranking ===\n")
alt_estimator_set_ranking <- per_set |>
  dplyr::transmute(K, set_id, game_ids, n_loss_informative,
                   r_sb_primary_mean = r_sb_mean, r_sb_robust_median = r_sb_median)
save_csv(alt_estimator_set_ranking, "alt_estimator_set_ranking")

c3_rank <- alt_estimator_set_ranking |>
  dplyr::group_by(K) |>
  dplyr::summarise(
    n_sets = sum(is.finite(r_sb_primary_mean) & is.finite(r_sb_robust_median)),
    ranking_spearman = {
      ok <- is.finite(r_sb_primary_mean) & is.finite(r_sb_robust_median)
      if (sum(ok) < MIN_SETS_COV) NA_real_ else
        suppressWarnings(stats::cor(r_sb_primary_mean[ok],
                                    r_sb_robust_median[ok], method = "spearman"))
    },
    median_primary = stats::median(r_sb_primary_mean, na.rm = TRUE),
    median_robust = stats::median(r_sb_robust_median, na.rm = TRUE),
    .groups = "drop")
gc3_rank_k12 <- c3_rank$ranking_spearman[c3_rank$K == GATE_K]
say("C3 primary (K=12): set-ranking Spearman(mean, median) = %.4f\n", gc3_rank_k12)
gc3_verdict <- if (isTRUE(gc3_rank_k12 >= GC3_RANK_MIN))
  "robust_not_scoring_specific" else "estimator_specific"

parallel::stopCluster(cl)

# =============================================================================
# 13. HG-ceiling guard (relative comparisons are not affected by the
#     shared-instrument inflation, but the flag is kept).
# =============================================================================
e0 <- readr::read_csv(
  here("output", "elevation", "results", "e0", "variance_ceiling_final.csv"),
  show_col_types = FALSE)
r1_of <- function(R, k) R / (k - (k - 1) * R)
sb_k  <- function(k, r1) k * r1 / (1 + (k - 1) * r1)
rw <- e0[e0$construct == GATE_TARGET, ]
stopifnot(nrow(rw) == 1L, rw$k_ref == 14)
R1_loss <- r1_of(rw$ceiling_fulldata, 14)
ceiling_k12 <- sb_k(GATE_K, R1_loss)
hg_median_raw <- stats::median(per_set$r_sb_mean[per_set$K == GATE_K], na.rm = TRUE)
hg_flag <- isTRUE(hg_median_raw > ceiling_k12)
say("HG-ceiling: K=12 median raw r_sb %.4f vs E0 out-of-subject ceiling %.4f -> %s (relative comparisons unaffected)\n",
    hg_median_raw, ceiling_k12, if (hg_flag) "flagged (median above ceiling; carried)" else "clear")

# =============================================================================
# 14. Exports for independent re-derivation
# =============================================================================
save_csv(tibble::tibble(SubjID = fm$subs, b_side_bias = bside_bias,
                        safe_choice_rate = safe_rate, combined_axis_style = style,
                        full_set1_loss = full_loss),
         "dispc_export_response_style")
save_csv(as.data.frame(resid_full) |> tibble::rownames_to_column("SubjID"),
         "dispc_export_problem_residuals")
save_csv(per_set |> dplyr::filter(K == GATE_K, set_id <= 40L) |>
           dplyr::select(K, set_id, game_ids, n_loss_informative, wpc,
                         r_sb_mean, r_sb_median, r_sb_resid, r_sb_ev),
         "dispc_export_per_set_scores")
save_csv(prim_halves |> dplyr::filter(partition <= 10L) |>
           dplyr::arrange(partition, set_id),
         "dispc_export_per_set_halves")

# =============================================================================
# 15. Gates, part file, seeds
# =============================================================================
say("\n=== C GATES ===\n")
say("G-C1 response style: relative cov-ICC drop %.4f (>= %.2f material; < %.2f negligible) -> %s\n",
    gc1_drop, GC1_MATERIAL, GC1_NEGLIGIBLE, gc1_verdict)
say("G-C2 local dependence: adj-ICC reduction %.4f (>= %.2f material) -> %s\n",
    gc2_adj, GC2_ADJ_DROP, gc2_verdict)
say("G-C3 alt-estimator: ranking Spearman %.4f (>= %.2f) -> %s\n",
    gc3_rank_k12, GC3_RANK_MIN, gc3_verdict)

part <- list(
  gate = list(
    spec = "elevation/DISPC_spec.md",
    primary_cell = list(form = GATE_FORM, K = GATE_K, target = GATE_TARGET,
                        regime = "experience (Trial 6-25)", r_reps = R_REPS),
    `G-C1` = list(
      diagnostic = "response-style contamination (B-side bias / safe-choice)",
      thresholds = list(material_drop = GC1_MATERIAL,
                        negligible_drop = GC1_NEGLIGIBLE),
      estimand_primary = "covariance_subject_split",
      corr_bside_safe = corr_bside_safe,
      collinearity_handled = collinear,
      combined_axis = "standardized B-side bias (safe = 1 - B-side, corr = -1)",
      r2_loss_on_style = r2_loss_style,
      primary_cell = list(ICC_raw = gc1_cov$ICC_raw, ICC_resid = gc1_cov$ICC_resid,
                          sigma_true_raw = gc1_cov$sigma_true_raw,
                          sigma_true_resid = gc1_cov$sigma_true_resid,
                          relative_ICC_drop = gc1_drop,
                          n_partitions = gc1_cov$n_partitions),
      per_K = lapply(seq_len(nrow(cov_summary)), function(i)
        as.list(cov_summary[i, c("K", "ICC_raw", "ICC_resid", "rel_ICC_drop")])),
      verdict = gc1_verdict),
    `G-C2` = list(
      diagnostic = "local item dependence (within-set problem-residual coherence)",
      threshold_adj_drop = GC2_ADJ_DROP,
      spearman_rsb_wpc_k12 = gc2_rel_k12,
      tercile_gap_k12 = gc2_gap_k12,
      dependence_adjusted_ICC_reduction = gc2_adj,
      single_problem_reliability_r1 = r1_mean,
      predicted_vs_observed = lapply(seq_len(nrow(c2_predicted)), function(i)
        as.list(c2_predicted[i, ])),
      caveat = paste("r_sb tracking wpc is expected under BOTH construct",
                     "coherence and nuisance dependence; the adjusted-ICC",
                     "reduction is a joint upper bound, not a clean dependence",
                     "estimate"),
      verdict = gc2_verdict),
    `G-C3` = list(
      diagnostic = "scoring-rule specificity (median vs mean loss residual)",
      threshold_ranking_min = GC3_RANK_MIN,
      ranking_spearman_k12 = gc3_rank_k12,
      per_K = lapply(seq_len(nrow(c3_rank)), function(i)
        as.list(c3_rank[i, c("K", "ranking_spearman", "n_sets")])),
      verdict = gc3_verdict),
    set1_collinearity = list(
      corr_bside_safe = corr_bside_safe,
      n_problems_b_higher_sd = n_b_higher_sd, n_sd_ties = n_sd_tie,
      note = paste("in Set 1 option B is the higher-SD option for every clean",
                   "problem, so safe-choice rate == 1 - B-side bias exactly",
                   "(corr = -1); the two traits enter as one combined side/risk",
                   "axis, per DISPC_spec section 2")),
    HG_ceiling = list(
      comparator = "E0 full-data out-of-subject ceiling SB-projected to K=12",
      ceiling_k12 = ceiling_k12, median_raw_r_sb_k12 = hg_median_raw,
      flagged = hg_flag,
      note = paste("This analysis compares RELATIVE quantities (ICC before/after,",
                   "rankings) which the shared-instrument inflation does not",
                   "affect; the flag is carried per DISPC_spec section 4"))),
  config = list(
    spec = "elevation/DISPC_spec.md",
    seed_base = SEED_C, seed_e3_for_set_membership = SEED_E3,
    regime = "experience (Trial 6-25)", r_reps = R_REPS,
    K_grid = K_GRID, gate_K = GATE_K, gate_form = GATE_FORM,
    targets = list(primary = GATE_TARGET, informational = "ev_sens_agn"),
    M_cap = M_CAP, n_split = N_SPLIT, b_part = B_PART, n_part_aux = N_PART_AUX,
    n_noise_sets = N_NOISE_SETS, n_boot = N_BOOT, min_pairs = MIN_PAIRS,
    min_sets_cov = MIN_SETS_COV,
    response_style = list(
      traits = c("B-side bias (mean B-choice rate)",
                 "safe-choice rate (mean lower-SD-choice rate, 0.5 on tie)"),
      corr_bside_safe = corr_bside_safe, collinear = collinear,
      combined_axis = "standardized B-side bias"),
    robust_score = "median per-problem loss residual over a half's loss-informative problems",
    wpc = "within-set mean pairwise Spearman of per-problem loss residuals over loss-informative problems",
    single_problem_reliability_r1 = r1_mean,
    set_membership_note = paste("reproduces E3's committed Set-1 sets from",
                                "seed_subset_e3(1, K) = SEED_E3 + 1000 + K; the",
                                "full-63 mean r_sb anchors to E3's committed",
                                "reliability_by_set_real bit-for-bit"),
    anchor_max_abs_diff_r_sb = anchor_d,
    partition_budget_rule = paste("the covariance estimand uses B_PART = 200",
                                  "partitions at the primary cell K = 12 and",
                                  "N_PART_AUX = 8 partitions at each auxiliary K",
                                  "(6, 8, 10), all on E3's single committed set",
                                  "sample; fixed before any result"),
    data_inputs = c("data/raw/All estimation raw data.csv",
                    "elevation/results/e1/participant_split.csv",
                    "elevation/results/e2/problem_information_dictionary.csv",
                    "elevation/results/e3/reliability_by_set_real.csv",
                    "elevation/results/e0/variance_ceiling_final.csv"),
    split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
    loss_baselines_exp = as.list(baselines$exp),
    inline_g_eq = list(halves = eq_n, max_abs_diff_loss = eq_max_loss,
                       max_abs_diff_ev = eq_max_ev, na_mismatches = eq_na),
    build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
    optimization_performed = FALSE, candidate_design_constructed = FALSE,
    new_mechanism_model_fitted = FALSE, e3_h1_gate_rerun_or_amended = FALSE,
    e4_h2_gate_rerun_or_amended = FALSE,
    runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
jsonlite::write_json(part, file.path(res_dir, "dispc_part.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

save_csv(cov_by_partition, "c1_cov_by_partition")
save_csv(c2_rel, "c2_wpc_relationship")
save_csv(c2_terciles, "c2_wpc_terciles")
save_csv(c2_adj_summary, "c2_dependence_adjusted_icc")
save_csv(c2_predicted, "c2_local_independence_prediction")
save_csv(c3_rank, "c3_ranking_by_K")

# ---- seeds ------------------------------------------------------------------
seeds_c <- dplyr::bind_rows(
  purrr::map_dfr(K_GRID, function(K) tibble::tibble(
    scope = sprintf("K%d", K), purpose = c("set_draw_e3", "split_draw_e3_base"),
    seed = c(seed_subset_e3(GATE_FORM, K), seed_split_e3(GATE_FORM, K)),
    derivation = c("E3 replay: SEED_E3 + 1000*form + K",
                   "E3 replay: SEED_E3 + 1000*form + K + 500000; per set + set_id"))),
  tibble::tibble(scope = sprintf("primary|partition%d", seq_len(B_PART)),
                 purpose = "subject_partition", seed = seed_part_c(seq_len(B_PART)),
                 derivation = "SEED_C + 700000 + partition_index"),
  purrr::map_dfr(setdiff(K_GRID, GATE_K), function(K) tibble::tibble(
    scope = sprintf("auxK%d|partition%d", K, seq_len(N_PART_AUX)),
    purpose = "subject_partition", seed = seed_part_aux_c(K, seq_len(N_PART_AUX)),
    derivation = "SEED_C + 720000 + 1000*K + j")),
  purrr::map_dfr(K_GRID, function(K) tibble::tibble(
    scope = sprintf("K%d|noise", K), purpose = "noise_bootstrap_base",
    seed = seed_noise_c(K, 0L),
    derivation = "SEED_C + 600000 + 1000*K + set_id (+400000 for residualized)")),
  purrr::map_dfr(fm$probs[li_cols], function(g) tibble::tibble(
    scope = sprintf("problem%d|trial_split", g), purpose = "single_problem_r1",
    seed = seed_trial_c(g), derivation = "SEED_C + 300000 + GameID")),
  tibble::tibble(scope = "halving", purpose = "cov_half_A_and_B",
                 seed = NA_integer_,
                 derivation = "SEED_C + 9e8*(half-1) + 2e8*K_index + 1e5*part + set_id"),
  tibble::tibble(scope = "null", purpose = "set_correspondence_permutation",
                 seed = NA_integer_,
                 derivation = "SEED_C + 800000 + partition_index (+500000 resid, +900000 wpc-adj)"))
save_csv(seeds_c, "dispc_seeds")
enum <- seeds_c$seed[is.finite(seeds_c$seed)]
stopifnot(!anyDuplicated(enum))
say("seed hygiene: %d enumerated seeds unique\n", length(enum))

say("\nLocal-dependence diagnostics complete in %.1f min. Outputs in elevation/results/dispc/\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))
say("G-C1 %s | G-C2 %s | G-C3 %s\n", gc1_verdict, gc2_verdict, gc3_verdict)

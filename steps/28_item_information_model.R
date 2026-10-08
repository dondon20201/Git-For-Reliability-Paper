# Step 28: joint item-information model
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e0/, e1/, e2/, e3b/ (steps 19-21, 24)
# Writes:  output/elevation/results/p2_irt/ (p2_gate.json, identifiability.csv, mechanism_regression.csv, confirmation_engine.csv, ...)
# Paper:   The reliability ceiling reflects true-score variance (A joint model of item information); Supplement
# Runtime: about 73 minutes with 12 parallel workers

# P2: exploratory joint item-response model of problem discrimination
# lambda_j and participant trait u_i. It tests the proposed pathway from
# information features to item discrimination to set reliability that the
# two-stage H2 test (E4, step 25) did not detect. The pre-specified H2 FAIL
# remains the confirmatory result, regardless of P2's outcome.
#
# All estimates use the 344 build subjects from step 20. Choice/outcome columns
# are dropped from the full-sample frame, excluding the 137 freeze and 205 test
# subjects' choices. No designs are optimized or constructed.
#
# The shared R/ functions (build_subject_data, fit_loss_baselines, estimate_half,
# load_cpc18, is_clean_row, clean_problem_features, add_scaled_features) are
# sourced without redefinition. E2 features (step 21), E3b outcomes (step 24),
# and E0 ceilings (step 18) are read from their frozen, md5-checked files.
# Before loading estimation, D1_SPOT checks exact agreement between this
# script's count matrix and the loss-informative trial construction through
# estimate_half for every build subject.
#
# Phase 1a fits the joint model and checks G0 (identifiability). If G0 fails,
# the script writes the identifiability results and stops. Phase 1b
# (G1/G2/G3 and HG1-HG4) runs only if G0 passes.
#
# P2_SMOKE=1 reduces loop counts and redirects outputs to a scratch directory.
# Smoke-run numbers are not research results.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
  library(mirt)
})

source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))
source(here("R", "step1_measures.R"))

SMOKE <- identical(Sys.getenv("P2_SMOKE"), "1")
res_dir <- if (SMOKE) {
  sd0 <- Sys.getenv("P2_SMOKE_DIR")
  if (nzchar(sd0)) sd0 else file.path(tempdir(), "p2_smoke")
} else here("output", "elevation", "results", "p2_irt")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
if (!SMOKE) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  invisible(df)
}
say <- function(...) cat(sprintf(...), sep = "")
t_start <- Sys.time()
if (SMOKE) say("*** SMOKE RUN: reduced constants, outputs to %s - NOT A RESULT ***\n",
               res_dir)

# ---- frozen constants (fixed BEFORE any result exists) ----
RAW_MD5 <- "38a26407a36c690be969f9d9d3fe3e5a"
MD5_IN <- c(
  e1_split      = "304b2b7ee67baedab0a8a24f9f093392",
  e2_dictionary = "dad7973dfacdf2426cc37840343a397c",
  e2_config     = "c338632d5437533365b6713580d55def",
  e3b_outcome   = "22b120e590cb6c0ec677c8292c25889d",
  e3b_subjsplit = "5434c077a1325f42badc10e26f649ba3",
  e3b_sim       = "ce7a2cbaf243bbcbffacc8fa52588d11",
  e3b_identity  = "f22607d0d7e83a49aa522d77aa4798a7",
  e3b_seeds     = "35a766e003cede7fe5f36c06e2bb40dc",
  e0_ceiling    = "b567ca580844e5d8abe37394fbeb4888")

SEED_P2 <- 20260724L
POS     <- 6:25                   # experience regime, frozen
R_REPS  <- length(POS)            # r = 20
NQ      <- 61L                    # Gauss-Hermite nodes (MML engine)
TARGETS <- c("loss_avoid_adj", "ev_sens_agn")
PRIMARY_TARGET <- "loss_avoid_adj"
FEATURES <- c("F1_target_predictor_spread", "F2_loss_exposure",
              "F3_predictor_collinearity", "F4_expected_information",
              "F5_choice_entropy")
DECOR_ORDER <- c("F1_target_predictor_spread", "F2_loss_exposure",
                 "F3_predictor_collinearity", "F5_choice_entropy")
F4_RESID <- "F4_expected_information_resid"

DEGEN_MAX_FRAC <- 0.10            # >= 10% degenerate items -> G0b FAIL
LAMBDA_ABS_MAX <- 6               # |lambda_j| > 6 flags degenerate
SE_MULT_MAX    <- 5               # SE_j > 5x median SE flags degenerate
NEAR_ZERO_VAR  <- 1e-4            # between-subject var of the cell proportion
COND_MAX       <- 30              # G0c design-matrix condition number
G0D_RECOV_MIN  <- 0.60            # corr(true lambda, est lambda)
G1_CV_MIN      <- 0.30
G3_MIN         <- 0.30
K_CALIB   <- c(8L, 10L, 12L)      # K = 6 reported but thin, non-gating
K_REPORT  <- c(6L, 8L, 10L, 12L)
BETA_PLANT <- 0.30                # planted feature -> loading slope (G0d)
LAM_BASE   <- 0.60
N_SIM_REP <- if (SMOKE) 6L else 200L
N_PERM    <- if (SMOKE) 6L else 200L
N_BOOT    <- if (SMOKE) 50L else 2000L
N_WORKERS <- if (SMOKE) 4L else 12L

seed_recovery <- function(rep) SEED_P2 + 1000000L + rep
seed_null_a   <- function(rep) SEED_P2 + 2000000L + rep
seed_null_b   <- function(rep) SEED_P2 + 3000000L + rep
seed_perm     <- function(rep) SEED_P2 + 4000000L + rep
seed_odnull   <- function(rep) SEED_P2 + 5000000L + rep
seed_boot     <- function(t_idx) SEED_P2 + 6000000L + 1000L * t_idx

# =============================================================================
# 1. Inputs and the build-only structural guarantee
# =============================================================================
raw_path <- here("data", "raw", "All estimation raw data.csv")
stopifnot(file.exists(raw_path))
raw_md5 <- unname(tools::md5sum(raw_path))
stopifnot(identical(raw_md5, RAW_MD5))
say("raw md5 verified: %s\n", raw_md5)

path_of <- c(
  e1_split      = here("output", "elevation", "results", "e1", "participant_split.csv"),
  e2_dictionary = here("output", "elevation", "results", "e2",
                       "problem_information_dictionary.csv"),
  e2_config     = here("output", "elevation", "results", "e2", "e2_feature_config.json"),
  e3b_outcome   = here("output", "elevation", "results", "e3b",
                       "outcome_reliability_real.csv"),
  e3b_subjsplit = here("output", "elevation", "results", "e3b",
                       "reliability_subject_split_real.csv"),
  e3b_sim       = here("output", "elevation", "results", "e3b",
                       "reliability_highprec_sim.csv"),
  e3b_identity  = here("output", "elevation", "results", "e3b", "e3b_set_identity.csv"),
  e3b_seeds     = here("output", "elevation", "results", "e3b", "e3b_seeds.csv"),
  e0_ceiling    = here("output", "elevation", "results", "e0",
                       "variance_ceiling_final.csv"))
in_md5 <- vapply(path_of, function(p) unname(tools::md5sum(p)), character(1))
for (k in names(MD5_IN))
  if (!identical(unname(in_md5[[k]]), unname(MD5_IN[[k]])))
    stop(sprintf("frozen input %s md5 mismatch: got %s expected %s",
                 k, in_md5[[k]], MD5_IN[[k]]))
say("frozen input md5s verified (%d files)\n", length(MD5_IN))

split_tbl <- readr::read_csv(path_of[["e1_split"]], show_col_types = FALSE)
stopifnot(identical(names(split_tbl),
                    c("SubjID", "Set", "split", "held_out_until_E7")))
build_ids <- split_tbl$SubjID[split_tbl$split == "build"]
held_out  <- split_tbl$SubjID[split_tbl$split != "build"]
n_build   <- length(build_ids)
stopifnot(n_build == 344L, sum(split_tbl$split == "freeze") == 137L,
          sum(split_tbl$split == "test") == 205L, !anyDuplicated(build_ids))
say("split: build %d / freeze %d (not read) / test %d (sealed)\n",
    n_build, 137L, 205L)

raw <- load_cpc18(raw_path)
stopifnot(all(cpc18_expected_columns %in% names(raw)))
choice_cols <- c("B", "Button", "Payoff", "Forgone", "RT", "Apay", "Bpay",
                 "Feedback")
design_all <- raw |> dplyr::select(-dplyr::any_of(choice_cols))
stopifnot(!any(choice_cols %in% names(design_all)))
build_raw <- raw |> dplyr::filter(SubjID %in% build_ids)
rm(raw); invisible(gc())
stopifnot(setequal(unique(build_raw$SubjID), build_ids),
          !any(build_raw$SubjID %in% held_out))                          # build only

feat_all <- clean_problem_features(design_all) |> add_scaled_features() |>
  dplyr::arrange(GameID)
stopifnot(nrow(feat_all) == 105L)
build_trials <- build_raw |>
  dplyr::filter(is_clean_row(build_raw)) |>
  dplyr::select(SubjID, Set, GameID, Trial, B) |>
  dplyr::left_join(
    feat_all |> dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev,
                              loss_informative, has_loss_b, ev_a, ev_b),
    by = "GameID")
stopifnot(!any(is.na(build_trials$z_ev)),
          dplyr::n_distinct(build_trials$SubjID) == 344L)

baselines <- fit_loss_baselines(build_trials)                            # FROZEN
subj_data <- build_subject_data(build_trials, baselines)                 # FROZEN
stopifnot(length(subj_data) == 344L)
ids <- names(subj_data)
subj_form <- build_trials |> dplyr::distinct(SubjID, Set) |>
  dplyr::mutate(o = match(as.character(SubjID), ids)) |> dplyr::arrange(o)
stopifnot(identical(as.character(subj_form$SubjID), ids))
say("frozen loss baselines (exp): b0=%.6f b1=%.6f\n",
    baselines$exp[1], baselines$exp[2])

E2 <- readr::read_csv(path_of[["e2_dictionary"]], show_col_types = FALSE) |>
  dplyr::arrange(GameID)
stopifnot(nrow(E2) == 105L, identical(E2$GameID, feat_all$GameID))
inS <- as.matrix(E2[, paste0("in_Set", 1:7)])
stopifnot(all(rowSums(inS) == 1L))                # forms are disjoint
E2$form <- max.col(inS)
say("item universe: %d clean problems, %d loss-informative; forms disjoint\n",
    nrow(E2), sum(E2$loss_informative))

e0 <- readr::read_csv(path_of[["e0_ceiling"]], show_col_types = FALSE)
r1_of <- function(R, k) R / (k - (k - 1) * R)
sb_k  <- function(k, r1) k * r1 / (1 + (k - 1) * r1)
R1 <- vapply(TARGETS, function(tg) {
  rw <- e0[e0$construct == tg, ]
  stopifnot(nrow(rw) == 1L, rw$k_ref == 14)
  r1_of(rw$ceiling_fulldata, 14)
}, numeric(1))
stopifnot(abs(round(R1["loss_avoid_adj"], 4) - 0.0787) < 1e-12,
          abs(round(R1["ev_sens_agn"], 5) - 0.02717) < 1e-12)
say("E0 anchors: r1 loss %.6f, ev %.6f\n", R1[1], R1[2])

# =============================================================================
# 2. Binomial count matrices + D1 spot check against the frozen estimate_half
# =============================================================================
build_counts <- function(games, mode) {
  Y <- matrix(NA_real_, length(ids), length(games),
              dimnames = list(ids, as.character(games)))
  for (i in seq_along(ids)) {
    sd <- subj_data[[ids[i]]]
    Bm <- sd$B_mat[, POS, drop = FALSE]
    flip <- if (mode == "noloss") sd$feat$has_loss_b else rep(FALSE, nrow(Bm))
    Bm2 <- Bm
    if (any(flip)) Bm2[flip, ] <- 1 - Bm[flip, , drop = FALSE]
    j <- match(as.character(sd$feat$GameID), colnames(Y))
    ok <- !is.na(j)
    Y[i, j[ok]] <- rowSums(Bm2)[ok]
  }
  Y
}
li_games  <- E2$GameID[E2$loss_informative]
all_games <- E2$GameID
if (SMOKE) {
  # Smoke mode also shrinks the ITEM universe so the item-factor fits take
  # seconds instead of tens of minutes. This exercises every code path; its
  # numbers are not results and it never writes to the results directory.
  li_games  <- li_games[seq_len(min(14L, length(li_games)))]
  all_games <- all_games[seq_len(min(30L, length(all_games)))]
  say("SMOKE: item universe shrunk to %d loss-informative / %d clean items\n",
      length(li_games), length(all_games))
}
Y_primary   <- build_counts(li_games, "noloss")
Y_secondary <- build_counts(all_games, "B")
say("primary counts: %d x %d, %.4f observed (cells/subject %d-%d)\n",
    nrow(Y_primary), ncol(Y_primary), mean(!is.na(Y_primary)),
    min(rowSums(!is.na(Y_primary))), max(rowSums(!is.na(Y_primary))))
say("secondary counts: %d x %d, %.4f observed (cells/subject %d-%d)\n",
    nrow(Y_secondary), ncol(Y_secondary), mean(!is.na(Y_secondary)),
    min(rowSums(!is.na(Y_secondary))), max(rowSums(!is.na(Y_secondary))))

spot_rows <- vector("list", length(ids)); bad_spot <- 0L
for (i in seq_along(ids)) {
  sd <- subj_data[[ids[i]]]
  rows <- seq_len(nrow(sd$feat))
  eh <- estimate_half(sd, rows, POS, "exp")                              # FROZEN
  li <- which(sd$feat$loss_informative)
  g <- as.character(sd$feat$GameID[li])
  keep_g <- g %in% colnames(Y_primary)          # all TRUE outside smoke mode
  ymat <- as.numeric(Y_primary[i, g[keep_g], drop = TRUE])
  nB <- rowSums(sd$B_mat[rows, POS, drop = FALSE])
  counts <- ifelse(sd$feat$has_loss_b[li], R_REPS - nB[li], nB[li])
  ok_counts <- identical(ymat, unname(as.numeric(counts[keep_g])))
  # The frozen index is defined over ALL of a subject's loss-informative
  # problems, so it is only reconstructible when the item universe is whole.
  full_universe <- all(keep_g)
  loss_from_mat <- if (length(li) == 0L || !full_universe) NA_real_ else
    (sum(ymat) - R_REPS * sum(sd$feat$p_base_exp[li])) / (R_REPS * length(li))
  ok_index <- if (!full_universe) NA else
    (is.na(eh[["loss_avoid_adj"]]) && is.na(loss_from_mat)) ||
    identical(eh[["loss_avoid_adj"]], loss_from_mat)
  if (!ok_counts || isFALSE(ok_index)) bad_spot <- bad_spot + 1L
  spot_rows[[i]] <- tibble::tibble(
    SubjID = ids[i], form = subj_form$Set[i], n_loss_items = length(li),
    counts_bit_identical = ok_counts, index_bit_identical = ok_index,
    loss_avoid_adj_frozen = eh[["loss_avoid_adj"]],
    loss_avoid_adj_from_counts = loss_from_mat)
}
d1_spotcheck <- dplyr::bind_rows(spot_rows)
say("D1_SPOT: %d subjects re-derived through frozen estimate_half; mismatches = %d\n",
    nrow(d1_spotcheck), bad_spot)
stopifnot(bad_spot == 0L)                                                # D1_SPOT

# =============================================================================
# 3. Per-item information features (frozen E2 values) + decorrelation
# =============================================================================
# The E2 dictionary is a PER-PROBLEM file; E4 consumed it at the SET level as
# sums/means over a set's problems. F4 and F5 decompose EXACTLY to the item
# (E4: F4 = sum of the per-item value; F5 = mean of the per-item value). F1, F2
# and F3 have no exact one-item decomposition of E4's set-level formula, so the
# per-item analogue below is a DEFINITION MADE IN THIS SCRIPT, labeled as such
# in the gate output (code_side_definitions):
#   F1 per item = squared deviation of the target predictor from its pool mean
#                 (the item's additive contribution to E4's set variance)
#   F2 per item = loss exposure magnitude eloss_a + eloss_b (E4's mean_eloss
#                 input). E4's set-level F2 was the COUNT of loss-informative
#                 problems, which is the constant 1 on the primary item
#                 universe (loss-informative by construction) and therefore
#                 carries no item-level variance.
#   F3 per item = max over the three predictor pairs of |zc_p * zc_q| on
#                 pool-centered z (the item's additive contribution to the
#                 cross-products that drive E4's set-level |corr|)
item_features <- function(target) {
  gid <- if (target == PRIMARY_TARGET) li_games else all_games
  D <- E2[match(gid, E2$GameID), ]
  zc <- scale(as.matrix(D[, c("z_ev", "z_sd", "z_eloss")]), TRUE, FALSE)
  f3 <- pmax(abs(zc[, 1] * zc[, 2]), abs(zc[, 1] * zc[, 3]),
             abs(zc[, 2] * zc[, 3]))
  tp <- if (target == PRIMARY_TARGET) D$d_ev_noloss else D$z_ev
  tibble::tibble(
    GameID = D$GameID, form = D$form, domain = D$domain,
    near_parity = D$near_parity, dominance = D$dominance,
    z_ev = D$z_ev, d_ev_noloss = D$d_ev_noloss,
    loss_informative = D$loss_informative,
    F1_target_predictor_spread = (tp - mean(tp))^2,
    F2_loss_exposure = D$eloss_a + D$eloss_b,
    F3_predictor_collinearity = as.numeric(f3),
    F4_expected_information = if (target == PRIMARY_TARGET)
      D$info_loss_avoid_adj else D$finfo_b_ev,
    F5_choice_entropy = D$choice_entropy_exp)
}
FEAT <- lapply(TARGETS, item_features); names(FEAT) <- TARGETS
for (tg in TARGETS) stopifnot(!anyNA(FEAT[[tg]][, FEATURES]))

decorrelate <- function(F) {
  X0 <- as.matrix(F[, DECOR_ORDER])
  r <- stats::lm.fit(cbind(1, X0), F[["F4_expected_information"]])$residuals
  out <- cbind(X0, r)
  colnames(out) <- c(DECOR_ORDER, F4_RESID)
  out
}
info_pc1 <- function(F) {
  M <- scale(as.matrix(F[, c("F2_loss_exposure", "F4_expected_information")]))
  as.numeric(M %*% eigen(stats::cov(M), symmetric = TRUE)$vectors[, 1])
}
zsc <- function(X) {
  mu <- colMeans(X); sg <- apply(X, 2, stats::sd)
  sg[!is.finite(sg) | sg == 0] <- 1
  sweep(sweep(X, 2, mu, "-"), 2, sg, "/")
}

# =============================================================================
# 4. Estimation engines
# =============================================================================
# The pre-specified model is
#   logit(E[y_ij / r]) = beta0_j + lambda_j * u_i + gamma * d_ev_noloss_ij.
# STRUCTURAL NOTE (reported, not worked around): d_ev_noloss_ij is a PROBLEM
# level constant (it does not vary across subjects within an item), so
# gamma * d_ev_noloss_j is exactly collinear with the free item easiness
# beta0_j. gamma is therefore not separately identified once beta0_j is free
# per item, and the model reduces algebraically to the one-factor (2PL) item
# model with the EV control absorbed into beta0_j. The reduction is forced by
# algebra, not chosen; it is quantified in identifiability.csv, and the
# E8-style common-intercept variant (which DOES identify gamma) is fitted
# alongside for reference. The SECONDARY model reduces the same way, with the
# identified loading lambda_j = kappa_j * z_ev_j; kappa_j = lambda_j / z_ev_j
# is reported with an instability flag where |z_ev_j| is near zero.

# PRIMARY engine: mirt with a custom Binomial(20) item. A cell is r = 20
# conditionally independent Bernoulli trials at one probability, i.e. a
# 21-category item with exactly two free parameters. This expresses the
# pre-specified binomial likelihood natively for mirt's EM. NOTE: the literal
# 20-column replicate expansion offered as an alternative is NOT usable on these
# data (many replicate columns are single-category and mirt refuses them).
bin_trace <- function(par, Theta, ncat, ...) {
  p <- stats::plogis(par[1] + par[2] * as.vector(Theta))
  vapply(0:(ncat - 1L), function(y) stats::dbinom(y, ncat - 1L, p),
         numeric(length(p)))
}
fit_mirt <- function(Y, se = TRUE) {
  dat <- as.data.frame(Y); names(dat) <- paste0("G", colnames(Y))
  it <- mirt::createItem(name = "binom20", par = c(d = 0, a = 0.3),
                         est = c(TRUE, TRUE), P = bin_trace,
                         lbound = c(-Inf, -Inf), ubound = c(Inf, Inf))
  t0 <- Sys.time()
  fit <- mirt::mirt(dat, 1, itemtype = rep("binom20", ncol(dat)),
                    customItems = list(binom20 = it), SE = se, verbose = FALSE,
                    technical = list(NCYCLES = 5000,
                                     customK = rep(R_REPS + 1L, ncol(dat))))
  cf <- mirt::coef(fit, printSE = TRUE)
  g <- function(row, nm) vapply(seq_len(ncol(Y)), function(j)
    unname(cf[[j]][row, nm]), numeric(1))
  list(a = g("par", "a"), d = g("par", "d"),
       se_a = if (se) g("SE", "a") else rep(NA_real_, ncol(Y)),
       se_d = if (se) g("SE", "d") else rep(NA_real_, ncol(Y)),
       converged = mirt::extract.mirt(fit, "converged"),
       logLik = mirt::extract.mirt(fit, "logLik"),
       n_cycles = tryCatch(mirt::extract.mirt(fit, "iterations"),
                           error = function(e) NA_integer_),
       secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       engine = "mirt::mirt custom Binomial(20) item, EM/MML")
}

# CONFIRMATION engine: an independent MML item-factor implementation.
# ENGINE NOTE (deviation, reported): the pre-specified first choice for
# confirmation is a Bayesian one-factor item model (brms/Stan). brms needs a
# C++ toolchain, which was unavailable, and rstanarm's precompiled models
# cannot express an item-specific loading on a latent trait (a product of two
# estimands), so the Stan path is unavailable in substance, not just in
# packaging. The pre-specified fallback is "a second frequentist item-factor
# engine"; this is one, written independently here (Golub-Welsch Gauss-Hermite
# quadrature, analytic gradient, BFGS, observed-information SEs) so the
# confirmation is a separate code path, not a second call into the same
# library. It is also the engine used for every simulation track, where a mirt
# fit per replicate is not tractable; the engine-agreement check on the real
# fit licenses that substitution.
gh_nodes <- function(nq) {
  i <- 1:(nq - 1L); J <- matrix(0, nq, nq)
  J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = rev(e$values), w = rev(e$vectors[1, ]^2))
}
GH <- gh_nodes(NQ)

fit_mml <- function(Y, se = FALSE, start_a = 0.3, gh = GH, r_reps = 20L) {
  n <- nrow(Y); J <- ncol(Y)
  obs <- !is.na(Y)
  Yo <- Y; Yo[!obs] <- 0
  No <- matrix(0, n, J); No[obs] <- r_reps
  nq <- length(gh$x)
  lnodes <- function(par) {
    d <- par[1:J]; a <- par[J + (1:J)]
    L <- matrix(0, n, nq)
    for (k in seq_len(nq)) {
      eta <- d + a * gh$x[k]
      L[, k] <- Yo %*% (-log1p(exp(-eta))) + (No - Yo) %*% (-log1p(exp(eta)))
    }
    L
  }
  negll <- function(par) {
    L <- lnodes(par); m <- apply(L, 1, max)
    -sum(m + log(rowSums(sweep(exp(L - m), 2, gh$w, "*"))))
  }
  grad <- function(par) {
    d <- par[1:J]; a <- par[J + (1:J)]
    L <- lnodes(par); m <- apply(L, 1, max)
    W <- sweep(exp(L - m), 2, gh$w, "*"); W <- W / rowSums(W)
    gd <- numeric(J); ga <- numeric(J)
    for (k in seq_len(nq)) {
      p <- stats::plogis(d + a * gh$x[k])
      cs <- colSums(W[, k] * (Yo - sweep(No, 2, p, "*")))
      gd <- gd + cs; ga <- ga + gh$x[k] * cs
    }
    -c(gd, ga)
  }
  pbar <- pmin(pmax(colSums(Yo) / pmax(colSums(No), 1), 0.02), 0.98)
  op <- stats::optim(c(stats::qlogis(pbar), rep(start_a, J)), negll, grad,
                     method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  out <- list(a = op$par[J + (1:J)], d = op$par[1:J], par = op$par,
              converged = op$convergence == 0L, logLik = -op$value,
              negll = negll, lnodes = lnodes, J = J, n = n,
              engine = "independent MML item-factor (Gauss-Hermite + BFGS)")
  if (se) {
    H <- numDeriv::hessian(negll, op$par)
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    V <- tryCatch(solve(H), error = function(e) NULL)
    out$hessian_min_eigen <- min(ev)
    out$hessian_cond <- max(ev) / min(ev)
    out$se_a <- if (is.null(V)) rep(NA_real_, J) else sqrt(pmax(diag(V)[J + (1:J)], 0))
    out$se_d <- if (is.null(V)) rep(NA_real_, J) else sqrt(pmax(diag(V)[1:J], 0))
  }
  out
}

# EAP trait scores from a fitted (d, a) pair.
eap_scores <- function(Y, d, a, gh = GH, r_reps = 20L) {
  obs <- !is.na(Y); Yo <- Y; Yo[!obs] <- 0
  No <- matrix(0, nrow(Y), ncol(Y)); No[obs] <- r_reps
  nq <- length(gh$x)
  L <- matrix(0, nrow(Y), nq)
  for (k in seq_len(nq)) {
    eta <- d + a * gh$x[k]
    L[, k] <- Yo %*% (-log1p(exp(-eta))) + (No - Yo) %*% (-log1p(exp(eta)))
  }
  m <- apply(L, 1, max)
  W <- sweep(exp(L - m), 2, gh$w, "*"); W <- W / rowSums(W)
  mu <- as.numeric(W %*% gh$x)
  sd <- sqrt(pmax(as.numeric(W %*% (gh$x^2)) - mu^2, 0))
  list(eap = mu, sd = sd)
}

# =============================================================================
# 5. Mechanism regression + calibration (shared by real and simulated tracks)
# =============================================================================
baseline_matrix <- function(F) {
  # Item-level baseline. E4's set-level baseline was problem count + broad
  # category; problem count has no item-level analogue (it is a set quantity),
  # so the item baseline is the broad category alone (domain plus the two
  # frozen E2 category flags). The substitution is reported.
  # Single-level factors are dropped rather than erroring: on the PRIMARY item
  # universe `domain` is constant ("mixed" for every loss-informative problem),
  # so the primary baseline rests on near_parity and dominance. If no category
  # varies, the baseline degenerates to the training-fold mean (handled by the
  # caller, which sees a zero-column matrix).
  d <- as.data.frame(F[, c("domain", "near_parity", "dominance")])
  d[] <- lapply(d, function(x) factor(as.character(x), exclude = NULL))
  keep <- names(d)[vapply(d, nlevels, integer(1)) >= 2L]
  if (length(keep) == 0L) return(matrix(numeric(0), nrow = nrow(F), ncol = 0L))
  stats::model.matrix(stats::reformulate(keep), data = d)[, -1, drop = FALSE]
}
mech_regression <- function(lambda, F, label) {
  keep <- is.finite(lambda)
  Xr <- decorrelate(F)[keep, , drop = FALSE]
  y <- lambda[keep]; form <- F$form[keep]
  B <- baseline_matrix(F)[keep, , drop = FALSE]
  Xs <- zsc(Xr)
  cond <- kappa(cbind(1, Xs), exact = TRUE)
  fit <- stats::lm.fit(cbind(1, Xs), y)
  bet <- fit$coefficients; bet[is.na(bet)] <- 0
  dfres <- max(length(y) - sum(!is.na(fit$coefficients)), 1L)
  s2 <- sum(fit$residuals^2) / dfres
  XtXi <- tryCatch(solve(crossprod(cbind(1, Xs))), error = function(e) NULL)
  sebet <- if (is.null(XtXi)) rep(NA_real_, length(bet)) else
    sqrt(pmax(diag(XtXi) * s2, 0))
  pc <- info_pc1(F)[keep]
  b_pc <- stats::lm.fit(cbind(1, scale(pc)), y)$coefficients[2]
  marg <- vapply(FEATURES, function(f)
    suppressWarnings(stats::cor(F[[f]][keep], y, method = "spearman")),
    numeric(1))
  cvp <- rep(NA_real_, length(y)); cvb <- rep(NA_real_, length(y))
  for (g in unique(form)) {
    tr <- form != g; te <- !tr
    if (sum(tr) < ncol(Xr) + 2L || sum(te) == 0L) next
    mu <- colMeans(Xr[tr, , drop = FALSE])
    sg <- apply(Xr[tr, , drop = FALSE], 2, stats::sd)
    sg[!is.finite(sg) | sg == 0] <- 1
    Ztr <- sweep(sweep(Xr[tr, , drop = FALSE], 2, mu, "-"), 2, sg, "/")
    Zte <- sweep(sweep(Xr[te, , drop = FALSE], 2, mu, "-"), 2, sg, "/")
    bb <- stats::lm.fit(cbind(1, Ztr), y[tr])$coefficients; bb[is.na(bb)] <- 0
    cvp[te] <- cbind(1, Zte) %*% bb
    kp <- apply(B[tr, , drop = FALSE], 2, stats::sd) > 0
    if (any(kp)) {
      cb <- stats::lm.fit(cbind(1, B[tr, kp, drop = FALSE]),
                          y[tr])$coefficients
      cb[is.na(cb)] <- 0
      cvb[te] <- cbind(1, B[te, kp, drop = FALSE]) %*% cb
    } else cvb[te] <- mean(y[tr])
  }
  sp <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5L || stats::sd(a[ok]) == 0 || stats::sd(b[ok]) == 0)
      return(NA_real_)
    suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
  }
  cvs <- sp(cvp, y); cvbs <- sp(cvb, y)
  list(label = label, n_items = length(y), cond = cond,
       beta = setNames(bet[-1], colnames(Xs)), se_beta = sebet[-1],
       beta_pc1 = unname(b_pc), marginal_spearman = marg,
       cv_spearman = cvs, cv_spearman_baseline = cvbs, cv_delta = cvs - cvbs,
       cv_pred = cvp, y = y, keep = keep)
}

# IRT-implied set reliability: one-factor test information at the population
# mean, I(S) = sum_j lambda_j^2 ; rho_hat(S) = I(S) / (I(S) + 1).
rho_hat_set <- function(lam, lam_names, game_ids_chr) {
  vapply(strsplit(game_ids_chr, ";", fixed = TRUE), function(v) {
    j <- match(v, lam_names); j <- j[!is.na(j)]
    I <- if (length(j) == 0L) 0 else sum(lam[j]^2, na.rm = TRUE)
    I / (I + 1)
  }, numeric(1))
}

# =============================================================================
# 6. PHASE 1a: fit the joint model on real BUILD data and evaluate G0
# =============================================================================
say("\n=== PHASE 1a: joint item-factor fit (PRIMARY = %s) ===\n", PRIMARY_TARGET)
FITS <- list(); MMLF <- list()
for (tg in TARGETS) {
  Y <- if (tg == PRIMARY_TARGET) Y_primary else Y_secondary
  say("fitting mirt [%s]: %d items ...\n", tg, ncol(Y))
  FITS[[tg]] <- fit_mirt(Y, se = TRUE)
  say("  mirt: converged=%s logLik=%.4f cycles=%s (%.1f min)\n",
      FITS[[tg]]$converged, FITS[[tg]]$logLik,
      paste(FITS[[tg]]$n_cycles, collapse = ","), FITS[[tg]]$secs / 60)
  say("fitting confirmation MML [%s] ...\n", tg)
  MMLF[[tg]] <- fit_mml(Y, se = TRUE)
  say("  MML: converged=%s logLik=%.4f minEigenH=%.6g\n",
      MMLF[[tg]]$converged, MMLF[[tg]]$logLik, MMLF[[tg]]$hessian_min_eigen)
}

# The EV control enters the primary model as d_ev_noloss and the secondary
# model as z_ev; both are problem-level constants, so each is tested against
# that target's own item-dummy design.
alias_rank <- function(gid, control) {
  D <- E2[match(gid, E2$GameID), ]
  ctrl <- D[[control]]
  stopifnot(all(is.finite(ctrl)))
  Dm <- diag(length(gid))
  c(qr(Dm)$rank, qr(cbind(Dm, ctrl))$rank)
}
al_p <- alias_rank(li_games, "d_ev_noloss")
al_s <- alias_rank(all_games, "z_ev")
say("gamma aliasing (PRIMARY): rank(items)=%d rank(items + d_ev_noloss)=%d -> gamma %s\n",
    al_p[1], al_p[2], if (al_p[1] == al_p[2]) "NOT identified" else "identified")

gamma_ref <- {
  fr <- do.call(rbind, lapply(seq_along(ids), function(i) {
    sd <- subj_data[[ids[i]]]
    li <- which(sd$feat$loss_informative)
    if (length(li) == 0L) return(NULL)
    nB <- rowSums(sd$B_mat[, POS, drop = FALSE])
    data.frame(y = ifelse(sd$feat$has_loss_b[li], R_REPS - nB[li], nB[li]) / R_REPS,
               d_ev_noloss = sd$feat$d_ev_noloss[li], w = R_REPS)
  }))
  unname(stats::coef(stats::glm(y ~ 1 + d_ev_noloss, family = stats::binomial(),
                                weights = fr$w, data = fr)))
}
say("E8-style common-intercept reference: intercept %.6f gamma %.6f\n",
    gamma_ref[1], gamma_ref[2])

lam <- FITS[[PRIMARY_TARGET]]$a; se_lam <- FITS[[PRIMARY_TARGET]]$se_a
names(lam) <- names(se_lam) <- as.character(li_games)
prop_var <- apply(Y_primary, 2, function(y) stats::var(y[!is.na(y)] / R_REPS))
degen_var <- prop_var < NEAR_ZERO_VAR
degen_se  <- is.finite(se_lam) & se_lam > SE_MULT_MAX * stats::median(se_lam, na.rm = TRUE)
degen_mag <- abs(lam) > LAMBDA_ABS_MAX
degen_nf  <- !is.finite(lam) | !is.finite(se_lam)
degenerate <- degen_var | degen_se | degen_mag | degen_nf
degen_frac <- mean(degenerate)

g0a <- isTRUE(FITS[[PRIMARY_TARGET]]$converged) && all(is.finite(lam)) &&
  all(is.finite(se_lam))
g0b <- all(is.finite(lam)) && all(is.finite(se_lam)) && degen_frac < DEGEN_MAX_FRAC
mech_real <- mech_regression(lam, FEAT[[PRIMARY_TARGET]], "real_primary")
g0c <- is.finite(mech_real$cond) && mech_real$cond < COND_MAX &&
  all(is.finite(mech_real$se_beta))

say("\nG0a convergence: converged=%s finite loadings=%s finite SEs=%s -> %s\n",
    FITS[[PRIMARY_TARGET]]$converged, all(is.finite(lam)),
    all(is.finite(se_lam)), if (g0a) "PASS" else "FAIL")
say("G0b degeneracy: %d/%d flagged (%.4f) [var %d, SE %d, |lambda|>%g %d] -> %s\n",
    sum(degenerate), length(lam), degen_frac, sum(degen_var), sum(degen_se),
    LAMBDA_ABS_MAX, sum(degen_mag), if (g0b) "PASS" else "FAIL")
say("  loadings min %.4f median %.4f max %.4f ; SEs min %.4f median %.4f max %.4f\n",
    min(lam), stats::median(lam), max(lam), min(se_lam),
    stats::median(se_lam), max(se_lam))
say("G0c conditioning: condition number %.4f (< %g) ; finite feature-beta SEs=%s -> %s\n",
    mech_real$cond, COND_MAX, all(is.finite(mech_real$se_beta)),
    if (g0c) "PASS" else "FAIL")

# ---- G0d: simulation recovery + null-A + null-B ---------------------------
OBS_P <- !is.na(Y_primary)
F_pri <- FEAT[[PRIMARY_TARGET]]
X_pri <- zsc(decorrelate(F_pri))
OD_RATIO <- vapply(seq_len(ncol(Y_primary)), function(j) {
  y <- Y_primary[!is.na(Y_primary[, j]), j] / R_REPS
  pb <- mean(y); stats::var(y) / (pb * (1 - pb) / R_REPS)
}, numeric(1))

sim_one <- function(mode, seed, obs, Xz, planted_col, beta_plant, lam_base,
                    r_reps) {
  set.seed(seed)
  n <- nrow(obs); J <- ncol(obs)
  lam_true <- switch(
    mode,
    recovery = lam_base + beta_plant * Xz[, planted_col] + stats::rnorm(J, 0, 0.10),
    null_a   = lam_base + stats::rnorm(J, 0, 0.30),
    null_b   = lam_base + beta_plant * Xz[, planted_col] + stats::rnorm(J, 0, 0.10))
  d_true <- stats::rnorm(J, 0, 0.5)
  u <- if (mode == "null_b") rep(0, n) else stats::rnorm(n, 0, 1)
  Y <- matrix(NA_real_, n, J, dimnames = dimnames(obs))
  for (j in seq_len(J)) {
    ii <- which(obs[, j])
    Y[ii, j] <- stats::rbinom(length(ii), r_reps,
                              stats::plogis(d_true[j] + lam_true[j] * u[ii]))
  }
  list(Y = Y, lam_true = lam_true)
}

say("\nG0d: simulation recovery + null-A + null-B (%d replicates each)\n", N_SIM_REP)
cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(stats); library(numDeriv)}); NULL }))
parallel::clusterExport(cl, c(
  "fit_mml", "sim_one", "R_REPS", "GH", "NQ", "decorrelate", "zsc",
  "DECOR_ORDER", "FEATURES", "F4_RESID", "info_pc1", "baseline_matrix",
  "mech_regression", "OBS_P", "X_pri", "F_pri", "BETA_PLANT", "LAM_BASE",
  "OD_RATIO", "SEED_P2", "Y_primary", "seed_odnull", "seed_perm", "lam"),
  envir = environment())

run_sim_block <- function(mode, seed_fun) {
  # force() both arguments: an unforced promise serialized to a worker would
  # try to resolve the caller's symbol there, where it does not exist.
  force(mode); force(seed_fun)
  parallel::parLapply(cl, seq_len(N_SIM_REP), function(rep) {
    s <- sim_one(mode, seed_fun(rep), OBS_P, X_pri, F4_RESID, BETA_PLANT,
                 LAM_BASE, R_REPS)
    f <- fit_mml(s$Y, se = FALSE)
    mr <- mech_regression(f$a, F_pri, paste0("sim_", mode))
    tibble::tibble(
      rep = rep, mode = mode, seed = seed_fun(rep),
      beta_planted_feature = unname(mr$beta[F4_RESID]),
      cv_spearman = mr$cv_spearman, cv_spearman_baseline = mr$cv_spearman_baseline,
      cv_delta = mr$cv_delta,
      corr_lambda_spearman = suppressWarnings(
        stats::cor(s$lam_true, f$a, method = "spearman")),
      corr_lambda_pearson = suppressWarnings(stats::cor(s$lam_true, f$a)),
      lam_est_median = stats::median(f$a), lam_est_absmax = max(abs(f$a)),
      converged = f$converged)
  })
}
sim_tbl <- dplyr::bind_rows(
  run_sim_block("recovery", seed_recovery),
  run_sim_block("null_a",   seed_null_a),
  run_sim_block("null_b",   seed_null_b))

q90 <- function(x) unname(stats::quantile(x, c(0.05, 0.95), na.rm = TRUE))
gv <- function(mode, col) sim_tbl[[col]][sim_tbl$mode == mode]
rec_beta_ci <- q90(gv("recovery", "beta_planted_feature"))
rec_corr <- stats::median(gv("recovery", "corr_lambda_spearman"), na.rm = TRUE)
nulla_ci <- q90(gv("null_a", "beta_planted_feature"))
nullb_ci <- q90(gv("null_b", "beta_planted_feature"))
nullb_lam <- stats::median(gv("null_b", "lam_est_median"), na.rm = TRUE)
g0d_recovery <- rec_beta_ci[1] <= BETA_PLANT && BETA_PLANT <= rec_beta_ci[2] &&
  rec_corr > G0D_RECOV_MIN
g0d_null_a <- nulla_ci[1] <= 0 && 0 <= nulla_ci[2]
g0d_null_b <- nullb_ci[1] <= 0 && 0 <= nullb_ci[2]
g0d <- g0d_recovery && g0d_null_a && g0d_null_b
say("G0d recovery: beta 90%% [%.4f, %.4f] vs planted %.2f ; median Spearman(true,est) %.4f -> %s\n",
    rec_beta_ci[1], rec_beta_ci[2], BETA_PLANT, rec_corr,
    if (g0d_recovery) "PASS" else "FAIL")
say("G0d null-A: beta 90%% [%.4f, %.4f] -> %s\n", nulla_ci[1], nulla_ci[2],
    if (g0d_null_a) "PASS" else "FAIL")
say("G0d null-B: beta 90%% [%.4f, %.4f] ; median lambda %.4f -> %s\n",
    nullb_ci[1], nullb_ci[2], nullb_lam, if (g0d_null_b) "PASS" else "FAIL")

# ---- Added overdispersion diagnostic (not a pre-declared gate) -----------
# The pre-declared Null-B simulates under the model's OWN Binomial(20)
# likelihood. The real cells are strongly overdispersed relative to
# Binomial(20), so an additional diagnostic destroys any common trait while
# preserving each item's observed marginal count distribution EXACTLY: for each
# item independently, the observed counts are permuted across the subjects who
# saw that item. Any loading the model then reports is manufactured from
# marginal overdispersion alone.
say("\nCode-added overdispersion null (item-wise subject permutation, %d reps)\n",
    N_PERM)
say("  observed overdispersion (obs var / Binomial(20) var of cell proportion): min %.2f median %.2f max %.2f\n",
    min(OD_RATIO), stats::median(OD_RATIO), max(OD_RATIO))
od_null <- dplyr::bind_rows(parallel::parLapply(cl, seq_len(N_PERM), function(rep) {
  set.seed(seed_odnull(rep))
  Y <- Y_primary
  for (j in seq_len(ncol(Y))) {
    ii <- which(!is.na(Y[, j]))
    Y[ii, j] <- Y[sample(ii), j]
  }
  f <- fit_mml(Y, se = FALSE)
  # Two K = 12 summaries. The "typical" one uses the median squared loading
  # (a fair stand-in for an arbitrary 12-item set); the "top" one takes the 12
  # largest and therefore capitalizes on chance, so it is reported only as the
  # best-case selection.
  I_typ <- 12 * stats::median(f$a^2)
  I_top <- sum(sort(f$a^2, decreasing = TRUE)[seq_len(min(12L, length(f$a)))])
  tibble::tibble(rep = rep, mode = "od_null_no_common_trait",
                 seed = seed_odnull(rep),
                 lam_median = stats::median(f$a), lam_absmax = max(abs(f$a)),
                 rho_hat_typical12 = I_typ / (I_typ + 1),
                 rho_hat_top12 = I_top / (I_top + 1), converged = f$converged)
}))
rho_typ12 <- function(a) { I <- 12 * stats::median(a^2); I / (I + 1) }
say("  permuted-null loadings: median %.4f vs real-fit median %.4f\n",
    stats::median(od_null$lam_median), stats::median(lam))
say("  implied rho_hat at K=12: permuted null typical %.4f / top %.4f ; real fit typical %.4f (E0 ceiling %.4f)\n",
    stats::median(od_null$rho_hat_typical12), stats::median(od_null$rho_hat_top12),
    rho_typ12(lam), sb_k(12, R1[[PRIMARY_TARGET]]))

# ---- G0 verdict -----------------------------------------------------------
g0 <- g0a && g0b && g0c && g0d
say("\n=== G0 VERDICT: %s (G0a %s / G0b %s / G0c %s / G0d %s) ===\n",
    if (g0) "PASS" else "FAIL", if (g0a) "PASS" else "FAIL",
    if (g0b) "PASS" else "FAIL", if (g0c) "PASS" else "FAIL",
    if (g0d) "PASS" else "FAIL")

# =============================================================================
# 7. Phase 1a outputs (written whatever the verdict)
# =============================================================================
mk_loadings <- function(tg) {
  gid <- if (tg == PRIMARY_TARGET) li_games else all_games
  Y <- if (tg == PRIMARY_TARGET) Y_primary else Y_secondary
  f <- FITS[[tg]]; m <- MMLF[[tg]]
  F <- FEAT[[tg]]
  pv <- apply(Y, 2, function(y) stats::var(y[!is.na(y)] / R_REPS))
  dse <- is.finite(f$se_a) & f$se_a > SE_MULT_MAX * stats::median(f$se_a, na.rm = TRUE)
  tb <- tibble::tibble(
    target = tg, GameID = gid, form = F$form,
    n_subjects_observed = colSums(!is.na(Y)),
    lambda = f$a, lambda_se = f$se_a, easiness_beta0 = f$d,
    easiness_se = f$se_d,
    lambda_confirmation_mml = m$a, lambda_se_confirmation_mml = m$se_a,
    between_subject_var_proportion = unname(pv),
    overdispersion_ratio = unname(pv / {
      pb <- colMeans(Y, na.rm = TRUE) / R_REPS; pb * (1 - pb) / R_REPS }),
    degenerate_near_zero_var = unname(pv < NEAR_ZERO_VAR),
    degenerate_large_se = unname(dse),
    degenerate_large_lambda = abs(f$a) > LAMBDA_ABS_MAX,
    degenerate_nonfinite = !is.finite(f$a) | !is.finite(f$se_a))
  tb$degenerate_any <- tb$degenerate_near_zero_var | tb$degenerate_large_se |
    tb$degenerate_large_lambda | tb$degenerate_nonfinite
  if (tg != PRIMARY_TARGET) {
    tb$kappa_j <- f$a / F$z_ev
    tb$kappa_unstable_small_z_ev <- abs(F$z_ev) < 0.1
  }
  dplyr::bind_cols(tb, F[, c(FEATURES, "domain", "near_parity", "dominance")])
}
loadings_primary   <- mk_loadings(PRIMARY_TARGET)
loadings_secondary <- mk_loadings("ev_sens_agn")
save_csv(loadings_primary, "loadings_primary")
save_csv(loadings_secondary, "loadings_secondary")
save_csv(d1_spotcheck, "d1_spotcheck")

identifiability <- dplyr::bind_rows(
  tibble::tibble(check = "G0a_convergence_primary", quantity = "mirt converged",
                 value = as.numeric(isTRUE(FITS[[PRIMARY_TARGET]]$converged)),
                 threshold = 1, passes = isTRUE(FITS[[PRIMARY_TARGET]]$converged),
                 detail = FITS[[PRIMARY_TARGET]]$engine),
  tibble::tibble(check = "G0a_finite_loadings", quantity = "n nonfinite lambda",
                 value = sum(!is.finite(lam)), threshold = 0,
                 passes = all(is.finite(lam)), detail = "primary target"),
  tibble::tibble(check = "G0a_finite_ses", quantity = "n nonfinite SE",
                 value = sum(!is.finite(se_lam)), threshold = 0,
                 passes = all(is.finite(se_lam)), detail = "observed information"),
  tibble::tibble(check = "G0b_degenerate_fraction", quantity = "fraction flagged",
                 value = degen_frac, threshold = DEGEN_MAX_FRAC,
                 passes = degen_frac < DEGEN_MAX_FRAC,
                 detail = sprintf("var %d, SE %d, |lambda|>%g %d, nonfinite %d of %d items",
                                  sum(degen_var), sum(degen_se), LAMBDA_ABS_MAX,
                                  sum(degen_mag), sum(degen_nf), length(lam))),
  tibble::tibble(check = "G0c_condition_number", quantity = "kappa(design)",
                 value = mech_real$cond, threshold = COND_MAX,
                 passes = mech_real$cond < COND_MAX,
                 detail = "decorrelated feature design, standardized, with intercept"),
  tibble::tibble(check = "G0c_finite_feature_betas", quantity = "n nonfinite beta SE",
                 value = sum(!is.finite(mech_real$se_beta)), threshold = 0,
                 passes = all(is.finite(mech_real$se_beta)), detail = "OLS on decorrelated features"),
  tibble::tibble(check = "G0d_recovery_beta_ci_lo", quantity = "beta 5th pct",
                 value = rec_beta_ci[1], threshold = BETA_PLANT,
                 passes = rec_beta_ci[1] <= BETA_PLANT,
                 detail = sprintf("planted %.2f, %d replicates", BETA_PLANT, N_SIM_REP)),
  tibble::tibble(check = "G0d_recovery_beta_ci_hi", quantity = "beta 95th pct",
                 value = rec_beta_ci[2], threshold = BETA_PLANT,
                 passes = rec_beta_ci[2] >= BETA_PLANT, detail = "90% interval covers planted"),
  tibble::tibble(check = "G0d_recovery_corr", quantity = "median Spearman(true,est lambda)",
                 value = rec_corr, threshold = G0D_RECOV_MIN,
                 passes = rec_corr > G0D_RECOV_MIN, detail = "recovery track"),
  tibble::tibble(check = "G0d_null_a_ci", quantity = "beta 90% interval includes 0",
                 value = as.numeric(g0d_null_a), threshold = 1, passes = g0d_null_a,
                 detail = sprintf("[%.4f, %.4f]", nulla_ci[1], nulla_ci[2])),
  tibble::tibble(check = "G0d_null_b_ci", quantity = "beta 90% interval includes 0",
                 value = as.numeric(g0d_null_b), threshold = 1, passes = g0d_null_b,
                 detail = sprintf("[%.4f, %.4f]; median lambda %.4f", nullb_ci[1],
                                  nullb_ci[2], nullb_lam)),
  tibble::tibble(check = "STRUCTURAL_gamma_aliasing",
                 quantity = "rank(items + d_ev_noloss) - rank(items)",
                 value = al_p[2] - al_p[1], threshold = 1,
                 passes = (al_p[2] - al_p[1]) == 1,
                 detail = sprintf("gamma not separately identified with free item easiness; E8-style common-intercept reference gamma = %.6f",
                                  gamma_ref[2])),
  tibble::tibble(check = "CODE_ADDED_overdispersion_null",
                 quantity = "median lambda under item-wise subject permutation",
                 value = stats::median(od_null$lam_median), threshold = 0,
                 passes = NA,
                 detail = sprintf("real-fit median lambda %.4f; implied rho_hat at K=12: permuted-null typical %.4f / top %.4f, real fit typical %.4f, E0 ceiling %.4f; diagnostic only, not a pre-declared gate",
                                  stats::median(lam),
                                  stats::median(od_null$rho_hat_typical12),
                                  stats::median(od_null$rho_hat_top12),
                                  rho_typ12(lam), sb_k(12, R1[[PRIMARY_TARGET]]))))
save_csv(identifiability, "identifiability")
save_csv(sim_tbl, "g0d_simulation_replicates")
save_csv(od_null, "overdispersion_null_diagnostic")

eap_p <- eap_scores(Y_primary, FITS[[PRIMARY_TARGET]]$d, FITS[[PRIMARY_TARGET]]$a)
per_subject_trait_sample <- tibble::tibble(
  target = PRIMARY_TARGET, SubjID = ids, form = subj_form$Set,
  n_items_observed = rowSums(!is.na(Y_primary)),
  trait_eap = eap_p$eap, trait_eap_sd = eap_p$sd)
save_csv(per_subject_trait_sample, "per_subject_trait_sample")

per_item_export <- dplyr::bind_rows(
  loadings_primary |> dplyr::select(target, GameID, form, lambda, lambda_se,
    easiness_beta0, dplyr::all_of(FEATURES), domain, near_parity, dominance,
    n_subjects_observed, degenerate_any),
  loadings_secondary |> dplyr::select(target, GameID, form, lambda, lambda_se,
    easiness_beta0, dplyr::all_of(FEATURES), domain, near_parity, dominance,
    n_subjects_observed, degenerate_any))
save_csv(per_item_export, "per_item_export")

# Likelihood-equivalence check (the binomial-count representation must be
# shown equivalent to the binomial likelihood). The
# primary engine's custom item evaluates the full dbinom density; the
# confirmation engine drops the binomial coefficient, which is a constant in
# the parameters. The two log-likelihoods therefore differ by sum over observed
# cells of log choose(r, y), up to the difference between the two quadrature
# rules (mirt integrates on an equally spaced grid, the confirmation engine on
# Golub-Welsch Gauss-Hermite nodes) and the two convergence criteria. The
# residual below is reported, not asserted to be zero; the substantive check is
# the loading agreement between the engines.
binom_coef_const <- function(Y) {
  y <- Y[!is.na(Y)]
  sum(lchoose(R_REPS, y))
}
ll_gap <- vapply(TARGETS, function(tg) {
  Y <- if (tg == PRIMARY_TARGET) Y_primary else Y_secondary
  (FITS[[tg]]$logLik - MMLF[[tg]]$logLik) - binom_coef_const(Y)
}, numeric(1))
say("engine likelihood residual (mirt - MML - log binomial coefficient): %s\n",
    paste(sprintf("%s %.6f (%.5f%% of |logLik|)", TARGETS, ll_gap,
                  100 * abs(ll_gap) / abs(vapply(TARGETS, function(t)
                    MMLF[[t]]$logLik, numeric(1)))), collapse = " ; "))

conf_rows <- lapply(TARGETS, function(tg) {
  f <- FITS[[tg]]; m <- MMLF[[tg]]
  Y <- if (tg == PRIMARY_TARGET) Y_primary else Y_secondary
  ok <- is.finite(f$a) & is.finite(m$a)
  mr_a <- mech_regression(f$a, FEAT[[tg]], "primary_engine")
  mr_b <- mech_regression(m$a, FEAT[[tg]], "confirmation_engine")
  tibble::tibble(
    target = tg,
    primary_engine = f$engine, confirmation_engine = m$engine,
    bayesian_engine_substituted = TRUE,
    substitution_reason = "brms/Stan needs a C++ toolchain absent on this machine (E8 finding); rstanarm's precompiled models cannot express an item-specific loading on a latent trait, so the spec's named fallback (a second frequentist item-factor engine) is used",
    n_items = sum(ok),
    loading_agreement_spearman = suppressWarnings(
      stats::cor(f$a[ok], m$a[ok], method = "spearman")),
    loading_agreement_pearson = suppressWarnings(stats::cor(f$a[ok], m$a[ok])),
    max_abs_loading_diff = max(abs(f$a[ok] - m$a[ok])),
    logLik_primary = f$logLik, logLik_confirmation = m$logLik,
    log_binomial_coefficient_constant = binom_coef_const(Y),
    logLik_equivalence_residual = unname(ll_gap[[tg]]),
    beta_sign_agreement_frac = mean(sign(mr_a$beta) == sign(mr_b$beta)),
    beta_max_abs_diff = max(abs(mr_a$beta - mr_b$beta)),
    cv_spearman_primary = mr_a$cv_spearman,
    cv_spearman_confirmation = mr_b$cv_spearman)
})
confirmation_engine <- dplyr::bind_rows(conf_rows)
save_csv(confirmation_engine, "confirmation_engine")

seeds_tbl <- dplyr::bind_rows(
  tibble::tibble(part = "base", purpose = "SEED_P2", seed = SEED_P2,
                 derivation = "fixed in script before any result"),
  tibble::tibble(part = "G0d", purpose = "recovery", seed = seed_recovery(1:N_SIM_REP),
                 derivation = "SEED_P2 + 1000000 + rep"),
  tibble::tibble(part = "G0d", purpose = "null_a", seed = seed_null_a(1:N_SIM_REP),
                 derivation = "SEED_P2 + 2000000 + rep"),
  tibble::tibble(part = "G0d", purpose = "null_b", seed = seed_null_b(1:N_SIM_REP),
                 derivation = "SEED_P2 + 3000000 + rep"),
  tibble::tibble(part = "HG2", purpose = "feature_permutation", seed = seed_perm(1:N_PERM),
                 derivation = "SEED_P2 + 4000000 + rep"),
  tibble::tibble(part = "diagnostic", purpose = "overdispersion_null",
                 seed = seed_odnull(1:N_PERM),
                 derivation = "SEED_P2 + 5000000 + rep"))
save_csv(seeds_tbl, "p2_seeds")

write_gate <- function(phase, extra = list()) {
  gate <- c(list(
    spec = "elevation/P2_spec.md", smoke_run = SMOKE,
    phase_reached = phase,
    exploratory = TRUE,
    frozen_h2_fail_stands = TRUE,
    build_only = TRUE, n_build = n_build, freeze_read = FALSE, test_sealed = TRUE,
    primary_target = PRIMARY_TARGET,
    G0 = list(
      verdict = if (g0) "PASS" else "FAIL",
      G0a_convergence = list(pass = g0a,
        mirt_converged = isTRUE(FITS[[PRIMARY_TARGET]]$converged),
        n_nonfinite_lambda = sum(!is.finite(lam)),
        n_nonfinite_se = sum(!is.finite(se_lam))),
      G0b_degeneracy = list(pass = g0b, degenerate_fraction = degen_frac,
        threshold = DEGEN_MAX_FRAC, n_items = length(lam),
        n_near_zero_var = sum(degen_var), n_large_se = sum(degen_se),
        n_large_lambda = sum(degen_mag),
        lambda_min = min(lam), lambda_median = stats::median(lam),
        lambda_max = max(lam)),
      G0c_conditioning = list(pass = g0c, condition_number = mech_real$cond,
        threshold = COND_MAX,
        feature_beta_ses_finite = all(is.finite(mech_real$se_beta))),
      G0d_simulation = list(pass = g0d,
        recovery_pass = g0d_recovery, recovery_beta_ci90 = rec_beta_ci,
        planted_beta = BETA_PLANT, recovery_corr_spearman = rec_corr,
        recovery_corr_threshold = G0D_RECOV_MIN,
        null_a_pass = g0d_null_a, null_a_beta_ci90 = nulla_ci,
        null_b_pass = g0d_null_b, null_b_beta_ci90 = nullb_ci,
        null_b_median_lambda = nullb_lam, n_replicates = N_SIM_REP)),
    structural_findings = list(
      gamma_not_identified = (al_p[1] == al_p[2]),
      gamma_note = "d_ev_noloss is a problem-level constant, so the spec's gamma is exactly aliased with the free item easiness beta0_j; the model reduces to the one-factor item model with the EV control absorbed. E8-style common-intercept reference gamma reported separately.",
      gamma_reference_common_intercept = gamma_ref[2],
      forms_disjoint = TRUE,
      n_forms = 7L,
      loss_informative_items_per_form = as.integer(table(FEAT[[PRIMARY_TARGET]]$form)),
      overdispersion_median = stats::median(OD_RATIO),
      overdispersion_note = "cell counts are overdispersed relative to Binomial(20); the one-factor binomial model has no uniqueness term, so subject-by-item variance can only be absorbed into lambda_j",
      code_added_permutation_null_median_lambda = stats::median(od_null$lam_median),
      code_added_permutation_null_rho_hat_top12 = stats::median(od_null$rho_hat_top12)),
    code_side_definitions = list(
      note = "F1, F2 and F3 have no exact one-item decomposition of E4's set-level formulas; the per-item analogues are Code-side definitions carried to planning review as an open item. F4 and F5 decompose exactly.",
      F1 = "squared deviation of the target predictor from its pool mean",
      F2 = "eloss_a + eloss_b (E4's mean_eloss input); E4's loss-informative COUNT is constant 1 on the primary item universe",
      F3 = "max over predictor pairs of |zc_p * zc_q| on pool-centered z",
      baseline = "item-level baseline is broad category only (domain + near_parity + dominance); E4's problem count has no item-level analogue")), extra)
  writeLines(jsonlite::toJSON(gate, auto_unbox = TRUE, pretty = TRUE,
                              digits = 15, null = "null"),
             file.path(res_dir, "p2_gate.json"))
  gate
}

cfg <- list(
  spec = "elevation/P2_spec.md",
  raw_md5 = raw_md5, input_md5 = as.list(in_md5),
  regime = "experience", trial_window = c(6L, 25L), r_reps = R_REPS,
  n_build = n_build, n_freeze_not_read = 137L, n_test_sealed = 205L,
  n_clean_problems = 105L, n_loss_informative = length(li_games),
  forms_disjoint = TRUE,
  primary_engine = "mirt::mirt custom Binomial(20) item (EM/MML), SE from the observed information",
  confirmation_engine = "independent MML item-factor (Golub-Welsch Gauss-Hermite, analytic gradient, BFGS)",
  bayesian_engine = "NOT RUN: brms/Stan needs a C++ toolchain absent on this machine (E8 finding); rstanarm cannot express item-specific loadings",
  simulation_engine = "confirmation MML engine (a mirt fit per replicate is not tractable); licensed by the engine-agreement check",
  quadrature_nodes = NQ, seed = SEED_P2,
  n_sim_replicates = N_SIM_REP, n_permutations = N_PERM, n_boot = N_BOOT,
  thresholds = list(degenerate_max_fraction = DEGEN_MAX_FRAC,
                    lambda_abs_max = LAMBDA_ABS_MAX, se_mult_max = SE_MULT_MAX,
                    near_zero_var = NEAR_ZERO_VAR, condition_max = COND_MAX,
                    g0d_recovery_min = G0D_RECOV_MIN, g1_cv_min = G1_CV_MIN,
                    g3_min = G3_MIN),
  float_precision_note = "CSVs are written by readr::write_csv (shortest round-trip representation); use utils::read.csv for bit-level checks (the E3b caveat)",
  no_freeze_or_test_read = TRUE, computes_optimization_or_designs = FALSE,
  exploratory = TRUE, frozen_h2_fail_stands = TRUE)
writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, digits = 15),
           file.path(res_dir, "p2_config.json"))

mech_rows <- function(mr, track, target) tibble::tibble(
  track = track, target = target, term = names(mr$beta),
  partial_beta_standardized = unname(mr$beta), partial_se = unname(mr$se_beta),
  marginal_spearman_raw_feature = unname(
    mr$marginal_spearman[match(sub("_resid$", "", names(mr$beta)), FEATURES)]),
  beta_pc1_info_factor = mr$beta_pc1,
  cv_spearman_grouped = mr$cv_spearman,
  cv_spearman_baseline = mr$cv_spearman_baseline,
  cv_delta_vs_baseline = mr$cv_delta,
  n_items = mr$n_items, condition_number = mr$cond)
mech_sec <- mech_regression(FITS[["ev_sens_agn"]]$a, FEAT[["ev_sens_agn"]],
                            "real_secondary")
mechanism_regression <- dplyr::bind_rows(
  mech_rows(mech_real, "real_build", PRIMARY_TARGET),
  mech_rows(mech_sec, "real_build", "ev_sens_agn"))

if (!g0) {
  say("\nG0 FAILED: writing the identifiability result and STOPPING before G1-G3 (spec section 4).\n")
  save_csv(mechanism_regression |>
             dplyr::mutate(phase_1a_only_not_a_gate_result = TRUE),
           "mechanism_regression")
  write_gate("1a_stopped_on_G0_fail",
             list(G1 = "not run (G0 FAIL)", G2 = "not run (G0 FAIL)",
                  G3 = "not run (G0 FAIL)",
                  honesty_guards = "not run (G0 FAIL)",
                  interpretation = "G0 FAIL: the joint mechanism is not identifiable on CPC18 BUILD at these budgets. A clean pre-declared methodological null that REINFORCES the measurement-limits contribution: even a principled joint model cannot extract the mechanism here. The frozen H2 FAIL (E4) stands; Study 1, E8 and P3 are untouched; the 205 test subjects remain sealed."))
  say("\nP2 complete (G0 FAIL branch) in %.1f min\n",
      as.numeric(difftime(Sys.time(), t_start, units = "mins")))
  parallel::stopCluster(cl)
  quit(save = "no", status = 0)
}

# =============================================================================
# 8. Phase 1b: G1 / G2 / G3 + HG1-HG4 (only if G0 passes)
# =============================================================================
say("\n=== PHASE 1b: G1 / G2 / G3 + HG1-HG4 ===\n")

# G1: mechanism on the clean simulation testbed (planted map), out-of-fold.
g1_cv <- stats::median(gv("recovery", "cv_spearman"), na.rm = TRUE)
g1_delta <- stats::median(gv("recovery", "cv_delta"), na.rm = TRUE)
g1_pass <- is.finite(g1_cv) && g1_cv >= G1_CV_MIN && is.finite(g1_delta) &&
  g1_delta > 0
say("G1 (simulation mechanism): median CV Spearman %.4f (>= %.2f), delta vs baseline %.4f (> 0) -> %s\n",
    g1_cv, G1_CV_MIN, g1_delta, if (g1_pass) "PASS" else "FAIL")

# G2: compare the same relationship on empirical build data with G1.
g2_cv <- mech_real$cv_spearman
g2_sign_agree <- sign(mech_real$beta[F4_RESID]) ==
  sign(stats::median(gv("recovery", "beta_planted_feature"), na.rm = TRUE))
say("G2 (real BUILD honesty check): CV Spearman %.4f ; sign agreement with G1 = %s\n",
    g2_cv, g2_sign_agree)

# G3: calibration of the IRT-implied set reliability against E3b.
e3b_sub <- readr::read_csv(path_of[["e3b_subjsplit"]], show_col_types = FALSE)
e3b_out <- readr::read_csv(path_of[["e3b_outcome"]], show_col_types = FALSE)
calib_rows <- list()
for (tg in TARGETS) {
  lamv <- if (tg == PRIMARY_TARGET) lam else
    setNames(FITS[["ev_sens_agn"]]$a, as.character(all_games))
  for (K in K_REPORT) {
    sub <- e3b_sub |> dplyr::filter(target == tg, K == !!K)
    if (nrow(sub) == 0L) next
    rh <- rho_hat_set(lamv, names(lamv), sub$game_ids)
    obs <- (sub$r_sb_A + sub$r_sb_B) / 2
    ok <- is.finite(rh) & is.finite(obs)
    ceil <- sb_k(K, R1[[tg]])
    set.seed(seed_boot(match(tg, TARGETS)) + K)
    bt <- replicate(N_BOOT, {
      i <- sample.int(sum(ok), replace = TRUE)
      suppressWarnings(stats::cor(rh[ok][i], obs[ok][i], method = "spearman"))
    })
    calib_rows[[length(calib_rows) + 1L]] <- tibble::tibble(
      target = tg, K = K, n_sets = sum(ok), gating = K %in% K_CALIB,
      spearman_rho_hat_vs_e3b = suppressWarnings(
        stats::cor(rh[ok], obs[ok], method = "spearman")),
      ci_lo = unname(stats::quantile(bt, 0.025, na.rm = TRUE)),
      ci_hi = unname(stats::quantile(bt, 0.975, na.rm = TRUE)),
      rho_hat_median = stats::median(rh[ok]),
      rho_hat_q05 = unname(stats::quantile(rh[ok], 0.05)),
      rho_hat_q95 = unname(stats::quantile(rh[ok], 0.95)),
      e3b_observed_median = stats::median(obs[ok]),
      rel_ceiling_K = ceil,
      frac_sets_above_ceiling = mean(rh[ok] > ceil),
      hg1b_ci_supported_exceedance =
        unname(stats::quantile(rh[ok], 0.05)) > ceil)
  }
}
calibration <- dplyr::bind_rows(calib_rows)
save_csv(calibration, "calibration")
g3_tbl <- calibration |> dplyr::filter(target == PRIMARY_TARGET, gating)
g3_pass <- nrow(g3_tbl) > 0 && all(g3_tbl$spearman_rho_hat_vs_e3b >= G3_MIN)
say("G3 (calibration): Spearman by K %s -> %s\n",
    paste(sprintf("K%d=%.4f", g3_tbl$K, g3_tbl$spearman_rho_hat_vs_e3b),
          collapse = " "), if (g3_pass) "PASS" else "FAIL")

# HG1a: cross-form anchor guard.
demo <- build_raw |> dplyr::distinct(SubjID, Set, Age, Gender) |>
  dplyr::group_by(Set) |>
  dplyr::summarise(n = dplyr::n(), age_mean = mean(Age, na.rm = TRUE),
                   age_sd = stats::sd(Age, na.rm = TRUE),
                   frac_female = mean(Gender == 1, na.rm = TRUE), .groups = "drop")
form_sign <- lapply(sort(unique(F_pri$form)), function(f) {
  idx <- which(F_pri$form == f)
  sub_ids <- which(subj_form$Set == f)
  Yf <- Y_primary[sub_ids, idx, drop = FALSE]
  ff <- tryCatch(fit_mml(Yf, se = FALSE), error = function(e) NULL)
  mspear <- if (is.null(ff) || length(idx) < 4L) NA_real_ else
    suppressWarnings(stats::cor(F_pri$F4_expected_information[idx], ff$a,
                                method = "spearman"))
  tibble::tibble(form = f, n_items = length(idx), n_subjects = length(sub_ids),
                 within_form_marginal_spearman_F4 = mspear,
                 pooled_marginal_spearman_F4 = unname(
                   mech_real$marginal_spearman["F4_expected_information"]),
                 sign_agrees_with_pooled = !is.na(mspear) &&
                   sign(mspear) == sign(unname(
                     mech_real$marginal_spearman["F4_expected_information"])))
}) |> dplyr::bind_rows()
crossform_sensitivity <- dplyr::left_join(
  form_sign, demo |> dplyr::rename(form = Set), by = "form")
save_csv(crossform_sensitivity, "crossform_sensitivity")
hg1a_ok <- mean(crossform_sensitivity$sign_agrees_with_pooled, na.rm = TRUE) > 0.5
say("HG1a: within-form sign agreement with pooled = %.3f -> %s\n",
    mean(crossform_sensitivity$sign_agrees_with_pooled, na.rm = TRUE),
    if (hg1a_ok) "ok" else "FLAGGED")

hg1b_flag <- any(calibration$hg1b_ci_supported_exceedance[
  calibration$target == PRIMARY_TARGET], na.rm = TRUE)
say("HG1b ceiling guard: CI-supported exceedance of the E0 ceiling = %s\n", hg1b_flag)

# HG2: null permutation of the feature rows across items.
perm_tbl <- dplyr::bind_rows(parallel::parLapply(cl, seq_len(N_PERM), function(rep) {
  set.seed(seed_perm(rep))
  Fp <- F_pri
  o <- sample.int(nrow(Fp))
  Fp[, FEATURES] <- Fp[o, FEATURES]
  mr <- mech_regression(lam, Fp, "hg2_permuted")
  tibble::tibble(rep = rep, seed = seed_perm(rep),
                 beta_permuted = unname(mr$beta[F4_RESID]),
                 cv_spearman_permuted = mr$cv_spearman)
}))
save_csv(perm_tbl, "null_permutation")
hg2_beta_ci <- q90(perm_tbl$beta_permuted)
hg2_cv_ci <- q90(perm_tbl$cv_spearman_permuted)
hg2_ok <- hg2_beta_ci[1] <= 0 && 0 <= hg2_beta_ci[2] &&
  hg2_cv_ci[1] <= 0 && 0 <= hg2_cv_ci[2]
say("HG2 permutation: beta 90%% [%.4f, %.4f] ; CV Spearman 90%% [%.4f, %.4f] -> %s\n",
    hg2_beta_ci[1], hg2_beta_ci[2], hg2_cv_ci[1], hg2_cv_ci[2],
    if (hg2_ok) "ok" else "FLAGGED")

# HG3: the real-track predictive correlation vs the E3b outcome ceiling.
hg3 <- lapply(K_REPORT, function(K) {
  o <- e3b_out |> dplyr::filter(target == PRIMARY_TARGET, K == !!K)
  cal <- calibration |> dplyr::filter(target == PRIMARY_TARGET, K == !!K)
  if (nrow(o) == 0L || nrow(cal) == 0L) return(NULL)
  tibble::tibble(K = K, e3b_outcome_reliability = o$outcome_reliability,
                 ceiling_sqrt = sqrt(max(o$outcome_reliability, 0)),
                 observed_predictive_spearman = cal$spearman_rho_hat_vs_e3b,
                 exceeds_ceiling = cal$spearman_rho_hat_vs_e3b >
                   sqrt(max(o$outcome_reliability, 0)))
}) |> dplyr::bind_rows()
save_csv(hg3, "outcome_ceiling_guard")
hg3_flag <- any(hg3$exceeds_ceiling, na.rm = TRUE)
say("HG3 outcome ceiling: exceedance = %s\n", hg3_flag)

hg4_ok <- stats::sd(lam) > 1e-6 && length(unique(round(lam, 8))) > 1L
say("HG4 non-degeneracy: sd(lambda) %.6f ; distinct loadings %d -> %s\n",
    stats::sd(lam), length(unique(round(lam, 8))), if (hg4_ok) "ok" else "FLAGGED")

save_csv(mechanism_regression, "mechanism_regression")
write_gate("1b_complete", list(
  G1 = list(pass = g1_pass, cv_spearman = g1_cv, threshold = G1_CV_MIN,
            cv_delta_vs_baseline = g1_delta),
  G2 = list(cv_spearman_real = g2_cv, sign_agrees_with_G1 = unname(g2_sign_agree),
            gating = FALSE),
  G3 = list(pass = g3_pass, threshold = G3_MIN,
            per_K = lapply(seq_len(nrow(g3_tbl)), function(i)
              as.list(g3_tbl[i, c("K", "spearman_rho_hat_vs_e3b", "ci_lo", "ci_hi")]))),
  honesty_guards = list(HG1a_ok = hg1a_ok, HG1b_leakage_flag = hg1b_flag,
                        HG2_ok = hg2_ok, HG3_leakage_flag = hg3_flag,
                        HG4_ok = hg4_ok)))

parallel::stopCluster(cl)
say("\nP2 complete in %.1f min\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))

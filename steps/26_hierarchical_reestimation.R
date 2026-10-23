# Step 26: hierarchical re-estimation and reliability by budget
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e0/, e1/ (steps 19-20)
# Writes:  output/elevation/results/e8/ (reliability_by_budget.csv, bootstrap, permutation and shrinkage tables)
# Paper:   The reliability ceiling reflects true-score variance (Better estimation does not lower the ceiling); Abstract
# Runtime: about 10 minutes with 12 parallel workers

# E8: Compare hierarchical partial pooling to frozen per-subject estimation on the
# same split draws. The main interest is loss_avoid_adj; ev_sens_agn is of secondary
# interest. The declared outcome relating to the discussion of the upper bound on true
# variance is a small change in reliability.
# The 344 build subjects from step 20 are used for all estimation. For the full sample
# data, the choice and outcome columns are removed including the choices of the 137
# freeze and 205 test subjects. Note that this does not attempt to optimize candidate
# designs, fit a mechanism model, or read in freeze/test outcomes.
# The following functions are loaded from their original files and used as is:
# R/step2_reliability.R build_subject_data, fit_loss_baselines,
# estimate_half, split_correlations, sb_correct,
# reliability_point
# R/step1_measures.R index definitions (used for provenance)
# R/cpc18_features.R load_cpc18, is_clean_row, clean_problem_features,
# add_scaled_features
# Partitioning and eligibility are drawn in the same way as for run_cell with lambda =
# NA in the draw.
# The raw code calls estimate_half directly, without an accelerator. G_EQ re-derives
# estimates for a committed sample of halves and checks that they match exactly before
# reliability computation. Both code paths use split_correlations, sb_correct, and
# reliability_point for interval and point estimation. The difference between the two
# code paths is in the vector of per-subject/per-half values passed in.
# When E8_SMOKE=1, fewer iterations are run and output is placed in a scratch folder
# instead of output/elevation/. The smoke run numbers are not part of the research.
# The constants for full runs are unchanged.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
  library(lme4)
})

source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))
source(here("R", "step1_measures.R"))

SMOKE <- identical(Sys.getenv("E8_SMOKE"), "1")

res_dir <- if (SMOKE) {
  # persists after exit so smoke mechanics can be inspected; never output/elevation/
  sd0 <- Sys.getenv("E8_SMOKE_DIR")
  if (nzchar(sd0)) sd0 else file.path(tempdir(), "e8_smoke")
} else here("output", "elevation", "results", "e8")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
if (!SMOKE) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  invisible(df)
}

t_start <- Sys.time()
say <- function(...) cat(sprintf(...), sep = "")
if (SMOKE) say("*** SMOKE RUN: reduced constants, outputs to %s - NOT A RESULT ***\n",
               res_dir)

# ---- frozen constants (fixed BEFORE any result exists) ----
RAW_MD5        <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_E8        <- 20260722L
REGIME         <- "exp"
POS            <- 6:25                  # experience regime, r = 20 (full pool)
R_REPS         <- length(POS)
B_GRID         <- c(6L, 8L, 10L, 12L)   # problems per half
TARGETS        <- c("loss_avoid_adj", "ev_sens_agn")   # primary, secondary
N_SPLITS       <- if (SMOKE) 6L else 100L              # N_SPLITS_HIER
MIN_PAIRS      <- 10L                   # frozen split_correlations default
N_BOOT         <- if (SMOKE) 40L else 1000L
GLMER_MAXFUN1  <- 1e5
GLMER_MAXFUN2  <- 1e6                   # refit-once budget on non-convergence
CONV_FAIL_MAX  <- 0.05                  # reported target
G1_DELTA_MIN   <- 0.05                  # SB units
G1_BUDGETS     <- c(8L, 10L, 12L)       # K = 6 reported but thin, non-gating
PRIMARY_TARGET <- "loss_avoid_adj"
PRIMARY_B      <- 12L
N_SPLITS_BAYES <- if (SMOKE) 2L else 20L
STAN_CHAINS    <- if (SMOKE) 2L else 4L
STAN_ITER      <- if (SMOKE) 500L else 2000L
STAN_ITER_RETRY<- if (SMOKE) 1000L else 4000L
G3_ABS_TOL     <- 0.05
RHAT_MAX       <- 1.01
ESS_MIN        <- 400
N_SAMPLE_SPLITS<- 3L                    # committed per-subject sample: splits 1-3
N_WORKERS      <- if (SMOKE) 6L else 12L
HG2_ABS_MAX    <- 0.02                  # 'materially non-zero' guard: |r| and
HG2_Z_MAX      <- 3                     # |z| must BOTH exceed to flag

# Seed derivations (all recorded in e8_seeds.csv). t_idx: 1 = loss_avoid_adj,
# 2 = ev_sens_agn; arm: 1 = raw, 2 = hier.
seed_draws <- function(B) SEED_E8 + 1000L * B
seed_boot  <- function(t_idx, B) SEED_E8 + 3000000L + 100000L * t_idx + 1000L * B
seed_perm  <- function(t_idx, B) SEED_E8 + 6000000L + 100000L * t_idx + 1000L * B
seed_relpt <- function(t_idx, B, arm)
  SEED_E8 + 9000000L + 1000000L * (arm - 1L) + 100000L * t_idx + 1000L * B
seed_stan  <- function(split, half) SEED_E8 + 12000000L + 10L * split + half

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
stopifnot(n_build == 344L, n_freeze == 137L, n_test == 205L)
stopifnot(!anyDuplicated(build_ids))
say("split: build %d / freeze %d (not read) / test %d (sealed)\n",
    n_build, n_freeze, n_test)

raw <- load_cpc18(raw_path)
stopifnot(all(cpc18_expected_columns %in% names(raw)))

# Structural guarantee: the full-sample frame keeps no
# choice/outcome column, so no freeze or test choice can reach any computation.
choice_cols <- c("B", "Button", "Payoff", "Forgone", "RT", "Apay", "Bpay",
                 "Feedback")
design_all <- raw |> dplyr::select(-dplyr::any_of(choice_cols))
stopifnot(!any(choice_cols %in% names(design_all)))

build_raw <- raw |> dplyr::filter(SubjID %in% build_ids)
rm(raw); invisible(gc())
stopifnot(setequal(unique(build_raw$SubjID), build_ids))
stopifnot(!any(build_raw$SubjID %in% held_out))

# Problem features from the choice-free frame (frozen loaders); no other
# results file is read (raw + E1 split + E0 anchors ONLY).
feat_all <- clean_problem_features(design_all) |> add_scaled_features()
feat_all <- feat_all |> dplyr::arrange(GameID)
stopifnot(nrow(feat_all) == 105L)
say("item universe: %d clean problems (frozen loaders)\n", nrow(feat_all))

build_trials <- build_raw |>
  dplyr::filter(is_clean_row(build_raw)) |>
  dplyr::select(SubjID, Set, GameID, Trial, B) |>
  dplyr::left_join(
    feat_all |> dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev,
                              loss_informative, has_loss_b, ev_a, ev_b),
    by = "GameID")
stopifnot(!any(is.na(build_trials$z_ev)))
stopifnot(dplyr::n_distinct(build_trials$SubjID) == 344L)

# Frozen build-only loss baselines and per-subject data structures.
baselines <- fit_loss_baselines(build_trials)
subj_data <- build_subject_data(build_trials, baselines)
stopifnot(length(subj_data) == 344L)
say("frozen loss baselines (exp): b0=%.6f b1=%.6f\n",
    baselines$exp[1], baselines$exp[2])

# =============================================================================
# 2. E0 ceiling anchors (read-only; HG1 guard reference + gap_closed)
# =============================================================================
e0 <- readr::read_csv(
  here("output", "elevation", "results", "e0", "variance_ceiling_final.csv"),
  show_col_types = FALSE)
ceil_row <- function(tg) e0[e0$construct == tg, ]
# Implied single-problem reliability r1 from the E0 full-data ceiling at
# k_ref = 14 via Spearman-Brown in problem count: r1 = R / (k - (k-1) R).
r1_of <- function(R, k) R / (k - (k - 1) * R)
sb_k  <- function(k, r1) k * r1 / (1 + (k - 1) * r1)
R1 <- vapply(TARGETS, function(tg) {
  rw <- ceil_row(tg)
  stopifnot(nrow(rw) == 1L, rw$k_ref == 14)
  r1 <- r1_of(rw$ceiling_fulldata, 14)
  # The ceiling file's k_0.7 / k_0.8 columns must reproduce under the same r1
  # (for loss: 27.31 / 46.82).
  k_needed <- function(R) R * (1 - r1) / (r1 * (1 - R))
  stopifnot(abs(k_needed(0.7) - rw$k_0.7) < 1e-9,
            abs(k_needed(0.8) - rw$k_0.8) < 1e-9)
  r1
}, numeric(1))
stopifnot(abs(round(R1["loss_avoid_adj"], 4) - 0.0787) < 1e-12)
stopifnot(abs(round(R1["ev_sens_agn"],   5) - 0.02717) < 1e-12)
say("E0 anchors: r1 loss %.6f, ev %.6f (k_0.7/k_0.8 reproduce exactly)\n",
    R1["loss_avoid_adj"], R1["ev_sens_agn"])

# =============================================================================
# 3. Eligibility + the committed split draws (frozen run_cell scheme,
#    lambda = NA, shared by BOTH arms; HG3 disjointness asserted)
# =============================================================================
# run_cell eligibility: a subject is eligible at budget B iff they saw >= 2B
# clean problems. The partition mirrors run_cell's lambda = NA branch
# draw-for-draw: one common pos per split (r = 20 = the full experience pool,
# so pos is the pool and consumes no RNG), then an independent
# sample(nrow(feat), 2B) per subject; h1 = first B, h2 = last B.
elig_of <- function(B) names(subj_data)[vapply(
  subj_data, function(sd) nrow(sd$feat) >= 2L * B, logical(1))]

draw_cell <- function(B) {
  elig <- elig_of(B)
  n <- length(elig)
  set.seed(seed_draws(B))
  pool <- 6:25
  draws <- vector("list", N_SPLITS)
  for (s in seq_len(N_SPLITS)) {
    pos <- if (R_REPS < length(pool)) sort(sample(pool, R_REPS)) else pool
    M <- matrix(0L, n, 2L * B)
    for (i in seq_len(n))
      M[i, ] <- sample(seq_len(nrow(subj_data[[elig[i]]]$feat)), 2L * B)
    stopifnot(all(vapply(seq_len(n),
                         function(i) anyDuplicated(M[i, ]) == 0L, logical(1))))
    draws[[s]] <- M
  }
  list(B = B, elig = elig, n = n, draws = draws)
}
cells <- lapply(B_GRID, draw_cell)
names(cells) <- as.character(B_GRID)
for (ce in cells)
  say("budget B=%2d: eligible N = %3d (%.3f of build)\n",
      ce$B, ce$n, ce$n / n_build)

# =============================================================================
# 4. Worker machinery: one job = one (budget, split): RAW estimates via the
#    FROZEN estimate_half + 4 hier glmer fits (2 halves x 2 targets)
# =============================================================================
fit_glmer_guarded <- function(form, dat, w, maxfun) {
  wrn <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      lme4::glmer(form, data = dat, family = stats::binomial(), weights = w,
                  control = lme4::glmerControl(
                    optimizer = "bobyqa", optCtrl = list(maxfun = maxfun),
                    calc.derivs = TRUE)),
      error = function(e) e),
    warning = function(cond) {
      wrn <<- c(wrn, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
  ok <- inherits(fit, "merMod") &&
    (is.null(fit@optinfo$conv$opt) || fit@optinfo$conv$opt == 0) &&
    !any(grepl("failed to converge",
               c(unlist(fit@optinfo$conv$lme4$messages), wrn),
               ignore.case = TRUE))
  list(fit = fit, converged = ok)
}

# Fit one half's hier model for one target; refit once on non-convergence.
# Returns per-eligible-subject estimates (NA when absent or failed) plus
# diagnostics. target 'loss': random intercept; 'ev': random slope on z_ev.
hier_fit_half <- function(frame, elig, target) {
  out <- list(est = rep(NA_real_, length(elig)), converged = FALSE,
              retried = FALSE, singular = NA, re_sd_target = NA_real_,
              re_sd_other = NA_real_, re_corr = NA_real_,
              fix_int = NA_real_, fix_slope = NA_real_,
              n_rows = nrow(frame), n_subj = dplyr::n_distinct(frame$SubjID))
  if (nrow(frame) == 0L) return(out)
  frame$SubjID <- factor(frame$SubjID, levels = elig)
  form <- if (target == "loss_avoid_adj")
    y ~ 1 + d_ev_noloss + (1 | SubjID) else
    y ~ 1 + z_ev + (1 + z_ev | SubjID)
  r1 <- fit_glmer_guarded(form, frame, frame$w, GLMER_MAXFUN1)
  res <- r1
  if (!r1$converged) {
    res <- fit_glmer_guarded(form, frame, frame$w, GLMER_MAXFUN2)
    out$retried <- TRUE
  }
  if (!res$converged) return(out)          # NA estimates, counted as failure
  fit <- res$fit
  out$converged <- TRUE
  out$singular <- lme4::isSingular(fit, tol = 1e-4)
  vc <- lme4::VarCorr(fit)$SubjID
  fx <- lme4::fixef(fit)
  re <- lme4::ranef(fit)$SubjID
  m <- match(rownames(re), elig)
  if (target == "loss_avoid_adj") {
    out$re_sd_target <- sqrt(unname(vc[1, 1]))
    out$fix_int <- unname(fx[1]); out$fix_slope <- unname(fx["d_ev_noloss"])
    out$est[m] <- re[["(Intercept)"]]
  } else {
    sds <- attr(vc, "stddev")
    out$re_sd_target <- unname(sds["z_ev"])
    out$re_sd_other  <- unname(sds["(Intercept)"])
    cr <- attr(vc, "correlation")
    out$re_corr <- if (!is.null(dim(cr))) unname(cr[2, 1]) else NA_real_
    out$fix_int <- unname(fx[1]); out$fix_slope <- unname(fx["z_ev"])
    out$est[m] <- re[["z_ev"]]
  }
  out
}

# Build the two aggregated per-half frames for one subject and half rows.
half_frames <- function(sid, rows) {
  sd <- subj_data[[sid]]
  f <- sd$feat[rows, ]
  nb <- if (length(rows) == 1L) sum(sd$B_mat[rows, POS]) else
    rowSums(sd$B_mat[rows, POS, drop = FALSE])
  ev <- data.frame(SubjID = sid, GameID = f$GameID,
                   y = nb / R_REPS, z_ev = f$z_ev, w = R_REPS)
  li <- which(f$loss_informative)
  ls <- if (length(li) == 0L) NULL else data.frame(
    SubjID = sid, GameID = f$GameID[li],
    y = ifelse(f$has_loss_b[li], R_REPS - nb[li], nb[li]) / R_REPS,
    d_ev_noloss = f$d_ev_noloss[li], w = R_REPS)
  list(ev = ev, loss = ls)
}

run_split_job <- function(job) {
  elig <- job$elig; B <- job$B; M <- job$M
  n <- length(elig)
  raw <- array(NA_real_, dim = c(n, 2L, 2L),
               dimnames = list(NULL, c("h1", "h2"), TARGETS))
  fr <- list(h1 = list(ev = vector("list", n), loss = vector("list", n)),
             h2 = list(ev = vector("list", n), loss = vector("list", n)))
  for (i in seq_len(n)) {
    sd <- subj_data[[elig[i]]]
    h1 <- M[i, seq_len(B)]; h2 <- M[i, (B + 1L):(2L * B)]
    e1 <- estimate_half(sd, h1, POS, REGIME)                      # FROZEN
    e2 <- estimate_half(sd, h2, POS, REGIME)                      # FROZEN
    raw[i, 1L, ] <- e1[TARGETS]; raw[i, 2L, ] <- e2[TARGETS]
    f1 <- half_frames(elig[i], h1); f2 <- half_frames(elig[i], h2)
    fr$h1$ev[[i]] <- f1$ev; fr$h1$loss[[i]] <- f1$loss
    fr$h2$ev[[i]] <- f2$ev; fr$h2$loss[[i]] <- f2$loss
  }
  ev1 <- dplyr::bind_rows(fr$h1$ev);  ev2 <- dplyr::bind_rows(fr$h2$ev)
  ls1 <- dplyr::bind_rows(fr$h1$loss); ls2 <- dplyr::bind_rows(fr$h2$loss)
  # HG3: the two halves share NO subject-problem cell (disjoint by draw;
  # asserted here on the actual fit frames, ev frame = the half's full set).
  overlap <- nrow(dplyr::inner_join(ev1[, c("SubjID", "GameID")],
                                    ev2[, c("SubjID", "GameID")],
                                    by = c("SubjID", "GameID")))
  hier <- list(
    loss = list(h1 = hier_fit_half(ls1, elig, "loss_avoid_adj"),
                h2 = hier_fit_half(ls2, elig, "loss_avoid_adj")),
    ev   = list(h1 = hier_fit_half(ev1, elig, "ev_sens_agn"),
                h2 = hier_fit_half(ev2, elig, "ev_sens_agn")))
  list(split = job$split, raw = raw, hier = hier, hg3_overlap = overlap)
}

say("\nstarting cluster (%d workers)\n", N_WORKERS)
cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(tidyverse); library(here); library(lme4)
  })
  source(here("R", "cpc18_features.R"))
  source(here("R", "step2_reliability.R"))
  source(here("R", "step1_measures.R"))
  NULL
}))
parallel::clusterExport(cl, c(
  "subj_data", "POS", "REGIME", "R_REPS", "TARGETS",
  "GLMER_MAXFUN1", "GLMER_MAXFUN2",
  "fit_glmer_guarded", "hier_fit_half", "half_frames", "run_split_job"),
  envir = environment())

# =============================================================================
# 5. Main loop: per budget, all splits; assemble estimate matrices
# =============================================================================
EST <- list()      # EST[[target]][[as.character(B)]]$raw$h1 etc: n x N_SPLITS
DIAG <- list()     # per-fit diagnostics rows
HG3_OVERLAPS <- integer(0)
for (ce in cells) {
  B <- ce$B
  jobs <- lapply(seq_len(N_SPLITS), function(s)
    list(split = s, B = B, elig = ce$elig, M = ce$draws[[s]]))
  t0 <- Sys.time()
  res <- parallel::parLapply(cl, jobs, run_split_job)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ord <- order(vapply(res, function(z) z$split, integer(1)))
  res <- res[ord]
  HG3_OVERLAPS <- c(HG3_OVERLAPS,
                    vapply(res, function(z) z$hg3_overlap, integer(1)))
  for (tg in TARGETS) {
    key <- if (tg == "loss_avoid_adj") "loss" else "ev"
    EST[[tg]][[as.character(B)]] <- list(
      raw = list(
        h1 = vapply(res, function(z) z$raw[, 1L, tg], numeric(ce$n)),
        h2 = vapply(res, function(z) z$raw[, 2L, tg], numeric(ce$n))),
      hier = list(
        h1 = vapply(res, function(z) z$hier[[key]]$h1$est, numeric(ce$n)),
        h2 = vapply(res, function(z) z$hier[[key]]$h2$est, numeric(ce$n))))
    for (s in seq_len(N_SPLITS)) for (hf in c("h1", "h2")) {
      hh <- res[[s]]$hier[[key]][[hf]]
      DIAG[[length(DIAG) + 1L]] <- tibble::tibble(
        target = tg, budget = B, split = s, half = hf,
        converged = hh$converged, retried = hh$retried,
        singular = hh$singular, re_sd_target = hh$re_sd_target,
        re_sd_other = hh$re_sd_other, re_corr = hh$re_corr,
        fix_int = hh$fix_int, fix_slope = hh$fix_slope,
        n_rows = hh$n_rows, n_subj_in_fit = hh$n_subj)
    }
  }
  say("budget B=%2d done in %.1f min (%d splits x [%d raw halves + 4 glmer])\n",
      B, el / 60, N_SPLITS, 2L * ce$n)
}
diag_tbl <- dplyr::bind_rows(DIAG)
stopifnot(all(HG3_OVERLAPS == 0L))                                    # HG3
say("HG3: zero subject-problem overlap between half fits in all %d splits\n",
    length(HG3_OVERLAPS))

# =============================================================================
# 6. Gate G_EQ: re-derive the RAW arm for the committed sample splits
#    through the frozen estimate_half, serially, and assert BIT-IDENTITY
# =============================================================================
say("\n=== G_EQ: frozen estimate_half re-derivation (splits 1-%d) ===\n",
    N_SAMPLE_SPLITS)
geq_n <- 0L; geq_bad <- 0L
for (ce in cells) {
  B <- ce$B
  for (s in seq_len(min(N_SAMPLE_SPLITS, N_SPLITS))) {
    M <- ce$draws[[s]]
    for (i in seq_len(ce$n)) {
      sd <- subj_data[[ce$elig[i]]]
      e1 <- estimate_half(sd, M[i, seq_len(B)], POS, REGIME)
      e2 <- estimate_half(sd, M[i, (B + 1L):(2L * B)], POS, REGIME)
      for (tg in TARGETS) {
        stored <- c(EST[[tg]][[as.character(B)]]$raw$h1[i, s],
                    EST[[tg]][[as.character(B)]]$raw$h2[i, s])
        fresh <- c(e1[[tg]], e2[[tg]])
        if (!identical(is.na(stored), is.na(fresh)) ||
            any(stored != fresh, na.rm = TRUE)) geq_bad <- geq_bad + 1L
        geq_n <- geq_n + 1L
      }
    }
  }
}
say("G_EQ: %d subject-split-target pairs re-derived; mismatches = %d\n",
    geq_n, geq_bad)
stopifnot(geq_bad == 0L)

# =============================================================================
# 7. Reliability, paired Delta bootstrap, HG1 + HG2 (frozen engine throughout)
# =============================================================================
# Per cell-target: point reliabilities via frozen split_correlations +
# sb_correct; per-arm 95% CI via the frozen reliability_point; the PAIRED
# subject-cluster bootstrap resamples subjects ONCE per draw and pushes the
# SAME index through both arms (bootstrap outside, splits inside, as frozen).
cell_inference <- function(job) {
  R1 <- job$R1; R2 <- job$R2; H1 <- job$H1; H2 <- job$H2
  n <- nrow(R1)
  rs_raw  <- split_correlations(R1, R2, MIN_PAIRS)                    # FROZEN
  rs_hier <- split_correlations(H1, H2, MIN_PAIRS)                    # FROZEN
  point_raw  <- mean(rs_raw,  na.rm = TRUE)
  point_hier <- mean(rs_hier, na.rm = TRUE)
  rp_raw  <- reliability_point(R1, R2, job$n_boot, job$seed_rp_raw,
                               MIN_PAIRS)                             # FROZEN
  rp_hier <- reliability_point(H1, H2, job$n_boot, job$seed_rp_hier,
                               MIN_PAIRS)                             # FROZEN
  set.seed(job$seed_boot)
  bt <- matrix(NA_real_, job$n_boot, 2L)
  for (b in seq_len(job$n_boot)) {
    idx <- sample.int(n, replace = TRUE)
    bt[b, 1L] <- sb_correct(mean(split_correlations(
      R1[idx, , drop = FALSE], R2[idx, , drop = FALSE], MIN_PAIRS),
      na.rm = TRUE))
    bt[b, 2L] <- sb_correct(mean(split_correlations(
      H1[idx, , drop = FALSE], H2[idx, , drop = FALSE], MIN_PAIRS),
      na.rm = TRUE))
  }
  dl <- bt[, 2L] - bt[, 1L]
  qs <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
  # HG2: fresh subject permutation per split, applied to the h2 vector,
  # breaking the true correspondence; frozen split_correlations on the
  # permuted pairing.
  set.seed(job$seed_perm)
  P2r <- R2; P2h <- H2
  for (s in seq_len(ncol(R2))) {
    pr <- sample.int(n); P2r[, s] <- R2[pr, s]
    ph <- sample.int(n); P2h[, s] <- H2[ph, s]
  }
  perm_raw  <- split_correlations(R1, P2r, MIN_PAIRS)
  perm_hier <- split_correlations(H1, P2h, MIN_PAIRS)
  list(
    rs_raw = rs_raw, rs_hier = rs_hier,
    point_raw = point_raw, point_hier = point_hier,
    rel_raw = sb_correct(point_raw), rel_hier = sb_correct(point_hier),
    rp_raw = rp_raw, rp_hier = rp_hier,
    boot_raw = bt[, 1L], boot_hier = bt[, 2L], boot_delta = dl,
    delta = sb_correct(point_hier) - sb_correct(point_raw),
    delta_ci90 = qs(dl, c(0.05, 0.95)), delta_ci95 = qs(dl, c(0.025, 0.975)),
    hier_ci90 = qs(bt[, 2L], c(0.05, 0.95)),
    perm = list(
      raw = list(mean = mean(perm_raw, na.rm = TRUE),
                 se = stats::sd(perm_raw, na.rm = TRUE) /
                   sqrt(sum(!is.na(perm_raw))),
                 n_used = sum(!is.na(perm_raw))),
      hier = list(mean = mean(perm_hier, na.rm = TRUE),
                  se = stats::sd(perm_hier, na.rm = TRUE) /
                    sqrt(sum(!is.na(perm_hier))),
                  n_used = sum(!is.na(perm_hier)))))
}

say("\n=== inference per (target, budget): frozen reliability + paired bootstrap (n_boot = %d) ===\n",
    N_BOOT)
inf_jobs <- list()
for (tg in TARGETS) for (ce in cells) {
  t_idx <- match(tg, TARGETS)
  e <- EST[[tg]][[as.character(ce$B)]]
  inf_jobs[[length(inf_jobs) + 1L]] <- list(
    target = tg, B = ce$B, n = ce$n,
    R1 = e$raw$h1, R2 = e$raw$h2, H1 = e$hier$h1, H2 = e$hier$h2,
    n_boot = N_BOOT,
    seed_boot = seed_boot(t_idx, ce$B), seed_perm = seed_perm(t_idx, ce$B),
    seed_rp_raw = seed_relpt(t_idx, ce$B, 1L),
    seed_rp_hier = seed_relpt(t_idx, ce$B, 2L))
}
parallel::clusterExport(cl, c("MIN_PAIRS", "cell_inference"),
                        envir = environment())
t0 <- Sys.time()
inf_res <- parallel::parLapply(cl, inf_jobs, cell_inference)
say("inference done in %.1f min\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

# =============================================================================
# 8. Assemble output tables
# =============================================================================
rel_rows <- list(); boot_rows <- list(); perm_rows <- list()
ceil_rows <- list(); split_rows <- list()
for (j in seq_along(inf_jobs)) {
  jb <- inf_jobs[[j]]; rz <- inf_res[[j]]
  tg <- jb$target; B <- jb$B
  ce <- cells[[as.character(B)]]
  ceiling_2B <- sb_k(2L * B, R1[tg])
  at_ceiling <- rz$rel_raw >= ceiling_2B
  gap_closed <- if (at_ceiling) NA_real_ else
    (rz$rel_hier - rz$rel_raw) / (ceiling_2B - rz$rel_raw)
  e <- EST[[tg]][[as.character(B)]]
  rel_rows[[j]] <- tibble::tibble(
    target = tg, budget_B = B, problems_total_2B = 2L * B,
    eligible_n = ce$n, eligible_frac_build = ce$n / n_build,
    n_splits = N_SPLITS,
    n_splits_used_raw = sum(!is.na(rz$rs_raw)),
    n_splits_used_hier = sum(!is.na(rz$rs_hier)),
    finite_frac_raw = mean(is.finite(e$raw$h1) & is.finite(e$raw$h2)),
    finite_frac_hier = mean(is.finite(e$hier$h1) & is.finite(e$hier$h2)),
    r_raw_uncorrected = rz$point_raw, rel_raw = rz$rel_raw,
    rel_raw_ci95_lo = rz$rp_raw$ci_lo_sb, rel_raw_ci95_hi = rz$rp_raw$ci_hi_sb,
    r_hier_uncorrected = rz$point_hier, rel_hier = rz$rel_hier,
    rel_hier_ci95_lo = rz$rp_hier$ci_lo_sb,
    rel_hier_ci95_hi = rz$rp_hier$ci_hi_sb,
    delta = rz$delta,
    delta_ci90_lo = rz$delta_ci90[1], delta_ci90_hi = rz$delta_ci90[2],
    delta_ci95_lo = rz$delta_ci95[1], delta_ci95_hi = rz$delta_ci95[2],
    rel_ceiling_2B = ceiling_2B, gap_closed = gap_closed,
    at_ceiling = at_ceiling, n_boot = N_BOOT)
  qb <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
  boot_rows[[j]] <- tibble::tibble(
    target = tg, budget_B = B, n_boot = N_BOOT,
    delta_point = rz$delta,
    boot_mean = mean(rz$boot_delta, na.rm = TRUE),
    boot_sd = stats::sd(rz$boot_delta, na.rm = TRUE),
    q025 = qb(rz$boot_delta, .025), q05 = qb(rz$boot_delta, .05),
    q25 = qb(rz$boot_delta, .25), q50 = qb(rz$boot_delta, .50),
    q75 = qb(rz$boot_delta, .75), q95 = qb(rz$boot_delta, .95),
    q975 = qb(rz$boot_delta, .975),
    frac_boot_delta_le_0 = mean(rz$boot_delta <= 0, na.rm = TRUE),
    frac_boot_hier_above_ceiling =
      mean(rz$boot_hier > ceiling_2B, na.rm = TRUE))
  for (arm in c("raw", "hier")) {
    pm <- rz$perm[[arm]]
    z <- pm$mean / pm$se
    perm_rows[[length(perm_rows) + 1L]] <- tibble::tibble(
      target = tg, budget_B = B, arm = arm,
      perm_r_mean = pm$mean, perm_r_se = pm$se, perm_z = z,
      perm_r_sb = sb_correct(pm$mean), n_splits_used = pm$n_used,
      materially_nonzero = abs(pm$mean) > HG2_ABS_MAX & abs(z) > HG2_Z_MAX)
  }
  ceil_rows[[j]] <- tibble::tibble(
    target = tg, budget_B = B, rel_hier = rz$rel_hier,
    rel_hier_ci90_lo = rz$hier_ci90[1], rel_hier_ci90_hi = rz$hier_ci90[2],
    rel_ceiling_2B = ceiling_2B,
    exceeds_point = rz$rel_hier > ceiling_2B,
    frac_boot_above = mean(rz$boot_hier > ceiling_2B, na.rm = TRUE),
    ci_supported_exceedance = rz$hier_ci90[1] > ceiling_2B)
  split_rows[[j]] <- tibble::tibble(
    target = tg, budget_B = B, split = seq_len(N_SPLITS),
    r_raw_split = rz$rs_raw, r_hier_split = rz$rs_hier,
    n_pairs_raw = vapply(seq_len(N_SPLITS), function(s)
      sum(is.finite(e$raw$h1[, s]) & is.finite(e$raw$h2[, s])), integer(1)),
    n_pairs_hier = vapply(seq_len(N_SPLITS), function(s)
      sum(is.finite(e$hier$h1[, s]) & is.finite(e$hier$h2[, s])), integer(1)))
}
reliability_by_budget <- dplyr::bind_rows(rel_rows) |>
  dplyr::arrange(target, budget_B)
delta_bootstrap <- dplyr::bind_rows(boot_rows) |>
  dplyr::arrange(target, budget_B)
null_permutation <- dplyr::bind_rows(perm_rows) |>
  dplyr::arrange(target, budget_B, arm)
ceiling_guard <- dplyr::bind_rows(ceil_rows) |>
  dplyr::arrange(target, budget_B)
per_split_correlations <- dplyr::bind_rows(split_rows) |>
  dplyr::arrange(target, budget_B, split)
print(as.data.frame(reliability_by_budget[, c(
  "target", "budget_B", "eligible_n", "rel_raw", "rel_hier", "delta",
  "delta_ci90_lo", "delta_ci90_hi", "rel_ceiling_2B", "gap_closed")]),
  digits = 4)

# =============================================================================
# 9. HG4 shrinkage diagnostics (+ convergence accounting)
# =============================================================================
# Cross-arm agreement and SD ratio per (split, half) over common finite
# subjects. NOTE (labeled in config): raw and hier live on different scales
# (probability-residual vs logit intercept; Spearman index vs logistic slope),
# so the SD ratio is a cross-scale collapse diagnostic, not a calibrated
# shrinkage estimate; the rank agreement (Spearman) is scale-free.
agree_rows <- list()
for (tg in TARGETS) for (ce in cells) {
  e <- EST[[tg]][[as.character(ce$B)]]
  rat <- c(); agr <- c()
  for (s in seq_len(N_SPLITS)) for (hf in c("h1", "h2")) {
    x <- e$raw[[hf]][, s]; y <- e$hier[[hf]][, s]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= MIN_PAIRS && stats::sd(x[ok]) > 0) {
      rat <- c(rat, stats::sd(y[ok]) / stats::sd(x[ok]))
      agr <- c(agr, suppressWarnings(
        stats::cor(x[ok], y[ok], method = "spearman")))
    }
  }
  agree_rows[[length(agree_rows) + 1L]] <- tibble::tibble(
    target = tg, budget_B = ce$B,
    sd_ratio_hier_over_raw_mean = mean(rat),
    spearman_raw_hier_mean = mean(agr, na.rm = TRUE),
    n_half_fits_compared = length(rat))
}
shrinkage_diagnostics <- diag_tbl |>
  dplyr::group_by(target, budget = budget) |>
  dplyr::summarise(
    n_fits = dplyr::n(),
    conv_fail_rate = mean(!converged),
    n_retried = sum(retried),
    singular_frac = mean(singular, na.rm = TRUE),
    re_sd_target_mean = mean(re_sd_target, na.rm = TRUE),
    re_sd_target_min = min(re_sd_target, na.rm = TRUE),
    re_sd_target_max = max(re_sd_target, na.rm = TRUE),
    re_sd_other_mean = mean(re_sd_other, na.rm = TRUE),
    re_corr_mean = mean(re_corr, na.rm = TRUE),
    fix_int_mean = mean(fix_int, na.rm = TRUE),
    fix_slope_mean = mean(fix_slope, na.rm = TRUE),
    re_sd_collapsed_frac = mean(re_sd_target < 1e-6, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::rename(budget_B = budget) |>
  dplyr::left_join(dplyr::bind_rows(agree_rows), by = c("target", "budget_B")) |>
  dplyr::arrange(target, budget_B)
print(as.data.frame(shrinkage_diagnostics[, c(
  "target", "budget_B", "conv_fail_rate", "singular_frac",
  "re_sd_target_mean", "sd_ratio_hier_over_raw_mean",
  "spearman_raw_hier_mean")]), digits = 4)

# =============================================================================
# 10. Committed per-subject sample (splits 1-3, both arms, both halves)
# =============================================================================
samp_rows <- list()
for (tg in TARGETS) for (ce in cells) {
  e <- EST[[tg]][[as.character(ce$B)]]
  for (s in seq_len(min(N_SAMPLE_SPLITS, N_SPLITS)))
    samp_rows[[length(samp_rows) + 1L]] <- tibble::tibble(
      target = tg, budget_B = ce$B, split = s, SubjID = ce$elig,
      raw_h1 = e$raw$h1[, s], raw_h2 = e$raw$h2[, s],
      hier_h1 = e$hier$h1[, s], hier_h2 = e$hier$h2[, s])
}
per_subject_estimates_sample <- dplyr::bind_rows(samp_rows)

# =============================================================================
# 11. G1 / G2 evaluation (pre-declared)
# =============================================================================
gate_eval <- reliability_by_budget |>
  dplyr::mutate(
    ci90_excludes_0 = delta_ci90_lo > 0 | delta_ci90_hi < 0,
    passes_g1_criteria = delta >= G1_DELTA_MIN & delta_ci90_lo > 0,
    in_g1_budgets = budget_B %in% G1_BUDGETS)
g1_tbl <- gate_eval |> dplyr::filter(target == PRIMARY_TARGET)
g1_pass <- any(g1_tbl$passes_g1_criteria & g1_tbl$in_g1_budgets)
g2_tbl <- gate_eval |> dplyr::filter(target == "ev_sens_agn")
g2_signal <- any(g2_tbl$passes_g1_criteria & g2_tbl$in_g1_budgets)
say("\nG1 (%s): %s\n", PRIMARY_TARGET, if (g1_pass) "PASS" else "FAIL")
say("G2 (ev_sens_agn, informational): %s\n",
    if (g2_signal) "signal present" else "no signal")

# =============================================================================
# 12. G3 Bayesian confirmation at the primary cell (rstanarm)
# =============================================================================
# ENGINE NOTE (deviation): the pre-declared engine is rstanarm::stan_glmer OR
# brms. brms needs a C++ toolchain for runtime Stan compilation; rstanarm's
# precompiled models are used instead.
# rstanarm cannot place a half-Normal prior on the random-intercept SD (its
# decov prior implies, for a single intercept, SD ~ Gamma(shape, scale)); the
# closest weakly-informative substitute Gamma(1, 2) = Exponential(mean 2) is
# used instead of half-Normal(0, 2) and the substitution is recorded in
# e8_config.json, with a prior-sensitivity spot check (scale 1 vs
# 2) saved alongside. Fixed effects keep the pre-declared Normal(0, 5).
say("\n=== G3: Bayesian confirmation (%s, B = %d, %d splits) ===\n",
    PRIMARY_TARGET, PRIMARY_B, N_SPLITS_BAYES)
ceP <- cells[[as.character(PRIMARY_B)]]
stan_frames <- list()
for (s in seq_len(N_SPLITS_BAYES)) {
  M <- ceP$draws[[s]]
  l1 <- list(); l2 <- list()
  for (i in seq_len(ceP$n)) {
    f1 <- half_frames(ceP$elig[i], M[i, seq_len(PRIMARY_B)])
    f2 <- half_frames(ceP$elig[i], M[i, (PRIMARY_B + 1L):(2L * PRIMARY_B)])
    l1[[i]] <- f1$loss; l2[[i]] <- f2$loss
  }
  stan_frames[[s]] <- list(h1 = dplyr::bind_rows(l1), h2 = dplyr::bind_rows(l2))
}
stan_job <- function(job) {
  suppressPackageStartupMessages(library(rstanarm))
  fr <- job$frame
  fr$SubjID <- factor(fr$SubjID, levels = job$elig)
  run_one <- function(iter) {
    fit <- stan_glmer(
      y ~ 1 + d_ev_noloss + (1 | SubjID), data = fr,
      family = stats::binomial(), weights = fr$w,
      prior = normal(0, 5), prior_intercept = normal(0, 5),
      prior_covariance = decov(regularization = 1, concentration = 1,
                               shape = job$decov_shape, scale = job$decov_scale),
      chains = job$chains, iter = iter, seed = job$seed,
      cores = 1L, refresh = 0)
    # as.array keeps rstanarm's pretty parameter names AND the chain
    # structure (fit$stanfit would expose raw Stan names and defeat the
    # variable filter below).
    dr <- posterior::as_draws_array(as.array(fit))
    vars <- posterior::variables(dr)
    keep <- grepl("^b\\[|^\\(Intercept\\)$|^d_ev_noloss$|^Sigma\\[", vars)
    ss <- posterior::summarise_draws(
      posterior::subset_draws(dr, variable = vars[keep]),
      "rhat", "ess_bulk")
    mm <- as.matrix(fit, regex_pars = "^b\\[\\(Intercept\\) SubjID:")
    ids <- sub("^b\\[\\(Intercept\\) SubjID:(.+)\\]$", "\\1", colnames(mm))
    est <- rep(NA_real_, length(job$elig))
    est[match(ids, job$elig)] <- colMeans(mm)
    list(est = est,
         max_rhat = max(ss$rhat, na.rm = TRUE),
         min_ess_bulk = min(ss$ess_bulk, na.rm = TRUE),
         divergences = rstan::get_num_divergent(fit$stanfit))
  }
  r <- run_one(job$iter)
  retried <- FALSE
  if (r$max_rhat >= job$rhat_max || r$min_ess_bulk <= job$ess_min ||
      r$divergences > 0) {
    r <- run_one(job$iter_retry); retried <- TRUE
  }
  c(r, list(retried = retried, split = job$split, half = job$half,
            decov_scale = job$decov_scale))
}
stan_jobs <- list()
for (s in seq_len(N_SPLITS_BAYES)) for (hf in 1:2) {
  stan_jobs[[length(stan_jobs) + 1L]] <- list(
    split = s, half = hf,
    frame = stan_frames[[s]][[c("h1", "h2")[hf]]], elig = ceP$elig,
    chains = STAN_CHAINS, iter = STAN_ITER, iter_retry = STAN_ITER_RETRY,
    seed = seed_stan(s, hf), rhat_max = RHAT_MAX, ess_min = ESS_MIN,
    decov_shape = 1, decov_scale = 2)
}
# Prior-sensitivity spot check: split 1, both halves, decov scale 1.
for (hf in 1:2) {
  stan_jobs[[length(stan_jobs) + 1L]] <- list(
    split = 1L, half = hf, frame = stan_frames[[1L]][[c("h1", "h2")[hf]]],
    elig = ceP$elig, chains = STAN_CHAINS, iter = STAN_ITER,
    iter_retry = STAN_ITER_RETRY, seed = seed_stan(1L, hf) + 5L,
    rhat_max = RHAT_MAX, ess_min = ESS_MIN,
    decov_shape = 1, decov_scale = 1)
}
t0 <- Sys.time()
stan_res <- parallel::parLapply(cl, stan_jobs, stan_job)
say("Stan fits done in %.1f min (%d fits incl. %d prior-sensitivity)\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), length(stan_res), 2L)

is_main <- vapply(stan_res, function(z) z$decov_scale == 2, logical(1))
main_res <- stan_res[is_main]
sens_res <- stan_res[!is_main]
EB1 <- matrix(NA_real_, ceP$n, N_SPLITS_BAYES)
EB2 <- matrix(NA_real_, ceP$n, N_SPLITS_BAYES)
for (z in main_res) {
  if (z$half == 1L) EB1[, z$split] <- z$est else EB2[, z$split] <- z$est
}
eP <- EST[[PRIMARY_TARGET]][[as.character(PRIMARY_B)]]
rs_raw_P <- split_correlations(eP$raw$h1, eP$raw$h2, MIN_PAIRS)
rs_glmer_P <- split_correlations(eP$hier$h1, eP$hier$h2, MIN_PAIRS)
rs_bayes <- split_correlations(EB1, EB2, MIN_PAIRS)                   # FROZEN
rel_raw_20 <- sb_correct(mean(rs_raw_P[seq_len(N_SPLITS_BAYES)], na.rm = TRUE))
rel_glmer_20 <- sb_correct(mean(rs_glmer_P[seq_len(N_SPLITS_BAYES)],
                                na.rm = TRUE))
rel_bayes_20 <- sb_correct(mean(rs_bayes, na.rm = TRUE))
delta_bayes <- rel_bayes_20 - rel_raw_20
delta_glmer_gate <- reliability_by_budget$delta[
  reliability_by_budget$target == PRIMARY_TARGET &
    reliability_by_budget$budget_B == PRIMARY_B]
delta_glmer_20 <- rel_glmer_20 - rel_raw_20
g3_diag_ok <- max(vapply(main_res, function(z) z$max_rhat, numeric(1))) <
  RHAT_MAX &&
  min(vapply(main_res, function(z) z$min_ess_bulk, numeric(1))) > ESS_MIN &&
  sum(vapply(main_res, function(z) z$divergences, numeric(1))) == 0
g3_pass <- g3_diag_ok &&
  sign(delta_bayes) == sign(delta_glmer_gate) &&
  abs(delta_bayes - delta_glmer_gate) < G3_ABS_TOL
# Prior sensitivity: Spearman between scale-2 and scale-1 estimates, split 1.
prior_sens <- vapply(1:2, function(hf) {
  a <- main_res[[which(vapply(main_res, function(z)
    z$split == 1L && z$half == hf, logical(1)))[1]]]$est
  b <- sens_res[[which(vapply(sens_res, function(z)
    z$half == hf, logical(1)))[1]]]$est
  ok <- is.finite(a) & is.finite(b)
  suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
}, numeric(1))
say("G3: delta_bayes %.4f vs delta_glmer(gate) %.4f (glmer same-20 %.4f); diag ok = %s -> %s\n",
    delta_bayes, delta_glmer_gate, delta_glmer_20, g3_diag_ok,
    if (g3_pass) "PASS" else "FAIL")
say("prior sensitivity (decov scale 2 vs 1, split 1): Spearman h1 %.5f h2 %.5f\n",
    prior_sens[1], prior_sens[2])

brms_confirmation <- dplyr::bind_rows(lapply(main_res, function(z)
  tibble::tibble(
    engine = "rstanarm::stan_glmer (precompiled; brms needs a C++ toolchain absent on this machine)",
    target = PRIMARY_TARGET, budget_B = PRIMARY_B,
    split = z$split, half = z$half, chains = STAN_CHAINS,
    iter = if (z$retried) STAN_ITER_RETRY else STAN_ITER,
    retried = z$retried, max_rhat = z$max_rhat,
    min_ess_bulk = z$min_ess_bulk, divergences = z$divergences))) |>
  dplyr::arrange(split, half) |>
  dplyr::mutate(
    n_splits_bayes = N_SPLITS_BAYES,
    rel_raw_same_splits = rel_raw_20,
    rel_hier_glmer_same_splits = rel_glmer_20,
    rel_hier_bayes = rel_bayes_20,
    delta_bayes = delta_bayes,
    delta_glmer_gate_100splits = delta_glmer_gate,
    delta_glmer_same_splits = delta_glmer_20,
    sign_match = sign(delta_bayes) == sign(delta_glmer_gate),
    abs_diff_vs_gate = abs(delta_bayes - delta_glmer_gate),
    g3_diagnostics_ok = g3_diag_ok, g3_pass = g3_pass,
    prior_sens_spearman_h1 = prior_sens[1],
    prior_sens_spearman_h2 = prior_sens[2])

parallel::stopCluster(cl)

# =============================================================================
# 13. Outputs
# =============================================================================
save_csv(reliability_by_budget, "reliability_by_budget")
save_csv(delta_bootstrap, "delta_bootstrap")
save_csv(null_permutation, "null_permutation")
save_csv(shrinkage_diagnostics, "shrinkage_diagnostics")
save_csv(ceiling_guard, "ceiling_guard")
save_csv(brms_confirmation, "brms_confirmation")
save_csv(per_subject_estimates_sample, "per_subject_estimates_sample")
save_csv(per_split_correlations, "per_split_correlations")

hg1_any <- any(ceiling_guard$ci_supported_exceedance)
hg2_any <- any(null_permutation$materially_nonzero)
hg4_conv_ok <- all(shrinkage_diagnostics$conv_fail_rate < CONV_FAIL_MAX)
hg4_collapse <- any(shrinkage_diagnostics$re_sd_collapsed_frac > 0.5,
                    na.rm = TRUE)

g1_gate <- list(
  spec = "elevation/E8_spec.md",
  smoke_run = SMOKE,
  primary = list(
    target = PRIMARY_TARGET, gating_budgets = G1_BUDGETS,
    delta_min = G1_DELTA_MIN, ci_level = 0.90,
    per_budget = lapply(seq_len(nrow(g1_tbl)), function(i) as.list(
      g1_tbl[i, c("budget_B", "delta", "delta_ci90_lo", "delta_ci90_hi",
                  "passes_g1_criteria", "in_g1_budgets")])),
    g1_verdict = if (g1_pass) "PASS" else "FAIL", g1_pass = g1_pass),
  secondary = list(
    target = "ev_sens_agn", informational_not_gating = TRUE,
    per_budget = lapply(seq_len(nrow(g2_tbl)), function(i) as.list(
      g2_tbl[i, c("budget_B", "delta", "delta_ci90_lo", "delta_ci90_hi",
                  "passes_g1_criteria", "in_g1_budgets")])),
    signal_present = g2_signal),
  g3_bayes = list(
    engine = "rstanarm::stan_glmer",
    engine_deviation = paste(
      "brms unavailable (no C++ toolchain on this machine); random-intercept",
      "SD prior is decov Gamma(1, 2) = Exponential(mean 2) instead of the",
      "spec's half-Normal(0, 2), reported with a prior-sensitivity check"),
    n_splits = N_SPLITS_BAYES, chains = STAN_CHAINS, iter = STAN_ITER,
    rel_raw_same_splits = rel_raw_20,
    rel_hier_bayes = rel_bayes_20, delta_bayes = delta_bayes,
    delta_glmer_gate = delta_glmer_gate,
    delta_glmer_same_splits = delta_glmer_20,
    sign_match = sign(delta_bayes) == sign(delta_glmer_gate),
    abs_diff_vs_gate = abs(delta_bayes - delta_glmer_gate),
    abs_tol = G3_ABS_TOL,
    max_rhat = max(vapply(main_res, function(z) z$max_rhat, numeric(1))),
    min_ess_bulk = min(vapply(main_res, function(z) z$min_ess_bulk,
                              numeric(1))),
    total_divergences = sum(vapply(main_res, function(z) z$divergences,
                                   numeric(1))),
    n_retried = sum(vapply(main_res, function(z) z$retried, logical(1))),
    diagnostics_ok = g3_diag_ok,
    prior_sensitivity_spearman = as.list(stats::setNames(prior_sens,
                                                         c("h1", "h2"))),
    g3_verdict = if (g3_pass) "PASS" else "FAIL", g3_pass = g3_pass),
  honesty_guards = list(
    HG1_ceiling = list(
      any_ci_supported_exceedance = hg1_any,
      flag_meaning = "leakage flag, never a gain",
      per_cell = lapply(seq_len(nrow(ceiling_guard)), function(i)
        as.list(ceiling_guard[i, ]))),
    HG2_null_permutation = list(
      any_materially_nonzero = hg2_any,
      criterion = sprintf("|mean perm r| > %.2f AND |z| > %g", HG2_ABS_MAX,
                          HG2_Z_MAX)),
    HG3_no_cross_half_leakage = list(
      subject_problem_overlaps_across_all_split_fits = sum(HG3_OVERLAPS),
      separate_model_calls_per_half = TRUE,
      hyperparameters_within_half = TRUE),
    HG4_shrinkage = list(
      all_conv_fail_rates_below_target = hg4_conv_ok,
      conv_fail_target = CONV_FAIL_MAX,
      any_cell_re_sd_collapsed = hg4_collapse)),
  g_eq = list(pairs_rederived = geq_n, mismatches = geq_bad,
              bit_identical = geq_bad == 0L),
  interpretation_predeclared = list(
    if_g1_pass = paste(
      "hierarchical re-estimation recovers reducible measurement noise; a",
      "constructive partial remedy; the population true-variance ceiling is",
      "untouched"),
    if_g1_fail = paste(
      "the frozen per-subject indices are already near-efficient at CPC18",
      "budgets; the wall is a true-variance ceiling, not an estimation",
      "artifact; equally publishable and arguably sharper for Study 1")),
  build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
  optimization_performed = FALSE, mechanism_model_fitted = FALSE,
  designs_constructed = FALSE)
jsonlite::write_json(g1_gate, file.path(res_dir, "g1_gate.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

# Seeds (vectors built OUTSIDE tibble() and length-checked, guarding against
# unintended recycling).
seed_rows <- list()
for (B in B_GRID) {
  seed_rows[[length(seed_rows) + 1L]] <- tibble::tibble(
    purpose = "split_partition_draws", target = NA_character_, budget_B = B,
    seed = seed_draws(B), derivation = "SEED_E8 + 1000*B")
  for (tg in TARGETS) {
    t_idx <- match(tg, TARGETS)
    sv <- c(seed_boot(t_idx, B), seed_perm(t_idx, B),
            seed_relpt(t_idx, B, 1L), seed_relpt(t_idx, B, 2L))
    stopifnot(length(sv) == 4L)
    seed_rows[[length(seed_rows) + 1L]] <- tibble::tibble(
      purpose = c("paired_delta_bootstrap", "hg2_permutation",
                  "reliability_point_raw", "reliability_point_hier"),
      target = tg, budget_B = B, seed = sv,
      derivation = c(
        "SEED_E8 + 3000000 + 100000*t_idx + 1000*B",
        "SEED_E8 + 6000000 + 100000*t_idx + 1000*B",
        "SEED_E8 + 9000000 + 0 + 100000*t_idx + 1000*B",
        "SEED_E8 + 9000000 + 1000000 + 100000*t_idx + 1000*B"))
  }
}
seed_rows[[length(seed_rows) + 1L]] <- tibble::tibble(
  purpose = "stan_fit_base", target = PRIMARY_TARGET, budget_B = PRIMARY_B,
  seed = seed_stan(0L, 0L),
  derivation = paste("SEED_E8 + 12000000 + 10*split + half; prior-sensitivity",
                     "fits use + 5"))
e8_seeds <- dplyr::bind_rows(seed_rows)
save_csv(e8_seeds, "e8_seeds")

cfg <- list(
  spec = "elevation/E8_spec.md", smoke_run = SMOKE,
  raw_md5 = raw_md5, seed_base = SEED_E8,
  regime = "experience (Trial 6-25)", r_reps = R_REPS,
  budgets_problems_per_half = B_GRID,
  targets = list(primary = PRIMARY_TARGET, secondary = "ev_sens_agn"),
  n_splits_hier = N_SPLITS, min_pairs = MIN_PAIRS, n_boot = N_BOOT,
  n_sample_splits_committed = N_SAMPLE_SPLITS,
  glmer = list(
    package_version = as.character(utils::packageVersion("lme4")),
    optimizer = "bobyqa", maxfun_first = GLMER_MAXFUN1,
    maxfun_refit = GLMER_MAXFUN2,
    convergence_rule = paste(
      "non-convergence (optimizer code != 0 or a 'failed to converge'",
      "message) triggers ONE refit with the longer budget; a fit that still",
      "fails contributes NA for every subject-half it covers, excluded from",
      "complete pairs and counted in conv_fail_rate"),
    response = "per subject-by-problem binomial proportion, weights = 20"),
  hier_models = list(
    loss_avoid_adj = "y ~ 1 + d_ev_noloss + (1 | SubjID); estimate = conditional mode of the subject random intercept",
    ev_sens_agn = "y ~ 1 + z_ev + (1 + z_ev | SubjID); estimate = conditional mode of the subject random slope on z_ev"),
  pairing = paste(
    "one committed draw stream per budget (seed_draws) generates every",
    "subject's h1/h2 rows once per split, mirroring the frozen run_cell",
    "lambda = NA branch (pos = full experience pool, consumes no RNG); BOTH",
    "arms consume identical indices, so per-split differences isolate the",
    "estimator"),
  raw_arm = paste(
    "the FROZEN estimate_half called directly per subject-half (no",
    "accelerator); G_EQ re-derives splits 1-3 of every cell and asserts",
    "bit-identity"),
  reliability_estimand = paste(
    "frozen split_correlations (Spearman, complete pairs, min_pairs 10) ->",
    "mean over splits -> frozen sb_correct; per-arm 95% CI from the frozen",
    "reliability_point; Delta CI from the paired subject-cluster bootstrap",
    "(one index vector per draw pushed through BOTH arms; bootstrap outside,",
    "splits inside)"),
  e0_anchors = list(
    source = "elevation/results/e0/variance_ceiling_final.csv (read-only)",
    r1_loss_avoid_adj = unname(R1["loss_avoid_adj"]),
    r1_ev_sens_agn = unname(R1["ev_sens_agn"]),
    derivation = paste(
      "r1 = R/(k - (k-1)R) at the LOCKED k_ref = 14 full-data ceiling;",
      "committed k_0.7/k_0.8 columns reproduce under the same r1 to < 1e-9;",
      "spec-quoted 0.0787 / 0.02717 are these values rounded"),
    rel_ceiling_2B = "2B * r1 / (1 + (2B - 1) * r1)"),
  hg2_permutation = paste(
    "fresh subject permutation of the h2 estimate vector per split (seeded,",
    "committed), frozen split_correlations on the broken pairing; flag if",
    sprintf("|mean r| > %.2f AND |z| > %g", HG2_ABS_MAX, HG2_Z_MAX)),
  hg4_note = paste(
    "raw and hier estimates live on different scales (probability-residual",
    "vs logit intercept; Spearman index vs logistic slope), so",
    "sd_ratio_hier_over_raw is a cross-scale collapse diagnostic, not a",
    "calibrated shrinkage estimate; spearman_raw_hier is scale-free"),
  bayes = list(
    engine = "rstanarm::stan_glmer",
    rstanarm_version = as.character(utils::packageVersion("rstanarm")),
    rstan_version = as.character(utils::packageVersion("rstan")),
    chains = STAN_CHAINS, iter = STAN_ITER, iter_retry = STAN_ITER_RETRY,
    priors = paste(
      "fixed effects Normal(0, 5) [spec]; random-intercept SD decov",
      "Gamma(shape 1, scale 2) = Exponential(mean 2) SUBSTITUTED for the",
      "spec's half-Normal(0, 2): brms cannot run without a C++ toolchain",
      "(no RTools on this machine) and rstanarm's decov cannot express",
      "half-Normal; prior-sensitivity spot check (scale 1 vs 2) committed in",
      "brms_confirmation.csv"),
    diagnostics_rule = paste(
      "per fit: max Rhat and min bulk-ESS over the subject intercepts b[...],",
      "both fixed effects and Sigma (lp__ and mean_PPD excluded), plus",
      "divergence count; a fit failing Rhat < 1.01 / ESS > 400 / 0",
      "divergences is refit once at iter_retry")),
  environment = list(
    r_version = R.version.string,
    lme4 = as.character(utils::packageVersion("lme4")),
    rstanarm = as.character(utils::packageVersion("rstanarm")),
    matrix = as.character(utils::packageVersion("Matrix")),
    n_workers = N_WORKERS,
    renv_snapshot_note = paste(
      "lme4 + rstanarm installed via renv::install into the project library;",
      "renv.lock was not updated in the same run;",
      "both packages were added to it",
      "afterwards with renv::snapshot()")),
  float_precision_note = paste(
    "CSVs written by readr::write_csv (shortest round-trippable",
    "representation); for bit-identity checks parse with utils::read.csv,",
    "not readr (the E3b vroom float-parse finding)"),
  data_inputs = c("data/raw/All estimation raw data.csv",
                  "elevation/results/e1/participant_split.csv",
                  "elevation/results/e0/variance_ceiling_final.csv"),
  split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
  loss_baselines_exp = as.list(baselines$exp),
  eligible_n = stats::setNames(
    lapply(cells, function(ce) ce$n), paste0("B", B_GRID)),
  no_freeze_or_test_choice_read = TRUE,
  optimization_performed = FALSE, mechanism_model_fitted = FALSE,
  runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
jsonlite::write_json(cfg, file.path(res_dir, "e8_config.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

say("\nE8 complete in %.1f min. Outputs in %s\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")), res_dir)
say("G1 (%s, Delta >= %.2f with 90%% CI excluding 0 at B in {8,10,12}): %s\n",
    PRIMARY_TARGET, G1_DELTA_MIN, if (g1_pass) "PASS" else "FAIL")
say("G3 (Bayesian confirmation): %s | HG1 exceedance: %s | HG2 nonzero: %s | HG4 conv ok: %s\n",
    if (g3_pass) "PASS" else "FAIL", hg1_any, hg2_any, hg4_conv_ok)

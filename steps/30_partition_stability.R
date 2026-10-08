# Step 30: stability of the subset ranking across problem partitions
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e0/, e1/, e2/, e3/ (steps 19-22)
# Writes:  output/elevation/results/dispb/ (b1_*.csv, dispb_part_b1.json, seeds)
# Paper:   Study 2 (Ruling out alternatives); Supplement (Study 2 robustness detail), through steps 33 and 36
# Runtime: about 2 hours with 12 parallel workers

# B1: test the stability of H1 with pre-declared gates G-B1a..d.
# The original E3 result used one gate cell, one Monte-Carlo set sample, and
# one noise-floor estimand. This check uses R_SEED independent set samples,
# two noise estimands, random subject partitions, and a null-permutation control.
# It re-scores the existing subsets without optimizing designs or fitting a
# new choice model. The original H1 PASS is retained, and any failed robustness
# gate is reported.
#
# Reliabilities use the 344 build subjects from step 20; Set-1 cells use its
# 63 build subjects. Choice/outcome columns are dropped from the full-sample
# frame at load, excluding the 137 freeze and 205 test subjects' choices.
#
# The shared functions are sourced without redefinition, checked by
# 33_partition_checks_verify.R. Scoring uses split_correlations and sb_correct;
# noise floors use reliability_point. The two target indices use the
# accelerator documented in steps 22 and 24. Step 33's G_EQ checks it against
# estimate_half to < 1e-12 with identical NA patterns, and a whole-set spot
# check verifies r_sb through the original engine.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})

source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))

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
# 0. Pre-declared constants (fixed BEFORE any result)
# =============================================================================
RAW_MD5      <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_B1      <- 20260722L
SEED_E3      <- 20260719L          # only to REPLAY E3's committed anchor sample
REGIME       <- "exp"
POS          <- 6:25               # experience regime, r = 20 (the full pool)
R_REPS       <- length(POS)
K_GRID       <- c(6L, 8L, 10L, 12L)
TARGETS      <- c("loss_avoid_adj", "ev_sens_agn")   # primary, secondary
M_CAP        <- 2000L
R_SEED       <- 8L                 # independent Monte-Carlo set samples
N_SPLIT_B1   <- 100L
N_NOISE_SETS <- 50L
N_BOOT       <- 1000L
B_PART_B1    <- 200L               # subject partitions at the primary cell
MIN_PAIRS    <- 10L                # frozen split_correlations default
GATE_FORM    <- 1L
GATE_K       <- 12L
GATE_TARGET  <- "loss_avoid_adj"
H1_SIGMA_MIN <- 0.05               # frozen E3 threshold; nothing re-tuned
H1_ICC_MIN   <- 0.50               # frozen E3 threshold; nothing re-tuned
N_HALF_A     <- 32L                # largest-remainder halves of 63
N_WORKERS    <- 12L

# Further constants, fixed HERE before any result exists and recorded in the
# part file. The design calls for the covariance estimand "for each K" and for
# a >= 200-partition distribution. Running 200 partitions at every (seed, K) is
# ~6400 partition-cells and is not affordable in one run, so the partition
# budget is allocated as follows, fixed in advance and reported in full:
#   - PRIMARY cell (K = 12): the full B_PART_B1 = 200 partitions. Partition p is
#     evaluated on set sample si = ((p - 1) mod R_SEED) + 1, so the 200-value
#     distribution marginalises over BOTH analyst choices (set sample and
#     subject partition) and every seed contributes exactly 25 partitions.
#   - Auxiliary K (6, 8, 10): N_PART_AUX partitions per (seed, K).
#   - Full map (G-B1d, reported not gated): N_PART_MAP partitions per cell.
N_PART_AUX <- 8L
N_PART_MAP <- 8L
# Sets whose r_sb is finite on BOTH subject halves enter the covariance
# estimand; a cell needs at least this many to be reported.
MIN_SETS_COV <- 30L

# ---- seed derivations (disjoint blocks, asserted) ---------------------------
# Pre-declared:
#   set-draw     = SEED_B1 + 1000*form + K + 10000*seed_index
#   split-draw   = set-draw + 500000            (per set: + set_id)
#   partition    = SEED_B1 + 700000 + partition_index
# Fixed in code (blocks chosen disjoint and asserted):
seed_set_b1   <- function(form, K, si) SEED_B1 + 1000L * form + K + 10000L * si
seed_split_b1 <- function(form, K, si) seed_set_b1(form, K, si) + 500000L
seed_noise_b1 <- function(form, K, si) SEED_B1 + 600000L + 1000L * form + K +
  10000L * si
seed_part_b1  <- function(p) SEED_B1 + 700000L + p
seed_part_aux <- function(K, si, j) SEED_B1 + 720000L + 1000L * K + 10L * si + j
seed_part_map <- function(form, K, j) SEED_B1 + 750000L + 1000L * form +
  10L * K + j
seed_perm_b1  <- function(p) SEED_B1 + 800000L + p
# Per-half halving seed. `part` is a partition INDEX unique within its family:
# 1..200 for the primary cell, 201..288 for auxiliary K, 1101..1808 for the map
# (1000 + 100*form + j). h is 1 (half A) or 2 (half B).
seed_half_b1 <- function(h, K, part, sid)
  SEED_B1 + 900000000L * (h - 1L) + 200000000L * match(K, K_GRID) +
  100000L * part + sid
# E3's committed derivations, replayed ONLY to reproduce its anchor sample.
seed_subset_e3 <- function(form, K) SEED_E3 + 1000L * form + K
seed_split_e3  <- function(form, K) SEED_E3 + 1000L * form + K + 500000L
seed_noise_e3  <- function(form, K) SEED_E3 + 600000L + 1000L * form + K

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
# 2. Per-form matrices (E3 construction, unchanged)
# =============================================================================
form_of <- build_trials |> dplyr::distinct(SubjID, Set)
form_problems <- lapply(1:7, function(S)
  sort(unique(build_trials$GameID[build_trials$Set == S])))
form_subjects <- lapply(1:7, function(S)
  as.character(sort(form_of$SubjID[form_of$Set == S])))
make_form_mats <- function(S) {
  probs <- form_problems[[S]]; subs <- form_subjects[[S]]
  NB <- matrix(NA_real_, length(subs), length(probs),
               dimnames = list(subs, as.character(probs)))
  for (sid in subs) {
    sd <- subj_data[[sid]]
    j <- match(sd$feat$GameID, probs)
    ok <- !is.na(j)
    NB[sid, j[ok]] <- rowSums(sd$B_mat[which(ok), POS, drop = FALSE])
  }
  fi <- feat_all[match(probs, feat_all$GameID), ]
  list(probs = probs, subs = subs, NB = NB, feat = fi,
       li = fi$loss_informative, has_loss_b = fi$has_loss_b,
       d_ev = fi$d_ev, p_base = fi$p_base_exp)
}
form_mats <- lapply(1:7, make_form_mats)
stopifnot(length(form_mats[[GATE_FORM]]$subs) == 63L)
stopifnot(!anyNA(form_mats[[GATE_FORM]]$NB))
say("Set 1: %d build subjects, %d clean problems, %d loss-informative\n",
    length(form_mats[[1]]$subs), length(form_mats[[1]]$probs),
    sum(form_mats[[1]]$li))

# =============================================================================
# 3. Accelerator (identical formulas to the frozen estimate_half; see header)
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

half_targets <- function(NBh, fm, cols, r) {
  k <- length(cols)
  n <- r * k
  R  <- ev_ranks(fm$d_ev[cols], r)
  sR <- sum(R); sR2 <- sum(R * R)
  n1  <- rowSums(NBh)
  Sxy <- as.vector(NBh %*% R) - n1 * sR / k
  Sxx <- n1 - n1 * n1 / n
  Syy <- r * (sR2 - sR * sR / k)
  ev <- rep(NA_real_, nrow(NBh))
  ok <- Sxx > 0 & Syy > 0
  if (any(ok)) ev[ok] <- Sxy[ok] / sqrt(Sxx[ok] * Syy)

  lic <- cols[fm$li[cols]]
  if (length(lic) == 0) {
    loss <- rep(NA_real_, nrow(NBh))
  } else {
    NBl <- NBh[, fm$li[cols], drop = FALSE]
    cnt <- matrix(0, nrow(NBl), ncol(NBl))
    for (j in seq_along(lic)) {
      cnt[, j] <- if (fm$has_loss_b[lic[j]]) r - NBl[, j] else NBl[, j]
    }
    loss <- (rowSums(cnt) - r * sum(fm$p_base[lic])) / (r * length(lic))
  }
  list(ev_sens_agn = ev, loss_avoid_adj = loss)
}

# rows = subject row indices into fm$NB. Subjects who did not see EVERY problem
# in S are dropped (the E3 eligibility rule, which reduces to "the whole form"
# wherever the common-set premise holds). The halving of S is COMMON to all
# subjects within a split, as in E3.
run_set_rows <- function(fm, cols, rows, n_split, seed, want_est = FALSE) {
  NBs <- fm$NB[rows, cols, drop = FALSE]
  keep <- which(!is.na(rowSums(NBs)))
  NBs <- NBs[keep, , drop = FALSE]
  ns <- length(keep); K <- length(cols); kh <- K %/% 2L
  E <- list()
  for (tg in TARGETS) E[[tg]] <- list(h1 = matrix(NA_real_, ns, n_split),
                                      h2 = matrix(NA_real_, ns, n_split))
  set.seed(seed)
  for (s in seq_len(n_split)) {
    take <- sample.int(K, 2L * kh)
    i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
    a <- half_targets(NBs[, i1, drop = FALSE], fm, cols[i1], R_REPS)
    b <- half_targets(NBs[, i2, drop = FALSE], fm, cols[i2], R_REPS)
    for (tg in TARGETS) { E[[tg]]$h1[, s] <- a[[tg]]; E[[tg]]$h2[, s] <- b[[tg]] }
  }
  out <- list(n_subjects = ns)
  for (tg in TARGETS) {
    rs <- split_correlations(E[[tg]]$h1, E[[tg]]$h2, MIN_PAIRS)        # FROZEN
    r_raw <- mean(rs, na.rm = TRUE)
    out[[tg]] <- list(
      r_raw = r_raw,
      r_sb = if (is.finite(r_raw)) sb_correct(r_raw) else NA_real_,    # FROZEN
      coverage = mean(!is.na(rs)),
      n_splits_used = sum(!is.na(rs)))
  }
  if (want_est) out$E <- E
  out
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

# =============================================================================
# 4. E3's committed anchor sample and the E0 ceiling anchors
# =============================================================================
e3_real <- readr::read_csv(
  here("output", "elevation", "results", "e3", "reliability_by_set_real.csv"),
  show_col_types = FALSE)
e3_s1 <- e3_real |> dplyr::filter(form == GATE_FORM)
stopifnot(nrow(e3_s1) == 2L * length(K_GRID) * M_CAP)

parse_ids <- function(s) as.integer(strsplit(s, ";", fixed = TRUE)[[1]])
anchor_sets <- lapply(K_GRID, function(K) {
  d <- e3_s1 |> dplyr::filter(target == TARGETS[1], K == !!K) |>
    dplyr::arrange(set_id)
  stopifnot(identical(as.integer(d$set_id), seq_len(M_CAP)))
  lapply(d$game_ids, function(s) {
    cc <- match(parse_ids(s), form_mats[[GATE_FORM]]$probs)
    stopifnot(!anyNA(cc), !anyDuplicated(cc))
    sort(cc)
  })
})
names(anchor_sets) <- as.character(K_GRID)
# Independent re-derivation from E3's committed SEEDS alone (the check used in
# 24_denoised_subset_reliability.R).
anchor_seed_match <- vapply(K_GRID, function(K) {
  rs <- draw_sets(length(form_mats[[GATE_FORM]]$probs), K, M_CAP,
                  seed_subset_e3(GATE_FORM, K))
  identical(lapply(rs[seq_len(M_CAP)], as.integer),
            lapply(anchor_sets[[as.character(K)]], as.integer))
}, logical(1))
stopifnot(all(anchor_seed_match))
say("anchor sample: E3's committed Set-1 sets re-derived from E3 seeds alone MATCH at all %d K\n",
    length(K_GRID))

e0 <- readr::read_csv(
  here("output", "elevation", "results", "e0", "variance_ceiling_final.csv"),
  show_col_types = FALSE)
r1_of <- function(R, k) R / (k - (k - 1) * R)
sb_k  <- function(k, r1) k * r1 / (1 + (k - 1) * r1)
R1 <- vapply(TARGETS, function(tg) {
  rw <- e0[e0$construct == tg, ]
  stopifnot(nrow(rw) == 1L, rw$k_ref == 14)
  r1 <- r1_of(rw$ceiling_fulldata, 14)
  k_needed <- function(R) R * (1 - r1) / (r1 * (1 - R))
  stopifnot(abs(k_needed(0.7) - rw$k_0.7) < 1e-9,
            abs(k_needed(0.8) - rw$k_0.8) < 1e-9)
  r1
}, numeric(1))
say("E0 ceiling anchors: r1 loss %.6f, ev %.6f (k_0.7 / k_0.8 reproduce exactly)\n",
    R1["loss_avoid_adj"], R1["ev_sens_agn"])
# HG-ceiling comparator, fixed here before any result: the E0 full-data ceiling
# SB-projected in PROBLEM COUNT to the cell's K, and E8's exceedance rule (a
# flag requires the frozen bootstrap CI LOWER bound to clear the ceiling, i.e.
# an exceedance "beyond bootstrap noise"). Point exceedance is reported too, as
# a descriptive fraction, because the maximum over 2000 sampled sets is a
# selection statistic and is expected to sit above a population average.
# The comparator is sb_k(K, R1[target]); it is applied in section 8.

# =============================================================================
# 5. Estimand (i): bootstrap-floor decomposition over R_SEED set samples
# =============================================================================
# seed_index 0 is E3's committed sample replayed with E3's own split and noise
# seeds. It is an ANCHOR: it must reproduce the frozen h1_gate numbers, and it
# is EXCLUDED from the G-B1a battery minimum (the battery is the R_SEED = 8
# pre-declared fresh samples).
say("\n=== B1 estimand (i): bootstrap-floor decomposition ===\n")
Z975 <- stats::qnorm(0.975)
si_grid <- 0:R_SEED

cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(tidyverse); library(here) })
  source(here("R", "step2_reliability.R"))
  NULL
}))
parallel::clusterExport(cl, c("form_mats", "half_targets", "ev_ranks",
                              "run_set_rows", "TARGETS", "R_REPS", "MIN_PAIRS",
                              "N_SPLIT_B1", "N_NOISE_SETS", "N_BOOT", "Z975",
                              "M_CAP", "GATE_FORM"),
                        envir = environment())

jobs_i <- list()
for (si in si_grid) for (K in K_GRID) {
  sets <- if (si == 0L) anchor_sets[[as.character(K)]] else
    draw_sets(length(form_mats[[GATE_FORM]]$probs), K, M_CAP,
              seed_set_b1(GATE_FORM, K, si))
  stopifnot(length(sets) == M_CAP)
  stopifnot(!anyDuplicated(vapply(sets, paste, character(1), collapse = ",")))
  jobs_i[[length(jobs_i) + 1L]] <- list(
    si = si, K = K, sets = sets,
    split_base = if (si == 0L) seed_split_e3(GATE_FORM, K) else
      seed_split_b1(GATE_FORM, K, si),
    noise_base = if (si == 0L) seed_noise_e3(GATE_FORM, K) else
      seed_noise_b1(GATE_FORM, K, si))
}
t0 <- Sys.time()
res_i <- parallel::parLapply(cl, jobs_i, function(job) {
  fm <- form_mats[[GATE_FORM]]
  rows <- seq_along(fm$subs)
  M <- length(job$sets)
  per <- vector("list", M)
  noise <- vector("list", min(N_NOISE_SETS, M))
  for (i in seq_len(M)) {
    want <- i <= N_NOISE_SETS
    r <- run_set_rows(fm, job$sets[[i]], rows, N_SPLIT_B1,
                      job$split_base + i, want_est = want)
    per[[i]] <- tibble::tibble(
      set_id = i,
      game_ids = paste(fm$probs[job$sets[[i]]], collapse = ";"),
      n_loss_informative = sum(fm$li[job$sets[[i]]]),
      n_subjects = r$n_subjects,
      r_sb_loss_avoid_adj = r$loss_avoid_adj$r_sb,
      r_sb_ev_sens_agn = r$ev_sens_agn$r_sb,
      cov_loss_avoid_adj = r$loss_avoid_adj$coverage,
      cov_ev_sens_agn = r$ev_sens_agn$coverage)
    if (want) {
      noise[[i]] <- purrr::map_dfr(TARGETS, function(tg) {
        rp <- reliability_point(r$E[[tg]]$h1, r$E[[tg]]$h2, N_BOOT,
                                job$noise_base + i, MIN_PAIRS)          # FROZEN
        tibble::tibble(set_id = i, target = tg, r_sb = rp$r_sb,
                       ci_lo_sb = rp$ci_lo_sb, ci_hi_sb = rp$ci_hi_sb,
                       sd_noise = (rp$ci_hi_sb - rp$ci_lo_sb) / (2 * Z975))
      })
    }
  }
  list(si = job$si, K = job$K, per = dplyr::bind_rows(per),
       noise = dplyr::bind_rows(noise))
})
say("estimand (i): %d cells in %.1f min\n", length(res_i),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

per_set_i <- purrr::map_dfr(res_i, function(z)
  z$per |> dplyr::mutate(form = GATE_FORM, seed_index = z$si, K = z$K,
                         .before = 1))
noise_i <- purrr::map_dfr(res_i, function(z)
  z$noise |> dplyr::mutate(form = GATE_FORM, seed_index = z$si, K = z$K,
                           .before = 1))

# ---- ANCHOR check: seed_index 0 must reproduce E3's committed per-set r_sb --
anchor_chk <- per_set_i |>
  dplyr::filter(seed_index == 0L) |>
  dplyr::select(K, set_id, r_sb_loss_avoid_adj, r_sb_ev_sens_agn) |>
  dplyr::left_join(
    e3_s1 |> dplyr::select(K, set_id, target, r_sb) |>
      tidyr::pivot_wider(names_from = target, values_from = r_sb,
                         names_prefix = "e3_"),
    by = c("K", "set_id"))
anchor_d_loss <- max(abs(anchor_chk$r_sb_loss_avoid_adj -
                           anchor_chk$e3_loss_avoid_adj), na.rm = TRUE)
anchor_d_ev <- max(abs(anchor_chk$r_sb_ev_sens_agn -
                         anchor_chk$e3_ev_sens_agn), na.rm = TRUE)
anchor_na_ok <- identical(is.na(anchor_chk$r_sb_loss_avoid_adj),
                          is.na(anchor_chk$e3_loss_avoid_adj))
say("anchor vs committed E3 per-set r_sb: max|d| loss = %.3g (NA pattern identical = %s), ev = %.3g\n",
    anchor_d_loss, anchor_na_ok, anchor_d_ev)
stopifnot(anchor_na_ok, anchor_d_loss < 1e-12)

q <- function(x, p) if (all(is.na(x))) NA_real_ else
  unname(stats::quantile(x, p, na.rm = TRUE))
long_i <- per_set_i |>
  tidyr::pivot_longer(c(r_sb_loss_avoid_adj, r_sb_ev_sens_agn),
                      names_to = "target", values_to = "r_sb") |>
  dplyr::mutate(target = sub("^r_sb_", "", target)) |>
  dplyr::left_join(
    per_set_i |>
      tidyr::pivot_longer(c(cov_loss_avoid_adj, cov_ev_sens_agn),
                          names_to = "target", values_to = "split_coverage") |>
      dplyr::mutate(target = sub("^cov_", "", target)) |>
      dplyr::select(seed_index, K, set_id, target, split_coverage),
    by = c("seed_index", "K", "set_id", "target"))

est_i <- long_i |>
  dplyr::group_by(form, seed_index, K, target) |>
  dplyr::summarise(
    M = dplyr::n(),
    n_sets_identified = sum(is.finite(r_sb)),
    set_coverage = mean(is.finite(r_sb)),
    split_coverage = mean(split_coverage),
    n_subjects_median = stats::median(n_subjects),
    sigma2_set_raw = stats::var(r_sb, na.rm = TRUE),
    r_sb_min = min(r_sb, na.rm = TRUE), r_sb_median = q(r_sb, .50),
    r_sb_max = max(r_sb, na.rm = TRUE),
    best_minus_median = max(r_sb, na.rm = TRUE) - q(r_sb, .50),
    .groups = "drop") |>
  dplyr::left_join(
    noise_i |> dplyr::group_by(seed_index, K, target) |>
      dplyr::summarise(sigma2_noise = mean(sd_noise^2, na.rm = TRUE),
                       n_noise_sets = sum(is.finite(sd_noise)), .groups = "drop"),
    by = c("seed_index", "K", "target")) |>
  dplyr::mutate(
    estimand = "bootstrap_floor",
    sigma_set_raw = sqrt(sigma2_set_raw),
    sigma_noise = sqrt(sigma2_noise),
    sigma2_set_true = pmax(0, sigma2_set_raw - sigma2_noise),
    sigma_set_true = sqrt(sigma2_set_true),
    ICC_set = sigma2_set_true / sigma2_set_raw,
    n_partitions = NA_integer_, is_anchor = seed_index == 0L)

gate_i <- est_i |> dplyr::filter(K == GATE_K, target == GATE_TARGET)
say("\nestimand (i) at the primary cell (Set %d, K %d, %s):\n",
    GATE_FORM, GATE_K, GATE_TARGET)
print(as.data.frame(gate_i |> dplyr::select(seed_index, sigma_set_raw,
                                            sigma_noise, sigma_set_true,
                                            ICC_set, set_coverage)))
anchor_row <- gate_i |> dplyr::filter(seed_index == 0L)
say("anchor sigma_set_true %.15f vs frozen h1_gate 0.187659511205571 (|d| = %.3g)\n",
    anchor_row$sigma_set_true, abs(anchor_row$sigma_set_true - 0.187659511205571))
say("anchor ICC_set      %.15f vs frozen h1_gate 0.625825410407600 (|d| = %.3g)\n",
    anchor_row$ICC_set, abs(anchor_row$ICC_set - 0.6258254104076))

# =============================================================================
# 6. Estimand (ii): covariance-based (subject-split) decomposition
# =============================================================================
say("\n=== B1 estimand (ii): covariance (subject-split) decomposition ===\n")
# One partition-cell: split the eligible subject rows into two disjoint halves,
# estimate each set's r_sb independently on each half (independent halvings),
# and return the per-set pair. Halves are a property of SUBJECTS, fixed across
# all sets within the partition (the rule used in step 24).
partition_rows <- function(n, seed, n_a) {
  set.seed(seed)
  p <- sample.int(n)
  list(A = sort(p[seq_len(n_a)]), B = sort(p[(n_a + 1L):n]))
}

cov_cell <- function(fm, sets, rows_all, seed_part, part_index, K, n_split) {
  n <- length(rows_all)
  n_a <- ceiling(n / 2)
  pr <- partition_rows(n, seed_part, n_a)
  rA <- rows_all[pr$A]; rB <- rows_all[pr$B]
  stopifnot(length(intersect(rA, rB)) == 0L)
  M <- length(sets)
  out <- matrix(NA_real_, M, 4L)
  nsub <- matrix(NA_integer_, M, 2L)
  for (i in seq_len(M)) {
    a <- run_set_rows(fm, sets[[i]], rA, n_split,
                      seed_half_b1(1L, K, part_index, i))
    b <- run_set_rows(fm, sets[[i]], rB, n_split,
                      seed_half_b1(2L, K, part_index, i))
    out[i, 1L] <- a$loss_avoid_adj$r_sb; out[i, 2L] <- b$loss_avoid_adj$r_sb
    out[i, 3L] <- a$ev_sens_agn$r_sb;    out[i, 4L] <- b$ev_sens_agn$r_sb
    nsub[i, ] <- c(a$n_subjects, b$n_subjects)
  }
  list(r = out, n_a = length(rA), n_b = length(rB), n_sub = nsub)
}

# Decomposition + the pre-declared null control (permute the set-to-set
# correspondence between the two subject halves; the covariance must collapse).
cov_stats <- function(xA, xB, seed_perm, n_perm = 20L) {
  ok <- is.finite(xA) & is.finite(xB)
  n_used <- sum(ok)
  if (n_used < MIN_SETS_COV)
    return(list(n_used = n_used, cov = NA_real_, var_a = NA_real_,
                var_b = NA_real_, icc = NA_real_, sigma_true = NA_real_,
                sigma_raw = NA_real_, rho = NA_real_,
                cov_null = NA_real_, sigma_true_null = NA_real_))
  a <- xA[ok]; b <- xB[ok]
  cv <- stats::cov(a, b)
  va <- stats::var(a); vb <- stats::var(b)
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

# ---- primary cell: the full B_PART_B1 partition distribution ---------------
prim_jobs <- lapply(seq_len(B_PART_B1), function(p) {
  si <- ((p - 1L) %% R_SEED) + 1L
  list(p = p, si = si, K = GATE_K,
       sets = draw_sets(length(form_mats[[GATE_FORM]]$probs), GATE_K, M_CAP,
                        seed_set_b1(GATE_FORM, GATE_K, si)),
       seed_part = seed_part_b1(p), part_index = p,
       seed_perm = seed_perm_b1(p))
})
parallel::clusterExport(cl, c("seed_half_b1", "partition_rows", "cov_cell",
                              "cov_stats", "MIN_SETS_COV", "K_GRID", "SEED_B1"),
                        envir = environment())
t0 <- Sys.time()
prim_res <- parallel::parLapply(cl, prim_jobs, function(job) {
  fm <- form_mats[[GATE_FORM]]
  cc <- cov_cell(fm, job$sets, seq_along(fm$subs), job$seed_part,
                 job$part_index, job$K, N_SPLIT_B1)
  sl <- cov_stats(cc$r[, 1L], cc$r[, 2L], job$seed_perm)
  se <- cov_stats(cc$r[, 3L], cc$r[, 4L], job$seed_perm + 500000L)
  list(p = job$p, si = job$si, n_a = cc$n_a, n_b = cc$n_b,
       loss = sl, ev = se,
       export = tibble::tibble(partition = job$p, seed_index = job$si,
                               K = job$K, set_id = seq_len(nrow(cc$r)),
                               r_sb_A_loss_avoid_adj = cc$r[, 1L],
                               r_sb_B_loss_avoid_adj = cc$r[, 2L],
                               r_sb_A_ev_sens_agn = cc$r[, 3L],
                               r_sb_B_ev_sens_agn = cc$r[, 4L]))
})
say("primary-cell covariance: %d partitions in %.1f min\n", B_PART_B1,
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

mk_part_row <- function(z, tg) {
  s <- if (tg == "loss_avoid_adj") z$loss else z$ev
  tibble::tibble(
    form = GATE_FORM, K = GATE_K, target = tg, partition = z$p,
    seed_index = z$si, n_subjects_A = z$n_a, n_subjects_B = z$n_b,
    n_sets_used = s$n_used, cov_set = s$cov, var_A = s$var_a, var_B = s$var_b,
    sigma_set_raw_cov = s$sigma_raw, sigma_set_true_cov = s$sigma_true,
    ICC_cov = s$icc, spearman_A_vs_B = s$rho,
    cov_set_null = s$cov_null, sigma_set_true_cov_null = s$sigma_true_null)
}
b1_partition_distribution <- dplyr::bind_rows(
  purrr::map_dfr(prim_res, mk_part_row, tg = "loss_avoid_adj"),
  purrr::map_dfr(prim_res, mk_part_row, tg = "ev_sens_agn")) |>
  dplyr::arrange(target, partition)

# ---- auxiliary K: N_PART_AUX partitions per (seed, K) ----------------------
aux_jobs <- list()
for (K in setdiff(K_GRID, GATE_K)) for (si in seq_len(R_SEED)) {
  sets <- draw_sets(length(form_mats[[GATE_FORM]]$probs), K, M_CAP,
                    seed_set_b1(GATE_FORM, K, si))
  for (j in seq_len(N_PART_AUX)) {
    aux_jobs[[length(aux_jobs) + 1L]] <- list(
      K = K, si = si, j = j, sets = sets,
      seed_part = seed_part_aux(K, si, j),
      part_index = 200L + 10L * si + j,
      seed_perm = seed_perm_b1(200L + 10L * si + j))
  }
}
t0 <- Sys.time()
aux_res <- parallel::parLapply(cl, aux_jobs, function(job) {
  fm <- form_mats[[GATE_FORM]]
  cc <- cov_cell(fm, job$sets, seq_along(fm$subs), job$seed_part,
                 job$part_index, job$K, N_SPLIT_B1)
  list(K = job$K, si = job$si, j = job$j, n_a = cc$n_a, n_b = cc$n_b,
       loss = cov_stats(cc$r[, 1L], cc$r[, 2L], job$seed_perm),
       ev = cov_stats(cc$r[, 3L], cc$r[, 4L], job$seed_perm + 500000L))
})
say("auxiliary-K covariance: %d partition-cells in %.1f min\n", length(aux_jobs),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

cov_long <- dplyr::bind_rows(
  b1_partition_distribution |>
    dplyr::select(K, seed_index, target, n_sets_used, cov_set, var_A, var_B,
                  sigma_set_raw_cov, sigma_set_true_cov, ICC_cov,
                  sigma_set_true_cov_null),
  purrr::map_dfr(c("loss_avoid_adj", "ev_sens_agn"), function(tg)
    purrr::map_dfr(aux_res, function(z) {
      s <- if (tg == "loss_avoid_adj") z$loss else z$ev
      tibble::tibble(K = z$K, seed_index = z$si, target = tg,
                     n_sets_used = s$n_used, cov_set = s$cov,
                     var_A = s$var_a, var_B = s$var_b,
                     sigma_set_raw_cov = s$sigma_raw,
                     sigma_set_true_cov = s$sigma_true, ICC_cov = s$icc,
                     sigma_set_true_cov_null = s$sigma_true_null)
    })))

est_ii <- cov_long |>
  dplyr::group_by(seed_index, K, target) |>
  dplyr::summarise(
    form = GATE_FORM, estimand = "covariance_subject_split",
    n_partitions = dplyr::n(),
    M = M_CAP, n_sets_identified = round(mean(n_sets_used)),
    set_coverage = mean(n_sets_used) / M_CAP,
    split_coverage = NA_real_, n_subjects_median = NA_real_,
    sigma2_set_raw = mean(sigma_set_raw_cov^2, na.rm = TRUE),
    sigma_set_raw = sqrt(mean(sigma_set_raw_cov^2, na.rm = TRUE)),
    sigma2_noise = NA_real_, sigma_noise = NA_real_,
    sigma2_set_true = mean(pmax(0, cov_set), na.rm = TRUE),
    sigma_set_true = mean(sigma_set_true_cov, na.rm = TRUE),
    ICC_set = mean(ICC_cov, na.rm = TRUE),
    sigma_set_true_null = mean(sigma_set_true_cov_null, na.rm = TRUE),
    r_sb_min = NA_real_, r_sb_median = NA_real_, r_sb_max = NA_real_,
    best_minus_median = NA_real_, n_noise_sets = NA_integer_,
    is_anchor = FALSE, .groups = "drop")

b1_stability_by_seed <- dplyr::bind_rows(
  est_i |> dplyr::mutate(sigma_set_true_null = NA_real_),
  est_ii) |>
  dplyr::select(form, K, target, seed_index, estimand, is_anchor, M,
                n_sets_identified, set_coverage, split_coverage,
                n_subjects_median, n_partitions, sigma_set_raw, sigma_noise,
                sigma_set_true, ICC_set, sigma2_set_raw, sigma2_noise,
                sigma2_set_true, sigma_set_true_null, r_sb_min, r_sb_median,
                r_sb_max, best_minus_median, n_noise_sets) |>
  dplyr::arrange(target, K, estimand, seed_index)
save_csv(b1_stability_by_seed, "b1_stability_by_seed")
save_csv(b1_partition_distribution, "b1_partition_distribution")

# =============================================================================
# 7. Full feasible-cell map under the covariance estimand (G-B1d, REPORTED)
# =============================================================================
# Sets are E3's committed sets for EVERY form and K (parsed from the committed
# per-set file), so the map is the SAME cell grid E3 reported 6/50 on, scored
# under the covariance estimand instead of the bootstrap floor.
say("\n=== B1 full-map under the covariance estimand ===\n")
map_cells <- e3_real |>
  dplyr::filter(target == TARGETS[1]) |>
  dplyr::distinct(form, K) |>
  dplyr::arrange(form, K)
say("feasible cells: %d (x %d targets = %d)\n", nrow(map_cells),
    length(TARGETS), nrow(map_cells) * length(TARGETS))

map_jobs <- list()
for (ci in seq_len(nrow(map_cells))) {
  S <- as.integer(map_cells$form[ci]); K <- as.integer(map_cells$K[ci])
  d <- e3_real |> dplyr::filter(form == S, K == !!K, target == TARGETS[1]) |>
    dplyr::arrange(set_id)
  cols <- lapply(d$game_ids, function(s) {
    cc <- match(parse_ids(s), form_mats[[S]]$probs)
    stopifnot(!anyNA(cc), !anyDuplicated(cc))
    sort(cc)
  })
  for (j in seq_len(N_PART_MAP)) {
    map_jobs[[length(map_jobs) + 1L]] <- list(
      form = S, K = K, j = j, sets = cols,
      seed_part = seed_part_map(S, K, j),
      part_index = 1000L + 100L * S + j,
      seed_perm = seed_perm_b1(1000L + 100L * S + j))
  }
}
t0 <- Sys.time()
map_res <- parallel::parLapply(cl, map_jobs, function(job) {
  fm <- form_mats[[job$form]]
  cc <- cov_cell(fm, job$sets, seq_along(fm$subs), job$seed_part,
                 job$part_index, job$K, N_SPLIT_B1)
  list(form = job$form, K = job$K, j = job$j, n_a = cc$n_a, n_b = cc$n_b,
       n_sets = length(job$sets),
       loss = cov_stats(cc$r[, 1L], cc$r[, 2L], job$seed_perm),
       ev = cov_stats(cc$r[, 3L], cc$r[, 4L], job$seed_perm + 500000L))
})
say("map: %d partition-cells in %.1f min\n", length(map_jobs),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
parallel::stopCluster(cl)

b1_map_cov <- purrr::map_dfr(TARGETS, function(tg)
  purrr::map_dfr(map_res, function(z) {
    s <- if (tg == "loss_avoid_adj") z$loss else z$ev
    tibble::tibble(form = z$form, K = z$K, target = tg, partition = z$j,
                   n_sets = z$n_sets, n_sets_used = s$n_used,
                   n_subjects_A = z$n_a, n_subjects_B = z$n_b,
                   cov_set = s$cov, sigma_set_raw_cov = s$sigma_raw,
                   sigma_set_true_cov = s$sigma_true, ICC_cov = s$icc,
                   spearman_A_vs_B = s$rho,
                   sigma_set_true_cov_null = s$sigma_true_null)
  })) |>
  dplyr::group_by(form, K, target) |>
  dplyr::summarise(
    # NOTE: every summary below reads the ORIGINAL per-partition columns, so
    # the quantiles and the evaluable count are computed BEFORE any name is
    # reused for its own mean (dplyr::summarise makes earlier results visible
    # to later expressions).
    n_partitions = dplyr::n(), n_sets = dplyr::first(n_sets),
    mean_n_sets_used = mean(n_sets_used),
    n_subjects_A = dplyr::first(n_subjects_A),
    n_subjects_B = dplyr::first(n_subjects_B),
    n_partitions_evaluable = sum(is.finite(ICC_cov)),
    sigma_set_true_cov_lo = q(sigma_set_true_cov, .05),
    sigma_set_true_cov_hi = q(sigma_set_true_cov, .95),
    sigma_set_raw_cov = mean(sigma_set_raw_cov, na.rm = TRUE),
    sigma_set_true_cov = mean(sigma_set_true_cov, na.rm = TRUE),
    ICC_cov = mean(ICC_cov, na.rm = TRUE),
    spearman_A_vs_B = mean(spearman_A_vs_B, na.rm = TRUE),
    sigma_set_true_cov_null = mean(sigma_set_true_cov_null, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(
    clears_sigma = is.finite(sigma_set_true_cov) &
      sigma_set_true_cov >= H1_SIGMA_MIN,
    clears_icc = is.finite(ICC_cov) & ICC_cov >= H1_ICC_MIN,
    clears_both = clears_sigma & clears_icc) |>
  dplyr::arrange(target, form, K)
save_csv(b1_map_cov, "b1_map_cov")

# E3's committed bootstrap-floor map, for the side-by-side count.
e3_map <- readr::read_csv(
  here("output", "elevation", "results", "e3", "heterogeneity_summary.csv"),
  show_col_types = FALSE)
e3_pass_both <- sum(e3_map$sigma_set_true >= H1_SIGMA_MIN &
                      e3_map$ICC_set >= H1_ICC_MIN, na.rm = TRUE)
say("map: covariance estimand clears both thresholds in %d of %d cells; E3 bootstrap floor %d of %d\n",
    sum(b1_map_cov$clears_both), nrow(b1_map_cov), e3_pass_both, nrow(e3_map))

# =============================================================================
# 8. HG-ceiling guard
# =============================================================================
# Two readings of "no recomputed reliability exceeds the E0 true-variance
# ceiling beyond bootstrap noise" are computed and BOTH reported. Neither is
# adjusted after the fact; which one fires is a result.
#   PER-SET (strict): the frozen reliability_point CI lower bound of an
#     individual set clears the ceiling (E8's exceedance rule, applied set by
#     set). Note in advance that this reading is in structural tension with H1
#     itself: a heterogeneity result asserts that some item sets are better
#     than the average, and the E0 ceiling is an average over item sets in a
#     DIFFERENT design (Study 1's pooled 686 subjects, each with their own
#     problems), so individual-set exceedance is expected under H1 and is not
#     by itself evidence of leakage.
#   PER-CELL (central): the cell's MEDIAN r_sb over the M sampled sets clears
#     the ceiling by more than the cell's own estimation noise (median minus
#     sigma_noise still above the ceiling). This is the reading that compares
#     like with like, because the E0 ceiling bounds the average reliability
#     attainable at K problems, not the best item set.
# Whatever fires is recorded and is NOT worked around here.
hg <- noise_i |>
  dplyr::mutate(ceiling_K = sb_k(K, unname(R1[target])),
                exceeds_point = r_sb > ceiling_K,
                exceeds_beyond_noise = ci_lo_sb > ceiling_K)
hg_set <- hg |>
  dplyr::group_by(target, K) |>
  dplyr::summarise(n_sets_with_ci = dplyr::n(),
                   ceiling_K = dplyr::first(ceiling_K),
                   frac_point_exceed_bootstrapped_sets =
                     mean(exceeds_point, na.rm = TRUE),
                   n_ci_supported_exceedance = sum(exceeds_beyond_noise,
                                                   na.rm = TRUE),
                   .groups = "drop")
hg_cell <- est_i |>
  dplyr::mutate(ceiling_K = sb_k(K, unname(R1[target])),
                median_exceeds_point = r_sb_median > ceiling_K,
                median_exceeds_beyond_noise =
                  (r_sb_median - sigma_noise) > ceiling_K) |>
  dplyr::select(form, K, target, seed_index, is_anchor, ceiling_K, r_sb_median,
                r_sb_max, sigma_noise, median_exceeds_point,
                median_exceeds_beyond_noise)
b1_hg_ceiling <- dplyr::full_join(
  hg_set,
  hg_cell |> dplyr::group_by(target, K) |>
    dplyr::summarise(n_cells = dplyr::n(),
                     median_r_sb_over_cells = stats::median(r_sb_median),
                     n_cells_median_exceeds_point =
                       sum(median_exceeds_point, na.rm = TRUE),
                     n_cells_median_exceeds_beyond_noise =
                       sum(median_exceeds_beyond_noise, na.rm = TRUE),
                     .groups = "drop"),
  by = c("target", "K"))
print(as.data.frame(b1_hg_ceiling))
hg_set_any  <- any(hg$exceeds_beyond_noise, na.rm = TRUE)
hg_cell_any <- any(hg_cell$median_exceeds_beyond_noise, na.rm = TRUE)
say("HG-ceiling (B1) per-set reading: %d of %d bootstrapped sets CI-clear the ceiling -> %s\n",
    sum(hg$exceeds_beyond_noise, na.rm = TRUE), nrow(hg),
    if (hg_set_any) "FLAGGED" else "clear")
say("HG-ceiling (B1) per-cell reading: %d of %d cells have a median above the ceiling beyond noise -> %s\n",
    sum(hg_cell$median_exceeds_beyond_noise, na.rm = TRUE), nrow(hg_cell),
    if (hg_cell_any) "FLAGGED" else "clear")
save_csv(b1_hg_ceiling, "b1_hg_ceiling")
save_csv(hg_cell, "b1_hg_ceiling_by_cell")

# =============================================================================
# 9. Exports for independent re-derivation
# =============================================================================
save_csv(per_set_i |> dplyr::arrange(seed_index, K, set_id),
         "b1_export_per_set_full")
save_csv(noise_i |> dplyr::arrange(seed_index, K, target, set_id),
         "b1_export_noise_floor")
EXPORT_PARTS <- 1:10       # fixed sample of partitions exported per set
save_csv(dplyr::bind_rows(lapply(prim_res[EXPORT_PARTS], function(z) z$export)),
         "b1_export_per_set_halves")

# =============================================================================
# 10. Gates G-B1a..d
# =============================================================================
batt <- b1_stability_by_seed |>
  dplyr::filter(form == GATE_FORM, K == GATE_K, target == GATE_TARGET,
                !is_anchor)
stopifnot(nrow(batt) == 2L * R_SEED)          # R_SEED seeds x 2 estimands
g_b1a_sigma_min <- min(batt$sigma_set_true, na.rm = TRUE)
g_b1a_icc_min   <- min(batt$ICC_set, na.rm = TRUE)
g_b1a <- isTRUE(g_b1a_sigma_min >= H1_SIGMA_MIN) &&
  isTRUE(g_b1a_icc_min >= H1_ICC_MIN)

prim_loss <- b1_partition_distribution |>
  dplyr::filter(target == GATE_TARGET)
p90 <- unname(stats::quantile(prim_loss$sigma_set_true_cov, c(0.05, 0.95),
                              na.rm = TRUE))
g_b1b <- isTRUE(p90[1] >= H1_SIGMA_MIN)

null_mean <- mean(prim_loss$sigma_set_true_cov_null, na.rm = TRUE)
real_mean <- mean(prim_loss$sigma_set_true_cov, na.rm = TRUE)
g_b1c <- isTRUE(null_mean < H1_SIGMA_MIN) && isTRUE(null_mean < 0.5 * real_mean)

g_b1d_n <- sum(b1_map_cov$clears_both)

say("\n=== B1 GATES ===\n")
say("G-B1a primary-cell stability: min sigma_set_true %.5f (>= %.2f) AND min ICC_set %.4f (>= %.2f) -> %s\n",
    g_b1a_sigma_min, H1_SIGMA_MIN, g_b1a_icc_min, H1_ICC_MIN,
    if (g_b1a) "PASS" else "FAIL")
say("G-B1b partition stability: 90%% interval of sigma_set_true_cov = [%.5f, %.5f]; lower bound >= %.2f -> %s\n",
    p90[1], p90[2], H1_SIGMA_MIN, if (g_b1b) "PASS" else "FAIL")
say("G-B1c null control: permuted sigma_set_true_cov mean %.5f vs real %.5f -> %s\n",
    null_mean, real_mean, if (g_b1c) "PASS" else "FAIL")
say("G-B1d map (REPORTED, not gated): %d of %d cells clear both under covariance; E3 bootstrap floor %d of %d\n",
    g_b1d_n, nrow(b1_map_cov), e3_pass_both, nrow(e3_map))

part <- list(
  gate = list(
    probe = "B1 heterogeneity stability",
    primary_cell = list(form = GATE_FORM, K = GATE_K, target = GATE_TARGET),
    thresholds = list(sigma_set_true_min = H1_SIGMA_MIN,
                      ICC_set_min = H1_ICC_MIN),
    frozen_e3_reference = list(
      sigma_set_true = 0.187659511205571, ICC_set = 0.6258254104076,
      best_minus_median = 0.195354087580587,
      source = "elevation/results/e3/h1_gate.json (read-only; NOT re-run)"),
    anchor_reproduction = list(
      seed_index = 0L,
      max_abs_diff_per_set_r_sb_loss_avoid_adj = anchor_d_loss,
      max_abs_diff_per_set_r_sb_ev_sens_agn = anchor_d_ev,
      na_pattern_identical = anchor_na_ok,
      sigma_set_true = anchor_row$sigma_set_true,
      ICC_set = anchor_row$ICC_set,
      excluded_from_gate_battery = TRUE),
    `G-B1a` = list(
      battery_size = nrow(batt), r_seed = R_SEED, n_estimands = 2L,
      min_sigma_set_true = g_b1a_sigma_min, min_ICC_set = g_b1a_icc_min,
      per_estimand = lapply(split(batt, batt$estimand), function(d)
        list(estimand = d$estimand[1], min_sigma_set_true = min(d$sigma_set_true),
             max_sigma_set_true = max(d$sigma_set_true),
             min_ICC_set = min(d$ICC_set), max_ICC_set = max(d$ICC_set))),
      verdict = if (g_b1a) "PASS" else "FAIL", pass = g_b1a),
    `G-B1b` = list(
      n_partitions = nrow(prim_loss),
      mean_sigma_set_true_cov = real_mean,
      interval_90 = list(lo = p90[1], hi = p90[2]),
      threshold = H1_SIGMA_MIN,
      verdict = if (g_b1b) "PASS" else "FAIL", pass = g_b1b),
    `G-B1c` = list(
      permuted_sigma_set_true_cov_mean = null_mean,
      real_sigma_set_true_cov_mean = real_mean,
      criterion = "permuted mean < 0.05 AND permuted mean < half the real mean",
      verdict = if (g_b1c) "PASS" else "FAIL", pass = g_b1c),
    `G-B1d` = list(
      reported_not_gated = TRUE,
      n_cells = nrow(b1_map_cov), n_cells_clearing_both = g_b1d_n,
      e3_bootstrap_floor_cells_clearing_both = e3_pass_both,
      e3_n_cells = nrow(e3_map)),
    HG_ceiling = list(
      comparator = paste("E0 full-data ceiling (k_ref 14) SB-projected in",
                         "problem count to the cell K"),
      ceiling_r1 = as.list(R1),
      per_set_reading = list(
        rule = "frozen reliability_point CI lower bound above the ceiling (the E8 rule), set by set",
        n_bootstrapped_sets = nrow(hg),
        n_ci_supported_exceedance = sum(hg$exceeds_beyond_noise, na.rm = TRUE),
        frac_point_exceedance = mean(hg$exceeds_point, na.rm = TRUE),
        flagged = hg_set_any,
        caveat = paste("individual-set exceedance is EXPECTED under H1: the E0",
                       "ceiling is an average over item sets in Study 1's",
                       "pooled design, and H1 asserts item sets differ")),
      per_cell_reading = list(
        rule = "cell median r_sb minus the cell sigma_noise still above the ceiling",
        n_cells = nrow(hg_cell),
        n_flagged = sum(hg_cell$median_exceeds_beyond_noise, na.rm = TRUE),
        n_median_above_point = sum(hg_cell$median_exceeds_point, na.rm = TRUE),
        flagged = hg_cell_any),
      flagged_any = hg_set_any || hg_cell_any,
      disposition = "reported; not worked around")),
  config = list(
    probe = "B1", seed_base = SEED_B1,
    regime = "experience (Trial 6-25)", r_reps = R_REPS,
    K_grid = K_GRID, targets = list(primary = TARGETS[1], secondary = TARGETS[2]),
    M_cap = M_CAP, r_seed = R_SEED, n_split_b1 = N_SPLIT_B1,
    n_noise_sets = N_NOISE_SETS, n_boot = N_BOOT, b_part_b1 = B_PART_B1,
    n_part_aux = N_PART_AUX, n_part_map = N_PART_MAP,
    min_sets_cov = MIN_SETS_COV, min_pairs = MIN_PAIRS,
    partition_budget_rule = paste(
      "DISPB_spec section 2 asks for the covariance estimand at every K and a",
      ">= 200-partition distribution; 200 partitions at every (seed, K) is",
      "~6400 partition-cells and is not affordable in one run. Fixed",
      "BEFORE any result: the full 200 partitions run at the PRIMARY cell",
      "(K = 12), with partition p evaluated on set sample ((p-1) mod 8) + 1 so",
      "every seed contributes exactly 25 partitions and the distribution",
      "marginalises over BOTH analyst choices; auxiliary K use 8 partitions per",
      "(seed, K); the full map uses 8 partitions per cell."),
    anchor_rule = paste(
      "seed_index 0 replays E3's committed Set-1 sets with E3's own split and",
      "noise seeds. It is an ANCHOR that must reproduce the frozen h1_gate and",
      "is EXCLUDED from the G-B1a battery minimum."),
    accelerator = paste(
      "The two targets are evaluated by a vectorised path using the SAME",
      "formulas as the frozen estimate_half, because estimate_half also runs",
      "the 4-parameter Firth fit for b_ev/b_risk/b_loss, which are not",
      "targets here. dispb_verify.R re-proves agreement to < 1e-12 with identical",
      "NA patterns before any result exists, and closes the loop on whole-set",
      "r_sb through the frozen engine on both the full sample and a half."),
    data_inputs = c("data/raw/All estimation raw data.csv",
                    "elevation/results/e1/participant_split.csv",
                    "elevation/results/e2/problem_information_dictionary.csv",
                    "elevation/results/e3/reliability_by_set_real.csv",
                    "elevation/results/e3/heterogeneity_summary.csv",
                    "elevation/results/e0/variance_ceiling_final.csv"),
    split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
    loss_baselines_exp = as.list(baselines$exp),
    build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
    optimization_performed = FALSE, candidate_design_constructed = FALSE,
    new_choice_model_fitted = FALSE, e3_h1_gate_rerun_or_amended = FALSE,
    runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
jsonlite::write_json(part, file.path(res_dir, "dispb_part_b1.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

seeds_b1 <- dplyr::bind_rows(
  purrr::map_dfr(si_grid, function(si) purrr::map_dfr(K_GRID, function(K) {
    sv <- c(if (si == 0L) seed_subset_e3(GATE_FORM, K) else
              seed_set_b1(GATE_FORM, K, si),
            if (si == 0L) seed_split_e3(GATE_FORM, K) else
              seed_split_b1(GATE_FORM, K, si),
            if (si == 0L) seed_noise_e3(GATE_FORM, K) else
              seed_noise_b1(GATE_FORM, K, si))
    tibble::tibble(
      probe = "B1", scope = sprintf("form%d|K%d|seed%d", GATE_FORM, K, si),
      purpose = c("set_draw", "split_draw_base", "noise_bootstrap_base"),
      seed = sv,
      derivation = if (si == 0L) c(
        "E3 replay: SEED_E3 + 1000*form + K",
        "E3 replay: SEED_E3 + 1000*form + K + 500000; per set = base + set_id",
        "E3 replay: SEED_E3 + 600000 + 1000*form + K; per set = base + set_id")
      else c(
        "SEED_B1 + 1000*form + K + 10000*seed_index",
        "SEED_B1 + 1000*form + K + 10000*seed_index + 500000; per set = base + set_id",
        "SEED_B1 + 600000 + 1000*form + K + 10000*seed_index; per set = base + set_id"))
  })),
  tibble::tibble(
    probe = "B1", scope = sprintf("primary|partition%d", seq_len(B_PART_B1)),
    purpose = "subject_partition", seed = seed_part_b1(seq_len(B_PART_B1)),
    derivation = "SEED_B1 + 700000 + partition_index"),
  purrr::map_dfr(setdiff(K_GRID, GATE_K), function(K)
    purrr::map_dfr(seq_len(R_SEED), function(si) tibble::tibble(
      probe = "B1", scope = sprintf("auxK%d|seed%d|partition%d", K, si,
                                    seq_len(N_PART_AUX)),
      purpose = "subject_partition",
      seed = seed_part_aux(K, si, seq_len(N_PART_AUX)),
      derivation = "SEED_B1 + 720000 + 1000*K + 10*seed_index + j"))),
  purrr::map_dfr(seq_len(nrow(map_cells)), function(ci) tibble::tibble(
    probe = "B1",
    scope = sprintf("map|form%d|K%d|partition%d", map_cells$form[ci],
                    map_cells$K[ci], seq_len(N_PART_MAP)),
    purpose = "subject_partition",
    seed = seed_part_map(as.integer(map_cells$form[ci]),
                         as.integer(map_cells$K[ci]), seq_len(N_PART_MAP)),
    derivation = "SEED_B1 + 750000 + 1000*form + 10*K + j")),
  tibble::tibble(
    probe = "B1", scope = "halving", purpose = "half_A_and_half_B_halvings",
    seed = NA_integer_,
    derivation = paste("SEED_B1 + 900000000*(half-1) + 200000000*K_index +",
                       "100000*partition_index + set_id; partition_index is",
                       "1..200 (primary), 200+10*seed_index+j (auxiliary K),",
                       "1000+100*form+j (map)")),
  tibble::tibble(
    probe = "B1", scope = "null", purpose = "set_correspondence_permutation",
    seed = NA_integer_,
    derivation = "SEED_B1 + 800000 + partition_index (+500000 for the secondary target)"))
readr::write_csv(seeds_b1, file.path(res_dir, "dispb_seeds_b1.csv"))

# Seed hygiene: every explicitly enumerated seed must be unique within B1.
enum <- seeds_b1$seed[is.finite(seeds_b1$seed)]
stopifnot(!anyDuplicated(enum))
half_enum <- unlist(lapply(c(1L, 2L), function(h)
  unlist(lapply(K_GRID, function(K)
    unlist(lapply(c(seq_len(B_PART_B1), 200L + 10L * rep(seq_len(R_SEED), each = N_PART_AUX) +
                      rep(seq_len(N_PART_AUX), R_SEED),
                    1000L + 100L * rep(1:7, each = N_PART_MAP) +
                      rep(seq_len(N_PART_MAP), 7)), function(p)
                        seed_half_b1(h, K, p, c(1L, M_CAP))))))))
stopifnot(!anyDuplicated(half_enum), max(half_enum) < .Machine$integer.max)
say("seed hygiene: %d enumerated seeds unique; halving-seed grid unique, max %d\n",
    length(enum), max(half_enum))

say("\nB1 complete in %.1f min. Outputs in elevation/results/dispb/\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))

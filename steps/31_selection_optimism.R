# Step 31: optimism of selecting the best subset
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e0/, e1/, e2/, e3/ (steps 19-22)
# Writes:  output/elevation/results/dispb/ (b2_*.csv, dispb_part_b2.json, seeds)
# Paper:   Study 2 (Ruling out alternatives); Supplement (Study 2 robustness detail), through steps 33 and 36
# Runtime: about 48 minutes with 12 parallel workers

# B2: quantify selection optimism with pre-declared gates G-B2a..c.
# E3's best-minus-median gap of 0.195 in r_sb used the same 63 subjects to
# select and score a set. Here sets are ranked on one subject half and scored
# on a disjoint half, over repeated partitions and in both directions.
#
# The analysis re-scores the committed E3 sets without optimizing or
# constructing a design or fitting a new choice model. The frozen E3 h1_gate
# and E4 h2_gate are neither rerun nor amended.
#
# Only the 63 Set-1 build subjects from step 20 are used. Choice/outcome
# columns are dropped from the full-sample frame at load, excluding the 137
# freeze and 205 test subjects' choices.
#
# Shared functions are sourced without redefinition. Each r_sb uses
# split_correlations and sb_correct; the strict ceiling reading uses
# reliability_point. Step 33 checks the accelerator from steps 22 and 24
# against estimate_half.

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
RAW_MD5     <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_B2     <- 20260723L
REGIME      <- "exp"
POS         <- 6:25
R_REPS      <- length(POS)
K_GRID_B2   <- c(8L, 10L, 12L)     # K = 12 primary
TARGETS     <- c("loss_avoid_adj", "ev_sens_agn")
M_CAP       <- 2000L
B_PART_B2   <- 200L
N_SPLIT_B2  <- 50L                 # bounds runtime; MC error reported
DECILE      <- 0.10                # top 10% by selection-half r_sb
MIN_PAIRS   <- 10L
GATE_FORM   <- 1L
GATE_K      <- 12L
GATE_TARGET <- "loss_avoid_adj"
N_HALF_A    <- 32L                 # largest-remainder halves of 63
N_WORKERS   <- 12L

# Further constants, fixed HERE before any result exists.
N_PERM_B2    <- 20L    # evaluation-half label permutations per partition
MIN_SETS_B2  <- 30L    # identifiable sets needed for a partition to be scored
NULL_ABS_TOL <- 0.01   # "within Monte-Carlo noise of 0" for G-B2c: the 90%
                       # interval of the permuted gain must cover 0 AND the
                       # mean magnitude must be below this absolute tolerance
N_CEIL_SETS  <- 50L    # sets given a frozen reliability_point CI on each half
N_BOOT_CEIL  <- 1000L  # matches E3's noise-floor N_BOOT
EXPORT_PARTS <- 1:5    # partitions exported per-set

# ---- seed derivations (blocks asserted disjoint) ----------------------------
# Pre-declared: partition seed = SEED_B2 + partition_index.
seed_part_b2 <- function(p) SEED_B2 + p
# Fixed in code, derived as in step 24 (disjoint per-half blocks, asserted
# below):
seed_half_b2 <- function(h, K, p, sid)
  SEED_B2 + 900000000L * (h - 1L) + 200000000L * match(K, K_GRID_B2) +
  100000L * p + sid
seed_perm_b2 <- function(K, p) SEED_B2 + 800000L + 1000L * match(K, K_GRID_B2) + p
seed_ceil_b2 <- function(h, sid) SEED_B2 + 600000L + 100000L * h + sid

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
stopifnot(all(split_tbl$held_out_until_E7 == (split_tbl$split == "test")))
say("split: build %d / freeze %d (not read) / test %d (SEALED)\n",
    n_build, n_freeze, n_test)

raw <- load_cpc18(raw_path)
choice_cols <- c("B", "Button", "Payoff", "Forgone", "RT", "Apay", "Bpay",
                 "Feedback")
design_all <- raw |> dplyr::select(-dplyr::any_of(choice_cols))
stopifnot(!any(choice_cols %in% names(design_all)))                    # build-only guard
build_raw <- raw |> dplyr::filter(SubjID %in% build_ids)
rm(raw); invisible(gc())
stopifnot(!any(build_raw$SubjID %in% held_out))                        # build-only guard

feat_all <- clean_problem_features(design_all) |> add_scaled_features()
feat_all <- feat_all |> dplyr::arrange(GameID)
stopifnot(nrow(feat_all) == 105L)
e2_dict <- readr::read_csv(
  here("output", "elevation", "results", "e2", "problem_information_dictionary.csv"),
  show_col_types = FALSE)
stopifnot(nrow(e2_dict) == 105L, setequal(e2_dict$GameID, feat_all$GameID))
stopifnot(all(abs(e2_dict$d_ev[order(e2_dict$GameID)] - feat_all$d_ev) < 1e-12))

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
stopifnot(abs(baselines$exp[1] - 0.4301295476692) < 1e-10,
          abs(baselines$exp[2] - 0.165107793793507) < 1e-10)
say("frozen loss baselines (exp): b0=%.13f b1=%.15f\n",
    baselines$exp[1], baselines$exp[2])

feat_all <- feat_all |>
  dplyr::mutate(
    d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a),
    p_base_exp  = ifelse(loss_informative,
                         stats::plogis(baselines$exp[1] +
                                         baselines$exp[2] * d_ev_noloss), NA))

s1_trials <- build_trials |> dplyr::filter(Set == GATE_FORM)
s1_subs <- as.character(sort(unique(s1_trials$SubjID)))
s1_probs <- sort(unique(s1_trials$GameID))
stopifnot(length(s1_subs) == 63L)
NB <- matrix(NA_real_, length(s1_subs), length(s1_probs),
             dimnames = list(s1_subs, as.character(s1_probs)))
for (sid in s1_subs) {
  sd <- subj_data[[sid]]
  j <- match(sd$feat$GameID, s1_probs)
  ok <- !is.na(j)
  NB[sid, j[ok]] <- rowSums(sd$B_mat[which(ok), POS, drop = FALSE])
}
stopifnot(!anyNA(NB))
fi <- feat_all[match(s1_probs, feat_all$GameID), ]
fm <- list(probs = s1_probs, subs = s1_subs, NB = NB, feat = fi,
           li = fi$loss_informative, has_loss_b = fi$has_loss_b,
           d_ev = fi$d_ev, p_base = fi$p_base_exp)
say("Set 1: %d build subjects, %d clean problems, %d loss-informative\n",
    length(s1_subs), length(s1_probs), sum(fm$li))

# =============================================================================
# 2. Accelerator (identical formulas to the frozen estimate_half; see header)
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

run_set_rows <- function(fm, cols, rows, n_split, seed, want_est = FALSE) {
  NBs <- fm$NB[rows, cols, drop = FALSE]
  ns <- nrow(NBs); K <- length(cols); kh <- K %/% 2L
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
    rs <- split_correlations(E[[tg]]$h1, E[[tg]]$h2, MIN_PAIRS)         # FROZEN
    r_raw <- mean(rs, na.rm = TRUE)
    out[[tg]] <- list(r_sb = if (is.finite(r_raw)) sb_correct(r_raw) else
      NA_real_)                                                         # FROZEN
  }
  if (want_est) out$E <- E
  out
}

# =============================================================================
# 3. The committed E3 item sets (the SAME sets E3 measured the naive gap on)
# =============================================================================
e3_real <- readr::read_csv(
  here("output", "elevation", "results", "e3", "reliability_by_set_real.csv"),
  show_col_types = FALSE)
e3_s1 <- e3_real |> dplyr::filter(form == GATE_FORM)
parse_ids <- function(s) as.integer(strsplit(s, ";", fixed = TRUE)[[1]])
sets_by_K <- lapply(K_GRID_B2, function(K) {
  d <- e3_s1 |> dplyr::filter(target == TARGETS[1], K == !!K) |>
    dplyr::arrange(set_id)
  stopifnot(identical(as.integer(d$set_id), seq_len(M_CAP)))
  chk <- e3_s1 |> dplyr::filter(target == TARGETS[2], K == !!K) |>
    dplyr::arrange(set_id)
  stopifnot(identical(d$game_ids, chk$game_ids))            # set identity
  lapply(d$game_ids, function(s) {
    cc <- match(parse_ids(s), fm$probs)
    stopifnot(!anyNA(cc), !anyDuplicated(cc))
    sort(cc)
  })
})
names(sets_by_K) <- as.character(K_GRID_B2)
say("B2 item sets: %d committed E3 Set-1 sets per K (K = %s)\n", M_CAP,
    paste(K_GRID_B2, collapse = ", "))

# E3's committed naive (full 63, in-sample) best-minus-median reference.
e3_het <- readr::read_csv(
  here("output", "elevation", "results", "e3", "heterogeneity_summary.csv"),
  show_col_types = FALSE)
naive_ref <- e3_het |>
  dplyr::filter(form == GATE_FORM, K %in% K_GRID_B2) |>
  dplyr::select(K, target, naive_best_minus_median = best_minus_median,
                naive_r_sb_median = r_sb_median, naive_r_sb_max = r_sb_max,
                naive_sigma_set_raw = sigma_set_raw)
print(as.data.frame(naive_ref))

# E0 ceiling anchors (HG-ceiling comparator; same construction as B1 and E8).
e0 <- readr::read_csv(
  here("output", "elevation", "results", "e0", "variance_ceiling_final.csv"),
  show_col_types = FALSE)
r1_of <- function(R, k) R / (k - (k - 1) * R)
sb_k  <- function(k, r1) k * r1 / (1 + (k - 1) * r1)
R1 <- vapply(TARGETS, function(tg) {
  rw <- e0[e0$construct == tg, ]
  stopifnot(nrow(rw) == 1L, rw$k_ref == 14)
  r1_of(rw$ceiling_fulldata, 14)
}, numeric(1))

# =============================================================================
# 4. Per-partition selection / evaluation
# =============================================================================
partition_rows <- function(n, seed, n_a) {
  set.seed(seed)
  p <- sample.int(n)
  list(A = sort(p[seq_len(n_a)]), B = sort(p[(n_a + 1L):n]))
}

# Selection statistics for ONE direction. `sel` ranks, `evl` scores.
gain_stats <- function(sel, evl, n_perm, seed) {
  ok <- is.finite(sel) & is.finite(evl)
  n_used <- sum(ok)
  na_out <- list(n_used = n_used, in_single = NA_real_, out_single = NA_real_,
                 in_decile = NA_real_, out_decile = NA_real_,
                 rho = NA_real_, n_top = NA_integer_,
                 null_single = NA_real_, null_decile = NA_real_)
  if (n_used < MIN_SETS_B2) return(na_out)
  a <- sel[ok]; b <- evl[ok]
  n_top <- max(1L, ceiling(DECILE * n_used))
  top <- order(a, decreasing = TRUE)[seq_len(n_top)]
  best <- which.max(a)
  in_single  <- a[best] - stats::median(a)
  out_single <- b[best] - stats::median(b)
  in_decile  <- mean(a[top]) - mean(a)
  out_decile <- mean(b[top]) - mean(b)
  set.seed(seed)
  nl <- vapply(seq_len(n_perm), function(z) {
    bp <- sample(b)
    c(bp[best] - stats::median(bp), mean(bp[top]) - mean(bp))
  }, numeric(2))
  list(n_used = n_used, in_single = in_single, out_single = out_single,
       in_decile = in_decile, out_decile = out_decile,
       rho = suppressWarnings(stats::cor(a, b, method = "spearman")),
       n_top = n_top,
       null_single = mean(nl[1, ]), null_decile = mean(nl[2, ]))
}

say("\n=== B2: %d partitions x %d K, N_SPLIT_B2 = %d ===\n",
    B_PART_B2, length(K_GRID_B2), N_SPLIT_B2)
cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(tidyverse); library(here) })
  source(here("R", "step2_reliability.R"))
  NULL
}))
parallel::clusterExport(cl, c("fm", "half_targets", "ev_ranks", "run_set_rows",
                              "partition_rows", "gain_stats", "seed_half_b2",
                              "seed_perm_b2", "TARGETS", "R_REPS", "MIN_PAIRS",
                              "N_SPLIT_B2", "N_HALF_A", "DECILE", "MIN_SETS_B2",
                              "N_PERM_B2", "K_GRID_B2", "SEED_B2"),
                        envir = environment())

jobs <- list()
for (K in K_GRID_B2) for (p in seq_len(B_PART_B2)) {
  jobs[[length(jobs) + 1L]] <- list(K = K, p = p, sets = sets_by_K[[as.character(K)]])
}
t0 <- Sys.time()
res <- parallel::parLapply(cl, jobs, function(job) {
  K <- job$K; p <- job$p
  pr <- partition_rows(length(fm$subs), SEED_B2 + p, N_HALF_A)
  M <- length(job$sets)
  rA <- matrix(NA_real_, M, 2L); rB <- matrix(NA_real_, M, 2L)
  for (i in seq_len(M)) {
    a <- run_set_rows(fm, job$sets[[i]], pr$A, N_SPLIT_B2,
                      seed_half_b2(1L, K, p, i))
    b <- run_set_rows(fm, job$sets[[i]], pr$B, N_SPLIT_B2,
                      seed_half_b2(2L, K, p, i))
    rA[i, ] <- c(a$loss_avoid_adj$r_sb, a$ev_sens_agn$r_sb)
    rB[i, ] <- c(b$loss_avoid_adj$r_sb, b$ev_sens_agn$r_sb)
  }
  rows <- list()
  for (ti in seq_along(TARGETS)) {
    tg <- TARGETS[ti]
    for (dir in c("A_selects_B_evaluates", "B_selects_A_evaluates")) {
      sel <- if (dir == "A_selects_B_evaluates") rA[, ti] else rB[, ti]
      evl <- if (dir == "A_selects_B_evaluates") rB[, ti] else rA[, ti]
      g <- gain_stats(sel, evl, N_PERM_B2,
                      seed_perm_b2(K, p) + 1000000L * ti +
                        500000L * (dir == "B_selects_A_evaluates"))
      rows[[length(rows) + 1L]] <- tibble::tibble(
        K = K, partition = p, target = tg, direction = dir,
        n_subjects_selection = if (dir == "A_selects_B_evaluates")
          length(pr$A) else length(pr$B),
        n_subjects_evaluation = if (dir == "A_selects_B_evaluates")
          length(pr$B) else length(pr$A),
        n_sets = M, n_sets_used = g$n_used, n_top_decile = g$n_top,
        in_sample_single_best_gain = g$in_single,
        out_sample_single_best_gain = g$out_single,
        in_sample_decile_gain = g$in_decile,
        out_sample_decile_gain = g$out_decile,
        optimism_single = g$in_single - g$out_single,
        optimism_decile = g$in_decile - g$out_decile,
        shrinkage_single = g$out_single / g$in_single,
        shrinkage_decile = g$out_decile / g$in_decile,
        spearman_sel_vs_eval = g$rho,
        null_out_sample_single_best_gain = g$null_single,
        null_out_sample_decile_gain = g$null_decile)
    }
  }
  list(K = K, p = p, gains = dplyr::bind_rows(rows),
       per_set = tibble::tibble(
         K = K, partition = p, set_id = seq_len(M),
         r_sb_A_loss_avoid_adj = rA[, 1L], r_sb_B_loss_avoid_adj = rB[, 1L],
         r_sb_A_ev_sens_agn = rA[, 2L], r_sb_B_ev_sens_agn = rB[, 2L]))
})
say("B2: %d partition-cells in %.1f min\n", length(jobs),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

b2_selection_gain <- purrr::map_dfr(res, function(z) z$gains) |>
  dplyr::arrange(target, K, partition, direction)
save_csv(b2_selection_gain, "b2_selection_gain")

per_set_all <- purrr::map_dfr(res, function(z) z$per_set)
save_csv(per_set_all |> dplyr::filter(partition %in% EXPORT_PARTS) |>
           dplyr::arrange(K, partition, set_id),
         "b2_export_per_set_halves")

# =============================================================================
# 5. Per-K summary with the 90% interval (the G-B2b statistic)
# =============================================================================
qi <- function(x, p) if (all(is.na(x))) NA_real_ else
  unname(stats::quantile(x, p, na.rm = TRUE))
b2_summary_by_K <- b2_selection_gain |>
  dplyr::group_by(target, K) |>
  dplyr::summarise(
    n_partition_directions = dplyr::n(),
    mean_n_sets_used = mean(n_sets_used, na.rm = TRUE),
    mean_in_sample_single = mean(in_sample_single_best_gain, na.rm = TRUE),
    mean_out_sample_single = mean(out_sample_single_best_gain, na.rm = TRUE),
    single_lo90 = qi(out_sample_single_best_gain, .05),
    single_hi90 = qi(out_sample_single_best_gain, .95),
    mean_in_sample_decile = mean(in_sample_decile_gain, na.rm = TRUE),
    mean_out_sample_decile = mean(out_sample_decile_gain, na.rm = TRUE),
    decile_lo90 = qi(out_sample_decile_gain, .05),
    decile_hi90 = qi(out_sample_decile_gain, .95),
    mean_optimism_single = mean(optimism_single, na.rm = TRUE),
    mean_optimism_decile = mean(optimism_decile, na.rm = TRUE),
    mean_shrinkage_single = mean(out_sample_single_best_gain, na.rm = TRUE) /
      mean(in_sample_single_best_gain, na.rm = TRUE),
    mean_shrinkage_decile = mean(out_sample_decile_gain, na.rm = TRUE) /
      mean(in_sample_decile_gain, na.rm = TRUE),
    mean_spearman_sel_vs_eval = mean(spearman_sel_vs_eval, na.rm = TRUE),
    mean_null_single = mean(null_out_sample_single_best_gain, na.rm = TRUE),
    null_single_lo90 = qi(null_out_sample_single_best_gain, .05),
    null_single_hi90 = qi(null_out_sample_single_best_gain, .95),
    mean_null_decile = mean(null_out_sample_decile_gain, na.rm = TRUE),
    null_decile_lo90 = qi(null_out_sample_decile_gain, .05),
    null_decile_hi90 = qi(null_out_sample_decile_gain, .95),
    mc_se_out_decile = stats::sd(out_sample_decile_gain, na.rm = TRUE) /
      sqrt(sum(is.finite(out_sample_decile_gain))),
    .groups = "drop") |>
  dplyr::left_join(naive_ref, by = c("K", "target")) |>
  dplyr::arrange(target, K)
print(as.data.frame(b2_summary_by_K |>
                      dplyr::filter(target == GATE_TARGET) |>
                      dplyr::select(K, mean_in_sample_single,
                                    mean_out_sample_single,
                                    mean_out_sample_decile, decile_lo90,
                                    decile_hi90, naive_best_minus_median)))
save_csv(b2_summary_by_K, "b2_summary_by_K")

# =============================================================================
# 6. HG-ceiling guard
# =============================================================================
# Strict (per-set, E8 rule) on a committed subsample of the PRIMARY cell's
# first partition, where the frozen reliability_point supplies a CI on each
# half; and the per-cell central reading over every partition. Both are
# reported; see the equivalent note in 30_partition_stability.R for why per-set
# exceedance is expected under H1 and is not by itself evidence of leakage.
say("\n=== HG-ceiling (B2) ===\n")
pr1 <- partition_rows(length(fm$subs), seed_part_b2(1L), N_HALF_A)
ceil_jobs <- lapply(seq_len(N_CEIL_SETS), function(i)
  list(i = i, cols = sets_by_K[[as.character(GATE_K)]][[i]]))
parallel::clusterExport(cl, c("pr1", "GATE_K", "N_BOOT_CEIL", "seed_ceil_b2",
                              "seed_half_b2"), envir = environment())
ceil_res <- dplyr::bind_rows(parallel::parLapply(cl, ceil_jobs, function(job) {
  a <- run_set_rows(fm, job$cols, pr1$A, N_SPLIT_B2,
                    seed_half_b2(1L, GATE_K, 1L, job$i), want_est = TRUE)
  b <- run_set_rows(fm, job$cols, pr1$B, N_SPLIT_B2,
                    seed_half_b2(2L, GATE_K, 1L, job$i), want_est = TRUE)
  purrr::map_dfr(TARGETS, function(tg) purrr::map_dfr(c("A", "B"), function(h) {
    E <- if (h == "A") a$E[[tg]] else b$E[[tg]]
    rp <- reliability_point(E$h1, E$h2, N_BOOT_CEIL,
                            seed_ceil_b2(if (h == "A") 1L else 2L, job$i),
                            MIN_PAIRS)                                  # FROZEN
    tibble::tibble(set_id = job$i, target = tg, subject_half = h,
                   K = GATE_K, r_sb = rp$r_sb, ci_lo_sb = rp$ci_lo_sb,
                   ci_hi_sb = rp$ci_hi_sb, n_subjects = rp$n_subjects)
  }))
}))
parallel::stopCluster(cl)

# PATH AGREEMENT. The ceiling pass recomputes partition 1 at K = 12
# independently of the production run, and scores it through the frozen
# reliability_point rather than the accelerator's own sb_correct(mean(rs)).
# The two must agree bit-for-bit on the primary target.
prod1 <- per_set_all |>
  dplyr::filter(K == GATE_K, partition == 1L, set_id <= N_CEIL_SETS)
path_d <- vapply(seq_len(nrow(ceil_res)), function(i) {
  r <- ceil_res[i, ]
  col <- paste0("r_sb_", r$subject_half, "_", r$target)
  v <- prod1[[col]][match(r$set_id, prod1$set_id)]
  if (is.na(v) && is.na(r$r_sb)) 0 else abs(v - r$r_sb)
}, numeric(1))
say("path agreement (B2): committed r_sb vs frozen reliability_point recompute, max |d| = %.3g over %d set-halves\n",
    max(path_d, na.rm = TRUE), length(path_d))
stopifnot(all(is.na(path_d) | path_d < 1e-12))

ceil_res <- ceil_res |>
  dplyr::mutate(ceiling_K = sb_k(K, unname(R1[target])),
                exceeds_point = r_sb > ceiling_K,
                exceeds_beyond_noise = ci_lo_sb > ceiling_K)
hg_set_any <- any(ceil_res$exceeds_beyond_noise, na.rm = TRUE)

hg_cell <- per_set_all |>
  tidyr::pivot_longer(dplyr::starts_with("r_sb_"),
                      names_to = "key", values_to = "r_sb") |>
  dplyr::mutate(subject_half = sub("^r_sb_([AB])_.*$", "\\1", key),
                target = sub("^r_sb_[AB]_", "", key)) |>
  dplyr::group_by(target, K, subject_half) |>
  dplyr::summarise(n_values = sum(is.finite(r_sb)),
                   median_r_sb = stats::median(r_sb, na.rm = TRUE),
                   frac_point_exceed = mean(r_sb > sb_k(K, unname(R1[target])),
                                            na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::mutate(ceiling_K = sb_k(K, unname(R1[target])),
                median_exceeds_point = median_r_sb > ceiling_K)
print(as.data.frame(hg_cell))
say("HG-ceiling (B2) per-set reading: %d of %d bootstrapped set-halves CI-clear the ceiling -> %s\n",
    sum(ceil_res$exceeds_beyond_noise, na.rm = TRUE), nrow(ceil_res),
    if (hg_set_any) "FLAGGED" else "clear")
save_csv(ceil_res, "b2_hg_ceiling_sets")
save_csv(hg_cell, "b2_hg_ceiling")

# =============================================================================
# 7. Gates G-B2a..c (evaluated at the primary cell: K = 12, loss_avoid_adj)
# =============================================================================
gk <- b2_summary_by_K |>
  dplyr::filter(K == GATE_K, target == GATE_TARGET)
stopifnot(nrow(gk) == 1L)

g_b2a <- isTRUE(gk$mean_out_sample_single < gk$mean_in_sample_single)
g_b2a_decile <- isTRUE(gk$mean_out_sample_decile < gk$mean_in_sample_decile)

g_b2b_lo <- gk$decile_lo90
g_b2b_hi <- gk$decile_hi90
g_b2b_positive <- isTRUE(g_b2b_lo > 0)
g_b2b_reading <- if (g_b2b_positive)
  paste("a modest, transferable composition advantage from direct empirical",
        "selection at build scale") else
  paste("selection does not reliably beat the average out of sample at these N;",
        "the best-minus-median gap is winner's curse")

null_cov0 <- isTRUE(gk$null_decile_lo90 <= 0 && gk$null_decile_hi90 >= 0)
null_small <- isTRUE(abs(gk$mean_null_decile) < NULL_ABS_TOL)
g_b2c <- null_cov0 && null_small

say("\n=== B2 GATES (Set %d, K %d, %s) ===\n", GATE_FORM, GATE_K, GATE_TARGET)
say("naive E3 in-sample best-minus-median (full 63) = %.6f\n",
    gk$naive_best_minus_median)
say("G-B2a optimism exists: out-of-sample single-best gain %.6f < in-sample %.6f -> %s (decile version: %s)\n",
    gk$mean_out_sample_single, gk$mean_in_sample_single,
    if (g_b2a) "PASS" else "FAIL", if (g_b2a_decile) "also lower" else "NOT lower")
say("G-B2b residual selectability: top-decile out-of-sample gain %.6f, 90%% interval [%.6f, %.6f] -> %s\n",
    gk$mean_out_sample_decile, g_b2b_lo, g_b2b_hi,
    if (g_b2b_positive) "lower bound > 0" else "interval covers 0")
say("       reading: %s\n", g_b2b_reading)
say("G-B2c null control: permuted top-decile out-of-sample gain %.6f, 90%% interval [%.6f, %.6f] -> %s\n",
    gk$mean_null_decile, gk$null_decile_lo90, gk$null_decile_hi90,
    if (g_b2c) "PASS" else "FAIL")
say("shrinkage: single-best %.4f, top-decile %.4f\n",
    gk$mean_shrinkage_single, gk$mean_shrinkage_decile)

part <- list(
  gate = list(
    probe = "B2 optimism / winner's curse",
    primary_cell = list(form = GATE_FORM, K = GATE_K, target = GATE_TARGET),
    naive_reference = list(
      e3_best_minus_median_full63 = gk$naive_best_minus_median,
      e3_r_sb_median = gk$naive_r_sb_median, e3_r_sb_max = gk$naive_r_sb_max,
      source = "elevation/results/e3/heterogeneity_summary.csv (read-only)"),
    `G-B2a` = list(
      criterion = "out-of-sample gain < in-sample apparent gain",
      in_sample_single_best_gain = gk$mean_in_sample_single,
      out_sample_single_best_gain = gk$mean_out_sample_single,
      optimism_single = gk$mean_optimism_single,
      shrinkage_single = gk$mean_shrinkage_single,
      in_sample_decile_gain = gk$mean_in_sample_decile,
      out_sample_decile_gain = gk$mean_out_sample_decile,
      optimism_decile = gk$mean_optimism_decile,
      shrinkage_decile = gk$mean_shrinkage_decile,
      decile_version_also_lower = g_b2a_decile,
      verdict = if (g_b2a) "PASS" else "FAIL", pass = g_b2a),
    `G-B2b` = list(
      rule = paste("interval-vs-zero rule fixed in the spec before results:",
                   "lower bound > 0 reports a modest transferable advantage,",
                   "an interval covering 0 reports winner's curse; EITHER",
                   "outcome is honest"),
      n_partition_directions = gk$n_partition_directions,
      top_decile_out_of_sample_gain = gk$mean_out_sample_decile,
      interval_90 = list(lo = g_b2b_lo, hi = g_b2b_hi),
      mc_se = gk$mc_se_out_decile,
      lower_bound_positive = g_b2b_positive,
      reading = g_b2b_reading),
    `G-B2c` = list(
      permuted_top_decile_gain = gk$mean_null_decile,
      interval_90 = list(lo = gk$null_decile_lo90, hi = gk$null_decile_hi90),
      permuted_single_best_gain = gk$mean_null_single,
      criterion = sprintf("90%% interval covers 0 AND |mean| < %.3f",
                          NULL_ABS_TOL),
      interval_covers_zero = null_cov0, mean_below_tolerance = null_small,
      verdict = if (g_b2c) "PASS" else "FAIL", pass = g_b2c),
    HG_ceiling = list(
      comparator = "E0 full-data ceiling SB-projected in problem count to K",
      ceiling_r1 = as.list(R1),
      per_set_reading = list(
        n_bootstrapped_set_halves = nrow(ceil_res),
        n_ci_supported_exceedance = sum(ceil_res$exceeds_beyond_noise,
                                        na.rm = TRUE),
        frac_point_exceedance = mean(ceil_res$exceeds_point, na.rm = TRUE),
        flagged = hg_set_any,
        caveat = paste("individual-set exceedance is EXPECTED under H1 and on",
                       "half samples, which carry more estimation noise than",
                       "the full 63")),
      per_cell_reading = lapply(seq_len(nrow(hg_cell)), function(i)
        as.list(hg_cell[i, ])),
      disposition = "reported; not worked around")),
  config = list(
    probe = "B2", seed_base = SEED_B2,
    regime = "experience (Trial 6-25)", r_reps = R_REPS,
    K_grid = K_GRID_B2, targets = list(primary = TARGETS[1],
                                       secondary = TARGETS[2]),
    M = M_CAP, b_part_b2 = B_PART_B2, n_split_b2 = N_SPLIT_B2,
    decile = DECILE, n_perm = N_PERM_B2, min_sets = MIN_SETS_B2,
    null_abs_tol = NULL_ABS_TOL, n_ceiling_sets = N_CEIL_SETS,
    n_boot_ceiling = N_BOOT_CEIL, half_sizes = c(N_HALF_A, 63L - N_HALF_A),
    item_sets = paste("the SAME committed E3 Set-1 sets (parsed from",
                      "elevation/results/e3/reliability_by_set_real.csv), so",
                      "the out-of-sample gains are directly comparable with",
                      "E3's in-sample best-minus-median on the full 63"),
    both_directions = TRUE,
    constants_fixed_by_code = paste(
      "The spec fixes B_PART_B2, N_SPLIT_B2 and the decile. N_PERM_B2 = 20,",
      "MIN_SETS_B2 = 30, NULL_ABS_TOL = 0.01, N_CEIL_SETS = 50,",
      "N_BOOT_CEIL = 1000 and EXPORT_PARTS = 1:5 were fixed in this script",
      "before any result existed and are recorded here."),
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
    e4_h2_gate_rerun_or_amended = FALSE,
    runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
jsonlite::write_json(part, file.path(res_dir, "dispb_part_b2.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

seeds_b2 <- dplyr::bind_rows(
  tibble::tibble(probe = "B2",
                 scope = sprintf("partition%d", seq_len(B_PART_B2)),
                 purpose = "subject_partition",
                 seed = seed_part_b2(seq_len(B_PART_B2)),
                 derivation = "SEED_B2 + partition_index"),
  purrr::map_dfr(K_GRID_B2, function(K) tibble::tibble(
    probe = "B2", scope = sprintf("K%d", K),
    purpose = c("halving_half_A_base", "halving_half_B_base",
                "evaluation_label_permutation_base"),
    seed = c(seed_half_b2(1L, K, 0L, 0L), seed_half_b2(2L, K, 0L, 0L),
             seed_perm_b2(K, 0L)),
    derivation = c(
      "SEED_B2 + 200000000*K_index; per unit = base + 100000*partition + set_id",
      "SEED_B2 + 900000000 + 200000000*K_index; per unit = base + 100000*partition + set_id",
      "SEED_B2 + 800000 + 1000*K_index; per unit = base + partition (+1000000*target_index +500000 for the reverse direction)"))),
  tibble::tibble(probe = "B2", scope = "ceiling",
                 purpose = c("reliability_point_half_A", "reliability_point_half_B"),
                 seed = c(seed_ceil_b2(1L, 0L), seed_ceil_b2(2L, 0L)),
                 derivation = "SEED_B2 + 600000 + 100000*half; per set = base + set_id"))
readr::write_csv(seeds_b2, file.path(res_dir, "dispb_seeds_b2.csv"))
stopifnot(!anyDuplicated(seeds_b2$seed[is.finite(seeds_b2$seed)]))
half_enum <- unlist(lapply(c(1L, 2L), function(h) unlist(lapply(K_GRID_B2,
  function(K) seed_half_b2(h, K, c(1L, B_PART_B2), c(1L, M_CAP))))))
stopifnot(!anyDuplicated(half_enum), max(half_enum) < .Machine$integer.max)
say("seed hygiene: %d enumerated seeds unique; halving-seed grid max %d\n",
    nrow(seeds_b2), max(half_enum))

say("\nB2 complete in %.1f min. Outputs in elevation/results/dispb/\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))

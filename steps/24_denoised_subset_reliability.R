# Step 24: high-precision outcome reliability of the problem subsets
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e1/, e2/, e3/ (steps 20-22)
# Writes:  output/elevation/results/e3b/ (outcome_reliability_real.csv, simulation tables, set identity, seeds)
# Paper:   Study 2 (The information mechanism)
# Runtime: about 6 hours with 12 parallel workers

# E3b: denoised / subject-split reliability outcome.
#
# EXPLORATORY. The frozen pre-specification and the frozen H2 gate result
# (step 25, FAIL) stand as the confirmatory record. Nothing here re-runs,
# re-scores, amends or replaces the frozen H2 gate. This script only BUILDS a
# lower-noise, honestly-evaluable reliability outcome for an exploratory
# re-test.
#
# It fits NO mechanism model, performs NO optimization, constructs NO candidate
# design, and reads NO freeze or test outcome.
#
# BUILD-ONLY. Part A uses only the 63 Set-1 BUILD subjects of the frozen
# participant split (step 20), partitioned into disjoint halves. The 137 FREEZE
# and 205 TEST subjects are NOT read. The guarantee is STRUCTURAL, as in steps
# 21-22: the full-sample frame has every choice / outcome column DROPPED at
# load, so only the build-restricted frame can carry a choice. Part B is
# SIMULATION under the frozen build generative model and is labeled as such
# everywhere it appears.
#
# ENGINE PROVENANCE. The frozen shared engines are SOURCED unchanged and none
# is redefined:
#   R/step2_reliability.R  build_subject_data, fit_loss_baselines,
#                          estimate_half, split_correlations, sb_correct,
#                          reliability_point
#   R/cpc18_features.R     load_cpc18, is_clean_row, clean_problem_features,
#                          add_scaled_features
#   R/simulate_recovery.R  draw_true_params, analytic_p_b, simulate_choices,
#                          boot_cor
# reliability_point (hence split_correlations and sb_correct, which it calls)
# produces EVERY reported per-set reliability in Part A. boot_cor produces every
# reported outcome-reliability CI. draw_true_params, simulate_choices and
# analytic_p_b produce every simulated quantity in Part B.
#
# ACCELERATOR (same device as step 22, not a redefinition). estimate_half
# always runs the 4-parameter Firth fit that yields b_ev / b_risk / b_loss,
# which are NOT E3b targets. Evaluating Part A through it directly would need
# ~2.0e8 calls (~140 h single-threaded at the 1.9-3.6 ms per call measured in
# step 22). `half_targets()` below evaluates ONLY the two targets,
# using the SAME formulas as the frozen estimate_half, vectorised over subjects.
# It is not a substitute definition: gate G_EQ re-runs the frozen estimate_half
# and asserts agreement to < 1e-12 with identical NA patterns BEFORE any result
# is produced, and the D1 spot-check re-derives whole sets end-to-end through the
# frozen engine on both subject halves. The definitions are byte-equivalent to
# E3's; the equivalence is re-proven here rather than inherited.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})

source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))
source(here("R", "simulate_recovery.R"))

res_dir <- here("output", "elevation", "results", "e3b")
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
# 0. Frozen constants (fixed BEFORE any result exists)
# =============================================================================
RAW_MD5     <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_E3B    <- 20260721L
REGIME      <- "exp"
POS         <- 6:25              # experience regime, r = 20 (the full pool)
R_REPS      <- length(POS)
K_GRID      <- c(6L, 8L, 10L, 12L)
TARGETS     <- c("loss_avoid_adj", "ev_sens_agn")   # primary, secondary
GATE_FORM   <- 1L                # Part A is Set 1 only
MIN_PAIRS   <- 10L               # frozen split_correlations default

# Part A
N_SPLIT_HALF <- 200L             # halvings of S per subject half
N_HALF_A     <- 32L              # largest-remainder halves of 63
N_HALF_B     <- 31L

# Part B (SIMULATION)
N_SIM_SUBJ  <- 2000L             # vs E3's 300
N_SPLIT_SIM <- 100L              # vs E3's 50
N_REPS_SIM  <- 20L
N_SIM_SUBJ_E3 <- 300L            # E3's setting, replayed in the noise study

# Further constants, fixed HERE before any result exists and recorded in
# e3b_config.json. N_BOOT_SET matches E3's noise-floor N_BOOT.
N_BOOT_SET  <- 1000L             # frozen reliability_point bootstrap, per set
N_BOOT_COR  <- 2000L             # frozen boot_cor default, outcome reliability
N_GEQ       <- 8000L             # G_EQ halves checked against estimate_half
N_SPOT      <- 12L               # sets re-derived end-to-end through the engine
SPOT_TOL_EV <- 5e-3              # ev_sens_agn rank-tie tolerance (as in step 22)
N_NOISE_SETS <- 50L              # sets in the Part B replicate noise study
N_NOISE_REPS <- 20L              # independent replicate simulations per set
N_WORKERS   <- 12L

# ---- seed derivations -------------------------------------------------------
# The BASE seeds are pre-declared. The per-set derivations below are
# deterministic functions of those, so every draw reproduces exactly under any
# execution order. Blocks are spaced so that the A stream and the B stream can
# never collide (asserted below), so A and B share no randomness.
SEED_HALF_A <- SEED_E3B + 1L
SEED_HALF_B <- SEED_E3B + 2L
seed_setA <- function(K, i) SEED_HALF_A + 1000000L + 10000L * K + i
seed_setB <- function(K, i) SEED_HALF_B + 9000000L + 10000L * K + i
seed_boot <- function(K, i, h) SEED_E3B + 20000000L + 100000L * h + 10000L * K + i
seed_simpar <- function(K) SEED_E3B + 100L + K
seed_simset <- function(K, i) SEED_E3B + 200L + 1000000L + 100000L * K + i
seed_noise  <- function(K, i, j, n) SEED_E3B + 40000000L + 4000000L * (n == N_SIM_SUBJ) +
  200000L * K + 100L * i + j

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
stopifnot(n_build == 344L, n_freeze == 137L, n_test == 205L)      # frozen split sizes
stopifnot(!anyDuplicated(build_ids))
stopifnot(all(split_tbl$held_out_until_E7 == (split_tbl$split == "test")))
say("split: build %d / freeze %d (E6, not read) / test %d (sealed to E7)\n",
    n_build, n_freeze, n_test)

raw <- load_cpc18(raw_path)
stopifnot(all(cpc18_expected_columns %in% names(raw)))

# Structural build-only guarantee (as in steps 21-22): the full-sample frame keeps no
# choice / outcome column, so no freeze or test choice can reach any computation.
choice_cols <- c("B", "Button", "Payoff", "Forgone", "RT", "Apay", "Bpay",
                 "Feedback")
design_all <- raw |> dplyr::select(-dplyr::any_of(choice_cols))
stopifnot(!any(choice_cols %in% names(design_all)))

build_raw <- raw |> dplyr::filter(SubjID %in% build_ids)
rm(raw); invisible(gc())
stopifnot(setequal(unique(build_raw$SubjID), build_ids))
stopifnot(!any(build_raw$SubjID %in% held_out))                   # build-only

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

# Frozen build-only loss baselines and per-subject data structures (as in steps 21-22).
baselines <- fit_loss_baselines(build_trials)
subj_data <- build_subject_data(build_trials, baselines)
stopifnot(length(subj_data) == 344L)
say("frozen loss baselines (exp): b0=%.6f b1=%.6f\n",
    baselines$exp[1], baselines$exp[2])

feat_all <- feat_all |>
  dplyr::mutate(
    d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a),
    p_base_exp  = ifelse(loss_informative,
                         stats::plogis(baselines$exp[1] +
                                         baselines$exp[2] * d_ev_noloss), NA))
chk <- subj_data[[1]]$feat
stopifnot(all(abs(chk$p_base_exp -
                    feat_all$p_base_exp[match(chk$GameID, feat_all$GameID)]) <
                1e-12, na.rm = TRUE))

# =============================================================================
# 2. Set 1 structure and the FIXED subject partition
# =============================================================================
s1_trials <- build_trials |> dplyr::filter(Set == GATE_FORM)
s1_subs   <- as.character(sort(unique(s1_trials$SubjID)))
s1_probs  <- sort(unique(s1_trials$GameID))
stopifnot(length(s1_subs) == 63L)                                  # step-20 split
say("Set 1 build subjects: %d; clean problems: %d\n",
    length(s1_subs), length(s1_probs))

# NB[subject, problem] = number of B choices over the r = 20 experience trials.
# POS is FIXED at the full experience pool, so NB is a constant and every half
# is a lookup. NA marks a problem that subject never saw.
NB <- matrix(NA_real_, length(s1_subs), length(s1_probs),
             dimnames = list(s1_subs, as.character(s1_probs)))
for (sid in s1_subs) {
  sd <- subj_data[[sid]]
  j <- match(sd$feat$GameID, s1_probs)
  ok <- !is.na(j)
  NB[sid, j[ok]] <- rowSums(sd$B_mat[which(ok), POS, drop = FALSE])
}
# Set 1 satisfies the E3 common-set premise exactly (one distinct clean set).
stopifnot(!anyNA(NB))
fi <- feat_all[match(s1_probs, feat_all$GameID), ]
fm <- list(probs = s1_probs, subs = s1_subs, NB = NB, feat = fi,
           li = fi$loss_informative, has_loss_b = fi$has_loss_b,
           d_ev = fi$d_ev, p_base = fi$p_base_exp)

# Deterministic partition, FIXED across all sets and both targets, so the halves
# are a property of SUBJECTS and not of sets.
set.seed(SEED_E3B)
perm  <- sample.int(length(s1_subs))
rowsA <- sort(perm[seq_len(N_HALF_A)])
rowsB <- sort(perm[(N_HALF_A + 1L):length(perm)])
stopifnot(length(rowsA) == N_HALF_A, length(rowsB) == N_HALF_B)
stopifnot(length(intersect(rowsA, rowsB)) == 0L)
stopifnot(setequal(c(rowsA, rowsB), seq_along(s1_subs)))
subs_A <- s1_subs[rowsA]; subs_B <- s1_subs[rowsB]
say("subject partition (seed %d): half A n=%d, half B n=%d, disjoint\n",
    SEED_E3B, length(rowsA), length(rowsB))

# The A and B random streams must not intersect (they share no randomness).
sA <- unlist(lapply(K_GRID, function(K) seed_setA(K, seq_len(2000L))))
sB <- unlist(lapply(K_GRID, function(K) seed_setB(K, seq_len(2000L))))
stopifnot(!anyDuplicated(sA), !anyDuplicated(sB),
          length(intersect(sA, sB)) == 0L)
rm(sA, sB)

# =============================================================================
# 3. The two targets, vectorised over subjects (formulas identical to the frozen
#    estimate_half; equivalence asserted at gate G_EQ below). Same definitions
#    as E3's accelerator.
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

# =============================================================================
# 4. GATE G_EQ: the vectorised path vs the FROZEN estimate_half
#    Re-proven here on E3b's own subject halves rather than inherited from E3.
# =============================================================================
say("\n=== G_EQ: vectorised targets vs frozen estimate_half (%d halves) ===\n",
    N_GEQ)
set.seed(SEED_E3B + 77L)
eq_n <- 0L; eq_max_ev <- 0; eq_max_loss <- 0; eq_na <- 0L
for (it in seq_len(N_GEQ)) {
  rows <- if (it %% 2L == 0L) rowsA else rowsB
  ri  <- rows[sample.int(length(rows), 1)]
  sid <- s1_subs[ri]
  sd  <- subj_data[[sid]]
  kh  <- sample(3:6, 1)
  cols <- sample(seq_along(fm$probs), kh)
  a <- estimate_half(sd, match(fm$probs[cols], sd$feat$GameID), POS, REGIME)
  b <- half_targets(fm$NB[ri, cols, drop = FALSE], fm, cols, R_REPS)
  for (nm in TARGETS) {
    x <- unname(a[[nm]]); y <- unname(b[[nm]])
    if (is.na(x) != is.na(y)) { eq_na <- eq_na + 1L; next }
    if (!is.na(x)) {
      d <- abs(x - y)
      if (nm == "ev_sens_agn") eq_max_ev <- max(eq_max_ev, d)
      else eq_max_loss <- max(eq_max_loss, d)
    }
  }
  eq_n <- eq_n + 1L
}
say("G_EQ: %d halves; max|d ev_sens_agn| = %.3g; max|d loss_avoid_adj| = %.3g; NA mismatches = %d\n",
    eq_n, eq_max_ev, eq_max_loss, eq_na)
stopifnot(eq_na == 0L, eq_max_ev < 1e-12, eq_max_loss < 1e-12)

# =============================================================================
# 5. Part A: real Set-1 subject-split reliabilities
# =============================================================================
# The item sets are REUSED from step 22's saved per-set file, parsed rather than
# redrawn, so set identity (D3) holds by construction. Section 10 additionally
# re-derives them from E3's seeds alone and asserts agreement.
e3_real <- readr::read_csv(
  here("output", "elevation", "results", "e3", "reliability_by_set_real.csv"),
  show_col_types = FALSE)
e3_s1 <- e3_real |> dplyr::filter(form == GATE_FORM)
stopifnot(nrow(e3_s1) == 2L * length(K_GRID) * 2000L)

# game_ids are identical across the two targets for a given (K, set_id);
# assert that, then key the set list on (K, set_id).
sets_real <- e3_s1 |>
  dplyr::filter(target == TARGETS[1]) |>
  dplyr::select(K, set_id, game_ids, n_subjects_e3 = n_subjects) |>
  dplyr::arrange(K, set_id)
chk_alt <- e3_s1 |>
  dplyr::filter(target == TARGETS[2]) |>
  dplyr::select(K, set_id, game_ids) |>
  dplyr::arrange(K, set_id)
stopifnot(identical(sets_real$game_ids, chk_alt$game_ids))       # D3
stopifnot(all(sets_real$n_subjects_e3 == 63L))

parse_ids <- function(s) as.integer(strsplit(s, ";", fixed = TRUE)[[1]])
sets_real_cols <- lapply(sets_real$game_ids, function(s) {
  g <- parse_ids(s)
  cc <- match(g, fm$probs)
  stopifnot(!anyNA(cc), !anyDuplicated(cc))
  cc
})
stopifnot(all(lengths(sets_real_cols) == sets_real$K))
say("Part A: %d committed E3 Set-1 sets reused (K = %s)\n",
    nrow(sets_real), paste(K_GRID, collapse = ", "))

# One set, one subject half: h1/h2 estimate matrices over N_SPLIT_HALF halvings of S,
# scored by the FROZEN reliability_point (which calls the frozen
# split_correlations and sb_correct).
run_half_sets <- function(fm, cols, rows, n_split, seed) {
  NBs <- fm$NB[rows, cols, drop = FALSE]
  ns <- length(rows); K <- length(cols); kh <- K %/% 2L
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
  E
}

say("\n=== Part A: subject-split reliabilities (frozen reliability_point, N_boot = %d) ===\n",
    N_BOOT_SET)
cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(tidyverse); library(here) })
  source(here("R", "step2_reliability.R"))
  source(here("R", "simulate_recovery.R"))
  NULL
}))
parallel::clusterExport(cl, c("fm", "half_targets", "ev_ranks", "run_half_sets",
                              "rowsA", "rowsB", "TARGETS", "R_REPS", "MIN_PAIRS",
                              "N_SPLIT_HALF", "N_BOOT_SET"),
                        envir = environment())

# Seeds are resolved HERE and carried in the job list, so the stream is
# identical under any execution order.
jobsA <- lapply(seq_len(nrow(sets_real)), function(i) {
  K <- sets_real$K[i]; sid <- sets_real$set_id[i]
  list(row = i, K = K, set_id = sid, cols = sets_real_cols[[i]],
       seed_A = seed_setA(K, sid), seed_B = seed_setB(K, sid),
       boot_A = seed_boot(K, sid, 1L), boot_B = seed_boot(K, sid, 2L))
})
t0 <- Sys.time()
resA <- parallel::parLapply(cl, jobsA, function(job) {
  EA <- run_half_sets(fm, job$cols, rowsA, N_SPLIT_HALF, job$seed_A)
  EB <- run_half_sets(fm, job$cols, rowsB, N_SPLIT_HALF, job$seed_B)
  purrr::map_dfr(TARGETS, function(tg) {
    pA <- reliability_point(EA[[tg]]$h1, EA[[tg]]$h2, N_BOOT_SET,
                            job$boot_A, MIN_PAIRS)            # FROZEN
    pB <- reliability_point(EB[[tg]]$h1, EB[[tg]]$h2, N_BOOT_SET,
                            job$boot_B, MIN_PAIRS)            # FROZEN
    tibble::tibble(
      row = job$row, K = job$K, set_id = job$set_id, target = tg,
      r_sb_A = pA$r_sb, ci_lo_sb_A = pA$ci_lo_sb, ci_hi_sb_A = pA$ci_hi_sb,
      r_sb_B = pB$r_sb, ci_lo_sb_B = pB$ci_lo_sb, ci_hi_sb_B = pB$ci_hi_sb,
      r_raw_A = pA$r_raw, r_raw_B = pB$r_raw,
      n_subj_A = pA$n_subjects, n_subj_B = pB$n_subjects,
      cov_A = pA$n_splits_used / pA$n_splits,
      cov_B = pB$n_splits_used / pB$n_splits,
      n_splits_used_A = pA$n_splits_used, n_splits_used_B = pB$n_splits_used)
  })
})
elA <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
say("Part A: %d set-targets scored in %.1f min\n",
    length(resA) * length(TARGETS), elA)

reliability_subject_split_real <- dplyr::bind_rows(resA) |>
  dplyr::left_join(sets_real |> dplyr::mutate(row = dplyr::row_number()) |>
                     dplyr::select(row, game_ids), by = "row") |>
  dplyr::mutate(
    form = GATE_FORM,
    n_loss_informative = vapply(sets_real_cols[row],
                                function(z) sum(fm$li[z]), integer(1))) |>
  dplyr::select(form, K, target, set_id, game_ids, n_loss_informative,
                r_sb_A, r_sb_B, n_subj_A, n_subj_B, cov_A, cov_B,
                r_raw_A, r_raw_B, ci_lo_sb_A, ci_hi_sb_A, ci_lo_sb_B,
                ci_hi_sb_B, n_splits_used_A, n_splits_used_B) |>
  dplyr::arrange(target, K, set_id)
stopifnot(all(reliability_subject_split_real$n_subj_A == N_HALF_A),
          all(reliability_subject_split_real$n_subj_B == N_HALF_B))

# =============================================================================
# 6. D1 spot-check: whole sets re-derived end to end through the FROZEN engine
# =============================================================================
# G_EQ proves per-value equivalence. This closes the loop at the level of a
# whole set's r_sb, on BOTH subject halves, through frozen estimate_half.
# Tolerances differ by target for the reason measured in step 22:
# loss_avoid_adj arithmetic is identical in both paths so every split's Spearman
# agrees EXACTLY (tolerance 0); ev_sens_agn values are built from integer choice
# counts, so subjects hold exactly tied values that split_correlations ranks, and
# a last-bit difference resolves a tie the other way. That is rank-tie
# sensitivity of the FROZEN estimator, not of the vectorised path.
say("\n=== D1 spot-check: frozen estimate_half re-derivation of %d whole sets ===\n",
    N_SPOT)
spot_idx <- which(sets_real$K == 12L)[seq_len(N_SPOT)]
spot_jobs <- lapply(spot_idx, function(i)
  list(set_id = sets_real$set_id[i], cols = sets_real_cols[[i]],
       seed_A = seed_setA(12L, sets_real$set_id[i]),
       seed_B = seed_setB(12L, sets_real$set_id[i])))
parallel::clusterExport(cl, c("subj_data", "s1_subs", "POS", "REGIME"),
                        envir = environment())
spot_detail <- dplyr::bind_rows(parallel::parLapply(cl, spot_jobs, function(job) {
  cols <- job$cols; K <- length(cols); kh <- K %/% 2L
  purrr::map_dfr(c("A", "B"), function(hh) {
    rows <- if (hh == "A") rowsA else rowsB
    seed <- if (hh == "A") job$seed_A else job$seed_B
    subs <- s1_subs[rows]
    Ef <- list()
    for (tg in TARGETS) Ef[[tg]] <- list(
      h1 = matrix(NA_real_, length(subs), N_SPLIT_HALF),
      h2 = matrix(NA_real_, length(subs), N_SPLIT_HALF))
    set.seed(seed)
    ties <- integer(N_SPLIT_HALF)
    for (s in seq_len(N_SPLIT_HALF)) {
      take <- sample.int(K, 2L * kh)
      i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
      for (u in seq_along(subs)) {
        sd <- subj_data[[subs[u]]]
        r1 <- estimate_half(sd, match(fm$probs[cols[i1]], sd$feat$GameID),
                            POS, REGIME)                       # FROZEN
        r2 <- estimate_half(sd, match(fm$probs[cols[i2]], sd$feat$GameID),
                            POS, REGIME)                       # FROZEN
        for (tg in TARGETS) { Ef[[tg]]$h1[u, s] <- r1[[tg]]
                              Ef[[tg]]$h2[u, s] <- r2[[tg]] }
      }
      v <- Ef$ev_sens_agn$h1[, s]; v <- v[is.finite(v)]
      ties[s] <- length(v) - length(unique(v))
    }
    Ev <- run_half_sets(fm, cols, rows, N_SPLIT_HALF, seed)
    purrr::map_dfr(TARGETS, function(tg) {
      rf <- split_correlations(Ef[[tg]]$h1, Ef[[tg]]$h2, MIN_PAIRS)   # FROZEN
      rv <- split_correlations(Ev[[tg]]$h1, Ev[[tg]]$h2, MIN_PAIRS)   # FROZEN
      frozen_sb <- sb_correct(mean(rf, na.rm = TRUE))                 # FROZEN
      vect_sb   <- sb_correct(mean(rv, na.rm = TRUE))                 # FROZEN
      tibble::tibble(
        set_id = job$set_id, subject_half = hh, target = tg,
        frozen_r_sb = frozen_sb, vectorised_r_sb = vect_sb,
        abs_diff = abs(frozen_sb - vect_sb),
        n_splits_differing = sum(abs(rf - rv) > 0, na.rm = TRUE),
        max_split_abs_diff = suppressWarnings(max(abs(rf - rv), na.rm = TRUE)),
        mean_tied_subjects_per_split = mean(ties))
    })
  })
}))
spot_summary <- spot_detail |>
  dplyr::group_by(target, subject_half) |>
  dplyr::summarise(n_sets = dplyr::n(),
                   max_abs_diff = max(abs_diff, na.rm = TRUE),
                   median_abs_diff = stats::median(abs_diff, na.rm = TRUE),
                   mean_splits_differing = mean(n_splits_differing),
                   mean_tied_subjects_per_split =
                     mean(mean_tied_subjects_per_split),
                   .groups = "drop")
print(as.data.frame(spot_summary))
spot_loss <- max(spot_summary$max_abs_diff[
  spot_summary$target == "loss_avoid_adj"], na.rm = TRUE)
spot_ev <- max(spot_summary$max_abs_diff[
  spot_summary$target == "ev_sens_agn"], na.rm = TRUE)
stopifnot(spot_loss == 0)                    # primary target: bit-identical
stopifnot(spot_ev < SPOT_TOL_EV)             # secondary: tie resolution only
say("D1 spot-check PASS: loss_avoid_adj max |d r_sb| = %.3g (exact); ev_sens_agn max = %.3g (< %.0e)\n",
    spot_loss, spot_ev, SPOT_TOL_EV)

# PATH AGREEMENT: the spot-check recomputes each set independently of the
# parallel production run. Assert the production r_sb_A / r_sb_B (produced by the
# frozen reliability_point inside parLapply) equal the spot-check's own frozen
# recomputation, so the saved numbers are tied to the frozen engine and not
# merely to a second copy of the fast path.
prodA <- reliability_subject_split_real |> dplyr::filter(K == 12L)
m <- match(paste(spot_detail$set_id, spot_detail$target),
           paste(prodA$set_id, prodA$target))
stopifnot(!anyNA(m))
prod_val <- ifelse(spot_detail$subject_half == "A", prodA$r_sb_A[m],
                   prodA$r_sb_B[m])
path_d <- abs(spot_detail$vectorised_r_sb - prod_val)
stopifnot(all(is.na(prod_val) == is.na(spot_detail$vectorised_r_sb)))
stopifnot(all(is.na(path_d) | path_d < 1e-12))
say("path agreement: committed r_sb vs independent recompute, max |d| = %.3g over %d set-halves\n",
    suppressWarnings(max(path_d, na.rm = TRUE)), nrow(spot_detail))

# =============================================================================
# 7. Outcome reliability: cross-half Spearman per (K, target)
# =============================================================================
# The honest denominator for a mechanism-model test: a model can at best predict
# half B to the extent half A already agrees with half B. NO disattenuation is
# applied; the SB projection column is auxiliary and labeled.
say("\n=== Outcome reliability (cross-half Spearman, frozen boot_cor) ===\n")
outcome_reliability_real <- purrr::map_dfr(TARGETS, function(tg) {
  purrr::map_dfr(K_GRID, function(Kv) {
    d <- reliability_subject_split_real[
      reliability_subject_split_real$target == tg &
        reliability_subject_split_real$K == Kv, ]
    ok <- is.finite(d$r_sb_A) & is.finite(d$r_sb_B)
    if (sum(ok) < MIN_PAIRS) {
      return(tibble::tibble(
        K = Kv, target = tg, n_sets = nrow(d), n_sets_used = sum(ok),
        outcome_reliability = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        outcome_rel_sb_projected = NA_real_,
        sd_r_sb_A = NA_real_, sd_r_sb_B = NA_real_,
        median_r_sb_A = NA_real_, median_r_sb_B = NA_real_,
        mean_abs_half_diff = NA_real_, n_boot = N_BOOT_COR))
    }
    set.seed(seed_boot(Kv, 0L, 3L))
    bc <- boot_cor(d$r_sb_A, d$r_sb_B, "spearman", N_BOOT_COR)   # FROZEN
    tibble::tibble(
      K = Kv, target = tg, n_sets = nrow(d), n_sets_used = bc$n_used,
      outcome_reliability = bc$r, ci_lo = bc$ci_lo, ci_hi = bc$ci_hi,
      # Auxiliary, clearly labeled: the Spearman-Brown projection of the
      # half-sample agreement to the FULL 63-subject outcome that step 25 modeled.
      # It is a length projection, not a disattenuation, and is never the
      # headline number.
      outcome_rel_sb_projected = sb_correct(bc$r),               # FROZEN
      sd_r_sb_A = stats::sd(d$r_sb_A, na.rm = TRUE),
      sd_r_sb_B = stats::sd(d$r_sb_B, na.rm = TRUE),
      median_r_sb_A = stats::median(d$r_sb_A, na.rm = TRUE),
      median_r_sb_B = stats::median(d$r_sb_B, na.rm = TRUE),
      mean_abs_half_diff = mean(abs(d$r_sb_A - d$r_sb_B), na.rm = TRUE),
      n_boot = N_BOOT_COR)
  })
})

# Context: how each half tracks step 22's FULL 63-subject outcome, which is what
# step 25 actually modeled.
# These correlations are INFLATED by shared subjects (each half is a subset of
# the 63) and are reported only as context, never as an out-of-subject result.
e3_full <- e3_s1 |> dplyr::select(K, target, set_id, r_sb_e3 = r_sb)
half_vs_e3 <- reliability_subject_split_real |>
  dplyr::left_join(e3_full, by = c("K", "target", "set_id")) |>
  dplyr::group_by(K, target) |>
  dplyr::summarise(
    spearman_A_vs_e3_full = {
      ok <- is.finite(r_sb_A) & is.finite(r_sb_e3)
      if (sum(ok) < MIN_PAIRS) NA_real_ else
        suppressWarnings(stats::cor(r_sb_A[ok], r_sb_e3[ok], method = "spearman"))
    },
    spearman_B_vs_e3_full = {
      ok <- is.finite(r_sb_B) & is.finite(r_sb_e3)
      if (sum(ok) < MIN_PAIRS) NA_real_ else
        suppressWarnings(stats::cor(r_sb_B[ok], r_sb_e3[ok], method = "spearman"))
    },
    n_e3_identifiable = sum(is.finite(r_sb_e3)),
    .groups = "drop")
outcome_reliability_real <- outcome_reliability_real |>
  dplyr::left_join(half_vs_e3, by = c("K", "target")) |>
  dplyr::arrange(target, K)
print(as.data.frame(outcome_reliability_real |>
                      dplyr::select(target, K, n_sets_used, outcome_reliability,
                                    ci_lo, ci_hi)))

# =============================================================================
# 8. Part B: SIMULATION high-precision reliabilities
#    Every number in this section is SIMULATED under the frozen generative
#    model (R/simulate_recovery.R) and the frozen uniform-range population (unchanged, deliberately not
#    calibrated to build, so this isolates the noise + precision effect).
#    It is never an empirical result.
# =============================================================================
say("\n=== Part B: SIMULATION high-precision reliabilities (N = %d subjects) ===\n",
    N_SIM_SUBJ)
stopifnot(N_REPS_SIM == R_REPS)
e3_sim <- readr::read_csv(
  here("output", "elevation", "results", "e3", "reliability_recovery_by_set_sim.csv"),
  show_col_types = FALSE)
e3_free <- e3_sim |> dplyr::filter(composition == "free_105")
sets_sim <- e3_free |>
  dplyr::filter(target == TARGETS[1]) |>
  dplyr::select(K, set_id, game_ids) |>
  dplyr::arrange(K, set_id)
chk_sim <- e3_free |>
  dplyr::filter(target == TARGETS[2]) |>
  dplyr::select(K, set_id, game_ids) |>
  dplyr::arrange(K, set_id)
stopifnot(identical(sets_sim$game_ids, chk_sim$game_ids))        # D3
stopifnot(nrow(sets_sim) == 500L * length(K_GRID))

uni_feat <- feat_all
uni_fm <- list(probs = uni_feat$GameID, feat = uni_feat,
               li = uni_feat$loss_informative,
               has_loss_b = uni_feat$has_loss_b, d_ev = uni_feat$d_ev,
               p_base = uni_feat$p_base_exp)
sets_sim_cols <- lapply(sets_sim$game_ids, function(s) {
  g <- parse_ids(s)
  cc <- match(g, uni_feat$GameID)
  stopifnot(!anyNA(cc), !anyDuplicated(cc))
  cc
})
stopifnot(all(lengths(sets_sim_cols) == sets_sim$K))

sim_one_set <- function(cols, pars, seed, n_split, n_reps) {
  set.seed(seed)
  f <- uni_feat[cols, ]
  K <- length(cols); kh <- K %/% 2L
  ns <- nrow(pars)
  NBm <- matrix(NA_real_, ns, K)
  for (i in seq_len(ns)) {
    ch <- simulate_choices(pars[i, ], f, n_reps)                 # FROZEN
    NBm[i, ] <- rowSums(matrix(ch$choice, ncol = n_reps, byrow = TRUE))
  }
  PT <- t(vapply(seq_len(ns), function(i)
    n_reps * analytic_p_b(pars[i, ], f), numeric(K)))            # FROZEN
  E <- list()
  for (tg in TARGETS) E[[tg]] <- list(h1 = matrix(NA_real_, ns, n_split),
                                      h2 = matrix(NA_real_, ns, n_split))
  for (s in seq_len(n_split)) {
    take <- sample.int(K, 2L * kh)
    i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
    a <- half_targets(NBm[, i1, drop = FALSE], uni_fm, cols[i1], n_reps)
    b <- half_targets(NBm[, i2, drop = FALSE], uni_fm, cols[i2], n_reps)
    for (tg in TARGETS) { E[[tg]]$h1[, s] <- a[[tg]]; E[[tg]]$h2[, s] <- b[[tg]] }
  }
  est_full  <- half_targets(NBm, uni_fm, cols, n_reps)
  true_full <- half_targets(PT, uni_fm, cols, n_reps)
  purrr::map_dfr(TARGETS, function(tg) {
    rs <- split_correlations(E[[tg]]$h1, E[[tg]]$h2, MIN_PAIRS)  # FROZEN
    r_raw <- mean(rs, na.rm = TRUE)
    x <- est_full[[tg]]; y <- true_full[[tg]]
    ok <- is.finite(x) & is.finite(y)
    rec <- if (sum(ok) < MIN_PAIRS || stats::sd(x[ok]) == 0 ||
               stats::sd(y[ok]) == 0) NA_real_ else
                 suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
    tibble::tibble(
      target = tg, r_raw_hp = r_raw,
      r_sb_hp = if (is.finite(r_raw)) sb_correct(r_raw) else NA_real_,  # FROZEN
      cov_hp = mean(!is.na(rs)),
      recovery_hp = rec, n_recovery_pairs = sum(ok))
  })
}

parallel::clusterExport(cl, c("uni_feat", "uni_fm", "sim_one_set",
                              "N_SPLIT_SIM", "N_REPS_SIM"),
                        envir = environment())
sim_rows <- list()
for (K in K_GRID) {
  set.seed(seed_simpar(K))
  pars <- draw_true_params(N_SIM_SUBJ)                           # FROZEN
  idx <- which(sets_sim$K == K)
  jobs <- lapply(idx, function(i)
    list(K = K, row = i, set_id = sets_sim$set_id[i], cols = sets_sim_cols[[i]],
         seed = seed_simset(K, sets_sim$set_id[i])))
  parallel::clusterExport(cl, c("pars"), envir = environment())
  t0 <- Sys.time()
  rr <- parallel::parLapply(cl, jobs, function(job)
    sim_one_set(job$cols, pars, job$seed, N_SPLIT_SIM, N_REPS_SIM) |>
      dplyr::mutate(K = job$K, row = job$row, set_id = job$set_id))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sim_rows[[length(sim_rows) + 1L]] <- dplyr::bind_rows(rr)
  say("  K %2d: %d sets, %.1f s\n", K, length(jobs), el)
}
reliability_highprec_sim <- dplyr::bind_rows(sim_rows) |>
  dplyr::mutate(
    game_ids = sets_sim$game_ids[row],
    n_loss_informative = vapply(sets_sim_cols[row],
                                function(z) sum(uni_feat$loss_informative[z]),
                                integer(1)),
    composition = "free_105", data_kind = "SIMULATION",
    n_sim_subjects = N_SIM_SUBJ, n_split = N_SPLIT_SIM, n_reps = N_REPS_SIM) |>
  dplyr::select(data_kind, composition, K, target, set_id, game_ids,
                n_loss_informative, r_sb_hp, r_raw_hp, cov_hp, recovery_hp,
                n_recovery_pairs, n_sim_subjects, n_split, n_reps) |>
  dplyr::arrange(target, K, set_id)
stopifnot(nrow(reliability_highprec_sim) == 2L * 500L * length(K_GRID))

# =============================================================================
# 9. Part B noise reduction vs E3, measured directly
# =============================================================================
# The across-set SD confounds true composition variance with estimation noise,
# so the estimation SD is measured DIRECTLY: a fixed subsample of sets is
# re-simulated N_NOISE_REPS times with independent seeds at BOTH E3's N = 300
# and E3b's N = 2000, and the SD across replicates is the estimation SD. The
# ratio of those SDs is the measured noise reduction; sqrt(2000/300) = 2.58 is
# the expected ratio.
say("\n=== Part B noise study: %d sets x %d replicates at N = %d and N = %d ===\n",
    N_NOISE_SETS, N_NOISE_REPS, N_SIM_SUBJ_E3, N_SIM_SUBJ)
noise_jobs <- list()
for (K in K_GRID) {
  idx <- which(sets_sim$K == K)[seq_len(N_NOISE_SETS)]
  for (n_sub in c(N_SIM_SUBJ_E3, N_SIM_SUBJ)) {
    for (jj in seq_len(N_NOISE_REPS)) {
      for (i in idx) {
        noise_jobs[[length(noise_jobs) + 1L]] <- list(
          K = K, set_id = sets_sim$set_id[i], cols = sets_sim_cols[[i]],
          n_sub = n_sub, rep = jj,
          seed = seed_noise(K, sets_sim$set_id[i], jj, n_sub),
          seed_par = SEED_E3B + 60000000L + 1000000L * (n_sub == N_SIM_SUBJ) +
            10000L * K + jj)
      }
    }
  }
}
t0 <- Sys.time()
noise_res <- parallel::parLapply(cl, noise_jobs, function(job) {
  set.seed(job$seed_par)
  pars <- draw_true_params(job$n_sub)                            # FROZEN
  sim_one_set(job$cols, pars, job$seed, N_SPLIT_SIM, N_REPS_SIM) |>
    dplyr::mutate(K = job$K, set_id = job$set_id, n_sim_subjects = job$n_sub,
                  replicate = job$rep)
})
say("noise study: %d runs in %.1f min\n", length(noise_jobs),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
e3b_sim_noise_replicates <- dplyr::bind_rows(noise_res) |>
  dplyr::mutate(data_kind = "SIMULATION") |>
  dplyr::select(data_kind, K, target, set_id, n_sim_subjects, replicate,
                r_sb_hp, cov_hp, recovery_hp) |>
  dplyr::arrange(target, K, set_id, n_sim_subjects, replicate)

# Column names are built from FIXED constants, never from data values, so the
# join keys cannot silently change with the simulated sample sizes.
noise_sd_long <- e3b_sim_noise_replicates |>
  dplyr::group_by(K, target, set_id, n_sim_subjects) |>
  dplyr::summarise(sd_rep = stats::sd(r_sb_hp, na.rm = TRUE),
                   n_finite = sum(is.finite(r_sb_hp)), .groups = "drop") |>
  dplyr::filter(n_finite >= 2L) |>
  dplyr::group_by(K, target, n_sim_subjects) |>
  dplyr::summarise(mean_estimation_sd = mean(sd_rep, na.rm = TRUE),
                   n_noise_sets = dplyr::n(), .groups = "drop")
noise_sd <- dplyr::full_join(
  noise_sd_long |>
    dplyr::filter(n_sim_subjects == N_SIM_SUBJ_E3) |>
    dplyr::select(K, target, estimation_sd_e3 = mean_estimation_sd,
                  n_noise_sets_e3 = n_noise_sets),
  noise_sd_long |>
    dplyr::filter(n_sim_subjects == N_SIM_SUBJ) |>
    dplyr::select(K, target, estimation_sd_hp = mean_estimation_sd,
                  n_noise_sets_hp = n_noise_sets),
  by = c("K", "target"))

e3_free_sd <- e3_free |>
  dplyr::group_by(K, target) |>
  dplyr::summarise(sd_across_sets_e3 = stats::sd(sim_r_sb, na.rm = TRUE),
                   median_r_sb_e3 = stats::median(sim_r_sb, na.rm = TRUE),
                   n_identifiable_e3 = sum(is.finite(sim_r_sb)),
                   .groups = "drop")

hp_vs_e3_pairs <- reliability_highprec_sim |>
  dplyr::select(K, target, set_id, r_sb_hp) |>
  dplyr::left_join(e3_free |> dplyr::select(K, target, set_id,
                                            r_sb_e3 = sim_r_sb),
                   by = c("K", "target", "set_id"))

highprec_vs_e3_sim <- hp_vs_e3_pairs |>
  dplyr::group_by(K, target) |>
  dplyr::summarise(
    n_sets = dplyr::n(),
    n_pairs_used = sum(is.finite(r_sb_hp) & is.finite(r_sb_e3)),
    sd_across_sets_hp = stats::sd(r_sb_hp, na.rm = TRUE),
    median_r_sb_hp = stats::median(r_sb_hp, na.rm = TRUE),
    n_identifiable_hp = sum(is.finite(r_sb_hp)),
    spearman_hp_vs_e3 = {
      ok <- is.finite(r_sb_hp) & is.finite(r_sb_e3)
      if (sum(ok) < MIN_PAIRS) NA_real_ else
        suppressWarnings(stats::cor(r_sb_hp[ok], r_sb_e3[ok],
                                    method = "spearman"))
    },
    mean_shift_hp_minus_e3 = mean(r_sb_hp - r_sb_e3, na.rm = TRUE),
    sd_shift_hp_minus_e3 = stats::sd(r_sb_hp - r_sb_e3, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::left_join(e3_free_sd, by = c("K", "target")) |>
  dplyr::left_join(noise_sd, by = c("K", "target")) |>
  dplyr::mutate(
    data_kind = "SIMULATION",
    n_sim_subjects_e3 = N_SIM_SUBJ_E3, n_sim_subjects_hp = N_SIM_SUBJ,
    noise_reduction_factor = estimation_sd_e3 / estimation_sd_hp,
    noise_reduction_expected = sqrt(N_SIM_SUBJ / N_SIM_SUBJ_E3),
    sd_across_sets_ratio = sd_across_sets_e3 / sd_across_sets_hp) |>
  dplyr::select(data_kind, K, target, n_sets, n_pairs_used,
                n_identifiable_e3, n_identifiable_hp,
                n_sim_subjects_e3, n_sim_subjects_hp,
                estimation_sd_e3, estimation_sd_hp,
                noise_reduction_factor, noise_reduction_expected,
                sd_across_sets_e3, sd_across_sets_hp, sd_across_sets_ratio,
                median_r_sb_e3, median_r_sb_hp, mean_shift_hp_minus_e3,
                sd_shift_hp_minus_e3, spearman_hp_vs_e3,
                n_noise_sets_e3, n_noise_sets_hp) |>
  dplyr::arrange(target, K)
print(as.data.frame(highprec_vs_e3_sim |>
                      dplyr::select(target, K, estimation_sd_e3,
                                    estimation_sd_hp,
                                    noise_reduction_factor,
                                    mean_shift_hp_minus_e3)))

# =============================================================================
# 10. D3 set identity, verified independently of the saved CSV
# =============================================================================
# The sets are parsed from step 22's saved files, so identity holds by
# construction. This re-derives them from E3's SEEDS alone (E3's draw_sets,
# replayed) and asserts the same membership, so the reuse is verified rather
# than merely asserted.
M_CAP_E3 <- 2000L
SEED_E3  <- 20260719L
draw_sets_e3 <- function(np, K, M, seed) {
  set.seed(seed)
  total <- choose(np, K)
  if (total <= M_CAP_E3) {
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
d3_rows <- purrr::map_dfr(K_GRID, function(K) {
  rs <- draw_sets_e3(length(fm$probs), K,
                     min(M_CAP_E3, choose(length(fm$probs), K)),
                     SEED_E3 + 1000L * GATE_FORM + K)
  rs <- rs[seq_len(min(M_CAP_E3, length(rs)))]
  mine <- sets_real_cols[sets_real$K == K]
  real_ok <- identical(lapply(rs, as.integer), lapply(mine, as.integer))
  ss <- draw_sets_e3(nrow(uni_feat), K, 500L, SEED_E3 + 700000L + K)
  ss <- ss[seq_len(min(500L, length(ss)))]
  mine_s <- sets_sim_cols[sets_sim$K == K]
  sim_ok <- identical(lapply(ss, as.integer), lapply(mine_s, as.integer))
  tibble::tibble(K = K, real_sets_match_e3_seeds = real_ok,
                 sim_sets_match_e3_seeds = sim_ok,
                 n_real_sets = length(mine), n_sim_sets = length(mine_s))
})
print(as.data.frame(d3_rows))
stopifnot(all(d3_rows$real_sets_match_e3_seeds),
          all(d3_rows$sim_sets_match_e3_seeds))
say("D3 set identity: real and simulated sets re-derived from E3's seeds alone MATCH\n")

parallel::stopCluster(cl)

# =============================================================================
# 11. Outputs
# =============================================================================
save_csv(reliability_subject_split_real, "reliability_subject_split_real")
save_csv(outcome_reliability_real, "outcome_reliability_real")
save_csv(reliability_highprec_sim, "reliability_highprec_sim")
save_csv(highprec_vs_e3_sim, "highprec_vs_e3_sim")
save_csv(e3b_sim_noise_replicates, "e3b_sim_noise_replicates")
save_csv(spot_detail, "e3b_d1_spotcheck")
save_csv(spot_summary, "e3b_d1_spotcheck_summary")
save_csv(d3_rows, "e3b_set_identity")

# Seed vectors are built BEFORE each tibble() call: inside tibble() the name `K`
# resolves to the COLUMN just created, which is recycled to the length of
# `purpose`, so calling seed_setA(K, 0L) in place would silently produce a
# length-4 seed vector.
seeds_tbl <- dplyr::bind_rows(
  tibble::tibble(
    part = "A", K = NA_integer_, purpose = "subject_partition",
    seed = SEED_E3B,
    derivation = "SEED_E3B; sample.int(63), first 32 = half A, rest = half B"),
  purrr::map_dfr(K_GRID, function(Kc) {
    sv <- c(seed_setA(Kc, 0L), seed_setB(Kc, 0L),
            seed_boot(Kc, 0L, 1L), seed_boot(Kc, 0L, 2L))
    stopifnot(length(sv) == 4L)
    tibble::tibble(
      part = "A", K = Kc,
      purpose = c("halving_half_A_base", "halving_half_B_base",
                  "bootstrap_half_A_base", "bootstrap_half_B_base"),
      seed = sv,
      derivation = c(
        "SEED_E3B + 1 + 1000000 + 10000*K; per-set seed = base + set_id",
        "SEED_E3B + 2 + 9000000 + 10000*K; per-set seed = base + set_id",
        "SEED_E3B + 20000000 + 100000*1 + 10000*K; per-set seed = base + set_id",
        "SEED_E3B + 20000000 + 100000*2 + 10000*K; per-set seed = base + set_id"))
  }),
  purrr::map_dfr(K_GRID, function(Kc) {
    sv <- seed_boot(Kc, 0L, 3L)
    stopifnot(length(sv) == 1L)
    tibble::tibble(part = "A", K = Kc,
                   purpose = "outcome_reliability_bootstrap", seed = sv,
                   derivation = "SEED_E3B + 20000000 + 100000*3 + 10000*K")
  }),
  purrr::map_dfr(K_GRID, function(Kc) {
    sv <- c(seed_simpar(Kc), seed_simset(Kc, 0L))
    stopifnot(length(sv) == 2L)
    tibble::tibble(
      part = "B", K = Kc,
      purpose = c("sim_param_draw", "sim_set_base"), seed = sv,
      derivation = c("SEED_E3B + 100 + K",
                     "SEED_E3B + 200 + 1000000 + 100000*K; per-set = base + set_id"))
  }),
  purrr::map_dfr(K_GRID, function(Kc) {
    sv <- c(seed_noise(Kc, 0L, 0L, N_SIM_SUBJ_E3),
            seed_noise(Kc, 0L, 0L, N_SIM_SUBJ),
            SEED_E3B + 60000000L + 10000L * Kc,
            SEED_E3B + 60000000L + 1000000L + 10000L * Kc)
    stopifnot(length(sv) == 4L)
    tibble::tibble(
      part = "B_noise", K = Kc,
      purpose = c("noise_choice_base_e3_n", "noise_choice_base_hp_n",
                  "noise_param_base_e3_n", "noise_param_base_hp_n"),
      seed = sv,
      derivation = c(
        "SEED_E3B + 40000000 + 200000*K; per-run = base + 100*set_id + replicate",
        "SEED_E3B + 44000000 + 200000*K; per-run = base + 100*set_id + replicate",
        "SEED_E3B + 60000000 + 10000*K; per-run = base + replicate",
        "SEED_E3B + 61000000 + 10000*K; per-run = base + replicate"))
  }),
  tibble::tibble(part = "gate", K = NA_integer_, purpose = "G_EQ",
                 seed = SEED_E3B + 77L, derivation = "SEED_E3B + 77"))
save_csv(seeds_tbl, "e3b_seeds")

cfg <- list(
  spec = "elevation/E3b_spec.md",
  status = "EXPLORATORY addendum; the frozen E1 pre-specification and the frozen H2 gate result (E4, FAIL) STAND and are unchanged",
  raw_md5 = raw_md5, seed_base = SEED_E3B,
  regime = "experience (Trial 6-25)", r_reps = R_REPS,
  targets = list(primary = "loss_avoid_adj", secondary = "ev_sens_agn"),
  K_grid = K_GRID, form = GATE_FORM,
  part_a = list(
    n_split_half = N_SPLIT_HALF, n_half_A = N_HALF_A, n_half_B = N_HALF_B,
    n_sets_per_K = 2000L, n_boot_per_set = N_BOOT_SET,
    n_boot_outcome_reliability = N_BOOT_COR,
    reliability_engine = "frozen reliability_point (calls frozen split_correlations + sb_correct)",
    outcome_reliability_engine = "frozen boot_cor (percentile bootstrap over SETS)",
    disattenuation_applied = FALSE),
  part_b = list(
    is_simulation = TRUE, n_sim_subjects = N_SIM_SUBJ,
    n_split = N_SPLIT_SIM, n_reps = N_REPS_SIM, n_sets_per_K = 500L,
    e3_n_sim_subjects = N_SIM_SUBJ_E3,
    noise_study = list(n_sets = N_NOISE_SETS, n_replicates = N_NOISE_REPS,
                       levels = c(N_SIM_SUBJ_E3, N_SIM_SUBJ)),
    population = "frozen uniform-range draw_true_params, UNCHANGED (calibration to build is deliberately not done here)"),
  constants_fixed_by_code = paste(
    "The spec fixes the seeds and the sample-size constants but leaves the",
    "bootstrap sizes, the G_EQ / spot-check sizes and the noise-study design",
    "open. N_BOOT_SET = 1000 (matching E3's noise-floor N_BOOT), N_BOOT_COR =",
    "2000 (the frozen boot_cor default), N_GEQ = 8000, N_SPOT = 12,",
    "N_NOISE_SETS = 50, N_NOISE_REPS = 20 were fixed in this script before any",
    "result existed and are recorded here."),
  seed_block_note = paste(
    "The spec fixes the BASE seeds (halving A = SEED_E3B + 1, halving B =",
    "SEED_E3B + 2). The per-set derivations are deterministic functions of",
    "those, placed in disjoint blocks so the A and B streams provably never",
    "intersect (asserted in-script)."),
  accelerator = paste(
    "The two targets are evaluated by a vectorised path using the SAME formulas",
    "as the frozen estimate_half, because estimate_half also runs the",
    "4-parameter Firth fit for b_ev/b_risk/b_loss, which are not E3b targets;",
    "Part A would otherwise need ~2.0e8 engine calls. Gate G_EQ and the D1",
    "spot-check verify agreement against the frozen engine on both subject",
    "halves. Equivalence is re-proven here, not inherited from E3."),
  g_eq_checks = eq_n, g_eq_max_abs_diff = max(eq_max_ev, eq_max_loss),
  g_eq_na_mismatches = eq_na,
  d1_spotcheck_max_abs_diff_loss_avoid_adj = spot_loss,
  d1_spotcheck_max_abs_diff_ev_sens_agn = spot_ev,
  d3_set_identity_verified_from_e3_seeds = all(d3_rows$real_sets_match_e3_seeds) &&
    all(d3_rows$sim_sets_match_e3_seeds),
  split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
  loss_baselines_exp = as.list(baselines$exp),
  data_inputs = c("data/raw/All estimation raw data.csv",
                  "elevation/results/e1/participant_split.csv",
                  "elevation/results/e2/problem_information_dictionary.csv",
                  "elevation/results/e3/reliability_by_set_real.csv",
                  "elevation/results/e3/reliability_recovery_by_set_sim.csv"),
  build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
  no_freeze_or_test_choice_read = TRUE,
  optimization_performed = FALSE, mechanism_model_fitted = FALSE,
  candidate_design_constructed = FALSE,
  h2_gate_rerun_or_amended = FALSE,
  runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
jsonlite::write_json(cfg, file.path(res_dir, "e3b_config.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

say("\nE3b complete in %.1f min. Outputs in elevation/results/e3b/\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))
say("Outcome reliability (cross-half Spearman), primary target:\n")
print(as.data.frame(outcome_reliability_real |>
                      dplyr::filter(target == "loss_avoid_adj") |>
                      dplyr::select(K, n_sets_used, outcome_reliability,
                                    ci_lo, ci_hi)))

# Step 33: verification and assembly of steps 30-32
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e1/, e2/, dispb/ (steps 20-21, 30-32); steps 30-33 (source text)
# Writes:  output/elevation/results/dispb/ (dispb_gate.json, dispb_verify.json, spot checks, engine provenance)
# Paper:   Study 2 and Supplement, through step 36
# Runtime: about 1 minute with 12 parallel workers

# Verify the scoring engine for steps 30-32 (gate D1), then assemble their
# outputs. B1 and B2 use a vectorised path for loss_avoid_adj and ev_sens_agn
# to avoid the unused four-parameter Firth fits for b_ev / b_risk / b_loss.
# G_EQ checks this path against estimate_half to < 1e-12 with identical NA
# patterns. Whole-set r_sb checks use the original engine on all 63 Set-1
# build subjects and on a disjoint subject half.
#
# When steps 30-32 have written their part files, combine the gate,
# configuration, and seed files. G_EQ uses the same seed on each run.
#
# Only the 344 build subjects from step 20 are read. Choice/outcome columns
# are dropped from the full-sample frame, excluding the 137 freeze and 205 test
# subjects' choices. The shared functions are sourced without redefinition
# here and in steps 30-32:
#   R/step2_reliability.R  build_subject_data, fit_loss_baselines,
#                          estimate_half, split_correlations, sb_correct,
#                          reliability_point
#   R/cpc18_features.R     load_cpc18, is_clean_row, clean_problem_features,
#                          add_scaled_features
# No designs are optimized or constructed, no new choice model is fitted,
# and the frozen E3 h1_gate and E4 h2_gate are neither rerun nor amended.

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
# 0. Constants (fixed before any result exists)
# =============================================================================
RAW_MD5   <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_B1   <- 20260722L
SEED_GEQ  <- SEED_B1 + 900000L      # disjoint from every production seed block
REGIME    <- "exp"
POS       <- 6:25                   # experience regime, r = 20 (the full pool)
R_REPS    <- length(POS)
TARGETS   <- c("loss_avoid_adj", "ev_sens_agn")   # primary, secondary
GATE_FORM <- 1L
GATE_K    <- 12L
MIN_PAIRS <- 10L                    # frozen split_correlations default
N_GEQ     <- 8000L                  # halves re-derived through estimate_half
N_SPOT    <- 8L                     # whole sets re-derived through the engine
N_SPLIT_SPOT <- 100L                # matches B1's N_SPLIT_B1
SPOT_TOL_EV  <- 5e-3                # ev_sens_agn rank-tie tolerance (E3/E3b)
N_HALF_A  <- 32L                    # largest-remainder halves of 63
N_WORKERS <- 12L

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
stopifnot(!any(choice_cols %in% names(design_all)))
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
# Must match the build-only loss baselines obtained in E3 (step 22).
stopifnot(abs(baselines$exp[1] - 0.4301295476692) < 1e-10,
          abs(baselines$exp[2] - 0.165107793793507) < 1e-10)
say("frozen loss baselines (exp): b0=%.13f b1=%.15f (match E3 constants)\n",
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

# =============================================================================
# 3. The accelerator, byte-equivalent to E3 / E3b (verified below, not assumed)
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

run_set_rows <- function(fm, cols, rows, n_split, seed) {
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
  list(n_subjects = ns, rows_kept = rows[keep], E = E)
}

draw_sets <- function(np, K, M, seed, m_cap = 2000L) {
  set.seed(seed)
  total <- choose(np, K)
  if (total <= m_cap) {
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
# 4. GATE G_EQ: the accelerator vs the FROZEN estimate_half
# =============================================================================
# Drawn across ALL seven forms and, on Set 1, across BOTH the full 63 subjects
# and the two disjoint subject halves that B1's covariance estimand and B2 use,
# so the sample covers every regime in which the accelerator is exercised.
say("\n=== G_EQ: accelerator vs frozen estimate_half (%d halves) ===\n", N_GEQ)
set.seed(SEED_B1)
perm1  <- sample.int(63L)
rows_A <- sort(perm1[seq_len(N_HALF_A)])
rows_B <- sort(perm1[(N_HALF_A + 1L):63L])
stopifnot(length(intersect(rows_A, rows_B)) == 0L,
          setequal(c(rows_A, rows_B), seq_len(63L)))

set.seed(SEED_GEQ)
eq_rows <- vector("list", N_GEQ)
eq_max_ev <- 0; eq_max_loss <- 0; eq_na <- 0L; eq_n <- 0L
for (it in seq_len(N_GEQ)) {
  S <- sample(1:7, 1)
  fm <- form_mats[[S]]
  if (S == GATE_FORM) {
    regime_rows <- switch(sample(1:3, 1), seq_len(63L), rows_A, rows_B)
  } else {
    regime_rows <- seq_along(fm$subs)
  }
  cand <- regime_rows[!is.na(rowSums(fm$NB[regime_rows, , drop = FALSE]))]
  if (!length(cand)) next
  ri  <- cand[sample.int(length(cand), 1)]
  sid <- fm$subs[ri]
  sd  <- subj_data[[sid]]
  kh   <- sample(3:6, 1)
  cols <- sample(seq_along(fm$probs), min(kh, length(fm$probs)))
  a <- estimate_half(sd, match(fm$probs[cols], sd$feat$GameID), POS, REGIME)
  b <- half_targets(fm$NB[ri, cols, drop = FALSE], fm, cols, R_REPS)
  drow <- c(ev_sens_agn = NA_real_, loss_avoid_adj = NA_real_)
  for (nm in TARGETS) {
    x <- unname(a[[nm]]); y <- unname(b[[nm]])
    if (is.na(x) != is.na(y)) { eq_na <- eq_na + 1L; next }
    if (!is.na(x)) {
      d <- abs(x - y)
      drow[nm] <- d
      if (nm == "ev_sens_agn") eq_max_ev <- max(eq_max_ev, d)
      else eq_max_loss <- max(eq_max_loss, d)
    }
  }
  eq_n <- eq_n + 1L
  eq_rows[[it]] <- tibble::tibble(
    check = it, form = S, subject_row = ri, k_half = length(cols),
    game_ids = paste(fm$probs[cols], collapse = ";"),
    frozen_ev_sens_agn = unname(a[["ev_sens_agn"]]),
    accel_ev_sens_agn = unname(b[["ev_sens_agn"]]),
    frozen_loss_avoid_adj = unname(a[["loss_avoid_adj"]]),
    accel_loss_avoid_adj = unname(b[["loss_avoid_adj"]]),
    abs_diff_ev_sens_agn = unname(drow["ev_sens_agn"]),
    abs_diff_loss_avoid_adj = unname(drow["loss_avoid_adj"]))
}
geq_tbl <- dplyr::bind_rows(eq_rows)
say("G_EQ: %d halves; max|d ev_sens_agn| = %.3g; max|d loss_avoid_adj| = %.3g; NA mismatches = %d\n",
    eq_n, eq_max_ev, eq_max_loss, eq_na)
stopifnot(eq_na == 0L, eq_max_ev < 1e-12, eq_max_loss < 1e-12)
save_csv(geq_tbl, "dispb_geq")

# =============================================================================
# 5. Whole-set spot check through the FROZEN engine (full 63 and one half)
# =============================================================================
# G_EQ proves per-value equivalence. This closes the loop at the level of a
# whole set's r_sb, on the full Set-1 sample AND on a disjoint subject half,
# entirely through frozen estimate_half + split_correlations + sb_correct.
# Tolerances differ by target for the reason E3 measured: loss_avoid_adj
# arithmetic is identical in both paths so every split's Spearman agrees
# EXACTLY (tolerance 0); ev_sens_agn values are built from integer choice
# counts, so subjects hold exactly tied values that split_correlations ranks,
# and a last-bit difference resolves a tie the other way. That is rank-tie
# sensitivity of the FROZEN estimator, not of the accelerator.
say("\n=== D1 spot check: %d whole sets re-derived through the frozen engine ===\n",
    N_SPOT)
fm1 <- form_mats[[GATE_FORM]]
spot_sets <- draw_sets(length(fm1$probs), GATE_K, 2000L,
                       SEED_B1 + 1000L * GATE_FORM + GATE_K + 10000L)
spot_jobs <- lapply(seq_len(N_SPOT), function(i)
  list(set_id = i, cols = spot_sets[[i]],
       seed = SEED_B1 + 1000L * GATE_FORM + GATE_K + 10000L + 500000L + i))

cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(tidyverse); library(here) })
  source(here("R", "step2_reliability.R"))
  NULL
}))
parallel::clusterExport(cl, c("fm1", "subj_data", "POS", "REGIME", "TARGETS",
                              "R_REPS", "MIN_PAIRS", "N_SPLIT_SPOT", "GATE_K",
                              "half_targets", "ev_ranks", "run_set_rows",
                              "rows_A"), envir = environment())
spot_detail <- dplyr::bind_rows(parallel::parLapply(cl, spot_jobs, function(job) {
  cols <- job$cols; kh <- GATE_K %/% 2L
  purrr::map_dfr(c("full63", "halfA"), function(scope) {
    rows <- if (scope == "full63") seq_along(fm1$subs) else rows_A
    subs <- fm1$subs[rows]
    Ef <- list()
    for (tg in TARGETS) Ef[[tg]] <- list(
      h1 = matrix(NA_real_, length(subs), N_SPLIT_SPOT),
      h2 = matrix(NA_real_, length(subs), N_SPLIT_SPOT))
    set.seed(job$seed)
    ties <- integer(N_SPLIT_SPOT)
    for (s in seq_len(N_SPLIT_SPOT)) {
      take <- sample.int(GATE_K, 2L * kh)
      i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
      for (u in seq_along(subs)) {
        sd <- subj_data[[subs[u]]]
        r1 <- estimate_half(sd, match(fm1$probs[cols[i1]], sd$feat$GameID),
                            POS, REGIME)                              # FROZEN
        r2 <- estimate_half(sd, match(fm1$probs[cols[i2]], sd$feat$GameID),
                            POS, REGIME)                              # FROZEN
        for (tg in TARGETS) { Ef[[tg]]$h1[u, s] <- r1[[tg]]
                              Ef[[tg]]$h2[u, s] <- r2[[tg]] }
      }
      v <- Ef$ev_sens_agn$h1[, s]; v <- v[is.finite(v)]
      ties[s] <- length(v) - length(unique(v))
    }
    Ev <- run_set_rows(fm1, cols, rows, N_SPLIT_SPOT, job$seed)$E
    purrr::map_dfr(TARGETS, function(tg) {
      rf <- split_correlations(Ef[[tg]]$h1, Ef[[tg]]$h2, MIN_PAIRS)    # FROZEN
      rv <- split_correlations(Ev[[tg]]$h1, Ev[[tg]]$h2, MIN_PAIRS)    # FROZEN
      tibble::tibble(
        set_id = job$set_id, scope = scope, target = tg,
        n_subjects = length(subs),
        frozen_r_sb = sb_correct(mean(rf, na.rm = TRUE)),              # FROZEN
        accel_r_sb  = sb_correct(mean(rv, na.rm = TRUE)),              # FROZEN
        abs_diff = abs(sb_correct(mean(rf, na.rm = TRUE)) -
                         sb_correct(mean(rv, na.rm = TRUE))),
        n_splits_differing = sum(abs(rf - rv) > 0, na.rm = TRUE),
        max_split_abs_diff = suppressWarnings(max(abs(rf - rv), na.rm = TRUE)),
        mean_tied_subjects_per_split = mean(ties))
    })
  })
}))
parallel::stopCluster(cl)

spot_summary <- spot_detail |>
  dplyr::group_by(target, scope) |>
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
stopifnot(spot_loss == 0)                 # primary target: bit-identical
stopifnot(spot_ev < SPOT_TOL_EV)          # secondary: tie resolution only
say("D1 spot check PASS: loss_avoid_adj max |d r_sb| = %.3g (exact); ev_sens_agn max = %.3g (< %.0e)\n",
    spot_loss, spot_ev, SPOT_TOL_EV)
save_csv(spot_detail, "dispb_spotcheck")
save_csv(spot_summary, "dispb_spotcheck_summary")

# =============================================================================
# 6. Source-level provenance: zero redefinitions of any frozen function
# =============================================================================
frozen_fns <- c("build_subject_data", "fit_loss_baselines", "firth_half_fit",
                "estimate_half", "run_cell", "split_correlations", "sb_correct",
                "reliability_point", "load_cpc18", "is_clean_row",
                "clean_problem_features", "add_scaled_features")
dispb_scripts <- c("33_partition_checks_verify.R", "30_partition_stability.R", "31_selection_optimism.R",
                   "32_leave_one_feature_out.R")
redef <- purrr::map_dfr(dispb_scripts, function(f) {
  p <- here("steps", f)
  if (!file.exists(p)) return(tibble::tibble(
    script = f, exists = FALSE, n_redefinitions = NA_integer_,
    redefined = NA_character_))
  txt <- readLines(p, warn = FALSE)
  hit <- vapply(frozen_fns, function(fn)
    any(grepl(paste0("^\\s*", fn, "\\s*(<-|=)\\s*function"), txt)), logical(1))
  tibble::tibble(script = f, exists = TRUE, n_redefinitions = sum(hit),
                 redefined = paste(frozen_fns[hit], collapse = ";"))
})
print(as.data.frame(redef))
stopifnot(all(redef$n_redefinitions[redef$exists] == 0L))
say("engine provenance: 0 redefinitions of %d frozen functions across %d scripts\n",
    length(frozen_fns), sum(redef$exists))
save_csv(redef, "dispb_engine_provenance")

verify <- list(
  spec = "elevation/DISPB_spec.md",
  gate = "D1 engine provenance",
  raw_md5 = raw_md5,
  split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
  loss_baselines_exp = as.list(baselines$exp),
  g_eq = list(halves_checked = eq_n, seed = SEED_GEQ,
              max_abs_diff_ev_sens_agn = eq_max_ev,
              max_abs_diff_loss_avoid_adj = eq_max_loss,
              na_mismatches = eq_na, tolerance = 1e-12,
              pass = eq_na == 0L && eq_max_ev < 1e-12 && eq_max_loss < 1e-12),
  spot_check = list(n_sets = N_SPOT, n_split = N_SPLIT_SPOT,
                    scopes = c("full63", "halfA"),
                    max_abs_diff_loss_avoid_adj = spot_loss,
                    max_abs_diff_ev_sens_agn = spot_ev,
                    tolerance_ev = SPOT_TOL_EV,
                    pass = spot_loss == 0 && spot_ev < SPOT_TOL_EV),
  frozen_function_redefinitions = sum(redef$n_redefinitions, na.rm = TRUE),
  build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
  optimization_performed = FALSE, candidate_design_constructed = FALSE,
  new_choice_model_fitted = FALSE,
  frozen_gates_rerun_or_amended = FALSE)
jsonlite::write_json(verify, file.path(res_dir, "dispb_verify.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

# =============================================================================
# 7. Assembly (only once all three probe part files exist)
# =============================================================================
parts <- c(B1 = "dispb_part_b1.json", B2 = "dispb_part_b2.json",
           B3 = "dispb_part_b3.json")
have <- file.exists(file.path(res_dir, parts))
if (all(have)) {
  say("\n=== assembling dispb_gate.json / dispb_config.json / dispb_seeds.csv ===\n")
  P <- lapply(parts, function(f)
    jsonlite::read_json(file.path(res_dir, f), simplifyVector = TRUE))
  gate <- list(
    spec = "elevation/DISPB_spec.md",
    status = paste("newly pre-declared robustness re-analysis of the frozen",
                   "Study 2 result; NOT part of the frozen E1 confirmatory arc"),
    frozen_results_untouched = list(
      study1_tag = "paper1-final-2026-07-17",
      e3_h1_verdict = "PASS (unchanged; h1_gate.json neither re-run nor amended)",
      e4_h2_verdict = "FAIL (unchanged; h2_gate.json neither re-run nor amended)"),
    d1_engine_provenance = verify$g_eq,
    B1 = P$B1$gate, B2 = P$B2$gate, B3 = P$B3$gate,
    honesty_guards = list(
      HG_ceiling = list(B1 = P$B1$gate$HG_ceiling, B2 = P$B2$gate$HG_ceiling),
      HG_frozen = list(
        e3_h1_gate_rerun_or_amended = FALSE,
        e4_h2_gate_rerun_or_amended = FALSE,
        b3_recomputes_reliability = FALSE,
        b1_b2_reliabilities_from_build_raw_with_frozen_engine = TRUE),
      HG_null = list(B1 = P$B1$gate$`G-B1c`, B2 = P$B2$gate$`G-B2c`)),
    build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
    optimization_performed = FALSE, candidate_design_constructed = FALSE,
    new_choice_model_fitted = FALSE)
  jsonlite::write_json(gate, file.path(res_dir, "dispb_gate.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null",
                       digits = NA)
  cfg <- list(
    spec = "elevation/DISPB_spec.md",
    raw_md5 = raw_md5,
    split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
    loss_baselines_exp = as.list(baselines$exp),
    d1_verification = verify,
    B1 = P$B1$config, B2 = P$B2$config, B3 = P$B3$config)
  jsonlite::write_json(cfg, file.path(res_dir, "dispb_config.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null",
                       digits = NA)
  seeds <- dplyr::bind_rows(lapply(c("b1", "b2", "b3"), function(p)
    readr::read_csv(file.path(res_dir, paste0("dispb_seeds_", p, ".csv")),
                    show_col_types = FALSE)))
  save_csv(seeds, "dispb_seeds")
  say("assembled: %d seed rows\n", nrow(seeds))
} else {
  say("\nassembly skipped: part files present = %s (run the three probes first)\n",
      paste(names(parts)[have], collapse = ","))
}

say("\ndispb_verify complete in %.1f min\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))

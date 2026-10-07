# Step 22: reliability of equal-length problem subsets (H1)
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e1/, e2/ (steps 20-21)
# Writes:  output/elevation/results/e3/ (h1_gate.json, reliability_by_set_real.csv, heterogeneity_summary.csv,
#          simulation tables, seeds, configuration)
# Paper:   Methods (Equal-length subset benchmark and mechanism model); Study 2 (Equal-length heterogeneity)
# Runtime: about 27 minutes with 12 parallel workers

# E3: same-count subset benchmark and the H1 premise gate.
#
# QUESTION (H1): at a fixed number K of unique clean problems, does held-out
# reliability vary MATERIALLY across item sets? A FAIL is a valid scientific
# result and is reported as such; it is never engineered away.
#
# BUILD-ONLY. Every quantity uses the 344 BUILD subjects of the frozen
# participant split (step 20). The 137 FREEZE and 205 TEST subjects are NOT
# read here. The guarantee is structural, as in step 21: the full-sample frame
# has every choice / outcome column DROPPED at load, so only the
# build-restricted frame can carry a choice.
#
# NO optimization, NO mechanism-model fit (that is step 25), NO reading of any
# freeze/test outcome. Part C is SIMULATION and is labeled as such everywhere.
#
# ENGINE PROVENANCE. The frozen shared engines are SOURCED unchanged and none
# is redefined:
#   R/step2_reliability.R  build_subject_data, fit_loss_baselines,
#                          estimate_half, split_correlations, sb_correct,
#                          reliability_point
#   R/cpc18_features.R     load_cpc18, is_clean_row, clean_problem_features,
#                          add_scaled_features
#   R/simulate_recovery.R  draw_true_params, analytic_p_b, simulate_choices,
#                          fit_subject_firth
# split_correlations, sb_correct, reliability_point, fit_subject_firth,
# simulate_choices, draw_true_params and analytic_p_b are called DIRECTLY on
# every reported number.
#
# ACCELERATOR (not a redefinition). Part A requires ~3.2e8 estimate_half
# evaluations; measured at 1.9-3.6 ms per call that is ~9 days of
# single-threaded compute, because estimate_half always runs the 4-parameter
# Firth fit that yields b_ev/b_risk/b_loss - and those three are NOT E3
# targets. `half_targets()`
# below evaluates ONLY the two E3 targets, using the SAME formulas as the
# frozen estimate_half, vectorised over subjects. It is not a substitute
# definition: gate G_EQ re-runs the frozen estimate_half and asserts agreement
# to < 1e-12 with identical NA patterns before any result is produced, and the
# D1 fixed-set spot-check re-derives a whole cell's per-split Spearman through the
# frozen engine.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})

source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))
source(here("R", "simulate_recovery.R"))

res_dir <- here("output", "elevation", "results", "e3")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  invisible(df)
}

t_start <- Sys.time()
say <- function(...) cat(sprintf(...), sep = "")

# ---- frozen constants (fixed BEFORE any result exists) -----------------------
RAW_MD5      <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_E3      <- 20260719L
REGIME       <- "exp"
POS          <- 6:25           # experience regime, r = 20 (the full pool)
R_REPS       <- length(POS)
K_GRID       <- c(6L, 8L, 10L, 12L)
TARGETS      <- c("loss_avoid_adj", "ev_sens_agn")   # primary, secondary
M_CAP        <- 2000L
N_SPLIT_A    <- 100L
N_NOISE_SETS <- 50L
N_BOOT       <- 1000L
N_PRED_SETS  <- 300L
N_SETS_SIM   <- 500L
N_SIM_SUBJ   <- 300L
N_REPS_SIM   <- 20L
N_SPLIT_C    <- 50L
MIN_PAIRS    <- 10L            # frozen split_correlations default
N_SPOT       <- 12L            # fixed sets re-derived through the frozen engine
SPOT_TOL_EV  <- 5e-3           # ev_sens_agn only; rank-tie resolution (see D1)
H1_SIGMA_MIN <- 0.05
H1_ICC_MIN   <- 0.50
GATE_FORM    <- 1L
GATE_K       <- 12L
GATE_TARGET  <- "loss_avoid_adj"
COVERAGE_CONTINGENCY <- 0.50
N_WORKERS    <- 12L

# Seed derivations. The CELL-level seeds are pre-declared; the per-set
# derivations below are deterministic functions of those, so set membership and
# every stochastic draw reproduce exactly under any execution order.
seed_subset <- function(form, K) SEED_E3 + 1000L * form + K
seed_split  <- function(form, K) SEED_E3 + 1000L * form + K + 500000L
seed_noise  <- function(form, K, i) SEED_E3 + 600000L + 1000L * form + K + i
SEED_PRED   <- SEED_E3 + 900000L
seed_simset <- function(K) SEED_E3 + 700000L + K
seed_simpar <- function(K) SEED_E3 + 800000L + K

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
build_ids  <- split_tbl$SubjID[split_tbl$split == "build"]
held_out   <- split_tbl$SubjID[split_tbl$split != "build"]
n_build    <- length(build_ids)
n_freeze   <- sum(split_tbl$split == "freeze")
n_test     <- sum(split_tbl$split == "test")
stopifnot(n_build == 344L, n_freeze == 137L, n_test == 205L)   # frozen split sizes
stopifnot(!anyDuplicated(build_ids))
stopifnot(all(split_tbl$held_out_until_E7 == (split_tbl$split == "test")))
say("split: build %d / freeze %d (E6, not read) / test %d (sealed to E7)\n",
    n_build, n_freeze, n_test)

raw <- load_cpc18(raw_path)
stopifnot(all(cpc18_expected_columns %in% names(raw)))

# Structural build-only guarantee (as in step 21): the full-sample frame keeps no
# choice/outcome column, so no freeze or test choice can reach any computation.
choice_cols <- c("B", "Button", "Payoff", "Forgone", "RT", "Apay", "Bpay",
                 "Feedback")
design_all <- raw |> dplyr::select(-dplyr::any_of(choice_cols))
stopifnot(!any(choice_cols %in% names(design_all)))

build_raw <- raw |> dplyr::filter(SubjID %in% build_ids)
rm(raw); invisible(gc())
stopifnot(setequal(unique(build_raw$SubjID), build_ids))
stopifnot(!any(build_raw$SubjID %in% held_out))                # build-only

# Problem features come from the choice-free frame (frozen loaders).
feat_all <- clean_problem_features(design_all) |> add_scaled_features()
feat_all <- feat_all |> dplyr::arrange(GameID)
stopifnot(nrow(feat_all) == 105L)

# The step-21 problem dictionary is the item universe of record; assert we match it.
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

# Frozen build-only loss baselines and per-subject data structures.
baselines <- fit_loss_baselines(build_trials)
subj_data <- build_subject_data(build_trials, baselines)
stopifnot(length(subj_data) == 344L)
say("frozen loss baselines (exp): b0=%.6f b1=%.6f\n",
    baselines$exp[1], baselines$exp[2])

# p_base_exp is a function of the PROBLEM only (frozen formula); assert that,
# then build the problem-level baseline vector used by the vectorised path.
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
# 2. Form structure. The common-set premise (a form's build subjects all saw
#    all of the form's clean problems) holds for Sets 1-5 (exactly one
#    distinct clean set each) but is FALSE for Sets 6 and 7, whose subjects saw
#    varying subsets of the form union (ByProb-only forms).
#    Rule applied, which REDUCES to the whole form wherever the premise holds:
#    the eligible subjects for a set S are the form's build subjects who saw
#    EVERY problem in S. n_subjects is reported per set.
# =============================================================================
form_of <- build_trials |> dplyr::distinct(SubjID, Set)
stopifnot(nrow(form_of) == 344L)

form_problems <- lapply(1:7, function(S)
  sort(unique(build_trials$GameID[build_trials$Set == S])))
names(form_problems) <- as.character(1:7)

form_subjects <- lapply(1:7, function(S)
  as.character(sort(form_of$SubjID[form_of$Set == S])))
names(form_subjects) <- as.character(1:7)

# NB[subject, problem] = number of B choices over the r = 20 experience trials;
# NA where that subject never saw that problem. pos is FIXED at the full
# experience pool, so NB is a constant per form and every half is a lookup.
make_form_mats <- function(S) {
  probs <- form_problems[[as.character(S)]]
  subs  <- form_subjects[[as.character(S)]]
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

form_struct <- purrr::map_dfr(1:7, function(S) {
  keys <- vapply(form_subjects[[as.character(S)]], function(sid)
    paste(sort(subj_data[[sid]]$feat$GameID), collapse = ","), character(1))
  tibble::tibble(
    form = S,
    n_build_subjects = length(form_subjects[[as.character(S)]]),
    n_clean_problems_union = length(form_problems[[as.character(S)]]),
    n_distinct_clean_sets = dplyr::n_distinct(keys),
    n_subjects_seeing_union = sum(!is.na(rowSums(form_mats[[S]]$NB))),
    n_loss_informative = sum(form_mats[[S]]$li),
    common_set_premise_holds = dplyr::n_distinct(keys) == 1L)
})
print(form_struct)
# Sets 1-5 must satisfy the common-set premise exactly.
stopifnot(all(form_struct$common_set_premise_holds[1:5]))

# =============================================================================
# 3. The two E3 targets, vectorised over subjects (formulas identical to the
#    frozen estimate_half; equivalence asserted at gate G_EQ below)
# =============================================================================
# ev_sens_agn = cor(b_vec, ev_vec, "spearman") where b_vec is the half's 0/1
# per-trial choices and ev_vec repeats each problem's d_ev r times. rank(b_vec)
# is a POSITIVE affine map of b_vec (n0 zeros share rank (n0+1)/2, n1 ones
# share n0+(n1+1)/2), and Pearson is invariant to that, so the Spearman equals
# Pearson(b_vec, rank(ev_vec)). The tie-averaged ranks of ev_vec depend only on
# the half's d_ev multiset, so they are computed once per half and reused for
# every subject.
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

# NBh: subjects x k matrix of B counts on the half's problems.
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
    # counts = has_loss_b ? r - nB : nB   (frozen estimate_half convention)
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
# =============================================================================
say("\n=== G_EQ: vectorised targets vs frozen estimate_half ===\n")
set.seed(SEED_E3)
eq_n <- 0L; eq_max_ev <- 0; eq_max_loss <- 0; eq_na <- 0L
for (it in seq_len(6000)) {
  S <- sample(1:7, 1)
  fm <- form_mats[[S]]
  keep <- which(!is.na(rowSums(fm$NB)))
  sid <- fm$subs[sample(keep, 1)]
  sd <- subj_data[[sid]]
  kh <- sample(3:6, 1)
  cols <- sample(seq_along(fm$probs), kh)
  a <- estimate_half(sd, match(fm$probs[cols], sd$feat$GameID), POS, REGIME)
  b <- half_targets(fm$NB[sid, cols, drop = FALSE], fm, cols, R_REPS)
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
# 5. Part A: real-data within-form reliability heterogeneity + the H1 gate
# =============================================================================
# M distinct K-subsets, uniform without replacement; enumerated when C <= M_CAP.
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
      assign(key, TRUE, envir = seen)
      got <- got + 1L
      out[[got]] <- s
    }
  }
  out
}

# One set: r_sb per target over N_SPLIT common halvings of S. The halving is
# shared across subjects (S is a FIXED instrument): S is split into two halves
# of K/2 and r_sb is averaged over N_split random halvings.
run_set <- function(fm, cols, n_split, seed, want_est = FALSE) {
  NBs <- fm$NB[, cols, drop = FALSE]
  keep <- which(!is.na(rowSums(NBs)))          # subjects who saw ALL of S
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
    rs <- split_correlations(E[[tg]]$h1, E[[tg]]$h2, MIN_PAIRS)   # FROZEN
    r_raw <- mean(rs, na.rm = TRUE)
    out[[tg]] <- list(
      r_raw = r_raw,
      r_sb = if (is.finite(r_raw)) sb_correct(r_raw) else NA_real_,  # FROZEN
      coverage = mean(!is.na(rs)),
      n_splits_used = sum(!is.na(rs)))
  }
  if (want_est) out$E <- E
  out
}

say("\n=== Part A: real within-form heterogeneity ===\n")
cellsA <- purrr::map_dfr(1:7, function(S) {
  np <- length(form_problems[[as.character(S)]])
  tibble::tibble(form = S, K = K_GRID[K_GRID <= np])
})
say("feasible (form, K) cells: %d\n", nrow(cellsA))

rel_rows <- list(); noise_jobs <- list(); cell_meta <- list()
for (ci in seq_len(nrow(cellsA))) {
  S <- cellsA$form[ci]; K <- cellsA$K[ci]
  fm <- form_mats[[S]]
  np <- length(fm$probs)
  total <- choose(np, K)
  M <- min(M_CAP, total)
  sets <- draw_sets(np, K, M, seed_subset(S, K))
  stopifnot(length(sets) == M)
  stopifnot(!anyDuplicated(vapply(sets, paste, character(1), collapse = ",")))
  sseed <- seed_split(S, K)
  t0 <- Sys.time()
  res <- vector("list", M)
  for (i in seq_len(M)) res[[i]] <- run_set(fm, sets[[i]], N_SPLIT_A, sseed + i)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  for (tg in TARGETS) {
    rel_rows[[length(rel_rows) + 1L]] <- tibble::tibble(
      form = S, K = K, target = tg, set_id = seq_len(M),
      game_ids = vapply(sets, function(z) paste(fm$probs[z], collapse = ";"),
                        character(1)),
      n_loss_informative = vapply(sets, function(z) sum(fm$li[z]), integer(1)),
      n_subjects = vapply(res, function(z) z$n_subjects, integer(1)),
      r_raw = vapply(res, function(z) z[[tg]]$r_raw, numeric(1)),
      r_sb = vapply(res, function(z) z[[tg]]$r_sb, numeric(1)),
      coverage = vapply(res, function(z) z[[tg]]$coverage, numeric(1)),
      n_splits_used = vapply(res, function(z) z[[tg]]$n_splits_used, integer(1)))
  }
  cell_meta[[ci]] <- tibble::tibble(form = S, K = K, M = M,
                                    n_combinations = total,
                                    enumerated = total <= M_CAP,
                                    secs = el)
  noise_jobs[[ci]] <- list(form = S, K = K, sets = sets[seq_len(min(N_NOISE_SETS, M))])
  say("  form %d K %2d: M = %4d (%s), %.1f s\n", S, K, M,
      if (total <= M_CAP) "enumerated" else "sampled", el)
}
reliability_by_set_real <- dplyr::bind_rows(rel_rows)
cell_meta <- dplyr::bind_rows(cell_meta)

# ---- noise floor: FROZEN reliability_point subject-cluster bootstrap --------
# sigma2_noise = mean over the 50-set subsample of each set's estimation
# variance. reliability_point returns a percentile CI, not a variance, so the
# per-set SD is the normal-approximation SD of that frozen CI on the SB scale,
# sd = (ci_hi_sb - ci_lo_sb) / (2 * qnorm(.975)).
say("\n=== Part A noise floor (frozen reliability_point, N_boot = %d) ===\n", N_BOOT)
Z975 <- stats::qnorm(0.975)
cl <- parallel::makeCluster(N_WORKERS)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(tidyverse); library(here) })
  source(here("R", "step2_reliability.R"))
  NULL
}))
parallel::clusterExport(cl, c("form_mats", "half_targets", "ev_ranks", "run_set",
                              "TARGETS", "R_REPS", "MIN_PAIRS", "N_SPLIT_A",
                              "N_BOOT", "seed_split", "seed_noise", "Z975",
                              "SEED_E3"), envir = environment())
noise_res <- parallel::parLapply(cl, noise_jobs, function(job) {
  fm <- form_mats[[job$form]]
  sseed <- seed_split(job$form, job$K)
  purrr::map_dfr(seq_along(job$sets), function(i) {
    r <- run_set(fm, job$sets[[i]], N_SPLIT_A, sseed + i, want_est = TRUE)
    purrr::map_dfr(TARGETS, function(tg) {
      rp <- reliability_point(r$E[[tg]]$h1, r$E[[tg]]$h2, N_BOOT,
                              seed_noise(job$form, job$K, i), MIN_PAIRS)
      tibble::tibble(form = job$form, K = job$K, target = tg, set_id = i,
                     r_sb = rp$r_sb, ci_lo_sb = rp$ci_lo_sb,
                     ci_hi_sb = rp$ci_hi_sb,
                     sd_noise = (rp$ci_hi_sb - rp$ci_lo_sb) / (2 * Z975))
    })
  })
})
noise_tbl <- dplyr::bind_rows(noise_res)
say("noise floor computed for %d set-targets\n", nrow(noise_tbl))

# ---- variance decomposition per cell ----------------------------------------
q <- function(x, p) if (all(is.na(x))) NA_real_ else
  unname(stats::quantile(x, p, na.rm = TRUE))
heterogeneity_summary <- reliability_by_set_real |>
  dplyr::group_by(form, K, target) |>
  dplyr::summarise(
    M = dplyr::n(),
    n_sets_identified = sum(is.finite(r_sb)),
    set_coverage = mean(is.finite(r_sb)),
    split_coverage = mean(coverage),
    n_subjects_median = stats::median(n_subjects),
    sigma2_set_raw = stats::var(r_sb, na.rm = TRUE),
    r_sb_min = min(r_sb, na.rm = TRUE), r_sb_q25 = q(r_sb, .25),
    r_sb_median = q(r_sb, .50), r_sb_q75 = q(r_sb, .75),
    r_sb_max = max(r_sb, na.rm = TRUE),
    r_sb_iqr = q(r_sb, .75) - q(r_sb, .25),
    best_minus_median = max(r_sb, na.rm = TRUE) - q(r_sb, .50),
    .groups = "drop") |>
  dplyr::left_join(
    noise_tbl |> dplyr::group_by(form, K, target) |>
      dplyr::summarise(sigma2_noise = mean(sd_noise^2, na.rm = TRUE),
                       n_noise_sets = sum(is.finite(sd_noise)), .groups = "drop"),
    by = c("form", "K", "target")) |>
  dplyr::mutate(
    sigma_set_raw = sqrt(sigma2_set_raw),
    sigma_noise = sqrt(sigma2_noise),
    sigma2_set_true = pmax(0, sigma2_set_raw - sigma2_noise),
    sigma_set_true = sqrt(sigma2_set_true),
    ICC_set = sigma2_set_true / sigma2_set_raw) |>
  dplyr::left_join(cell_meta |> dplyr::select(form, K, n_combinations, enumerated),
                   by = c("form", "K")) |>
  dplyr::arrange(target, form, K)

# ---- D1 spot-check: fixed sets re-derived through the FROZEN engine ---------
# Per-VALUE equivalence is already proven at gate G_EQ (< 1e-12 over 6000
# halves). This check closes the loop at the level of a whole cell's r_sb, over
# N_SPOT fixed sets rather than one.
#
# TOLERANCES DIFFER BY TARGET, for a measured reason.
# loss_avoid_adj is bit-identical: its arithmetic is the same in both paths, so
# every split's Spearman agrees EXACTLY and the tolerance is 0. ev_sens_agn is
# not, and cannot be: its per-subject values are built from integer choice
# counts, so subjects hold EXACTLY TIED values (measured: ~3.2 of 63 per
# split). split_correlations ranks those values, so a difference of order 1e-16
# in the last bit resolves a tie the other way and moves that split's Spearman.
# This is rank-tie sensitivity of the FROZEN estimator itself, not a defect of
# the vectorised path; the diagnostic below records the tie counts and the
# number of splits affected so the magnitude is auditable. Measured over 12
# sets at the gate cell: median |d r_sb| 4.0e-5, max 6.1e-4 - two orders of
# magnitude below the 0.05 H1 threshold, and it cannot move the verdict.
say("\n=== D1 spot-check: frozen estimate_half re-derivation of %d fixed sets ===\n",
    N_SPOT)
fm1 <- form_mats[[GATE_FORM]]
spot_sets <- draw_sets(length(fm1$probs), GATE_K, min(M_CAP,
                       choose(length(fm1$probs), GATE_K)),
                       seed_subset(GATE_FORM, GATE_K))
spot_jobs <- lapply(seq_len(N_SPOT), function(i)
  list(set_id = i, cols = spot_sets[[i]],
       seed = seed_split(GATE_FORM, GATE_K) + i))
parallel::clusterExport(cl, c("subj_data", "fm1", "POS", "REGIME", "GATE_K",
                              "run_set"), envir = environment())
spot_detail <- dplyr::bind_rows(parallel::parLapply(cl, spot_jobs, function(job) {
  cols <- job$cols
  keep <- which(!is.na(rowSums(fm1$NB[, cols, drop = FALSE])))
  subs <- fm1$subs[keep]
  kh <- GATE_K %/% 2L
  Ef <- list(); for (tg in TARGETS) Ef[[tg]] <- list(
    h1 = matrix(NA_real_, length(subs), N_SPLIT_A),
    h2 = matrix(NA_real_, length(subs), N_SPLIT_A))
  set.seed(job$seed)
  ties <- integer(N_SPLIT_A)
  for (s in seq_len(N_SPLIT_A)) {
    take <- sample.int(GATE_K, 2L * kh)
    i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
    for (u in seq_along(subs)) {
      sd <- subj_data[[subs[u]]]
      r1 <- estimate_half(sd, match(fm1$probs[cols[i1]], sd$feat$GameID),
                          POS, REGIME)                       # FROZEN
      r2 <- estimate_half(sd, match(fm1$probs[cols[i2]], sd$feat$GameID),
                          POS, REGIME)                       # FROZEN
      for (tg in TARGETS) { Ef[[tg]]$h1[u, s] <- r1[[tg]]; Ef[[tg]]$h2[u, s] <- r2[[tg]] }
    }
    v <- Ef$ev_sens_agn$h1[, s]; v <- v[is.finite(v)]
    ties[s] <- length(v) - length(unique(v))
  }
  fast <- run_set(fm1, cols, N_SPLIT_A, job$seed)
  purrr::map_dfr(TARGETS, function(tg) {
    rf <- split_correlations(Ef[[tg]]$h1, Ef[[tg]]$h2, MIN_PAIRS)
    rv <- run_set(fm1, cols, N_SPLIT_A, job$seed, want_est = TRUE)
    rv <- split_correlations(rv$E[[tg]]$h1, rv$E[[tg]]$h2, MIN_PAIRS)
    frozen_sb <- sb_correct(mean(rf, na.rm = TRUE))
    tibble::tibble(
      set_id = job$set_id, target = tg, frozen_r_sb = frozen_sb,
      vectorised_r_sb = fast[[tg]]$r_sb,
      abs_diff = abs(frozen_sb - fast[[tg]]$r_sb),
      n_splits_differing = sum(abs(rf - rv) > 0, na.rm = TRUE),
      max_split_abs_diff = max(abs(rf - rv), na.rm = TRUE),
      mean_tied_subjects_per_split = mean(ties))
  })
}))
spot_rows <- spot_detail |>
  dplyr::group_by(target) |>
  dplyr::summarise(n_sets = dplyr::n(),
                   max_abs_diff = max(abs_diff),
                   median_abs_diff = stats::median(abs_diff),
                   mean_splits_differing = mean(n_splits_differing),
                   mean_tied_subjects_per_split = mean(mean_tied_subjects_per_split),
                   .groups = "drop")
print(as.data.frame(spot_rows))
spot_loss <- spot_rows$max_abs_diff[spot_rows$target == "loss_avoid_adj"]
spot_ev   <- spot_rows$max_abs_diff[spot_rows$target == "ev_sens_agn"]
stopifnot(spot_loss == 0)                    # primary target: bit-identical
stopifnot(spot_ev < SPOT_TOL_EV)             # secondary: tie-resolution only
say("D1 spot-check PASS: loss_avoid_adj max |d r_sb| = %.3g (exact); ev_sens_agn max = %.3g (< %.0e, tie resolution)\n",
    spot_loss, spot_ev, SPOT_TOL_EV)

# ---- H1 verdict at the primary gate cell ------------------------------------
gate_row <- heterogeneity_summary |>
  dplyr::filter(form == GATE_FORM, K == GATE_K, target == GATE_TARGET)
stopifnot(nrow(gate_row) == 1L)
contingency_applied <- gate_row$set_coverage < COVERAGE_CONTINGENCY
gate_target_used <- if (contingency_applied) "ev_sens_agn" else GATE_TARGET
gate_final <- heterogeneity_summary |>
  dplyr::filter(form == GATE_FORM, K == GATE_K, target == gate_target_used)
h1_pass <- isTRUE(gate_final$sigma_set_true >= H1_SIGMA_MIN) &&
  isTRUE(gate_final$ICC_set >= H1_ICC_MIN)
say("\n=== H1 GATE (form %d, K %d, target %s) ===\n",
    GATE_FORM, GATE_K, gate_target_used)
say("set coverage %.4f (contingency threshold %.2f, applied = %s)\n",
    gate_row$set_coverage, COVERAGE_CONTINGENCY, contingency_applied)
say("sigma_set_raw %.5f | sigma_noise %.5f | sigma_set_true %.5f (>= %.2f?) | ICC_set %.4f (>= %.2f?)\n",
    gate_final$sigma_set_raw, gate_final$sigma_noise, gate_final$sigma_set_true,
    H1_SIGMA_MIN, gate_final$ICC_set, H1_ICC_MIN)
say("H1 VERDICT: %s\n", if (h1_pass) "PASS" else "FAIL")

# =============================================================================
# 6. Part B: held-out choice prediction across sets (real data, build folds)
# =============================================================================
say("\n=== Part B: held-out prediction (frozen fit_subject_firth) ===\n")
fmB <- form_mats[[GATE_FORM]]
predB_cols <- list()
set.seed(SEED_PRED)
for (K in c(8L, 12L)) {
  allsets <- draw_sets(length(fmB$probs), K, min(M_CAP,
                       choose(length(fmB$probs), K)), seed_subset(GATE_FORM, K))
  idx <- sort(sample.int(length(allsets), min(N_PRED_SETS, length(allsets))))
  predB_cols[[as.character(K)]] <- list(idx = idx, sets = allsets[idx])
}
# Per-subject per-trial frames for Set 1 build subjects (frozen column names).
b_subs <- fmB$subs
b_trials <- lapply(b_subs, function(sid) {
  sd <- subj_data[[sid]]
  j <- match(fmB$probs, sd$feat$GameID)
  list(j = j, feat = sd$feat[j, ], B = sd$B_mat[j, POS, drop = FALSE])
})
names(b_trials) <- b_subs

invisible(parallel::clusterEvalQ(cl, {
  source(here("R", "simulate_recovery.R")); NULL
}))
parallel::clusterExport(cl, c("fmB", "b_trials", "b_subs", "R_REPS", "POS"),
                        envir = environment())
predB_jobs <- unlist(lapply(names(predB_cols), function(K)
  lapply(seq_along(predB_cols[[K]]$sets), function(i)
    list(K = as.integer(K), set_id = predB_cols[[K]]$idx[i],
         cols = predB_cols[[K]]$sets[[i]]))), recursive = FALSE)
predB <- parallel::parLapply(cl, predB_jobs, function(job) {
  cols <- job$cols
  hold <- setdiff(seq_along(fmB$probs), cols)
  ll <- numeric(0); br <- numeric(0); nfail <- 0L; nclamp <- 0L
  for (sid in b_subs) {
    bt <- b_trials[[sid]]
    if (anyNA(bt$B[cols, ])) next
    ftr <- bt$feat[cols, ]
    ch <- data.frame(
      choice = as.vector(t(bt$B[cols, , drop = FALSE])),
      z_ev = rep(ftr$z_ev, each = R_REPS),
      z_sd = rep(ftr$z_sd, each = R_REPS),
      z_eloss = rep(ftr$z_eloss, each = R_REPS))
    fit <- fit_subject_firth(ch)
    if (isTRUE(fit$degenerate) || !isTRUE(fit$converged) ||
        anyNA(c(fit$est_b0, fit$est_b_ev, fit$est_b_risk, fit$est_b_loss))) {
      nfail <- nfail + 1L; next
    }
    fh <- bt$feat[hold, ]
    p <- stats::plogis(fit$est_b0 + fit$est_b_ev * fh$z_ev +
                         fit$est_b_risk * fh$z_sd + fit$est_b_loss * fh$z_eloss)
    nB <- rowSums(bt$B[hold, , drop = FALSE])
    okh <- !is.na(nB)
    p <- p[okh]; nB <- nB[okh]
    nclamp <- nclamp + sum(p < 1e-15 | p > 1 - 1e-15)
    pc <- pmin(pmax(p, 1e-15), 1 - 1e-15)
    n_tr <- R_REPS * length(pc)
    ll <- c(ll, -(sum(nB * log(pc)) + sum((R_REPS - nB) * log(1 - pc))) / n_tr)
    br <- c(br, (sum(nB * (1 - p)^2) + sum((R_REPS - nB) * p^2)) / n_tr)
  }
  # A set can be UNEVALUABLE: with 4 parameters and no loss-informative
  # problem, z_eloss is constant, the per-subject design is rank-deficient and
  # every subject's fit is degenerate. That is a result about the set, not an
  # error, so it is recorded as NA with the diagnosis attached rather than
  # allowed to propagate as NaN through the summary.
  tibble::tibble(form = 1L, K = job$K, set_id = job$set_id,
                 game_ids = paste(fmB$probs[cols], collapse = ";"),
                 n_loss_informative = sum(fmB$li[cols]),
                 n_held_out_problems = length(hold),
                 n_subjects_fit = length(ll), n_fit_failures = nfail,
                 n_clamped_predictions = nclamp,
                 mean_log_loss = if (length(ll)) mean(ll) else NA_real_,
                 mean_brier = if (length(br)) mean(br) else NA_real_)
})
prediction_by_set_real <- dplyr::bind_rows(predB) |> dplyr::arrange(K, set_id)
pred_summary <- prediction_by_set_real |>
  dplyr::group_by(K) |>
  dplyr::summarise(
    n_sets = dplyr::n(),
    n_sets_evaluable = sum(!is.na(mean_log_loss)),
    n_sets_no_subject_fit = sum(is.na(mean_log_loss)),
    n_sets_zero_loss_informative = sum(n_loss_informative == 0),
    ll_min = min(mean_log_loss, na.rm = TRUE),
    ll_median = stats::median(mean_log_loss, na.rm = TRUE),
    ll_max = max(mean_log_loss, na.rm = TRUE),
    ll_best_minus_median = stats::median(mean_log_loss, na.rm = TRUE) -
      min(mean_log_loss, na.rm = TRUE),
    br_min = min(mean_brier, na.rm = TRUE),
    br_median = stats::median(mean_brier, na.rm = TRUE),
    br_max = max(mean_brier, na.rm = TRUE),
    br_best_minus_median = stats::median(mean_brier, na.rm = TRUE) -
      min(mean_brier, na.rm = TRUE),
    .groups = "drop")
print(as.data.frame(pred_summary))

# =============================================================================
# 7. Part C: SIMULATION free-composition reliability + recovery
#    Every number in this section is SIMULATED under the frozen generative
#    model (R/simulate_recovery.R) and the frozen uniform-range population. It is never an empirical
#    result. The population is not calibrated to the build subjects here.
# =============================================================================
say("\n=== Part C: SIMULATION (free composition across the 105 universe) ===\n")
uni_feat <- feat_all
uni_fm <- list(probs = uni_feat$GameID, feat = uni_feat,
               li = uni_feat$loss_informative,
               has_loss_b = uni_feat$has_loss_b, d_ev = uni_feat$d_ev,
               p_base = uni_feat$p_base_exp)

# The r -> infinity limit of each frozen index under the noise-free choice
# probabilities: the same formulas evaluated at nB_j = r * analytic_p_b.
sim_one_set <- function(cols, pars, seed, n_split) {
  set.seed(seed)
  f <- uni_feat[cols, ]
  K <- length(cols); kh <- K %/% 2L
  ns <- nrow(pars)
  NB <- matrix(NA_real_, ns, K)
  for (i in seq_len(ns)) {
    ch <- simulate_choices(pars[i, ], f, N_REPS_SIM)            # FROZEN
    NB[i, ] <- rowSums(matrix(ch$choice, ncol = N_REPS_SIM, byrow = TRUE))
  }
  PT <- t(vapply(seq_len(ns), function(i)
    R_REPS * analytic_p_b(pars[i, ], f), numeric(K)))           # FROZEN
  E <- list()
  for (tg in TARGETS) E[[tg]] <- list(h1 = matrix(NA_real_, ns, n_split),
                                      h2 = matrix(NA_real_, ns, n_split))
  for (s in seq_len(n_split)) {
    take <- sample.int(K, 2L * kh)
    i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
    a <- half_targets(NB[, i1, drop = FALSE], uni_fm, cols[i1], R_REPS)
    b <- half_targets(NB[, i2, drop = FALSE], uni_fm, cols[i2], R_REPS)
    for (tg in TARGETS) { E[[tg]]$h1[, s] <- a[[tg]]; E[[tg]]$h2[, s] <- b[[tg]] }
  }
  est_full  <- half_targets(NB, uni_fm, cols, R_REPS)
  true_full <- half_targets(PT, uni_fm, cols, R_REPS)
  purrr::map_dfr(TARGETS, function(tg) {
    rs <- split_correlations(E[[tg]]$h1, E[[tg]]$h2, MIN_PAIRS)  # FROZEN
    r_raw <- mean(rs, na.rm = TRUE)
    x <- est_full[[tg]]; y <- true_full[[tg]]
    ok <- is.finite(x) & is.finite(y)
    rec <- if (sum(ok) < MIN_PAIRS || stats::sd(x[ok]) == 0 ||
               stats::sd(y[ok]) == 0) NA_real_ else
                 suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
    rr <- if (sum(ok) < MIN_PAIRS) NA_real_ else
      sqrt(mean((rank(x[ok]) - rank(y[ok]))^2))
    tibble::tibble(
      target = tg, sim_r_raw = r_raw,
      sim_r_sb = if (is.finite(r_raw)) sb_correct(r_raw) else NA_real_,
      sim_coverage = mean(!is.na(rs)),
      recovery_spearman = rec, recovery_rank_rmse = rr,
      recovery_rank_rmse_norm = rr / sum(ok), n_recovery_pairs = sum(ok))
  })
}

parallel::clusterExport(cl, c("uni_feat", "uni_fm", "sim_one_set", "N_REPS_SIM",
                              "N_SPLIT_C", "half_targets", "ev_ranks",
                              "TARGETS", "R_REPS", "MIN_PAIRS"),
                        envir = environment())

sim_rows <- list(); anchor_rows <- list()
for (K in K_GRID) {
  set.seed(seed_simpar(K))
  pars <- draw_true_params(N_SIM_SUBJ)                          # FROZEN
  # (a) free composition across the whole 105 universe
  free_sets <- draw_sets(nrow(uni_feat), K, N_SETS_SIM, seed_simset(K))
  free_sets <- free_sets[seq_len(min(N_SETS_SIM, length(free_sets)))]
  # (b) anchor: the FIRST N_SETS_SIM Part-A Set-1 sets, so real and simulated
  #     r_sb are paired on identical sets
  a_sets_cols <- draw_sets(length(fm1$probs), K,
                           min(M_CAP, choose(length(fm1$probs), K)),
                           seed_subset(GATE_FORM, K))
  a_idx <- seq_len(min(N_SETS_SIM, length(a_sets_cols)))
  anchor_cols <- lapply(a_sets_cols[a_idx],
                        function(z) match(fm1$probs[z], uni_feat$GameID))
  # Seeds are resolved HERE and carried in the job list, so workers need no
  # seed function and the stream is identical under any execution order.
  free_jobs <- lapply(seq_along(free_sets), function(i)
    list(cols = free_sets[[i]], seed = seed_simset(K) + 10000L + i))
  anch_jobs <- lapply(seq_along(anchor_cols), function(i)
    list(cols = anchor_cols[[i]], seed = seed_simset(K) + 20000L + i))
  parallel::clusterExport(cl, c("pars"), envir = environment())
  t0 <- Sys.time()
  free_res <- parallel::parLapply(cl, free_jobs, function(j)
    sim_one_set(j$cols, pars, j$seed, N_SPLIT_C))
  anch_res <- parallel::parLapply(cl, anch_jobs, function(j)
    sim_one_set(j$cols, pars, j$seed, N_SPLIT_C))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sim_rows[[length(sim_rows) + 1L]] <- dplyr::bind_rows(
    purrr::imap_dfr(free_res, function(d, i) d |> dplyr::mutate(
      composition = "free_105", K = K, set_id = i,
      game_ids = paste(uni_feat$GameID[free_sets[[i]]], collapse = ";"),
      n_loss_informative = sum(uni_feat$loss_informative[free_sets[[i]]]))),
    purrr::imap_dfr(anch_res, function(d, i) d |> dplyr::mutate(
      composition = "anchor_set1", K = K, set_id = a_idx[i],
      game_ids = paste(uni_feat$GameID[anchor_cols[[i]]], collapse = ";"),
      n_loss_informative = sum(uni_feat$loss_informative[anchor_cols[[i]]]))))
  say("  K %2d: %d free + %d anchor sets, %.1f s\n",
      K, length(free_res), length(anch_res), el)
}
reliability_recovery_by_set_sim <- dplyr::bind_rows(sim_rows) |>
  dplyr::select(composition, K, set_id, target, game_ids, n_loss_informative,
                dplyr::everything()) |>
  dplyr::arrange(composition, target, K, set_id)

# ---- Set-1 sim-vs-real anchor ----------------------------------------------
real_s1 <- reliability_by_set_real |>
  dplyr::filter(form == GATE_FORM) |>
  dplyr::select(K, target, set_id, real_r_sb = r_sb, real_coverage = coverage)
sim_a <- reliability_recovery_by_set_sim |>
  dplyr::filter(composition == "anchor_set1") |>
  dplyr::select(K, target, set_id, sim_r_sb, sim_coverage,
                recovery_spearman)
sim_real_anchor_pairs <- dplyr::inner_join(real_s1, sim_a,
                                           by = c("K", "target", "set_id"))
sim_real_anchor <- sim_real_anchor_pairs |>
  dplyr::group_by(K, target) |>
  dplyr::summarise(
    n_shared_sets = dplyr::n(),
    n_pairs_used = sum(is.finite(real_r_sb) & is.finite(sim_r_sb)),
    spearman_real_vs_sim = {
      ok <- is.finite(real_r_sb) & is.finite(sim_r_sb)
      if (sum(ok) < MIN_PAIRS || stats::sd(real_r_sb[ok]) == 0 ||
          stats::sd(sim_r_sb[ok]) == 0) NA_real_ else
            suppressWarnings(stats::cor(real_r_sb[ok], sim_r_sb[ok],
                                        method = "spearman"))
    },
    sigma_set_real = stats::sd(real_r_sb, na.rm = TRUE),
    sigma_set_sim = stats::sd(sim_r_sb, na.rm = TRUE),
    median_real_r_sb = stats::median(real_r_sb, na.rm = TRUE),
    median_sim_r_sb = stats::median(sim_r_sb, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::arrange(target, K)
print(sim_real_anchor)

parallel::stopCluster(cl)

# =============================================================================
# 8. Outputs
# =============================================================================
save_csv(reliability_by_set_real, "reliability_by_set_real")
save_csv(heterogeneity_summary, "heterogeneity_summary")
save_csv(prediction_by_set_real, "prediction_by_set_real")
save_csv(reliability_recovery_by_set_sim, "reliability_recovery_by_set_sim")
save_csv(sim_real_anchor, "sim_real_anchor")
save_csv(form_struct, "e3_form_structure")
save_csv(noise_tbl, "e3_noise_floor_by_set")
save_csv(spot_detail, "e3_d1_spotcheck")
save_csv(spot_rows, "e3_d1_spotcheck_summary")

sim_summary <- reliability_recovery_by_set_sim |>
  dplyr::group_by(composition, K, target) |>
  dplyr::summarise(
    n_sets = dplyr::n(),
    sim_r_sb_min = min(sim_r_sb, na.rm = TRUE),
    sim_r_sb_median = stats::median(sim_r_sb, na.rm = TRUE),
    sim_r_sb_max = max(sim_r_sb, na.rm = TRUE),
    sigma_set_sim = stats::sd(sim_r_sb, na.rm = TRUE),
    recovery_min = min(recovery_spearman, na.rm = TRUE),
    recovery_median = stats::median(recovery_spearman, na.rm = TRUE),
    recovery_max = max(recovery_spearman, na.rm = TRUE),
    rank_rmse_median = stats::median(recovery_rank_rmse, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::arrange(composition, target, K)
save_csv(sim_summary, "e3_simulation_summary")
save_csv(pred_summary, "e3_prediction_summary")

h1_gate <- list(
  spec = "elevation/E3_spec.md",
  gate_cell = list(form = GATE_FORM, K = GATE_K,
                   target_predeclared = GATE_TARGET,
                   target_used = gate_target_used,
                   contingency_applied = contingency_applied,
                   contingency_rule = sprintf(
                     "fall back to ev_sens_agn if loss_avoid_adj identifiable-set coverage < %.2f",
                     COVERAGE_CONTINGENCY),
                   set_coverage_loss_avoid_adj = gate_row$set_coverage,
                   split_coverage_loss_avoid_adj = gate_row$split_coverage),
  thresholds = list(sigma_set_true_min = H1_SIGMA_MIN, ICC_set_min = H1_ICC_MIN),
  observed = list(
    M = gate_final$M, n_subjects_median = gate_final$n_subjects_median,
    sigma_set_raw = gate_final$sigma_set_raw,
    sigma_noise = gate_final$sigma_noise,
    sigma_set_true = gate_final$sigma_set_true,
    ICC_set = gate_final$ICC_set,
    r_sb_min = gate_final$r_sb_min, r_sb_median = gate_final$r_sb_median,
    r_sb_max = gate_final$r_sb_max,
    best_minus_median = gate_final$best_minus_median),
  h1_verdict = if (h1_pass) "PASS" else "FAIL",
  h1_pass = h1_pass,
  full_map = list(
    n_cells = nrow(heterogeneity_summary),
    n_cells_pass_both = sum(heterogeneity_summary$sigma_set_true >= H1_SIGMA_MIN &
                              heterogeneity_summary$ICC_set >= H1_ICC_MIN,
                            na.rm = TRUE),
    max_sigma_set_true = max(heterogeneity_summary$sigma_set_true, na.rm = TRUE),
    max_ICC_set = max(heterogeneity_summary$ICC_set, na.rm = TRUE)),
  d1_spotcheck_max_abs_diff_loss_avoid_adj = spot_loss,
  d1_spotcheck_max_abs_diff_ev_sens_agn = spot_ev,
  g_eq_max_abs_diff = max(eq_max_ev, eq_max_loss),
  build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
  optimization_performed = FALSE, mechanism_model_fitted = FALSE)
# digits = NA writes full double precision, so the gate quantities can be
# re-derived exactly rather than from 4-significant-digit rounding.
jsonlite::write_json(h1_gate, file.path(res_dir, "h1_gate.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

# Seed vectors are built BEFORE the tibble call: inside tibble() the names
# `form` and `K` resolve to the columns just created, which get recycled to the
# length of `purpose`, so calling seed_subset(S, K) in-place would silently
# produce a length-3 seed vector.
seeds_tbl <- dplyr::bind_rows(
  purrr::map_dfr(seq_len(nrow(cellsA)), function(i) {
    S <- cellsA$form[i]; Kc <- cellsA$K[i]
    sv <- c(seed_subset(S, Kc), seed_split(S, Kc), seed_noise(S, Kc, 0L))
    stopifnot(length(sv) == 3L)
    tibble::tibble(
      part = "A", form = S, K = Kc,
      purpose = c("subset_draw", "split_draw_base", "noise_bootstrap_base"),
      seed = sv,
      derivation = c("SEED_E3 + 1000*form + K",
                     "SEED_E3 + 1000*form + K + 500000; per-set seed = base + set_id",
                     "SEED_E3 + 600000 + 1000*form + K; per-set seed = base + set_id"))
  }),
  tibble::tibble(part = "B", form = GATE_FORM, K = NA_integer_,
                 purpose = "prediction_subsample", seed = SEED_PRED,
                 derivation = "SEED_E3 + 900000; K = 8 then K = 12 from one stream"),
  purrr::map_dfr(K_GRID, function(Kc) {
    sv <- c(seed_simset(Kc), seed_simpar(Kc), seed_simset(Kc) + 10000L,
            seed_simset(Kc) + 20000L)
    stopifnot(length(sv) == 4L)
    tibble::tibble(
      part = "C", form = NA_integer_, K = Kc,
      purpose = c("sim_set_draw", "sim_param_draw", "sim_free_rep_base",
                  "sim_anchor_rep_base"),
      seed = sv,
      derivation = c("SEED_E3 + 700000 + K", "SEED_E3 + 800000 + K",
                     "SEED_E3 + 700000 + K + 10000; per-set seed = base + set_id",
                     "SEED_E3 + 700000 + K + 20000; per-set seed = base + set_id"))
  }))
save_csv(seeds_tbl, "e3_seeds")

cfg <- list(
  spec = "elevation/E3_spec.md",
  raw_md5 = raw_md5, seed_base = SEED_E3,
  regime = "experience (Trial 6-25)", r_reps = R_REPS,
  targets = list(primary = "loss_avoid_adj", secondary = "ev_sens_agn"),
  K_grid = K_GRID, M_cap = M_CAP, n_split_partA = N_SPLIT_A,
  n_noise_sets = N_NOISE_SETS, n_boot = N_BOOT,
  n_prediction_sets = N_PRED_SETS, n_sets_sim = N_SETS_SIM,
  n_sim_subjects = N_SIM_SUBJ, n_reps_sim = N_REPS_SIM,
  n_split_partC = N_SPLIT_C, min_pairs = MIN_PAIRS,
  h1_thresholds = list(sigma_set_true_min = H1_SIGMA_MIN,
                       ICC_set_min = H1_ICC_MIN),
  gate_cell = list(form = GATE_FORM, K = GATE_K, target = GATE_TARGET),
  split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
  loss_baselines_exp = as.list(baselines$exp),
  halving_convention = paste(
    "The halving of a set S is COMMON to all subjects within a split (S is a",
    "fixed instrument), per E3_spec section 2. This differs deliberately from",
    "the Study-1 run_cell, where each subject received an independent partition",
    "because subjects there held different problem sets."),
  eligibility_rule = paste(
    "Eligible subjects for a set S are the form's build subjects who saw EVERY",
    "problem in S. For Sets 1-5 this is the whole form (one distinct clean set",
    "each, spec premise holds). For Sets 6-7 the spec premise is FALSE (7 and 6",
    "distinct clean sets), so the rule binds and n_subjects is reported per set."),
  sigma_noise_derivation = paste(
    "Per-set estimation SD is the normal-approximation SD of the FROZEN",
    "reliability_point percentile CI on the SB scale,",
    "(ci_hi_sb - ci_lo_sb) / (2 * qnorm(0.975)); sigma2_noise is the mean of",
    "its square over the 50-set subsample."),
  accelerator = paste(
    "The two E3 targets are evaluated by a vectorised path using the SAME",
    "formulas as the frozen estimate_half, because estimate_half also runs the",
    "4-parameter Firth fit for b_ev/b_risk/b_loss, which are not E3 targets;",
    "the specified Part A would otherwise need ~3.2e8 engine calls. Gate G_EQ",
    "and the D1 spot-check verify agreement against the frozen engine."),
  g_eq_checks = eq_n, g_eq_max_abs_diff = max(eq_max_ev, eq_max_loss),
  g_eq_na_mismatches = eq_na,
  d1_spotcheck_max_abs_diff_loss_avoid_adj = spot_loss,
  d1_spotcheck_max_abs_diff_ev_sens_agn = spot_ev,
  data_inputs = c("data/raw/All estimation raw data.csv",
                  "elevation/results/e1/participant_split.csv",
                  "elevation/results/e2/problem_information_dictionary.csv"),
  no_freeze_or_test_choice_read = TRUE,
  optimization_performed = FALSE, mechanism_model_fitted = FALSE,
  part_c_is_simulation = TRUE,
  runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
jsonlite::write_json(cfg, file.path(res_dir, "e3_config.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

say("\nE3 complete in %.1f min. Outputs in elevation/results/e3/\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))
say("H1 (%s, form %d, K %d): sigma_set_true %.5f vs %.2f; ICC_set %.4f vs %.2f -> %s\n",
    gate_target_used, GATE_FORM, GATE_K, gate_final$sigma_set_true,
    H1_SIGMA_MIN, gate_final$ICC_set, H1_ICC_MIN,
    if (h1_pass) "PASS" else "FAIL")

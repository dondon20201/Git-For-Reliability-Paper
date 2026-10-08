# Step 29: detectability of the mechanism and the attenuation bound (simulation)
# Reads:   output/elevation/results/e1/, e2/, e3/, e3b/, e4/ (steps 20-25); definitions read from steps 22 and 25
# Writes:  output/elevation/results/sens/ (attenuation_bound.csv, detection and robustness tables, configuration)
# Paper:   Study 2; Supplement (Mechanism sensitivity and the attenuation bound), through step 36
# Runtime: about 7.3 hours with 12 parallel workers

# SENS: simulation calibration of mechanism detectability and attenuation.
# Every planted-mechanism result is labeled SIMULATION. The H2 FAIL from
# 25_mechanism_model.R remains unchanged; this script neither fits a new
# empirical mechanism nor optimizes or constructs candidate designs.
#
# Inputs are build-only outputs from E2 (step 21), E3 (22), E3b (24), and E4
# (25), plus the frozen simulation and scoring functions. No raw participant
# data is loaded. The step 20 split file is read only to check the 344/137/205
# group counts; freeze and test choices are excluded.
#
# R/step2_reliability.R and R/simulate_recovery.R are sourced unchanged.
# The needed functions and constants from steps 22 and 25 are extracted by
# parsing and evaluating their named top-level definitions (section 0),
# preserving the original function bodies.
#
# Before any planted-data gate is evaluated, G_E4EQ reproduces the committed
# E4 cv_performance and pooled statistics on the real form-1 inputs with E4's
# seeds. G_EQ checks the simulated scoring path against estimate_half. G_GEN
# checks the planted generator against analytic_p_b under the declared
# parameter substitution and against the original generator at beta = 0.
#
# The frozen E4 gate uses repeated set-level 70/30 holdout within K and a
# problem-level cluster bootstrap. Its empirical anchor is form 1; there is no
# leave-one-form-out fold. The planted data follow that same design: 63
# subjects, 24 clean problems, and r = 20 experience trials.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})

source(here("R", "step2_reliability.R"))
source(here("R", "simulate_recovery.R"))

res_dir <- here("output", "elevation", "results", "sens")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  invisible(df)
}
say <- function(...) cat(sprintf(...), sep = "")
t_start <- Sys.time()

# =============================================================================
# 0. Frozen-definition extractor
# =============================================================================
# Parses a committed frozen script and evaluates EXACTLY the requested
# top-level `name <- value` definitions, in file order, into the global
# environment. The evaluated bodies are byte-identical to the frozen file (no
# retyping); nothing else in the file is executed. This is the sourcing
# mechanism for engines that live inside committed analysis scripts rather
# than R/ libraries.
extract_frozen <- function(path, names_wanted) {
  exprs <- parse(file = path, keep.source = FALSE)
  found <- character(0)
  for (e in exprs) {
    if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
        is.symbol(e[[2]]) && as.character(e[[2]]) %in% names_wanted) {
      eval(e, envir = globalenv())
      found <- c(found, as.character(e[[2]]))
    }
  }
  missing <- setdiff(names_wanted, found)
  if (length(missing)) stop("frozen definitions not found in ", path, ": ",
                            paste(missing, collapse = ", "))
  tibble::tibble(source_file = path, definition = found)
}

p_e3 <- here("steps", "22_problem_subset_reliability.R")
p_e4 <- here("steps", "25_mechanism_model.R")

prov_e3 <- extract_frozen(p_e3, c("REGIME", "POS", "R_REPS", "TARGETS",
                                  "MIN_PAIRS", "ev_ranks", "half_targets",
                                  "run_set"))
TARGETS_E3 <- TARGETS
prov_e4 <- extract_frozen(p_e4, c("SEED_E4", "N_REP", "TRAIN_FRAC", "N_BOOT",
                                  "RIDGE_GRID", "PLS_NCOMP", "INNER_FOLDS",
                                  "TARGETS", "K_GRID", "K_GATE", "FEATURES",
                                  "MODELS", "GATE_SPEARMAN_MIN",
                                  "GATE_BOOT_SUPPORT", "MIN_CELL_N",
                                  "CARRIER_MIN_BETA",
                                  "abs_cor", "ridge_fit", "simpls_fit",
                                  "inner_folds", "safe_lsq", "fit_predict",
                                  "sp", "cv_cell", "set_problem_matrix",
                                  "bootstrap_cell_group", "ci_lo", "ci_hi"))
stopifnot(identical(TARGETS, TARGETS_E3))          # same frozen target pair
stopifnot(identical(K_GATE, c(8L, 10L, 12L)))
stopifnot(GATE_SPEARMAN_MIN == 0.30, GATE_BOOT_SUPPORT == 0.90)
stopifnot(N_BOOT == 2000L, N_REP == 200L, R_REPS == 20L, MIN_PAIRS == 10L)
PRIMARY <- TARGETS[1]                              # loss_avoid_adj
stopifnot(identical(PRIMARY, "loss_avoid_adj"))
frozen_provenance <- dplyr::bind_rows(prov_e3, prov_e4)
say("frozen definitions extracted: %d from E3, %d from E4\n",
    nrow(prov_e3), nrow(prov_e4))

# =============================================================================
# 1. SENS constants (fixed HERE, before any result exists)
# =============================================================================
SEED_SENS   <- 20260722L
BETA_GRID   <- c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40)   # FIXED
N_SIM_MAIN  <- 100L      # replicates per beta, main detection curve
N_SIM_COLL  <- 40L       # replicates per beta per collinearity scenario
N_SIM_NEG   <- 100L      # negative-control replicates (beta = 0.40, carrier
                         # orthogonalized against the true one)
N_SIM_RATE  <- 40L       # per beta in {0, 0.15}: uncalibrated-rate sensitivity
COLL_TARGETS <- c(0.5, 0.7, 0.9, 0.95)                     # FIXED
BETA_NEG    <- 0.40
BETA_POS    <- 0.40      # positive control = top of the grid
K_SENS      <- K_GATE    # gate budgets; K = 6 is not part of the frozen gate
N_SPLIT_SENS <- 100L     # E3's committed N_SPLIT_A: the empirical split noise
M_SETS      <- 2000L     # the committed E3 form-1 sets per K, used in full
N_HP_SUBJ   <- 2000L     # high-precision arm subjects (near-true reliability)
N_HP_REPS   <- 2L        # high-precision replicates per beta
N_SPLIT_HP  <- 40L
M_SETS_HP_812 <- 600L    # K in {8,10}: first 600 committed sets in the hp arm
M_SETS_FORM <- 200L      # per form in the seven-form generalization arm
N_CAL       <- 100000L   # frozen-population draws for the rate calibration
N_GEQ       <- 2000L     # scorer-equivalence checks vs frozen estimate_half
N_GGEN      <- 4000L     # generator-equivalence checks vs frozen analytic_p_b
E4EQ_TOL    <- 1e-10
MU_L        <- -0.75     # mean of the FROZEN draw_true_params b_loss range
                         # U(-2.0, 0.5); asserted against a large draw below
QUANT_GRAD  <- 0.20      # top/bottom quintile for the induced-change gradient
MATTERS_SB  <- 0.05      # pre-declared "matters" threshold at K = 12
POWER_BAR   <- 0.80      # G-power / G-detect-floor power bar
NEG_MAX_P   <- 0.20      # G-controls: negative control "not detected" bound
N_WORKERS   <- 12L

# Detection-event definitions (FIXED before results; the frozen-gate clause
# structure applied to the planted track, which is the only track a planted
# world has; the G-real "sign consistency with the simulation" clause maps to
# sign consistency with the PLANTED truth, which speaks only through F4):
#  C1 (G-sim structure): gate statistic = min(pooled, repeat-mean) out-of-fold
#     Spearman of FULL >= 0.30 at EVERY gate K with the cluster-bootstrap CI
#     lower bound > 0, AND the pooled-across-gate-K statistic >= 0.30 with
#     pooled CI lower bound > 0.
#  C2 (G-beyond): at EVERY gate K, (delta vs BASE_category > 0 in Spearman OR
#     in R^2) AND delta bootstrap support >= 0.90.
#  C3 (G-real sign structure vs planted truth): at EVERY gate K, F4 is a
#     carrier (|standardized full-sample FULL beta| >= CARRIER_MIN_BETA) AND
#     its coefficient is positive (the planted direction).
#  P(gate) = P(C1 and C2 and C3).
#  P(threshold) = P(gate statistic >= 0.30 at every gate K AND pooled), the
#     literal "out-of-fold Spearman exceeds the frozen E4 threshold".
#  P(sign) = P(F4 coefficient > 0 at every gate K).
# Controls (G-controls, FIXED): positive control detected iff
# P(gate | beta = 0.40, true carrier) >= 0.80; negative control not detected
# iff P(gate | beta = 0.40, orthogonalized carrier) <= 0.20. beta = 0 rows are
# the baseline (composition-channel) rate and are reported as such.

# ---- seed derivations (all < 2^31 - 1; resolved in job lists) ---------------
seed_cal    <- SEED_SENS + 11L
seed_geq    <- SEED_SENS + 22L
seed_ggen   <- SEED_SENS + 33L
seed_par    <- function(arm, rep)      100000000L + arm * 1000000L + rep * 100L
seed_choice <- function(arm, rep)      100000000L + arm * 1000000L + rep * 100L + 1L
seed_negcar <- function(arm, rep)      100000000L + arm * 1000000L + rep * 100L + 2L
seed_split  <- function(arm, rep, K, set_id)
  300000000L + arm * 3000000L + rep * 30000L + ((K - 8L) %/% 2L) * 10000L + set_id
seed_cv     <- function(arm, rep, K, variant)
  700000000L + variant * 200000000L + arm * 3000000L + rep * 30000L + K
seed_boot_g <- function(arm, rep, variant)
  700000000L + variant * 200000000L + arm * 3000000L + rep * 30000L + 99L

# Arm registry (each arm = one generator configuration).
arm_tbl <- dplyr::bind_rows(
  tibble::tibble(arm = seq_along(BETA_GRID), family = "main",
                 beta = BETA_GRID, carrier = "canonical",
                 coll_target = NA_real_, calibrated = TRUE, n_rep = N_SIM_MAIN),
  tibble::tibble(arm = 8L, family = "negative_control", beta = BETA_NEG,
                 carrier = "orthogonalized_random", coll_target = NA_real_,
                 calibrated = TRUE, n_rep = N_SIM_NEG),
  tidyr::expand_grid(sc = seq_along(COLL_TARGETS), bi = seq_along(BETA_GRID)) |>
    dplyr::transmute(arm = 8L + (sc - 1L) * length(BETA_GRID) + bi,
                     family = "collinearity", beta = BETA_GRID[bi],
                     carrier = paste0("scenario_", COLL_TARGETS[sc]),
                     coll_target = COLL_TARGETS[sc], calibrated = TRUE,
                     n_rep = N_SIM_COLL),
  tibble::tibble(arm = 37L:38L, family = "rate_sensitivity",
                 beta = c(0, 0.15), carrier = "canonical",
                 coll_target = NA_real_, calibrated = FALSE, n_rep = N_SIM_RATE),
  tibble::tibble(arm = 40L + seq_along(BETA_GRID), family = "high_precision",
                 beta = BETA_GRID, carrier = "canonical",
                 coll_target = NA_real_, calibrated = TRUE, n_rep = N_HP_REPS))
stopifnot(!anyDuplicated(arm_tbl$arm))

# =============================================================================
# 2. Inputs: committed build-only outputs (utils::read.csv for float exactness)
# =============================================================================
p_dict <- here("output", "elevation", "results", "e2", "problem_information_dictionary.csv")
p_real <- here("output", "elevation", "results", "e3", "reliability_by_set_real.csv")
p_sf   <- here("output", "elevation", "results", "e4", "set_features.csv")
p_cv   <- here("output", "elevation", "results", "e4", "cv_performance.csv")
p_cvp  <- here("output", "elevation", "results", "e4", "cv_performance_pooled_K.csv")
p_orel <- here("output", "elevation", "results", "e3b", "outcome_reliability_real.csv")
p_e3cf <- here("output", "elevation", "results", "e3", "e3_config.json")
p_splt <- here("output", "elevation", "results", "e1", "participant_split.csv")

input_md5 <- vapply(c(e2_dictionary = p_dict, e3_real = p_real,
                      e4_set_features = p_sf, e4_cv = p_cv, e4_cv_pooled = p_cvp,
                      e3b_outcome_reliability = p_orel, e3_config = p_e3cf,
                      e1_split = p_splt),
                    function(f) unname(tools::md5sum(f)), character(1))

DICT <- utils::read.csv(p_dict, check.names = FALSE)
REAL <- utils::read.csv(p_real, check.names = FALSE)
SF   <- utils::read.csv(p_sf,   check.names = FALSE)
CVE4 <- utils::read.csv(p_cv,   check.names = FALSE)
CVP4 <- utils::read.csv(p_cvp,  check.names = FALSE)
OREL <- utils::read.csv(p_orel, check.names = FALSE)
E3CF <- jsonlite::read_json(p_e3cf)
SPLT <- utils::read.csv(p_splt, check.names = FALSE)

# Split counts asserted; no participant CHOICE is read anywhere (this
# script touches no raw file; the split file carries IDs and labels only).
stopifnot(sum(SPLT$split == "build") == 344L,
          sum(SPLT$split == "freeze") == 137L,
          sum(SPLT$split == "test") == 205L)
stopifnot(all(SPLT$held_out_until_E7 == (SPLT$split == "test")))
say("split counts asserted: build 344 / freeze 137 (not read) / test 205 (sealed)\n")

stopifnot(nrow(DICT) == 105L, !anyDuplicated(DICT$GameID))

# The frozen build loss baseline (E3 config) reproduces the dictionary's
# p_noloss_ref_exp on loss-informative rows: the scorer's baseline is tied to
# committed constants, not re-fit here.
b0 <- E3CF$loss_baselines_exp[["(Intercept)"]]
b1 <- E3CF$loss_baselines_exp[["d_ev_noloss"]]
li_rows <- DICT$loss_informative == "TRUE" | DICT$loss_informative == TRUE
pb_chk <- plogis(b0 + b1 * DICT$d_ev_noloss[li_rows])
stopifnot(max(abs(pb_chk - DICT$p_noloss_ref_exp[li_rows])) < 1e-12)
say("frozen loss baseline reproduces dictionary p_noloss_ref_exp (max dev < 1e-12)\n")

# =============================================================================
# 3. Form tables, committed sets, carriers
# =============================================================================
form_col <- function(f) paste0("in_Set", f)
form_items <- function(f) {
  d <- DICT[DICT[[form_col(f)]] == "TRUE" | DICT[[form_col(f)]] == TRUE, ]
  d <- d[order(d$GameID), ]
  d$loss_informative <- d$loss_informative == "TRUE" | d$loss_informative == TRUE
  d$has_loss_b <- d$has_loss_b == "TRUE" | d$has_loss_b == TRUE
  d
}
D1F <- form_items(1L)
stopifnot(nrow(D1F) == 24L, sum(D1F$loss_informative) == 6L)

# fm template for the frozen run_set: probs / NB / li / has_loss_b / d_ev /
# p_base. NB is filled per replicate with SIMULATED counts.
make_fm <- function(D, NB) list(probs = D$GameID, NB = NB,
                                li = D$loss_informative,
                                has_loss_b = D$has_loss_b,
                                d_ev = D$d_ev, p_base = D$p_noloss_ref_exp)

# Committed E3 form-1 sets at the gate budgets (identical across targets).
parse_ids <- function(s) as.integer(strsplit(s, ";", fixed = TRUE)[[1]])
sets_by_K <- lapply(stats::setNames(K_SENS, K_SENS), function(K) {
  a <- REAL[REAL$form == 1L & REAL$K == K & REAL$target == PRIMARY, ]
  b <- REAL[REAL$form == 1L & REAL$K == K & REAL$target == TARGETS[2], ]
  a <- a[order(a$set_id), ]; b <- b[order(b$set_id), ]
  stopifnot(nrow(a) == M_SETS, identical(a$game_ids, b$game_ids))
  cols <- lapply(a$game_ids, function(s) {
    cc <- match(parse_ids(s), D1F$GameID)
    stopifnot(!anyNA(cc), !anyDuplicated(cc))
    sort(cc)
  })
  list(set_id = a$set_id, game_ids = a$game_ids, cols = cols,
       r_sb_real = a$r_sb)
})

# ---- E4 feature aggregation, re-derived and verified against the committed
# set_features.csv (the aggregation is deterministic and re-derivable) --------
sens_features <- function(D, cols_list) {
  z_ev <- D$z_ev; z_sd <- D$z_sd; z_eloss <- D$z_eloss
  li <- D$loss_informative; dnl <- D$d_ev_noloss
  ent <- D$choice_entropy_exp; ila <- D$info_loss_avoid_adj
  purrr::map_dfr(cols_list, function(i) {
    livals <- dnl[i][li[i]]; livals <- livals[is.finite(livals)]
    cc <- c(abs_cor(z_ev[i], z_sd[i]), abs_cor(z_ev[i], z_eloss[i]),
            abs_cor(z_sd[i], z_eloss[i]))
    tibble::tibble(
      F1_target_predictor_variance =
        if (length(livals) >= 2L) stats::var(livals) else NA_real_,
      F2_loss_exposure = sum(li[i]),
      F3_predictor_collinearity =
        if (all(is.na(cc))) NA_real_ else max(cc, na.rm = TRUE),
      F4_expected_information = sum(ila[i]),
      F5_choice_entropy = mean(ent[i]))
  })
}
FEAT1 <- lapply(stats::setNames(as.character(K_SENS), K_SENS), function(Kc) {
  K <- as.integer(Kc)
  ft <- sens_features(D1F, sets_by_K[[Kc]]$cols)
  cm <- SF[SF$track == "real_form1" & SF$target == PRIMARY & SF$K == K, ]
  cm <- cm[order(cm$set_id), ]
  stopifnot(nrow(cm) == M_SETS,
            identical(cm$set_id, sets_by_K[[Kc]]$set_id))
  for (f in FEATURES) {
    stopifnot(identical(is.na(ft[[f]]), is.na(cm[[f]])))
    ok <- !is.na(ft[[f]])
    stopifnot(max(abs(ft[[f]][ok] - cm[[f]][ok])) < 1e-9)
  }
  dplyr::bind_cols(tibble::tibble(set_id = sets_by_K[[Kc]]$set_id,
                                  game_ids = sets_by_K[[Kc]]$game_ids), ft)
})
say("E4 feature aggregation re-derived: matches committed set_features at K = %s\n",
    paste(K_SENS, collapse = ", "))

# ---- ALTERNATIVE per-item F1/F2/F3 definitions ("alternative_P2_erratum" in
# the outputs), aggregated to the set level as means of the per-item analogues
sens_features_alt <- function(D, cols_list) {
  li <- D$loss_informative
  pool_li_mean <- mean(D$d_ev_noloss[li])
  zc <- scale(as.matrix(D[, c("z_ev", "z_sd", "z_eloss")]), TRUE, FALSE)
  f3i <- pmax(abs(zc[, 1] * zc[, 2]), abs(zc[, 1] * zc[, 3]),
              abs(zc[, 2] * zc[, 3]))
  f2i <- D$eloss_a + D$eloss_b
  ent <- D$choice_entropy_exp; ila <- D$info_loss_avoid_adj
  purrr::map_dfr(cols_list, function(i) {
    ii <- i[li[i]]
    tibble::tibble(
      F1_target_predictor_variance = if (length(ii) >= 1L)
        mean((D$d_ev_noloss[ii] - pool_li_mean)^2) else NA_real_,
      F2_loss_exposure = mean(f2i[i]),
      F3_predictor_collinearity = mean(f3i[i]),
      F4_expected_information = sum(ila[i]),
      F5_choice_entropy = mean(ent[i]))
  })
}
FEAT1_ALT <- lapply(stats::setNames(as.character(K_SENS), K_SENS), function(Kc) {
  dplyr::bind_cols(
    tibble::tibble(set_id = sets_by_K[[Kc]]$set_id,
                   game_ids = sets_by_K[[Kc]]$game_ids),
    sens_features_alt(D1F, sets_by_K[[Kc]]$cols))
})

# ---- carriers ---------------------------------------------------------------
# Canonical planted carrier g: per-item expected information for the primary
# target (info_loss_avoid_adj), centered and scaled over the form's
# loss-informative pool; 0 on non-loss-informative items (the primary target
# scores only loss-informative items). Mean 0 over the li pool preserves the
# average trait dispersion and the marginal choice rates.
make_carrier <- function(D) {
  g <- rep(0, nrow(D))
  li <- D$loss_informative
  v <- D$info_loss_avoid_adj[li]
  g[li] <- (v - mean(v)) / stats::sd(v)
  g
}
g_can <- make_carrier(D1F)

# Collinearity scenario carriers: g_c = z_c + c0 on the li pool (unit item-SD
# by construction), c0 solved by bisection so the SET-level Pearson
# corr(F2, sum of g_c over the set's li items) over the 2000 committed K = 12
# sets hits the scenario target. The generator plants THROUGH g_c and the
# analysis F4 column is REPLACED by the set-sum of g_c: the scenario world's
# information feature has the constructed collinearity with loss exposure.
coll_corr <- function(g, Kc = "12") {
  f4p <- vapply(sets_by_K[[Kc]]$cols, function(i)
    sum(g[i][D1F$loss_informative[i]]), numeric(1))
  f2 <- FEAT1[[Kc]]$F2_loss_exposure
  stats::cor(f2, f4p)
}
z_c <- g_can                      # already unit-SD, mean 0 on the li pool
solve_c0 <- function(target) {
  lo <- -60; hi <- 60
  for (it in 1:200) {
    mid <- (lo + hi) / 2
    if (coll_corr(z_c + ifelse(D1F$loss_informative, mid, 0)) < target)
      lo <- mid else hi <- mid
  }
  mid <- (lo + hi) / 2
  g <- z_c + ifelse(D1F$loss_informative, mid, 0)
  list(c0 = mid, g = g, achieved = coll_corr(g))
}
coll_carriers <- lapply(stats::setNames(COLL_TARGETS, COLL_TARGETS), solve_c0)
for (ct in names(coll_carriers))
  stopifnot(abs(coll_carriers[[ct]]$achieved - as.numeric(ct)) < 1e-4)
corr_f2f4_observed <- vapply(as.character(K_SENS), function(Kc)
  stats::cor(FEAT1[[Kc]]$F2_loss_exposure, FEAT1[[Kc]]$F4_expected_information),
  numeric(1))
say("observed set-level corr(F2, F4): %s; scenario c0 solved for %s\n",
    paste(sprintf("K%s %.3f", as.character(K_SENS), corr_f2f4_observed),
          collapse = ", "),
    paste(names(coll_carriers), collapse = ", "))

# =============================================================================
# 4. Rate calibration to the committed build empirical choice rates
# =============================================================================
# Per-problem offset c_j so that the frozen-population average P(choose B)
# equals the committed E2 build empirical rate p_choose_B_exp. Evaluated on a
# large frozen draw_true_params sample; deterministic given seed_cal. The
# planted term preserves the calibration in expectation (mean multiplier 1).
set.seed(seed_cal)
cal_pars <- draw_true_params(N_CAL)                              # FROZEN
stopifnot(abs(mean(cal_pars$b_loss) - MU_L) < 0.01)              # MU_L check
eta_all <- with(cal_pars,
  outer(b0, rep(1, nrow(DICT))) +
  outer(b_ev, DICT$z_ev) + outer(b_risk, DICT$z_sd) +
  outer(b_loss, DICT$z_eloss))
p_emp <- DICT$p_choose_B_exp
cal_offset <- vapply(seq_len(nrow(DICT)), function(j) {
  f <- function(cc) mean(plogis(eta_all[, j] + cc)) - p_emp[j]
  stats::uniroot(f, c(-12, 12), tol = 1e-12)$root
}, numeric(1))
cal_achieved <- vapply(seq_len(nrow(DICT)), function(j)
  mean(plogis(eta_all[, j] + cal_offset[j])), numeric(1))
stopifnot(max(abs(cal_achieved - p_emp)) < 1e-9)
rm(eta_all); invisible(gc())
say("rate calibration: 105 offsets solved; max |achieved - empirical| < 1e-9\n")

# =============================================================================
# 5. Planted generator (frozen model + declared discrimination term)
# =============================================================================
# For subject i, problem j:
#   eta_ij = b0_i + b_ev_i z_ev_j + b_risk_i z_sd_j
#            + (MU_L + (1 + beta g_j) d_i) z_eloss_j + calibrated c_j
# with d_i = b_loss_i - MU_L. beta = 0 and c = 0 reduce EXACTLY to the frozen
# analytic_p_b (gate G_GEN). NB counts are Binomial(r = 20, p_ij): 20
# independent Bernoulli trials, the frozen simulate_choices noise model.
gen_p <- function(pars, D, beta, g, offs) {
  d_i <- pars$b_loss - MU_L
  eta <- outer(pars$b0, rep(1, nrow(D))) +
    outer(pars$b_ev, D$z_ev) + outer(pars$b_risk, D$z_sd) +
    (MU_L + outer(d_i, 1 + beta * g)) * outer(rep(1, nrow(pars)), D$z_eloss)
  plogis(sweep(eta, 2, offs, "+"))
}
gen_NB <- function(pars, D, beta, g, offs, seed) {
  P <- gen_p(pars, D, beta, g, offs)
  set.seed(seed)
  matrix(stats::rbinom(length(P), R_REPS, as.vector(P)), nrow = nrow(P))
}

# ---- gate G_GEN: generator equals the frozen analytic_p_b under the declared
# substitution (b_loss := MU_L + (1 + beta g_j) d_i, plus the declared offset)
say("\n=== G_GEN: planted generator vs frozen analytic_p_b (%d checks) ===\n",
    N_GGEN)
set.seed(seed_ggen)
gg_max <- 0
for (it in seq_len(N_GGEN)) {
  i <- sample.int(50L, 1L); j <- sample.int(24L, 1L)
  beta <- sample(c(0, BETA_GRID), 1L)
  cal <- sample(c(TRUE, FALSE), 1L)
  pars <- cal_pars[i, ]
  offs <- if (cal) cal_offset[match(D1F$GameID, DICT$GameID)] else rep(0, 24L)
  p_mine <- gen_p(pars, D1F, beta, g_can, offs)[1, j]
  pars_sub <- pars
  pars_sub$b_loss <- MU_L + (1 + beta * g_can[j]) * (pars$b_loss - MU_L)
  frow <- D1F[j, c("GameID", "z_ev", "z_sd", "z_eloss")]
  p_frozen <- plogis(qlogis(analytic_p_b(pars_sub, frow)) + offs[j])  # FROZEN
  gg_max <- max(gg_max, abs(p_mine - p_frozen))
}
stopifnot(gg_max < 1e-12)
say("G_GEN PASS: max |p_planted - p_frozen_substituted| = %.3g\n", gg_max)

# ---- gate G_EQ: scoring path for simulated counts equals the frozen
# estimate_half (re-proven here on planted data, not inherited) --------------
say("=== G_EQ: half_targets on simulated counts vs frozen estimate_half (%d halves) ===\n",
    N_GEQ)
set.seed(seed_geq)
offs1 <- cal_offset[match(D1F$GameID, DICT$GameID)]
eq_max <- 0; eq_na <- 0L
feat_sd <- D1F[, c("GameID", "z_ev", "z_sd", "z_eloss", "d_ev",
                   "loss_informative", "has_loss_b")]
feat_sd$p_base_exp <- D1F$p_noloss_ref_exp
for (it in seq_len(N_GEQ)) {
  beta <- sample(c(0, BETA_GRID), 1L)
  pars <- cal_pars[sample.int(200L, 1L), ]
  P <- gen_p(pars, D1F, beta, g_can, offs1)[1, ]
  B_mat <- matrix(0L, 24L, 25L)
  B_mat[, POS] <- matrix(stats::rbinom(24L * R_REPS, 1L, rep(P, R_REPS)),
                         nrow = 24L)
  sdx <- list(feat = feat_sd, B_mat = B_mat)
  kh <- sample(3:6, 1L)
  cols <- sample.int(24L, kh)
  a <- estimate_half(sdx, cols, POS, REGIME)                     # FROZEN
  NBh <- matrix(rowSums(B_mat[cols, POS, drop = FALSE]), nrow = 1L)
  fm1 <- make_fm(D1F, NULL)
  b <- half_targets(NBh, fm1, cols, R_REPS)                      # FROZEN (E3)
  for (nm in TARGETS) {
    x <- unname(a[[nm]]); y <- unname(b[[nm]])
    if (is.na(x) != is.na(y)) eq_na <- eq_na + 1L
    else if (!is.na(x)) eq_max <- max(eq_max, abs(x - y))
  }
}
stopifnot(eq_na == 0L, eq_max < 1e-12)
say("G_EQ PASS: max |d| = %.3g, NA mismatches = %d\n", eq_max, eq_na)

# =============================================================================
# 6. The frozen E4 gate harness (extracted machinery, verbatim call sequence)
# =============================================================================
build_cells <- function(feat_by_K, r_sb_by_K) {
  out <- list()
  for (Kc in names(feat_by_K)) {
    d <- feat_by_K[[Kc]]
    d$r_sb <- r_sb_by_K[[Kc]]
    keep <- is.finite(d$r_sb) &
      Reduce(`&`, lapply(FEATURES, function(f) is.finite(d[[f]])))
    out[[Kc]] <- d[keep, , drop = FALSE]
  }
  out
}

# Mirrors E4's per-cell flow: cv_cell + set_problem_matrix +
# bootstrap_cell_group + the per-K and pooled-across-gate-K statistics + the
# full-sample standardized FULL coefficients. All computation inside is done
# by the extracted frozen functions.
run_gate <- function(cells, pool, seeds_cv, seed_boot) {
  cvs <- list(); grp <- list()
  for (Kc in names(cells)) {
    d <- cells[[Kc]]
    if (nrow(d) < MIN_CELL_N) next
    cv <- cv_cell(d, seeds_cv[[Kc]])                             # FROZEN
    if (is.null(cv)) next
    W <- set_problem_matrix(d, pool)                             # FROZEN
    grp[[Kc]] <- list(y = cv$y, oof = cv$oof, W = W)
    cvs[[Kc]] <- list(cv = cv, d = d)
  }
  if (!length(grp)) return(NULL)
  if (!length(intersect(as.character(K_GATE), names(cvs)))) return(NULL)
  B <- bootstrap_cell_group(grp, pool, seed_boot)                # FROZEN
  per_K <- purrr::map_dfr(names(cvs), function(Kc) {
    cv <- cvs[[Kc]]$cv; d <- cvs[[Kc]]$d
    pooled <- vapply(MODELS, function(m) sp(cv$oof[, m], cv$y), numeric(1))
    rmean  <- colMeans(cv$rep_sp, na.rm = TRUE)[MODELS]
    r2mean <- colMeans(cv$rep_r2, na.rm = TRUE)[MODELS]
    Xs <- scale(as.matrix(d[, FEATURES])); Xs[!is.finite(Xs)] <- 0
    bfull <- safe_lsq(cbind(1, Xs), d$r_sb)                      # FROZEN
    ki <- match(Kc, dimnames(B)[[2]])
    dfb <- B[, ki, "FULL"] - B[, ki, "BASE_category"]
    tibble::tibble(
      K = as.integer(Kc), n_sets_used = cv$n,
      oof_spearman_pooled = pooled["FULL"],
      oof_spearman_repeatmean = rmean["FULL"],
      gate_stat = min(pooled["FULL"], rmean["FULL"]),
      boot_ci_lo = ci_lo(B[, ki, "FULL"]), boot_ci_hi = ci_hi(B[, ki, "FULL"]),
      base_pooled = pooled["BASE_category"],
      base_repeatmean = rmean["BASE_category"],
      delta_vs_base_category = pooled["FULL"] - pooled["BASE_category"],
      delta_repeatmean = rmean["FULL"] - rmean["BASE_category"],
      delta_r2 = r2mean["FULL"] - r2mean["BASE_category"],
      delta_boot_ci_lo = ci_lo(dfb), delta_boot_ci_hi = ci_hi(dfb),
      delta_boot_support = mean(dfb > 0, na.rm = TRUE),
      beta_F1 = bfull[2], beta_F2 = bfull[3], beta_F3 = bfull[4],
      beta_F4 = bfull[5], beta_F5 = bfull[6])
  })
  ks <- intersect(as.character(K_GATE), names(cvs))
  pooled_obs <- mean(per_K$oof_spearman_pooled[match(ks, as.character(per_K$K))])
  pooled_rm  <- mean(per_K$oof_spearman_repeatmean[match(ks, as.character(per_K$K))])
  draws_full <- rowMeans(B[, ks, "FULL", drop = FALSE])
  draws_base <- rowMeans(B[, ks, "BASE_category", drop = FALSE])
  pooled <- tibble::tibble(
    K_pooled = paste(ks, collapse = ";"),
    pooled_gate_stat = min(pooled_obs, pooled_rm),
    pooled_oof_spearman = pooled_obs, pooled_oof_repeatmean = pooled_rm,
    pooled_ci_lo = ci_lo(draws_full), pooled_ci_hi = ci_hi(draws_full),
    pooled_delta = pooled_obs -
      mean(per_K$base_pooled[match(ks, as.character(per_K$K))]),
    pooled_delta_support = mean(draws_full - draws_base > 0, na.rm = TRUE))
  list(per_K = per_K, pooled = pooled)
}

empty_pooled <- function() tibble::tibble(
  K_pooled = NA_character_, pooled_gate_stat = NA_real_,
  pooled_oof_spearman = NA_real_, pooled_oof_repeatmean = NA_real_,
  pooled_ci_lo = NA_real_, pooled_ci_hi = NA_real_,
  pooled_delta = NA_real_, pooled_delta_support = NA_real_)

gate_clauses <- function(gres) {
  if (is.null(gres)) return(tibble::tibble(
    cells_complete = FALSE, c1 = FALSE, c2 = FALSE, c3 = FALSE,
    pass_gate = FALSE, pass_threshold = FALSE, sign_correct = FALSE))
  pk <- gres$per_K[gres$per_K$K %in% K_GATE, ]
  complete <- nrow(pk) == length(K_GATE)
  c1_perK <- complete && all(pk$gate_stat >= GATE_SPEARMAN_MIN) &&
    all(pk$boot_ci_lo > 0)
  c1_pool <- isTRUE(gres$pooled$pooled_gate_stat >= GATE_SPEARMAN_MIN) &&
    isTRUE(gres$pooled$pooled_ci_lo > 0)
  c1 <- c1_perK && c1_pool
  c2 <- complete && all((pk$delta_vs_base_category > 0 | pk$delta_r2 > 0) &
                          pk$delta_boot_support >= GATE_BOOT_SUPPORT)
  c3 <- complete && all(abs(pk$beta_F4) >= CARRIER_MIN_BETA & pk$beta_F4 > 0)
  thr <- complete && all(pk$gate_stat >= GATE_SPEARMAN_MIN) &&
    isTRUE(gres$pooled$pooled_gate_stat >= GATE_SPEARMAN_MIN)
  sgn <- complete && all(pk$beta_F4 > 0)
  tibble::tibble(cells_complete = complete, c1 = c1, c2 = c2, c3 = c3,
                 pass_gate = c1 && c2 && c3, pass_threshold = thr,
                 sign_correct = sgn)
}

# ---- gate G_E4EQ: the harness reproduces the committed frozen-gate numbers --
# Rebuild the real_form1 / loss_avoid_adj cells from the committed
# set_features.csv rows (file order), run the harness with E4's own seed
# derivations (track_index 1, target_index 1), and assert every reproduced
# statistic equals the committed cv_performance / pooled values.
say("\n=== G_E4EQ: extracted gate machinery vs committed E4 numbers ===\n")
e4_cells <- list(); pool_form1 <- sort(unique(unlist(
  lapply(SF$game_ids[SF$track == "real_form1"], parse_ids))))
stopifnot(identical(pool_form1, D1F$GameID))   # planted arms use the same pool
for (K in K_GRID) {
  d <- SF[SF$track == "real_form1" & SF$target == PRIMARY & SF$K == K, ]
  keep <- is.finite(d$r_sb) &
    Reduce(`&`, lapply(FEATURES, function(f) is.finite(d[[f]])))
  e4_cells[[as.character(K)]] <- d[keep, , drop = FALSE]
}
e4_seeds <- lapply(stats::setNames(as.character(K_GRID), K_GRID), function(Kc)
  SEED_E4 + 1000L * 1L + 100L * 1L + as.integer(Kc))
g4 <- run_gate(e4_cells, pool_form1, e4_seeds,
               SEED_E4 + 500000L + 1000L * 1L + 100L * 1L)
e4eq_rows <- purrr::map_dfr(K_GRID, function(K) {
  mine <- g4$per_K[g4$per_K$K == K, ]
  ref <- CVE4[CVE4$track == "real_form1" & CVE4$target == PRIMARY &
                CVE4$K == K & CVE4$model == "FULL", ]
  tibble::tibble(
    K = K, quantity = c("n_sets_used", "oof_spearman_pooled",
                        "oof_spearman_repeatmean", "boot_ci_lo", "boot_ci_hi",
                        "delta_vs_base_category", "delta_boot_support",
                        "delta_r2", "delta_repeatmean"),
    committed = c(ref$n_sets_used, ref$oof_spearman_pooled,
                  ref$oof_spearman_repeatmean, ref$boot_ci_lo, ref$boot_ci_hi,
                  ref$delta_vs_base_category, ref$delta_boot_support,
                  ref$delta_r2_vs_base_category,
                  ref$delta_repeatmean_vs_base_category),
    reproduced = c(mine$n_sets_used, mine$oof_spearman_pooled,
                   mine$oof_spearman_repeatmean, mine$boot_ci_lo,
                   mine$boot_ci_hi, mine$delta_vs_base_category,
                   mine$delta_boot_support, mine$delta_r2,
                   mine$delta_repeatmean))
})
refp <- CVP4[CVP4$track == "real_form1" & CVP4$target == PRIMARY &
               CVP4$model == "FULL", ]
e4eq_rows <- dplyr::bind_rows(e4eq_rows, tibble::tibble(
  K = NA_integer_,
  quantity = c("pooled_oof_spearman", "pooled_oof_repeatmean",
               "pooled_ci_lo", "pooled_ci_hi", "pooled_delta_support"),
  committed = c(refp$pooled_oof_spearman, refp$pooled_oof_spearman_repeatmean,
                refp$pooled_boot_ci_lo, refp$pooled_boot_ci_hi,
                refp$pooled_delta_support),
  reproduced = c(g4$pooled$pooled_oof_spearman, g4$pooled$pooled_oof_repeatmean,
                 g4$pooled$pooled_ci_lo, g4$pooled$pooled_ci_hi,
                 g4$pooled$pooled_delta_support)))
e4eq_rows$abs_diff <- abs(e4eq_rows$committed - e4eq_rows$reproduced)
stopifnot(max(e4eq_rows$abs_diff) < E4EQ_TOL)
say("G_E4EQ PASS: %d statistics reproduced, max |diff| = %.3g\n",
    nrow(e4eq_rows), max(e4eq_rows$abs_diff))
save_csv(e4eq_rows, "sens_e4eq_check")

# Committed real-data reference values for the attenuation exhibit.
real_ref <- g4$per_K |> dplyr::mutate(track = "real_form1_committed")

# =============================================================================
# 7. Replicate engine
# =============================================================================
one_replicate <- function(arm_row, rep_i, want_r_sb = FALSE) {
  arm <- arm_row$arm
  g <- switch(arm_row$carrier,
              canonical = g_can,
              orthogonalized_random = {
                set.seed(seed_negcar(arm, rep_i))
                li <- D1F$loss_informative
                z <- stats::rnorm(sum(li))
                r <- stats::residuals(stats::lm(z ~ g_can[li]))
                h <- rep(0, nrow(D1F)); h[li] <- r / stats::sd(r)
                h
              },
              coll_carriers[[as.character(arm_row$coll_target)]]$g)
  offs <- if (arm_row$calibrated) offs1 else rep(0, 24L)
  n_subj <- if (arm_row$family == "high_precision") N_HP_SUBJ else 63L
  set.seed(seed_par(arm, rep_i))
  pars <- draw_true_params(n_subj)                               # FROZEN
  NB <- gen_NB(pars, D1F, arm_row$beta, g, offs, seed_choice(arm, rep_i))
  fm <- make_fm(D1F, NB)
  n_split <- if (arm_row$family == "high_precision") N_SPLIT_HP else N_SPLIT_SENS
  Ks <- if (arm_row$family == "high_precision") as.character(c(8L, 10L, 12L)) else
    as.character(K_SENS)
  r_sb_by_K <- lapply(stats::setNames(Ks, Ks), function(Kc) {
    K <- as.integer(Kc)
    sb <- sets_by_K[[Kc]]
    n_use <- if (arm_row$family == "high_precision" && K %in% c(8L, 10L))
      M_SETS_HP_812 else length(sb$cols)
    vapply(seq_len(n_use), function(i)
      run_set(fm, sb$cols[[i]], n_split,                          # FROZEN (E3)
              seed_split(arm, rep_i, K, sb$set_id[i]))[[PRIMARY]]$r_sb,
      numeric(1))
  })
  if (arm_row$family == "high_precision") {
    ach <- colMeans(gen_p(pars, D1F, arm_row$beta, g, offs))
    return(list(arm = arm, rep = rep_i, r_sb_by_K = r_sb_by_K,
                achieved_rate = ach))
  }
  feat_use <- FEAT1
  if (arm_row$family == "collinearity") {
    for (Kc in names(feat_use)) {
      f4p <- vapply(sets_by_K[[Kc]]$cols, function(i)
        sum(g[i][D1F$loss_informative[i]]), numeric(1))
      feat_use[[Kc]]$F4_expected_information <- f4p
    }
  }
  cells <- build_cells(feat_use,
                       lapply(r_sb_by_K, function(v)
                         c(v, rep(NA_real_, M_SETS - length(v)))[1:M_SETS]))
  seeds_cvv <- lapply(stats::setNames(as.character(K_SENS), K_SENS),
                      function(Kc) seed_cv(arm, rep_i, as.integer(Kc), 0L))
  gres <- run_gate(cells, D1F$GameID, seeds_cvv, seed_boot_g(arm, rep_i, 0L))
  out <- gate_row(gres, arm_row, rep_i, "canonical_features")
  if (want_r_sb) list(row = out, r_sb_by_K = r_sb_by_K) else list(row = out)
}

gate_row <- function(gres, arm_row, rep_i, variant) {
  head_tbl <- tibble::tibble(
    arm = arm_row$arm, family = arm_row$family, beta = arm_row$beta,
    carrier = arm_row$carrier, coll_target = arm_row$coll_target,
    calibrated = arm_row$calibrated, rep = rep_i, variant = variant)
  cl <- gate_clauses(gres)
  if (is.null(gres)) return(dplyr::bind_cols(head_tbl, cl, empty_pooled()))
  dplyr::bind_cols(
    head_tbl, cl, gres$pooled,
    gres$per_K |>
      dplyr::select(K, n_sets_used, gate_stat, boot_ci_lo,
                    delta_vs_base_category, delta_r2, delta_boot_support,
                    beta_F4) |>
      tidyr::pivot_wider(names_from = K,
                         values_from = c(n_sets_used, gate_stat, boot_ci_lo,
                                         delta_vs_base_category, delta_r2,
                                         delta_boot_support, beta_F4)))
}

regate_alt <- function(arm_row, rep_i, r_sb_by_K) {
  cells <- build_cells(FEAT1_ALT, r_sb_by_K)
  seeds_cvv <- lapply(stats::setNames(as.character(K_SENS), K_SENS),
                      function(Kc) seed_cv(arm_row$arm, rep_i, as.integer(Kc), 1L))
  gres <- run_gate(cells, D1F$GameID, seeds_cvv,
                   seed_boot_g(arm_row$arm, rep_i, 1L))
  gate_row(gres, arm_row, rep_i, "alternative_features")
}

# =============================================================================
# 8. Production runs (parallel; seeds resolved in the job list)
# =============================================================================
cl <- parallel::makeCluster(min(N_WORKERS, max(1L, parallel::detectCores() - 2L)))
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(tidyverse); library(here) })
  source(here("R", "step2_reliability.R"))
  source(here("R", "simulate_recovery.R"))
  NULL
}))
parallel::clusterExport(cl, c(
  "REGIME", "POS", "R_REPS", "TARGETS", "MIN_PAIRS", "ev_ranks", "half_targets",
  "run_set", "SEED_E4", "N_REP", "TRAIN_FRAC", "N_BOOT", "RIDGE_GRID",
  "PLS_NCOMP", "INNER_FOLDS", "K_GRID", "K_GATE", "FEATURES", "MODELS",
  "GATE_SPEARMAN_MIN", "GATE_BOOT_SUPPORT", "MIN_CELL_N", "CARRIER_MIN_BETA",
  "abs_cor", "ridge_fit", "simpls_fit", "inner_folds", "safe_lsq",
  "fit_predict", "sp", "cv_cell", "set_problem_matrix", "bootstrap_cell_group",
  "ci_lo", "ci_hi", "PRIMARY", "MU_L", "D1F", "sets_by_K", "FEAT1", "FEAT1_ALT",
  "make_fm", "gen_p", "gen_NB", "g_can", "coll_carriers", "offs1", "cal_offset",
  "build_cells", "run_gate", "gate_clauses", "empty_pooled", "gate_row",
  "one_replicate", "regate_alt",
  "seed_par", "seed_choice", "seed_negcar", "seed_split", "seed_cv",
  "seed_boot_g", "N_SPLIT_SENS", "N_SPLIT_HP", "N_HP_SUBJ", "M_SETS",
  "M_SETS_HP_812", "K_SENS", "parse_ids", "draw_true_params", "analytic_p_b"),
  envir = environment())

sim_arms <- arm_tbl[arm_tbl$family != "high_precision", ]
jobs <- purrr::map_dfr(seq_len(nrow(sim_arms)), function(i)
  tibble::tibble(idx = i, rep = seq_len(sim_arms$n_rep[i]),
                 keep_rsb = sim_arms$family[i] == "main"))
parallel::clusterExport(cl, c("jobs", "sim_arms"), envir = environment())
say("\n=== production: %d simulation replicates over %d arms (workers: %d) ===\n",
    nrow(jobs), nrow(sim_arms), length(cl))
t0 <- Sys.time()
res <- parallel::parLapply(cl, seq_len(nrow(jobs)), function(ji) {
  jb <- jobs[ji, ]
  ar <- sim_arms[jb$idx, ]
  r <- one_replicate(ar, jb$rep, want_r_sb = jb$keep_rsb)
  alt <- if (jb$keep_rsb) regate_alt(ar, jb$rep, r$r_sb_by_K) else NULL
  list(row = r$row, alt = alt)
})
say("production done in %.1f min\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
detail <- dplyr::bind_rows(lapply(res, `[[`, "row"))
detail_alt <- dplyr::bind_rows(lapply(res, `[[`, "alt"))
rm(res); invisible(gc())

# ---- high-precision arm -----------------------------------------------------
hp_arms <- arm_tbl[arm_tbl$family == "high_precision", ]
hp_jobs <- purrr::map_dfr(seq_len(nrow(hp_arms)), function(i)
  tibble::tibble(idx = i, rep = seq_len(hp_arms$n_rep[i])))
parallel::clusterExport(cl, c("hp_jobs", "hp_arms"), envir = environment())
say("\n=== high-precision arm: %d runs (N = %d simulated subjects) ===\n",
    nrow(hp_jobs), N_HP_SUBJ)
t0 <- Sys.time()
hp_res <- parallel::parLapply(cl, seq_len(nrow(hp_jobs)), function(ji) {
  jb <- hp_jobs[ji, ]
  one_replicate(hp_arms[jb$idx, ], jb$rep)
})
say("high-precision done in %.1f min\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

# ---- seven-form generalization (hp rep 1 machinery, K* per form) ------------
# Forms 2-7 at the largest committed within-form budget <= 12 (form pools cap
# K); first M_SETS_FORM committed sets; N = 2000 simulated subjects with
# full-exposure within form (forms 6-7 violate the common-set premise on real
# data and are flagged). Carrier per form: centered/scaled item information
# over that form's loss-informative pool. SIMULATION.
form_K <- vapply(2:7, function(f)
  max(REAL$K[REAL$form == f & REAL$target == PRIMARY]), integer(1))
seven_jobs <- purrr::map_dfr(seq_along(BETA_GRID), function(bi)
  tibble::tibble(form = 2:7, K = form_K, beta = BETA_GRID[bi], beta_idx = bi))
say("\n=== seven-form generalization: forms 2-7 at K* = %s, %d runs ===\n",
    paste(form_K, collapse = "/"), nrow(seven_jobs))
parallel::clusterExport(cl, c("DICT", "REAL", "form_items", "form_col",
                              "make_carrier", "sens_features", "M_SETS_FORM",
                              "BETA_GRID", "QUANT_GRAD", "seven_jobs"),
                        envir = environment())
t0 <- Sys.time()
seven_res <- parallel::parLapply(cl, seq_len(nrow(seven_jobs)), function(ji) {
  jb <- seven_jobs[ji, ]
  Df <- form_items(jb$form)
  gf <- make_carrier(Df)
  offf <- cal_offset[match(Df$GameID, DICT$GameID)]
  a <- REAL[REAL$form == jb$form & REAL$K == jb$K & REAL$target == PRIMARY, ]
  a <- a[order(a$set_id), ][seq_len(min(M_SETS_FORM, sum(REAL$form == jb$form &
        REAL$K == jb$K & REAL$target == PRIMARY))), ]
  colsf <- lapply(a$game_ids, function(s) sort(match(parse_ids(s), Df$GameID)))
  arm_id <- 60L + jb$form
  set.seed(seed_par(arm_id, jb$beta_idx))
  pars <- draw_true_params(N_HP_SUBJ)                            # FROZEN
  NB <- gen_NB(pars, Df, jb$beta, gf, offf, seed_choice(arm_id, jb$beta_idx))
  fmf <- make_fm(Df, NB)
  r_sb <- vapply(seq_along(colsf), function(i)
    run_set(fmf, colsf[[i]], N_SPLIT_HP,                          # FROZEN (E3)
            seed_split(arm_id, jb$beta_idx, 8L, a$set_id[i]))[[PRIMARY]]$r_sb,
    numeric(1))
  f4 <- sens_features(Df, colsf)$F4_expected_information
  ok <- is.finite(r_sb)
  qs <- stats::quantile(f4[ok], c(QUANT_GRAD, 1 - QUANT_GRAD))
  tibble::tibble(
    form = jb$form, K_star = jb$K, beta = jb$beta,
    n_sets = length(r_sb), n_finite = sum(ok),
    n_loss_informative_pool = sum(Df$loss_informative),
    mean_r_sb = mean(r_sb[ok]),
    grad_f4_quintile = mean(r_sb[ok][f4[ok] >= qs[2]]) -
      mean(r_sb[ok][f4[ok] <= qs[1]]),
    spearman_f4_r_sb = suppressWarnings(
      stats::cor(f4[ok], r_sb[ok], method = "spearman")),
    premise_note = ifelse(jb$form %in% 6:7,
                          "common-set premise FALSE on real data; simulated full exposure",
                          "common-set premise holds"))
})
seven_form <- dplyr::bind_rows(seven_res) |> dplyr::arrange(form, beta)
say("seven-form done in %.1f min\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

parallel::stopCluster(cl)

# =============================================================================
# 9. High-precision aggregation: r_true, induced change, calibration drift
# =============================================================================
hp_beta <- vapply(hp_res, function(r) hp_arms$beta[hp_arms$arm == r$arm],
                  numeric(1))
hp_setlevel <- purrr::map_dfr(seq_along(hp_res), function(i) {
  r <- hp_res[[i]]
  purrr::map_dfr(names(r$r_sb_by_K), function(Kc) {
    v <- r$r_sb_by_K[[Kc]]
    tibble::tibble(data_kind = "SIMULATION", beta = hp_beta[i], hp_rep = r$rep,
                   K = as.integer(Kc), set_id = sets_by_K[[Kc]]$set_id[seq_along(v)],
                   r_sb_hp = v)
  })
})
hp_mean <- hp_setlevel |>
  dplyr::group_by(beta, K, set_id) |>
  dplyr::summarise(r_sb_hp = mean(r_sb_hp, na.rm = TRUE), .groups = "drop")

r_true_tbl <- purrr::map_dfr(BETA_GRID, function(b) {
  purrr::map_dfr(c(8L, 10L, 12L), function(K) {
    h <- hp_mean[hp_mean$beta == b & hp_mean$K == K, ]
    ft <- FEAT1[[as.character(K)]]
    f4 <- ft$F4_expected_information[match(h$set_id, ft$set_id)]
    ok <- is.finite(h$r_sb_hp)
    per_rep <- vapply(1:N_HP_REPS, function(rp) {
      hh <- hp_setlevel[hp_setlevel$beta == b & hp_setlevel$K == K &
                          hp_setlevel$hp_rep == rp, ]
      f4r <- ft$F4_expected_information[match(hh$set_id, ft$set_id)]
      okr <- is.finite(hh$r_sb_hp)
      suppressWarnings(stats::cor(f4r[okr], hh$r_sb_hp[okr], method = "spearman"))
    }, numeric(1))
    tibble::tibble(
      beta = b, K = K, n_sets = nrow(h), n_finite = sum(ok),
      r_true_spearman = suppressWarnings(
        stats::cor(f4[ok], h$r_sb_hp[ok], method = "spearman")),
      r_true_pearson = suppressWarnings(stats::cor(f4[ok], h$r_sb_hp[ok])),
      r_true_rep_range = max(per_rep) - min(per_rep))
  })
})

grad_tbl <- purrr::map_dfr(BETA_GRID, function(b) {
  h <- hp_mean[hp_mean$beta == b & hp_mean$K == 12L, ]
  ft <- FEAT1[["12"]]
  f4 <- ft$F4_expected_information[match(h$set_id, ft$set_id)]
  ok <- is.finite(h$r_sb_hp)
  qs <- stats::quantile(f4[ok], c(QUANT_GRAD, 1 - QUANT_GRAD))
  tibble::tibble(beta = b,
                 grad_f4_quintile_K12 = mean(h$r_sb_hp[ok][f4[ok] >= qs[2]]) -
                   mean(h$r_sb_hp[ok][f4[ok] <= qs[1]]),
                 mean_r_sb_K12 = mean(h$r_sb_hp[ok]))
})
grad0 <- grad_tbl$grad_f4_quintile_K12[grad_tbl$beta == 0]
mean0 <- grad_tbl$mean_r_sb_K12[grad_tbl$beta == 0]
grad_tbl <- grad_tbl |>
  dplyr::mutate(induced_reliability_change_K12 = grad_f4_quintile_K12 - grad0,
                mean_r_sb_shift_vs_beta0 = mean_r_sb_K12 - mean0)

cal_drift <- purrr::map_dfr(seq_along(hp_res), function(i) {
  r <- hp_res[[i]]
  tibble::tibble(beta = hp_beta[i], hp_rep = r$rep,
                 max_abs_rate_drift = max(abs(r$achieved_rate -
                   p_emp[match(D1F$GameID, DICT$GameID)])),
                 mean_abs_rate_drift = mean(abs(r$achieved_rate -
                   p_emp[match(D1F$GameID, DICT$GameID)])))
})

# =============================================================================
# 10. Detection curves and pre-declared gate evaluation
# =============================================================================
bin_ci <- function(x, n) {
  lo <- if (x == 0) 0 else stats::qbeta(0.025, x, n - x + 1)
  hi <- if (x == n) 1 else stats::qbeta(0.975, x + 1, n - x)
  c(lo, hi)
}
curve_from <- function(d) {
  d |>
    dplyr::group_by(beta) |>
    dplyr::summarise(
      n_sim = dplyr::n(),
      p_gate = mean(pass_gate), p_threshold = mean(pass_threshold),
      p_sign = mean(sign_correct),
      p_c1 = mean(c1), p_c2 = mean(c2), p_c3 = mean(c3),
      mc_se_p_gate = sqrt(p_gate * (1 - p_gate) / n_sim),
      p_gate_ci_lo = bin_ci(sum(pass_gate), n_sim)[1],
      p_gate_ci_hi = bin_ci(sum(pass_gate), n_sim)[2],
      mean_pooled_gate_stat = mean(pooled_gate_stat, na.rm = TRUE),
      sd_pooled_gate_stat = stats::sd(pooled_gate_stat, na.rm = TRUE),
      mean_gate_stat_K12 = mean(gate_stat_12, na.rm = TRUE),
      mean_delta_K12 = mean(delta_vs_base_category_12, na.rm = TRUE),
      mean_support_K12 = mean(delta_boot_support_12, na.rm = TRUE),
      mean_beta_F4_K12 = mean(beta_F4_12, na.rm = TRUE),
      .groups = "drop")
}

main_curve <- curve_from(detail |> dplyr::filter(family == "main"))
sensitivity_detection_by_beta <- main_curve |>
  dplyr::left_join(r_true_tbl |> dplyr::filter(K == 12L) |>
                     dplyr::select(beta, r_true_spearman_K12 = r_true_spearman,
                                   r_true_pearson_K12 = r_true_pearson,
                                   r_true_rep_range), by = "beta") |>
  dplyr::left_join(grad_tbl, by = "beta") |>
  dplyr::mutate(data_kind = "SIMULATION", .before = 1)
save_csv(sensitivity_detection_by_beta, "sensitivity_detection_by_beta")

neg_curve <- curve_from(detail |> dplyr::filter(family == "negative_control"))
rate_curve <- curve_from(detail |> dplyr::filter(family == "rate_sensitivity"))
alt_curve <- curve_from(detail_alt)

coll_curve <- detail |>
  dplyr::filter(family == "collinearity") |>
  dplyr::group_by(coll_target, beta) |>
  dplyr::summarise(
    n_sim = dplyr::n(),
    p_gate = mean(pass_gate), p_threshold = mean(pass_threshold),
    p_sign = mean(sign_correct), p_c2_beyond = mean(c2),
    mc_se_p_gate = sqrt(p_gate * (1 - p_gate) / n_sim),
    mean_pooled_gate_stat = mean(pooled_gate_stat, na.rm = TRUE),
    mean_delta_K12 = mean(delta_vs_base_category_12, na.rm = TRUE),
    mean_support_K12 = mean(delta_boot_support_12, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(
    data_kind = "SIMULATION", .before = 1) |>
  dplyr::mutate(
    achieved_corr_K12 = vapply(coll_target, function(ct)
      coll_carriers[[as.character(ct)]]$achieved, numeric(1)),
    carrier_c0 = vapply(coll_target, function(ct)
      coll_carriers[[as.character(ct)]]$c0, numeric(1)),
    observed_corr_K12_reference = corr_f2f4_observed[["12"]])
save_csv(coll_curve, "collinearity_scenarios")

# ---- pre-declared gates -----------------------------------------------------
matters_idx <- which(grad_tbl$induced_reliability_change_K12 >= MATTERS_SB)
beta_matters <- if (length(matters_idx)) min(BETA_GRID[matters_idx]) else NA_real_
floor_idx <- which(main_curve$p_gate >= POWER_BAR & main_curve$beta > 0)
beta_floor <- if (length(floor_idx)) min(main_curve$beta[floor_idx]) else NA_real_
floor_idx_alt <- which(alt_curve$p_gate >= POWER_BAR & alt_curve$beta > 0)
beta_floor_alt <- if (length(floor_idx_alt))
  min(alt_curve$beta[floor_idx_alt]) else NA_real_

p_gate_at <- function(curve, b) {
  v <- curve$p_gate[curve$beta == b]
  if (length(v)) v else NA_real_
}
p_pos <- p_gate_at(main_curve, BETA_POS)
p_neg <- neg_curve$p_gate[1]
p_base0 <- p_gate_at(main_curve, 0)

g_controls_pos <- isTRUE(p_pos >= POWER_BAR)
g_controls_neg <- isTRUE(p_neg <= NEG_MAX_P)
g_controls_pass <- g_controls_pos && g_controls_neg

g_power_main <- !is.na(beta_matters) &&
  isTRUE(p_gate_at(main_curve, beta_matters) >= POWER_BAR)
g_power_alt <- !is.na(beta_matters) &&
  isTRUE(p_gate_at(alt_curve, beta_matters) >= POWER_BAR)
g_power_pass <- g_power_main && g_power_alt

detect_floor_verdict <- if (!g_controls_pass) {
  "VOID (G-controls failed; detection probabilities cannot be interpreted as mechanism detection)"
} else if (is.na(beta_floor)) {
  "DETECTABILITY-limited (no grid effect reaches power 0.80; minimum detectable effect exceeds the grid)"
} else if (!is.na(beta_matters) && beta_floor <= beta_matters) {
  "INFORMATIVE (minimum detectable effect at or below the design-relevant effect)"
} else {
  "DETECTABILITY-limited (minimum detectable effect exceeds the design-relevant effect)"
}

# =============================================================================
# 11. Attenuation bound (committed E3b outcome noise; no generator involved)
# =============================================================================
orl <- OREL[OREL$target == PRIMARY, ]
att_K <- purrr::map_dfr(c(6L, 8L, 10L, 12L), function(K) {
  o <- orl[orl$K == K, ]
  rr <- real_ref[real_ref$K == K, ]
  tibble::tibble(
    component = "per_budget", K = K,
    rho_o_out_of_subject = o$outcome_reliability,
    rho_o_ci_lo = o$ci_lo, rho_o_ci_hi = o$ci_hi,
    rho_o_sb_projected_full63 = o$outcome_rel_sb_projected,
    max_observable_spearman = sqrt(o$outcome_reliability),
    max_observable_sb_projected = sqrt(o$outcome_rel_sb_projected),
    committed_real_gate_stat = if (nrow(rr)) rr$gate_stat else NA_real_,
    committed_real_oof_pooled = if (nrow(rr)) rr$oof_spearman_pooled else NA_real_,
    real_exceeds_out_of_subject_ceiling =
      if (nrow(rr)) rr$gate_stat > sqrt(o$outcome_reliability) else NA,
    note = "committed real_form1 association shares subjects and overlapping items across sets, so it can exceed the out-of-subject ceiling; that excess is shared-noise fit, not transferable mechanism signal")
})
o12 <- orl[orl$K == 12L, ]
att_req <- tibble::tibble(
  component = "required_r_true_K12", K = 12L,
  target_oof = c(0.20, 0.30, 0.20, 0.30),
  rho_basis = c("out_of_subject", "out_of_subject", "sb_projected_full63",
                "sb_projected_full63"),
  r_true_required = c(0.20 / sqrt(o12$outcome_reliability),
                      0.30 / sqrt(o12$outcome_reliability),
                      0.20 / sqrt(o12$outcome_rel_sb_projected),
                      0.30 / sqrt(o12$outcome_rel_sb_projected)))
att_beta <- r_true_tbl |>
  dplyr::filter(K == 12L) |>
  dplyr::transmute(component = "planted_r_true_vs_realized", K = K, beta,
                   r_true_spearman = r_true_spearman,
                   predicted_ceiling_out_of_subject =
                     r_true_spearman * sqrt(o12$outcome_reliability),
                   predicted_ceiling_sb_projected =
                     r_true_spearman * sqrt(o12$outcome_rel_sb_projected)) |>
  dplyr::left_join(main_curve |>
                     dplyr::select(beta, realized_mean_gate_stat_K12 =
                                     mean_gate_stat_K12), by = "beta")
att_ctx <- tibble::tibble(
  component = "context", K = NA_integer_,
  note = c(
    "honest real-data set-level association: the committed E4 real_form1 gate statistics above (item-overlap-inflated; see per_budget rows)",
    "P2 real per-ITEM mechanism association, different granularity: grouped-CV Spearman -0.2975, inside its own permutation null"))
attenuation_bound <- dplyr::bind_rows(att_K, att_req, att_beta, att_ctx)
save_csv(attenuation_bound, "attenuation_bound")

# =============================================================================
# 12. Feature-definition robustness (canonical vs alternative definitions)
# =============================================================================
# Real-data arm: the frozen gate machinery on the COMMITTED real form-1 r_sb
# with the alternative features, PAIRED with the committed canonical numbers by
# reusing E4's own seed derivations (identical fold sequences; only the feature
# columns differ). Sim-track arm: same pairing on the committed sim_free105
# r_sb. Simulated arm: the alt-features re-gate of every main-curve replicate.
alt_cells_real <- lapply(stats::setNames(as.character(K_GRID), K_GRID),
                         function(Kc) {
  K <- as.integer(Kc)
  cm <- SF[SF$track == "real_form1" & SF$target == PRIMARY & SF$K == K, ]
  cm <- cm[order(cm$set_id), ]
  if (K %in% K_SENS) {
    fa <- FEAT1_ALT[[Kc]]
    stopifnot(identical(cm$set_id, fa$set_id))
  } else {
    sa <- REAL[REAL$form == 1L & REAL$K == K & REAL$target == PRIMARY, ]
    sa <- sa[order(sa$set_id), ]
    colsK <- lapply(sa$game_ids, function(s) sort(match(parse_ids(s), D1F$GameID)))
    fa <- dplyr::bind_cols(tibble::tibble(set_id = sa$set_id,
                                          game_ids = sa$game_ids),
                           sens_features_alt(D1F, colsK))
  }
  d <- fa
  d$r_sb <- cm$r_sb[match(d$set_id, cm$set_id)]
  keep <- is.finite(d$r_sb) &
    Reduce(`&`, lapply(FEATURES, function(f) is.finite(d[[f]])))
  d[keep, , drop = FALSE]
})
g_alt_real <- run_gate(alt_cells_real, D1F$GameID, e4_seeds,
                       SEED_E4 + 500000L + 1000L * 1L + 100L * 1L)

D105 <- DICT[order(DICT$GameID), ]
D105$loss_informative <- D105$loss_informative == "TRUE" | D105$loss_informative == TRUE
tracks_sim <- sort(unique(SF$track))
ti_sim <- match("sim_free105", tracks_sim)
sim_cells_can <- list(); sim_cells_alt <- list()
for (K in K_GRID) {
  cm <- SF[SF$track == "sim_free105" & SF$target == PRIMARY & SF$K == K, ]
  cm <- cm[order(cm$set_id), ]
  colsK <- lapply(cm$game_ids, function(s) sort(match(parse_ids(s), D105$GameID)))
  fa <- dplyr::bind_cols(tibble::tibble(set_id = cm$set_id,
                                        game_ids = cm$game_ids),
                         sens_features_alt(D105, colsK))
  fa$r_sb <- cm$r_sb
  dcan <- cm
  keepc <- is.finite(dcan$r_sb) &
    Reduce(`&`, lapply(FEATURES, function(f) is.finite(dcan[[f]])))
  keepa <- is.finite(fa$r_sb) &
    Reduce(`&`, lapply(FEATURES, function(f) is.finite(fa[[f]])))
  sim_cells_can[[as.character(K)]] <- dcan[keepc, , drop = FALSE]
  sim_cells_alt[[as.character(K)]] <- fa[keepa, , drop = FALSE]
}
sim_seeds <- lapply(stats::setNames(as.character(K_GRID), K_GRID), function(Kc)
  SEED_E4 + 1000L * ti_sim + 100L * 1L + as.integer(Kc))
sim_pool <- sort(unique(unlist(lapply(SF$game_ids[SF$track == "sim_free105"],
                                      parse_ids))))
g_sim_can <- run_gate(sim_cells_can, sim_pool, sim_seeds,
                      SEED_E4 + 500000L + 1000L * ti_sim + 100L * 1L)
g_sim_alt <- run_gate(sim_cells_alt, sim_pool, sim_seeds,
                      SEED_E4 + 500000L + 1000L * ti_sim + 100L * 1L)
# Fidelity: the canonical sim-track rerun must match committed E4 numbers too.
for (K in K_GRID) {
  mine <- g_sim_can$per_K[g_sim_can$per_K$K == K, ]
  ref <- CVE4[CVE4$track == "sim_free105" & CVE4$target == PRIMARY &
                CVE4$K == K & CVE4$model == "FULL", ]
  stopifnot(abs(mine$oof_spearman_pooled - ref$oof_spearman_pooled) < E4EQ_TOL,
            abs(mine$delta_boot_support - ref$delta_boot_support) < E4EQ_TOL)
}
say("sim_free105 canonical rerun matches committed E4 (fidelity extended)\n")

fd_real <- purrr::map_dfr(list(
  list(def = "canonical_committed", g = g4),
  list(def = "alternative_P2_erratum", g = g_alt_real)), function(x) {
    x$g$per_K |>
      dplyr::transmute(arm = "real_form1_committed_r_sb", definition = x$def,
                       K, n_sets_used, gate_stat, delta_vs_base_category,
                       delta_boot_support, beta_F4,
                       threshold_met = gate_stat >= GATE_SPEARMAN_MIN,
                       beats_base = delta_vs_base_category > 0)
  })
fd_sim <- purrr::map_dfr(list(
  list(def = "canonical_committed", g = g_sim_can),
  list(def = "alternative_P2_erratum", g = g_sim_alt)), function(x) {
    x$g$per_K |>
      dplyr::transmute(arm = "sim_free105_committed_r_sb (SIMULATION)",
                       definition = x$def,
                       K, n_sets_used, gate_stat, delta_vs_base_category,
                       delta_boot_support, beta_F4,
                       threshold_met = gate_stat >= GATE_SPEARMAN_MIN,
                       beats_base = delta_vs_base_category > 0)
  })
fd_curve <- dplyr::bind_rows(
  main_curve |> dplyr::mutate(definition = "canonical_committed"),
  alt_curve |> dplyr::mutate(definition = "alternative_P2_erratum")) |>
  dplyr::transmute(arm = "planted_detection_curve (SIMULATION)", definition,
                   beta, n_sim, p_gate, p_threshold, p_sign, mc_se_p_gate)
fd_floor <- tibble::tibble(
  arm = "minimum_detectable_effect", definition = c("canonical_committed",
                                                    "alternative_P2_erratum"),
  beta_floor_power80 = c(beta_floor, beta_floor_alt))
feature_definition_robustness <- dplyr::bind_rows(
  fd_real |> dplyr::mutate(section = "real_data_gate", .before = 1),
  fd_sim |> dplyr::mutate(section = "sim_track_gate", .before = 1),
  fd_curve |> dplyr::mutate(section = "detection_curve", .before = 1),
  fd_floor |> dplyr::mutate(section = "detection_floor", .before = 1))
save_csv(feature_definition_robustness, "feature_definition_robustness")

# =============================================================================
# 13. Remaining exports, seeds, config
# =============================================================================
detail_out <- dplyr::bind_rows(detail, detail_alt) |>
  dplyr::mutate(data_kind = "SIMULATION", .before = 1) |>
  dplyr::arrange(family, coll_target, beta, variant, rep)
save_csv(detail_out, "sens_replicate_detail")
save_csv(hp_setlevel |> dplyr::filter(K == 12L), "sens_hp_set_reliability_K12")
save_csv(seven_form |> dplyr::mutate(data_kind = "SIMULATION", .before = 1),
         "sens_sevenform_generalization")
save_csv(dplyr::bind_rows(
  neg_curve |> dplyr::mutate(arm = "negative_control (SIMULATION)", .before = 1),
  rate_curve |> dplyr::mutate(arm = "uncalibrated_rates (SIMULATION)", .before = 1)),
  "sens_control_and_rate_curves")
save_csv(tibble::tibble(GameID = DICT$GameID, empirical_p_choose_B = p_emp,
                        calibration_offset = cal_offset,
                        achieved_rate = cal_achieved),
         "sens_generator_calibration")
save_csv(cal_drift, "sens_hp_rate_drift")
save_csv(frozen_provenance, "sens_frozen_provenance")

seeds_tbl <- dplyr::bind_rows(
  tibble::tibble(purpose = "base", scope = "all", seed = SEED_SENS,
                 derivation = "SEED_SENS"),
  tibble::tibble(purpose = c("rate_calibration", "G_EQ", "G_GEN"),
                 scope = "gates",
                 seed = c(seed_cal, seed_geq, seed_ggen),
                 derivation = c("SEED_SENS + 11", "SEED_SENS + 22",
                                "SEED_SENS + 33")),
  tibble::tibble(
    purpose = c("params", "choices", "negative_carrier", "set_halvings",
                "gate_cv", "gate_bootstrap"),
    scope = "per (arm, rep); resolved in the job list",
    seed = c(seed_par(1L, 1L), seed_choice(1L, 1L), seed_negcar(8L, 1L),
             seed_split(1L, 1L, 8L, 1L), seed_cv(1L, 1L, 8L, 0L),
             seed_boot_g(1L, 1L, 0L)),
    derivation = c(
      "1e8 + arm*1e6 + rep*100 (value shown: arm 1, rep 1)",
      "1e8 + arm*1e6 + rep*100 + 1",
      "1e8 + arm*1e6 + rep*100 + 2 (negative-control arm)",
      "3e8 + arm*3e6 + rep*3e4 + ((K-8)/2)*1e4 + set_id",
      "7e8 + variant*2e8 + arm*3e6 + rep*3e4 + K (variant 0 canonical, 1 alt)",
      "7e8 + variant*2e8 + arm*3e6 + rep*3e4 + 99")),
  tibble::tibble(purpose = "e4_reproduction", scope = "G_E4EQ + committed-data arms",
                 seed = c(e4_seeds[["12"]], SEED_E4 + 500000L + 1100L),
                 derivation = c("SEED_E4 + 1000*track_index + 100*target_index + K (E4's own derivation)",
                                "SEED_E4 + 500000 + 1000*track_index + 100*target_index")))
save_csv(seeds_tbl, "sens_seeds")

cfg <- list(
  spec = "elevation/SENS_spec.md",
  status = paste("SIMULATION detectability calibration; the frozen H2 gate",
                 "result (E4 = FAIL) STANDS and is not re-run,",
                 "amended or replaced"),
  seed_base = SEED_SENS,
  beta_grid = BETA_GRID,
  n_sim = list(main = N_SIM_MAIN, collinearity = N_SIM_COLL,
               negative_control = N_SIM_NEG, rate_sensitivity = N_SIM_RATE,
               high_precision_reps = N_HP_REPS),
  design = list(
    form = 1L, n_subjects = 63L, n_problems = 24L, n_loss_informative = 6L,
    r_reps = R_REPS, n_split = N_SPLIT_SENS, m_sets_per_K = M_SETS,
    K_gate = K_GATE,
    seven_form_note = paste(
      "The frozen E4 gate's empirical anchor is form 1 only; the planted-data",
      "gate replicates that design exactly. The seven-form structure enters",
      "through the committed per-form set machinery: a high-precision",
      "generalization arm scores forms 2-7 at their largest committed budget",
      "(sens_sevenform_generalization.csv). Forms 6-7 violate the common-set",
      "premise on real data and are flagged; their simulated subjects have",
      "full within-form exposure.")),
  generator = list(
    base = "frozen draw_true_params + Binomial(20, p) per problem, p from the frozen analytic_p_b linear predictor (G_GEN-proven)",
    planted_term = "b_loss_ij = MU_L + (1 + beta * g_j) * (b_loss_i - MU_L); MU_L = -0.75 = mean of the frozen U(-2.0, 0.5) b_loss range",
    carrier = "per-item info_loss_avoid_adj centered/scaled over the form's loss-informative pool; 0 on non-loss-informative items",
    rate_calibration = "per-problem offset so the frozen-population mean P(choose B) equals the committed E2 build empirical rate (p_choose_B_exp); solved once, seed_cal; sensitivity arm runs uncalibrated at beta in {0, 0.15}",
    choice_rates_preserved = "exactly at beta = 0 by construction; drift at beta > 0 measured in sens_hp_rate_drift.csv",
    heterogeneity_preserved = "carrier mean 0 over the loss-informative pool, so the average trait dispersion is unchanged",
    noise_model_caveat = paste(
      "Binomial trial noise; P2 measured ~8.5x overdispersion in real per-item",
      "counts, so planted-data detection probabilities are OPTIMISTIC (upper",
      "bounds on detectability). The attenuation bound uses the E3b MEASURED",
      "outcome noise and needs no generator; conclusions lean on it wherever",
      "the two disagree.")),
  detection_event = list(
    C1 = "gate_stat = min(pooled, repeat-mean) oof Spearman >= 0.30 AND boot CI lower > 0, at every gate K and pooled",
    C2 = "(delta vs BASE_category > 0 in Spearman or R^2) AND delta boot support >= 0.90, at every gate K",
    C3 = "F4 carrier (|standardized FULL beta| >= 0.05) with positive sign at every gate K (the planted direction)",
    gate = "C1 and C2 and C3",
    p_threshold = "gate_stat >= 0.30 at every gate K and pooled",
    p_sign = "F4 coefficient positive at every gate K",
    mapping_note = paste(
      "The frozen gate's G-sim/G-real split dissolves in a planted world with",
      "one data source; C1/C2 carry the G-sim + G-beyond structure and C3",
      "carries the G-real sign-consistency structure against the planted",
      "truth, which speaks only through F4. Fixed before results."),
    seed_policy = "fresh derived seeds per replicate (the detection probability marginalizes over fold randomness)"),
  matters_threshold = list(
    definition = "top-vs-bottom F4 quintile difference in high-precision K=12 set reliability, minus the same at beta = 0 (the planted increment over the composition channel)",
    threshold_sb = MATTERS_SB, quantile = QUANT_GRAD),
  gates = list(
    g_controls = list(positive_rule = "P(gate | beta = 0.40 canonical) >= 0.80",
                      negative_rule = "P(gate | beta = 0.40 orthogonalized carrier) <= 0.20",
                      p_positive = p_pos, p_negative = p_neg,
                      p_gate_beta0_baseline = p_base0,
                      positive_pass = g_controls_pos,
                      negative_pass = g_controls_neg,
                      pass = g_controls_pass),
    g_power = list(rule = "P(gate | beta_matters) >= 0.80 under canonical AND alternative feature definitions",
                   beta_matters = beta_matters,
                   p_at_matters_canonical = if (!is.na(beta_matters))
                     p_gate_at(main_curve, beta_matters) else NA_real_,
                   p_at_matters_alternative = if (!is.na(beta_matters))
                     p_gate_at(alt_curve, beta_matters) else NA_real_,
                   pass = g_power_pass),
    g_detect_floor = list(beta_floor_canonical = beta_floor,
                          beta_floor_alternative = beta_floor_alt,
                          verdict = detect_floor_verdict)),
  observed_corr_F2_F4 = as.list(corr_f2f4_observed),
  collinearity_scenarios = lapply(coll_carriers, function(x)
    list(c0 = x$c0, achieved_corr_K12 = x$achieved)),
  equivalence_gates = list(
    G_GEN = list(checks = N_GGEN, max_abs_diff = gg_max),
    G_EQ = list(checks = N_GEQ, max_abs_diff = eq_max, na_mismatches = eq_na),
    G_E4EQ = list(statistics = nrow(e4eq_rows),
                  max_abs_diff = max(e4eq_rows$abs_diff), tolerance = E4EQ_TOL)),
  alternative_features = "P2 erratum s3 RATIFIED per-item analogues, set level = mean of per-item values (F1: squared deviation of d_ev_noloss from the pool loss-informative mean; F2: eloss_a + eloss_b; F3: max pairwise |zc_p * zc_q| on pool-centered z); F4/F5 unchanged (exact decompositions)",
  spec_wording_note = "SENS_spec says 'grouped-by-form cross-validation gate'; the frozen E4 gate is repeated set-level holdout with a problem-cluster bootstrap (E4_spec s4 'grouped; cluster-aware'), anchored on form 1; run as frozen, wording difference reported",
  input_files = as.list(input_md5),
  build_only = TRUE, freeze_read = FALSE, test_read = FALSE,
  raw_data_read = FALSE, participant_choices_read = FALSE,
  optimization_performed = FALSE, real_data_mechanism_fitted = FALSE,
  h2_gate_rerun_or_amended = FALSE, candidate_designs_constructed = FALSE,
  runtime_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
jsonlite::write_json(cfg, file.path(res_dir, "sens_config.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null",
                     digits = NA)

say("\n=== SENS complete in %.1f min ===\n",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")))
say("main detection curve (SIMULATION):\n")
print(as.data.frame(sensitivity_detection_by_beta |>
  dplyr::select(beta, n_sim, p_gate, p_threshold, p_sign,
                induced_reliability_change_K12, r_true_spearman_K12)))
say("beta_matters = %s; beta_floor = %s (alt %s)\n",
    format(beta_matters), format(beta_floor), format(beta_floor_alt))
say("G-controls: positive %.3f (pass %s), negative %.3f (pass %s)\n",
    p_pos, g_controls_pos, p_neg, g_controls_neg)
say("G-power pass: %s | detect-floor verdict: %s\n",
    g_power_pass, detect_floor_verdict)

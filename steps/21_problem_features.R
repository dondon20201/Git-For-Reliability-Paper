# Step 21: problem-information dictionary
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e1/participant_split.csv (step 20)
# Writes:  output/elevation/results/e2/ (problem_information_dictionary.csv, e2_feature_config.json, audits)
# Paper:   Methods (Problem-information dictionary)
# Runtime: a few seconds

# E2: deterministic item-level problem-information dictionary.
# Each field is a closed-form function of the raw data restricted to build
# subjects. There is no random sampling or bootstrap, and no reliability,
# recovery, prediction, optimization, or mechanism model is fitted.
#
# Choice-dependent fields use the 344 build subjects from step 20. The 137
# freeze and 205 test subjects' choices are excluded; test choices remain
# reserved for E7. Choice/outcome columns are dropped from the full-sample
# frame at load, so only build_trials contains choices.
#
# The loaders in R/cpc18_features.R keep the structure and predictor fields
# consistent with earlier steps.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

source(here("R", "cpc18_features.R"))

res_dir <- here("output", "elevation", "results", "e2")
log_dir <- here("output", "elevation", "logs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(res_dir, paste0(name, ".csv")))
  invisible(df)
}

RAW_MD5 <- "38a26407a36c690be969f9d9d3fe3e5a"
EXP_LO <- 6L; EXP_HI <- 25L      # experience regime (primary)
DESC_LO <- 1L; DESC_HI <- 5L     # description regime (secondary generalization)

# ---- inputs -----------------------------------------------------------------
raw_path <- here("data", "raw", "All estimation raw data.csv")
stopifnot(file.exists(raw_path))
raw_md5 <- unname(tools::md5sum(raw_path))
stopifnot(identical(raw_md5, RAW_MD5))              # verified data only
cat(sprintf("raw md5 verified: %s\n", raw_md5))

split_path <- here("output", "elevation", "results", "e1", "participant_split.csv")
stopifnot(file.exists(split_path))
split_tbl <- readr::read_csv(split_path, show_col_types = FALSE)
stopifnot(identical(names(split_tbl),
                    c("SubjID", "Set", "split", "held_out_until_E7")))

build_ids  <- split_tbl$SubjID[split_tbl$split == "build"]
n_build    <- length(build_ids)
n_freeze   <- sum(split_tbl$split == "freeze")   # counted, choices never read
n_test     <- sum(split_tbl$split == "test")     # counted, choices never read
stopifnot(n_build == 344L, n_freeze == 137L, n_test == 205L)
stopifnot(!anyDuplicated(build_ids))
# the held-out flag must partition exactly as the E1 split set it
stopifnot(all(split_tbl$held_out_until_E7 == (split_tbl$split == "test")))
cat(sprintf("split read: build %d / freeze %d / test %d\n",
            n_build, n_freeze, n_test))

# ---- load raw, then structurally separate design from choices ---------------
raw <- load_cpc18(raw_path)
stopifnot(all(cpc18_expected_columns %in% names(raw)))

# Every per-trial choice / outcome column. Dropped from the full-sample frame so
# that no freeze or test choice can reach any computation even by accident.
choice_cols <- c("B", "Button", "Payoff", "Forgone", "RT", "Apay", "Bpay",
                 "Feedback")

# design_all: ALL subjects, but choice-free (problem descriptions + form ids).
design_all <- raw |> dplyr::select(-dplyr::any_of(choice_cols))
stopifnot(!any(choice_cols %in% names(design_all)))

# build_trials: the ONLY frame carrying choices, restricted to build subjects
# BEFORE any choice-based field is computed.
build_trials <- raw |> dplyr::filter(SubjID %in% build_ids)
rm(raw); invisible(gc())

stopifnot(setequal(unique(build_trials$SubjID), build_ids))
stopifnot(length(unique(build_trials$SubjID)) == 344L)
# no freeze/test subject in the choice-bearing frame
held_out_ids <- split_tbl$SubjID[split_tbl$split != "build"]
stopifnot(!any(build_trials$SubjID %in% held_out_ids))

# ---- clean-problem universe (rows of the dictionary) ------------------------
clean_all   <- design_all |> dplyr::filter(is_clean_row(design_all))
clean_build <- build_trials |> dplyr::filter(is_clean_row(build_trials))

universe_ids   <- sort(unique(clean_build$GameID))   # clean AND build-visible
full_clean_ids <- sort(unique(clean_all$GameID))     # clean, whole sample
unseen_by_build <- setdiff(full_clean_ids, universe_ids)

stopifnot(length(universe_ids) == 105L)              # expected clean set size
stopifnot(setequal(universe_ids, full_clean_ids))    # build sees the full set
stopifnot(length(unseen_by_build) == 0L)
cat(sprintf("universe: %d clean problems; %d unseen by build\n",
            length(universe_ids), length(unseen_by_build)))

# ---- B. structure + C. predictor (shared loaders; subject-independent) ------
feat <- clean_problem_features(design_all) |> add_scaled_features()
scales <- attr(feat, "scales")
stopifnot(setequal(feat$GameID, universe_ids))
feat <- feat |> dplyr::arrange(GameID)

# near_parity: reproduces near_parity_sensitivity from step 04, medians over the
# clean-problem universe.
med_abs_d_ev <- stats::median(abs(feat$d_ev))
med_abs_d_sd <- stats::median(abs(feat$d_sd))

# domain: literal on the four outcome values (same convention as
# has_loss, which ignores probability).
outcome_min <- pmin(feat$Ha, feat$La, feat$Hb, feat$Lb)
outcome_max <- pmax(feat$Ha, feat$La, feat$Hb, feat$Lb)

# First-order stochastic dominance between two two-outcome options. F is a step
# function with jumps at the outcome values, so comparing the CDFs at the
# distinct outcome values is exactly sufficient (F is constant between jumps and
# both CDFs are 0 to the left of the smallest outcome).
fosd <- function(h_a, p_a, l_a, h_b, p_b, l_b, tol = 1e-12) {
  cdf <- function(t, h, p, l) (l <= t) * (1 - p) + (h <= t) * p
  ts <- unique(c(h_a, l_a, h_b, l_b))
  fa <- vapply(ts, cdf, numeric(1), h = h_a, p = p_a, l = l_a)
  fb <- vapply(ts, cdf, numeric(1), h = h_b, p = p_b, l = l_b)
  d <- fa - fb
  a_dom <- all(d <= tol) && any(d < -tol)   # F_A <= F_B everywhere, strict once
  b_dom <- all(d >= -tol) && any(d > tol)
  if (a_dom) "A" else if (b_dom) "B" else "none"
}
dominance <- vapply(seq_len(nrow(feat)), function(i) {
  fosd(feat$Ha[i], feat$pHa[i], feat$La[i],
       feat$Hb[i], feat$pHb[i], feat$Lb[i])
}, character(1))

feat <- feat |>
  dplyr::mutate(
    domain = dplyr::case_when(
      outcome_min >= 0 ~ "gain",
      outcome_max <= 0 ~ "loss",
      TRUE ~ "mixed"
    ),
    near_parity = abs(d_ev) <= med_abs_d_ev & abs(d_sd) >= med_abs_d_sd,
    dominance = dominance,
    abs_d_ev = abs(d_ev), abs_d_sd = abs(d_sd), abs_d_eloss = abs(d_eloss),
    # EV(no-loss option) - EV(loss option); same convention as R/step1_measures.R.
    d_ev_noloss = dplyr::if_else(loss_informative,
                                 dplyr::if_else(has_loss_b, ev_a - ev_b,
                                                ev_b - ev_a),
                                 NA_real_)
  )

# ---- A. identity and provenance ---------------------------------------------
# Form membership is a design fact, so it is read from the choice-free frame.
set_membership <- clean_all |>
  dplyr::distinct(GameID, Set) |>
  dplyr::mutate(present = TRUE) |>
  tidyr::pivot_wider(names_from = Set, values_from = present,
                     names_prefix = "in_Set", values_fill = FALSE) |>
  dplyr::select(GameID, dplyr::num_range("in_Set", 1:7)) |>
  dplyr::arrange(GameID)
set_cols <- paste0("in_Set", 1:7)
stopifnot(all(set_cols %in% names(set_membership)))
set_membership$n_forms <- rowSums(as.matrix(set_membership[set_cols]))

# The build-restricted form map must agree (every form is in the build split).
set_membership_build <- clean_build |>
  dplyr::distinct(GameID, Set) |>
  dplyr::mutate(present = TRUE) |>
  tidyr::pivot_wider(names_from = Set, values_from = present,
                     names_prefix = "in_Set", values_fill = FALSE) |>
  dplyr::select(GameID, dplyr::num_range("in_Set", 1:7)) |>
  dplyr::arrange(GameID)
stopifnot(identical(as.data.frame(set_membership[c("GameID", set_cols)]),
                    as.data.frame(set_membership_build[c("GameID", set_cols)])))

build_counts <- clean_build |>
  dplyr::group_by(GameID) |>
  dplyr::summarise(
    n_build_subjects = dplyr::n_distinct(SubjID),
    n_build_trials_exp = sum(Trial >= EXP_LO & Trial <= EXP_HI),
    n_build_trials_desc = sum(Trial >= DESC_LO & Trial <= DESC_HI),
    .groups = "drop"
  ) |>
  dplyr::arrange(GameID)

# ---- D. entropy (BUILD choices; deterministic) ------------------------------
# Binary entropy in bits; 0 at p in {0, 1}.
binary_entropy <- function(p) {
  ifelse(is.na(p) | p <= 0 | p >= 1, ifelse(is.na(p), NA_real_, 0),
         -p * log2(p) - (1 - p) * log2(1 - p))
}

clean_build_feat <- clean_build |>
  dplyr::select(SubjID, GameID, Trial, B) |>
  dplyr::left_join(
    feat |> dplyr::select(GameID, z_ev, z_sd, z_eloss, loss_informative,
                          has_loss_b, d_ev_noloss),
    by = "GameID"
  )
stopifnot(!any(is.na(clean_build_feat$z_ev)))

exp_trials  <- clean_build_feat |>
  dplyr::filter(Trial >= EXP_LO, Trial <= EXP_HI)
desc_trials <- clean_build_feat |>
  dplyr::filter(Trial >= DESC_LO, Trial <= DESC_HI)

p_exp <- exp_trials |>
  dplyr::group_by(GameID) |>
  dplyr::summarise(p_choose_B_exp = mean(B), .groups = "drop")
p_desc <- desc_trials |>
  dplyr::group_by(GameID) |>
  dplyr::summarise(p_choose_B_desc = mean(B), .groups = "drop")

# chose_noloss on loss-informative problems (same convention as R/step1_measures.R).
p_nl_exp <- exp_trials |>
  dplyr::filter(loss_informative) |>
  dplyr::mutate(chose_noloss = ifelse(has_loss_b, 1 - B, B)) |>
  dplyr::group_by(GameID) |>
  dplyr::summarise(p_choose_noloss_exp = mean(chose_noloss), .groups = "drop")

# ---- BUILD reference operating point (deterministic; exported) --------------
# theta_bar_build: pooled logistic on ALL build clean experience trials, one
# Bernoulli row per trial. A population-average reference for EXPECTED
# information only; it is not a per-subject estimate and carries no reliability
# claim. Pooling ~344 subjects x problems x 20 trials avoids the per-subject
# quasi-separation flagged in E0 (step 18).
fit_ref <- stats::glm(B ~ z_ev + z_sd + z_eloss, family = stats::binomial(),
                      data = exp_trials)
theta_bar_build <- stats::coef(fit_ref)
stopifnot(fit_ref$converged, all(is.finite(theta_bar_build)))

# Build loss baseline: the build-only analogue of fit_loss_baselines()$exp in
# R/step2_reliability.R, used ONLY to weight the loss information atom. It does
# not replace that baseline or the reliability engine.
loss_exp_trials <- exp_trials |>
  dplyr::filter(loss_informative) |>
  dplyr::mutate(chose_noloss = ifelse(has_loss_b, 1 - B, B))
fit_loss <- stats::glm(chose_noloss ~ d_ev_noloss, family = stats::binomial(),
                       data = loss_exp_trials)
coef_loss <- stats::coef(fit_loss)
stopifnot(fit_loss$converged, all(is.finite(coef_loss)))

cat(sprintf("theta_bar_build (n=%d): %s\n", nrow(exp_trials),
            paste(sprintf("%s=%.6f", names(theta_bar_build), theta_bar_build),
                  collapse = ", ")))
cat(sprintf("build loss baseline (n=%d): %s\n", nrow(loss_exp_trials),
            paste(sprintf("%s=%.6f", names(coef_loss), coef_loss),
                  collapse = ", ")))

# ---- E. Fisher information at the build reference operating point -----------
dict <- feat |>
  dplyr::left_join(set_membership, by = "GameID") |>
  dplyr::left_join(build_counts, by = "GameID") |>
  dplyr::left_join(p_exp, by = "GameID") |>
  dplyr::left_join(p_desc, by = "GameID") |>
  dplyr::left_join(p_nl_exp, by = "GameID") |>
  dplyr::mutate(
    choice_entropy_exp = binary_entropy(p_choose_B_exp),
    choice_entropy_desc = binary_entropy(p_choose_B_desc),
    # p_ref = plogis(x . theta_bar_build), x = (1, z_ev, z_sd, z_eloss)
    p_ref_exp = stats::plogis(theta_bar_build[["(Intercept)"]] +
                              theta_bar_build[["z_ev"]] * z_ev +
                              theta_bar_build[["z_sd"]] * z_sd +
                              theta_bar_build[["z_eloss"]] * z_eloss),
    fisher_w_exp = p_ref_exp * (1 - p_ref_exp),
    finfo_b_ev = fisher_w_exp * z_ev^2,
    finfo_b_risk = fisher_w_exp * z_sd^2,
    finfo_b_loss = fisher_w_exp * z_eloss^2,
    p_noloss_ref_exp = dplyr::if_else(
      loss_informative,
      stats::plogis(coef_loss[["(Intercept)"]] +
                    coef_loss[["d_ev_noloss"]] * d_ev_noloss),
      NA_real_),
    w_noloss_exp = p_noloss_ref_exp * (1 - p_noloss_ref_exp),
    # Design proxy for how much a problem informs the PRIMARY target. A problem
    # that is not loss-informative contributes ZERO loss information (not an
    # undefined quantity), so set-level sums stay well defined; the genuinely
    # undefined fields (d_ev_noloss, p_choose_noloss_exp, w_noloss_exp) stay NA.
    info_loss_avoid_adj = dplyr::if_else(loss_informative, w_noloss_exp, 0)
  )

dict <- dict |>
  dplyr::select(
    # A. identity and provenance
    GameID, dplyr::all_of(set_cols), n_forms,
    n_build_subjects, n_build_trials_exp, n_build_trials_desc,
    # B. structure
    Ha, pHa, La, Hb, pHb, Lb,
    ev_a, sd_a, eloss_a, has_loss_a, ev_b, sd_b, eloss_b, has_loss_b,
    domain, loss_informative, near_parity, dominance,
    # C. predictor
    d_ev, d_sd, d_eloss, z_ev, z_sd, z_eloss,
    abs_d_ev, abs_d_sd, abs_d_eloss, d_ev_noloss,
    # D. entropy (build choices)
    p_choose_B_exp, choice_entropy_exp, p_choose_noloss_exp,
    p_choose_B_desc, choice_entropy_desc,
    # E. Fisher information at the build reference operating point
    p_ref_exp, fisher_w_exp, finfo_b_ev, finfo_b_risk, finfo_b_loss,
    p_noloss_ref_exp, w_noloss_exp, info_loss_avoid_adj
  ) |>
  dplyr::arrange(GameID)

# ---- in-script assertions (pre-declared gates) ------------------------------
# G1 universe and shape
stopifnot(nrow(dict) == 105L)
stopifnot(!anyDuplicated(dict$GameID))
stopifnot(setequal(dict$GameID, full_clean_ids))
stopifnot(all(dict$n_build_subjects >= 1L))
stopifnot(all(dict$n_forms >= 1L))
# every experience-regime trial count is 20 per build subject who saw the item
stopifnot(all(dict$n_build_trials_exp == dict$n_build_subjects * 20L))
stopifnot(all(dict$n_build_trials_desc == dict$n_build_subjects * 5L))
# G3 structure/predictor internal consistency
stopifnot(all(abs(dict$d_ev - (dict$ev_b - dict$ev_a)) < 1e-12))
stopifnot(all(abs(dict$d_sd - (dict$sd_b - dict$sd_a)) < 1e-12))
stopifnot(all(abs(dict$d_eloss - (dict$eloss_b - dict$eloss_a)) < 1e-12))
stopifnot(all(abs(dict$z_ev - dict$d_ev / scales[["d_ev"]]) < 1e-12))
stopifnot(all(abs(dict$z_sd - dict$d_sd / scales[["d_sd"]]) < 1e-12))
stopifnot(all(abs(dict$z_eloss - dict$d_eloss / scales[["d_eloss"]]) < 1e-12))
stopifnot(identical(dict$loss_informative,
                    xor(dict$has_loss_a, dict$has_loss_b)))
stopifnot(all(dict$domain %in% c("gain", "loss", "mixed")))
stopifnot(all(dict$dominance %in% c("A", "B", "none")))
# G4 entropy / proportion
stopifnot(all(dict$p_choose_B_exp >= 0 & dict$p_choose_B_exp <= 1))
stopifnot(all(dict$choice_entropy_exp >= 0 & dict$choice_entropy_exp <= 1))
stopifnot(all(dict$choice_entropy_desc >= 0 & dict$choice_entropy_desc <= 1))
stopifnot(all(is.finite(dict$p_choose_B_exp)))
# loss-only fields are present exactly on loss-informative rows
stopifnot(identical(!is.na(dict$d_ev_noloss), dict$loss_informative))
stopifnot(identical(!is.na(dict$p_choose_noloss_exp), dict$loss_informative))
stopifnot(identical(!is.na(dict$w_noloss_exp), dict$loss_informative))
stopifnot(all(!is.na(dict$info_loss_avoid_adj)))
stopifnot(all(dict$info_loss_avoid_adj[!dict$loss_informative] == 0))
# G5 Fisher weights are proper Bernoulli variances
stopifnot(all(dict$fisher_w_exp > 0 & dict$fisher_w_exp <= 0.25))
stopifnot(all(dict$w_noloss_exp > 0 & dict$w_noloss_exp <= 0.25,
              na.rm = TRUE))
# G7 rows are per-problem: no subject identifier and no freeze/test column, so
# the output carries no subject-level leakage vector. (n_build_subjects is an
# aggregate count over build subjects only, not a subject-level field.)
stopifnot(!("SubjID" %in% names(dict)))
stopifnot(!any(grepl("freeze|test", names(dict), ignore.case = TRUE)))
stopifnot(identical(anyDuplicated(dict$GameID), 0L))   # one row per problem

cat(sprintf("dictionary: %d rows x %d columns; all assertions passed\n",
            nrow(dict), ncol(dict)))

# ---- outputs ----------------------------------------------------------------
save_csv(dict, "problem_information_dictionary")

build_reference <- dplyr::bind_rows(
  tibble::tibble(
    fit = "theta_bar_build",
    formula = "B ~ z_ev + z_sd + z_eloss",
    regime = sprintf("experience (Trial %d-%d)", EXP_LO, EXP_HI),
    term = names(theta_bar_build),
    estimate = as.numeric(theta_bar_build),
    n_trials = nrow(exp_trials),
    n_subjects = dplyr::n_distinct(exp_trials$SubjID),
    n_problems = dplyr::n_distinct(exp_trials$GameID),
    converged = fit_ref$converged
  ),
  tibble::tibble(
    fit = "build_loss_baseline",
    formula = "chose_noloss ~ d_ev_noloss",
    regime = sprintf("experience (Trial %d-%d), loss-informative", EXP_LO, EXP_HI),
    term = names(coef_loss),
    estimate = as.numeric(coef_loss),
    n_trials = nrow(loss_exp_trials),
    n_subjects = dplyr::n_distinct(loss_exp_trials$SubjID),
    n_problems = dplyr::n_distinct(loss_exp_trials$GameID),
    converged = fit_loss$converged
  )
)
save_csv(build_reference, "e2_build_reference")

per_set_clean <- vapply(set_cols, function(cc) sum(dict[[cc]]), integer(1))
universe_audit <- dplyr::bind_rows(
  tibble::tibble(
    check = c("n_clean_problems", "n_clean_problems_build_visible",
              "n_clean_problems_unseen_by_build",
              "n_build_subjects", "n_freeze_subjects", "n_test_subjects",
              "n_build_clean_trials_exp", "n_build_clean_trials_desc",
              "median_abs_d_ev", "median_abs_d_sd",
              "scale_sd_d_ev", "scale_sd_d_sd", "scale_sd_d_eloss",
              "n_loss_informative", "n_near_parity", "n_dominated"),
    value = c(length(full_clean_ids), length(universe_ids),
              length(unseen_by_build),
              n_build, n_freeze, n_test,
              nrow(exp_trials), nrow(desc_trials),
              med_abs_d_ev, med_abs_d_sd,
              scales[["d_ev"]], scales[["d_sd"]], scales[["d_eloss"]],
              sum(dict$loss_informative), sum(dict$near_parity),
              sum(dict$dominance != "none")),
    note = c("LotNumA==1 & LotNumB==1 & Amb==0",
             "clean GameIDs seen by >= 1 build subject",
             "must be 0",
             "E1 split, choices read", "E1 split, choices NEVER read",
             "E1 split, sealed until E7 - choices NEVER read",
             sprintf("Trial %d-%d", EXP_LO, EXP_HI),
             sprintf("Trial %d-%d", DESC_LO, DESC_HI),
             "near_parity median over the clean universe",
             "near_parity median over the clean universe",
             "frozen add_scaled_features convention",
             "frozen add_scaled_features convention",
             "frozen add_scaled_features convention",
             "xor(has_loss_a, has_loss_b)",
             "reproduces Study-1 near_parity_sensitivity",
             "first-order stochastic dominance, either direction")
  ),
  tibble::tibble(
    check = paste0("n_clean_problems_", set_cols),
    value = as.numeric(per_set_clean),
    note = "clean problems appearing in this fixed form"
  )
)
save_csv(universe_audit, "e2_universe_audit")

cfg <- list(
  spec = "elevation/E2_spec.md",
  raw_md5 = raw_md5,
  regimes = list(experience = c(EXP_LO, EXP_HI), description = c(DESC_LO, DESC_HI)),
  clean_rule = "LotNumA == 1 & LotNumB == 1 & Amb == 0",
  n_clean_problems = length(universe_ids),
  n_clean_problems_unseen_by_build = length(unseen_by_build),
  scaling_constants = list(d_ev = scales[["d_ev"]], d_sd = scales[["d_sd"]],
                           d_eloss = scales[["d_eloss"]]),
  near_parity_medians = list(abs_d_ev = med_abs_d_ev, abs_d_sd = med_abs_d_sd),
  split_counts = list(build = n_build, freeze = n_freeze, test = n_test),
  theta_bar_build = as.list(theta_bar_build),
  theta_bar_build_n_trials = nrow(exp_trials),
  build_loss_baseline = as.list(coef_loss),
  build_loss_baseline_n_trials = nrow(loss_exp_trials),
  info_loss_avoid_adj_convention =
    "0 on non-loss-informative rows (contributes no loss information, so set-level sums are defined); d_ev_noloss / p_choose_noloss_exp / w_noloss_exp stay NA there",
  data_inputs = c("data/raw/All estimation raw data.csv",
                  "elevation/results/e1/participant_split.csv"),
  input_provenance = paste(
    "Choice-based fields and both reference fits read ONLY the 344 build",
    "subjects. The full-sample frame has every choice/outcome column dropped at",
    "load, so no freeze (137) or test (205) choice can enter any field. The 205",
    "test subjects stay sealed until E7."),
  no_freeze_or_test_choice_read = TRUE,
  computes_reliability_recovery_or_prediction = FALSE,
  deterministic = TRUE,
  rng_used = FALSE,
  seed = NULL
)
jsonlite::write_json(cfg, file.path(res_dir, "e2_feature_config.json"),
                     pretty = TRUE, auto_unbox = TRUE, null = "null")

cat("\nE2 dictionary written to elevation/results/e2/\n")
cat(sprintf("PASS: %d clean problems, all build-visible; build-only (%d subjects); deterministic (no RNG).\n",
            nrow(dict), n_build))
cat("PASS: no reliability / recovery / prediction computed; test group sealed.\n")

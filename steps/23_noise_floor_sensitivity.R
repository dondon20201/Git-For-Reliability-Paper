# Step 23: sensitivity of the subset noise floor to the split draw
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e1/, e3/ (steps 20, 22)
# Writes:  output/elevation/results/e3/e3_split_noise_by_set.csv, output/elevation/results/e3/e3_noise_floor_sensitivity.csv
# Paper:   Study 2 (input to the mechanism model in step 25)
# Runtime: about 20 seconds with 12 parallel workers

# E3 noise-floor SENSITIVITY (diagnostic, reported alongside the
# pre-declared gate; it does NOT replace it).
#
# The pre-declared gate defines sigma2_noise via the reliability_point
# subject-cluster bootstrap, taken to capture subject-sampling + split
# noise for a single set. That bootstrap resamples SUBJECTS while holding the
# 100 halvings FIXED, so the split-draw component of a set's r_sb is not in it.
# Any missing noise inflates sigma2_set_true = sigma2_set_raw - sigma2_noise,
# which biases the H1 gate TOWARD passing. This script measures the split-draw
# component directly by recomputing r_sb for the same sets under independent
# split seeds, and reports what the gate would be with the augmented floor.
#
# RESULT: the split-draw component is small (sigma2 = 1.6e-4 at the gate cell,
# 0.8% of the bootstrap noise floor), because r_sb already averages over 100
# halvings. The H1 verdict is unchanged (ICC 0.6258 -> 0.6230), so this
# omission does not carry the PASS.
#
# The OPPOSITE bias is the material one and is NOT corrected here: the
# subject-cluster bootstrap
# measures how much a set's r_sb moves when SUBJECTS are resampled, but all
# sets within a form are scored on the SAME subjects, so that component is
# largely COMMON across sets and does not contribute to across-set variance.
# Subtracting it therefore over-subtracts, which is why 32 of 50 cells return
# sigma2_set_raw - sigma2_noise < 0 (clamped to 0). Reported sigma_set_true is
# consequently a LOWER bound and the H1 PASS is conservative.
#
# This script re-uses the vectorised target path by copy (it is a standalone
# diagnostic, not an independent re-derivation) and produces no headline number;
# split_correlations / sb_correct from R/step2_reliability.R are called unchanged.

suppressPackageStartupMessages({ library(tidyverse); library(here); library(parallel) })
source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))

RAW_MD5 <- "38a26407a36c690be969f9d9d3fe3e5a"
SEED_E3 <- 20260719L; POS <- 6:25; R_REPS <- 20L
TARGETS <- c("loss_avoid_adj", "ev_sens_agn")
N_SPLIT_A <- 100L; MIN_PAIRS <- 10L; M_CAP <- 2000L
GATE_FORM <- 1L; GATE_K <- 12L
N_SETS_DIAG <- 50L      # the same 50-set subsample the noise floor uses
N_REPLICATES <- 8L      # independent split seeds per set

raw_path <- here("data", "raw", "All estimation raw data.csv")
stopifnot(identical(unname(tools::md5sum(raw_path)), RAW_MD5))
split_tbl <- readr::read_csv(here("output", "elevation", "results", "e1",
                                 "participant_split.csv"), show_col_types = FALSE)
build_ids <- split_tbl$SubjID[split_tbl$split == "build"]
stopifnot(length(build_ids) == 344L)

raw <- load_cpc18(raw_path)
feat_all <- clean_problem_features(raw) |> add_scaled_features() |> dplyr::arrange(GameID)
build_raw <- raw |> dplyr::filter(SubjID %in% build_ids); rm(raw); invisible(gc())
build_trials <- build_raw |>
  dplyr::filter(is_clean_row(build_raw)) |>
  dplyr::select(SubjID, Set, GameID, Trial, B) |>
  dplyr::left_join(feat_all |> dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev,
                                             loss_informative, has_loss_b, ev_a, ev_b),
                   by = "GameID")
baselines <- fit_loss_baselines(build_trials)
subj_data <- build_subject_data(build_trials, baselines)
feat_all <- feat_all |>
  dplyr::mutate(d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a),
                p_base_exp = ifelse(loss_informative,
                                    stats::plogis(baselines$exp[1] +
                                                    baselines$exp[2] * d_ev_noloss), NA))

probs <- sort(unique(build_trials$GameID[build_trials$Set == GATE_FORM]))
subs <- as.character(sort(unique(build_trials$SubjID[build_trials$Set == GATE_FORM])))
NB <- matrix(NA_real_, length(subs), length(probs), dimnames = list(subs, NULL))
for (sid in subs) {
  sd <- subj_data[[sid]]; j <- match(sd$feat$GameID, probs); ok <- !is.na(j)
  NB[sid, j[ok]] <- rowSums(sd$B_mat[which(ok), POS, drop = FALSE])
}
fi <- feat_all[match(probs, feat_all$GameID), ]
fm <- list(probs = probs, li = fi$loss_informative, has_loss_b = fi$has_loss_b,
           d_ev = fi$d_ev, p_base = fi$p_base_exp)

ev_ranks <- function(d_ev, r) {
  k <- length(d_ev); o <- order(d_ev); dv <- d_ev[o]
  grp <- cumsum(c(TRUE, diff(dv) != 0)); sizes <- tabulate(grp)
  starts <- cumsum(c(0, sizes[-length(sizes)])) * r
  Rg <- starts + (sizes * r + 1) / 2
  R <- numeric(k); R[o] <- Rg[grp]; R
}
half_targets <- function(NBh, fm, cols, r) {
  k <- length(cols); n <- r * k
  R <- ev_ranks(fm$d_ev[cols], r); sR <- sum(R); sR2 <- sum(R * R)
  n1 <- rowSums(NBh)
  Sxy <- as.vector(NBh %*% R) - n1 * sR / k
  Sxx <- n1 - n1 * n1 / n; Syy <- r * (sR2 - sR * sR / k)
  ev <- rep(NA_real_, nrow(NBh)); ok <- Sxx > 0 & Syy > 0
  if (any(ok)) ev[ok] <- Sxy[ok] / sqrt(Sxx[ok] * Syy)
  lic <- cols[fm$li[cols]]
  if (length(lic) == 0) loss <- rep(NA_real_, nrow(NBh)) else {
    NBl <- NBh[, fm$li[cols], drop = FALSE]
    cnt <- matrix(0, nrow(NBl), ncol(NBl))
    for (j in seq_along(lic)) cnt[, j] <- if (fm$has_loss_b[lic[j]]) r - NBl[, j] else NBl[, j]
    loss <- (rowSums(cnt) - r * sum(fm$p_base[lic])) / (r * length(lic))
  }
  list(ev_sens_agn = ev, loss_avoid_adj = loss)
}
run_set <- function(cols, seed) {
  NBs <- NB[, cols, drop = FALSE]
  keep <- which(!is.na(rowSums(NBs))); NBs <- NBs[keep, , drop = FALSE]
  K <- length(cols); kh <- K %/% 2L
  E <- list(); for (tg in TARGETS) E[[tg]] <- list(
    h1 = matrix(NA_real_, length(keep), N_SPLIT_A),
    h2 = matrix(NA_real_, length(keep), N_SPLIT_A))
  set.seed(seed)
  for (s in seq_len(N_SPLIT_A)) {
    take <- sample.int(K, 2L * kh)
    i1 <- take[seq_len(kh)]; i2 <- take[(kh + 1L):(2L * kh)]
    a <- half_targets(NBs[, i1, drop = FALSE], fm, cols[i1], R_REPS)
    b <- half_targets(NBs[, i2, drop = FALSE], fm, cols[i2], R_REPS)
    for (tg in TARGETS) { E[[tg]]$h1[, s] <- a[[tg]]; E[[tg]]$h2[, s] <- b[[tg]] }
  }
  vapply(TARGETS, function(tg) {
    rr <- mean(split_correlations(E[[tg]]$h1, E[[tg]]$h2, MIN_PAIRS), na.rm = TRUE)
    if (is.finite(rr)) sb_correct(rr) else NA_real_
  }, numeric(1))
}
draw_sets <- function(np, K, M, seed) {
  set.seed(seed); total <- choose(np, K)
  if (total <= M_CAP) { cmb <- utils::combn(np, K)
    return(lapply(seq_len(ncol(cmb)), function(i) sort(cmb[, i]))) }
  seen <- new.env(hash = TRUE, parent = emptyenv()); out <- vector("list", M); got <- 0L
  while (got < M) { s <- sort(sample.int(np, K)); key <- paste(s, collapse = ",")
    if (is.null(seen[[key]])) { assign(key, TRUE, envir = seen); got <- got + 1L; out[[got]] <- s } }
  out
}

out_all <- list()
for (K in c(8L, 12L)) {
  sets <- draw_sets(length(probs), K, min(M_CAP, choose(length(probs), K)),
                    SEED_E3 + 1000L * GATE_FORM + K)
  sets_sub <- sets[seq_len(min(N_SETS_DIAG, length(sets)))]
  base_seed <- SEED_E3 + 1000L * GATE_FORM + K + 500000L
  cl <- makeCluster(12)
  clusterEvalQ(cl, { suppressPackageStartupMessages({library(tidyverse); library(here)})
                     source(here("R", "step2_reliability.R")); NULL })
  clusterExport(cl, c("NB", "fm", "half_targets", "ev_ranks", "run_set", "TARGETS",
                      "R_REPS", "MIN_PAIRS", "N_SPLIT_A", "base_seed", "N_REPLICATES", "sets_sub"),
                envir = environment())
  res <- parLapply(cl, seq_along(sets_sub), function(i) {
    m <- t(vapply(seq_len(N_REPLICATES), function(rep)
      run_set(sets_sub[[i]], base_seed + i + rep * 7919L), numeric(length(TARGETS))))
    colnames(m) <- TARGETS
    tibble::tibble(set_id = i,
                   target = TARGETS,
                   mean_r_sb = apply(m, 2, mean, na.rm = TRUE),
                   var_across_split_seeds = apply(m, 2, stats::var, na.rm = TRUE),
                   n_finite = apply(m, 2, function(z) sum(is.finite(z))))
  })
  stopCluster(cl)
  out_all[[length(out_all) + 1L]] <- dplyr::bind_rows(res) |> dplyr::mutate(K = K)
  cat(sprintf("K=%d done\n", K))
}
diag_tbl <- dplyr::bind_rows(out_all)

het <- readr::read_csv(here("output", "elevation", "results", "e3",
                            "heterogeneity_summary.csv"), show_col_types = FALSE)
summ <- diag_tbl |>
  dplyr::group_by(K, target) |>
  dplyr::summarise(n_sets = dplyr::n(),
                   sigma2_split_draw = mean(var_across_split_seeds, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::left_join(het |> dplyr::filter(form == GATE_FORM) |>
                     dplyr::select(K, target, sigma2_set_raw, sigma2_noise,
                                   sigma2_set_true, ICC_set),
                   by = c("K", "target")) |>
  dplyr::mutate(
    sigma_split_draw = sqrt(sigma2_split_draw),
    sigma2_noise_augmented = sigma2_noise + sigma2_split_draw,
    sigma2_set_true_augmented = pmax(0, sigma2_set_raw - sigma2_noise_augmented),
    sigma_set_true_augmented = sqrt(sigma2_set_true_augmented),
    ICC_set_augmented = sigma2_set_true_augmented / sigma2_set_raw,
    h1_pass_predeclared = sigma_set_true_augmented >= 0 &
      sqrt(sigma2_set_true) >= 0.05 & ICC_set >= 0.50,
    h1_pass_augmented = sigma_set_true_augmented >= 0.05 & ICC_set_augmented >= 0.50)
print(as.data.frame(summ), digits = 5)
readr::write_csv(diag_tbl, here("output", "elevation", "results", "e3",
                                "e3_split_noise_by_set.csv"))
readr::write_csv(summ, here("output", "elevation", "results", "e3",
                            "e3_noise_floor_sensitivity.csv"))
cat("\nwritten: e3_split_noise_by_set.csv, e3_noise_floor_sensitivity.csv\n")

# Step 38: within-form robustness of the variance decomposition and split-half reliability
# Reads:   data/raw/All estimation raw data.csv; output/elevation/results/e0/ (step 18);
#          output/results/step2/, step4/ (steps 05, 11); steps/18_variance_ceiling.R (functions)
# Writes:  output/elevation/results/audit/ (b1_between_form_anova.csv, b2_within_form_decomposition.csv,
#          b3_within_form_splithalf.csv, reproduction checks, audit_config.json)
# Paper:   Results (pooling across forms); Study 2 (Locating the ceiling); Supplement (Within-form robustness)
# Runtime: about 50 minutes with 12 parallel workers

# Within-form robustness of the reliability results (adds results only; changes
# no method or threshold).
#
# The CPC18 estimation data (md5-verified) were collected
# in seven fixed problem forms (raw column `Set`, 1-7, constant within
# participant). Every reliability the manuscript reports pools participants
# across forms. This script asks whether that pooling matters, by re-reading the
# frozen engines and recomputing within-form counterparts:
#
#   B1  Between-form share. Per-subject full-data point estimates for the five
#       constructs (frozen compute_subject_point_var, B_inner 300, the E0 inner
#       seeds from step 18) -> one-way ANOVA by form: F, p, form means and Ns,
#       and the between-form share of observed variance (SS_between /
#       SS_total). Kruskal-Wallis is added as a labelled rank-based supplement.
#   B2  Within-form variance decomposition for ev_sens_agn and loss_avoid_adj:
#       s2_within = pooled within-form variance of the points (df n - 7),
#       s2_err unchanged (mean within-subject sampling variance),
#       s2_true_within = max(0, s2_within - s2_err); implied reliability at 14
#       problems and k for reliability .7, each with k_ref = 14 (the frozen
#       median) and with k_ref = harmonic mean of per-subject problem counts.
#       The pooled values are recomputed with the frozen
#       cluster_bootstrap_ceiling() and ASSERTED to reproduce the step 18
#       variance_ceiling.csv exactly. Within-form
#       intervals: 2000-draw subject bootstrap stratified by form, percentile 95%.
#   B3  Within-form split-half. Strict grid (step 05, axis i_exp,
#       k 2-7, r 20) and expanded grid (step 11, axis i_exp_ext, k 2-13,
#       r 20), all five constructs. run_cell() is re-run with the original cell
#       seeds and reliability_point() with the original bootstrap seeds; the
#       pooled rows are ASSERTED to reproduce reliability_by_cell_construct.csv
#       (step 05) and reliability_extended.csv (step 11) exactly. Within-form
#       statistic per split:
#       n-weighted mean of per-form Spearman correlations over forms with at
#       least 10 complete pairs (weights = complete pairs), averaged over
#       splits, Spearman-Brown corrected. Interval: 2000-draw subject bootstrap
#       stratified by form, splits held fixed, percentile 95%.
#
# Frozen engines are SOURCED, never redefined: R/cpc18_features.R,
# R/step2_reliability.R, steps/18_variance_ceiling.R (with
# E0_LIB_ONLY <- TRUE), and, for the expanded pool's features only,
# R/lotshape_reconstruction.R (the same reconstruction step 11 uses).
#
# Pre-declared interpretation (fixed before results existed):
# a within-form estimate is the reliability of one fixed instrument. If a
# construct's s2_true_within interval includes zero, the claim that it carries
# trustworthy between-person variance is not supported. A within-form
# split-half value that differs from the pooled one by more than the pooled
# interval's half-width (Spearman-Brown scale) is reported as material. This
# script reports; it decides and applies nothing.
#
# Seeds: base 20261003L. Offsets (documented in audit_config.json):
#   +1          B2 stratified bootstrap
#   +10000 + j  B3 strict-grid within-form bootstrap, task j (task order below)
#   +20000 + j  B3 expanded-grid within-form bootstrap, task j
# The E0 per-subject recomputation reuses the E0 inner seeds from step 18
# (733000 + subject index); the pooled B3 re-runs reuse the original cell and
# bootstrap seeds from steps 05 and 11. Nothing else is stochastic.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(parallel)
})
E0_LIB_ONLY <- TRUE
source(here("R", "cpc18_features.R"))
source(here("R", "step2_reliability.R"))
source(here("R", "lotshape_reconstruction.R"))
source(here("steps", "18_variance_ceiling.R"))  # helpers only

out_dir <- here("output", "elevation", "results", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(df, name) {
  readr::write_csv(df, file.path(out_dir, paste0(name, ".csv")), na = "NA")
  invisible(df)
}
say <- function(...) cat(sprintf(...), sep = "")
mins_since <- function(t0) round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
t_start <- Sys.time()

# ---- constants ---------------------------------------------------------------
SEED_BASE <- 20261003L
SEED_B2 <- SEED_BASE + 1L
SEED_B3_STRICT <- SEED_BASE + 10000L
SEED_B3_EXT <- SEED_BASE + 20000L
N_BOOT_AUDIT <- 2000L
MIN_PAIRS_FORM <- 10L
B_INNER <- 300L                      # frozen E0 value
E0_SEED_INNER_BASE <- 733000L        # frozen E0 value
E0_SEED_OUTER <- 734001L             # frozen E0 value
E0_N_OUT <- 2000L                    # frozen E0 value
N_SPLITS <- 200L                     # frozen value (steps 05 and 11)
N_BOOT <- 2000L                      # frozen value (steps 05 and 11)
K_EVAL <- 14                         # "implied reliability at 14 problems"
RAW_MD5 <- "38a26407a36c690be969f9d9d3fe3e5a"
# Workers do not affect any result (every stochastic task sets its own seed);
# AUDIT_WORKERS lets the audit run beside another job without oversubscribing.
N_WORKERS <- if (nzchar(Sys.getenv("AUDIT_WORKERS"))) as.integer(Sys.getenv("AUDIT_WORKERS")) else
  max(1L, min(12L, parallel::detectCores() - 2L))

# exact numeric reader for saved result CSVs (utils::read.csv round-trips the
# shortest-representation doubles readr writes; the readr/vroom float-parse
# caveat found at E3b is avoided)
read_exact <- function(path) utils::read.csv(path, check.names = FALSE,
                                             stringsAsFactors = FALSE)

# ---- data (frozen preparation, identical to steps 05 / 18) -------------------
raw_path <- here("data", "raw", "All estimation raw data.csv")
stopifnot(identical(unname(tools::md5sum(raw_path)), RAW_MD5))
raw <- load_cpc18(raw_path)
stopifnot(dplyr::n_distinct(raw$SubjID) == 686L,
          dplyr::n_distinct(raw$GameID) == 210L)

# Form = raw Set column, constant within participant (asserted).
form_tbl <- raw |>
  dplyr::distinct(SubjID, Set) |>
  dplyr::arrange(SubjID)
stopifnot(!anyDuplicated(form_tbl$SubjID),             # one Set per subject
          nrow(form_tbl) == 686L,
          all(form_tbl$Set %in% 1:7),
          identical(sort(unique(as.integer(form_tbl$Set))), 1:7))
form_of <- setNames(as.integer(form_tbl$Set), as.character(form_tbl$SubjID))
say("form assertion PASS: Set constant within each of %d participants; forms 1-7 present\n",
    nrow(form_tbl))

features <- add_scaled_features(clean_problem_features(raw))
trials <- raw |>
  dplyr::filter(is_clean_row(raw)) |>
  dplyr::select(SubjID, GameID, Trial, B) |>
  dplyr::left_join(features |>
                     dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev,
                                   loss_informative, has_loss_b, ev_a, ev_b),
                   by = "GameID")
baselines <- fit_loss_baselines(trials)
subj_data <- build_subject_data(trials, baselines)
stopifnot(length(subj_data) == 686L)
k_subj <- vapply(subj_data, function(sd) nrow(sd$feat), integer(1))
k_ref_frozen <- stats::median(k_subj)
stopifnot(k_ref_frozen == 14)
k_harm <- length(k_subj) / sum(1 / k_subj)
say("clean pool: %d subjects; k_ref (median) = %g; harmonic mean of problem counts = %.6f\n",
    length(subj_data), k_ref_frozen, k_harm)

form_sizes <- tibble::tibble(form = 1:7,
                             n_subjects = as.integer(table(factor(form_of, levels = 1:7))))
save_csv(form_sizes, "form_sizes")

# ---- cluster -------------------------------------------------------------------
cl <- makeCluster(N_WORKERS)
invisible(clusterCall(cl, function(lp, root) {
  .libPaths(lp)
  E0_LIB_ONLY <<- TRUE
  suppressPackageStartupMessages(library(tidyverse))
  source(file.path(root, "R", "step2_reliability.R"))
  source(file.path(root, "steps", "18_variance_ceiling.R"))
  NULL
}, .libPaths(), here::here()))
say("workers: %d\n", N_WORKERS)

# =============================================================================
# B1 / B2 input: E0 per-subject full-data points + sampling variances
# =============================================================================
t0 <- Sys.time()
ids <- names(subj_data)
seeds_inner <- E0_SEED_INNER_BASE + seq_along(ids)
committed_seeds <- read_exact(here("output", "elevation", "results", "e0", "seeds_inner.csv"))
stopifnot(identical(as.character(committed_seeds$SubjID), ids),
          identical(as.integer(committed_seeds$inner_seed), as.integer(seeds_inner)))
clusterExport(cl, c("subj_data", "seeds_inner", "B_INNER"), envir = environment())
per_subj <- parLapplyLB(cl, seq_along(ids), function(i) {
  compute_subject_point_var(subj_data[[i]], B_INNER, seeds_inner[i])
})
names(per_subj) <- ids
points_mat <- t(vapply(per_subj, function(x) x$point[E0_CONSTRUCTS],
                       numeric(length(E0_CONSTRUCTS))))
vars_mat <- t(vapply(per_subj, function(x) x$samp_var[E0_CONSTRUCTS],
                     numeric(length(E0_CONSTRUCTS))))
colnames(points_mat) <- colnames(vars_mat) <- E0_CONSTRUCTS
say("E0 per-subject recomputation: %.2f min\n", mins_since(t0))

# Reproduction check against the step 18 per-subject file (exact bits).
pse <- read_exact(here("output", "elevation", "results", "e0", "per_subject_estimates.csv"))
e0_repro <- purrr::map_dfr(E0_CONSTRUCTS, function(cn) {
  d <- pse[pse$construct == cn, ]
  stopifnot(identical(as.character(d$SubjID), ids))
  tibble::tibble(
    check = paste0("per_subject_estimates.csv point/samp_var, ", cn),
    n = nrow(d),
    point_identical = identical(unname(points_mat[, cn]), d$point),
    samp_var_identical = identical(unname(vars_mat[, cn]), d$samp_var),
    point_max_abs_diff = max(abs(points_mat[, cn] - d$point), na.rm = TRUE),
    samp_var_max_abs_diff = max(abs(vars_mat[, cn] - d$samp_var), na.rm = TRUE))
})
save_csv(e0_repro, "e0_per_subject_reproduction")
say("E0 per-subject reproduction (point and samp_var bit-identical): %d of %d constructs\n",
    sum(e0_repro$point_identical & e0_repro$samp_var_identical), nrow(e0_repro))

# =============================================================================
# B1 between-form share (one-way ANOVA by form; KW supplement)
# =============================================================================
form_vec <- unname(form_of[ids])
b1_rows <- list(); b1_means <- list()
for (cn in E0_CONSTRUCTS) {
  ok <- is.finite(points_mat[, cn])
  y <- unname(points_mat[ok, cn]); f <- factor(form_vec[ok], levels = 1:7)
  fit <- stats::lm(y ~ f)
  av <- stats::anova(fit)
  ss_between <- av[["Sum Sq"]][1]; ss_within <- av[["Sum Sq"]][2]
  ss_total <- sum((y - mean(y))^2)
  stopifnot(isTRUE(all.equal(ss_between + ss_within, ss_total)))
  kw <- stats::kruskal.test(y, f)
  b1_rows[[cn]] <- tibble::tibble(
    construct = cn, n_subjects = length(y),
    F = av[["F value"]][1], df_between = av[["Df"]][1], df_within = av[["Df"]][2],
    p = av[["Pr(>F)"]][1],
    ss_between = ss_between, ss_within = ss_within, ss_total = ss_total,
    between_form_share = ss_between / ss_total,
    grand_mean = mean(y),
    kw_chisq_supplement = unname(kw$statistic), kw_p_supplement = kw$p.value)
  b1_means[[cn]] <- tibble::tibble(construct = cn, form = 1:7) |>
    dplyr::left_join(tibble::tibble(form = as.integer(as.character(f)), y = y) |>
                       dplyr::group_by(form) |>
                       dplyr::summarise(n = dplyr::n(), mean = mean(y),
                                        median = stats::median(y), sd = stats::sd(y),
                                        .groups = "drop"),
                     by = "form")
}
b1 <- dplyr::bind_rows(b1_rows)
save_csv(b1, "b1_between_form_anova")
save_csv(dplyr::bind_rows(b1_means), "b1_form_means")
say("\n== B1 between-form share ==\n"); print(as.data.frame(b1[, c("construct", "n_subjects", "F", "p", "between_form_share", "kw_p_supplement")]), digits = 5)

# =============================================================================
# B2 within-form decomposition (ev_sens_agn, loss_avoid_adj)
# =============================================================================
# Pooled: frozen cluster_bootstrap_ceiling(); must reproduce variance_ceiling.csv
t0 <- Sys.time()
vc_pooled <- cluster_bootstrap_ceiling(points_mat, vars_mat, k_ref_frozen,
                                       E0_N_OUT, E0_SEED_OUTER)
vc_committed <- read_exact(here("output", "elevation", "results", "e0", "variance_ceiling.csv"))
num_cols <- setdiff(names(vc_committed), c("construct", "optimization_target", "target_flag"))
stopifnot(identical(vc_pooled$construct, vc_committed$construct))
vc_repro <- purrr::map_dfr(num_cols, function(nm) {
  a <- as.numeric(vc_pooled[[nm]]); b <- as.numeric(vc_committed[[nm]])
  tibble::tibble(check = paste0("variance_ceiling.csv column ", nm),
                 n = length(a), identical = identical(a, b),
                 max_abs_diff = suppressWarnings(max(abs(a - b)[is.finite(a - b)], 0)))
})
vc_repro <- dplyr::bind_rows(vc_repro, tibble::tibble(
  check = "variance_ceiling.csv column optimization_target", n = nrow(vc_committed),
  identical = identical(as.logical(vc_pooled$optimization_target),
                        as.logical(vc_committed$optimization_target)),
  max_abs_diff = 0))
save_csv(vc_repro, "b2_pooled_reproduction")
stopifnot(all(vc_repro$identical))   # gate: pooled values must reproduce exactly
say("pooled variance_ceiling.csv reproduced exactly (%d columns x 5 constructs; %.2f min)\n",
    nrow(vc_repro), mins_since(t0))

B2_CONSTRUCTS <- c("ev_sens_agn", "loss_avoid_adj")
# s2_within with df n - G, G = number of forms present among finite subjects
within_decomp <- function(p, v, f, kk) {
  ok <- is.finite(p) & is.finite(v)
  p <- p[ok]; v <- v[ok]; f <- f[ok]; kk <- kk[ok]
  n <- length(p); G <- length(unique(f))
  dev <- p - stats::ave(p, f)
  s2_within <- sum(dev^2) / (n - G)
  s2_obs <- stats::var(p)
  s2_err <- mean(v)
  s2_diff_w <- s2_within - s2_err
  s2_true_w <- max(0, s2_diff_w)
  s2_true_p <- max(0, s2_obs - s2_err)
  kh <- length(kk) / sum(1 / kk)
  rel_at <- function(s2t, kref) if (s2t <= 0) 0 else s2t / (s2t + s2_err * (kref / K_EVAL))
  k07 <- function(s2t, kref) if (s2t <= 0) Inf else s2_err * kref / (s2t * (1 / 0.7 - 1))
  c(n = n, n_forms = G, s2_obs = s2_obs, s2_within = s2_within, s2_err = s2_err,
    s2_diff_within = s2_diff_w, s2_true_within = s2_true_w,
    ceiling_within = if (s2_within > 0) s2_true_w / s2_within else NA_real_,
    k_harm = kh,
    rel14_within_kref14 = rel_at(s2_true_w, k_ref_frozen),
    rel14_within_kharm = rel_at(s2_true_w, kh),
    k07_within_kref14 = k07(s2_true_w, k_ref_frozen),
    k07_within_kharm = k07(s2_true_w, kh),
    s2_true_pooled = s2_true_p,
    rel14_pooled_kref14 = rel_at(s2_true_p, k_ref_frozen),
    rel14_pooled_kharm = rel_at(s2_true_p, kh),
    k07_pooled_kref14 = k07(s2_true_p, k_ref_frozen),
    k07_pooled_kharm = k07(s2_true_p, kh))
}

t0 <- Sys.time()
strata <- split(seq_along(ids), form_vec)
b2_point <- lapply(B2_CONSTRUCTS, function(cn)
  within_decomp(points_mat[, cn], vars_mat[, cn], form_vec, k_subj))
names(b2_point) <- B2_CONSTRUCTS
b2_boot <- array(NA_real_, c(N_BOOT_AUDIT, length(B2_CONSTRUCTS), length(b2_point[[1]])),
                 dimnames = list(NULL, B2_CONSTRUCTS, names(b2_point[[1]])))
set.seed(SEED_B2)
for (b in seq_len(N_BOOT_AUDIT)) {
  idx <- unlist(lapply(strata, function(s) s[sample.int(length(s), replace = TRUE)]),
                use.names = FALSE)
  for (cn in B2_CONSTRUCTS) {
    b2_boot[b, cn, ] <- within_decomp(points_mat[idx, cn], vars_mat[idx, cn],
                                      form_vec[idx], k_subj[idx])
  }
}
qci <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
b2 <- purrr::map_dfr(B2_CONSTRUCTS, function(cn) {
  pt <- b2_point[[cn]]
  vc <- vc_pooled[vc_pooled$construct == cn, ]
  # the pooled point branch of within_decomp must equal the frozen decomposition
  stopifnot(identical(unname(pt[["s2_obs"]]), vc$sigma2_obs),
            identical(unname(pt[["s2_err"]]), vc$sigma2_err),
            identical(unname(pt[["s2_true_pooled"]]), vc$sigma2_true),
            isTRUE(all.equal(unname(pt[["k07_pooled_kref14"]]), vc$k_0.7, tolerance = 1e-12)))
  ci <- function(key) qci(b2_boot[, cn, key])
  out <- tibble::tibble(
    construct = cn, n_subjects = pt[["n"]], n_forms = pt[["n_forms"]],
    k_ref_median = k_ref_frozen, k_ref_harmonic = pt[["k_harm"]],
    pooled_sigma2_obs = vc$sigma2_obs, pooled_sigma2_err = vc$sigma2_err,
    pooled_sigma2_true = vc$sigma2_true,
    pooled_sigma2_true_lo = vc$sigma2_true_lo, pooled_sigma2_true_hi = vc$sigma2_true_hi,
    pooled_ceiling = vc$ceiling_fulldata,
    pooled_k07_kref14 = vc$k_0.7, pooled_k07_kref14_lo = vc$k_0.7_lo,
    pooled_k07_kref14_hi = vc$k_0.7_hi,
    pooled_rel14_kref14 = pt[["rel14_pooled_kref14"]],
    pooled_rel14_kharm = pt[["rel14_pooled_kharm"]],
    pooled_k07_kharm = pt[["k07_pooled_kharm"]],
    sigma2_within = pt[["s2_within"]],
    sigma2_within_lo = ci("s2_within")[1], sigma2_within_hi = ci("s2_within")[2],
    sigma2_err = pt[["s2_err"]],
    sigma2_err_lo = ci("s2_err")[1], sigma2_err_hi = ci("s2_err")[2],
    sigma2_diff_within = pt[["s2_diff_within"]],
    sigma2_diff_within_lo = ci("s2_diff_within")[1],
    sigma2_diff_within_hi = ci("s2_diff_within")[2],
    sigma2_true_within = pt[["s2_true_within"]],
    sigma2_true_within_lo = ci("s2_true_within")[1],
    sigma2_true_within_hi = ci("s2_true_within")[2],
    ceiling_within = pt[["ceiling_within"]],
    ceiling_within_lo = ci("ceiling_within")[1], ceiling_within_hi = ci("ceiling_within")[2],
    rel14_within_kref14 = pt[["rel14_within_kref14"]],
    rel14_within_kref14_lo = ci("rel14_within_kref14")[1],
    rel14_within_kref14_hi = ci("rel14_within_kref14")[2],
    rel14_within_kharm = pt[["rel14_within_kharm"]],
    rel14_within_kharm_lo = ci("rel14_within_kharm")[1],
    rel14_within_kharm_hi = ci("rel14_within_kharm")[2],
    k07_within_kref14 = pt[["k07_within_kref14"]],
    k07_within_kref14_lo = ci("k07_within_kref14")[1],
    k07_within_kref14_hi = ci("k07_within_kref14")[2],
    k07_within_kharm = pt[["k07_within_kharm"]],
    k07_within_kharm_lo = ci("k07_within_kharm")[1],
    k07_within_kharm_hi = ci("k07_within_kharm")[2],
    n_boot = N_BOOT_AUDIT, bootstrap = "subject resample stratified by form",
    boot_seed = SEED_B2)
  out$sigma2_true_within_ci_excludes_zero <- out$sigma2_diff_within_lo > 0
  out
})
save_csv(b2, "b2_within_form_decomposition")
say("\n== B2 within-form decomposition (%.2f min) ==\n", mins_since(t0))
print(as.data.frame(b2[, c("construct", "pooled_sigma2_true", "sigma2_within", "sigma2_err",
                           "sigma2_true_within", "sigma2_true_within_lo",
                           "sigma2_true_within_hi", "sigma2_true_within_ci_excludes_zero")]),
      digits = 5)

# =============================================================================
# B3 within-form split-half on the strict and expanded grids
# =============================================================================
# Expanded-pool subject data, prepared exactly as in 11_expanded_pool_reliability.R
param_cols <- c("Ha", "pHa", "La", "LotShapeA", "LotNumA",
                "Hb", "pHb", "Lb", "LotShapeB", "LotNumB", "Amb", "Corr")
problems_exp <- raw |>
  dplyr::distinct(GameID, dplyr::across(dplyr::all_of(param_cols))) |>
  dplyr::arrange(GameID) |>
  dplyr::filter(Amb == 0)
feat_exp <- expanded_problem_features(problems_exp,
                                      reconstruct_all_problems(problems_exp))
trials_ext <- raw |>
  dplyr::semi_join(problems_exp, by = "GameID") |>
  dplyr::select(SubjID, GameID, Trial, B) |>
  dplyr::left_join(feat_exp |>
                     dplyr::select(GameID, z_ev, z_sd, z_eloss, d_ev, loss_informative,
                                   has_loss_b, ev_a, ev_b), by = "GameID")
subj_data_ext <- build_subject_data(trials_ext, fit_loss_baselines(trials_ext))
stopifnot(length(subj_data_ext) == 686L)

seeds_strict <- read_exact(here("output", "results", "step2", "seeds_cells.csv"))
boot_strict <- read_exact(here("output", "results", "step2", "seeds_bootstrap.csv"))
seeds_ext <- read_exact(here("output", "results", "step4", "seeds_cells_extended.csv"))
boot_ext <- read_exact(here("output", "results", "step4", "seeds_bootstrap_extended.csv"))
rel_strict <- read_exact(here("output", "results", "step2", "reliability_by_cell_construct.csv"))
rel_ext <- read_exact(here("output", "results", "step4", "reliability_extended.csv"))

grid_specs <- c(
  lapply(2:7, function(k) {
    id <- sprintf("i_exp_k%d_r20", k)
    list(grid = "strict", data = "clean", axis = "i_exp", regime = "exp", k = k,
         r = 20, lambda = NA, cell_id = id,
         seed = seeds_strict$seed[seeds_strict$cell_id == id])
  }),
  lapply(2:13, function(k) {
    id <- sprintf("ext_exp_k%d_r20", k)
    list(grid = "expanded", data = "expanded", axis = "i_exp_ext", regime = "exp", k = k,
         r = 20, lambda = NA, cell_id = id,
         seed = seeds_ext$seed[seeds_ext$cell_id == id])
  }))
stopifnot(all(vapply(grid_specs, function(s) length(s$seed) == 1L, logical(1))))

t0 <- Sys.time()
clusterExport(cl, c("subj_data_ext", "N_SPLITS"), envir = environment())
cells <- parLapplyLB(cl, grid_specs, function(spec) {
  sdat <- if (spec$data == "clean") subj_data else subj_data_ext
  run_cell(spec, sdat, N_SPLITS, seed = spec$seed)
})
names(cells) <- vapply(cells, function(x) x$spec$cell_id, character(1))
say("\nB3 run_cell re-run (%d cells): %.2f min\n", length(cells), mins_since(t0))

# tasks in a fixed, documented order: strict cells then expanded cells, each
# cell's constructs in step2_constructs order
tasks <- purrr::list_flatten(lapply(cells, function(cell) {
  lapply(step2_constructs, function(con) {
    sp <- cell$spec
    bs <- if (sp$grid == "strict") boot_strict else boot_ext
    list(grid = sp$grid, cell_id = sp$cell_id, k = sp$k, construct = con,
         E1 = cell$est[, , 1, con], E2 = cell$est[, , 2, con],
         forms = unname(form_of[dimnames(cell$est)[[1]]]),
         boot_seed_committed = bs$seed[bs$cell_id == sp$cell_id & bs$construct == con])
  })
}))
j_strict <- 0L; j_ext <- 0L
for (j in seq_along(tasks)) {
  if (tasks[[j]]$grid == "strict") {
    j_strict <- j_strict + 1L; tasks[[j]]$audit_seed <- SEED_B3_STRICT + j_strict
  } else {
    j_ext <- j_ext + 1L; tasks[[j]]$audit_seed <- SEED_B3_EXT + j_ext
  }
  stopifnot(length(tasks[[j]]$boot_seed_committed) == 1L, !anyNA(tasks[[j]]$forms))
}

# New audit statistic (not a frozen engine): per split, n-weighted mean over
# forms of the frozen split_correlations() value, forms with >= 10 complete
# pairs only (split_correlations returns NA below min_pairs), weights = the
# form's complete pairs in that split; averaged over usable splits.
within_form_mean_r <- function(E1, E2, forms, min_pairs) {
  fl <- sort(unique(forms))
  num <- numeric(ncol(E1)); den <- numeric(ncol(E1)); nf <- integer(ncol(E1))
  for (g in fl) {
    rows <- which(forms == g)
    e1 <- E1[rows, , drop = FALSE]; e2 <- E2[rows, , drop = FALSE]
    rs <- split_correlations(e1, e2, min_pairs)
    np <- colSums(is.finite(e1) & is.finite(e2))
    use <- !is.na(rs)
    num[use] <- num[use] + np[use] * rs[use]
    den[use] <- den[use] + np[use]
    nf[use] <- nf[use] + 1L
  }
  per_split <- ifelse(den > 0, num / den, NA_real_)
  c(r = mean(per_split, na.rm = TRUE), n_splits_used = sum(!is.na(per_split)),
    mean_forms_used = mean(nf))
}
clusterExport(cl, c("within_form_mean_r", "N_BOOT", "N_BOOT_AUDIT",
                    "MIN_PAIRS_FORM"), envir = environment())

t0 <- Sys.time()
b3_list <- parLapplyLB(cl, tasks, function(t) {
  ident <- any(is.finite(t$E1) & is.finite(t$E2))
  if (!ident) {
    pooled <- tibble::tibble(
      r_raw = NA_real_, ci_lo_raw = NA_real_, ci_hi_raw = NA_real_,
      r_sb = NA_real_, ci_lo_sb = NA_real_, ci_hi_sb = NA_real_,
      n_subjects = nrow(t$E1), n_splits_used = 0L, n_splits = ncol(t$E1),
      finite_frac = 0, n_boot = N_BOOT)
    wf <- c(r = NA_real_, n_splits_used = 0, mean_forms_used = 0)
    wf_ci <- c(NA_real_, NA_real_)
  } else {
    pooled <- reliability_point(t$E1, t$E2, N_BOOT, seed = t$boot_seed_committed)
    wf <- within_form_mean_r(t$E1, t$E2, t$forms, MIN_PAIRS_FORM)
    strata <- split(seq_len(nrow(t$E1)), t$forms)
    set.seed(t$audit_seed)
    boots <- replicate(N_BOOT_AUDIT, {
      idx <- unlist(lapply(strata, function(s) s[sample.int(length(s), replace = TRUE)]),
                    use.names = FALSE)
      within_form_mean_r(t$E1[idx, , drop = FALSE], t$E2[idx, , drop = FALSE],
                         t$forms[idx], MIN_PAIRS_FORM)[["r"]]
    })
    wf_ci <- stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  }
  fn <- table(factor(t$forms, levels = 1:7))
  tibble::tibble(grid = t$grid, cell_id = t$cell_id, k = t$k, construct = t$construct) |>
    dplyr::bind_cols(pooled) |>
    dplyr::bind_cols(tibble::tibble(
      within_r_raw = wf[["r"]], within_ci_lo_raw = wf_ci[1], within_ci_hi_raw = wf_ci[2],
      within_r_sb = sb_correct(wf[["r"]]), within_ci_lo_sb = sb_correct(wf_ci[1]),
      within_ci_hi_sb = sb_correct(wf_ci[2]),
      within_n_splits_used = wf[["n_splits_used"]],
      within_mean_forms_used = wf[["mean_forms_used"]],
      eligible_by_form = paste(sprintf("%d:%d", 1:7, as.integer(fn)), collapse = " "),
      boot_seed_committed = t$boot_seed_committed, audit_seed = t$audit_seed,
      n_boot_audit = N_BOOT_AUDIT))
})
b3 <- dplyr::bind_rows(b3_list)
say("B3 pooled + within-form bootstraps (%d tasks): %.2f min\n", length(tasks), mins_since(t0))

# Pooled reproduction: every recomputed pooled row must equal the saved row.
repro_cols <- c("r_raw", "ci_lo_raw", "ci_hi_raw", "r_sb", "ci_lo_sb", "ci_hi_sb",
                "n_subjects", "n_splits_used", "n_splits", "finite_frac", "n_boot")
b3_repro <- purrr::map_dfr(seq_len(nrow(b3)), function(i) {
  row <- b3[i, ]
  ref <- if (row$grid == "strict") rel_strict else rel_ext
  cr <- ref[ref$cell_id == row$cell_id & ref$construct == row$construct, ]
  stopifnot(nrow(cr) == 1L)
  same <- vapply(repro_cols, function(nm) {
    a <- as.numeric(row[[nm]]); b <- as.numeric(cr[[nm]])
    identical(a, b) || (is.na(a) && is.na(b))
  }, logical(1))
  tibble::tibble(grid = row$grid, cell_id = row$cell_id, construct = row$construct,
                 committed_file = if (row$grid == "strict")
                   "results/step2/reliability_by_cell_construct.csv" else
                   "results/step4/reliability_extended.csv",
                 n_columns_checked = length(repro_cols), all_identical = all(same),
                 differing_columns = paste(repro_cols[!same], collapse = ";"))
})
save_csv(b3_repro, "b3_pooled_reproduction")
stopifnot(all(b3_repro$all_identical))
say("B3 pooled reproduction: %d of %d cell-construct rows bit-identical to the committed files\n",
    sum(b3_repro$all_identical), nrow(b3_repro))

b3 <- b3 |>
  dplyr::mutate(
    pooled_half_width_sb = (ci_hi_sb - ci_lo_sb) / 2,
    diff_within_minus_pooled_sb = within_r_sb - r_sb,
    material = abs(diff_within_minus_pooled_sb) > pooled_half_width_sb,
    within_ci_excludes_pooled_point = within_ci_lo_sb > r_sb | within_ci_hi_sb < r_sb)
save_csv(b3, "b3_within_form_splithalf")
say("\n== B3 (SB scale) ==\n")
print(as.data.frame(b3[, c("grid", "k", "construct", "r_sb", "ci_lo_sb", "ci_hi_sb",
                           "within_r_sb", "within_ci_lo_sb", "within_ci_hi_sb", "material")]),
      digits = 4)

stopCluster(cl)

# ---- reproduction summary + config ---------------------------------------------
save_csv(dplyr::bind_rows(
  e0_repro |> dplyr::transmute(check, n, identical = point_identical & samp_var_identical,
                               max_abs_diff = pmax(point_max_abs_diff, samp_var_max_abs_diff)),
  vc_repro,
  b3_repro |> dplyr::transmute(check = paste(committed_file, cell_id, construct),
                               n = n_columns_checked, identical = all_identical,
                               max_abs_diff = NA_real_)),
  "reproduction_checks")

in_files <- c(
  raw = raw_path,
  e0_per_subject_estimates = here("output", "elevation", "results", "e0", "per_subject_estimates.csv"),
  e0_seeds_inner = here("output", "elevation", "results", "e0", "seeds_inner.csv"),
  e0_variance_ceiling = here("output", "elevation", "results", "e0", "variance_ceiling.csv"),
  step2_seeds_cells = here("output", "results", "step2", "seeds_cells.csv"),
  step2_seeds_bootstrap = here("output", "results", "step2", "seeds_bootstrap.csv"),
  step2_reliability = here("output", "results", "step2", "reliability_by_cell_construct.csv"),
  step4_seeds_cells = here("output", "results", "step4", "seeds_cells_extended.csv"),
  step4_seeds_bootstrap = here("output", "results", "step4", "seeds_bootstrap_extended.csv"),
  step4_reliability = here("output", "results", "step4", "reliability_extended.csv"))
jsonlite::write_json(list(
  seed_base = SEED_BASE,
  seed_offsets = list(b2_stratified_bootstrap = 1L,
                      b3_strict_within_form_bootstrap = "10000 + task index j (strict tasks in cell order k2..k7, constructs in step2_constructs order)",
                      b3_expanded_within_form_bootstrap = "20000 + task index j (expanded tasks in cell order k2..k13, constructs in step2_constructs order)"),
  reused_committed_seeds = list(e0_inner = "733000 + subject index (elevation/results/e0/seeds_inner.csv)",
                                e0_outer = E0_SEED_OUTER,
                                b3_cells = "results/step2/seeds_cells.csv, results/step4/seeds_cells_extended.csv",
                                b3_pooled_bootstrap = "results/step2/seeds_bootstrap.csv, results/step4/seeds_bootstrap_extended.csv"),
  constants = list(B_inner = B_INNER, e0_n_out = E0_N_OUT, n_splits = N_SPLITS,
                   n_boot_pooled = N_BOOT, n_boot_audit = N_BOOT_AUDIT,
                   min_pairs_per_form = MIN_PAIRS_FORM, k_eval = K_EVAL,
                   k_ref_median = k_ref_frozen, k_ref_harmonic = k_harm,
                   workers = N_WORKERS),
  form_definition = "raw column Set (1-7), asserted constant within participant",
  frozen_engines_sourced = c("R/cpc18_features.R", "R/step2_reliability.R",
                             "R/lotshape_reconstruction.R",
                             "steps/18_variance_ceiling.R (E0_LIB_ONLY)"),
  input_md5 = as.list(vapply(in_files, function(f) unname(tools::md5sum(f)), character(1))),
  wall_minutes = mins_since(t_start),
  R_version = R.version.string),
  file.path(out_dir, "audit_config.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

say("\naudit_within_form complete in %.2f min. Outputs in elevation/results/audit/\n",
    mins_since(t_start))

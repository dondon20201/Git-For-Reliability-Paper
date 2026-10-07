# Step 16: reliability on a fixed participant cohort
# Reads:   data/raw/All estimation raw data.csv; output/results/step2/, step4/ (steps 05, 09-11)
# Writes:  output/results/step2/fixed_cohort_reliability.csv, output/results/step2/fixed_cohort_membership.csv
# Paper:   Supplement (Fixed-cohort robustness)
# Runtime: about 22 minutes with 12 parallel workers

# Fixed-cohort robustness of the reliability versus problem-count curves.
#
# REAL DATA. The original curves (steps 05 and 11) let the SUBJECT set shrink
# as k grows (run_cell keeps subjects with >= 2k problems), so the count
# effect could partly reflect who is measured rather than budget. This script
# re-runs the SAME engine (run_cell / reliability_point — no
# re-implementation) with the SAME saved cell and bootstrap seeds on a
# CONSTANT participant cohort per axis: subjects eligible at the largest
# budget of that axis, held fixed at every k. Generalizes step 07 (whose
# fixed_pool_check.csv is untouched) to all five constructs, both pools, and
# the description regime.
#
# Reproduction anchor (pre-declared): at each axis's k_max the fixed cohort
# EQUALS the original eligible set, so with the saved seeds the fixed-cohort
# r_sb must equal the original r_sb to < 1e-9 for every construct. In-script
# stopifnot; halt on any miss.
#
# Only the two fixed_cohort_* CSVs are written; the earlier reliability
# tables, trials-needed tables, step 07 outputs and all figures are untouched.

library(tidyverse)
library(here)
library(parallel)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "fetch_official_cpc18.R"))
source(here("R", "lotshape_reconstruction.R"))
source(here("R", "step2_reliability.R"))

res2_dir <- here("output", "results", "step2")
res4_dir <- here("output", "results", "step4")

N_SPLITS <- 200   # pre-declared
N_BOOT <- 2000    # pre-declared

# ---- saved reference files --------------------------------------------------
committed_strict <- read_csv(file.path(res2_dir, "reliability_by_cell_construct.csv"),
                             show_col_types = FALSE)
committed_ext <- read_csv(file.path(res4_dir, "reliability_extended.csv"),
                          show_col_types = FALSE)
seeds_strict <- read_csv(file.path(res2_dir, "seeds_cells.csv"),
                         show_col_types = FALSE)
seeds_ext <- read_csv(file.path(res4_dir, "seeds_cells_extended.csv"),
                      show_col_types = FALSE)
boot_seeds <- bind_rows(
  read_csv(file.path(res2_dir, "seeds_bootstrap.csv"), show_col_types = FALSE),
  read_csv(file.path(res4_dir, "seeds_bootstrap_extended.csv"),
           show_col_types = FALSE)
)

# ---- raw data (md5-verified fetch, as in every step) -------------------------
fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))
raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))

# ---- STRICT pool subject data (exact step 05 path) ---------------------------
features <- add_scaled_features(clean_problem_features(raw))
trials_strict <- raw |>
  filter(is_clean_row(raw)) |>
  select(SubjID, GameID, Trial, B) |>
  left_join(features |>
              select(GameID, z_ev, z_sd, z_eloss, d_ev, loss_informative,
                     has_loss_b, ev_a, ev_b), by = "GameID")
baselines_strict <- fit_loss_baselines(trials_strict)
subj_data_strict <- build_subject_data(trials_strict, baselines_strict)
cat(sprintf("strict pool: %d subjects\n", length(subj_data_strict)))

# ---- EXPANDED pool subject data (exact step 11 path, gates included) ---------
gate <- read_csv(file.path(res4_dir, "reconstruction_validation.csv"),
                 show_col_types = FALSE)
stopifnot(all(gate$pass))
verdict_a4 <- read_csv(file.path(res4_dir, "expanded_recovery_verdict.csv"),
                       show_col_types = FALSE)
stopifnot(all(verdict_a4$pass))

param_cols <- c("Ha", "pHa", "La", "LotShapeA", "LotNumA",
                "Hb", "pHb", "Lb", "LotShapeB", "LotNumB", "Amb", "Corr")
problems_exp <- raw |>
  distinct(GameID, across(all_of(param_cols))) |>
  arrange(GameID) |>
  filter(Amb == 0)
feat_exp <- expanded_problem_features(problems_exp,
                                      reconstruct_all_problems(problems_exp))
trials_ext <- raw |>
  semi_join(problems_exp, by = "GameID") |>
  select(SubjID, GameID, Trial, B) |>
  left_join(feat_exp |>
              select(GameID, z_ev, z_sd, z_eloss, d_ev, loss_informative,
                     has_loss_b, ev_a, ev_b), by = "GameID")
baselines_ext <- fit_loss_baselines(trials_ext)
subj_data_ext <- build_subject_data(trials_ext, baselines_ext)
cat(sprintf("expanded pool: %d subjects\n", length(subj_data_ext)))

# ---- FIXED COHORTS: eligible at the axis's largest budget, constant across k.
# Subsetting by name preserves the original subj_data order, so at k_max the
# eligible vector inside run_cell is identical to the original run's and the
# saved seeds reproduce it draw for draw. ---------------------------------------
cohort_strict <- names(subj_data_strict)[
  vapply(subj_data_strict, function(s) nrow(s$feat), integer(1)) >= 14]
cohort_ext <- names(subj_data_ext)[
  vapply(subj_data_ext, function(s) nrow(s$feat), integer(1)) >= 26]
stopifnot(length(cohort_strict) > 0, length(cohort_ext) > 0)
sd_fixed_strict <- subj_data_strict[cohort_strict]
sd_fixed_ext <- subj_data_ext[cohort_ext]
cat(sprintf("fixed cohorts: strict >=14 clean -> %d subjects; expanded >=26 -> %d subjects\n",
            length(cohort_strict), length(cohort_ext)))

# Stable hash of the sorted SubjID vector (md5 of the newline-joined ids,
# written in binary mode so the bytes are platform-independent).
cohort_hash <- function(ids) {
  f <- tempfile()
  con <- file(f, "wb")
  writeLines(sort(ids), con, sep = "\n")
  close(con)
  h <- unname(tools::md5sum(f))
  unlink(f)
  h
}
hash_strict <- cohort_hash(cohort_strict)
hash_ext <- cohort_hash(cohort_ext)

# ---- cell specs: all reported problem-count cells, seeds bound from the
# saved seeds CSVs by cell_id (never re-derived). Strict iii_lambda and
# ii_reps are out of scope (problem-count axes only). --------------------------
bind_seed <- function(s, seeds_df) {
  hit <- seeds_df$seed[seeds_df$cell_id == s$cell_id]
  stopifnot(length(hit) == 1L)
  s$seed <- hit
  s
}
strict_specs <- c(
  lapply(2:7, function(k) list(axis = "i_exp", regime = "exp", k = k, r = 20,
                               lambda = NA, pool = "strict",
                               cell_id = sprintf("i_exp_k%d_r20", k))),
  lapply(2:7, function(k) list(axis = "i_desc", regime = "desc", k = k, r = 5,
                               lambda = NA, pool = "strict",
                               cell_id = sprintf("i_desc_k%d_r5", k)))
)
strict_specs <- lapply(strict_specs, bind_seed, seeds_df = seeds_strict)

realized <- read_csv(file.path(res4_dir, "realized_grid.csv"),
                     show_col_types = FALSE)
ext_specs <- lapply(realized$k, function(k) {
  list(axis = "i_exp_ext", regime = "exp", k = k, r = 20, lambda = NA,
       pool = "expanded", cell_id = sprintf("ext_exp_k%d_r20", k))
})
ext_specs <- lapply(ext_specs, bind_seed, seeds_df = seeds_ext)

all_specs <- c(strict_specs, ext_specs)
stopifnot(length(all_specs) == 12 + nrow(realized))

# ---- run every cell on its fixed cohort (same engine, saved seeds).
# Bootstrap CIs computed in-worker with the saved (cell_id, construct)
# bootstrap seed, mirroring the NA guard of steps 05 and 11 exactly. -----------
n_workers <- max(1, min(12, detectCores() - 2))
cat("workers:", n_workers, "\n")
cl <- makeCluster(n_workers)
invisible(clusterCall(cl, function(lp, root) {
  .libPaths(lp)
  source(file.path(root, "R", "step2_reliability.R"))
  NULL
}, .libPaths(), here::here()))
clusterExport(cl, c("sd_fixed_strict", "sd_fixed_ext", "N_SPLITS", "N_BOOT",
                    "boot_seeds"), envir = environment())

t0 <- Sys.time()
cell_rows <- parLapplyLB(cl, all_specs, function(spec) {
  sd_pool <- if (spec$pool == "strict") sd_fixed_strict else sd_fixed_ext
  cohort <- names(sd_pool)
  out <- run_cell(spec, sd_pool, N_SPLITS, seed = spec$seed)
  # COHORT CONSTANCY (pre-declared): every fixed-cohort member has
  # >= 2*k_max >= 2k problems, so the eligible set must be the whole cohort.
  stopifnot(out$n_eligible == length(cohort),
            identical(dimnames(out$est)[[1]], cohort))
  rows <- lapply(step2_constructs, function(con) {
    E1 <- out$est[, , 1, con]
    E2 <- out$est[, , 2, con]
    bs <- boot_seeds$seed[boot_seeds$cell_id == spec$cell_id &
                            boot_seeds$construct == con]
    stopifnot(length(bs) == 1L)
    res <- if (!any(is.finite(E1) & is.finite(E2))) {
      tibble::tibble(
        r_raw = NA_real_, ci_lo_raw = NA_real_, ci_hi_raw = NA_real_,
        r_sb = NA_real_, ci_lo_sb = NA_real_, ci_hi_sb = NA_real_,
        n_subjects = nrow(E1), n_splits_used = 0L, n_splits = ncol(E1),
        finite_frac = 0, n_boot = N_BOOT)
    } else {
      reliability_point(E1, E2, N_BOOT, seed = bs)
    }
    dplyr::bind_cols(
      tibble::tibble(pool = spec$pool, cell_id = spec$cell_id,
                     axis = spec$axis, regime = spec$regime, k = spec$k,
                     construct = con), res)
  })
  dplyr::bind_rows(rows)
})
stopCluster(cl)
fixed <- bind_rows(cell_rows)
cat("estimation + bootstrap wall time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
stopifnot(nrow(fixed) == length(all_specs) * length(step2_constructs))

# ---- join the original variable-cohort estimates -----------------------------
committed <- bind_rows(
  committed_strict |> mutate(pool = "strict"),
  committed_ext |> mutate(pool = "expanded")
) |> select(pool, cell_id, construct, r_sb_variable = r_sb)
fixed <- fixed |>
  left_join(committed, by = c("pool", "cell_id", "construct")) |>
  mutate(total_problems = 2L * k,
         d_fixed_minus_variable = r_sb - r_sb_variable)
# every reported cell-construct must have an original counterpart, and an NA
# fixed-cohort estimate must pair with an NA original one (k < 4 cells)
stopifnot(!any(is.na(fixed$r_sb) & is.finite(fixed$r_sb_variable)),
          sum(is.finite(fixed$r_sb)) ==
            sum(is.finite(fixed$r_sb) & is.finite(fixed$r_sb_variable)))

# ---- REPRODUCTION ANCHOR (pre-declared): at each axis's k_max the
# fixed-cohort r_sb equals the original r_sb to < 1e-9 for every construct -----
anchor_cells <- c("i_exp_k7_r20", "i_desc_k7_r5",
                  sprintf("ext_exp_k%d_r20", max(realized$k)))
anc <- fixed |> filter(cell_id %in% anchor_cells)
stopifnot(nrow(anc) == 3 * length(step2_constructs),
          all(is.finite(anc$r_sb)), all(is.finite(anc$r_sb_variable)))
max_anchor_dev <- max(abs(anc$r_sb - anc$r_sb_variable))
cat(sprintf("reproduction anchor: %d cell-constructs, max |dev| = %.3g\n",
            nrow(anc), max_anchor_dev))
if (!(max_anchor_dev < 1e-9)) {
  stop(sprintf("REPRODUCTION ANCHOR FAILED: max |fixed r_sb - committed r_sb| at k_max = %.3g >= 1e-9; the re-run is not faithful — STOP",
               max_anchor_dev))
}

# ---- COHORT CONSTANCY summary (pre-declared) ----------------------------------
constancy <- fixed |>
  group_by(pool, axis) |>
  summarise(n_levels = n_distinct(n_subjects), n = unique(n_subjects),
            .groups = "drop")
print(constancy)
stopifnot(all(constancy$n_levels == 1),
          constancy$n[constancy$pool == "strict"] == length(cohort_strict),
          constancy$n[constancy$pool == "expanded"] == length(cohort_ext))

# ---- GATE G1 (pre-declared interpretation gate) -------------------------------
# loss_avoid_adj on the fixed cohort: monotone-nondecreasing in k on the
# STRICT experience axis up to tolerance 0.02, AND reaches SB 0.5 within the
# strict observed grid.
la <- fixed |>
  filter(axis == "i_exp", construct == "loss_avoid_adj") |>
  arrange(k)
stopifnot(nrow(la) == 6, all(is.finite(la$r_sb)))
steps <- diff(la$r_sb)
g1_mono <- all(steps >= -0.02)
g1_cross <- any(la$r_sb >= 0.5)
cat("\n---- GATE G1 (pre-declared) ----\n")
cat(sprintf("loss_avoid_adj fixed-cohort r_sb (strict i_exp, k=2..7): %s\n",
            paste(sprintf("%.4f", la$r_sb), collapse = ", ")))
cat(sprintf("min consecutive step = %+.4f (bound: >= -0.02) -> monotone %s\n",
            min(steps), ifelse(g1_mono, "OK", "VIOLATED")))
cat(sprintf("SB 0.5 crossing within strict grid: %s (first at total problems = %s)\n",
            ifelse(g1_cross, "YES", "NO"),
            ifelse(g1_cross, as.character(min(la$total_problems[la$r_sb >= 0.5])), "none")))
cat(sprintf("G1 %s\n", ifelse(g1_mono && g1_cross, "PASSES", "FAILS")))

# ---- outputs -----------------------------------------------------------------
out_reliability <- fixed |>
  select(pool, axis, regime, k, total_problems, construct, n_subjects,
         r_raw, r_sb, ci_lo_sb, ci_hi_sb, finite_frac, n_splits_used,
         r_sb_variable, d_fixed_minus_variable)
write_csv(out_reliability, file.path(res2_dir, "fixed_cohort_reliability.csv"))

# Cohort composition by Set (leading SubjID digit, as in step 07;
# verified against the Set column of the raw file).
subj_info <- raw |> distinct(SubjID, Set)
mk_membership <- function(pool, axes, definition, ids, h) {
  tibble(SubjID = as.numeric(ids)) |>
    left_join(subj_info, by = "SubjID") |>
    count(Set, name = "n_in_set") |>
    mutate(pool = pool, axes = axes, cohort_definition = definition,
           n_cohort = length(ids), cohort_hash = h, .before = 1)
}
membership <- bind_rows(
  mk_membership("strict", "i_exp;i_desc",
                ">= 14 clean problems (eligible at k_max = 7)",
                cohort_strict, hash_strict),
  mk_membership("expanded", "i_exp_ext",
                ">= 26 Amb==0 problems (eligible at k_max = 13)",
                cohort_ext, hash_ext)
)
write_csv(membership, file.path(res2_dir, "fixed_cohort_membership.csv"))

cat("\nDone: fixed_cohort_reliability.csv, fixed_cohort_membership.csv\n")

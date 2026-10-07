# Step 36: summary cells for Study 2, the ceiling section and the Supplement
# Reads:   output/elevation/results/dispb/, dispc/, sens/, p3_frey/, e4/ (steps 25, 27, 29, 33, 35)
# Writes:  output/elevation/results/dispb/b_cells.csv, dispc/c_cells.csv, sens/sens_cells.csv,
#          e4/feature_collinearity.csv, p3_frey/provenance.csv; adds summary rows to p3_frey/modelpar_reliability.csv
# Paper:   Study 2; The reliability ceiling reflects true-score variance; Supplement; Abstract
# Runtime: a few seconds

# Surface already-computed values into small "manuscript cell" CSVs.
#
# PURE SURFACING. Nothing here is estimated, re-fit, or re-derived from data:
# every value is COPIED out of an existing gate JSON or results CSV written by
# an earlier step, and the source field is recorded alongside it. The only
# arithmetic performed is (a) a median and (b) a count over the existing
# per-parameter retest values in p3_frey/modelpar_reliability.csv, both
# self-checked below.
#
# Fidelity policy: values are written at the FULL precision of the source
# (shortest exact round-trip representation), never at a rounded precision.
# Each write is self-checked against an expected value at the precision that
# value is quoted. The manuscript renders these through num(..., digits = 2),
# so full precision in the cell file and rounded precision in the prose agree
# by construction.
#
# p3_frey/modelpar_reliability.csv is modified additively (2 summary rows +
# 1 column; existing cells stay byte-identical).
#
# Deterministic - no RNG.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
  library(here)
})

er <- function(...) here("output", "elevation", "results", ...)

# ---- helpers ---------------------------------------------------------------

# Shortest decimal string that reads back as EXACTLY the committed double.
fmt_exact <- function(x) {
  stopifnot(length(x) == 1L, is.numeric(x), is.finite(x))
  if (x == round(x) && abs(x) < 1e15) {
    return(format(x, scientific = FALSE, trim = TRUE))
  }
  for (d in 1:17) {
    s <- formatC(x, format = "g", digits = d)
    if (identical(as.numeric(s), x)) return(s)
  }
  stop("fmt_exact: no exact round-trip for ", x)
}

# Self-check a surfaced value against its expected value, at the precision the
# expected value is quoted. A mismatch is a hard stop, never a nudge.
checked <- function(x, expected_str, what) {
  nd <- nchar(sub("^[^.]*\\.?", "", expected_str))
  got <- formatC(as.numeric(x), format = "f", digits = nd)
  if (!identical(got, formatC(as.numeric(expected_str), format = "f", digits = nd))) {
    stop("cell self-check FAILED for ", what, ": committed source gives ", got,
         " at ", nd, " dp, expected ", expected_str, call. = FALSE)
  }
  fmt_exact(as.numeric(x))
}

jget <- function(obj, path) {
  for (p in path) {
    stopifnot(!is.null(obj[[p]]))
    obj <- obj[[p]]
  }
  stopifnot(length(obj) == 1L)
  obj[[1]]
}

write_cells <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write_csv(df, path, na = "NA")
  # re-read and assert every value survived the round trip verbatim
  rr <- read_csv(path, col_types = cols(.default = col_character()))
  stopifnot(identical(rr$cell, df$cell), identical(rr$value, df$value))
  invisible(path)
}

# ---- 1. e4/feature_collinearity.csv ----------------------------------------
# Source: dispb_gate.json B3.feature_collinearity, real_form1 / loss_avoid_adj
# / K = 12 Spearman between F2 (loss exposure) and F4 (expected information).

b <- fromJSON(er("dispb", "dispb_gate.json"), simplifyVector = FALSE)

fc_rows <- b$B3$feature_collinearity$real_track_primary_target
pick_fc <- function(K, field) {
  hit <- Filter(function(r) r$track == "real_form1" &&
                  r$target == "loss_avoid_adj" && r$K == K, fc_rows)
  stopifnot(length(hit) == 1L)
  hit[[1]][[field]]
}

f2f4_spearman_k12 <- checked(pick_fc(12, "spearman_F2_F4"), "0.90", "F2_F4_real")

write_cells(tibble(
  cell    = "F2_F4_real",
  value   = f2f4_spearman_k12,
  spearman = f2f4_spearman_k12,
  source  = paste0("dispb_gate.json B3.feature_collinearity ",
                   "real_form1 loss_avoid_adj K12 spearman_F2_F4")
), er("e4", "feature_collinearity.csv"))

# ---- 2. dispb/b_cells.csv ---------------------------------------------------

b_spec <- list(
  list("b1_sigma_set_true_cov", "0.2214",
       c("B1", "G-B1b", "mean_sigma_set_true_cov")),
  list("b1_sigma_cov_lo", "0.1922", c("B1", "G-B1b", "interval_90", "lo")),
  list("b1_sigma_cov_hi", "0.2469", c("B1", "G-B1b", "interval_90", "hi")),
  list("b1_n_partitions", "200",    c("B1", "G-B1b", "n_partitions")),
  list("b1_icc_cov_lo", "0.6265",
       c("B1", "G-B1a", "per_estimand", "covariance_subject_split", "min_ICC_set")),
  list("b1_icc_cov_hi", "0.7154",
       c("B1", "G-B1a", "per_estimand", "covariance_subject_split", "max_ICC_set")),
  list("b1_null_perm", "0.0164",
       c("B1", "G-B1c", "permuted_sigma_set_true_cov_mean")),
  list("b2_single_best_insample", "0.2366",
       c("B2", "G-B2a", "in_sample_single_best_gain")),
  list("b2_single_best_oos", "0.0883",
       c("B2", "G-B2a", "out_sample_single_best_gain")),
  list("b2_decile_oos", "0.1534",
       c("B2", "G-B2b", "top_decile_out_of_sample_gain")),
  list("b2_decile_oos_lo", "0.0452", c("B2", "G-B2b", "interval_90", "lo")),
  list("b2_decile_oos_hi", "0.2476", c("B2", "G-B2b", "interval_90", "hi")),
  list("b3_delta_dropF4", "0.0279",
       c("B3", "G-B3b", "max_abs_delta_gate_metric"))
)

b_cells <- bind_rows(lapply(b_spec, function(s) {
  tibble(cell = s[[1]],
         value = checked(jget(b, s[[3]]), s[[2]], s[[1]]),
         source = paste0("dispb_gate.json ", paste(s[[3]], collapse = ".")))
}))

# the two F2/F4 Pearson collinearity cells come from the same B3 block
b_cells <- bind_rows(
  b_cells,
  tibble(
    cell = c("b3_corr_f2f4_k6_pearson", "b3_corr_f2f4_k12_pearson"),
    value = c(checked(pick_fc(6, "pearson_F2_F4"), "0.7970",
                      "b3_corr_f2f4_k6_pearson"),
              checked(pick_fc(12, "pearson_F2_F4"), "0.9138",
                      "b3_corr_f2f4_k12_pearson")),
    source = paste0("dispb_gate.json B3.feature_collinearity real_form1 ",
                    "loss_avoid_adj K", c(6, 12), " pearson_F2_F4")
  )
) |> arrange(cell)

write_cells(b_cells, er("dispb", "b_cells.csv"))

# ---- 3. dispc/c_cells.csv ---------------------------------------------------

cj <- fromJSON(er("dispc", "dispc_gate.json"), simplifyVector = FALSE)

# G-C2.predicted_vs_observed is a list of per-K records; select K = 12.
pvo <- cj[["G-C2"]]$predicted_vs_observed
pick_k12 <- function(field) {
  hit <- Filter(function(r) as.character(r$K) == "12", pvo)
  stopifnot(length(hit) == 1L)
  hit[[1]][[field]]
}
# G-C3.per_K is a list of per-K records; select K = 6.
perK <- cj[["G-C3"]]$per_K
pick_c3 <- function(K, field) {
  hit <- Filter(function(r) as.character(r$K) == as.character(K), perK)
  stopifnot(length(hit) == 1L)
  hit[[1]][[field]]
}

c_spec <- list(
  list("c1_r2_loss_on_style", "0.6900", c("G-C1", "r2_loss_on_style")),
  list("c1_icc_raw", "0.6458", c("G-C1", "primary_cell", "ICC_raw")),
  list("c1_icc_resid", "0.3954", c("G-C1", "primary_cell", "ICC_resid")),
  list("c1_rel_icc_drop", "0.3876",
       c("G-C1", "primary_cell", "relative_ICC_drop")),
  list("c1_sigma_raw", "0.2186", c("G-C1", "primary_cell", "sigma_true_raw")),
  list("c1_sigma_resid", "0.2240",
       c("G-C1", "primary_cell", "sigma_true_resid")),
  list("c2_dep_adj_icc_reduction", "0.5548",
       c("G-C2", "dependence_adjusted_ICC_reduction")),
  list("c2_r1", "0.9374", c("G-C2", "single_problem_reliability_r1")),
  list("c2_spearman_rsb_wpc", "0.8091", c("G-C2", "spearman_rsb_wpc_k12")),
  list("c3_ranking_spearman_k12", "0.9982", c("G-C3", "ranking_spearman_k12"))
)

c_cells <- bind_rows(lapply(c_spec, function(s) {
  tibble(cell = s[[1]],
         value = checked(jget(cj, s[[3]]), s[[2]], s[[1]]),
         source = paste0("dispc_gate.json ", paste(s[[3]], collapse = ".")))
}))

c_cells <- bind_rows(
  c_cells,
  tibble(
    cell = c("c2_predicted_r_sb_k12", "c2_observed_r_sb_k12",
             "c3_ranking_spearman_k6"),
    value = c(checked(pick_k12("predicted_r_sb"), "0.9794",
                      "c2_predicted_r_sb_k12"),
              checked(pick_k12("observed_r_sb"), "0.4802",
                      "c2_observed_r_sb_k12"),
              checked(pick_c3(6, "ranking_spearman"), "0.9997",
                      "c3_ranking_spearman_k6")),
    source = c("dispc_gate.json G-C2.predicted_vs_observed[K=12].predicted_r_sb",
               "dispc_gate.json G-C2.predicted_vs_observed[K=12].observed_r_sb",
               "dispc_gate.json G-C3.per_K[K=6].ranking_spearman")
  )
) |> arrange(cell)

write_cells(c_cells, er("dispc", "c_cells.csv"))

# ---- 4. sens/sens_cells.csv -------------------------------------------------
# Source: the SENS attenuation-bound output (29_mechanism_sensitivity.R),
# K = 12 rows. Cell names differ from the source column names; the mapping is
# recorded in the `source` column. Nothing is recomputed.

ab <- read_csv(er("sens", "attenuation_bound.csv"),
               col_types = cols(.default = col_character()))

ab_cell <- function(comp, K, col, extra = NULL) {
  d <- ab |> filter(component == comp, K == !!as.character(K))
  if (!is.null(extra)) {
    for (nm in names(extra)) d <- d[d[[nm]] == extra[[nm]], ]
  }
  stopifnot(nrow(d) == 1L)
  d[[col]]
}

sens_cells <- tibble(
  cell = c("atten_max_observable_sb_k12", "atten_max_observable_raw_k12",
           "atten_rtrue_for_030_sb", "real_anchor_k12",
           "real_anchor_ceiling_k12"),
  value = c(
    checked(ab_cell("per_budget", 12, "max_observable_sb_projected"),
            "0.8025", "atten_max_observable_sb_k12"),
    checked(ab_cell("per_budget", 12, "max_observable_spearman"),
            "0.6891", "atten_max_observable_raw_k12"),
    checked(ab_cell("required_r_true_K12", 12, "r_true_required",
                    list(target_oof = "0.3", rho_basis = "sb_projected_full63")),
            "0.3739", "atten_rtrue_for_030_sb"),
    checked(ab_cell("per_budget", 12, "committed_real_oof_pooled"),
            "0.7464", "real_anchor_k12"),
    checked(ab_cell("per_budget", 12, "max_observable_spearman"),
            "0.6891", "real_anchor_ceiling_k12")
  ),
  source = c(
    "sens/attenuation_bound.csv per_budget K=12 max_observable_sb_projected",
    "sens/attenuation_bound.csv per_budget K=12 max_observable_spearman",
    paste0("sens/attenuation_bound.csv required_r_true_K12 K=12 ",
           "target_oof=0.3 rho_basis=sb_projected_full63 r_true_required"),
    paste0("sens/attenuation_bound.csv per_budget K=12 ",
           "committed_real_oof_pooled (E4 real-anchor association at K12)"),
    paste0("sens/attenuation_bound.csv per_budget K=12 ",
           "max_observable_spearman = sqrt(rho_o_out_of_subject), the ",
           "association ceiling the outcome noise permits")
  )
)

# the ceiling cell is by construction sqrt of the committed out-of-subject
# outcome reliability; assert that identity against the same committed row
stopifnot(isTRUE(all.equal(
  as.numeric(sens_cells$value[sens_cells$cell == "real_anchor_ceiling_k12"]),
  sqrt(as.numeric(ab_cell("per_budget", 12, "rho_o_out_of_subject"))),
  tolerance = 1e-12
)))

write_cells(sens_cells, er("sens", "sens_cells.csv"))

# ---- 5. p3_frey/provenance.csv ---------------------------------------------
# Surfaced from the committed p3_frey/provenance_manifest.csv `_n` row, which
# records the analysed Basel-Berlin main-sample size.

pm <- read_csv(er("p3_frey", "provenance_manifest.csv"),
               col_types = cols(.default = col_character()))
n_row <- pm$role[pm$file == "_n"]
stopifnot(length(n_row) == 1L)
n_main <- sub("^main N=([0-9]+).*$", "\\1", n_row)
stopifnot(grepl("^[0-9]+$", n_main))
n_main <- checked(n_main, "1507", "p3 sample n_subjects")

# cross-check against the committed per-measure split-half sample size
dd <- read_csv(er("p3_frey", "dfd_dfe_reliability.csv"),
               col_types = cols(.default = col_character()))
stopifnot(any(dd$splithalf_n == n_main))

prov <- tibble(
  item = "sample",
  n_subjects = n_main,
  source = paste0("p3_frey/provenance_manifest.csv row `_n` (main N), ",
                  "cross-checked against dfd_dfe_reliability.csv splithalf_n")
)
dir.create(dirname(er("p3_frey", "provenance.csv")), showWarnings = FALSE,
           recursive = TRUE)
write_csv(prov, er("p3_frey", "provenance.csv"), na = "NA")

# ---- 6. p3_frey/modelpar_reliability.csv - append two summary rows ----------
# Additive: the 50 committed per-parameter rows keep every existing cell
# byte-identical; one column (n_params) and two summary rows are added.

mp_path <- er("p3_frey", "modelpar_reliability.csv")
mp <- read_csv(mp_path, col_types = cols(.default = col_character()))
# idempotent: strip any summary rows/column a previous run added, then rebuild
mp <- mp |> filter(!parameter %in% c("gate_median", "gate_family"))
mp$n_params <- NULL
mp_before <- mp
stopifnot(nrow(mp) == 50L)

gate <- mp |> filter(gate_param == "TRUE")
gate_median <- checked(median(as.numeric(gate$retest_spearman)), "0.3443",
                       "p3 gate_median retest_spearman")
gate_n <- checked(nrow(gate), "30", "p3 gate_family n_params")

mp_out <- mp |>
  mutate(n_params = NA_character_) |>
  bind_rows(tibble(
    parameter = c("gate_median", "gate_family"),
    task = NA_character_, model = NA_character_, param = NA_character_,
    param_class = c(
      "summary: median retest_spearman over the gate_param rows",
      "summary: count of gate_param rows"),
    gate_param = NA_character_,
    n_retest_pairs = NA_character_,
    retest_spearman = c(gate_median, NA_character_),
    n_params = c(gate_n, gate_n)
  ))

write_csv(mp_out, mp_path, na = "NA")

# every pre-existing cell must be byte-identical after the rewrite
mp_after <- read_csv(mp_path, col_types = cols(.default = col_character()))
stopifnot(nrow(mp_after) == 52L)
for (nm in names(mp_before)) {
  stopifnot(identical(mp_before[[nm]], head(mp_after[[nm]], 50L)))
}

cat("e11 surfacing OK: 5 cells files written; modelpar summary rows appended;",
    "every value copied from a committed source and self-checked\n")

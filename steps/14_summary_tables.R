# Step 14: trials-needed table and claim-status table
# Reads:   data/claim_prose.csv; output/results/step0/ to step4/ (steps 01-12)
# Writes:  output/results/step5/trials_needed_combined.csv, output/results/step5/claim_status_skeleton.csv
# Paper:   Results (trials-needed values and thresholds); Appendix and Supplement (claim-status tables)
# Runtime: a few seconds

# Combine trials-needed values by construct, pool, regime and threshold
# (.5/.7/.8), and assemble one row per quantitative claim.
#
# Values come from saved result CSVs. Prior-study values use the byte-identical
# copy of data/prior_studies.csv written by step 08. Operations are limited
# to filtering, finding the first saved k at each reliability threshold, and
# converting units (total problems = 2 * k; total trials = total problems * r).
# The threshold rule matches step 05 and is checked against the pre-declared
# answers in output/results/step4/extended_interpretation.csv. There is no
# random sampling or new estimation.
#
# Interpretation text from data/claim_prose.csv is joined by claim_id. The run
# fails if a claim lacks prose or a prose row names no claim. That text contains
# no quantitative values; each claim_value comes from its source CSV.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(here)
})

res_dir <- here("output", "results", "step5")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

read_chr <- function(...) read_csv(here("output", ...), col_types = cols(.default = col_character()))

needed   <- read_chr("results", "step2", "trials_needed_table.csv")
ext      <- read_chr("results", "step4", "reliability_extended.csv")
interp   <- read_chr("results", "step4", "extended_interpretation.csv")
pool_sum <- read_chr("results", "step4", "expanded_pool_summary.csv")

pool_clean_lbl <- paste0("clean (", pool_sum$n_problems_clean_step2,
                         " two-outcome no-ambiguity problems)")
pool_ext_lbl   <- paste0("expanded (", pool_sum$n_problems_expanded,
                         " Amb==0 problems, reconstruction-validated)")

# ============================================================================
# Combined trials-needed table
# ============================================================================

# Clean pool: the step 05 table, carried verbatim (every column
# string-identical; 'not reached' entries preserved exactly), tagged with pool
# and regime.
combined_clean <- needed |>
  mutate(
    pool   = pool_clean_lbl,
    regime = if_else(axis == "i_desc", "description (r=5)", "experience (r=20)"),
    source_file = "results/step2/trials_needed_table.csv"
  )

# Expanded pool: assembled from output/results/step4/reliability_extended.csv
# using the same rule that built the step 05 table: the first realized k whose
# saved r_sb value meets the threshold; never extrapolate beyond the
# realized grid (k <= 13).
constructs <- c("b_ev", "b_loss", "b_risk", "ev_sens_agn", "loss_avoid_adj")
thresholds <- c("0.5", "0.7", "0.8")

ext_verified <- ext |>
  filter(axis == "i_exp_ext", regime == "exp", status == "VERIFIED") |>
  mutate(k_int = as.integer(k), r_num = as.numeric(r_sb)) |>
  arrange(construct, k_int)
stopifnot(identical(sort(unique(ext_verified$k_int[ext_verified$construct == "loss_avoid_adj"])), 2:13))

ext_rows <- crossing(construct = constructs, threshold = thresholds) |>
  rowwise() |>
  do({
    co  <- .$construct; thr <- .$threshold
    d   <- ext_verified |> filter(construct == co, r_num >= as.numeric(thr))
    if (nrow(d) == 0) {
      tibble(
        axis = "i_exp_ext",
        axis_label = "problems (experience, r=20, expanded pool)",
        construct = co, threshold = thr,
        budget = NA_character_, total_problems = NA_character_,
        total_trials = NA_character_, r_sb = NA_character_,
        ci_lo_sb = NA_character_, ci_hi_sb = NA_character_,
        n_subjects = NA_character_,
        entry = "not reached within realized grid (k<=13)"
      )
    } else {
      hit <- d |> slice_min(k_int, n = 1)
      tibble(
        axis = "i_exp_ext",
        axis_label = "problems (experience, r=20, expanded pool)",
        construct = co, threshold = thr,
        budget = hit$k,                                    # verbatim string
        total_problems = as.character(2L * hit$k_int),     # unit conversion
        total_trials   = as.character(2L * hit$k_int * 20L),
        r_sb = hit$r_sb, ci_lo_sb = hit$ci_lo_sb, ci_hi_sb = hit$ci_hi_sb,
        n_subjects = hit$n_subjects,
        entry = "reached"
      )
    }
  }) |>
  ungroup() |>
  mutate(pool = pool_ext_lbl, regime = "experience (r=20)",
         source_file = "results/step4/reliability_extended.csv")

# Cross-assert the expanded crossings against the PRE-DECLARED answers in
# output/results/step4/extended_interpretation.csv (step 11).
ans <- function(q) interp$answer[interp$question == q]
stopifnot(
  ext_rows$budget[ext_rows$construct == "loss_avoid_adj" & ext_rows$threshold == "0.5"] ==
    ans("loss_avoid_adj crosses SB .5 at k (problems/half)"),
  ext_rows$entry[ext_rows$construct == "loss_avoid_adj" & ext_rows$threshold == "0.7"] ==
    "not reached within realized grid (k<=13)",
  all(vapply(constructs, function(co)
    ans(paste0(co, " reaches SB .7 within realized grid (k<=13)")) == "no" &&
      ext_rows$entry[ext_rows$construct == co & ext_rows$threshold == "0.7"] ==
        "not reached within realized grid (k<=13)",
    logical(1)))
)
# Nothing reaches .8 on the expanded grid either.
stopifnot(all(ext_rows$entry[ext_rows$threshold == "0.8"] ==
                "not reached within realized grid (k<=13)"))

combined <- bind_rows(combined_clean, ext_rows) |>
  select(pool, regime, axis, axis_label, construct, threshold, budget,
         total_problems, total_trials, r_sb, ci_lo_sb, ci_hi_sb, n_subjects,
         entry, source_file)

write_csv(combined, file.path(res_dir, "trials_needed_combined.csv"), na = "NA")

# ---- Verbatim assertions on the file as written to disk
reread <- read_csv(file.path(res_dir, "trials_needed_combined.csv"),
                   col_types = cols(.default = col_character()))
# (a) Clean rows byte-equal to the step 05 table, column by column.
needed_again <- read_chr("results", "step2", "trials_needed_table.csv")
clean_back <- reread |>
  filter(source_file == "results/step2/trials_needed_table.csv") |>
  select(all_of(names(needed_again)))
stopifnot(identical(as.data.frame(clean_back), as.data.frame(needed_again)))
# (b) Every 'reached' expanded row's r_sb / CI strings identical to the source.
ext_again <- read_chr("results", "step4", "reliability_extended.csv")
for (i in which(reread$source_file == "results/step4/reliability_extended.csv" &
                reread$entry == "reached")) {
  src <- ext_again |> filter(k == reread$budget[i], construct == reread$construct[i])
  stopifnot(nrow(src) == 1L,
            identical(src$r_sb,     reread$r_sb[i]),
            identical(src$ci_lo_sb, reread$ci_lo_sb[i]),
            identical(src$ci_hi_sb, reread$ci_hi_sb[i]),
            identical(src$n_subjects, reread$n_subjects[i]))
}
n_reached_ext <- sum(reread$entry == "reached" &
                       reread$source_file == "results/step4/reliability_extended.csv")

# ============================================================================
# Claim-status skeleton (prose joined from data/claim_prose.csv)
# ============================================================================
# Every claim_value is pulled programmatically from its source CSV; the
# claim() helper takes the already-loaded data frame cell, so a hand-typed
# value cannot enter. Statuses: VERIFIED (computed here from md5-verified real
# data), REPRODUCED (matched a named external source), SIMULATION
# (from simulated data; recovery gates), EXPLORATORY (labeled exploratory).

s0_rec  <- read_chr("results", "step0", "recovery_verdict.csv")
s0_rate <- read_chr("results", "step0", "sanity_b_rate.csv")
s0_mono <- read_chr("results", "step0", "sanity_ev_monotonicity_bins.csv")
s1_agr  <- read_chr("results", "step1", "agreement_bar_verdict.csv")
s1a     <- read_chr("results", "step1_addendum", "cross_half_agreement.csv")
s2_fix  <- read_chr("results", "step2", "fixed_pool_check.csv")
s2_rel  <- read_chr("results", "step2", "reliability_by_cell_construct.csv")
studies <- read_chr("results", "step3", "prior_studies_table.csv")
s4_gate <- read_chr("results", "step4", "reconstruction_validation.csv")
s4_off  <- read_chr("results", "step4", "official_code_verification.csv")
s4_rec  <- read_chr("results", "step4", "expanded_recovery_verdict.csv")
c13_chk <- read_chr("results", "step4", "choices13k_aggregate_checks.csv")
c13_exp <- read_chr("results", "step4", "choices13k_exploratory_loss.csv")

claims <- list()
claim <- function(claim_id, claim_value, source_file, source_key, status) {
  stopifnot(length(claim_value) == 1L, !is.na(claim_value))
  tibble(claim_id = claim_id, claim_value = claim_value,
         source_file = source_file, source_key = source_key,
         status = status, prose = "")
}

# ---- step0_ claims (step 01): provenance sanity + recovery gate (recovery is SIMULATION)
claims$s0_rate <- claim(
  "step0_overall_b_rate_all_trials_near_.5",
  s0_rate$b_rate[s0_rate$scope == "all trials"],
  "results/step0/sanity_b_rate.csv", "scope=all trials; col b_rate", "VERIFIED")
claims$s0_mono_lo <- claim(
  "step0_ev_monotonicity_b_rate_lowest_ev_bin",
  s0_mono$b_rate[s0_mono$ev_bin == "1"],
  "results/step0/sanity_ev_monotonicity_bins.csv", "ev_bin=1; col b_rate", "VERIFIED")
claims$s0_mono_hi <- claim(
  "step0_ev_monotonicity_b_rate_highest_ev_bin",
  s0_mono$b_rate[s0_mono$ev_bin == "10"],
  "results/step0/sanity_ev_monotonicity_bins.csv", "ev_bin=10; col b_rate", "VERIFIED")
for (p in s0_rec$param) {
  claims[[paste0("s0_rec_", p)]] <- claim(
    paste0("step0_recovery_spearman_", p, "_at_14x25_(gate_>=.7_",
           if (s0_rec$recoverable_at_14x25[s0_rec$param == p] == "TRUE") "PASS" else "FAIL", ")"),
    s0_rec$spearman_r[s0_rec$param == p],
    "results/step0/recovery_verdict.csv", paste0("param=", p, "; col spearman_r"),
    "SIMULATION")
}

# ---- step1_ claims (step 04): agreement verdicts (real data; FAILs are results)
for (i in seq_len(nrow(s1_agr))) {
  claims[[paste0("s1_agr_", i)]] <- claim(
    paste0("step1_agreement_spearman_", s1_agr$agnostic[i], "_vs_",
           s1_agr$model_based[i], "_(bar_",
           if (s1_agr$pass[i] == "TRUE") "PASS" else "FAIL", ")"),
    s1_agr$r[i],
    "results/step1/agreement_bar_verdict.csv",
    paste0("agnostic=", s1_agr$agnostic[i], "; col r"),
    "VERIFIED")
}

# ---- step1_addendum_ claims (step 06): cross-half vs same-half agreement
for (i in seq_len(nrow(s1a))) {
  claims[[paste0("s1a_", i)]] <- claim(
    paste0("step1_addendum_", s1a$measure[i], "_agreement_", s1a$pair[i]),
    s1a$rho[i],
    "results/step1_addendum/cross_half_agreement.csv",
    paste0("pair=", s1a$pair[i], "; measure=", s1a$measure[i], "; col rho"),
    "VERIFIED")
}

# ---- step2_ claims (steps 05, 07): the wall, the .5 crossings, the k=5 artifact
wall_rows <- needed |> filter(axis == "i_exp", threshold == "0.7")
stopifnot(nrow(wall_rows) == 5L,
          all(wall_rows$entry == "not reached within CPC18 budgets"))
claims$s2_wall <- claim(
  "step2_wall_no_construct_reaches_SB.7_clean_pool_(all_5_constructs)",
  unique(wall_rows$entry),
  "results/step2/trials_needed_table.csv",
  "axis=i_exp; threshold=0.7; col entry (identical for all constructs)",
  "VERIFIED")
row_c <- needed |> filter(axis == "i_exp", construct == "loss_avoid_adj", threshold == "0.5")
claims$s2_cross_total <- claim(
  "step2_loss_avoid_adj_SB.5_crossing_clean_total_problems",
  row_c$total_problems,
  "results/step2/trials_needed_table.csv",
  "axis=i_exp; construct=loss_avoid_adj; threshold=0.5; col total_problems",
  "VERIFIED")
claims$s2_cross_rsb <- claim(
  "step2_loss_avoid_adj_SB_at_clean_crossing",
  row_c$r_sb,
  "results/step2/trials_needed_table.csv",
  "axis=i_exp; construct=loss_avoid_adj; threshold=0.5; col r_sb",
  "VERIFIED")
row_d <- needed |> filter(axis == "i_desc", construct == "loss_avoid_adj", threshold == "0.5")
claims$s2_desc_total <- claim(
  "step2_loss_avoid_adj_SB.5_crossing_description_total_problems",
  row_d$total_problems,
  "results/step2/trials_needed_table.csv",
  "axis=i_desc; construct=loss_avoid_adj; threshold=0.5; col total_problems",
  "VERIFIED")
# k=5 artifact: on the fixed problem pool the curves are monotone in k;
# computed (not typed) from the step 07 fixed-pool CSV.
mono <- function(co) {
  v <- s2_fix |> filter(construct == co) |>
    mutate(k = as.integer(k), r = as.numeric(r_sb)) |> arrange(k)
  all(diff(v$r) > 0)
}
stopifnot(mono("ev_sens_agn"), mono("loss_avoid_adj"))
claims$s2_k5 <- claim(
  "step2_k5_dip_is_set_composition_artifact_(fixed_pool_SB_monotone_k2-7)",
  as.character(mono("ev_sens_agn") && mono("loss_avoid_adj")),
  "results/step2/fixed_pool_check.csv",
  "computed: diff(r_sb) > 0 over k=2..7 for ev_sens_agn AND loss_avoid_adj",
  "VERIFIED")

# ---- step3_ claims (step 08): prior-study facts (verbatim, carrying the table's verification labels)
for (i in seq_len(nrow(studies))) {
  claims[[paste0("s3_", studies$study_id[i])]] <- claim(
    paste0("step3_prior_study_", studies$study_id[i], "_",
           str_replace_all(studies$reliability_type[i], "[^A-Za-z0-9.]+", "_")),
    studies$reliability_value[i],
    "results/step3/prior_studies_table.csv",
    paste0("study_id=", studies$study_id[i], "; col reliability_value"),
    studies$verification[i])
}

# ---- step4_ claims (steps 09-11): reconstruction gate, recovery re-gate, wall at 26, 12->24, b_ev rise
stopifnot(all(s4_gate$pass == "TRUE"), all(s4_off$verified == "TRUE"))
claims$s4_gate <- claim(
  "step4_reconstruction_exact_match_official_CPC18_implementation_(gate_PASS)",
  s4_gate$pass[s4_gate$check == "1_exact_match_official_implementation"],
  "results/step4/reconstruction_validation.csv",
  "check=1_exact_match_official_implementation; col pass",
  "REPRODUCED")
claims$s4_md5 <- claim(
  "step4_official_code_and_data_artifacts_all_md5_verified",
  as.character(all(s4_off$verified == "TRUE")),
  "results/step4/official_code_verification.csv",
  "computed: all(verified == TRUE) over 8 pinned artifacts",
  "VERIFIED")
for (p in s4_rec$param) {
  claims[[paste0("s4_rec_", p)]] <- claim(
    paste0("step4_recovery_regate_spearman_", p, "_expanded_design_(gate_>=.7_",
           if (s4_rec$pass[s4_rec$param == p] == "TRUE") "PASS" else "FAIL", ")"),
    s4_rec$spearman_r[s4_rec$param == p],
    "results/step4/expanded_recovery_verdict.csv",
    paste0("param=", p, "; col spearman_r"),
    "SIMULATION")
}
for (i in seq_len(nrow(interp))) {
  claims[[paste0("s4_interp_", i)]] <- claim(
    paste0("step4_expanded_", str_replace_all(interp$question[i], "[^A-Za-z0-9.<=]+", "_")),
    interp$answer[i],
    "results/step4/extended_interpretation.csv",
    paste0("question=", interp$question[i], "; col answer"),
    "VERIFIED")
}
bev13 <- ext |> filter(cell_id == "ext_exp_k13_r20", construct == "b_ev")
claims$s4_bev_ext <- claim(
  "step4_b_ev_SB_expanded_pool_at_k13_(rise_endpoint,_still_<.7)",
  bev13$r_sb,
  "results/step4/reliability_extended.csv",
  "cell_id=ext_exp_k13_r20; construct=b_ev; col r_sb",
  "VERIFIED")
bev7 <- s2_rel |> filter(cell_id == "i_exp_k7_r20", construct == "b_ev")
claims$s4_bev_clean <- claim(
  "step4_b_ev_SB_clean_pool_at_k7_(rise_baseline)",
  bev7$r_sb,
  "results/step2/reliability_by_cell_construct.csv",
  "cell_id=i_exp_k7_r20; construct=b_ev; col r_sb",
  "VERIFIED")

# ---- step4_choices13k claims (step 12; aggregate-only; exploratory rows labeled EXPLORATORY)
for (i in seq_len(nrow(c13_chk))) {
  claims[[paste0("c13_", i)]] <- claim(
    paste0("step4_choices13k_", str_replace_all(
      str_extract(c13_chk$check[i], "^B2\\.\\d+ [^:,]+"), "[^A-Za-z0-9.]+", "_"),
      "_(aggregate_only,_no_per-subject_claims)"),
    c13_chk$value[i],
    "results/step4/choices13k_aggregate_checks.csv",
    paste0("check=", c13_chk$check[i], "; col value"),
    c13_chk$status[i])
}
for (i in seq_len(nrow(c13_exp))) {
  claims[[paste0("c13x_", i)]] <- claim(
    paste0("step4_choices13k_exploratory_loss_direction_",
           str_replace_all(c13_exp$dataset[i], "[^A-Za-z0-9]+", "_"),
           "_(direction_only,_aggregate_only)"),
    c13_exp$direction_at_ev_parity[i],
    "results/step4/choices13k_exploratory_loss.csv",
    paste0("dataset=", c13_exp$dataset[i], "; col direction_at_ev_parity"),
    "EXPLORATORY")
}

skeleton <- bind_rows(claims)
stopifnot(!any(duplicated(skeleton$claim_id)), all(skeleton$prose == ""))

# Join the interpretation prose: one non-empty prose row per claim, no orphans
# in either direction.
claim_prose <- read_csv(here("data", "claim_prose.csv"),
                        col_types = cols(.default = col_character()),
                        na = character(0))
stopifnot(identical(names(claim_prose), c("claim_id", "prose")),
          !anyDuplicated(claim_prose$claim_id))
no_prose <- setdiff(skeleton$claim_id, claim_prose$claim_id[nzchar(claim_prose$prose)])
no_claim <- setdiff(claim_prose$claim_id, skeleton$claim_id)
if (length(no_prose) > 0) stop("step5b: claim(s) without prose: ",
                               paste(no_prose, collapse = "; "))
if (length(no_claim) > 0) stop("step5b: prose row(s) without a claim: ",
                               paste(no_claim, collapse = "; "))
skeleton$prose <- claim_prose$prose[match(skeleton$claim_id, claim_prose$claim_id)]
stopifnot(!anyNA(skeleton$prose), all(nzchar(skeleton$prose)))

write_csv(skeleton, file.path(res_dir, "claim_status_skeleton.csv"), na = "")

# ---- Verification pass: re-read the skeleton from disk and re-read every
# source file fresh; assert each claim_value appears verbatim among the source
# file's cell values (computed booleans are re-derived instead).
sk <- read_csv(file.path(res_dir, "claim_status_skeleton.csv"),
               col_types = cols(.default = col_character())) |>
  mutate(prose = if_else(is.na(prose), "", prose))
stopifnot(nrow(sk) == nrow(skeleton),
          identical(sk$prose,
                    claim_prose$prose[match(sk$claim_id, claim_prose$claim_id)]))
computed_ids <- c(
  "step2_k5_dip_is_set_composition_artifact_(fixed_pool_SB_monotone_k2-7)",
  "step4_official_code_and_data_artifacts_all_md5_verified")
src_cache <- list()
for (i in seq_len(nrow(sk))) {
  f <- sk$source_file[i]
  if (is.null(src_cache[[f]])) {
    src_cache[[f]] <- read_csv(here("output", f), col_types = cols(.default = col_character()))
  }
  if (sk$claim_id[i] %in% computed_ids) next  # re-derived via stopifnot above
  cells <- unlist(src_cache[[f]], use.names = FALSE)
  stopifnot(sk$claim_value[i] %in% cells)
}
cat("step5b OK:", nrow(combined), "combined trials-needed rows (",
    n_reached_ext, "expanded 'reached');", nrow(sk),
    "skeleton claims, all values re-verified against sources; prose joined",
    "from data/claim_prose.csv\n")

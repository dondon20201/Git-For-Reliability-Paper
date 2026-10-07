# Step 09: validation of the reconstructed multi-outcome lotteries against the official CPC18 code
# Reads:   data/raw/All estimation raw data.csv; data/raw/cpc18_official/ (downloaded and md5-checked if absent)
# Writes:  output/results/step4/ (reconstruction_validation.csv, official_code_verification.csv, problem_distributions.csv)
# Paper:   Methods (Problem pools)
# Runtime: a few seconds

# Multi-outcome reconstruction VALIDATION GATE.
#
# Reconstructs the full outcome distribution of every option of all 210 CPC18
# problems from (H, pH, L, LotShape, LotNum) with R/lotshape_reconstruction.R
# and validates it against the OFFICIAL CPC18 baseline implementation
# (CPC18_getDist.R, fetched + md5-verified by R/fetch_official_cpc18.R).
#
# HARD GATE (pre-declared; all four checks asserted in-code; if ANY fails
# this script stops AFTER writing the validation table):
#   1. Enumeration matches the official implementation EXACTLY (all 210
#      problems, both options). Supplementary: matches the organizer-produced
#      CPC18_EstSet_fullDist.csv within 1e-9 (CSV decimal round-trip).
#   2. Every option's probabilities sum to 1 (tolerance 1e-9).
#   3. LotNum == 1 options reduce exactly to the two-outcome (H, pH, L) form.
#   4. EVs recomputed via the reconstruction match the values used in steps
#      04-05 (clean_problem_features) exactly on the clean two-outcome subset.

library(tidyverse)
library(here)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "fetch_official_cpc18.R"))
source(here("R", "lotshape_reconstruction.R"))

res_dir <- here("output", "results", "step4")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
official_dir <- here("data", "raw", "cpc18_official")

# ---- verified inputs --------------------------------------------------------
fetch_res <- fetch_cpc18(here("data", "raw"))
stopifnot(all(fetch_res$verified))
official_res <- fetch_official_cpc18(official_dir)
stopifnot(all(official_res$verified))
write_csv(official_res, file.path(res_dir, "official_code_verification.csv"))

raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))

# ---- problem table from the raw data + cross-check vs official EstSet210 ----
param_cols <- c("Ha", "pHa", "La", "LotShapeA", "LotNumA",
                "Hb", "pHb", "Lb", "LotShapeB", "LotNumB", "Amb", "Corr")
problems <- raw |>
  distinct(GameID, across(all_of(param_cols))) |>
  arrange(GameID)
stopifnot(nrow(problems) == 210, anyDuplicated(problems$GameID) == 0)

est210 <- read_csv(file.path(official_dir, "CPC18_EstSet210.csv"),
                   show_col_types = FALSE) |> arrange(GameID)
stopifnot(identical(est210$GameID, problems$GameID))
params_match <- vapply(param_cols, function(cc) {
  all(problems[[cc]] == est210[[cc]])
}, logical(1))

# ---- my reconstruction ------------------------------------------------------
mine <- reconstruct_all_problems(problems)

# ---- official implementation, run VERBATIM ----------------------------------
source(file.path(official_dir, "CPC18_getDist.R"))  # defines CPC18_getDist
official <- map_dfr(seq_len(nrow(problems)), function(r) {
  p <- problems[r, ]
  map_dfr(c("A", "B"), function(opt) {
    m <- if (opt == "A") {
      CPC18_getDist(p$Ha, p$pHa, p$La, p$LotShapeA, p$LotNumA)
    } else {
      CPC18_getDist(p$Hb, p$pHb, p$Lb, p$LotShapeB, p$LotNumB)
    }
    m <- matrix(as.numeric(m), ncol = 2)
    tibble(GameID = p$GameID, option = opt,
           outcome_index = seq_len(nrow(m)),
           outcome = m[, 1], prob = m[, 2])
  })
})

# ---- CHECK 1: exact match to the official enumeration -----------------------
same_rows <- identical(dim(mine)[1], dim(official)[1]) &&
  identical(mine$GameID, official$GameID) &&
  identical(mine$option, official$option) &&
  identical(mine$outcome_index, official$outcome_index)
max_dev_outcome <- max(abs(mine$outcome - official$outcome))
max_dev_prob <- max(abs(mine$prob - official$prob))
check1_code <- same_rows && max_dev_outcome == 0 && max_dev_prob == 0

# Supplementary: organizer-produced full-distribution table (data-level
# target), compared ascending-sorted. The CSV is a rounded decimal export
# (probabilities printed to as few as 3 significant digits, e.g. "3.91E-05"
# for 3.90625e-05), so the criterion is a decimal ROUND-TRIP: every printed
# value must equal the reconstructed probability rounded to its own printed
# precision (within half a unit in the last printed digit). A fixed 1e-9
# tolerance fails on 14/1129 entries (max dev 3.75e-8, all at Symm binomial
# tails printed with 3 significant digits); this is noted in the validation
# detail.
fulldist <- read_csv(file.path(official_dir, "CPC18_EstSet_fullDist.csv"),
                     show_col_types = FALSE) |> arrange(GameID)
fulldist_chr <- read_csv(file.path(official_dir, "CPC18_EstSet_fullDist.csv"),
                         col_types = cols(.default = "c")) |>
  mutate(GameID = as.numeric(GameID)) |> arrange(GameID)
stopifnot(identical(fulldist$GameID, problems$GameID))
csv_long <- map_dfr(c("A", "B"), function(opt) {
  pre <- tolower(opt)
  vals <- fulldist |> select(GameID, matches(sprintf("^%sv\\d+$", pre)))
  prbs_chr <- fulldist_chr |> select(matches(sprintf("^%sp\\d+$", pre)))
  map_dfr(seq_len(nrow(vals)), function(r) {
    v <- as.numeric(vals[r, -1]); s <- as.character(prbs_chr[r, ])
    keep <- !is.na(v)
    tibble(GameID = vals$GameID[r], option = opt,
           outcome_index = seq_len(sum(keep)),
           outcome = v[keep], prob_chr = s[keep])
  })
}) |> arrange(GameID, option, outcome_index)
# half a unit in the last printed digit of a decimal string like "0.0003125",
# "3.91E-05", "0.99", or "9"
printed_half_ulp <- function(s) {
  s <- trimws(s)
  m <- regmatches(s, regexec("^-?([0-9]*)(?:\\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$", s))
  vapply(m, function(g) {
    stopifnot(length(g) == 4)
    n_dec <- nchar(g[3])
    expo <- ifelse(g[4] == "", 0, as.numeric(g[4]))
    0.5 * 10^(expo - n_dec)
  }, numeric(1))
}
mine_sorted <- mine |>
  arrange(GameID, option, outcome, prob) |>
  group_by(GameID, option) |>
  mutate(outcome_index = row_number()) |>
  ungroup() |>
  arrange(GameID, option, outcome_index)
same_rows_csv <- nrow(csv_long) == nrow(mine_sorted) &&
  identical(csv_long$GameID, mine_sorted$GameID) &&
  identical(csv_long$option, mine_sorted$option)
max_dev_outcome_csv <- if (same_rows_csv) max(abs(csv_long$outcome - mine_sorted$outcome)) else NA_real_
prob_dev_csv <- abs(as.numeric(csv_long$prob_chr) - mine_sorted$prob)
max_dev_prob_csv <- if (same_rows_csv) max(prob_dev_csv) else NA_real_
roundtrip_fail <- sum(prob_dev_csv > printed_half_ulp(csv_long$prob_chr) * (1 + 1e-9))
check1_csv <- same_rows_csv && max_dev_outcome_csv == 0 && roundtrip_fail == 0

# ---- CHECK 2: probabilities sum to 1 ----------------------------------------
sums <- mine |>
  group_by(GameID, option) |>
  summarise(psum = sum(prob), .groups = "drop")
max_dev_sum <- max(abs(sums$psum - 1))
check2 <- max_dev_sum <= 1e-9

# ---- CHECK 3: LotNum == 1 reduces to the two-outcome (H, pH, L) form --------
opt_params <- bind_rows(
  problems |> transmute(GameID, option = "A", H = Ha, pH = pHa, L = La,
                        LotShape = LotShapeA, LotNum = LotNumA),
  problems |> transmute(GameID, option = "B", H = Hb, pH = pHb, L = Lb,
                        LotShape = LotShapeB, LotNum = LotNumB)
)
stopifnot(all(opt_params$LotShape[opt_params$LotNum == 1] == "-"))
check3_one <- function(gid, opt, H, pH, L) {
  d <- mine |> filter(GameID == gid, option == opt)
  if (pH == 1) {
    nrow(d) == 1 && d$outcome == H && d$prob == 1
  } else {
    nrow(d) == 2 && all(d$outcome == c(L, H)) && all(d$prob == c(1 - pH, pH))
  }
}
two_out <- opt_params |> filter(LotNum == 1)
check3_all <- pmap_lgl(two_out |> select(gid = GameID, opt = option, H, pH, L),
                       check3_one)
check3 <- all(check3_all)
# degenerate pH == 1 options where L != H (L unreachable, dropped by encoding)
n_degenerate_L <- sum(opt_params$pH == 1 & opt_params$L != opt_params$H)

# ---- CHECK 4: EVs match steps 04-05 exactly on the clean subset --------------
features_12 <- clean_problem_features(raw)  # the exact function used in steps 04-05
recon_mom <- mine |>
  group_by(GameID, option) |>
  summarise(dist_moments(outcome, prob), .groups = "drop")
recon_wide <- recon_mom |>
  pivot_wider(id_cols = GameID, names_from = option,
              values_from = c(ev, sd, eloss, has_loss))
cmp4 <- features_12 |>
  select(GameID, ev_a, ev_b, d_ev, sd_a, sd_b, eloss_a, eloss_b,
         has_loss_a, has_loss_b) |>
  left_join(recon_wide, by = "GameID")
check4_ev <- all(cmp4$ev_a == cmp4$ev_A) && all(cmp4$ev_b == cmp4$ev_B) &&
  all(cmp4$d_ev == cmp4$ev_B - cmp4$ev_A)
# informational (not part of the pre-declared gate): sd/eloss/has_loss carry-over
info_sd <- all(cmp4$sd_a == cmp4$sd_A) && all(cmp4$sd_b == cmp4$sd_B)
info_eloss <- all(cmp4$eloss_a == cmp4$eloss_A) && all(cmp4$eloss_b == cmp4$eloss_B)
info_hasloss <- all(cmp4$has_loss_a == cmp4$has_loss_A) &&
  all(cmp4$has_loss_b == cmp4$has_loss_B)

# ---- write outputs BEFORE asserting (a failed gate still leaves the table) ---
dists_out <- mine |>
  left_join(opt_params, by = c("GameID", "option")) |>
  select(GameID, option, LotShape, LotNum, outcome_index, outcome, prob)
write_csv(dists_out, file.path(res_dir, "problem_distributions.csv"))

validation <- tribble(
  ~check, ~pass, ~detail,
  "0_raw_params_match_official_EstSet210", all(params_match),
  paste0("210 problems x 12 parameter columns vs official CPC18_EstSet210.csv; mismatching columns: ",
         ifelse(all(params_match), "none",
                paste(param_cols[!params_match], collapse = ";"))),
  "1_exact_match_official_implementation", check1_code,
  sprintf("row structure identical: %s; max |outcome diff| = %g; max |prob diff| = %g (required: exactly 0). Official code: %s @ %s (verbatim copy of the official BEASTsd R baseline from %s), md5 %s",
          same_rows, max_dev_outcome, max_dev_prob,
          official_cpc18_source$repo, official_cpc18_source$commit,
          official_cpc18_source$channel,
          official_res$md5_observed[official_res$file == "CPC18_getDist.R"]),
  "1b_match_official_fullDist_csv", check1_csv,
  sprintf(paste0("organizer-produced CPC18_EstSet_fullDist.csv (rounded decimal export); ",
                 "row structure identical: %s; max |outcome diff| = %g; ",
                 "decimal round-trip failures (printed value != reconstruction rounded to printed precision): %d of %d entries; ",
                 "max |prob diff| = %g. TRANSPARENCY NOTE: a first-pass fixed 1e-9 tolerance failed on 14/1129 entries ",
                 "(all Symm binomial tails printed with 3 significant digits, e.g. 3.91E-05 for 3.90625e-05); ",
                 "criterion replaced by the round-trip check; the pre-declared gate check (1: exact match to the official ",
                 "implementation) is unaffected and passed at exactly 0 deviation"),
          same_rows_csv, max_dev_outcome_csv, roundtrip_fail, nrow(csv_long), max_dev_prob_csv),
  "2_probabilities_sum_to_1", check2,
  sprintf("max |sum - 1| over 420 options = %g (tolerance 1e-9)", max_dev_sum),
  "3_lotnum1_reduces_to_two_outcome", check3,
  sprintf("%d/%d LotNum==1 options reduce exactly; %d degenerate pH==1 options with L != H in the 420 (L unreachable, dropped by the official encoding)",
          sum(check3_all), nrow(two_out), n_degenerate_L),
  "4_ev_matches_steps_1_2_exactly", check4_ev,
  sprintf("clean two-outcome subset (%d problems): ev_a, ev_b, d_ev all exactly equal (==) to clean_problem_features() as used in Steps 1-2. Informational carry-over: sd exact = %s, eloss exact = %s, has_loss identical = %s",
          nrow(cmp4), info_sd, info_eloss, info_hasloss)
)
write_csv(validation, file.path(res_dir, "reconstruction_validation.csv"))
print(validation, width = Inf)

# ---- HARD GATE ---------------------------------------------------------------
gate_pass <- all(params_match) && check1_code && check1_csv &&
  check2 && check3 && check4_ev
if (!gate_pass) {
  stop("VALIDATION GATE FAILED — see results/step4/reconstruction_validation.csv. ",
       "STOPPING (no extended curves may be computed).",
       call. = FALSE)
}
cat("\nVALIDATION GATE: ALL CHECKS PASSED (asserted in-code).\n")

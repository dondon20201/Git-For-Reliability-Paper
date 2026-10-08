# Step 27: external re-analysis of the Frey et al. (2017) Basel-Berlin data
# Reads:   data/raw/Frey Data/ (see data/README.md)
# Writes:  output/elevation/results/p3_frey/ (dfd_dfe_reliability.csv, modelpar_reliability.csv,
#          reliability_by_tradition.csv, provenance_manifest.csv, ...)
# Paper:   The reliability ceiling reflects true-score variance (The ceiling replicates in independent data and over time); Methods; Abstract
# Runtime: about 40 seconds, single process

# P3 external replication on the Frey et al. (2017) Basel-Berlin Risk Study.
#
# Uses only the public Frey data (OSF rce7g); no CPC18 data, no CPC18 participant
# split and no CPC18-derived result is read. The analysis plan was fixed before this
# script ran. All reliability statistics are rank-based (Spearman). Split-half uses the frozen
# Spearman-Brown form 2r/(1+r). Test-retest (temporal, 6-month) and split-half (internal)
# are distinguished: every output column carries a retest_ or splithalf_ prefix.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lme4)
  library(jsonlite)
})

t_start <- Sys.time()
set.seed(20260719)

# ---------------------------------------------------------------------------
# 0. Frozen configuration (fixed before any result)
# ---------------------------------------------------------------------------
CFG <- list(
  spec = "elevation/P3_spec.md",
  osf_project = "rce7g (Basel-Berlin Risk Study, Frey et al. 2017, Science Advances)",
  seed_p3 = 20260719,
  seed_boot = 20260720,
  seed_perm = 20260721,
  seed_sim = 20260722,
  n_boot = 2000,
  n_perm = 1000,
  n_sim_ceiling = 20,
  boot_ci_levels = c(0.90, 0.95),
  sb_form = "2r/(1+r) applied to a Spearman half correlation (frozen Study-1 form)",
  linkage_rule = "6-digit core id = drop the trailing 2-digit occasion suffix (main ids end 01, retest ids end 02); retest links to main on the core",
  representative_set = "the 53 cbalpha measures; traditions from grouping$group; the 37 with grouping$category == risk fall in {beh, prop, freq} and form the risk set",
  factor_conversion = "8 shipped measures (FTND, GABS, PG, DAST, NUM, CAREa, CAREs, CAREw) arrive as ordered factors with numeric level labels (the frequency screens are 0/1 in the shipped measures object); converted via as.numeric(as.character(.)), asserted NA-preserving; rank statistics are invariant to this labeling",
  general_factor = list(
    measures = "the 37 risk-category (beh/prop/freq) cbalpha measures",
    method = "first principal component of the pairwise-complete Pearson correlation matrix of T1-standardized measures (all 1507 T1 subjects)",
    orientation = "sign fixed so the mean loading over the six DOSPERT risk-taking scores (Deth..Dsoc) is positive",
    scoring = "score = sum(w_j * z_j) / sum(|w_j|) over measures available for the subject; z uses T1 mean/sd at BOTH time points; T2 scored with T1 loadings",
    min_available = 30,
    score_round_decimals = 12
  ),
  dfd_halves = "odd vs even presentation position (gamble_ind): R_sh1 = mean R at positions 1,3,5,7; R_sh2 = positions 2,4,6,8; recovered exactly from dfd_perprob and asserted against Frey's shipped R_sh1/R_sh2 in main AND retest_berlin",
  dfe_halves = "same odd/even gamble_ind rule applied to dfe_perprob (R and Rexp); half value = mean over non-NA problems in the half; a half needs >= 2 valid problems else NA",
  dfe_min_valid_per_half = 2,
  retest_min_pairs = 10,
  degenerate_flag = "prop of linked T1 subjects at the modal value > 0.90 (Frey trivial-measure caveat); flag only, never a gate",
  degenerate_prop_modal = 0.90,
  material_below_margin = 0.15,
  modelpar_columns = "all modelpars.csv columns except partid and *_DEV / *_BIC fit indices",
  modelpar_gate_family = "tasks ABCT/MPL/CCT, parameter suffix in {OS, OSg, OSl, PW, PWg, PWl, LA} (outcome sensitivity, probability weighting, loss aversion); CS and BART reported but not gated",
  modelpar_gate_rule = "G-param PASS if strictly more than half of the gate-family parameters have 6-month test-retest Spearman < 0.60 (report each)",
  hier_model = "glmer(risky ~ 1 + domain + (1 | partid), family = binomial); risky = shipped R indicator; per-subject estimate = random-intercept conditional mode (BLUP), rounded to 12 decimals before correlation and export",
  hier_optimizer = "default bobyqa path; maxfun 1e5, one pre-declared refit at 1e6 on non-convergence",
  blup_round_decimals = 12,
  perm_flag_rule = "|mean permuted r| > 0.02 AND |z| > 3 (z = mean / (sd/sqrt(n_perm))); same criterion as E8",
  ceiling_sim = "SIMULATION, labeled: draw u_i ~ N(0, sd_hat) from the full-8-problem glmer fit, simulate choices on each subject's actual domain design, recompute raw and hier split-half; 20 sims; model-implied attainable reliability at this thin design",
  frey_benchmarks = list(beh = 0.46, prop = 0.68, freq = 0.65, general_factor_bifactor = 0.85,
                         source = "restated from P3_spec.md section 4 (Frey et al. 2017 published values)"),
  gates = list(
    g_repro = "beh mean retest in [0.36, 0.56]; prop mean and freq mean each > beh mean and in [0.52, 0.78]; general-factor retest > 0.75; on FAIL stop before Part B",
    g_wall = "DfD (R) and DfE (R and Rexp) each: split-half SB < 0.70 AND 6-month retest < 0.35, and both statistics at least 0.15 below the general-factor retest",
    g_param = "strictly more than half of the 30 gate-family model parameters have retest < 0.60",
    g_nolift = "PASS (E8 null replicated) unless Delta >= 0.05 AND the 90 percent bootstrap CI on Delta excludes 0; a material lift is reported as a real external difference"
  )
)

base_dir <- here::here("data", "raw", "Frey Data")
out_dir <- here::here("output", "elevation", "results", "p3_frey")
SMOKE <- Sys.getenv("P3_SMOKE") == "1"
if (SMOKE) {
  # smoke run (mechanics only): reduced constants, outputs to scratch, never to the results directory
  CFG$n_boot <- 100; CFG$n_perm <- 100; CFG$n_sim_ceiling <- 2
  out_dir <- file.path(Sys.getenv("P3_SMOKE_DIR", tempdir()), "p3_frey")
  cat("SMOKE RUN: reduced constants; outputs to", out_dir, "\n")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sb <- function(r) 2 * r / (1 + r)
spearman <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
n_pairs <- function(x, y) sum(!is.na(x) & !is.na(y))
core6 <- function(id) substr(as.character(id), 1, nchar(as.character(id)) - 2)

# ---------------------------------------------------------------------------
# 1. Guards + load (stop on any structural mismatch)
# ---------------------------------------------------------------------------
stopifnot(system2("git", c("check-ignore", shQuote(file.path("data", "raw", "Frey Data"))),
                  stdout = NULL, stderr = NULL) == 0)
cat("guard: data/raw/Frey Data is gitignored\n")

load_measures <- function(sample) {
  e <- new.env()
  load(file.path(base_dir, sample, "measures", "measures.Rdata"), envir = e)
  meas <- get("measures", envir = e)
  cbn <- names(get("cbalpha", envir = e))
  # 8 shipped measures (FTND, GABS, PG, DAST, NUM, CAREa, CAREs, CAREw) are ordered factors
  # with numeric level labels; convert by label (rank-identical, label-preserving), never silently
  for (v in cbn) {
    if (is.factor(meas[[v]])) {
      conv <- as.numeric(as.character(meas[[v]]))
      stopifnot(sum(is.na(conv)) == sum(is.na(meas[[v]])))
      meas[[v]] <- conv
    }
  }
  list(measures = meas, grouping = get("grouping", envir = e),
       participants = get("participants", envir = e), cbalpha = get("cbalpha", envir = e))
}
main <- load_measures("main")
berl <- load_measures("retest_berlin")

read_frey <- function(sample, ...) read.csv(file.path(base_dir, sample, ...))
dfd_pp_m <- read_frey("main", "dfd", "dfd_perpers.csv")
dfd_prob_m <- read_frey("main", "dfd", "dfd_perprob.csv")
dfe_pp_m <- read_frey("main", "dfe", "dfe_perpers.csv")
dfe_prob_m <- read_frey("main", "dfe", "dfe_perprob.csv")
mp_m <- read_frey("main", "modelpars", "modelpars.csv")
part_m <- read_frey("main", "participants", "participants.csv")
dfd_pp_b <- read_frey("retest_berlin", "dfd", "dfd_perpers.csv")
dfd_prob_b <- read_frey("retest_berlin", "dfd", "dfd_perprob.csv")
dfe_pp_b <- read_frey("retest_berlin", "dfe", "dfe_perpers.csv")
dfe_prob_b <- read_frey("retest_berlin", "dfe", "dfe_perprob.csv")
mp_b <- read_frey("retest_berlin", "modelpars", "modelpars.csv")
part_b <- read_frey("retest_berlin", "participants", "participants.csv")
part_bas <- read_frey("retest_basel", "participants", "participants.csv")

stopifnot(nrow(main$measures) == 1507, length(main$participants) == 1507,
          nrow(part_m) == 1507, nrow(dfd_pp_m) == 1507, nrow(dfe_pp_m) == 1507, nrow(mp_m) == 1507)
stopifnot(nrow(berl$measures) == 109, length(berl$participants) == 109,
          nrow(part_b) == 109, nrow(dfd_pp_b) == 109, nrow(dfe_pp_b) == 109, nrow(mp_b) == 109)
stopifnot(nrow(part_bas) == 64)
stopifnot(identical(colnames(main$measures), colnames(berl$measures)))
stopifnot(identical(as.integer(main$participants), as.integer(dfd_pp_m$partid)))
stopifnot(identical(as.integer(berl$participants), as.integer(dfd_pp_b$partid)))

cb_names <- names(main$cbalpha)
stopifnot(length(cb_names) == 53, all(cb_names %in% colnames(main$measures)))
grp <- main$grouping |>
  mutate(across(c(label, category, group, task), as.character))
gmap <- grp[match(cb_names, grp$label), c("label", "category", "group", "task")]
gmap$label <- cb_names
stopifnot(!anyNA(gmap$group))
risk_set <- cb_names[gmap$category == "risk"]
stopifnot(length(risk_set) == 37, all(gmap$group[gmap$category == "risk"] %in% c("beh", "prop", "freq")))

# linkage (fixed rule; see CFG$linkage_rule)
m_ids <- as.character(main$participants)
b_ids <- as.character(berl$participants)
stopifnot(all(nchar(m_ids) == 8), all(nchar(b_ids) == 8),
          all(substr(m_ids, 7, 8) == "01"), all(substr(b_ids, 7, 8) == "02"))
linked_core <- intersect(core6(b_ids), core6(m_ids))
nonlink_b <- b_ids[!(core6(b_ids) %in% core6(m_ids))]
stopifnot(length(linked_core) == 108, identical(nonlink_b, "65040102"))
bas_link_n <- sum(core6(as.character(part_bas$partid)) %in% core6(m_ids))
cat("linkage: 108 of 109 Berlin retest subjects link; non-linking id 65040102; Basel", bas_link_n, "of 64 link (provenance only)\n")

# identity asserts: representative task indices equal shipped per-person files
idx_dfd <- match(dfd_pp_m$partid, main$participants)
stopifnot(isTRUE(all.equal(main$measures$DFDr[idx_dfd], dfd_pp_m$R)))
idx_dfe <- match(dfe_pp_m$partid, main$participants)
stopifnot(isTRUE(all.equal(main$measures$DFEr[idx_dfe], dfe_pp_m$R)),
          isTRUE(all.equal(main$measures$DFEre[idx_dfe], dfe_pp_m$Rexp)),
          isTRUE(all.equal(main$measures$DFEss[idx_dfe], dfe_pp_m$samples)))
cat("guard: DFDr == dfd R, DFEr == dfe R, DFEre == dfe Rexp, DFEss == dfe samples\n")

# DfD half-assignment recovery: odd/even gamble_ind reproduces Frey's shipped halves exactly
half_means <- function(prob, valcol, min_valid = 1) {
  prob |>
    mutate(half = if_else(gamble_ind %% 2 == 1, "sh1", "sh2"), val = .data[[valcol]]) |>
    summarise(m = if (sum(!is.na(val)) >= min_valid) mean(val, na.rm = TRUE) else NA_real_,
              .by = c(partid, half)) |>
    pivot_wider(names_from = half, values_from = m)
}
hm_m <- half_means(dfd_prob_m, "R")
hm_m <- hm_m[match(dfd_pp_m$partid, hm_m$partid), ]
stopifnot(isTRUE(all.equal(hm_m$sh1, dfd_pp_m$R_sh1)), isTRUE(all.equal(hm_m$sh2, dfd_pp_m$R_sh2)))
hm_b <- half_means(dfd_prob_b, "R")
hm_b <- hm_b[match(dfd_pp_b$partid, hm_b$partid), ]
stopifnot(isTRUE(all.equal(hm_b$sh1, dfd_pp_b$R_sh1)), isTRUE(all.equal(hm_b$sh2, dfd_pp_b$R_sh2)))
cat("guard: odd/even gamble_ind halves reproduce Frey's R_sh1/R_sh2 exactly (main and retest_berlin)\n")

# ---------------------------------------------------------------------------
# 2. Provenance manifest
# ---------------------------------------------------------------------------
used_files <- c(
  file.path("main", c("measures/measures.Rdata", "dfd/dfd_perpers.csv", "dfd/dfd_perprob.csv",
                      "dfe/dfe_perpers.csv", "dfe/dfe_perprob.csv", "modelpars/modelpars.csv",
                      "participants/participants.csv")),
  file.path("retest_berlin", c("measures/measures.Rdata", "dfd/dfd_perpers.csv", "dfd/dfd_perprob.csv",
                               "dfe/dfe_perpers.csv", "dfe/dfe_perprob.csv", "modelpars/modelpars.csv",
                               "participants/participants.csv")),
  "retest_basel/participants/participants.csv",
  "Codebook.pdf"
)
manifest <- tibble(
  file = used_files,
  role = c(rep("read by script (main, N=1507)", 7),
           rep("read by script (retest_berlin, n=109, 6-month)", 7),
           "read by script (retest_basel, n=64, 3-month; provenance count only, not analyzed)",
           "documentation (not read by script)"),
  md5 = unname(tools::md5sum(file.path(base_dir, used_files)))
)
manifest_meta <- tibble(
  file = c("_source", "_linkage", "_n"),
  role = c(paste0("OSF project ", CFG$osf_project, "; local copy at data/raw/Frey Data/ (gitignored, never committed)"),
           paste0(CFG$linkage_rule, "; Berlin links 108 of 109 (non-linking id 65040102); Basel links ", bas_link_n, " of 64"),
           "main N=1507; retest_berlin n=109 (6-month, primary retest); retest_basel n=64 (3-month, secondary, not analyzed here)"),
  md5 = NA_character_
)
readr::write_csv(bind_rows(manifest_meta, manifest), file.path(out_dir, "provenance_manifest.csv"))
cat("provenance manifest written (", nrow(manifest), "files md5-recorded )\n")

# ---------------------------------------------------------------------------
# 3. Part A: per-measure and per-tradition 6-month test-retest + general factor
# ---------------------------------------------------------------------------
m_t1 <- main$measures |> mutate(partid = as.integer(main$participants), core = core6(partid), .before = 1)
m_t2 <- berl$measures |> mutate(partid = as.integer(berl$participants), core = core6(partid), .before = 1)
linked_t1 <- m_t1 |> filter(core %in% linked_core)
linked_t2 <- m_t2 |> filter(core %in% linked_core)
linked_t2 <- linked_t2[match(linked_t1$core, linked_t2$core), ]
stopifnot(nrow(linked_t1) == 108, identical(linked_t1$core, linked_t2$core))

prop_modal <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  max(table(x)) / length(x)
}
alpha_t1 <- unlist(main$cbalpha[1, cb_names])
rel_meas <- lapply(cb_names, function(v) {
  x <- linked_t1[[v]]; y <- linked_t2[[v]]
  n <- n_pairs(x, y)
  tibble(
    measure = v,
    tradition = gmap$group[gmap$label == v],
    category = gmap$category[gmap$label == v],
    task = gmap$task[gmap$label == v],
    cronbach_alpha_t1_shipped = unname(alpha_t1[v]),
    n_retest_pairs = n,
    retest_spearman = if (n >= CFG$retest_min_pairs) spearman(x, y) else NA_real_,
    t1_sd_full = sd(m_t1[[v]], na.rm = TRUE),
    t1_sd_linked = sd(x, na.rm = TRUE),
    t2_sd_linked = sd(y, na.rm = TRUE),
    t1_prop_modal_linked = prop_modal(x),
    flag_degenerate = prop_modal(x) > CFG$degenerate_prop_modal
  )
}) |> bind_rows()

# split-half columns for the DfD/DfE representative indices (internal, same testing occasion; filled in Part B)
rel_meas$splithalf_spearman_sb_t1 <- NA_real_
rel_meas$splithalf_n_t1 <- NA_integer_

# A2: general risk factor (frozen operationalization; see CFG$general_factor)
t1_mu <- sapply(risk_set, function(v) mean(m_t1[[v]], na.rm = TRUE))
t1_sdv <- sapply(risk_set, function(v) sd(m_t1[[v]], na.rm = TRUE))
stopifnot(all(t1_sdv > 0))
z_t1 <- scale(as.matrix(m_t1[, risk_set]), center = t1_mu, scale = t1_sdv)
cmat <- cor(z_t1, use = "pairwise.complete.obs")
stopifnot(!anyNA(cmat))
eig <- eigen(cmat, symmetric = TRUE)
load1 <- eig$vectors[, 1]
names(load1) <- risk_set
dospert_rt <- c("Deth", "Dinv", "Dgam", "Dhea", "Drec", "Dsoc")
if (mean(load1[dospert_rt]) < 0) load1 <- -load1
gf_score <- function(zmat, w, min_avail) {
  apply(zmat, 1, function(z) {
    ok <- !is.na(z)
    if (sum(ok) < min_avail) return(NA_real_)
    sum(w[ok] * z[ok]) / sum(abs(w[ok]))
  })
}
score_t1 <- round(gf_score(z_t1, load1, CFG$general_factor$min_available), CFG$general_factor$score_round_decimals)
z_t2 <- scale(as.matrix(m_t2[, risk_set]), center = t1_mu, scale = t1_sdv)
score_t2 <- round(gf_score(z_t2, load1, CFG$general_factor$min_available), CFG$general_factor$score_round_decimals)
gf_t1_linked <- score_t1[match(linked_t1$partid, m_t1$partid)]
gf_t2_linked <- score_t2[match(linked_t2$partid, m_t2$partid)]
gf_n <- n_pairs(gf_t1_linked, gf_t2_linked)
gf_retest <- spearman(gf_t1_linked, gf_t2_linked)
var_expl <- eig$values[1] / sum(eig$values)
cat(sprintf("general factor: first-PC variance share %.3f; retest Spearman %.4f on n=%d linked\n", var_expl, gf_retest, gf_n))

readr::write_csv(tibble(measure = risk_set, tradition = gmap$group[match(risk_set, gmap$label)],
                        loading_pc1 = load1[risk_set],
                        note = "loadings of the first PC (pairwise-complete Pearson on T1-standardized measures, N=1507); T2 scored with these loadings and T1 standardization"),
                 file.path(out_dir, "general_factor.csv"))

# per-tradition aggregation (risk traditions gate; covar traditions reported for completeness)
trad <- rel_meas |>
  summarise(n_measures = n(), n_with_retest = sum(!is.na(retest_spearman)),
            retest_spearman_mean = mean(retest_spearman, na.rm = TRUE),
            retest_spearman_median = median(retest_spearman, na.rm = TRUE),
            cronbach_alpha_t1_mean = mean(cronbach_alpha_t1_shipped, na.rm = TRUE),
            .by = tradition) |>
  mutate(risk_tradition = tradition %in% c("beh", "prop", "freq"),
         frey_published_retest_mean = case_when(tradition == "beh" ~ CFG$frey_benchmarks$beh,
                                                tradition == "prop" ~ CFG$frey_benchmarks$prop,
                                                tradition == "freq" ~ CFG$frey_benchmarks$freq,
                                                .default = NA_real_)) |>
  arrange(desc(risk_tradition), tradition)
gf_row <- tibble(tradition = "general_factor", n_measures = length(risk_set), n_with_retest = gf_n,
                 retest_spearman_mean = gf_retest, retest_spearman_median = gf_retest,
                 cronbach_alpha_t1_mean = NA_real_, risk_tradition = TRUE,
                 frey_published_retest_mean = CFG$frey_benchmarks$general_factor_bifactor)
trad_out <- bind_rows(trad, gf_row) |>
  mutate(note = if_else(tradition == "general_factor",
                        "simple first-PC general risk factor over the 37 risk measures; Frey benchmark is the published bifactor R, a different (more powerful) aggregate",
                        "mean/median of per-measure 6-month test-retest Spearman over linked subjects"))
readr::write_csv(trad_out, file.path(out_dir, "reliability_by_tradition.csv"))

beh_mean <- trad$retest_spearman_mean[trad$tradition == "beh"]
prop_mean <- trad$retest_spearman_mean[trad$tradition == "prop"]
freq_mean <- trad$retest_spearman_mean[trad$tradition == "freq"]
cat(sprintf("tradition retest means: beh %.4f | prop %.4f | freq %.4f (Frey: 0.46 / 0.68 / 0.65)\n",
            beh_mean, prop_mean, freq_mean))

# G-repro gate (must PASS before Part B is computed)
g_repro <- list(
  beh_mean = beh_mean, prop_mean = prop_mean, freq_mean = freq_mean,
  general_factor_retest = gf_retest,
  beh_in_band = beh_mean >= 0.36 && beh_mean <= 0.56,
  prop_ok = prop_mean > beh_mean && prop_mean >= 0.52 && prop_mean <= 0.78,
  freq_ok = freq_mean > beh_mean && freq_mean >= 0.52 && freq_mean <= 0.78,
  gf_ok = gf_retest > 0.75
)
g_repro$pass <- with(g_repro, beh_in_band && prop_ok && freq_ok && gf_ok)
cat("G-repro:", if (g_repro$pass) "PASS" else "FAIL", "\n")
if (!g_repro$pass) {
  write_json(list(g_repro = g_repro, verdict = "G-repro FAIL: harmonization/linkage suspect; Part B NOT run per the pre-declared gate"),
             file.path(out_dir, "p3_gate.json"), pretty = TRUE, auto_unbox = TRUE, digits = 10)
  readr::write_csv(rel_meas, file.path(out_dir, "reliability_by_measure.csv"))
  stop("G-repro FAILED; stopping before Part B per P3_spec.md section 4")
}

# ---------------------------------------------------------------------------
# 4. Part B1: DfD and DfE single-index reliability (split-half SB + 6-month retest)
# ---------------------------------------------------------------------------
dfd_core_m <- core6(dfd_pp_m$partid); dfd_core_b <- core6(dfd_pp_b$partid)
dfe_core_m <- core6(dfe_pp_m$partid); dfe_core_b <- core6(dfe_pp_b$partid)

dfe_hm_m_R <- half_means(dfe_prob_m, "R", CFG$dfe_min_valid_per_half)
dfe_hm_m_R <- dfe_hm_m_R[match(dfe_pp_m$partid, dfe_hm_m_R$partid), ]
dfe_hm_m_Rexp <- half_means(dfe_prob_m, "Rexp", CFG$dfe_min_valid_per_half)
dfe_hm_m_Rexp <- dfe_hm_m_Rexp[match(dfe_pp_m$partid, dfe_hm_m_Rexp$partid), ]
dfe_hm_b_R <- half_means(dfe_prob_b, "R", CFG$dfe_min_valid_per_half)
dfe_hm_b_R <- dfe_hm_b_R[match(dfe_pp_b$partid, dfe_hm_b_R$partid), ]
dfe_hm_b_Rexp <- half_means(dfe_prob_b, "Rexp", CFG$dfe_min_valid_per_half)
dfe_hm_b_Rexp <- dfe_hm_b_Rexp[match(dfe_pp_b$partid, dfe_hm_b_Rexp$partid), ]

index_row <- function(index, task, sh1, sh2, t1, t1_core, t2, t2_core) {
  r_half <- spearman(sh1, sh2)
  i1 <- t1_core %in% linked_core
  tt1 <- t1[i1][match(linked_core, t1_core[i1])]
  i2 <- t2_core %in% linked_core
  tt2 <- t2[i2][match(linked_core, t2_core[i2])]
  tibble(index = index, task = task,
         splithalf_n = n_pairs(sh1, sh2), splithalf_spearman_raw = r_half, splithalf_spearman_sb = sb(r_half),
         retest_n = n_pairs(tt1, tt2), retest_spearman = spearman(tt1, tt2))
}
dfd_dfe <- bind_rows(
  index_row("dfd_R (= DFDr, proportion risky, description)", "DfD",
            dfd_pp_m$R_sh1, dfd_pp_m$R_sh2, dfd_pp_m$R, dfd_core_m, dfd_pp_b$R, dfd_core_b),
  index_row("dfe_R (= DFEr, proportion risky, experience)", "DfE",
            dfe_hm_m_R$sh1, dfe_hm_m_R$sh2, dfe_pp_m$R, dfe_core_m, dfe_pp_b$R, dfe_core_b),
  index_row("dfe_Rexp (= DFEre, proportion risky by experienced values; cbalpha-representative)", "DfE",
            dfe_hm_m_Rexp$sh1, dfe_hm_m_Rexp$sh2, dfe_pp_m$Rexp, dfe_core_m, dfe_pp_b$Rexp, dfe_core_b)
) |>
  mutate(general_factor_retest = gf_retest,
         note = "splithalf = internal consistency at T1 (SB-corrected Spearman over odd/even problem halves); retest = 6-month temporal stability (Spearman, linked subjects); the two estimands are distinct and never conflated")
readr::write_csv(dfd_dfe, file.path(out_dir, "dfd_dfe_reliability.csv"))
print(as.data.frame(dfd_dfe[, 1:7]), row.names = FALSE)

# backfill split-half columns into reliability_by_measure for the representative indices
fill_sh <- function(meas, row_i) {
  rel_meas$splithalf_spearman_sb_t1[rel_meas$measure == meas] <<- dfd_dfe$splithalf_spearman_sb[row_i]
  rel_meas$splithalf_n_t1[rel_meas$measure == meas] <<- dfd_dfe$splithalf_n[row_i]
}
fill_sh("DFDr", 1); fill_sh("DFEr", 2); fill_sh("DFEre", 3)
readr::write_csv(rel_meas, file.path(out_dir, "reliability_by_measure.csv"))

# ---------------------------------------------------------------------------
# 5. Part B2: 6-month test-retest of model-based parameters
# ---------------------------------------------------------------------------
par_cols <- setdiff(names(mp_m), "partid")
par_cols <- par_cols[!grepl("_(DEV|BIC)$", par_cols)]
gate_suffix <- c("OS", "OSg", "OSl", "PW", "PWg", "PWl", "LA")
mp_core_m <- core6(mp_m$partid); mp_core_b <- core6(mp_b$partid)
i1 <- mp_core_m %in% linked_core
i2 <- mp_core_b %in% linked_core
mp1 <- mp_m[i1, ][match(linked_core, mp_core_m[i1]), ]
mp2 <- mp_b[i2, ][match(linked_core, mp_core_b[i2]), ]
modelpar_rel <- lapply(par_cols, function(v) {
  prefix <- sub("_.*$", "", v)
  task <- sub("(EU|3P|4P|6P|RN|1|3)$", "", prefix)
  model <- substring(prefix, nchar(task) + 1)
  suffix <- sub("^.*_", "", v)
  tibble(parameter = v, task = task, model = model, param = suffix,
         param_class = case_when(suffix %in% c("OS", "OSg", "OSl") ~ "outcome sensitivity",
                                 suffix %in% c("PW", "PWg", "PWl") ~ "probability weighting",
                                 suffix == "LA" ~ "loss aversion",
                                 suffix == "CS" ~ "choice sensitivity",
                                 .default = "other"),
         gate_param = task %in% c("ABCT", "MPL", "CCT") & suffix %in% gate_suffix,
         n_retest_pairs = n_pairs(mp1[[v]], mp2[[v]]),
         retest_spearman = if (n_pairs(mp1[[v]], mp2[[v]]) >= CFG$retest_min_pairs) spearman(mp1[[v]], mp2[[v]]) else NA_real_)
}) |> bind_rows()
readr::write_csv(modelpar_rel, file.path(out_dir, "modelpar_reliability.csv"))
gate_pars <- modelpar_rel |> filter(gate_param)
cat(sprintf("B2: %d gate-family parameters (ABCT/MPL/CCT x OS/PW/LA family); %d of %d with retest < 0.60; median retest %.3f\n",
            nrow(gate_pars), sum(gate_pars$retest_spearman < 0.60, na.rm = TRUE), nrow(gate_pars),
            median(gate_pars$retest_spearman, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 6. Part B3: E8 (step 26) parallel - hierarchical partial pooling of the DfD index (T1, N=1507)
# ---------------------------------------------------------------------------
dfd_prob_m <- dfd_prob_m |> mutate(half = if_else(gamble_ind %% 2 == 1, "odd", "even"))
fit_half <- function(dat) {
  ctrl <- glmerControl(optCtrl = list(maxfun = 1e5))
  fit <- suppressMessages(glmer(R ~ 1 + domain + (1 | partid), data = dat, family = binomial, control = ctrl))
  conv <- length(fit@optinfo$conv$lme4$messages) == 0
  if (!conv) {
    ctrl <- glmerControl(optCtrl = list(maxfun = 1e6))
    fit <- suppressMessages(glmer(R ~ 1 + domain + (1 | partid), data = dat, family = binomial, control = ctrl))
    conv <- length(fit@optinfo$conv$lme4$messages) == 0
  }
  list(fit = fit, converged = conv, singular = isSingular(fit))
}
blups <- function(fit, ids) {
  re <- ranef(fit)$partid
  round(re[match(as.character(ids), rownames(re)), "(Intercept)"], CFG$blup_round_decimals)
}
h_odd <- fit_half(dfd_prob_m |> filter(half == "odd"))
h_even <- fit_half(dfd_prob_m |> filter(half == "even"))
blup_odd <- blups(h_odd$fit, dfd_pp_m$partid)
blup_even <- blups(h_even$fit, dfd_pp_m$partid)
stopifnot(!anyNA(blup_odd), !anyNA(blup_even))

r_raw <- spearman(dfd_pp_m$R_sh1, dfd_pp_m$R_sh2)
r_hier <- spearman(blup_odd, blup_even)
rel_raw <- sb(r_raw); rel_hier <- sb(r_hier)
delta <- rel_hier - rel_raw
cat(sprintf("B3: rel_raw(SB) %.4f | rel_hier(SB) %.4f | Delta %+.4f (N=%d)\n", rel_raw, rel_hier, delta, nrow(dfd_pp_m)))

# subject-level bootstrap over frozen per-subject half estimates (no refit; B=2000)
set.seed(CFG$seed_boot)
N <- nrow(dfd_pp_m)
boot <- replicate(CFG$n_boot, {
  ii <- sample.int(N, N, replace = TRUE)
  rr <- sb(spearman(dfd_pp_m$R_sh1[ii], dfd_pp_m$R_sh2[ii]))
  rh <- sb(spearman(blup_odd[ii], blup_even[ii]))
  c(raw = rr, hier = rh, delta = rh - rr)
})
qs <- function(x, lv) quantile(x, probs = c((1 - lv) / 2, 1 - (1 - lv) / 2), names = FALSE, na.rm = TRUE)
ci_d90 <- qs(boot["delta", ], 0.90); ci_d95 <- qs(boot["delta", ], 0.95)
ci_r90 <- qs(boot["raw", ], 0.90); ci_h90 <- qs(boot["hier", ], 0.90)

# null permutation: break the within-subject half link (same shuffle both arms)
set.seed(CFG$seed_perm)
perm <- replicate(CFG$n_perm, {
  jj <- sample.int(N, N, replace = FALSE)
  c(raw = spearman(dfd_pp_m$R_sh1, dfd_pp_m$R_sh2[jj]), hier = spearman(blup_odd, blup_even[jj]))
})
perm_stats <- function(v) {
  z <- mean(v) / (sd(v) / sqrt(length(v)))
  c(mean = mean(v), sd = sd(v), z = z, flag = as.numeric(abs(mean(v)) > 0.02 && abs(z) > 3))
}
perm_raw <- perm_stats(perm["raw", ]); perm_hier <- perm_stats(perm["hier", ])

# ceiling sanity (SIMULATION, labeled): model-implied attainable reliability at 4+4 binary trials
h_full <- fit_half(dfd_prob_m)
sd_hat <- sqrt(unname(VarCorr(h_full$fit)$partid[1, 1]))
fe <- fixef(h_full$fit)
design <- dfd_prob_m |> select(partid, gamble_ind, domain, half)
set.seed(CFG$seed_sim)
sim_res <- vapply(seq_len(CFG$n_sim_ceiling), function(s) {
  u <- rnorm(N, 0, sd_hat)
  names(u) <- as.character(dfd_pp_m$partid)
  eta <- fe[["(Intercept)"]] + fe[["domainloss"]] * (design$domain == "loss") + u[as.character(design$partid)]
  ysim <- rbinom(nrow(design), 1, plogis(eta))
  dsim <- design |> mutate(R = ysim)
  raw1 <- tapply(dsim$R[dsim$half == "odd"], dsim$partid[dsim$half == "odd"], mean)
  raw2 <- tapply(dsim$R[dsim$half == "even"], dsim$partid[dsim$half == "even"], mean)
  raw2 <- raw2[names(raw1)]
  f1 <- fit_half(dsim |> filter(half == "odd")); f2 <- fit_half(dsim |> filter(half == "even"))
  ids <- names(raw1)
  bs1 <- blups(f1$fit, ids); bs2 <- blups(f2$fit, ids)
  c(raw = sb(spearman(raw1, raw2)), hier = sb(spearman(bs1, bs2)),
    conv = as.numeric(f1$converged && f2$converged))
}, numeric(3))
sim_ok <- sim_res["conv", ] == 1
stopifnot(sum(sim_ok) >= 1)
cat(sprintf("ceiling SIMULATION: %d/%d sims converged; raw %.4f [%.4f, %.4f]; hier %.4f [%.4f, %.4f]\n",
            sum(sim_ok), CFG$n_sim_ceiling,
            mean(sim_res["raw", sim_ok]), min(sim_res["raw", sim_ok]), max(sim_res["raw", sim_ok]),
            mean(sim_res["hier", sim_ok]), min(sim_res["hier", sim_ok]), max(sim_res["hier", sim_ok])))

hier_out <- tribble(
  ~metric, ~value,
  "n_subjects", N,
  "n_problems_per_half", 4,
  "splithalf_spearman_raw_uncorrected", r_raw,
  "splithalf_rel_raw_sb", rel_raw,
  "splithalf_spearman_hier_uncorrected", r_hier,
  "splithalf_rel_hier_sb", rel_hier,
  "delta_hier_minus_raw", delta,
  "delta_ci90_lo", ci_d90[1], "delta_ci90_hi", ci_d90[2],
  "delta_ci95_lo", ci_d95[1], "delta_ci95_hi", ci_d95[2],
  "rel_raw_ci90_lo", ci_r90[1], "rel_raw_ci90_hi", ci_r90[2],
  "rel_hier_ci90_lo", ci_h90[1], "rel_hier_ci90_hi", ci_h90[2],
  "n_boot", CFG$n_boot,
  "perm_raw_mean_r", perm_raw[["mean"]], "perm_raw_z", perm_raw[["z"]], "perm_raw_flag", perm_raw[["flag"]],
  "perm_hier_mean_r", perm_hier[["mean"]], "perm_hier_z", perm_hier[["z"]], "perm_hier_flag", perm_hier[["flag"]],
  "n_perm", CFG$n_perm,
  "rankagree_blup_vs_raw_odd", spearman(blup_odd, dfd_pp_m$R_sh1),
  "rankagree_blup_vs_raw_even", spearman(blup_even, dfd_pp_m$R_sh2),
  "n_distinct_raw_half_values", length(unique(dfd_pp_m$R_sh1)),
  "n_distinct_blup_odd_values", length(unique(blup_odd)),
  "glmer_converged_odd", as.numeric(h_odd$converged), "glmer_converged_even", as.numeric(h_even$converged),
  "glmer_singular_odd", as.numeric(h_odd$singular), "glmer_singular_even", as.numeric(h_even$singular),
  "full_fit_random_intercept_sd", sd_hat,
  "full_fit_intercept", unname(fe[["(Intercept)"]]),
  "full_fit_domainloss", unname(fe[["domainloss"]]),
  "sim_ceiling_n_used", sum(sim_ok),
  "sim_ceiling_raw_mean", mean(sim_res["raw", sim_ok]),
  "sim_ceiling_raw_min", min(sim_res["raw", sim_ok]), "sim_ceiling_raw_max", max(sim_res["raw", sim_ok]),
  "sim_ceiling_hier_mean", mean(sim_res["hier", sim_ok]),
  "sim_ceiling_hier_min", min(sim_res["hier", sim_ok]), "sim_ceiling_hier_max", max(sim_res["hier", sim_ok]),
  "sim_ceiling_delta_mean", mean(sim_res["hier", sim_ok] - sim_res["raw", sim_ok])
) |>
  mutate(note = "DfD T1 split-half in SB units; halves = odd/even gamble_ind (Frey's own R_sh1/R_sh2 rule); sim_ceiling rows are SIMULATION (model-implied attainable reliability at this thin 8-problem design), not empirical estimates")
readr::write_csv(hier_out, file.path(out_dir, "dfd_hier_reestimation.csv"))

# ---------------------------------------------------------------------------
# 7. Gates
# ---------------------------------------------------------------------------
wall_rows <- dfd_dfe |>
  mutate(sh_low = splithalf_spearman_sb < 0.70,
         rt_low = retest_spearman < 0.35,
         materially_below_gf = (gf_retest - retest_spearman) >= CFG$material_below_margin &
                               (gf_retest - splithalf_spearman_sb) >= CFG$material_below_margin)
g_wall <- list(
  per_index = wall_rows |> select(index, splithalf_spearman_sb, retest_spearman, sh_low, rt_low, materially_below_gf),
  pass = all(wall_rows$sh_low & wall_rows$rt_low & wall_rows$materially_below_gf)
)
g_param <- list(
  n_gate_params = nrow(gate_pars),
  n_below_060 = sum(gate_pars$retest_spearman < 0.60, na.rm = TRUE),
  frac_below_060 = mean(gate_pars$retest_spearman < 0.60, na.rm = TRUE),
  median_retest = median(gate_pars$retest_spearman, na.rm = TRUE),
  max_retest = max(gate_pars$retest_spearman, na.rm = TRUE),
  pass = mean(gate_pars$retest_spearman < 0.60, na.rm = TRUE) > 0.5
)
lift_material <- (delta >= 0.05) && (ci_d90[1] > 0)
g_nolift <- list(
  delta = delta, delta_ci90 = ci_d90, delta_ci95 = ci_d95,
  material_lift_observed = lift_material,
  pass = !lift_material,
  note = "PASS means the E8 null replicates externally (no material lift); a material lift would be reported as a real difference from E8"
)
gates <- list(
  g_repro = g_repro,
  g_wall = g_wall,
  g_param = g_param,
  g_nolift = g_nolift,
  b4_contrast = list(general_factor_retest = gf_retest,
                     beh_mean_retest = beh_mean,
                     dfd_R_retest = dfd_dfe$retest_spearman[1],
                     dfe_R_retest = dfd_dfe$retest_spearman[2],
                     dfe_Rexp_retest = dfd_dfe$retest_spearman[3],
                     note = "aggregate reliable vs single measures noisy: Study 1's decomposition, external"),
  verdict = list(
    g_repro = if (g_repro$pass) "PASS" else "FAIL",
    g_wall = if (g_wall$pass) "PASS" else "FAIL",
    g_param = if (g_param$pass) "PASS" else "FAIL",
    g_nolift = if (g_nolift$pass) "PASS (E8 null replicated)" else "FAIL (material lift; real external difference from E8)"
  )
)
write_json(gates, file.path(out_dir, "p3_gate.json"), pretty = TRUE, auto_unbox = TRUE, digits = 10)
cat("gate verdicts:", paste(names(gates$verdict), unlist(gates$verdict), sep = "=", collapse = " | "), "\n")

# ---------------------------------------------------------------------------
# 8. Per-subject export (every reliability can be re-derived from this file)
# ---------------------------------------------------------------------------
core_all <- core6(m_t1$partid)
linked <- core_all %in% linked_core
t2_at <- function(df_t2, core_t2, col) {
  out <- rep(NA_real_, length(core_all))
  hit <- match(core_all, core_t2)
  out[!is.na(hit)] <- df_t2[[col]][hit[!is.na(hit)]]
  out
}
exp_df <- tibble(partid_main = m_t1$partid, core6 = core_all, linked_berlin = linked)
for (v in cb_names) {
  exp_df[[paste0(v, "_t1")]] <- m_t1[[v]]
  exp_df[[paste0(v, "_t2")]] <- t2_at(m_t2, m_t2$core, v)
}
exp_df <- exp_df |>
  mutate(dfd_R_t1 = dfd_pp_m$R[match(partid_main, dfd_pp_m$partid)],
         dfd_R_sh1_t1 = dfd_pp_m$R_sh1[match(partid_main, dfd_pp_m$partid)],
         dfd_R_sh2_t1 = dfd_pp_m$R_sh2[match(partid_main, dfd_pp_m$partid)],
         dfd_R_t2 = t2_at(dfd_pp_b, core6(dfd_pp_b$partid), "R"),
         dfd_R_sh1_t2 = t2_at(dfd_pp_b, core6(dfd_pp_b$partid), "R_sh1"),
         dfd_R_sh2_t2 = t2_at(dfd_pp_b, core6(dfd_pp_b$partid), "R_sh2"),
         dfe_R_t1 = dfe_pp_m$R[match(partid_main, dfe_pp_m$partid)],
         dfe_R_sh1_t1 = dfe_hm_m_R$sh1[match(partid_main, dfe_hm_m_R$partid)],
         dfe_R_sh2_t1 = dfe_hm_m_R$sh2[match(partid_main, dfe_hm_m_R$partid)],
         dfe_Rexp_t1 = dfe_pp_m$Rexp[match(partid_main, dfe_pp_m$partid)],
         dfe_Rexp_sh1_t1 = dfe_hm_m_Rexp$sh1[match(partid_main, dfe_hm_m_Rexp$partid)],
         dfe_Rexp_sh2_t1 = dfe_hm_m_Rexp$sh2[match(partid_main, dfe_hm_m_Rexp$partid)],
         dfe_R_t2 = t2_at(dfe_pp_b, core6(dfe_pp_b$partid), "R"),
         dfe_Rexp_t2 = t2_at(dfe_pp_b, core6(dfe_pp_b$partid), "Rexp"),
         hier_blup_odd_t1 = blup_odd[match(partid_main, dfd_pp_m$partid)],
         hier_blup_even_t1 = blup_even[match(partid_main, dfd_pp_m$partid)],
         gf_score_t1 = score_t1[match(partid_main, m_t1$partid)],
         gf_score_t2 = t2_at(tibble(v = score_t2), core6(m_t2$partid), "v"))
for (v in par_cols) {
  exp_df[[paste0("par_", v, "_t1")]] <- mp_m[[v]][match(exp_df$partid_main, mp_m$partid)]
  exp_df[[paste0("par_", v, "_t2")]] <- t2_at(mp_b, core6(mp_b$partid), v)
}
readr::write_csv(exp_df, file.path(out_dir, "per_subject_export.csv"))
cat("per-subject export:", nrow(exp_df), "rows x", ncol(exp_df), "cols\n")

# ---------------------------------------------------------------------------
# 9. Config + in-run re-derivation from the written export
# ---------------------------------------------------------------------------
CFG$r_version <- R.version.string
CFG$lme4_version <- as.character(packageVersion("lme4"))
write_json(CFG, file.path(out_dir, "p3_config.json"), pretty = TRUE, auto_unbox = TRUE, digits = 10)

chk <- utils::read.csv(file.path(out_dir, "per_subject_export.csv"))
tol <- 1e-9
rederive_ok <- c(
  gf = abs(spearman(chk$gf_score_t1[chk$linked_berlin], chk$gf_score_t2[chk$linked_berlin]) - gf_retest) < tol,
  dfd_sh = abs(sb(spearman(chk$dfd_R_sh1_t1, chk$dfd_R_sh2_t1)) - dfd_dfe$splithalf_spearman_sb[1]) < tol,
  dfd_rt = abs(spearman(chk$dfd_R_t1[chk$linked_berlin], chk$dfd_R_t2[chk$linked_berlin]) - dfd_dfe$retest_spearman[1]) < tol,
  dfe_sh_Rexp = abs(sb(spearman(chk$dfe_Rexp_sh1_t1, chk$dfe_Rexp_sh2_t1)) - dfd_dfe$splithalf_spearman_sb[3]) < tol,
  hier = abs(sb(spearman(chk$hier_blup_odd_t1, chk$hier_blup_even_t1)) - rel_hier) < tol,
  beh_mean = abs(mean(sapply(gmap$label[gmap$group == "beh"], function(v)
    spearman(chk[[paste0(v, "_t1")]][chk$linked_berlin], chk[[paste0(v, "_t2")]][chk$linked_berlin]))) - beh_mean) < tol,
  par_example = abs(spearman(chk$par_ABCT4P_LA_t1[chk$linked_berlin], chk$par_ABCT4P_LA_t2[chk$linked_berlin]) -
    modelpar_rel$retest_spearman[modelpar_rel$parameter == "ABCT4P_LA"]) < tol
)
cat("re-derivation from committed export:", paste(names(rederive_ok), rederive_ok, sep = "=", collapse = " "), "\n")
stopifnot(all(rederive_ok))

cat(sprintf("DONE p3_frey_replication in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

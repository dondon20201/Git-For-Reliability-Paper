# Step 13: Figure 1: reliability curves for both problem pools
# Reads:   data/prior_studies.csv; output/results/step2/, step3/, step4/ (steps 05, 08, 10, 11)
# Writes:  output/results/step5/ (fig_reconciliation_v2_data.csv, fig_reconciliation_v2.png/.pdf, caption)
# Paper:   Results, Figure 1; the grid maxima quoted in the text
# Runtime: a few seconds

# Reconciliation figure v2: BOTH pools.
#
# ZERO estimation: every curve value is READ, string-verbatim, from the two
# saved reliability CSVs (step 05 for the clean pool; step 11 for the expanded
# pool). Every prior-study fact is READ from data/prior_studies.csv.
# The only numeric operations are (a) total problems = 2 * k for x-axis
# coordinates and (b) rounding for labels. Fully deterministic — no RNG.
# The figure data CSV carries the curve values string-verbatim (asserted) and
# the caption's key phrases are machine-checked.
#
# The step 08 figure (output/results/step3/fig_reconciliation.*) is RETAINED
# (supplement, clean-pool-only) — asserted below.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(here)
})

res_dir <- here("output", "results", "step5")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

clean_path   <- here("output", "results", "step2", "reliability_by_cell_construct.csv")
ext_path     <- here("output", "results", "step4", "reliability_extended.csv")
needed_path  <- here("output", "results", "step2", "trials_needed_table.csv")
interp_path  <- here("output", "results", "step4", "extended_interpretation.csv")
pool_path    <- here("output", "results", "step4", "expanded_pool_summary.csv")
studies_path <- here("data", "prior_studies.csv")

# Read EVERYTHING as character so the saved values flow into the figure-data
# CSV string-verbatim (numeric copies are made only for plotting coordinates).
clean_chr <- read_csv(clean_path,   col_types = cols(.default = col_character()))
ext_chr   <- read_csv(ext_path,     col_types = cols(.default = col_character()))
needed    <- read_csv(needed_path,  col_types = cols(.default = col_character()))
interp    <- read_csv(interp_path,  col_types = cols(.default = col_character()))
pool_sum  <- read_csv(pool_path,    col_types = cols(.default = col_character()))
studies   <- read_csv(studies_path, col_types = cols(.default = col_character()))

# ---- Input assertions (fail loudly if inputs are not the expected ones)
stopifnot(nrow(studies) == 7L)
expected_in_figure <- c("mkrtchian2023_loss", "mkrtchian2023_risk", "walasek_stewart2021")
stopifnot(identical(
  sort(studies$study_id[studies$in_figure == "yes"]),
  sort(expected_in_figure)
))
stopifnot(pool_sum$n_problems_clean_step2 == "105",
          pool_sum$n_problems_expanded == "182")

constructs_headline <- c("ev_sens_agn", "loss_avoid_adj")
constructs_grey     <- c("b_ev", "b_risk", "b_loss")

curve_clean <- clean_chr |>
  filter(axis == "i_exp", regime == "exp") |>
  mutate(pool = "clean") |>
  select(pool, cell_id, construct, k, r_sb, ci_lo_sb, ci_hi_sb, n_subjects, status)

curve_ext <- ext_chr |>
  filter(axis == "i_exp_ext", regime == "exp") |>
  mutate(pool = "expanded") |>
  select(pool, cell_id, construct, k, r_sb, ci_lo_sb, ci_hi_sb, n_subjects, status)

stopifnot(identical(sort(unique(as.integer(curve_clean$k))), 2:7))
stopifnot(identical(sort(unique(as.integer(curve_ext$k))), 2:13))

plotted_chr <- bind_rows(curve_clean, curve_ext) |> filter(status == "VERIFIED")

# Expected VERIFIED counts: clean — headline 6 k-points x 2, grey (k>=4) 4 x 3
# = 24; expanded — headline 12 x 2, grey (k>=4) 10 x 3 = 54. Total 78.
counts <- plotted_chr |> count(pool, construct)
stopifnot(
  all(counts$n[counts$pool == "clean"    & counts$construct %in% constructs_headline] == 6L),
  all(counts$n[counts$pool == "clean"    & counts$construct %in% constructs_grey]     == 4L),
  all(counts$n[counts$pool == "expanded" & counts$construct %in% constructs_headline] == 12L),
  all(counts$n[counts$pool == "expanded" & counts$construct %in% constructs_grey]     == 10L),
  nrow(plotted_chr) == 78L,
  !any(is.na(as.numeric(plotted_chr$r_sb)))
)

# Numeric copies for plotting coordinates only.
plotted <- plotted_chr |>
  mutate(
    k              = as.integer(k),
    total_problems = 2L * k,     # x-axis unit: total problems (2 x k per split half)
    r_num          = as.numeric(r_sb),
    lo_num         = as.numeric(ci_lo_sb),
    hi_num         = as.numeric(ci_hi_sb),
    layer          = if_else(construct %in% constructs_headline, "headline", "cautionary_grey")
  )

# Observed grid ends at 26 total problems (expanded pool, k = 13).
# HARD assertion: no plotted x beyond 26.
grid_end <- 26L
stopifnot(max(plotted$total_problems) == grid_end,
          max(plotted$total_problems[plotted$pool == "clean"]) == 14L,
          all(plotted$total_problems <= grid_end))

# ---- Crossing values for the caption (READ from the saved CSVs, never typed)
row_cross_clean <- needed |>
  filter(axis == "i_exp", construct == "loss_avoid_adj", threshold == "0.5")
stopifnot(row_cross_clean$entry == "reached")
cross_clean_total <- row_cross_clean$total_problems            # "12"

row_cross_ext <- interp |>
  filter(question == "loss_avoid_adj crosses SB .5 at k (problems/half)")
stopifnot(nrow(row_cross_ext) == 1L)
cross_ext_k     <- as.integer(row_cross_ext$answer)            # 12
cross_ext_total <- as.character(2L * cross_ext_k)              # "24"

# Cross-check the pre-declared crossings against the curve CSVs themselves.
first_k_at <- function(df, constr, thr) {
  d <- df |> filter(construct == constr) |>
    mutate(k = as.integer(k), r = as.numeric(r_sb)) |> arrange(k)
  ks <- d$k[!is.na(d$r) & d$r >= thr]
  if (length(ks) == 0) NA_integer_ else min(ks)
}
stopifnot(first_k_at(curve_clean |> filter(status == "VERIFIED"), "loss_avoid_adj", 0.5) ==
            as.integer(row_cross_clean$budget))
stopifnot(first_k_at(curve_ext |> filter(status == "VERIFIED"), "loss_avoid_adj", 0.5) ==
            cross_ext_k)
# The wall: nothing reaches SB .7 on either pool anywhere on the observed grid.
for (co in c(constructs_headline, constructs_grey)) {
  stopifnot(is.na(first_k_at(plotted_chr |> filter(pool == "clean"), co, 0.7)))
  stopifnot(is.na(first_k_at(plotted_chr |> filter(pool == "expanded"), co, 0.7)))
}

# ---- Study placements (x-axis only — no y-position on our curves)
fig_studies <- studies |> filter(in_figure == "yes")
mk <- fig_studies |> filter(str_starts(study_id, "mkrtchian"))
ws <- fig_studies |> filter(study_id == "walasek_stewart2021")
stopifnot(all(mk$design_points_equiv == "120"), ws$design_points_equiv == "64-100")

mk_loss <- mk |> filter(study_id == "mkrtchian2023_loss")
mk_risk <- mk |> filter(study_id == "mkrtchian2023_risk")

# The prior-study overlays (Mkrtchian marker, Walasek band) are not drawn on
# this figure. The mk_loss / mk_risk / ws reads above are kept because the
# study rows stay in fig_reconciliation_v2_data.csv and the caption draft
# references them.
lab_grey <- paste0(
  "Firth b_ev / b_risk / b_loss: near-0 at clean-pool budgets;\n",
  "b_ev climbs on the expanded pool (dashed) but stays below .7;\n",
  "k < 4 per half structurally unidentified"
)

lab_mix <- paste0(
  "problem-mix effect: the loss-avoidance index crosses r = .5 at ",
  cross_clean_total, " problems (clean, solid)\nbut ",
  cross_ext_total, " problems (expanded, dashed)"
)

# ---- Figure
x_max   <- 160
col_map <- c(ev_sens_agn = "#0072B2", loss_avoid_adj = "#D55E00")
lbl_map <- c(ev_sens_agn = "EV sensitivity (model-agnostic index, SB)",
             loss_avoid_adj = "Loss avoidance, EV-adjusted (SB)")
n_clean <- pool_sum$n_problems_clean_step2
n_ext   <- pool_sum$n_problems_expanded
lt_map  <- c(clean = "solid", expanded = "42")
lt_lbl  <- c(clean    = paste0("clean pool (", n_clean, " two-outcome, no-ambiguity problems)"),
             expanded = paste0("expanded pool (", n_ext, " reconstructed Amb==0 problems)"))

headline  <- plotted |> filter(layer == "headline")
greylayer <- plotted |> filter(layer == "cautionary_grey")

p <- ggplot() +
  # Pre-declared reliability thresholds
  geom_hline(yintercept = c(0.5, 0.7), linetype = "dotted", colour = "grey45",
             linewidth = 0.4) +
  # Headline curves: CI ribbons + straight segments between observed points
  # only (no smoothing, no extrapolation); pools distinguished by linetype
  geom_ribbon(data = headline,
              aes(x = total_problems, ymin = lo_num, ymax = hi_num,
                  fill = construct, group = interaction(construct, pool)),
              alpha = 0.13, colour = NA) +
  geom_line(data = headline,
            aes(total_problems, r_num, colour = construct, linetype = pool),
            linewidth = 0.8) +
  geom_point(data = headline,
             aes(total_problems, r_num, colour = construct,
                 shape = pool), size = 1.7) +
  # Cautionary grey layer: Firth model parameters, both pools
  geom_line(data = greylayer,
            aes(total_problems, r_num, linetype = pool,
                group = interaction(construct, pool)),
            colour = "grey55", linewidth = 0.45) +
  geom_point(data = greylayer,
             aes(total_problems, r_num, shape = pool,
                 group = interaction(construct, pool)),
             colour = "grey55", size = 1.0) +
  # Annotations (curves only)
  annotate("text", x = 9, y = -0.09, label = lab_grey, size = 2.7,
           colour = "grey35", lineheight = 1.05) +
  annotate("text", x = 10.5, y = 0.88, label = lab_mix, size = 2.8,
           colour = "grey10", lineheight = 1.08) +
  annotate("text", x = 4.15, y = c(0.515, 0.715), hjust = 0, vjust = -0.4,
           label = c("r = .5", "r = .7"), size = 2.7, colour = "grey45") +
  scale_x_log10(breaks = c(4, 6, 8, 10, 14, 20, 26),
                limits = c(4, grid_end)) +
  scale_colour_manual(values = col_map, labels = lbl_map, name = NULL) +
  scale_fill_manual(values = col_map, labels = lbl_map, name = NULL) +
  scale_linetype_manual(values = lt_map, labels = lt_lbl, name = NULL) +
  scale_shape_manual(values = c(clean = 16, expanded = 17), labels = lt_lbl,
                     name = NULL) +
  coord_cartesian(ylim = c(-0.15, 1.0)) +
  # Explicit legend order: without it ggplot2 stacks the two legend blocks in
  # a run-dependent order. Pool key (linetype + shape) first, construct key
  # (colour + fill) second. Merged guides must share an order or they split.
  guides(linetype = guide_legend(nrow = 2, order = 1),
         shape = guide_legend(nrow = 2, order = 1),
         colour = guide_legend(nrow = 2, order = 2),
         fill = guide_legend(order = 2)) +
  labs(
    x = "Total problems (2 × k per split half; log scale)",
    y = "Split-half reliability (Spearman rank, Spearman–Brown corrected)",
    title = "CPC18 split-half reliability versus distinct-problem count in the clean and expanded pools",
    subtitle = paste0(
      "Experience regime, r = 20 trials/problem; ribbons = 95% subject-cluster ",
      "bootstrap CIs. Solid = clean pool, dashed = expanded pool."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 8.5, colour = "grey30")
  )

ggsave(file.path(res_dir, "fig_reconciliation_v2.png"), p,
       width = 10, height = 7, dpi = 300)
ggsave(file.path(res_dir, "fig_reconciliation_v2.pdf"), p,
       width = 10, height = 7, device = cairo_pdf)

# ---- Caption (draft); machine-checked below
caption_text <- paste0(
  "Figure. Prior studies mapped onto the CPC18 reliability-versus-problem-count ",
  "curves for BOTH problem pools (experience regime, r = 20 trials per problem; ",
  "split-half Spearman with Spearman–Brown correction; ribbons are 95% ",
  "subject-cluster bootstrap CIs). Solid curves: the clean pool (", n_clean,
  " two-outcome, no-ambiguity problems; observed grid 4–14 total problems). ",
  "Dashed curves: the expanded pool (", n_ext, " Amb==0 problems, multi-outcome ",
  "lotteries reconstructed and validated exactly against the official CPC18 ",
  "implementation; observed grid 4–", grid_end, " total problems). Coloured ",
  "curves: EV-adjusted loss avoidance (loss_avoid_adj) and the model-agnostic ",
  "EV-sensitivity index (ev_sens_agn); grey curves: the Firth model parameters ",
  "(b_ev, b_risk, b_loss), structurally unidentified below 4 problems per half. ",
  "All curves are drawn only over their observed grids, with straight segments ",
  "between observed points — no smoothing and no extrapolation; the shaded ",
  "region beyond ", grid_end, " problems contains no curve values from either ",
  "pool. No construct reaches r = .7 on either pool at any observed budget. ",
  "The problem-mix effect is visible within the figure itself: loss_avoid_adj ",
  "crosses r = .5 at ", cross_clean_total, " total problems on the clean pool ",
  "but only at ", cross_ext_total, " total problems on the expanded pool — the ",
  "same construct, estimator, and subjects, at matched problem counts, yield ",
  "systematically different reliability under a different problem mix. ",
  "Prior studies are placed on the x-axis only, at their ",
  "unique-informative-design-points equivalent, and are given no position on ",
  "our curves: Mkrtchian et al. (2023; ~120 individually calibrated offers, ",
  "120 modelled trials; loss aversion test–retest ICC(A,1) = ",
  mk_loss$reliability_value, " [", mk_loss$ci_95, "], risk ICC(A,1) = ",
  mk_risk$reliability_value, " [", mk_risk$ci_95, "]) and Walasek & Stewart ",
  "(2021; 64–100 unique gambles per set; simulation-based “extremely poor” ",
  "recoverability of λ). The mapping is suggestive, not probative: problem ",
  "count appears necessary but not sufficient. Walasek & Stewart’s 64–100-gamble ",
  "sets exceed both pools’ observed grids and λ still fails — via a λ–bias ",
  "tradeoff tied to task format and via choice-set dependence tied to problem ",
  "mix (the same dependence our two pools display directly). Mkrtchian et ",
  "al.’s high reliability sits at both a larger calibrated design AND a ",
  "different task format; we do not adjudicate between these. ",
  "Estimation-regime caveat: Mkrtchian et al.’s ICCs come from per-occasion ",
  "hierarchical Bayesian fits (partial pooling across subjects), whereas our ",
  "curves are split-half Spearman correlations of independent per-subject ",
  "fits — the comparison is direction-level only, never numeric."
)
writeLines(caption_text, file.path(res_dir, "fig_reconciliation_v2_caption.txt"))

# Machine checks (as in step 08, plus the mix-effect sentence with both crossings)
stopifnot(
  str_detect(caption_text, "suggestive"),
  str_detect(caption_text, "task format"),
  str_detect(caption_text, "problem mix"),
  str_detect(caption_text, "hierarchical Bayesian"),
  str_detect(caption_text, "direction-level"),
  str_detect(caption_text, "necessary but not sufficient"),
  str_detect(caption_text, fixed(paste0(
    "crosses r = .5 at ", cross_clean_total, " total problems on the clean pool"))),
  str_detect(caption_text, fixed(paste0(
    "at ", cross_ext_total, " total problems on the expanded pool")))
)

# ---- Figure data (string-verbatim)
curve_out <- plotted_chr |>
  mutate(
    element        = "curve_point",
    layer          = if_else(construct %in% constructs_headline,
                             "headline", "cautionary_grey"),
    total_problems = as.character(2L * as.integer(k)),
    x_min          = NA_character_,
    x_max          = NA_character_,
    label_value    = NA_character_,
    label_stat     = NA_character_,
    source_file    = if_else(pool == "clean",
                             "results/step2/reliability_by_cell_construct.csv",
                             "results/step4/reliability_extended.csv"),
    source_key     = cell_id
  ) |>
  select(element, layer, pool, construct, k, total_problems,
         r_sb, ci_lo_sb, ci_hi_sb, status,
         x_min, x_max, label_value, label_stat, source_file, source_key)

study_out <- fig_studies |>
  mutate(
    element        = if_else(study_id == "walasek_stewart2021",
                             "study_band", "study_marker"),
    layer          = "prior_study_x_only",
    pool           = NA_character_,
    k              = NA_character_,
    total_problems = NA_character_,
    r_sb           = NA_character_,
    ci_lo_sb       = NA_character_,
    ci_hi_sb       = NA_character_,
    status         = verification,
    x_min          = if_else(study_id == "walasek_stewart2021", "64", "120"),
    x_max          = if_else(study_id == "walasek_stewart2021", "100", "120"),
    label_value    = reliability_value,
    label_stat     = reliability_type,
    source_file    = "results/step3/prior_studies_table.csv",
    source_key     = study_id
  ) |>
  select(element, layer, pool, construct, k, total_problems,
         r_sb, ci_lo_sb, ci_hi_sb, status,
         x_min, x_max, label_value, label_stat, source_file, source_key)

region_out <- tibble(
  element = "beyond_grid_region", layer = "annotation",
  pool = NA_character_, construct = NA_character_, k = NA_character_,
  total_problems = NA_character_, r_sb = NA_character_,
  ci_lo_sb = NA_character_, ci_hi_sb = NA_character_,
  status = paste0("no curve values beyond ", grid_end, " problems (either pool)"),
  x_min = as.character(grid_end), x_max = as.character(x_max),
  label_value = NA_character_, label_stat = NA_character_,
  source_file = NA_character_, source_key = NA_character_
)

figure_data <- bind_rows(curve_out, study_out, region_out)

# Integrity: r_sb / CI strings in figure_data must be IDENTICAL to the
# source-CSV strings (never converted to numeric and back).
check <- plotted_chr |>
  left_join(figure_data |> filter(element == "curve_point"),
            by = c("cell_id" = "source_key", "construct", "pool"),
            suffix = c(".src", ".out"))
stopifnot(
  nrow(check) == 78L,
  identical(check$r_sb.src, check$r_sb.out),
  identical(check$ci_lo_sb.src, check$ci_lo_sb.out),
  identical(check$ci_hi_sb.src, check$ci_hi_sb.out)
)
# Re-read the sources fresh and re-verify against the output file on disk.
write_csv(figure_data, file.path(res_dir, "fig_reconciliation_v2_data.csv"), na = "NA")
reread <- read_csv(file.path(res_dir, "fig_reconciliation_v2_data.csv"),
                   col_types = cols(.default = col_character()))
src_again <- bind_rows(
  read_csv(clean_path, col_types = cols(.default = col_character())) |>
    filter(axis == "i_exp", regime == "exp", status == "VERIFIED") |>
    mutate(pool = "clean"),
  read_csv(ext_path, col_types = cols(.default = col_character())) |>
    filter(axis == "i_exp_ext", regime == "exp", status == "VERIFIED") |>
    mutate(pool = "expanded")
) |> select(cell_id, pool, construct, r_sb, ci_lo_sb, ci_hi_sb)
check2 <- reread |>
  filter(element == "curve_point") |>
  inner_join(src_again, by = c("source_key" = "cell_id", "pool", "construct"),
             suffix = c(".out", ".src"))
stopifnot(
  nrow(check2) == 78L,
  identical(check2$r_sb.out, check2$r_sb.src),
  identical(check2$ci_lo_sb.out, check2$ci_lo_sb.src),
  identical(check2$ci_hi_sb.out, check2$ci_hi_sb.src)
)

# ---- The step 08 figure is NOT deleted (kept as clean-pool-only supplement)
stopifnot(
  file.exists(here("output", "results", "step3", "fig_reconciliation.png")),
  file.exists(here("output", "results", "step3", "fig_reconciliation.pdf")),
  file.exists(here("output", "results", "step3", "figure_data.csv"))
)

cat("step5a OK: figure v2 written;", nrow(plotted_chr),
    "curve points string-verbatim from both pools; grid ends at",
    grid_end, "problems; all assertions passed\n")

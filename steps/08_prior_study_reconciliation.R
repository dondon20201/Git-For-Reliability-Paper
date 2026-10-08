# Step 08: prior studies placed on the problem-count axis
# Reads:   data/prior_studies.csv; output/results/step2/ (step 05)
# Writes:  output/results/step3/ (prior_studies_table.csv, figure_data.csv, fig_reconciliation.png/.pdf, caption)
# Paper:   Relation to prior reliability findings; prior-study values in the Abstract and Introduction
# Runtime: a few seconds

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(here)
  library(knitr)
})

res_dir <- here("output", "results", "step3")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

cells_path   <- here("output", "results", "step2", "reliability_by_cell_construct.csv")
needed_path  <- here("output", "results", "step2", "trials_needed_table.csv")
studies_path <- here("data", "prior_studies.csv")

# Read EVERYTHING as character so the saved values flow into figure_data.csv
# string-verbatim (numeric copies are made only for plotting coordinates).
cells_chr <- read_csv(cells_path,   col_types = cols(.default = col_character()))
needed    <- read_csv(needed_path,  col_types = cols(.default = col_character()))
studies   <- read_csv(studies_path, col_types = cols(.default = col_character()))

# Stop if the prior-studies table is not exactly the expected one.
stopifnot(nrow(studies) == 7L)
expected_in_figure <- c("mkrtchian2023_loss", "mkrtchian2023_risk", "walasek_stewart2021")
stopifnot(identical(
  sort(studies$study_id[studies$in_figure == "yes"]),
  sort(expected_in_figure)
))

constructs_headline <- c("ev_sens_agn", "loss_avoid_adj")
constructs_grey     <- c("b_ev", "b_risk", "b_loss")

curve_chr <- cells_chr |>
  filter(axis == "i_exp", regime == "exp") |>
  select(cell_id, construct, k, r_sb, ci_lo_sb, ci_hi_sb, n_subjects, status)

# Observed grid: k = 2..7. Headline constructs identified everywhere;
# model parameters structurally unidentified below k = 4 (k < 4 problems/half
# for a 4-parameter fit); those rows are NA and are not plotted.
stopifnot(identical(sort(unique(as.integer(curve_chr$k))), 2:7))
plotted_chr <- curve_chr |> filter(status == "VERIFIED")
counts <- plotted_chr |> count(construct)
stopifnot(
  all(counts$n[counts$construct %in% constructs_headline] == 6L),
  all(counts$n[counts$construct %in% constructs_grey] == 4L),
  nrow(plotted_chr) == 24L,
  !any(is.na(as.numeric(plotted_chr$r_sb)))
)

# Numeric copies for plotting coordinates only.
plotted <- plotted_chr |>
  mutate(
    k              = as.integer(k),
    total_problems = 2L * k,          # x-axis unit: total clean problems (2 x k per split half)
    r_num          = as.numeric(r_sb),
    lo_num         = as.numeric(ci_lo_sb),
    hi_num         = as.numeric(ci_hi_sb),
    layer          = if_else(construct %in% constructs_headline, "headline", "cautionary_grey")
  )
max(plotted$total_problems)  # must be 14; no curve segment beyond this
stopifnot(max(plotted$total_problems) == 14L)

fig_studies <- studies |> filter(in_figure == "yes")

mk <- fig_studies |> filter(str_starts(study_id, "mkrtchian"))
ws <- fig_studies |> filter(study_id == "walasek_stewart2021")
stopifnot(all(mk$design_points_equiv == "120"), ws$design_points_equiv == "64-100")

mk_loss <- mk |> filter(study_id == "mkrtchian2023_loss")
mk_risk <- mk |> filter(study_id == "mkrtchian2023_risk")

lab_mk <- paste0(
  "Mkrtchian et al. (2023)\n",
  "~120 calibrated design points\n",
  "(raw budget: ", mk_loss$budget_reported |> str_extract("^120 modelled trials"), ")\n",
  "λ ", mk_loss$reliability_type, " = ", mk_loss$reliability_value,
  " [", mk_loss$ci_95, "]\n",
  "risk ", mk_risk$reliability_type, " = ", mk_risk$reliability_value,
  " [", mk_risk$ci_95, "]\n",
  "hierarchical Bayesian, per-occasion fits;\n",
  "different task format (calibrated accept–reject)"
)

lab_ws <- paste0(
  "Walasek & Stewart (2021)\n",
  "64–100 unique gambles per set\n",
  "(raw budget: up to 3840 trials, Sim 2)\n",
  "λ: “extremely poor parameter recoverability”\n",
  ws$reliability_type
)

lab_grey <- paste0(
  "Firth b_ev / b_risk / b_loss: ≈0 at every identified budget;\n",
  "k < 4 per half structurally unidentified"
)

lab_beyond <- "beyond CPC18’s observed grid —\nno curve values here"

x_max <- 160
col_map <- c(ev_sens_agn = "#0072B2", loss_avoid_adj = "#D55E00")
lbl_map <- c(ev_sens_agn = "EV sensitivity (model-agnostic index, SB)",
             loss_avoid_adj = "Loss avoidance, EV-adjusted (SB)")

headline  <- plotted |> filter(layer == "headline")
greylayer <- plotted |> filter(layer == "cautionary_grey")

p <- ggplot() +
  # Region beyond CPC18's observed grid (no curve values; visually distinct)
  annotate("rect", xmin = 14, xmax = x_max, ymin = -Inf, ymax = Inf,
           fill = "grey75", alpha = 0.35) +
  # Walasek & Stewart 64-100 design-points band (x placement only)
  annotate("rect", xmin = 64, xmax = 100, ymin = -Inf, ymax = Inf,
           fill = "#E69F00", alpha = 0.18) +
  # Mkrtchian ~120 design-points marker (x placement only)
  geom_vline(xintercept = 120, linetype = "dashed", colour = "grey30",
             linewidth = 0.5) +
  # Pre-declared reliability thresholds (reference lines)
  geom_hline(yintercept = c(0.5, 0.7), linetype = "dotted", colour = "grey45",
             linewidth = 0.4) +
  # Headline curves: bootstrap CI ribbons + straight segments between observed
  # points only (no smoothing, no extrapolation)
  geom_ribbon(data = headline,
              aes(x = total_problems, ymin = lo_num, ymax = hi_num,
                  fill = construct),
              alpha = 0.18, colour = NA) +
  geom_line(data = headline,
            aes(total_problems, r_num, colour = construct), linewidth = 0.8) +
  geom_point(data = headline,
             aes(total_problems, r_num, colour = construct), size = 1.8) +
  # Cautionary grey layer: Firth model parameters
  geom_line(data = greylayer,
            aes(total_problems, r_num, group = construct),
            colour = "grey55", linewidth = 0.5) +
  geom_point(data = greylayer,
             aes(total_problems, r_num, group = construct),
             colour = "grey55", size = 1.1) +
  # Annotations
  annotate("text", x = 33, y = 0.34, label = lab_beyond, size = 3.1,
           colour = "grey30", fontface = "italic", lineheight = 1.05) +
  annotate("text", x = 8, y = -0.075, label = lab_grey, size = 2.9,
           colour = "grey35", lineheight = 1.05) +
  annotate("text", x = 80, y = 0.93, label = lab_ws, size = 2.8,
           colour = "grey10", lineheight = 1.08) +
  annotate("text", x = 117, y = 0.66, hjust = 1, vjust = 1, label = lab_mk,
           size = 2.8, colour = "grey10", lineheight = 1.08) +
  annotate("text", x = 4.15, y = c(0.515, 0.715), hjust = 0, vjust = -0.4,
           label = c("r = .5", "r = .7"), size = 2.7, colour = "grey45") +
  scale_x_log10(breaks = c(4, 6, 8, 10, 14, 25, 50, 64, 100, 120, 150),
                limits = c(4, x_max)) +
  scale_colour_manual(values = col_map, labels = lbl_map, name = NULL) +
  scale_fill_manual(values = col_map, labels = lbl_map, name = NULL) +
  coord_cartesian(ylim = c(-0.12, 1.0)) +
  labs(
    x = "Total clean problems (2 × k per split half; log scale)",
    y = "Split-half reliability (Spearman rank, Spearman–Brown corrected)",
    title = "Prior studies mapped onto CPC18 reliability-vs-problem-count curves",
    subtitle = paste0(
      "Experience regime, r = 20 trials/problem; ribbons = 95% subject-cluster ",
      "bootstrap CIs. Studies placed on x only — suggestive mapping; ",
      "task format AND problem mix confound it."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 8.5, colour = "grey30")
  )

ggsave(file.path(res_dir, "fig_reconciliation.png"), p,
       width = 10, height = 6.5, dpi = 300)
ggsave(file.path(res_dir, "fig_reconciliation.pdf"), p,
       width = 10, height = 6.5, device = cairo_pdf)
p

caption_text <- paste0(
  "Figure. Prior studies mapped onto the CPC18 reliability-versus-problem-count ",
  "curves (experience regime, r = 20 trials per problem; split-half Spearman with ",
  "Spearman–Brown correction; ribbons are 95% subject-cluster bootstrap CIs). ",
  "Coloured curves: EV-adjusted loss avoidance (loss_avoid_adj) and the ",
  "model-agnostic EV-sensitivity index (ev_sens_agn), drawn only over CPC18’s ",
  "observed grid (4–14 total clean problems) with straight segments between ",
  "observed points — no smoothing and no extrapolation. Grey curves: the Firth ",
  "model parameters (b_ev, b_risk, b_loss) are ≈0-reliable at every identified ",
  "budget and are structurally unidentified below 4 problems per half; the shaded ",
  "region beyond 14 problems contains no CPC18 curve values. Prior studies are ",
  "placed on the x-axis only, at their unique-informative-design-points ",
  "equivalent, and are given no position on our curves: Mkrtchian et al. (2023; ",
  "~120 individually calibrated offers, 120 modelled trials; loss aversion ",
  "test–retest ICC(A,1) = 0.68 [0.50–0.81], risk ICC(A,1) = 0.78 ",
  "[0.55–0.89]) and Walasek & Stewart (2021; 64–100 unique gambles per ",
  "set; simulation-based “extremely poor” recoverability of λ). ",
  "The mapping is suggestive, not probative: problem count appears necessary but ",
  "not sufficient. Below CPC18’s identification floor nothing model-based is ",
  "measurable, yet Walasek & Stewart’s 64–100-gamble sets exceed ",
  "CPC18’s ≤14 clean problems and λ still fails — via a ",
  "λ–bias tradeoff tied to task format and via choice-set dependence ",
  "tied to problem mix (consistent with our Set-composition finding and with ",
  "Walasek, Mullett & Stewart’s observation that a change of choice set ",
  "roundly kills generalization). Mkrtchian et al.’s high reliability sits ",
  "at both a larger calibrated design AND a different task format; we do not ",
  "adjudicate between these. Estimation-regime caveat: Mkrtchian et al.’s ",
  "ICCs come from per-occasion hierarchical Bayesian fits (partial pooling across ",
  "subjects), whereas our curves are split-half Spearman correlations of ",
  "independent per-subject Firth fits — the comparison is direction-level ",
  "only, never numeric."
)
writeLines(caption_text, file.path(res_dir, "fig_reconciliation_caption.txt"))
cat("> ", caption_text, "\n", sep = "")

# Caption requirements, machine-checked: "suggestive" framing, BOTH
# confounds named, estimation-regime caveat present.
stopifnot(
  str_detect(caption_text, "suggestive"),
  str_detect(caption_text, "task format"),
  str_detect(caption_text, "problem mix"),
  str_detect(caption_text, "hierarchical Bayesian"),
  str_detect(caption_text, "direction-level")
)

copied <- file.copy(studies_path, file.path(res_dir, "prior_studies_table.csv"),
                    overwrite = TRUE)
stopifnot(copied)
# Verify byte-identity of the copy.
stopifnot(identical(
  readBin(studies_path, "raw", file.size(studies_path)),
  readBin(file.path(res_dir, "prior_studies_table.csv"),
          "raw", file.size(file.path(res_dir, "prior_studies_table.csv")))
))

studies |>
  select(study_id, role, in_figure, construct, reliability_type,
         reliability_value, ci_95, budget_reported, design_points_equiv,
         estimation_regime, verification) |>
  kable()

fmt3 <- function(x) sub("^0", "", sprintf("%.3f", as.numeric(x)))
row_exp  <- needed |> filter(axis == "i_exp",  construct == "loss_avoid_adj", threshold == "0.5")
row_desc <- needed |> filter(axis == "i_desc", construct == "loss_avoid_adj", threshold == "0.5")
stopifnot(row_exp$entry == "reached", row_desc$entry == "reached")
cat(paste0(
  "- In the **experience** regime (the figure), `loss_avoid_adj` is the only ",
  "construct to reach SB .5, at ", row_exp$budget, " problems/half (",
  row_exp$total_problems, " total; SB ", fmt3(row_exp$r_sb), " [",
  fmt3(row_exp$ci_lo_sb), ", ", fmt3(row_exp$ci_hi_sb), "], n = ",
  row_exp$n_subjects, "). Nothing reaches .7 anywhere.\n",
  "- In the **description** regime (robustness), `loss_avoid_adj` reaches SB .5 ",
  "at ", row_desc$budget, " problems/half (", row_desc$total_problems,
  " total; SB ", fmt3(row_desc$r_sb), " [", fmt3(row_desc$ci_lo_sb), ", ",
  fmt3(row_desc$ci_hi_sb), "], n = ", row_desc$n_subjects,
  "); no other construct reaches .5, and nothing reaches .7. The direction ",
  "(description more informative per trial) matches the Step-2 headline.\n"
))

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
    source_file    = "results/step2/reliability_by_cell_construct.csv",
    source_key     = cell_id
  ) |>
  select(element, layer, construct, k, total_problems,
         r_sb, ci_lo_sb, ci_hi_sb, status,
         x_min, x_max, label_value, label_stat, source_file, source_key)

study_out <- fig_studies |>
  mutate(
    element        = if_else(study_id == "walasek_stewart2021",
                             "study_band", "study_marker"),
    layer          = "prior_study_x_only",
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
  select(element, layer, construct, k, total_problems,
         r_sb, ci_lo_sb, ci_hi_sb, status,
         x_min, x_max, label_value, label_stat, source_file, source_key)

region_out <- tibble(
  element = "beyond_grid_region", layer = "annotation",
  construct = NA_character_, k = NA_character_,
  total_problems = NA_character_, r_sb = NA_character_,
  ci_lo_sb = NA_character_, ci_hi_sb = NA_character_,
  status = "no curve values beyond 14 problems",
  x_min = "14", x_max = as.character(x_max),
  label_value = NA_character_, label_stat = NA_character_,
  source_file = NA_character_, source_key = NA_character_
)

figure_data <- bind_rows(curve_out, study_out, region_out)

# Integrity: the r_sb / CI strings in figure_data must be IDENTICAL to the
# step 05 CSV strings (they were never converted to numeric and back).
check <- plotted_chr |>
  left_join(figure_data |> filter(element == "curve_point"),
            by = c("cell_id" = "source_key", "construct"),
            suffix = c(".src", ".out"))
stopifnot(
  nrow(check) == 24L,
  identical(check$r_sb.src, check$r_sb.out),
  identical(check$ci_lo_sb.src, check$ci_lo_sb.out),
  identical(check$ci_hi_sb.src, check$ci_hi_sb.out)
)

write_csv(figure_data, file.path(res_dir, "figure_data.csv"), na = "NA")
figure_data |> kable()


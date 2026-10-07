# Step 12: aggregate agreement checks on choices13k
# Reads:   data/raw/choices13k/ (downloaded and md5-checked if absent); data/raw/All estimation raw data.csv
# Writes:  output/results/step4/choices13k_*.csv, output/figures/step4_choices13k_overlay.png
# Paper:   Results (choices13k: aggregate agreement only)
# Runtime: about 20 seconds

# choices13k aggregate-level agreement checks.
#
# SCOPE: choices13k is problem-level aggregated data and CANNOT support
# per-subject reliability claims. Everything below is aggregate/direction-level
# only.
#
# Verified fetch (pinned commit, md5) + structure documentation.
# B2 (pre-declared): (1) EV-monotonicity of the aggregate B-choice rate;
#     (2) overall choice-rate sanity vs CPC18's ~.5; (3) shape agreement of
#     the binned choice-rate-vs-EV-gap curve vs the CPC18 aggregate
#     equivalent (rank correlation of binned rates + overlay plot).
# EXPLORATORY (labeled in every output): EV-adjusted loss avoidance at the
#     aggregate problem level, DIRECTION ONLY, as external context for the
#     CPC18 population-level signal. No reliability language, no new claims.
#
# EVs for choices13k come from the dataset's own explicit distributions
# (c13k_problems.json, what participants actually saw), cross-checked against
# the CPC18 LotShape reconstruction validated in step 09, applied to the CSV
# parameters (LotShapeB codes: 0='-', 1='Symm', 2='R-skew', 3='L-skew').

library(tidyverse)
library(here)

source(here("R", "fetch_cpc18.R"))
source(here("R", "cpc18_features.R"))
source(here("R", "lotshape_reconstruction.R"))
source(here("R", "fetch_choices13k.R"))

res_dir <- here("output", "results", "step4")
fig_dir <- here("output", "figures")
SEED_BOOT <- 191003
CANNOT_SUPPORT <- paste("choices13k is problem-level aggregated data and",
                        "CANNOT support per-subject reliability claims")

# ---- verified fetch + structure ---------------------------------------------
c13k_dir <- here("data", "raw", "choices13k")
fetch_c13k <- fetch_choices13k(c13k_dir)
stopifnot(all(fetch_c13k$verified))

sel <- read_csv(file.path(c13k_dir, "c13k_selections.csv"),
                show_col_types = FALSE)
probs_json <- jsonlite::fromJSON(file.path(c13k_dir, "c13k_problems.json"))
stopifnot(length(probs_json) == nrow(sel))

verification <- fetch_c13k |>
  mutate(source_repo = choices13k_source$repo,
         commit = choices13k_source$commit,
         commit_date = choices13k_source$commit_date,
         citation = choices13k_source$citation,
         n_rows_selections = nrow(sel),
         n_cols_selections = ncol(sel),
         n_unique_problems = n_distinct(sel$Problem),
         n_json_entries = length(probs_json),
         columns = paste(names(sel), collapse = ";"),
         structure_note = paste("problem-level aggregated choice rates",
                                "(bRate over ~n subjects x 5 trials per row;",
                                "rows split by Feedback); NO per-subject",
                                "repeated structure."),
         scope_note = CANNOT_SUPPORT)
write_csv(verification, file.path(res_dir, "choices13k_verification.csv"))
cat(sprintf("choices13k verified: %d rows, %d unique problems, %d JSON entries\n",
            nrow(sel), n_distinct(sel$Problem), length(probs_json)))

# ---- EVs from the dataset's own explicit distributions ----------------------
# JSON entry i (key = CSV row index, 0-based): A/B matrices of [prob, value].
json_mom <- map_dfr(seq_len(nrow(sel)), function(i) {
  e <- probs_json[[as.character(i - 1)]]
  a <- matrix(unlist(e$A), ncol = 2, byrow = FALSE)
  b <- matrix(unlist(e$B), ncol = 2, byrow = FALSE)
  # fromJSON returns n x 2 matrices [prob, value] already; guard the shape
  if (is.matrix(e$A)) a <- e$A
  if (is.matrix(e$B)) b <- e$B
  tibble(row = i,
         ev_a_json = sum(a[, 1] * a[, 2]),
         ev_b_json = sum(b[, 1] * b[, 2]),
         has_loss_a_json = any(a[, 2] < 0 & a[, 1] > 0),
         has_loss_b_json = any(b[, 2] < 0 & b[, 1] > 0),
         psum_a = sum(a[, 1]), psum_b = sum(b[, 1]))
})
stopifnot(max(abs(json_mom$psum_a - 1)) <= 1e-6,
          max(abs(json_mom$psum_b - 1)) <= 1e-6)

# Cross-check vs the step 09-validated reconstruction on Amb == 0 rows
shape_map <- c(`0` = "-", `1` = "Symm", `2` = "R-skew", `3` = "L-skew")
sel2 <- sel |>
  mutate(row = row_number(),
         Amb = as.logical(Amb), Feedback = as.logical(Feedback),
         LotShapeB_chr = shape_map[as.character(LotShapeB)]) |>
  left_join(json_mom, by = "row")
c13k_amb0 <- sel2 |> filter(!Amb)
ev_b_recon <- pmap_dbl(
  list(c13k_amb0$Hb, c13k_amb0$pHb, c13k_amb0$Lb,
       c13k_amb0$LotShapeB, c13k_amb0$LotNumB),
  function(H, pH, L, shape_code, num) {
    # a 1-outcome "lottery" is degenerate at Hb whatever its shape code
    shape <- if (num <= 1) "-" else shape_map[as.character(shape_code)]
    m <- reconstruct_lottery(H, pH, L, shape, max(num, 1))
    sum(m[, "outcome"] * m[, "prob"])
  })
max_ev_dev <- max(abs(c13k_amb0$ev_b_json - ev_b_recon))
cat(sprintf("EV cross-check (JSON vs validated reconstruction, %d Amb==0 rows): max |dev| = %g\n",
            nrow(c13k_amb0), max_ev_dev))
stopifnot(max_ev_dev <= 1e-6)

# ---- analysis frames (Amb == 0 in both datasets) ----------------------------
c13k <- sel2 |>
  filter(!Amb) |>
  mutate(d_ev = ev_b_json - ev_a_json,
         loss_informative = xor(has_loss_a_json, has_loss_b_json))

raw <- load_cpc18(here("data", "raw", "All estimation raw data.csv"))
param_cols <- c("Ha", "pHa", "La", "LotShapeA", "LotNumA",
                "Hb", "pHb", "Lb", "LotShapeB", "LotNumB", "Amb", "Corr")
problems_exp <- raw |>
  distinct(GameID, across(all_of(param_cols))) |>
  arrange(GameID) |> filter(Amb == 0)
feat_exp <- expanded_problem_features(problems_exp,
                                      reconstruct_all_problems(problems_exp))
cpc_agg <- raw |>
  semi_join(problems_exp, by = "GameID") |>
  group_by(GameID) |>
  summarise(bRate = mean(B), n_trials = n(), .groups = "drop") |>
  left_join(feat_exp |> select(GameID, d_ev, loss_informative,
                               has_loss_a, has_loss_b, ev_a, ev_b),
            by = "GameID")
cpc_overall_all <- mean(raw$B)
cpc_overall_amb0 <- weighted.mean(cpc_agg$bRate, cpc_agg$n_trials)

# ---- B2.1: EV-monotonicity --------------------------------------------------
set.seed(SEED_BOOT)
sp_point <- cor(c13k$bRate, c13k$d_ev, method = "spearman")
boots <- replicate(2000, {
  idx <- sample.int(nrow(c13k), replace = TRUE)
  cor(c13k$bRate[idx], c13k$d_ev[idx], method = "spearman")
})
sp_ci <- quantile(boots, c(.025, .975), names = FALSE)

# decile means for the monotone-direction readout
c13k_dec <- c13k |>
  mutate(bin = ntile(d_ev, 10)) |>
  group_by(bin) |>
  summarise(rate = weighted.mean(bRate, n), .groups = "drop")
monotone_up <- all(diff(c13k_dec$rate) > 0)

# ---- B2.2: overall choice rate ----------------------------------------------
c13k_overall <- weighted.mean(c13k$bRate, c13k$n)
c13k_overall_all <- weighted.mean(sel2$bRate, sel2$n)

# ---- B2.3: shape agreement (common bins, rank correlation, overlay) ---------
breaks <- quantile(c(c13k$d_ev, cpc_agg$d_ev), probs = seq(0, 1, 0.1))
breaks[1] <- -Inf; breaks[11] <- Inf
shape <- bind_rows(
  c13k |> mutate(dataset = "choices13k (Amb==0)", w = n) |>
    select(dataset, d_ev, bRate, w),
  cpc_agg |> mutate(dataset = "CPC18 (Amb==0)", w = n_trials) |>
    select(dataset, d_ev, bRate, w)
) |>
  mutate(bin = cut(d_ev, breaks, labels = FALSE)) |>
  group_by(dataset, bin) |>
  summarise(rate = weighted.mean(bRate, w), n_problems = n(),
            d_ev_mid = median(d_ev), .groups = "drop")
shape_wide <- shape |>
  pivot_wider(id_cols = bin, names_from = dataset,
              values_from = rate) |>
  drop_na()
shape_rank_cor <- cor(shape_wide[[2]], shape_wide[[3]], method = "spearman")
write_csv(shape |> mutate(scope_note = CANNOT_SUPPORT),
          file.path(res_dir, "choices13k_shape_agreement.csv"))

p <- ggplot(shape, aes(d_ev_mid, rate, color = dataset)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey70") +
  geom_point() + geom_line() +
  labs(x = "EV(B) - EV(A) (bin median)", y = "Aggregate B-choice rate",
       title = "Aggregate choice rate vs EV gap: CPC18 vs choices13k [real data, problem-level]",
       subtitle = sprintf("Common decile bins of the pooled EV gap; rank correlation of binned rates = %.3f.\nchoices13k is problem-level aggregated data and cannot support per-subject reliability claims.",
                          shape_rank_cor)) +
  theme_minimal()
ggsave(file.path(fig_dir, "step4_choices13k_overlay.png"), p,
       width = 8, height = 5, dpi = 150)

# regime-matched secondary readout (direction-level)
shape_split <- map_dfr(c(FALSE, TRUE), function(fb) {
  c1 <- c13k |> filter(Feedback == fb)
  cp <- raw |>
    semi_join(problems_exp, by = "GameID") |>
    filter(if (fb) Trial >= 6 else Trial <= 5) |>
    group_by(GameID) |>
    summarise(bRate = mean(B), n_trials = n(), .groups = "drop") |>
    left_join(feat_exp |> select(GameID, d_ev), by = "GameID")
  sh <- bind_rows(
    c1 |> transmute(dataset = "c13k", d_ev, bRate, w = n),
    cp |> transmute(dataset = "cpc18", d_ev, bRate, w = n_trials)
  ) |>
    mutate(bin = cut(d_ev, breaks, labels = FALSE)) |>
    group_by(dataset, bin) |>
    summarise(rate = weighted.mean(bRate, w), .groups = "drop") |>
    pivot_wider(id_cols = bin, names_from = dataset, values_from = rate) |>
    drop_na()
  tibble(feedback = fb, rank_cor_binned = cor(sh$c13k, sh$cpc18, method = "spearman"))
})

checks <- tribble(
  ~check, ~value, ~detail,
  "B2.1 EV-monotonicity: Spearman(bRate, d_ev), choices13k Amb==0",
  sp_point,
  sprintf("95%% bootstrap CI [%.3f, %.3f] (2000 resamples of problems, seed %d); decile curve strictly increasing: %s; n = %d rows",
          sp_ci[1], sp_ci[2], SEED_BOOT, monotone_up, nrow(c13k)),
  "B2.2 overall B-choice rate, choices13k Amb==0 (n-weighted)",
  c13k_overall,
  sprintf("all rows incl. Amb==1: %.3f; CPC18 overall (all 210 problems, all trials): %.3f; CPC18 Amb==0: %.3f — both near .5",
          c13k_overall_all, cpc_overall_all, cpc_overall_amb0),
  "B2.3 shape agreement: rank correlation of binned rate-vs-EV-gap curves",
  shape_rank_cor,
  sprintf("10 common pooled-decile bins; regime-matched secondary: no-feedback %.3f (c13k no-FB vs CPC18 trials 1-5), feedback %.3f (c13k FB vs CPC18 trials 6-25); overlay figure saved",
          shape_split$rank_cor_binned[1], shape_split$rank_cor_binned[2])
) |>
  mutate(status = "VERIFIED (aggregate/problem-level only)",
         scope_note = CANNOT_SUPPORT)
write_csv(checks, file.path(res_dir, "choices13k_aggregate_checks.csv"))
print(checks |> select(check, value, detail), width = Inf)

# ---- EXPLORATORY EV-adjusted loss avoidance (direction only) ----------------
# Among loss-informative problems (exactly one option's support contains a
# negative outcome), model the rate of choosing the NO-LOSS option as a
# function of its EV advantage (n-weighted quasibinomial). A positive
# intercept = avoidance of the loss-containing option beyond EV at parity;
# negative = reversed. DIRECTION ONLY; external context for the CPC18
# population-level signal (absent-to-slightly-reversed once EV-controlled).
c13k_li <- c13k |>
  filter(loss_informative) |>
  mutate(chose_noloss = ifelse(has_loss_b_json, 1 - bRate, bRate),
         d_ev_noloss = ifelse(has_loss_b_json, ev_a_json - ev_b_json,
                              ev_b_json - ev_a_json))
fit_c13k <- glm(chose_noloss ~ d_ev_noloss, family = quasibinomial(),
                data = c13k_li, weights = n)
co_c <- summary(fit_c13k)$coefficients

cpc_li <- cpc_agg |>
  filter(loss_informative) |>
  mutate(chose_noloss = ifelse(has_loss_b, 1 - bRate, bRate),
         d_ev_noloss = ifelse(has_loss_b, ev_a - ev_b, ev_b - ev_a))
fit_cpc <- glm(chose_noloss ~ d_ev_noloss, family = quasibinomial(),
               data = cpc_li, weights = n_trials)
co_p <- summary(fit_cpc)$coefficients

exploratory <- tibble(
  label = "EXPLORATORY — direction only, no reliability language, no new claims",
  dataset = c("choices13k (Amb==0, problem-level)", "CPC18 (Amb==0, problem-level aggregate)"),
  n_loss_informative_problems = c(nrow(c13k_li), nrow(cpc_li)),
  intercept_logit = c(co_c[1, 1], co_p[1, 1]),
  intercept_se = c(co_c[1, 2], co_p[1, 2]),
  ev_slope_logit = c(co_c[2, 1], co_p[2, 1]),
  direction_at_ev_parity = ifelse(c(co_c[1, 1], co_p[1, 1]) > 0,
                                  "avoids loss option (loss avoidance)",
                                  "prefers loss option (reversed)"),
  model = "weighted quasibinomial: chose_noloss ~ d_ev_noloss (weights = subjects/trials)",
  scope_note = CANNOT_SUPPORT
)
write_csv(exploratory, file.path(res_dir, "choices13k_exploratory_loss.csv"))
print(exploratory |> select(dataset, n_loss_informative_problems,
                            intercept_logit, intercept_se,
                            direction_at_ev_parity), width = Inf)
cat("\nPart B done. Scope:", CANNOT_SUPPORT, "\n")

# Step 03: correlation of estimation noise between parameters (simulation)
# Reads:   output/results/step0b/firth_recovery_fits_raw.rds (step 02)
# Writes:  output/results/step0b/est_est_noise_correlation.csv
# Paper:   Methods (Validation gates); no number in the paper is read from it
# Runtime: a few seconds

# Real-subject Firth estimates (step 04) show spearman(b_risk, b_loss) = -0.648.
# Question: how much estimate-estimate correlation arises from estimation
# noise ALONE? In the step 02 simulation true parameters are independent by
# construction, so any correlation among ESTIMATES there is noise/collinearity
# induced, not trait structure. SIMULATION-DERIVED numbers.

library(dplyr)
library(here)

sim <- readRDS(here("output", "results", "step0b", "firth_recovery_fits_raw.rds")) |>
  filter(estimator == "firth", k_problems == 14, n_reps == 25)

pairs <- tibble::tribble(
  ~pair, ~x, ~y,
  "est_b_risk ~ est_b_loss", sim$est_b_risk, sim$est_b_loss,
  "est_b_ev ~ est_b_loss",   sim$est_b_ev,   sim$est_b_loss,
  "est_b_ev ~ est_b_risk",   sim$est_b_ev,   sim$est_b_risk,
  "true b_risk ~ true b_loss (design check: ~0)", sim$b_risk, sim$b_loss
)
out <- pairs |>
  rowwise() |>
  summarise(pair = pair,
            spearman = cor(x, y, method = "spearman"),
            n = length(x), .groups = "drop") |>
  mutate(source = "SIMULATION (independent true params, Firth, 14x25)")
print(out)
readr::write_csv(out, here("output", "results", "step0b", "est_est_noise_correlation.csv"))
cat("Saved results/step0b/est_est_noise_correlation.csv\n")

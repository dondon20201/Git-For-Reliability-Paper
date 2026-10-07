# Load and featurize the CPC18 estimation raw data (fetched by R/fetch_cpc18.R).
#
# Encodings follow the official data dictionary PDF shipped with Zenodo record
# 845873: an option pays a draw from its "high lottery" with probability pH and
# the low payoff L otherwise; with LotNum == 1 the high lottery is degenerate
# at H, so the option is a simple two-outcome gamble (H with pH, else L).
# The clean subset for the main analysis is LotNumA==1 & LotNumB==1 & Amb==0.

cpc18_expected_columns <- c(
  "SubjID", "Location", "Gender", "Age", "Set", "Condition", "GameID",
  "Ha", "pHa", "La", "LotShapeA", "LotNumA",
  "Hb", "pHb", "Lb", "LotShapeB", "LotNumB",
  "Amb", "Corr", "Order", "Trial", "Button", "B",
  "Payoff", "Forgone", "RT", "Apay", "Bpay", "Feedback", "block"
)

load_cpc18 <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

# TRUE for rows in the clean two-outcome, no-ambiguity subset.
is_clean_row <- function(df) {
  df$LotNumA == 1 & df$LotNumB == 1 & df$Amb == 0
}

# One row per clean problem, with description-based moments of each option and
# the difference features used by the choice model:
#   d_ev    = EV(B) - EV(A)
#   d_sd    = SD(B) - SD(A)      (variance / risk difference)
#   d_eloss = E[loss magnitude](B) - E[loss magnitude](A),
#             where E[loss magnitude] = sum over outcomes of p * max(-x, 0)
# loss_informative = exactly one option contains a negative outcome.
clean_problem_features <- function(df) {
  probs <- df |>
    dplyr::filter(is_clean_row(df)) |>
    dplyr::distinct(GameID, Ha, pHa, La, Hb, pHb, Lb)
  stopifnot(anyDuplicated(probs$GameID) == 0)

  two_outcome_moments <- function(h, p, l) {
    ev <- h * p + l * (1 - p)
    v <- p * h^2 + (1 - p) * l^2 - ev^2
    tibble::tibble(
      ev = ev,
      sd = sqrt(pmax(v, 0)),
      eloss = pmax(-h, 0) * p + pmax(-l, 0) * (1 - p),
      has_loss = h < 0 | l < 0
    )
  }
  a <- two_outcome_moments(probs$Ha, probs$pHa, probs$La)
  b <- two_outcome_moments(probs$Hb, probs$pHb, probs$Lb)

  probs |>
    dplyr::mutate(
      ev_a = a$ev, sd_a = a$sd, eloss_a = a$eloss, has_loss_a = a$has_loss,
      ev_b = b$ev, sd_b = b$sd, eloss_b = b$eloss, has_loss_b = b$has_loss,
      d_ev = ev_b - ev_a,
      d_sd = sd_b - sd_a,
      d_eloss = eloss_b - eloss_a,
      loss_informative = xor(has_loss_a, has_loss_b)
    )
}

# Scale the difference features by their SDs across the clean problem set, so
# model coefficients read as logits per problem-set SD. The scaling constants
# are fixed dataset-level constants shared by all (real or simulated) subjects;
# they rescale coefficients uniformly and cannot affect correlations.
add_scaled_features <- function(features) {
  scales <- c(
    d_ev = stats::sd(features$d_ev),
    d_sd = stats::sd(features$d_sd),
    d_eloss = stats::sd(features$d_eloss)
  )
  out <- features |>
    dplyr::mutate(
      z_ev = d_ev / scales[["d_ev"]],
      z_sd = d_sd / scales[["d_sd"]],
      z_eloss = d_eloss / scales[["d_eloss"]]
    )
  attr(out, "scales") <- scales
  out
}

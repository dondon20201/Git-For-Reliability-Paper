# Multi-outcome lottery reconstruction for CPC18.
#
# An option in CPC18 is (H, pH, L, LotShape, LotNum): with probability pH the
# payoff is a draw from the "high lottery" (LotNum outcomes shaped around H),
# otherwise the payoff is L. Encoding per the official CPC18 baseline
# implementation (CPC18_getDist; fetched and md5-verified by
# R/fetch_official_cpc18.R) and the Zenodo data-dictionary PDF:
#   LotShape "-"      : degenerate high lottery at H (LotNum == 1); the option
#                       is {(L, 1-pH), (H, pH)}, or {(H, 1)} when pH == 1
#                       (L is never paid and is dropped by the official code).
#   LotShape "Symm"   : LotNum outcomes H - k/2 + i, i = 0..k with k = LotNum-1,
#                       probabilities pH * Binom(i; k, 1/2).
#   LotShape "R-skew" : outcomes H - (LotNum+1) + 2^i, i = 1..LotNum,
#                       probabilities pH / 2^i, the last one doubled.
#   LotShape "L-skew" : outcomes H + (LotNum+1) - 2^i, i = 1..LotNum,
#                       probabilities pH / 2^i, the last one doubled.
# In the multi-outcome case L is then merged in: if L equals an existing
# outcome its probability gains (1-pH), otherwise the outcome (L, 1-pH) is
# appended (only when pH < 1), and rows are sorted ascending by outcome.
# The validation gate (step 09) asserts this
# function matches the official implementation exactly on all 210 problems.

# One option -> matrix with columns (outcome, prob), in the official row order
# ("-" branch: L then H; multi-outcome branches: sorted ascending by outcome).
reconstruct_lottery <- function(H, pH, L, LotShape, LotNum) {
  stopifnot(LotShape %in% c("-", "Symm", "L-skew", "R-skew"))
  if (LotShape == "-") {
    stopifnot(LotNum == 1)
    if (pH == 1) {
      m <- cbind(outcome = H, prob = pH)
    } else {
      m <- cbind(outcome = c(L, H), prob = c(1 - pH, pH))
    }
    return(m)
  }
  stopifnot(LotNum >= 2)
  if (LotShape == "Symm") {
    k <- LotNum - 1
    i <- 0:k
    outcome <- H - k / 2 + i
    prob <- pH * stats::dbinom(i, k, 0.5)
  } else {
    C <- if (LotShape == "R-skew") -1 - LotNum else 1 + LotNum
    distsign <- if (LotShape == "R-skew") 1 else -1
    i <- seq_len(LotNum)
    outcome <- H + C + distsign * 2^i
    prob <- pH / 2^i
    prob[LotNum] <- prob[LotNum] * 2
  }
  # merge the low outcome L into the high-lottery distribution
  j <- match(L, outcome)
  if (!is.na(j)) {
    prob[j] <- prob[j] + (1 - pH)
  } else if (pH < 1) {
    outcome <- c(outcome, L)
    prob <- c(prob, 1 - pH)
  }
  ord <- order(outcome)
  cbind(outcome = outcome[ord], prob = prob[ord])
}

# All problems -> long tibble (GameID, option, outcome_index, outcome, prob).
# `problems` needs GameID, Ha, pHa, La, LotShapeA, LotNumA, Hb, pHb, Lb,
# LotShapeB, LotNumB (one row per problem).
reconstruct_all_problems <- function(problems) {
  purrr::map_dfr(seq_len(nrow(problems)), function(r) {
    p <- problems[r, ]
    purrr::map_dfr(c("A", "B"), function(opt) {
      m <- if (opt == "A") {
        reconstruct_lottery(p$Ha, p$pHa, p$La, p$LotShapeA, p$LotNumA)
      } else {
        reconstruct_lottery(p$Hb, p$pHb, p$Lb, p$LotShapeB, p$LotNumB)
      }
      tibble::tibble(
        GameID = p$GameID, option = opt,
        outcome_index = seq_len(nrow(m)),
        outcome = m[, "outcome"], prob = m[, "prob"]
      )
    })
  })
}

# Moments of one reconstructed distribution. Definitions mirror
# R/cpc18_features.R exactly (E[loss magnitude] = sum p * max(-x, 0);
# has_loss = any outcome in the support is negative).
dist_moments <- function(outcome, prob) {
  ev <- sum(outcome * prob)
  v <- sum(prob * outcome^2) - ev^2
  tibble::tibble(
    ev = ev,
    sd = sqrt(pmax(v, 0)),
    eloss = sum(prob * pmax(-outcome, 0)),
    has_loss = any(outcome < 0),
    n_outcomes = length(outcome)
  )
}

# Expanded-pool problem features (all Amb == 0 problems, LotNum unrestricted),
# from reconstructed distributions. Same columns and conventions as
# clean_problem_features() + add_scaled_features() in R/cpc18_features.R, with
# the moments computed on the full reconstructed support and the scaling SDs
# taken across the EXPANDED problem set (documented dataset-level constants;
# monotone-linear per feature, so they cannot affect rank statistics).
expanded_problem_features <- function(problems, dists) {
  mom <- dists |>
    dplyr::group_by(GameID, option) |>
    dplyr::summarise(dist_moments(outcome, prob), .groups = "drop")
  a <- mom |> dplyr::filter(option == "A")
  b <- mom |> dplyr::filter(option == "B")
  stopifnot(identical(a$GameID, problems$GameID),
            identical(b$GameID, problems$GameID))

  out <- problems |>
    dplyr::mutate(
      ev_a = a$ev, sd_a = a$sd, eloss_a = a$eloss, has_loss_a = a$has_loss,
      ev_b = b$ev, sd_b = b$sd, eloss_b = b$eloss, has_loss_b = b$has_loss,
      d_ev = ev_b - ev_a,
      d_sd = sd_b - sd_a,
      d_eloss = eloss_b - eloss_a,
      loss_informative = xor(has_loss_a, has_loss_b)
    )
  scales <- c(
    d_ev = stats::sd(out$d_ev),
    d_sd = stats::sd(out$d_sd),
    d_eloss = stats::sd(out$d_eloss)
  )
  out <- out |>
    dplyr::mutate(
      z_ev = d_ev / scales[["d_ev"]],
      z_sd = d_sd / scales[["d_sd"]],
      z_eloss = d_eloss / scales[["d_eloss"]]
    )
  attr(out, "scales") <- scales
  out
}

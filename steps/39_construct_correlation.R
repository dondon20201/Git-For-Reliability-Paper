# Step 39: rank correlation between loss avoidance and EV sensitivity
# Reads:   output/elevation/results/e0/per_subject_estimates.csv (step 18)
# Writes:  output/elevation/results/audit/construct_correlation.csv
# Paper:   Methods (Constructs)
# Runtime: a few seconds

# Compute the association between loss_avoid_adj and ev_sens_agn across
# participants with both full-data point estimates. Input is the step 18
# E0 table for the experience regime (r = 20, clean pool).
#
# Report Spearman correlation as the primary statistic and Pearson as a
# secondary statistic, with n and a 95% percentile participant-bootstrap
# interval for each. This adds one results file and changes no method or
# threshold in earlier steps.
#
# Bootstrap seed: 20261005, recorded in the output. There is no other
# random sampling.

suppressPackageStartupMessages({
  library(here)
})

SEED_BOOT <- 20261005L
N_BOOT <- 2000L
X_CONSTRUCT <- "loss_avoid_adj"
Y_CONSTRUCT <- "ev_sens_agn"

in_file <- here("output", "elevation", "results", "e0", "per_subject_estimates.csv")
out_file <- here("output", "elevation", "results", "audit", "construct_correlation.csv")

# exact numeric reader (utils::read.csv round-trips the shortest-representation
# doubles readr writes; same reader as step 38)
pse <- utils::read.csv(in_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("SubjID", "construct", "point") %in% names(pse)),
          !anyDuplicated(pse[c("SubjID", "construct")]))

pick <- function(cn) {
  d <- pse[pse$construct == cn, c("SubjID", "point")]
  stopifnot(nrow(d) == 686L)
  names(d)[2] <- cn
  d
}
wide <- merge(pick(X_CONSTRUCT), pick(Y_CONSTRUCT), by = "SubjID")
stopifnot(nrow(wide) == 686L)
ok <- stats::complete.cases(wide[c(X_CONSTRUCT, Y_CONSTRUCT)])
wide <- wide[ok, ]
x <- wide[[X_CONSTRUCT]]
y <- wide[[Y_CONSTRUCT]]
n <- length(x)

spearman <- stats::cor(x, y, method = "spearman")
pearson <- stats::cor(x, y, method = "pearson")

set.seed(SEED_BOOT)
boots <- t(replicate(N_BOOT, {
  i <- sample.int(n, n, replace = TRUE)
  c(stats::cor(x[i], y[i], method = "spearman"),
    stats::cor(x[i], y[i], method = "pearson"))
}))
ci <- function(v) stats::quantile(v, c(0.025, 0.975), names = FALSE)
ci_s <- ci(boots[, 1])
ci_p <- ci(boots[, 2])

out <- data.frame(
  pair = paste0(X_CONSTRUCT, "~", Y_CONSTRUCT),
  x = X_CONSTRUCT,
  y = Y_CONSTRUCT,
  n = n,
  n_excluded_missing = sum(!ok),
  spearman = spearman,
  spearman_ci_lo = ci_s[1],
  spearman_ci_hi = ci_s[2],
  pearson = pearson,
  pearson_ci_lo = ci_p[1],
  pearson_ci_hi = ci_p[2],
  n_boot = N_BOOT,
  boot_seed = SEED_BOOT,
  bootstrap = "participant resample, percentile 95%",
  source = "elevation/results/e0/per_subject_estimates.csv (E0 full-data points)",
  stringsAsFactors = FALSE
)
readr::write_csv(out, out_file, na = "NA")

# verification pass: re-read from disk and recompute the point values
back <- utils::read.csv(out_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(back) == 1L, back$n == n,
          isTRUE(all.equal(back$spearman, spearman, tolerance = 0)),
          isTRUE(all.equal(back$pearson, pearson, tolerance = 0)))

cat(sprintf("audit_construct_correlation OK: n = %d; Spearman %.4f [%.4f, %.4f]; Pearson %.4f [%.4f, %.4f]\n",
            n, spearman, ci_s[1], ci_s[2], pearson, ci_p[1], ci_p[2]))

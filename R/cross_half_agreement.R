# Cross-half vs same-half agreement between an agnostic index and a
# model-based coefficient (used by step 06). REAL DATA ONLY.
# Reuses the split machinery in R/step2_reliability.R (run_cell /
# split_correlations) unchanged; this file only adds the pairing metric.
#
# Cross-half: agnostic estimated on one half, model coefficient on the OTHER
# half, so shared sampling noise is removed; both directions (agn h1 ~ mod h2
# and agn h2 ~ mod h1) are computed and averaged. Same-half (agn h1 ~ mod h1,
# agn h2 ~ mod h2) is the contrast: it retains the correlated sampling error
# that the full-data agreement in step 04 (computed on all trials at once)
# contained.
# Per-split Spearman on complete pairs (min 10), averaged over splits and
# directions; subject-cluster bootstrap OUTSIDE, splits INSIDE (same nesting
# as step 05). NO Spearman-Brown correction and NO disattenuation anywhere:
# these are raw agreement correlations between half-budget measures.

# A1/A2: subject x split matrices of the agnostic index on half 1 / half 2;
# M1/M2: same for the model coefficient. Rows aligned (same subjects).
cross_half_point <- function(A1, A2, M1, M2, n_boot, seed, min_pairs = 10) {
  set.seed(seed)
  series <- function(a1, a2, m1, m2) {
    list(
      cross_d1 = split_correlations(a1, m2, min_pairs),  # agn h1 ~ mod h2
      cross_d2 = split_correlations(a2, m1, min_pairs),  # agn h2 ~ mod h1
      same_d1  = split_correlations(a1, m1, min_pairs),  # agn h1 ~ mod h1
      same_d2  = split_correlations(a2, m2, min_pairs)   # agn h2 ~ mod h2
    )
  }
  pt <- series(A1, A2, M1, M2)
  n_subj <- nrow(A1)

  boots <- replicate(n_boot, {
    idx <- sample.int(n_subj, replace = TRUE)
    b <- series(A1[idx, , drop = FALSE], A2[idx, , drop = FALSE],
                M1[idx, , drop = FALSE], M2[idx, , drop = FALSE])
    c(cross = mean(c(b$cross_d1, b$cross_d2), na.rm = TRUE),
      same  = mean(c(b$same_d1,  b$same_d2),  na.rm = TRUE))
  })

  row_for <- function(name, v1, v2, boot_row) {
    vals <- c(v1, v2)
    ci <- stats::quantile(boot_row, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    tibble::tibble(
      measure = name,
      rho = mean(vals, na.rm = TRUE),
      ci_lo = ci[1], ci_hi = ci[2],
      rho_dir1 = mean(v1, na.rm = TRUE), rho_dir2 = mean(v2, na.rm = TRUE),
      n_subjects = n_subj,
      n_split_dirs_used = sum(!is.na(vals)), n_split_dirs = length(vals),
      n_boot = n_boot
    )
  }
  dplyr::bind_rows(
    row_for("cross_half", pt$cross_d1, pt$cross_d2, boots["cross", ]),
    row_for("same_half",  pt$same_d1,  pt$same_d2,  boots["same", ])
  )
}

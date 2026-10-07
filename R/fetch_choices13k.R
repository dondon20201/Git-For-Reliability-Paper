# Fetch and md5-verify the choices13k dataset (secondary, robustness only).
#
# Source: https://github.com/jcpeterson/choices13k (Peterson, Bourgin,
# Agrawal, Reichman & Griffiths, 2021, Science 372:1209-1214), pinned by
# commit hash; files land in data/raw/choices13k/ (gitignored, like CPC18).
#
# Structure warning: choices13k is PROBLEM-LEVEL
# AGGREGATED data (one bRate per problem x feedback row; no per-subject
# repeated structure). It CANNOT support per-subject reliability claims and is
# used only for aggregate-level agreement checks.

choices13k_source <- list(
  repo   = "https://github.com/jcpeterson/choices13k",
  commit = "821ae7e88386b508ebb46fae76fac63cb62ec876",
  commit_date = "2026-03-23",
  citation = "Peterson et al. (2021), Science 372:1209-1214, doi:10.1126/science.abe2629"
)

choices13k_files <- data.frame(
  filename = c("c13k_selections.csv", "c13k_problems.json", "README.md"),
  url = sprintf("https://raw.githubusercontent.com/jcpeterson/choices13k/%s/%s",
                choices13k_source$commit,
                c("c13k_selections.csv", "c13k_problems.json", "README.md")),
  md5_expected = c(
    "65fd730812014d115c2cd0579f3122f9",
    "cc86e4e845dc2f48ff07c13b4663c20b",
    "da9da2743206c7a1bc8517a84bb5614d"
  ),
  stringsAsFactors = FALSE
)

# Download (if absent or failing verification) and md5-check every file.
# Errors out on mismatch: never analyze an unverified file.
fetch_choices13k <- function(dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout")
  options(timeout = 600)
  on.exit(options(timeout = old_timeout), add = TRUE)

  out <- lapply(seq_len(nrow(choices13k_files)), function(i) {
    f <- choices13k_files[i, ]
    dest <- file.path(dest_dir, f$filename)
    action <- "cached"
    if (!file.exists(dest) || unname(tools::md5sum(dest)) != f$md5_expected) {
      download.file(f$url, dest, mode = "wb", quiet = TRUE)
      action <- "downloaded"
    }
    md5_observed <- unname(tools::md5sum(dest))
    if (md5_observed != f$md5_expected) {
      stop(sprintf(
        paste0("md5 MISMATCH for choices13k file '%s': expected %s, observed %s. ",
               "STOPPING: do not analyze an unverified file."),
        f$filename, f$md5_expected, md5_observed
      ), call. = FALSE)
    }
    data.frame(
      file = f$filename, action = action, url = f$url,
      md5_expected = f$md5_expected, md5_observed = md5_observed,
      verified = TRUE, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

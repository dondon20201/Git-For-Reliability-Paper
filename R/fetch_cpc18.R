# Fetch and md5-verify the CPC18 estimation dataset (Zenodo record 845873, CC-BY).
#
# Files:
#   - "All estimation raw data.csv" (~47.5 MB): trial-level choices,
#     686 subjects x 30 problems x 25 repetitions (expected).
#   - "Dictionary for the choice prediction competition 2018 raw data.pdf":
#     the authoritative column/encoding reference. Read it before computing EVs.
#
# Raw data lands in data/raw/ (gitignored, read-only after download).
# On md5 mismatch the fetch STOPS with an error: never analyze an unverified file.

cpc18_files <- data.frame(
  filename = c(
    "All estimation raw data.csv",
    "Dictionary for the choice prediction competition 2018 raw data.pdf"
  ),
  url = c(
    "https://zenodo.org/records/845873/files/All%20estimation%20raw%20data.csv?download=1",
    "https://zenodo.org/records/845873/files/Dictionary%20for%20the%20choice%20prediction%20competition%202018%20raw%20data.pdf?download=1"
  ),
  md5_expected = c(
    "38a26407a36c690be969f9d9d3fe3e5a",
    "3980a6da798b4883167ebfcc7a2c52c5"
  ),
  stringsAsFactors = FALSE
)

# Download (if not already present and verified) and md5-check both files.
# Returns a data.frame with one row per file: action taken and md5 result.
# Errors out on any mismatch so downstream code cannot run on unverified data.
fetch_cpc18 <- function(dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout")
  options(timeout = 600)
  on.exit(options(timeout = old_timeout), add = TRUE)

  out <- lapply(seq_len(nrow(cpc18_files)), function(i) {
    f <- cpc18_files[i, ]
    dest <- file.path(dest_dir, f$filename)
    action <- "cached"
    if (!file.exists(dest) || unname(tools::md5sum(dest)) != f$md5_expected) {
      download.file(f$url, dest, mode = "wb", quiet = TRUE)
      action <- "downloaded"
    }
    md5_observed <- unname(tools::md5sum(dest))
    if (md5_observed != f$md5_expected) {
      stop(sprintf(
        paste0("md5 MISMATCH for '%s': expected %s, observed %s. ",
               "STOPPING: do not analyze an unverified file."),
        f$filename, f$md5_expected, md5_observed
      ), call. = FALSE)
    }
    data.frame(
      file = f$filename, action = action,
      md5_expected = f$md5_expected, md5_observed = md5_observed,
      verified = TRUE, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

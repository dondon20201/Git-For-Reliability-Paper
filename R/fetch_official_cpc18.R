# Fetch and md5-verify the OFFICIAL CPC18 baseline (BEASTsd) implementation
# used by the reconstruction validation in step 09.
#
# Provenance: the canonical distribution channel is the zip downloads linked at
# https://cpc-18.com/baseline-models-and-source-code/, which are hosted on a
# login-gated Google Sites page and cannot be fetched non-interactively. The
# files below are a verbatim copy of the official "CPC18 BEASTsd R" zip,
# committed to GitHub in April 2018 (contemporaneous with the competition) by
# an independent lab, pinned here by commit hash. A second, independent mirror
# (a Python port of CPC18_getDist in a different repository) is fetched as a
# cross-check that the algorithm is the same across mirrors; the official
# organizer-produced full-distribution table for the 210-problem estimation
# set (CPC18_EstSet_fullDist.csv) provides a third, data-level target.
# On md5 mismatch the fetch STOPS: never validate against unverified code.

official_cpc18_source <- list(
  channel   = "https://cpc-18.com/baseline-models-and-source-code/",
  repo      = "https://github.com/naecker-lab/cpc18",
  commit    = "b9caaca5874e3e388aea6580491ac639370b5bbd",
  commit_date = "2018-04-27",
  subdir    = "src/CPC18 BEASTsd R",
  crosscheck_repo   = "https://github.com/appgpk/choice_prediction",
  crosscheck_commit = "c2e7c571e7f0d7ca44c121265ae26a2d9df7e533"
)

official_cpc18_files <- data.frame(
  filename = c(
    "CPC18_getDist.R",
    "CPC18_get_pBetter.R",
    "distSample.R",
    "CPC18_BEASTsd_simulation.R",
    "CPC18_BEASTsd_pred.R",
    "CPC18_EstSet210.csv",
    "CPC18_EstSet_fullDist.csv",
    "GetDistribution_appgpk.py"
  ),
  url = c(
    sprintf("https://raw.githubusercontent.com/naecker-lab/cpc18/%s/src/CPC18%%20BEASTsd%%20R/%s",
            official_cpc18_source$commit,
            c("CPC18_getDist.R", "CPC18_get_pBetter.R", "distSample.R",
              "CPC18_BEASTsd_simulation.R", "CPC18_BEASTsd_pred.R",
              "CPC18_EstSet210.csv", "CPC18_EstSet_fullDist.csv")),
    sprintf("https://raw.githubusercontent.com/appgpk/choice_prediction/%s/GetDistribution.py",
            official_cpc18_source$crosscheck_commit)
  ),
  md5_expected = c(
    "338d7590949e11b598b5a8a081e4b79f",
    "b7a70b8b15977cccfc71dc0aeb60eaba",
    "389378a1febe93a1dbadced1a9fd57fa",
    "c042c61d6149c5d939d4c9fd72a3bff9",
    "cd8d3430e9d473ee956c01173e550543",
    "6c0817761e072ea4d642d7ff5909f839",
    "f09bb9b51fcb2bb641a1b10a9ad6bad8",
    "deab1ce2b5d2db8c18311201f396ecc3"
  ),
  stringsAsFactors = FALSE
)

# Download (if absent or failing verification) and md5-check every file.
# Returns one row per file with the action taken and md5 result; errors out on
# any mismatch so the validation gate cannot run against unverified code.
fetch_official_cpc18 <- function(dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout")
  options(timeout = 600)
  on.exit(options(timeout = old_timeout), add = TRUE)

  out <- lapply(seq_len(nrow(official_cpc18_files)), function(i) {
    f <- official_cpc18_files[i, ]
    dest <- file.path(dest_dir, f$filename)
    action <- "cached"
    if (!file.exists(dest) || unname(tools::md5sum(dest)) != f$md5_expected) {
      download.file(f$url, dest, mode = "wb", quiet = TRUE)
      action <- "downloaded"
    }
    md5_observed <- unname(tools::md5sum(dest))
    if (md5_observed != f$md5_expected) {
      stop(sprintf(
        paste0("md5 MISMATCH for official file '%s': expected %s, observed %s. ",
               "STOPPING: do not validate against unverified code."),
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

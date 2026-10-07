# Compare every file in expected/ with the file of the same name in output/.
#
# Usage (from the repository root): Rscript check.R
#
# CSV files are compared cell by cell and JSON files leaf by leaf, by position
# (CSV: same columns, same rows in the same order; JSON: same key path).
#   - A cell that is a number in expected/ must be a number in output/ and agree
#     within 1e-10 * max(1, |expected|). Inf, -Inf and NaN must match exactly.
#   - A missing cell (NA, empty, null) must not have become a number, and a
#     number must not have become missing.
#   - TRUE/FALSE flags must match exactly.
#   - All other text is ignored.
# PNG files are not checked. One line is printed per file, then a summary; the
# script exits with an error if any file differs or is missing.

suppressPackageStartupMessages(library(jsonlite))

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
}
root <- script_dir()
exp_dir <- file.path(root, "expected")
out_dir <- file.path(root, "output")
TOL <- 1e-10

is_missing <- function(x) x %in% c("NA", "", "null")
is_flag <- function(x) x %in% c("TRUE", "FALSE")
as_num <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  v[x == "NaN"] <- NaN
  v
}
is_num <- function(x) !is.na(as_num(x)) | x == "NaN"

# Compare two character vectors of cells; returns list(n_num, n_diff, max_diff, where)
compare_cells <- function(e, o, labels) {
  n_num <- 0L; n_diff <- 0L; max_diff <- 0; where <- ""
  en <- is_num(e); on <- is_num(o)
  ev <- as_num(e); ov <- as_num(o)
  for (i in seq_along(e)) {
    bad <- FALSE; d <- 0
    if (en[i]) {
      n_num <- n_num + 1L
      if (!on[i]) {
        bad <- TRUE; d <- Inf
      } else if (is.nan(ev[i]) || is.infinite(ev[i])) {
        bad <- !identical(ev[i], ov[i]); if (bad) d <- Inf
      } else if (!is.finite(ov[i])) {
        bad <- TRUE; d <- Inf
      } else {
        d <- abs(ev[i] - ov[i])
        bad <- d > TOL * max(1, abs(ev[i]))
      }
    } else if (is_missing(e[i])) {
      bad <- on[i]; if (bad) d <- Inf
    } else if (is_flag(e[i])) {
      bad <- !identical(e[i], o[i]); if (bad) d <- Inf
    }
    if (bad) {
      n_diff <- n_diff + 1L
      if (d >= max_diff) { max_diff <- d; where <- labels[i] }
    } else if (en[i] && d > max_diff && n_diff == 0L) {
      max_diff <- d
    }
  }
  list(n_num = n_num, n_diff = n_diff, max_diff = max_diff, where = where)
}

read_cells_csv <- function(path) {
  utils::read.csv(path, colClasses = "character", na.strings = character(0),
                  check.names = FALSE, stringsAsFactors = FALSE,
                  encoding = "UTF-8", strip.white = FALSE)
}

check_csv <- function(e_path, o_path) {
  e <- read_cells_csv(e_path); o <- read_cells_csv(o_path)
  if (!identical(names(e), names(o)))
    return(list(ok = FALSE, msg = "column names differ"))
  if (nrow(e) != nrow(o))
    return(list(ok = FALSE, msg = sprintf("%d rows expected, %d found", nrow(e), nrow(o))))
  labels <- as.vector(outer(seq_len(nrow(e)), names(e),
                            function(r, c) paste0("row ", r, ", column ", c)))
  r <- compare_cells(unlist(e, use.names = FALSE), unlist(o, use.names = FALSE), labels)
  list(ok = r$n_diff == 0L, r = r)
}

# Flatten parsed JSON into a named character vector of leaves (path -> text).
flatten_json <- function(x, path = "") {
  if (is.null(x)) return(stats::setNames("null", path))
  if (is.list(x)) {
    if (length(x) == 0L) return(character(0))
    nm <- names(x)
    if (is.null(nm)) nm <- paste0("[", seq_along(x), "]")
    out <- character(0)
    for (i in seq_along(x)) {
      key <- if (startsWith(nm[i], "[")) paste0(path, nm[i]) else
        if (nzchar(path)) paste0(path, ".", nm[i]) else nm[i]
      out <- c(out, flatten_json(x[[i]], key))
    }
    return(out)
  }
  if (length(x) != 1L) {
    return(unlist(lapply(seq_along(x), function(i) flatten_json(x[i], paste0(path, "[", i, "]")))))
  }
  txt <- if (is.logical(x)) { if (is.na(x)) "NA" else if (x) "TRUE" else "FALSE" }
         else if (is.numeric(x)) format(x, digits = 17) else paste0("§", x)
  stats::setNames(txt, path)
}

check_json <- function(e_path, o_path) {
  e <- flatten_json(fromJSON(e_path, simplifyVector = FALSE))
  o <- flatten_json(fromJSON(o_path, simplifyVector = FALSE))
  o_al <- o[names(e)]
  o_al[is.na(o_al)] <- "§<absent>"
  r <- compare_cells(unname(e), unname(o_al), names(e))
  list(ok = r$n_diff == 0L, r = r)
}

files <- sort(list.files(exp_dir, recursive = TRUE))
files <- files[!grepl("\\.png$", files, ignore.case = TRUE)]
if (length(files) == 0L) stop("no files in expected/")

n_ok <- 0L; n_bad <- 0L
width <- max(nchar(files))
for (f in files) {
  o_path <- file.path(out_dir, f)
  if (!file.exists(o_path)) {
    cat(sprintf("%-*s  MISSING  (not in output/)\n", width, f)); n_bad <- n_bad + 1L; next
  }
  res <- tryCatch(
    if (grepl("\\.json$", f)) check_json(file.path(exp_dir, f), o_path)
    else check_csv(file.path(exp_dir, f), o_path),
    error = function(err) list(ok = FALSE, msg = paste("could not compare:", conditionMessage(err))))
  if (!is.null(res$msg)) {
    cat(sprintf("%-*s  DIFFERS  (%s)\n", width, f, res$msg)); n_bad <- n_bad + 1L
  } else if (res$ok) {
    cat(sprintf("%-*s  OK       (%d numbers, largest difference %.3g)\n",
                width, f, res$r$n_num, res$r$max_diff)); n_ok <- n_ok + 1L
  } else {
    cat(sprintf("%-*s  DIFFERS  (%d of %d cells; largest difference %.3g at %s)\n",
                width, f, res$r$n_diff, res$r$n_num, res$r$max_diff, res$r$where))
    n_bad <- n_bad + 1L
  }
}
cat(sprintf("\n%d files checked: %d OK, %d differ or missing\n", length(files), n_ok, n_bad))
if (n_bad > 0L) quit(status = 1L)

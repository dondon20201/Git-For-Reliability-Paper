# Run every analysis step in order with Rscript, then check.R.
#
# Usage (from the repository root):
#   Rscript run_all.R          all steps, then check.R
#   Rscript run_all.R 1 10     steps 01 to 10 only
#   Rscript run_all.R 22       step 22 to the last step, then check.R
#
# Each step runs in its own fresh R process. Its console output goes to
# output/logs/<step>.log. The run stops at the first step that fails.
# check.R runs at the end whenever the range includes the last step.

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
}
root <- script_dir()
setwd(root)

steps <- sort(list.files("steps", pattern = "^[0-9]{2}_.*\\.R$"))
if (length(steps) == 0L) stop("no step files found in steps/")
nums <- as.integer(substr(steps, 1, 2))

args <- commandArgs(trailingOnly = TRUE)
from <- if (length(args) >= 1L) as.integer(args[1]) else min(nums)
to   <- if (length(args) >= 2L) as.integer(args[2]) else max(nums)
if (is.na(from) || is.na(to) || from > to) stop("usage: Rscript run_all.R [first_step] [last_step]")
todo <- steps[nums >= from & nums <= to]

rscript <- file.path(R.home("bin"), "Rscript")
log_dir <- file.path("output", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
timing_file <- file.path(log_dir, "run_times.csv")
if (!file.exists(timing_file))
  writeLines("step,status,start,end,minutes", timing_file)

t_all <- Sys.time()
for (s in todo) {
  log <- file.path(log_dir, sub("\\.R$", ".log", s))
  t0 <- Sys.time()
  cat(sprintf("[%s] %-40s ", format(t0, "%Y-%m-%d %H:%M:%S"), s))
  status <- system2(rscript, file.path("steps", s), stdout = log, stderr = log)
  t1 <- Sys.time()
  mins <- as.numeric(difftime(t1, t0, units = "mins"))
  # A step that draws a plot without saving it leaves Rplots.pdf in the
  # working directory; keep it with the step's log instead.
  if (file.exists("Rplots.pdf"))
    file.rename("Rplots.pdf", file.path(log_dir, sub("\\.R$", "_Rplots.pdf", s)))
  ok <- identical(as.integer(status), 0L)
  cat(sprintf("%s  %.1f min\n", if (ok) "OK" else "FAILED", mins))
  cat(sprintf("%s,%s,%s,%s,%.2f\n", s, if (ok) "OK" else "FAILED",
              format(t0, "%Y-%m-%d %H:%M:%S"), format(t1, "%Y-%m-%d %H:%M:%S"), mins),
      file = timing_file, append = TRUE)
  if (!ok) {
    cat("Stopped: ", s, " failed (exit status ", status, "). See ", log, "\n", sep = "")
    quit(status = 1L)
  }
}
cat(sprintf("Steps %02d to %02d finished in %.1f hours.\n", from, to,
            as.numeric(difftime(Sys.time(), t_all, units = "hours"))))

if (to == max(nums)) {
  cat("\nRunning check.R\n")
  status <- system2(rscript, "check.R")
  quit(status = as.integer(status))
} else {
  cat("Run the remaining steps, then Rscript check.R\n")
}

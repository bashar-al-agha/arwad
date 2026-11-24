# Phase 2, CVPAT: where aw_assess_cvpat() spends its time, on the model and
# settings of pls_cvpat_rounds.R (10 folds, 10 repetitions, earliest
# antecedents, 2,000 bootstrap resamples, seed 123). Medians of five runs of
# (a) the cross-validated losses alone, (b) the whole call and (c) the whole
# call with 2 resamples; the cross-validation's share is a / b and the
# bootstrap's (b - c) / b. Shares only: absolute times depend on the
# processor's state (REPORT.md, section 7).
# Usage (from the repository root, with ARWAD_GOOSE set):
#   Rscript bench/phase2/pls_cvpat_split.R
# Writes bench/phase2/pls_cvpat_split.csv; REPORT.md, section 3f, uses it.

suppressPackageStartupMessages({
  library(arwad)
  library(seminr)
})
options(arwad.cache_dir = file.path(tempdir(), "arwad-bench-cache"))
stopifnot(aw_toolchain()$ok)

mm <- constructs(
  composite("QUAL", multi_items("qual_", 1:8), weights = mode_B),
  composite("PERF", multi_items("perf_", 1:5), weights = mode_B),
  composite("CSOR", multi_items("csor_", 1:5), weights = mode_B),
  composite("ATTR", multi_items("attr_", 1:3), weights = mode_B),
  composite("COMP", multi_items("comp_", 1:3)),
  composite("LIKE", multi_items("like_", 1:3)),
  composite("CUSA", single_item("cusa")),
  composite("CUSL", multi_items("cusl_", 1:3)))
sm <- relationships(
  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE")),
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = "CUSA", to = "CUSL"))
model <- suppressMessages(estimate_pls(corp_rep_data, mm, sm, maxIt = 3000L))
seed <- 123L

# Not timed: builds the PLS engine and the bootstrap unit.
invisible(aw_assess_cvpat(model, nboot = 20L, seed = seed, technique = predict_EA, noFolds = 10L, reps = 2L))
med <- function(f) stats::median(replicate(5L, system.time(f())[["elapsed"]]))
t_cv <- med(function() {
  set.seed(seed)
  stopifnot(!is.null(arwad:::aw_cvpat_native_losses(model, predict_EA, 10L, 10L, 1L)))
})
t_all <- med(function() aw_assess_cvpat(model, nboot = 2000L, seed = seed, technique = predict_EA, noFolds = 10L, reps = 10L))
t_b2 <- med(function() aw_assess_cvpat(model, nboot = 2L, seed = seed, technique = predict_EA, noFolds = 10L, reps = 10L))
res <- data.frame(part = c("cross-validated losses", "whole call, 2000 resamples", "whole call, 2 resamples"),
                  median_s = round(c(t_cv, t_all, t_b2), 3), stringsAsFactors = FALSE)
print(res, row.names = FALSE)
cat(sprintf("shares of the whole call: cross-validation %.0f%%, bootstrap of 1998 resamples %.0f%%\n",
            100 * t_cv / t_all, 100 * (t_all - t_b2) / t_all))
cat(R.version.string, "| arwad", as.character(packageVersion("arwad")), "| seminr", as.character(packageVersion("seminr")), "|",
    utils::sessionInfo()$running, "\n")
write.csv(res, file.path("bench", "phase2", "pls_cvpat_split.csv"), row.names = FALSE)

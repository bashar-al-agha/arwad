# Phase 2, PLSpredict: seminr::predict_pls() against aw_predict_pls() on the
# Primer's extended corporate reputation model (four Mode B drivers, n =
# 344, maxIt = 3000), with the settings of the project lead's protocol: 10
# folds, 10 repetitions, earliest antecedents; and leave-one-out
# cross-validation of the same model. Three rounds, seminr first in each,
# memory freed between runs; seminr runs sequentially (cores = NULL, its
# default) and on all cores, Arwad on one core and on all cores. Every Arwad
# result is checked against seminr's of the same round (same seed, so the
# same folds). Figures: seminr's fastest run against Arwad's median
# (ADR-0006).
# Usage (from the repository root, with ARWAD_GOOSE set):
#   Rscript bench/phase2/pls_predict_rounds.R
# Writes bench/phase2/pls_predict_rounds.csv; REPORT.md interprets it.

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
cores <- parallel::detectCores()
set.seed(seed)
stopifnot(aw_predict_pls(model, technique = predict_EA, noFolds = 10L, reps = 2L)$native)  # builds the engine and the averaging unit (not timed)

same <- function(a, s) isTRUE(a$native) && identical(unclass(a)[names(s)], unclass(s), num.eq = FALSE)
time_run <- function(f) {
  set.seed(seed)
  t <- system.time(out <- f())[["elapsed"]]
  list(t = t, out = out)
}
designs <- list(kfold_10x10 = list(noFolds = 10L, reps = 10L), loocv = list(noFolds = NULL, reps = NULL))
res <- NULL
for (r in 1:3) for (d in names(designs)) {
  a <- designs[[d]]
  s1 <- time_run(function() suppressWarnings(predict_pls(model, technique = predict_EA, noFolds = a$noFolds, reps = a$reps)))
  invisible(gc())
  sc <- time_run(function() suppressWarnings(predict_pls(model, technique = predict_EA, noFolds = a$noFolds, reps = a$reps, cores = cores)))
  samesc <- identical(unclass(sc$out), unclass(s1$out), num.eq = FALSE)
  sc$out <- NULL
  invisible(gc())
  a1 <- time_run(function() aw_predict_pls(model, technique = predict_EA, noFolds = a$noFolds, reps = a$reps))
  same1 <- same(a1$out, s1$out)
  a1$out <- NULL
  invisible(gc())
  ac <- time_run(function() aw_predict_pls(model, technique = predict_EA, noFolds = a$noFolds, reps = a$reps, cores = cores))
  samec <- same(ac$out, s1$out)
  res <- rbind(res, data.frame(round = r, design = d, seminr_s = round(s1$t, 2), seminr_all_cores_s = round(sc$t, 2),
                               arwad_s = round(a1$t, 2), arwad_all_cores_s = round(ac$t, 2), identical = samesc && same1 && samec,
                               stringsAsFactors = FALSE))
  rm(s1, sc, a1, ac)
  invisible(gc())
  print(res[nrow(res), ], row.names = FALSE)
}
for (d in names(designs)) {
  x <- res[res$design == d, ]
  fast <- min(x$seminr_s, x$seminr_all_cores_s)
  m1 <- stats::median(x$arwad_s)
  mc <- stats::median(x$arwad_all_cores_s)
  cat(sprintf(paste0("\n%s, three rounds: seminr fastest %.2f s (sequential median %.2f, %d cores median %.2f); ",
                     "Arwad median %.2f s on one core, %.2f s on %d cores; seminr fastest / Arwad median: %.1f (one core), %.1f (%d cores); ",
                     "identical results in every run: %s\n"),
              d, fast, stats::median(x$seminr_s), cores, stats::median(x$seminr_all_cores_s), m1, mc, cores, fast / m1, fast / mc,
              cores, all(x$identical)))
}
cat(R.version.string, "| arwad", as.character(packageVersion("arwad")), "| seminr", as.character(packageVersion("seminr")), "|",
    utils::sessionInfo()$running, "\n")
write.csv(res, file.path("bench", "phase2", "pls_predict_rounds.csv"), row.names = FALSE)

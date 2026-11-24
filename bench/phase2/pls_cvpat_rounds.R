# Phase 2, CVPAT: seminrExtras::assess_cvpat() against aw_assess_cvpat() on
# the Primer's extended corporate reputation model (four Mode B drivers, n =
# 344, maxIt = 3000), with the settings of the project lead's protocol: 10
# folds, 10 repetitions, earliest antecedents; 2,000 bootstrap resamples for
# the t-tests and seed 123 (both seminrExtras' defaults). Three rounds,
# seminrExtras first in each, memory freed between runs; seminrExtras runs
# sequentially (cores = NULL, its default) and on all cores, Arwad on one core
# and on all cores. Every result is checked against seminrExtras' sequential
# result of the same round (same seed, so the same folds and resamples).
# Figures: seminrExtras' fastest run against Arwad's median (ADR-0006).
# Usage (from the repository root, with ARWAD_GOOSE set and seminrExtras
# installed):
#   Rscript bench/phase2/pls_cvpat_rounds.R
# Writes bench/phase2/pls_cvpat_rounds.csv; REPORT.md interprets it.

suppressPackageStartupMessages({
  library(arwad)
  library(seminr)
  library(seminrExtras)
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
nboot <- 2000L
folds <- 10L
reps <- 10L
cores <- parallel::detectCores()

# Not timed: Arwad's native path takes this model and these settings on the
# folds of this seed (aw_assess_cvpat() would otherwise hand the call to
# seminrExtras), which also builds the engine; a short run builds the compiled
# bootstrap unit, which must have no bailouts to GNU R.
set.seed(seed)
stopifnot(!is.null(arwad:::aw_cvpat_native_losses(model, predict_EA, folds, reps, 1L)))
invisible(aw_assess_cvpat(model, nboot = 20L, seed = seed, technique = predict_EA, noFolds = folds, reps = 2L))
cat(utils::capture.output(info <- aw_explain(arwad:::aw_cvpat_boot_unit())), sep = "\n")
stopifnot(sum(info$bailouts) == 0L)

# Both tables, their attributes and every column (the t-statistics and
# p-values included), bit for bit.
same <- function(a, s) {
  identical(names(a), names(s)) && all(vapply(names(s), function(nm)
    identical(attributes(a[[nm]]), attributes(s[[nm]])) &&
      identical(unclass(a[[nm]]), unclass(s[[nm]]), num.eq = FALSE), TRUE))
}
time_run <- function(f) {
  t <- system.time(out <- suppressWarnings(suppressMessages(f())))[["elapsed"]]
  list(t = t, out = out)
}
res <- NULL
for (r in 1:3) {
  s1 <- time_run(function() assess_cvpat(model, nboot = nboot, seed = seed, technique = predict_EA,
                                         noFolds = folds, reps = reps))
  invisible(gc())
  sc <- time_run(function() assess_cvpat(model, nboot = nboot, seed = seed, technique = predict_EA,
                                         noFolds = folds, reps = reps, cores = cores))
  samesc <- same(sc$out, s1$out)
  sc$out <- NULL
  invisible(gc())
  a1 <- time_run(function() aw_assess_cvpat(model, nboot = nboot, seed = seed, technique = predict_EA,
                                            noFolds = folds, reps = reps))
  same1 <- same(a1$out, s1$out)
  a1$out <- NULL
  invisible(gc())
  ac <- time_run(function() aw_assess_cvpat(model, nboot = nboot, seed = seed, technique = predict_EA,
                                            noFolds = folds, reps = reps, cores = cores))
  samec <- same(ac$out, s1$out)
  if (r == 1L) print(ac$out)
  ac$out <- NULL
  if (!(samesc && same1 && samec)) {
    cat(sprintf("round %d: seminrExtras on %d cores identical: %s; Arwad one core: %s; Arwad %d cores: %s\n",
                r, cores, samesc, same1, cores, samec))
  }
  res <- rbind(res, data.frame(round = r, seminrextras_s = round(s1$t, 2), seminrextras_all_cores_s = round(sc$t, 2),
                               arwad_s = round(a1$t, 2), arwad_all_cores_s = round(ac$t, 2),
                               identical = samesc && same1 && samec, stringsAsFactors = FALSE))
  rm(s1, sc, a1, ac)
  invisible(gc())
  print(res[nrow(res), ], row.names = FALSE)
}
fast <- min(res$seminrextras_s, res$seminrextras_all_cores_s)
m1 <- stats::median(res$arwad_s)
mc <- stats::median(res$arwad_all_cores_s)
cat(sprintf(paste0("\nCVPAT, %d folds x %d repetitions, %d resamples, three rounds: seminrExtras fastest %.2f s ",
                   "(sequential median %.2f, %d cores median %.2f); Arwad median %.2f s on one core, %.2f s on %d cores; ",
                   "seminrExtras fastest / Arwad median: %.1f (one core), %.1f (%d cores); identical results in every run: %s\n"),
            folds, reps, nboot, fast, stats::median(res$seminrextras_s), cores, stats::median(res$seminrextras_all_cores_s),
            m1, mc, cores, fast / m1, fast / mc, cores, all(res$identical)))
cat(R.version.string, "| arwad", as.character(packageVersion("arwad")), "| seminr", as.character(packageVersion("seminr")),
    "| seminrExtras", as.character(packageVersion("seminrExtras")), "|", utils::sessionInfo()$running, "\n")
write.csv(res, file.path("bench", "phase2", "pls_cvpat_rounds.csv"), row.names = FALSE)

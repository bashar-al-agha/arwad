# Phase 2, the PLS-SEM engine with a two-stage interaction term: 10,000
# resamples in three rounds, seminr first in each round, on the moderation
# model of pls_moderation_bench.R (same model, seed and settings). Single
# rounds of this model varied by up to 1.7 times between runs on the
# reference laptop (REPORT.md §3d), so the 10,000-resample figures come from
# here: seminr's fastest round against Arwad's median, as for the ECSI model
# (ADR-0006). Memory is freed between runs; every run is checked against
# seminr's arrays of the same round.
# Usage (from the repository root, with ARWAD_GOOSE set):
#   Rscript bench/phase2/pls_moderation_rounds.R
# Writes bench/phase2/pls_moderation_rounds.csv; REPORT.md interprets it.

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
  composite("SC", multi_items("switch_", 1:4)),
  composite("CUSL", multi_items("cusl_", 1:3)),
  interaction_term(iv = "CUSA", moderator = "SC", method = two_stage))
sm <- relationships(
  paths(from = c("QUAL", "PERF", "CSOR", "ATTR"), to = c("COMP", "LIKE")),
  paths(from = c("COMP", "LIKE"), to = c("CUSA", "CUSL")),
  paths(from = c("CUSA", "SC", "CUSA*SC"), to = "CUSL"))
model <- suppressMessages(estimate_pls(corp_rep_data, mm, sm, maxIt = 3000L))
B <- 10000L
seed <- 123L
cores <- parallel::detectCores()
exact <- c("boot_paths", "boot_loadings", "boot_weights", "boot_total_paths", "boots")
stopifnot(suppressMessages(aw_boot_pls(model, nboot = 10L, seed = seed))$native)  # builds the engine (not timed)

res <- NULL
for (r in 1:3) {
  ts <- system.time(s <- suppressWarnings(suppressMessages(bootstrap_model(model, nboot = B, cores = NULL, seed = seed))))[["elapsed"]]
  ref <- s[exact]
  rm(s); invisible(gc())
  t1 <- system.time(a <- suppressMessages(aw_boot_pls(model, nboot = B, seed = seed, cores = 1L)))[["elapsed"]]
  same1 <- a$native && all(vapply(exact, function(f) identical(a[[f]], ref[[f]], num.eq = FALSE), NA))
  rm(a); invisible(gc())
  tc <- system.time(a <- suppressMessages(aw_boot_pls(model, nboot = B, seed = seed, cores = cores)))[["elapsed"]]
  samec <- a$native && all(vapply(exact, function(f) identical(a[[f]], ref[[f]], num.eq = FALSE), NA))
  rm(a, ref); invisible(gc())
  res <- rbind(res, data.frame(round = r, seminr_all_cores_s = round(ts, 2), arwad_s = round(t1, 2),
                               arwad_all_cores_s = round(tc, 2), identical = same1 && samec))
  print(res[r, ], row.names = FALSE)
}
fast <- min(res$seminr_all_cores_s)
m1 <- stats::median(res$arwad_s)
mc <- stats::median(res$arwad_all_cores_s)
cat(sprintf(paste0("\n10,000 resamples, three rounds: seminr (all cores) fastest %.1f s, median %.1f s; ",
                   "Arwad median %.1f s on one core, %.1f s on %d cores; seminr fastest / Arwad median: %.1f (%d cores), %.1f (one core); ",
                   "identical arrays in every run: %s\n"),
            fast, stats::median(res$seminr_all_cores_s), m1, mc, cores, fast / mc, cores, fast / m1, all(res$identical)))
cat(R.version.string, "| arwad", as.character(packageVersion("arwad")), "| seminr", as.character(packageVersion("seminr")), "|",
    utils::sessionInfo()$running, "\n")
write.csv(res, file.path("bench", "phase2", "pls_moderation_rounds.csv"), row.names = FALSE)

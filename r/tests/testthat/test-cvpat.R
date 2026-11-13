# aw_assess_cvpat() and aw_assess_cvpat_compare() against seminrExtras (ADR-0008).

skip_if_old_seminr <- function() testthat::skip_if_not(arwad:::aw_seminr_260(), "seminr before 2.6.0")

# Exact structure, attributes (in pairlist order) and every column, T and P
# included, bit for bit.
expect_same_cvpat_table <- function(got, want, label) {
  expect_identical(attributes(got), attributes(want), label = paste(label, "attributes"))
  expect_true(identical(unclass(got), unclass(want), num.eq = FALSE), label = paste(label, "identical"))
}

expect_same_cvpat <- function(got, want, label) {
  expect_identical(names(got), names(want), label = paste(label, "names"))
  expect_same_cvpat_table(got$CVPAT_compare_LM, want$CVPAT_compare_LM, paste(label, "LM"))
  expect_same_cvpat_table(got$CVPAT_compare_IA, want$CVPAT_compare_IA, paste(label, "IA"))
}

test_that("aw_assess_cvpat() matches seminrExtras::assess_cvpat() on at least 50 models of composites", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_not_installed("seminrExtras")
  skip_if_old_seminr()
  fits <- aw_pls_fit_models()
  comp <- Filter(function(m) !any(m$mmMatrix[, "type"] == "C"), fits)
  expect_gte(length(comp), 50L)
  for (k in seq_along(comp)) {
    nm <- names(comp)[k]
    tech <- if (k %% 2L == 1L) seminr::predict_EA else seminr::predict_DA
    nf <- c(10L, 5L, 7L)[k %% 3L + 1L]
    tt <- if (k %% 4L == 0L) "greater" else "two.sided"
    # reps = NULL (the default) on two models in three, reps = 1 on the third.
    rp <- if (k %% 3L == 0L) 1L else NULL
    w_want <- testthat::capture_warnings(want <- suppressMessages(
      seminrExtras::assess_cvpat(comp[[nm]], testtype = tt, nboot = 20L, seed = k,
                                 technique = tech, noFolds = nf, reps = rp, cores = 1L)))
    sw <- .Random.seed
    w_got <- testthat::capture_warnings(got <- aw_assess_cvpat(comp[[nm]], testtype = tt, nboot = 20L, seed = k,
                                                               technique = tech, noFolds = nf, reps = rp, cores = 1L))
    expect_same_cvpat(got, want, nm)
    expect_identical(w_got, w_want, label = paste(nm, "warnings"))
    expect_identical(.Random.seed, sw, label = paste(nm, "random-number state"))
  }
  # Verify the compiled bootstrap kernel has zero bailouts.
  u <- arwad:::aw_cvpat_boot_unit()
  expect_s3_class(u, "arwad_function")
  utils::capture.output(info <- aw_explain(u))
  expect_identical(sum(info$bailouts), 0L)
})

test_that("aw_assess_cvpat() handles repetitions, leave-one-out, single indicator, cores, seed = NULL, and testtype = 'less'", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_not_installed("seminrExtras")
  skip_if_old_seminr()
  fits <- aw_pls_fit_models()
  ecsi <- fits[["mobi_ecsi_core"]]
  for (cfg in list(list(seminr::predict_EA, 10L, 2L, 1L, "two.sided"),
                   list(seminr::predict_DA, 5L, 3L, 1L, "greater"),
                   list(seminr::predict_EA, 5L, 2L, 2L, "two.sided"))) {
    want <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat(ecsi, testtype = cfg[[5L]], nboot = 25L, seed = 42L,
                                 technique = cfg[[1L]], noFolds = cfg[[2L]], reps = cfg[[3L]], cores = cfg[[4L]])))
    sw <- .Random.seed
    got <- aw_assess_cvpat(ecsi, testtype = cfg[[5L]], nboot = 25L, seed = 42L,
                           technique = cfg[[1L]], noFolds = cfg[[2L]], reps = cfg[[3L]], cores = cfg[[4L]])
    lab <- sprintf("reps=%d cores=%d %s", cfg[[3L]], cfg[[4L]], cfg[[5L]])
    expect_same_cvpat(got, want, lab)
    expect_identical(.Random.seed, sw, label = paste(lab, "random-number state"))
  }
  # seed = NULL consumes the caller's RNG stream without resetting.
  want_null <- suppressWarnings(suppressMessages(
    seminrExtras::assess_cvpat(ecsi, nboot = 20L, seed = 101L, noFolds = 5L, reps = 1L)))
  sw_null <- .Random.seed
  set.seed(101L)
  got_null <- aw_assess_cvpat(ecsi, nboot = 20L, seed = NULL, noFolds = 5L, reps = 1L)
  expect_same_cvpat(got_null, want_null, "seed = NULL")
  expect_identical(.Random.seed, sw_null)
  # Leave-one-out (noFolds = NULL) on the smallest composite model.
  comp <- Filter(function(m) !any(m$mmMatrix[, "type"] == "C"), fits)
  small <- comp[[which.min(vapply(comp, function(m) nrow(m$data), 1L))]]
  want_loo <- suppressWarnings(suppressMessages(
    seminrExtras::assess_cvpat(small, nboot = 15L, seed = 7L, noFolds = NULL, reps = 1L)))
  sw_loo <- .Random.seed
  got_loo <- aw_assess_cvpat(small, nboot = 15L, seed = 7L, noFolds = NULL, reps = 1L)
  expect_same_cvpat(got_loo, want_loo, "leave-one-out")
  expect_identical(.Random.seed, sw_loo)
  # Single-indicator endogenous construct.
  mm1 <- seminr::constructs(seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
                            seminr::composite("Satisfaction", seminr::single_item("CUSA1")))
  m1 <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm1, seminr::relationships(seminr::paths(from = "Image", to = "Satisfaction"))))
  for (R in c(1L, 2L)) {
    want1 <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat(m1, nboot = 20L, seed = 8L, noFolds = 5L, reps = R)))
    sw1 <- .Random.seed
    got1 <- aw_assess_cvpat(m1, nboot = 20L, seed = 8L, noFolds = 5L, reps = R)
    expect_same_cvpat(got1, want1, paste("single indicator reps", R))
    expect_identical(.Random.seed, sw1, label = paste("single indicator reps", R, "random-number state"))
  }
  # testtype = "less" (which crashes in seminrExtras 1.0.4 due to an unbound
  # variable): verify identical cols 1-4 and lower-tail t p-values.
  got_gt <- aw_assess_cvpat(ecsi, testtype = "greater", nboot = 30L, seed = 55L, noFolds = 5L, reps = 1L)
  got_lt <- aw_assess_cvpat(ecsi, testtype = "less", nboot = 30L, seed = 55L, noFolds = 5L, reps = 1L)
  for (nm in c("CVPAT_compare_LM", "CVPAT_compare_IA")) {
    expect_identical(got_lt[[nm]][, 1:4], got_gt[[nm]][, 1:4])
    expect_equal(unname(got_lt[[nm]][, 5L] + got_gt[[nm]][, 5L]), rep(1, nrow(got_lt[[nm]])), tolerance = 1e-6)
  }
})

test_that("aw_assess_cvpat_compare() matches seminrExtras::assess_cvpat_compare() and supports testtype = 'less'", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_not_installed("seminrExtras")
  skip_if_old_seminr()
  mm <- seminr::constructs(
    seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
    seminr::composite("Expectation", seminr::multi_items("CUEX", 1:3)),
    seminr::composite("Value", seminr::multi_items("PERV", 1:2)),
    seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3))
  )
  sm1 <- seminr::relationships(
    seminr::paths(from = "Image", to = c("Expectation", "Value", "Satisfaction")),
    seminr::paths(from = "Expectation", to = c("Value", "Satisfaction")),
    seminr::paths(from = "Value", to = "Satisfaction")
  )
  sm2 <- seminr::relationships(
    seminr::paths(from = "Image", to = c("Expectation", "Value")),
    seminr::paths(from = "Expectation", to = "Satisfaction"),
    seminr::paths(from = "Value", to = "Satisfaction")
  )
  m1 <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm, sm1))
  m2 <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm, sm2))
  for (cfg in list(list("two.sided", 1L, 1L, seminr::predict_DA, 5L), list("greater", 2L, 1L, seminr::predict_DA, 5L),
                   list("two.sided", 2L, 2L, seminr::predict_DA, 5L), list("two.sided", NULL, 1L, seminr::predict_EA, 10L))) {
    want <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat_compare(m1, m2, testtype = cfg[[1L]], nboot = 25L, seed = 17L,
                                         technique = cfg[[4L]], noFolds = cfg[[5L]], reps = cfg[[2L]], cores = cfg[[3L]])))
    sw <- .Random.seed
    got <- aw_assess_cvpat_compare(m1, m2, testtype = cfg[[1L]], nboot = 25L, seed = 17L,
                                   technique = cfg[[4L]], noFolds = cfg[[5L]], reps = cfg[[2L]], cores = cfg[[3L]])
    lab <- paste("compare", cfg[[1L]], "reps", if (is.null(cfg[[2L]])) "NULL" else cfg[[2L]], "cores", cfg[[3L]],
                 if (identical(cfg[[4L]], seminr::predict_EA)) "EA" else "DA", "folds", cfg[[5L]])
    expect_same_cvpat_table(got, want, lab)
    expect_identical(.Random.seed, sw, label = paste(lab, "random-number state"))
  }
  # testtype = "less" in aw_assess_cvpat_compare().
  gt <- aw_assess_cvpat_compare(m1, m2, testtype = "greater", nboot = 25L, seed = 17L, noFolds = 5L, reps = 1L)
  lt <- aw_assess_cvpat_compare(m1, m2, testtype = "less", nboot = 25L, seed = 17L, noFolds = 5L, reps = 1L)
  expect_identical(lt[, 1:4], gt[, 1:4])
  expect_equal(unname(lt[, 5L] + gt[, 5L]), rep(1, nrow(lt)), tolerance = 1e-6)
  # Reordered smMatrix target order, Mode A vs Mode B on endogenous construct,
  # and reordered indicator rows in mmMatrix all match seminrExtras.
  sm2_rev <- seminr::relationships(
    seminr::paths(from = "Image", to = c("Satisfaction", "Value", "Expectation")),
    seminr::paths(from = "Expectation", to = c("Satisfaction", "Value")),
    seminr::paths(from = "Value", to = "Satisfaction")
  )
  mm_modeB <- seminr::constructs(
    seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
    seminr::composite("Expectation", seminr::multi_items("CUEX", 1:3)),
    seminr::composite("Value", seminr::multi_items("PERV", 1:2), weights = seminr::mode_B),
    seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3))
  )
  mm_reord <- seminr::constructs(
    seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
    seminr::composite("Expectation", seminr::multi_items("CUEX", 1:3)),
    seminr::composite("Value", seminr::multi_items("PERV", 2:1)),
    seminr::composite("Satisfaction", seminr::multi_items("CUSA", 3:1))
  )
  alts <- list(
    rev_targets = suppressMessages(seminr::estimate_pls(seminr::mobi, mm, sm2_rev)),
    modeB_endo  = suppressMessages(seminr::estimate_pls(seminr::mobi, mm_modeB, sm2)),
    reord_items = suppressMessages(seminr::estimate_pls(seminr::mobi, mm_reord, sm2))
  )
  for (nm in names(alts)) {
    want_a <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat_compare(m1, alts[[nm]], nboot = 20L, seed = 33L, noFolds = 5L, reps = 1L)))
    sw_a <- .Random.seed
    got_a <- aw_assess_cvpat_compare(m1, alts[[nm]], nboot = 20L, seed = 33L, noFolds = 5L, reps = 1L)
    expect_same_cvpat_table(got_a, want_a, nm)
    expect_identical(.Random.seed, sw_a, label = paste(nm, "random-number state"))
  }
  # Different sample sizes error with "not all arguments have the same length"
  # after fold estimation, matching seminrExtras's RNG state.
  m_short <- suppressMessages(seminr::estimate_pls(seminr::mobi[1:200, ], mm, sm2))
  expect_error(suppressMessages(seminrExtras::assess_cvpat_compare(m1, m_short, nboot = 10L, seed = 19L, noFolds = 5L)),
               "not all arguments have the same length")
  sw_short <- .Random.seed
  expect_error(aw_assess_cvpat_compare(m1, m_short, nboot = 10L, seed = 19L, noFolds = 5L),
               "not all arguments have the same length")
  expect_identical(.Random.seed, sw_short)
  # Mismatched endogenous constructs or indicators error before seeding RNG.
  sm_diff <- seminr::relationships(
    seminr::paths(from = c("Image", "Expectation"), to = "Value"),
    seminr::paths(from = "Value", to = "Satisfaction")
  )
  m_diff <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm, sm_diff))
  set.seed(88L)
  s0 <- .Random.seed
  expect_error(aw_assess_cvpat_compare(m1, m_diff, seed = 123L),
               "CVPAT can only be applied to models with identical endogenous constructs and measures")
  expect_identical(.Random.seed, s0)
  mm_diff_items <- seminr::constructs(
    seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
    seminr::composite("Expectation", seminr::multi_items("CUEX", 1:2)),
    seminr::composite("Value", seminr::multi_items("PERV", 1:2)),
    seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3))
  )
  m_diff_items <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm_diff_items, sm1))
  expect_error(aw_assess_cvpat_compare(m1, m_diff_items, seed = 123L),
               "CVPAT can only be applied to models with identical endogenous constructs and measures")
  expect_identical(.Random.seed, s0)
})

test_that("CVPAT validates arguments, warns on unsupported objects, and delegates or falls back cleanly", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  fits <- aw_pls_fit_models()
  ecsi <- fits[["mobi_ecsi_core"]]
  # Non-SEMinR object warns and returns NULL; assess_cvpat seeds RNG before checking.
  set.seed(123L)
  s_seeded <- .Random.seed
  set.seed(1L)
  expect_warning(res <- aw_assess_cvpat(list(), seed = 123L), "assess_cvpat only works with SEMinR models")
  expect_null(res)
  expect_identical(.Random.seed, s_seeded)
  set.seed(1L)
  s_unseeded <- .Random.seed
  expect_warning(res_cmp <- aw_assess_cvpat_compare(list(), ecsi, seed = 123L), "assess_cvpat_compare only works with SEMinR models")
  expect_null(res_cmp)
  expect_identical(.Random.seed, s_unseeded)
  # Higher-order construct model warns and returns NULL.
  hoc_model <- ecsi
  hoc_model$hoc <- TRUE
  expect_warning(res_hoc <- aw_assess_cvpat(hoc_model, seed = 123L),
                 "There is no published solution for applying PLSpredict to higher-order models")
  expect_null(res_hoc)
  expect_warning(res_hoc_cmp <- aw_assess_cvpat_compare(hoc_model, ecsi, seed = 123L),
                 "There is no published solution for applying PLSpredict to higher-order models")
  expect_null(res_hoc_cmp)
  # Invalid nboot and reps arguments.
  for (bad_b in list(-1, "abc", c(10, 20), NA_real_)) {
    expect_error(aw_assess_cvpat(ecsi, nboot = bad_b), "invalid 'times' argument")
    expect_error(aw_assess_cvpat_compare(ecsi, ecsi, nboot = bad_b), "invalid 'times' argument")
  }
  for (bad_r in list(0, -1, 1.5, "abc", c(1, 2), NA_real_)) {
    expect_error(aw_assess_cvpat(ecsi, reps = bad_r), "reps must be NULL or a whole number of at least 1")
    expect_error(aw_assess_cvpat_compare(ecsi, ecsi, reps = bad_r), "reps must be NULL or a whole number of at least 1")
  }
  # nboot = 0 or 1 warns on near-zero bootstrap variance and sets T and P to NA.
  for (b01 in c(0L, 1L)) {
    w <- testthat::capture_warnings(res_b <- aw_assess_cvpat(ecsi, nboot = b01, noFolds = 5L))
    expect_gt(length(w), 0L)
    expect_match(w, "Bootstrap variance near zero; t-statistic set to NA")
    expect_true(all(is.na(res_b$CVPAT_compare_LM[, c("Boot T value", "Boot P Value")])))
  }
  # Standalone fallback loss helper matches native losses on composite models
  # (including models with missing_value = -99 mean replacement) and works on
  # interaction models without requiring seminrExtras.
  set.seed(19L)
  l_nat <- arwad:::aw_cvpat_native_losses(ecsi, seminr::predict_DA, 5L, 1L, 1L, need_lm = TRUE)
  set.seed(19L)
  l_sem <- arwad:::aw_cvpat_seminr_losses(ecsi, seminr::predict_DA, 5L, 1L, 1L, need_lm = TRUE)
  expect_identical(l_sem, l_nat)
  mobi_miss <- seminr::mobi
  mobi_miss[1L, "IMAG1"] <- -99
  mobi_miss[2L, "CUSA1"] <- -99
  m_miss <- suppressMessages(seminr::estimate_pls(mobi_miss, ecsi$measurement_model, ecsi$structural_model, missing_value = -99))
  miss_lt <- aw_assess_cvpat(m_miss, testtype = "less", nboot = 25L, seed = 29L, noFolds = 5L, reps = 1L)
  expect_false(anyNA(miss_lt$CVPAT_compare_IA))
  expect_false(anyNA(miss_lt$CVPAT_compare_LM))
  ints <- aw_pls_fit_interaction_models()
  if (requireNamespace("seminrExtras", quietly = TRUE)) {
    miss_gt <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat(m_miss, testtype = "greater", nboot = 25L, seed = 29L, noFolds = 5L, reps = 1L)))
    for (nm in c("CVPAT_compare_LM", "CVPAT_compare_IA")) {
      expect_identical(miss_lt[[nm]][, 1:3], miss_gt[[nm]][, 1:3])
      expect_identical(miss_lt[[nm]][, 4L], miss_gt[[nm]][, 4L])
      expect_equal(unname(miss_lt[[nm]][, 5L] + miss_gt[[nm]][, 5L]), rep(1, nrow(miss_lt[[nm]])), tolerance = 1e-6)
    }
    want_miss <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat(m_miss, nboot = 15L, seed = 12L, noFolds = 5L, reps = 1L)))
    got_miss <- aw_assess_cvpat(m_miss, nboot = 15L, seed = 12L, noFolds = 5L, reps = 1L)
    expect_same_cvpat(got_miss, want_miss, "missing_value = -99")
    want_int <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat(ints[[1L]], nboot = 15L, seed = 12L, noFolds = 5L, reps = 1L)))
    sw_int <- .Random.seed
    got_int <- aw_assess_cvpat(ints[[1L]], nboot = 15L, seed = 12L, noFolds = 5L, reps = 1L)
    expect_same_cvpat(got_int, want_int, "interaction model delegated")
    expect_identical(.Random.seed, sw_int, label = "interaction model delegated random-number state")
    # Reflective constructs (consistent PLS) are predicted by seminr's own algorithm: delegated.
    refl <- fits[["mobi_ecsi_reflective"]]
    want_refl <- suppressWarnings(suppressMessages(
      seminrExtras::assess_cvpat(refl, nboot = 15L, seed = 13L, noFolds = 5L)))
    sw_refl <- .Random.seed
    got_refl <- aw_assess_cvpat(refl, nboot = 15L, seed = 13L, noFolds = 5L)
    expect_same_cvpat(got_refl, want_refl, "reflective model delegated")
    expect_identical(.Random.seed, sw_refl, label = "reflective model delegated random-number state")
  }
})

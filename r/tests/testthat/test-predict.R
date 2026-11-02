# aw_predict_pls() against seminr::predict_pls(), bit for bit (ADR-0008).

# Every part of a predict_pls_model, signed zeros included, and the class.
expect_same_prediction <- function(got, want, label, native = TRUE) {
  expect_identical(got$native, native, label = paste(label, "native"))
  expect_identical(class(got), class(want), label = paste(label, "class"))
  expect_true(identical(unclass(got)[names(want)], unclass(want), num.eq = FALSE), label = paste(label, "identical to predict_pls()"))
}

# aw_predict_pls() is native from seminr 2.6.0 on and passes every call to
# predict_pls() under earlier versions, whose predict_pls() can loop forever
# with predict_EA on some structural models (fixed in 2.6.0); the tests need
# 2.6.0.
skip_if_old_seminr <- function() testthat::skip_if_not(arwad:::aw_seminr_260(), "seminr before 2.6.0")

# seminr::predict_pls() and aw_predict_pls() from the same seed; the random-
# number state after each call must agree too.
aw_predict_both <- function(model, seed, ...) {
  set.seed(seed)
  want <- suppressWarnings(suppressMessages(seminr::predict_pls(model, ...)))
  sw <- .Random.seed
  set.seed(seed)
  got <- aw_predict_pls(model, ...)
  list(got = got, want = want, same_rng = identical(.Random.seed, sw))
}

test_that("aw_predict_pls() equals seminr::predict_pls() on at least 50 models of composites", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  fits <- aw_pls_fit_models()
  comp <- Filter(function(m) !any(m$mmMatrix[, "type"] == "C"), fits)
  expect_gte(length(comp), 50L)
  for (k in seq_along(comp)) {
    nm <- names(comp)[k]
    tech <- if (k %% 2L == 1L) seminr::predict_EA else seminr::predict_DA
    r <- aw_predict_both(comp[[nm]], k, technique = tech, noFolds = c(10L, 5L, 7L)[k %% 3L + 1L])
    expect_same_prediction(r$got, r$want, nm)
    expect_true(r$same_rng, label = paste(nm, "random-number state"))
  }
  # The comparisons above are with native code, not with GNU R running the
  # kernel after a static fallback or a bailout.
  u <- arwad:::aw_pls_unit()
  expect_s3_class(u, "arwad_function")
  utils::capture.output(info <- aw_explain(u))
  expect_identical(sum(info$bailouts), 0L)
})

test_that("aw_predict_pls() repeats the cross-validation and runs leave-one-out as predict_pls() does", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  fits <- aw_pls_fit_models()
  ecsi <- fits[["mobi_ecsi_core"]]
  for (cfg in list(list(seminr::predict_EA, 10L, 3L), list(seminr::predict_DA, 5L, 2L), list(seminr::predict_EA, 10L, 1L))) {
    r <- aw_predict_both(ecsi, 42L, technique = cfg[[1L]], noFolds = cfg[[2L]], reps = cfg[[3L]])
    expect_same_prediction(r$got, r$want, paste("reps", cfg[[3L]]))
    expect_true(r$same_rng, label = paste("reps", cfg[[3L]], "random-number state"))
  }
  # Leave-one-out (noFolds = NULL) on the smallest model of composites.
  comp <- Filter(function(m) !any(m$mmMatrix[, "type"] == "C"), fits)
  small <- comp[[which.min(vapply(comp, function(m) nrow(m$data), 1L))]]
  r <- aw_predict_both(small, 7L, technique = seminr::predict_EA)
  expect_same_prediction(r$got, r$want, "leave-one-out")
  expect_true(r$same_rng)
  # A single predicted indicator: seminr's matrices drop to vectors, which
  # keep or lose names depending on reps.
  mm1 <- seminr::constructs(seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
                            seminr::composite("Satisfaction", seminr::single_item("CUSA1")))
  m1 <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm1, seminr::relationships(seminr::paths(from = "Image", to = "Satisfaction"))))
  for (R in list(NULL, 2L)) {
    r <- aw_predict_both(m1, 8L, technique = seminr::predict_DA, noFolds = 5L, reps = R)
    expect_same_prediction(r$got, r$want, paste("single indicator, reps", format(R)))
  }
  # Data with row names of their own.
  mob <- seminr::mobi
  rownames(mob) <- paste0("case", seq_len(nrow(mob)))
  m <- suppressMessages(seminr::estimate_pls(mob, ecsi$measurement_model, ecsi$structural_model))
  r <- aw_predict_both(m, 3L, technique = seminr::predict_DA, noFolds = 10L)
  expect_same_prediction(r$got, r$want, "row names")
  expect_identical(utils::capture.output(summary(r$got)), utils::capture.output(summary(r$want)))
})

test_that("aw_predict_pls() gives the same result on several cores", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  ecsi <- aw_pls_fit_models()[["mobi_ecsi_core"]]
  set.seed(5)
  one <- aw_predict_pls(ecsi, technique = seminr::predict_EA, noFolds = 10L, reps = 2L)
  set.seed(5)
  two <- aw_predict_pls(ecsi, technique = seminr::predict_EA, noFolds = 10L, reps = 2L, cores = 2L)
  expect_identical(two$native, one$native)
  expect_true(identical(two, one, num.eq = FALSE))
})

test_that("aw_predict_pls() hands consistent PLS and interaction terms to seminr", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  fits <- aw_pls_fit_models()
  plsc <- Filter(function(m) any(m$mmMatrix[, "type"] == "C"), fits)
  bad <- attr(fits, "inadmissible")
  # A consistent-PLS model whose full-sample solution is admissible; seminr
  # warns about inadmissible folds, or stops when every fold is.
  nm <- names(plsc)[!bad[names(plsc)]][1L]
  set.seed(1)
  want <- tryCatch(suppressWarnings(suppressMessages(seminr::predict_pls(plsc[[nm]], technique = seminr::predict_EA, noFolds = 10L))),
                   error = function(e) e)
  set.seed(1)
  got <- tryCatch(suppressWarnings(suppressMessages(aw_predict_pls(plsc[[nm]], technique = seminr::predict_EA, noFolds = 10L))),
                  error = function(e) e)
  if (inherits(want, "error")) {
    expect_identical(conditionMessage(got), conditionMessage(want))
  } else {
    expect_same_prediction(got, want, nm, native = FALSE)
  }
  ints <- aw_pls_fit_interaction_models()
  r <- aw_predict_both(ints[[1L]], 2L, technique = seminr::predict_DA, noFolds = 5L)
  expect_same_prediction(r$got, r$want, names(ints)[1L], native = FALSE)
  expect_true(r$same_rng)
})

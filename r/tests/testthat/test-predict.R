# aw_predict_pls() against seminr::predict_pls(), bit for bit (ADR-0008).

# Every part of a predict_pls_model, signed zeros included, and the class.
expect_same_prediction <- function(got, want, label, native = TRUE) {
  expect_identical(got$native, native, label = paste(label, "native"))
  expect_identical(class(got), class(want), label = paste(label, "class"))
  expect_true(identical(unclass(got)[names(want)], unclass(want), num.eq = FALSE), label = paste(label, "identical to predict_pls()"))
  if (native) {
    expect_type(got$q2_predict, "double")
    expect_identical(aw_q2_predict(got), got$q2_predict, label = paste(label, "aw_q2_predict(got)"))
  } else {
    expect_null(got$q2_predict, label = paste(label, "q2_predict null when delegated"))
  }
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
    if (cfg[[3L]] < 2L) next
    # The repetitions' averages, against single runs from the same stream; a
    # failure names the forms of averaging that would have matched.
    set.seed(42L)
    one <- lapply(seq_len(cfg[[3L]]), function(i) suppressWarnings(suppressMessages(
      seminr::predict_pls(ecsi, technique = cfg[[1L]], noFolds = cfg[[2L]]))))
    for (el in c("PLS_out_of_sample", "PLS_in_sample", "lm_out_of_sample", "lm_in_sample")) {
      A <- simplify2array(lapply(one, function(z) unname(z$items[[el]])))
      w <- unname(r$want$items[[el]])
      why <- sprintf("mean() %s, rowMeans() %s, sum() / reps %s, Reduce(`+`) / reps %s", identical(apply(A, c(1L, 2L), mean), w),
                     identical(rowMeans(A, dims = 2L), w), identical(apply(A, c(1L, 2L), sum) / dim(A)[3L], w),
                     identical(Reduce(`+`, lapply(seq_len(dim(A)[3L]), function(i) A[, , i])) / dim(A)[3L], w))
      expect_identical(unname(r$got$items[[el]]), w, label = paste("reps", cfg[[3L]], el, "average; seminr's equals", why))
    }
  }
  # An exogenous construct (Quality) that does not reach every endogenous
  # one (Expectation): the benchmark under earliest antecedents still
  # regresses each endogenous construct on the indicators of every
  # exogenous construct.
  mmr <- seminr::constructs(seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
                            seminr::composite("Expectation", seminr::multi_items("CUEX", 1:3)),
                            seminr::composite("Quality", seminr::multi_items("PERQ", 1:7)),
                            seminr::composite("Value", seminr::multi_items("PERV", 1:2)),
                            seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3)))
  smr <- seminr::relationships(seminr::paths(from = "Image", to = "Expectation"), seminr::paths(from = "Quality", to = "Value"),
                               seminr::paths(from = c("Expectation", "Value"), to = "Satisfaction"))
  mr <- suppressMessages(seminr::estimate_pls(seminr::mobi, mmr, smr))
  for (tech in list(seminr::predict_EA, seminr::predict_DA)) {
    r <- aw_predict_both(mr, 11L, technique = tech, noFolds = 10L)
    expect_same_prediction(r$got, r$want, "partial reach")
    expect_true(r$same_rng)
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

test_that("the repetitions' average is mean() of each element, bit for bit, compiled", {
  skip_if_no_toolchain()
  # Signed zeros, extremes, small and repeating values, infinities; NA and
  # NaN only in elements of their own (R leaves to the platform which of
  # the two a mix of them gives). No element's finite values sum past the
  # largest double: there GNU R on macOS arm64 returned a mean one unit in
  # the last place from the correctly rounded one, which Arwad's runtime
  # returns (ADR-0008); predictions never come near that range. A
  # failure lists each differing element: its values, then mean()'s result
  # and the unit's.
  odd <- c(-0, 0, Inf, -Inf, 1e300, 1e-300, 1 / 3, 2.5)
  set.seed(3)
  for (R in c(1L, 2L, 3L, 10L)) {
    x <- lapply(seq_len(R), function(r) matrix(c(sample(odd), stats::rnorm(32, 4, 2)), 8, 5))
    x[[1L]][1L, 1L] <- -0
    x[[1L]][3L, 3L] <- NA
    x[[R]][4L, 4L] <- NaN
    if (R > 1L) x[[R]][2L, 2L] <- 1e300
    A <- simplify2array(x)
    want <- apply(A, c(1L, 2L), mean)
    got <- arwad:::aw_pls_rep_mean(x)
    same <- vapply(seq_along(want), function(i) identical(got[i], want[i], num.eq = FALSE), NA)
    dim(A) <- c(length(want), R)
    why <- vapply(which(!same), function(i) sprintf("[%d] %s: mean() %a, compiled %a", i, paste(sprintf("%a", A[i, ]), collapse = " "),
                                                    want[i], got[i]), "")
    expect_true(all(same), label = paste(R, "repetitions;", paste(why, collapse = "; ")))
  }
  u <- arwad:::aw_pls_rep_mean_unit()
  expect_s3_class(u, "arwad_function")
  utils::capture.output(info <- aw_explain(u))
  expect_identical(sum(info$bailouts), 0L)
})

# Training sets of a few models: each model re-estimated by seminr on 90% of
# its cases, in shuffled order, with its own settings, and the held-out rest.
aw_predict_folds <- function() {
  fits <- aw_pls_fit_models()
  comp <- Filter(function(m) !any(m$mmMatrix[, "type"] == "C"), fits)
  pick <- unique(c(intersect(c("mobi_ecsi_core", "mobi_ecsi_full", "corp_rep_extended_modeB"), names(comp)), names(comp)[c(5L, 20L, 40L)]))
  lapply(stats::setNames(pick, pick), function(nm) {
    m <- comp[[nm]]
    n <- nrow(m$data)
    set.seed(9)
    perm <- sample(n)
    tr <- m$data[perm[seq_len(n - n %/% 10L)], , drop = FALSE]
    list(model = suppressMessages(seminr::estimate_pls(tr, m$measurement_model, m$structural_model, maxIt = m$settings$maxIt,
                                                        stopCriterion = m$settings$stopCriterion)),
         train = tr, test = m$data[perm[-seq_len(n - n %/% 10L)], , drop = FALSE])
  })
}

# These two tests pin down, platform by platform, the arithmetic that
# aw_predict_pls() takes from seminr; a failure names the alternative forms
# that would have matched.
test_that("fold models standardise with the means and standard deviations estimate_pls() stores", {
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  for (nm in names(folds <- aw_predict_folds())) {
    mt <- folds[[nm]]$model
    items <- as.vector(mt$mmVariables)
    X <- as.matrix(folds[[nm]]$train[, items])
    sc <- arwad:::aw_pls_train_scale(X)
    md <- unname(mt$meanData[items])
    sd <- unname(mt$sdData[items])
    mc <- colMeans(X)
    mm <- apply(X, 2L, mean)
    why <- sprintf("means as colMeans() %s, mean() %s; sds centred by colMeans() %s, by mean() %s, sd() %s",
                   identical(unname(mc), md), identical(unname(mm), md),
                   identical(unname(sqrt(colSums(sweep(X, 2L, mc)^2) / (nrow(X) - 1L))), sd),
                   identical(unname(sqrt(colSums(sweep(X, 2L, mm)^2) / (nrow(X) - 1L))), sd), identical(unname(apply(X, 2L, stats::sd)), sd))
    expect_identical(unname(sc$mu), md, label = paste(nm, "training means;", why))
    expect_identical(unname(sc$s), sd, label = paste(nm, "training sds;", why))
  }
})

test_that("a fold model predicts held-out cases as predict() does", {
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  for (nm in names(folds <- aw_predict_folds())) {
    mt <- folds[[nm]]$model
    spec <- arwad:::aw_pls_spec(mt)
    p <- length(spec$items)
    W <- mt$outer_weights[spec$items, spec$cons, drop = FALSE]
    o <- c(W[cbind(seq_len(p), spec$blk)], as.vector(mt$outer_loadings[spec$items, spec$cons]), as.vector(mt$path_coef[spec$cons, spec$cons]))
    ord <- arwad:::aw_pls_order(spec$pre)
    Xte <- as.matrix(folds[[nm]]$test[, spec$items])
    mu <- mt$meanData[spec$items]
    s <- mt$sdData[spec$items]
    for (ea in c(TRUE, FALSE)) {
      pr <- predict(mt, testData = folds[[nm]]$test, technique = if (ea) seminr::predict_EA else seminr::predict_DA)
      got <- arwad:::aw_pls_fold_predict(o, spec, Xte, mu, s, ea, ord)
      want <- unname(as.matrix(pr$predicted_composite_scores)[, spec$cons])
      # The composite scores as one product with the zero weights and one
      # construct at a time; the chain over the antecedents, over the whole
      # path column and by explicit products, without BLAS. A BLAS like
      # OpenBLAS tells these apart; R's reference BLAS does not.
      Z <- t((t(Xte) - mu) / s)
      S <- Z %*% W
      S1 <- sapply(seq_along(spec$cons), function(j) { k <- which(spec$blk == j); Z[, k, drop = FALSE] %*% W[k, j] })
      B <- mt$path_coef[spec$cons, spec$cons]
      chain <- function(f) {
        P <- S1
        for (j in ord) if (any(spec$pre[, j] > 0)) P[, j] <- f(if (ea) P else S1, j, which(spec$pre[, j] > 0))
        unname(P)
      }
      endo <- which(colSums(spec$pre > 0) > 0)
      why <- sprintf(paste("composite scores one construct at a time %s, as Z %%*%% W %s; endogenous scores over the antecedents %s,",
                           "over the path column %s, by explicit products %s"),
                     identical(unname(S1[, -endo]), want[, -endo]), identical(unname(S[, -endo]), want[, -endo]),
                     identical(chain(function(X, j, a) X[, a, drop = FALSE] %*% B[a, j]), want),
                     identical(chain(function(X, j, a) X %*% B[, j]), want),
                     identical(chain(function(X, j, a) { v <- 0; for (i in a) v <- v + X[, i] * B[i, j]; v }), want))
      lab <- paste(nm, if (ea) "EA" else "DA")
      expect_identical(unname(got$scores), want, label = paste(lab, "scores;", why))
      expect_identical(unname(got$items), unname(as.matrix(pr$predicted_items)[, spec$items]), label = paste(lab, "items;", why))
    }
  }
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

# Reference Q2_predict from the fold partition drawn by sample(n): each held-
# out case is benchmarked against the indicator means of its training fold,
# averaged across repetitions as PLS_out_of_sample is.
aw_ref_q2_predict <- function(model, pls_oos, seed, noFolds = NULL, reps = NULL) {
  n <- nrow(model$data)
  K <- if (is.null(noFolds)) n else as.integer(noFolds)
  R <- if (is.null(reps)) 1L else as.integer(reps)
  endo <- unique(model$smMatrix[, "target"])
  mm <- model$mmMatrix[model$mmMatrix[, "construct"] %in% endo, , drop = FALSE]
  yitems <- unname(unlist(lapply(endo, function(co) mm[mm[, "construct"] == co, "measurement"])))
  Y <- as.matrix(model$data[, yitems, drop = FALSE])
  set.seed(seed)
  tps <- split(seq_len(n), cut(seq(1, n), breaks = K, labels = FALSE))
  perms <- if (arwad:::aw_seminr_260()) lapply(seq_len(R), function(r) sample(n)) else rep(list(sample(n)), R)
  ia_reps <- lapply(seq_len(R), function(r) {
    perm <- perms[[r]]
    ia <- matrix(0, n, length(yitems))
    for (f in seq_len(K)) {
      te <- perm[tps[[f]]]
      tr <- perm[-tps[[f]]]
      mu <- colMeans(Y[tr, , drop = FALSE])
      ia[te, ] <- rep(mu, each = length(te))
    }
    ia
  })
  ia_oos <- apply(simplify2array(ia_reps), c(1L, 2L), mean)
  P <- as.matrix(pls_oos)
  stats::setNames(1 - colSums((Y - P)^2) / colSums((Y - ia_oos)^2), yitems)
}

test_that("aw_predict_pls() and aw_q2_predict() compute Q2_predict from training-fold means", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  fits <- aw_pls_fit_models()
  ecsi <- fits[["mobi_ecsi_core"]]
  for (cfg in list(list(seminr::predict_DA, 10L, NULL), list(seminr::predict_EA, 10L, 1L), list(seminr::predict_EA, 5L, 3L))) {
    set.seed(77L)
    pr <- aw_predict_pls(ecsi, technique = cfg[[1L]], noFolds = cfg[[2L]], reps = cfg[[3L]])
    sw <- .Random.seed
    want <- aw_ref_q2_predict(ecsi, pr$items$PLS_out_of_sample, 77L, noFolds = cfg[[2L]], reps = cfg[[3L]])
    expect_identical(pr$q2_predict, want, label = paste("q2_predict reps", format(cfg[[3L]])))
    # Extracting from a native predict_pls_model does not touch the RNG state.
    set.seed(99L)
    s0 <- .Random.seed
    expect_identical(aw_q2_predict(pr), want)
    expect_identical(.Random.seed, s0)
    # Computing directly from a pls_model matches and advances the RNG state identically.
    set.seed(77L)
    q2 <- aw_q2_predict(ecsi, technique = cfg[[1L]], noFolds = cfg[[2L]], reps = cfg[[3L]])
    expect_identical(q2, want)
    expect_identical(.Random.seed, sw)
  }
  # Leave-one-out (noFolds = NULL).
  comp <- Filter(function(m) !any(m$mmMatrix[, "type"] == "C"), fits)
  small <- comp[[which.min(vapply(comp, function(m) nrow(m$data), 1L))]]
  set.seed(19L)
  pr_loo <- aw_predict_pls(small, technique = seminr::predict_DA, noFolds = NULL)
  expect_identical(pr_loo$q2_predict, aw_ref_q2_predict(small, pr_loo$items$PLS_out_of_sample, 19L, noFolds = NULL))
  # Single endogenous indicator (where items$PLS_out_of_sample drops to a vector).
  mm1 <- seminr::constructs(seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
                            seminr::composite("Satisfaction", seminr::single_item("CUSA1")))
  m1 <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm1, seminr::relationships(seminr::paths(from = "Image", to = "Satisfaction"))))
  for (R in list(NULL, 2L)) {
    set.seed(23L)
    pr1 <- aw_predict_pls(m1, technique = seminr::predict_DA, noFolds = 5L, reps = R)
    expect_identical(pr1$q2_predict, aw_ref_q2_predict(m1, pr1$items$PLS_out_of_sample, 23L, noFolds = 5L, reps = R))
    expect_identical(names(pr1$q2_predict), "CUSA1")
  }
  # A foreign predict_pls_model (no training-fold means stored) errors clearly;
  # passing a delegated pls_model directly to aw_q2_predict() computes Q2_predict via seminr folds.
  set.seed(31L)
  foreign <- suppressWarnings(suppressMessages(seminr::predict_pls(ecsi, noFolds = 5L)))
  expect_error(aw_q2_predict(foreign), "training-fold indicator means")
  expect_error(aw_q2_predict(list()), "pls_model")
  ints <- aw_pls_fit_interaction_models()
  set.seed(41L)
  q2_int <- aw_q2_predict(ints[[1L]], technique = seminr::predict_DA, noFolds = 5L, reps = 2L)
  set.seed(41L)
  pr_int <- aw_predict_pls(ints[[1L]], technique = seminr::predict_DA, noFolds = 5L, reps = 2L)
  expect_identical(q2_int, aw_ref_q2_predict(ints[[1L]], pr_int$items$PLS_out_of_sample, 41L, noFolds = 5L, reps = 2L))
  # Higher-order construct models emit seminr's message and return NULL in aw_q2_predict().
  mm_hoc <- seminr::constructs(
    seminr::composite("Image", seminr::multi_items("IMAG", 1:3)),
    seminr::composite("Value", seminr::multi_items("PERV", 1:2)),
    seminr::higher_composite("Higher", c("Image", "Value"), method = seminr::two_stage),
    seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3))
  )
  sm_hoc <- seminr::relationships(seminr::paths(from = "Higher", to = "Satisfaction"))
  m_hoc <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm_hoc, sm_hoc))
  expect_message(q2_hoc <- aw_q2_predict(m_hoc, noFolds = 5L), "higher-order")
  expect_null(q2_hoc)
})

test_that("aw_q2_predict() matches seminr's fold draws across versions, including seminr 2.5.0 with reps > 1", {
  skip_if_no_toolchain()
  skip_if_not_installed("seminr")
  mm <- seminr::constructs(
    seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
    seminr::composite("Value", seminr::multi_items("PERV", 1:2)),
    seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3))
  )
  sm <- seminr::relationships(
    seminr::paths(from = "Image", to = c("Value", "Satisfaction")),
    seminr::paths(from = "Value", to = "Satisfaction")
  )
  m <- suppressMessages(seminr::estimate_pls(seminr::mobi, mm, sm))
  for (R in list(NULL, 1L, 2L)) {
    set.seed(53L)
    pr <- aw_predict_pls(m, technique = seminr::predict_DA, noFolds = 5L, reps = R)
    sw <- .Random.seed
    set.seed(53L)
    q2 <- aw_q2_predict(m, technique = seminr::predict_DA, noFolds = 5L, reps = R)
    expect_identical(.Random.seed, sw)
    expect_identical(q2, aw_ref_q2_predict(m, pr$items$PLS_out_of_sample, 53L, noFolds = 5L, reps = R))
  }
})

test_that("aw_q2_predict() redraws a delegated model's folds from the same state in a fresh session", {
  skip_if_not_installed("seminr")
  skip_if_old_seminr()
  # Reflective constructs go to seminr::predict_pls(); aw_q2_predict() then
  # draws the folds again for the training-fold means. Without a
  # .Random.seed, both draws must come from one initialisation, not two.
  m <- aw_pls_fit_models()[["mobi_ecsi_reflective"]]
  n <- nrow(m$data)
  keep <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (!is.null(keep)) assign(".Random.seed", keep, envir = globalenv()), add = TRUE)
  if (!is.null(keep)) rm(".Random.seed", envir = globalenv())
  perms <- new.env()
  perms$p <- list()
  record <- function(v) perms$p[[length(perms$p) + 1L]] <- v
  invisible(utils::capture.output(suppressMessages(
    trace("sample", print = FALSE, where = asNamespace("arwad"),
          exit = bquote(if (length(x) == 1L && x == .(n)) .(record)(returnValue()))))))
  on.exit(suppressMessages(untrace("sample", where = asNamespace("arwad"))), add = TRUE)
  # The folds are random here, and seminr 2.6.0 gives NA rows for a fold
  # whose consistent-PLS solution is inadmissible, so Q2 itself may be NA;
  # the test is about the folds.
  q2 <- suppressWarnings(aw_q2_predict(m, noFolds = 5L))
  expect_length(perms$p, 2L)
  expect_identical(perms$p[[2L]], perms$p[[1L]], label = "redrawn folds")
  expect_true(is.double(q2) && !is.null(names(q2)))
})

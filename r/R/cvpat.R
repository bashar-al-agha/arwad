# Cross-validated predictive ability test (CVPAT; Liengaard et al., 2021;
# Sharma et al., 2023), as seminrExtras::assess_cvpat() and
# seminrExtras::assess_cvpat_compare() compute it. Each rule below was matched
# to seminrExtras 1.0.4 and seminr 2.6.0 by black-box tests on their outputs
# and on the values they pass to base R's mean(), var(), pt() and sample()
# (help pages, results and those values only; docs/adr/0008-pls-predict.md):
# - when seed is not NULL, set.seed(seed) is called once at entry; for reps > 1,
#   predict_pls() is run reps times in sequence on the continuous RNG stream
#   (in assess_cvpat_compare(), all reps of established_model run first, then
#   all reps of alternative_model);
# - per-observation loss for endogenous construct c is the unweighted mean of
#   its squared out-of-sample indicator residuals, and "Overall" is the
#   unweighted mean across the J endogenous constructs; with reps > 1 the
#   construct losses are summed over the repetitions in order and divided by
#   their number, and "Overall" is taken after that;
# - the indicator-average (IA) benchmark in assess_cvpat() centers each
#   predicted indicator by its full-sample mean on seminr_model$data;
# - after the predictions, bootstrap resampling runs across 2*(J + 1) slots
#   [IA Overall, LM Overall, IA endo_1..J, LM endo_1..J] in assess_cvpat() and
#   (J + 1) slots [Overall, endo_1..J] in assess_cvpat_compare(); each slot
#   draws 2*nboot resamples of size n with replacement and uses every second
#   resample to compute the bootstrap standard error of the mean loss
#   difference D = L_bench - L_pls (or L_alt - L_base): the mean of the centred
#   D - mean(D) on that resample, whose sd() is the standard error (centring
#   changes the last bits of the standard error, so it is part of the rule).

# The compiled bootstrap kernel for one CVPAT slot: draws two resamples of size
# n per replicate from R's RNG stream and averages D, passed centred, on the
# second resample with GNU R's extended-precision mean().
aw_cvpat_boot_kernel <- function(D, B) {
  n <- length(D)
  out <- numeric(B)
  for (b in seq_len(B)) {
    s1 <- sample.int(n, n, replace = TRUE)
    s2 <- sample.int(n, n, replace = TRUE)
    out[b] <- mean(D[s2])
  }
  out
}

aw_cvpat_boot_unit <- function() {
  if (is.null(aw_state$cvpat_boot)) aw_state$cvpat_boot <- aw_compile(aw_cvpat_boot_kernel, name = "aw_cvpat_boot")
  aw_state$cvpat_boot
}

# Per-observation construct losses (n by J) from an n by q residual matrix E;
# col_idx is a list of column positions in E for each endogenous construct in
# endo.
aw_cvpat_loss_mat <- function(E, endo, col_idx) {
  n <- nrow(E)
  J <- length(endo)
  L <- matrix(0, n, J, dimnames = list(NULL, endo))
  for (j in seq_len(J)) L[, j] <- rowMeans(E[, col_idx[[j]], drop = FALSE]^2)
  L
}

# The construct losses with the "Overall" column (n by J + 1), taken after any
# averaging over repetitions.
aw_cvpat_overall <- function(L) cbind(L, Overall = rowMeans(L))

# The mean loss per row of a table: colMeans() for the constructs and mean()
# for "Overall", as seminrExtras computes them. The two differ in the last
# bit where R has no extended precision (macOS arm64), which the 15-digit
# columns of assess_cvpat_compare() show.
aw_cvpat_col_means <- function(L) {
  J <- ncol(L) - 1L
  c(colMeans(L[, seq_len(J), drop = FALSE]), Overall = mean(L[, J + 1L]))
}

# Native cross-validated per-observation loss matrices (PLS and optionally LM)
# averaged over R sequential single-rep runs on the caller's RNG stream, using
# one batched aw_pls_run_cores() pass per fold size. Returns NULL when the
# model or arguments are outside the native PLSpredict subset (without
# advancing .Random.seed).
aw_cvpat_native_losses <- function(model, technique, noFolds, R, cores, need_lm = TRUE) {
  cv <- aw_pls_cv_fits(model, technique, noFolds, R, cores, need_lm = need_lm)
  if (is.null(cv)) return(NULL)
  spec <- cv$spec
  ea <- cv$ea
  ord <- cv$ord
  n <- cv$n
  k <- cv$k
  endo <- cv$endo
  yi <- cv$yi
  q <- cv$q
  sets <- cv$sets
  perms <- cv$perms
  tps <- cv$tps
  fits <- cv$fits
  X <- spec$X
  col_idx <- lapply(cv$itc, function(cols) match(cols, yi))
  Xyi <- X[, yi, drop = FALSE]
  LP <- matrix(0, n, length(endo))
  LL <- if (need_lm) LP
  for (r in seq_len(R)) {
    oos <- matrix(NA_real_, n, q)
    oos_lm <- if (need_lm) oos
    for (f in seq_len(k)) {
      te <- perms[[r]][tps[[f]]]
      tr <- perms[[r]][-tps[[f]]]
      Xtr <- X[tr, , drop = FALSE]
      Xte <- X[te, , drop = FALSE]
      sc <- aw_pls_train_scale(Xtr)
      o <- fits[[(r - 1L) * k + f]]
      pout <- aw_pls_fold_predict(o, spec, Xte, sc$mu, sc$s, ea, ord)
      oos[te, ] <- pout$items[, yi]
      if (need_lm) {
        for (st in sets) {
          A <- cbind(1, Xtr[, st$x, drop = FALSE])
          fit <- .lm.fit(A, Xtr[, st$y, drop = FALSE])
          if (fit$rank < ncol(A)) return(cv$restore())
          Ate <- cbind(1, Xte[, st$x, drop = FALSE])
          b <- as.matrix(fit$coefficients)
          for (g in st$by) {
            cols <- match(st$y[g], yi)
            oos_lm[te, cols] <- Ate %*% b[, g, drop = FALSE]
          }
        }
      }
    }
    LP <- LP + aw_cvpat_loss_mat(Xyi - (oos + 0), endo, col_idx)
    if (need_lm) LL <- LL + aw_cvpat_loss_mat(Xyi - (oos_lm + 0), endo, col_idx)
  }
  LP <- aw_cvpat_overall(LP / R)
  if (need_lm) LL <- aw_cvpat_overall(LL / R)
  list(LP = LP, LL = LL, endo = endo, items = spec$items[yi], col_idx = col_idx, Xyi = Xyi)
}

# Fallback loss calculation through seminr::predict_pls() when seminrExtras is
# not installed or when testtype = "less" (which crashes in seminrExtras 1.0.4).
aw_cvpat_seminr_losses <- function(model, technique, noFolds, R, cores, need_lm = TRUE) {
  if (!inherits(model, "pls_model") && inherits(model, "cfa_model")) {
    stop("Prediction is not available for CFA models estimated with estimate_cfa(): a CFA has no structural model to predict from.",
         call. = FALSE)
  }
  if (is.null(model$data) || is.null(model$smMatrix)) {
    stop("'data' must be of a vector type, was 'NULL'", call. = FALSE)
  }
  endo <- unique(model$smMatrix[, "target"])
  itc <- lapply(endo, function(c) model$mmMatrix[model$mmMatrix[, "construct"] == c, "measurement"])
  items <- unlist(itc, use.names = FALSE)
  col_idx <- lapply(itc, function(cols) match(cols, items))
  Xyi <- unname(as.matrix(model$data[, items, drop = FALSE]))
  storage.mode(Xyi) <- "double"
  n <- nrow(Xyi)
  LP <- matrix(0, n, length(endo))
  LL <- if (need_lm) LP
  for (r in seq_len(R)) {
    p <- seminr::predict_pls(model, technique = technique, noFolds = noFolds, cores = cores)
    ep <- as.matrix(p$items$PLS_out_of_sample_residuals)
    if (anyNA(ep)) stop("'mu' must be a single number", call. = FALSE)
    LP <- LP + aw_cvpat_loss_mat(ep, endo, col_idx)
    if (need_lm) {
      el <- as.matrix(p$items$lm_out_of_sample_residuals)
      if (anyNA(el)) stop("'mu' must be a single number", call. = FALSE)
      LL <- LL + aw_cvpat_loss_mat(el, endo, col_idx)
    }
  }
  LP <- aw_cvpat_overall(LP / R)
  if (need_lm) LL <- aw_cvpat_overall(LL / R)
  list(LP = LP, LL = LL, endo = endo, items = items, col_idx = col_idx, Xyi = Xyi)
}

# Evaluates one CVPAT bootstrap slot on D = L_other - L_model (length n).
aw_cvpat_slot <- function(D, nboot, testtype) {
  n <- length(D)
  bm <- if (nboot >= 1L) aw_cvpat_boot_unit()(as.double(D - mean(D)), nboot) else numeric(0)
  s_bm <- if (nboot >= 2L) stats::sd(bm) else NA_real_
  if (is.na(s_bm) || s_bm < .Machine$double.eps) {
    warning("Bootstrap variance near zero; t-statistic set to NA", call. = FALSE)
    return(c(NA_real_, NA_real_))
  }
  tv <- mean(D) / s_bm
  pv <- switch(testtype,
    two.sided = 2 * stats::pt(-abs(tv), df = n - 1L),
    greater   = stats::pt(tv, df = n - 1L, lower.tail = FALSE),
    less      = stats::pt(tv, df = n - 1L, lower.tail = TRUE)
  )
  c(tv, pv)
}

# Validates nboot and reps as seminrExtras 1.0.4 does.
aw_cvpat_check_args <- function(nboot, reps) {
  if (!is.null(reps)) {
    if (!is.numeric(reps) || length(reps) != 1L || !is.finite(reps) || reps < 1 || reps != round(reps)) {
      stop("reps must be NULL or a whole number of at least 1.", call. = FALSE)
    }
  }
  if (is.null(nboot) || !is.numeric(nboot) || length(nboot) != 1L || is.na(nboot) || !is.finite(nboot) || nboot < 0) {
    stop("invalid 'times' argument", call. = FALSE)
  }
}

#' Cross-validated predictive ability test (CVPAT) for a PLS-SEM model, compiled
#'
#' A compiled replacement for `seminrExtras::assess_cvpat()` (Liengaard et al.,
#' 2021; Sharma et al., 2023). It compares the out-of-sample predictive loss of
#' a PLS-SEM model estimated by `seminr::estimate_pls()` against two benchmarks:
#' the linear-model (LM) benchmark and the indicator-average (IA) benchmark.
#' Cross-validated predictions are estimated with the compiled PLS-SEM engine
#' of [aw_predict_pls()] and the bootstrap t-statistics are computed with a
#' compiled resampling kernel, drawing from the same random-number stream and
#' returning the exact `table_output` matrices of `seminrExtras::assess_cvpat()`.
#'
#' All three `testtype` hypotheses (`"two.sided"`, `"greater"`, and `"less"`)
#' are supported (`seminrExtras` 1.0.4 errors on `testtype = "less"`). Models
#' outside the native PLSpredict subset (consistent PLS, interaction terms,
#' rank-deficient LM folds, or `seminr` before 2.6.0) are handed to
#' `seminrExtras::assess_cvpat()` when `seminrExtras` is installed and
#' `testtype != "less"`, or evaluated through `seminr::predict_pls()` otherwise.
#'
#' @param seminr_model A model returned by `seminr::estimate_pls()`.
#' @param testtype `"two.sided"` (the default), `"greater"` (one-sided test
#'   that the benchmark loss exceeds the PLS loss), or `"less"`.
#' @param nboot The number of bootstrap resamples for the t-statistic (default
#'   `2000`).
#' @param seed Random seed called once via `set.seed(seed)` before prediction
#'   (default `123`; `NULL` to use the current `.Random.seed`).
#' @param technique `seminr::predict_DA` (direct antecedents, the default) or
#'   `seminr::predict_EA` (earliest antecedents).
#' @param noFolds The number of cross-validation folds; `NULL` for
#'   leave-one-out.
#' @param reps The number of cross-validation repetitions; `NULL` for one.
#' @param cores The number of worker processes for fold estimation; `NULL` for
#'   one.
#' @return A list with two `(J + 1) x 5` matrices of class `table_output`:
#'   `CVPAT_compare_LM` and `CVPAT_compare_IA`, with rows for each endogenous
#'   construct and `"Overall"`, and columns `"PLS Loss"`, `"LM Loss"` (or
#'   `"IA Loss"`), `"Diff"`, `"Boot T value"`, and `"Boot P Value"`. Returns
#'   `NULL` with a warning when `seminr_model` is not a SEMinR model or is a
#'   higher-order model.
#' @references
#' Liengaard, B. D., Sharma, P. N., Hult, G. T. M., Jensen, M. B., Sarstedt,
#' M., Hair, J. F., & Ringle, C. M. (2021). Prediction: Coveted, yet forsaken?
#' Introducing a cross-validated predictive ability test in partial least
#' squares path modeling. *Decision Sciences, 52*(2), 362-392.
#'
#' Sharma, P. N., Liengaard, B. D., Hair, J. F., Sarstedt, M., & Ringle, C. M.
#' (2023). Predictive model assessment and selection in composite-based
#' modeling using PLS-SEM: Extensions and guidelines for using CVPAT. *European
#' Journal of Marketing, 57*(6), 1662-1677.
#' @export
#' @examplesIf requireNamespace("seminr", quietly = TRUE)
#' library(seminr)
#' mm <- constructs(composite("Image", multi_items("IMAG", 1:5)),
#'                  composite("Value", multi_items("PERV", 1:2)),
#'                  composite("Satisfaction", multi_items("CUSA", 1:3)))
#' sm <- relationships(paths(from = "Image", to = c("Value", "Satisfaction")),
#'                     paths(from = "Value", to = "Satisfaction"))
#' model <- estimate_pls(mobi, mm, sm)
#' aw_assess_cvpat(model, noFolds = 10, nboot = 200, seed = 123)
aw_assess_cvpat <- function(seminr_model, testtype = "two.sided", nboot = 2000, seed = 123,
                            technique = seminr::predict_DA, noFolds = NULL, reps = NULL, cores = NULL) {
  if (!requireNamespace("seminr", quietly = TRUE)) stop("aw_assess_cvpat() needs the seminr package")
  testtype <- match.arg(testtype, c("two.sided", "less", "greater"))
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  if (!is.null(seed)) set.seed(seed)
  if (!inherits(seminr_model, "seminr_model")) {
    warning("assess_cvpat only works with SEMinR models.", call. = FALSE)
    return(NULL)
  }
  if (isTRUE(seminr_model$hoc)) {
    warning("There is no published solution for applying PLSpredict to higher-order models.", call. = FALSE)
    return(NULL)
  }
  aw_cvpat_check_args(nboot, reps)
  R <- if (is.null(reps)) 1L else as.integer(reps)
  B <- as.integer(nboot)
  losses <- aw_cvpat_native_losses(seminr_model, technique, noFolds, R, cores, need_lm = TRUE)
  if (is.null(losses)) {
    if (testtype != "less" && requireNamespace("seminrExtras", quietly = TRUE)) {
      if (had) assign(".Random.seed", old, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
      return(seminrExtras::assess_cvpat(seminr_model = seminr_model, testtype = testtype, nboot = nboot, seed = seed,
                                        technique = technique, noFolds = noFolds, reps = reps, cores = cores))
    }
    losses <- aw_cvpat_seminr_losses(seminr_model, technique, noFolds, R, cores, need_lm = TRUE)
  }
  LP <- losses$LP
  LL <- losses$LL
  endo <- losses$endo
  J <- length(endo)
  Xyi <- losses$Xyi
  n <- nrow(Xyi)
  EI <- Xyi - matrix(colMeans(Xyi), n, ncol(Xyi), byrow = TRUE)
  LI <- aw_cvpat_overall(aw_cvpat_loss_mat(EI, endo, losses$col_idx))
  # Slot order in seminrExtras 1.0.4:
  # 1: IA Overall, 2: LM Overall, 3..(2+J): IA endo_1..J, (3+J)..(2+2J): LM endo_1..J
  tp_ia <- matrix(NA_real_, J + 1L, 2L)
  tp_lm <- matrix(NA_real_, J + 1L, 2L)
  tp_ia[J + 1L, ] <- aw_cvpat_slot(LI[, J + 1L] - LP[, J + 1L], B, testtype)
  tp_lm[J + 1L, ] <- aw_cvpat_slot(LL[, J + 1L] - LP[, J + 1L], B, testtype)
  for (j in seq_len(J)) tp_ia[j, ] <- aw_cvpat_slot(LI[, j] - LP[, j], B, testtype)
  for (j in seq_len(J)) tp_lm[j, ] <- aw_cvpat_slot(LL[, j] - LP[, j], B, testtype)
  mk_tab <- function(LB, tp, bench_col) {
    c1 <- aw_cvpat_col_means(LP)
    c2 <- aw_cvpat_col_means(LB)
    c3 <- c1 - c2
    m <- cbind(
      as.numeric(format(c1)),
      as.numeric(format(c2)),
      as.numeric(format(c3)),
      as.numeric(as.character(tp[, 1L])),
      as.numeric(as.character(tp[, 2L]))
    )
    dimnames(m) <- list(c(endo, "Overall"), c("PLS Loss", bench_col, "Diff", "Boot T value", "Boot P Value"))
    attr(m, "comment") <- "CVPAT as per Sharma et al. (2023)."
    class(m) <- c("table_output", "matrix", "array")
    m
  }
  list(
    CVPAT_compare_LM = mk_tab(LL, tp_lm, "LM Loss"),
    CVPAT_compare_IA = mk_tab(LI, tp_ia, "IA Loss")
  )
}

#' Compare two PLS-SEM models with CVPAT, compiled
#'
#' A compiled replacement for `seminrExtras::assess_cvpat_compare()` (Liengaard
#' et al., 2021; Sharma et al., 2023). It compares the cross-validated
#' out-of-sample predictive loss of an `established_model` and an
#' `alternative_model` that share identical endogenous constructs and
#' indicators, estimating every training fold with the compiled PLS-SEM engine
#' of [aw_predict_pls()] and computing the bootstrap t-statistics with a
#' compiled resampling kernel from the same random-number stream.
#'
#' All three `testtype` hypotheses (`"two.sided"`, `"greater"`, and `"less"`)
#' are supported (`seminrExtras` 1.0.4 errors on `testtype = "less"`).
#'
#' @param established_model Baseline model returned by `seminr::estimate_pls()`.
#' @param alternative_model Competing model returned by `seminr::estimate_pls()`,
#'   with the same endogenous constructs and indicators as `established_model`.
#' @param testtype `"two.sided"` (the default), `"greater"` (one-sided test
#'   that the alternative model loss exceeds the base model loss), or `"less"`.
#' @param nboot The number of bootstrap resamples for the t-statistic (default
#'   `2000`).
#' @param seed Random seed called once via `set.seed(seed)` before prediction
#'   (default `123`; `NULL` to use the current `.Random.seed`).
#' @param technique `seminr::predict_DA` (direct antecedents, the default) or
#'   `seminr::predict_EA` (earliest antecedents).
#' @param noFolds The number of cross-validation folds; `NULL` for
#'   leave-one-out.
#' @param reps The number of cross-validation repetitions; `NULL` for one.
#' @param cores The number of worker processes for fold estimation; `NULL` for
#'   one.
#' @return A `(J + 1) x 5` matrix of class `table_output` with rows for each
#'   endogenous construct and `"Overall"`, and columns `"Base Model Loss"`,
#'   `"Alt Model Loss"`, `"Diff"`, `"Boot T value"`, and `"Boot P Value"`.
#'   Returns `NULL` with a warning when either model is not a SEMinR model or
#'   is a higher-order model.
#' @export
#' @examplesIf requireNamespace("seminr", quietly = TRUE)
#' library(seminr)
#' mm1 <- constructs(composite("Image", multi_items("IMAG", 1:5)),
#'                   composite("Value", multi_items("PERV", 1:2)),
#'                   composite("Satisfaction", multi_items("CUSA", 1:3)))
#' sm1 <- relationships(paths(from = "Image", to = c("Value", "Satisfaction")),
#'                      paths(from = "Value", to = "Satisfaction"))
#' mm2 <- constructs(composite("Image", multi_items("IMAG", 1:5)),
#'                   composite("Expectation", multi_items("CUEX", 1:3)),
#'                   composite("Value", multi_items("PERV", 1:2)),
#'                   composite("Satisfaction", multi_items("CUSA", 1:3)))
#' sm2 <- relationships(paths(from = c("Image", "Expectation"), to = c("Value", "Satisfaction")),
#'                      paths(from = "Value", to = "Satisfaction"))
#' m1 <- estimate_pls(mobi, mm1, sm1)
#' m2 <- estimate_pls(mobi, mm2, sm2)
#' aw_assess_cvpat_compare(m1, m2, noFolds = 10, nboot = 200, seed = 123)
aw_assess_cvpat_compare <- function(established_model, alternative_model, testtype = "two.sided",
                                    nboot = 2000, seed = 123, technique = seminr::predict_DA,
                                    noFolds = NULL, reps = NULL, cores = NULL) {
  if (!requireNamespace("seminr", quietly = TRUE)) stop("aw_assess_cvpat_compare() needs the seminr package")
  if (!inherits(established_model, "seminr_model") || !inherits(alternative_model, "seminr_model")) {
    warning("assess_cvpat_compare only works with SEMinR models.", call. = FALSE)
    return(NULL)
  }
  if (isTRUE(established_model$hoc) || isTRUE(alternative_model$hoc)) {
    warning("There is no published solution for applying PLSpredict to higher-order models.", call. = FALSE)
    return(NULL)
  }
  testtype <- match.arg(testtype, c("two.sided", "less", "greater"))
  endo1 <- unique(established_model$smMatrix[, "target"])
  endo2 <- unique(alternative_model$smMatrix[, "target"])
  mm1 <- established_model$mmMatrix[established_model$mmMatrix[, "construct"] %in% endo1, , drop = FALSE]
  mm2 <- alternative_model$mmMatrix[alternative_model$mmMatrix[, "construct"] %in% endo2, , drop = FALSE]
  if (!setequal(endo1, endo2) || !setequal(mm1[, "measurement"], mm2[, "measurement"])) {
    stop("CVPAT can only be applied to models with identical endogenous constructs and measures", call. = FALSE)
  }
  aw_cvpat_check_args(nboot, reps)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  if (!is.null(seed)) set.seed(seed)
  R <- if (is.null(reps)) 1L else as.integer(reps)
  B <- as.integer(nboot)
  l1 <- aw_cvpat_native_losses(established_model, technique, noFolds, R, cores, need_lm = FALSE)
  l2 <- if (!is.null(l1)) aw_cvpat_native_losses(alternative_model, technique, noFolds, R, cores, need_lm = FALSE)
  if (is.null(l1) || is.null(l2)) {
    if (testtype != "less" && requireNamespace("seminrExtras", quietly = TRUE)) {
      if (had) assign(".Random.seed", old, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
      return(seminrExtras::assess_cvpat_compare(established_model = established_model,
                                                alternative_model = alternative_model,
                                                testtype = testtype, nboot = nboot, seed = seed,
                                                technique = technique, noFolds = noFolds,
                                                reps = reps, cores = cores))
    }
    if (!is.null(seed)) set.seed(seed)
    else if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    l1 <- aw_cvpat_seminr_losses(established_model, technique, noFolds, R, cores, need_lm = FALSE)
    l2 <- aw_cvpat_seminr_losses(alternative_model, technique, noFolds, R, cores, need_lm = FALSE)
  }
  endo <- l1$endo
  J <- length(endo)
  L_base <- l1$LP
  L_alt <- l2$LP[, c(match(endo, l2$endo), J + 1L), drop = FALSE]
  if (nrow(L_base) != nrow(L_alt)) {
    stop("not all arguments have the same length", call. = FALSE)
  }
  # Slot order in seminrExtras 1.0.4: 1: Overall, 2..(1+J): endo_1..J
  tp <- matrix(NA_real_, J + 1L, 2L)
  tp[J + 1L, ] <- aw_cvpat_slot(L_alt[, J + 1L] - L_base[, J + 1L], B, testtype)
  for (j in seq_len(J)) tp[j, ] <- aw_cvpat_slot(L_alt[, j] - L_base[, j], B, testtype)
  c1 <- aw_cvpat_col_means(L_base)
  c2 <- aw_cvpat_col_means(L_alt)
  c3 <- c1 - c2
  m <- cbind(
    as.numeric(as.character(c1)),
    as.numeric(as.character(c2)),
    as.numeric(as.character(c3)),
    as.numeric(as.character(tp[, 1L])),
    as.numeric(as.character(tp[, 2L]))
  )
  dimnames(m) <- list(c(endo, "Overall"), c("Base Model Loss", "Alt Model Loss", "Diff", "Boot T value", "Boot P Value"))
  attr(m, "comment") <- paste0(
    "CVPAT as per Sharma, Liengaard, Hair, Sarstedt, & Ringle, (2023).\n",
    "  Both models under comparison have identical endogenous constructs with identical measurement models.\n",
    "  Purely exogenous constructs can differ in regards to their relationships with both nomological\n",
    "  partners and measurement indicators."
  )
  class(m) <- c("table_output", "matrix", "array")
  m
}

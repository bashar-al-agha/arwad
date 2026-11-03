# PLSpredict: cross-validated predictions of the indicators of a model's
# endogenous constructs (Shmueli, Ray, Velasquez Estrada and Chatla, 2016;
# Shmueli et al., 2019), as seminr::predict_pls() computes them. Each rule
# below was matched to seminr 2.6.0 by black-box tests on its outputs (help
# pages and results only; docs/adr/0008-pls-predict.md):
# - each repetition draws sample(n) and cuts the shuffled cases into the k
#   consecutive folds of cut(seq(1, n), breaks = k) (leave-one-out: k = n,
#   shuffled as well); the training set of fold f is the shuffled data
#   without the fold, in shuffled order, estimated with the model's own
#   maxIt and stop criterion;
# - a fold model predicts a case from its indicators standardised with the
#   training means and the standard deviations of scale(): each construct's
#   composite scores, one matrix-vector product over its own indicators,
#   then each endogenous construct, after its antecedents, one over its
#   antecedents (in the model's construct order), from their predicted
#   scores (earliest antecedents) or their composite scores (direct
#   antecedents), then the indicators from the scores and the loadings,
#   scaled back;
# - the linear-model benchmark regresses the indicators of each endogenous
#   construct (.lm.fit() with an intercept, raw training data) on the
#   indicators of every exogenous construct (earliest antecedents, whether
#   it reaches the construct or not) or of the construct's direct
#   antecedents, and predicts with X %*% b, one construct at a time, in
#   sample too;
# - the shapes of those products are seminr's: an optimised BLAS (OpenBLAS
#   on Linux) rounds one wide product, or a product with zero weights,
#   differently from the narrower products, where R's reference BLAS gives
#   the same result for all of them;
# - a case's in-sample prediction is the sum, in long double, of the
#   predictions of the k - 1 fold models it trained, over k - 1; the item
#   predictions of the repetitions are averaged with mean() (rowMeans()
#   differs on Linux), and the composite predictions are those of the last
#   repetition.

# The constructs in an order in which each one follows its antecedents
# (pre[i, j] > 0 when construct i is an antecedent of j); NULL for a
# structural model with a cycle.
aw_pls_order <- function(pre) {
  J <- nrow(pre)
  done <- logical(J)
  out <- integer(0)
  while (length(out) < J) {
    ready <- which(!done & colSums(pre[!done, , drop = FALSE] > 0) == 0)
    if (!length(ready)) return(NULL)
    out <- c(out, ready)
    done[ready] <- TRUE
  }
  out
}

# The training means and standard deviations a fold model standardises
# with, as estimate_pls() stores them (meanData, sdData).
aw_pls_train_scale <- function(Xtr) {
  sc <- scale(Xtr)
  list(mu = attr(sc, "scaled:center"), s = attr(sc, "scaled:scale"))
}

# The average of the repetitions' predictions, element by element: x is a
# list of equal-sized matrices.
aw_pls_rep_mean <- function(x) apply(simplify2array(x), c(1L, 2L), mean)

# One fold model's predictions for the cases in the rows of Xs: the
# indicators, scaled back with the training means mu and standard
# deviations s, and the construct scores with the endogenous ones predicted.
# o is the model's row from the PLS kernel.
aw_pls_fold_predict <- function(o, spec, Xs, mu, s, ea, ord) {
  p <- length(spec$items)
  J <- length(spec$cons)
  w <- o[seq_len(p)]
  L <- matrix(o[p + seq_len(p * J)], p, J)
  B <- matrix(o[p + p * J + seq_len(J * J)], J, J)
  Z <- t((t(Xs) - mu) / s)
  # One matrix-vector product per construct, over its own indicators, and
  # per endogenous construct, over its antecedents: an optimised BLAS
  # (OpenBLAS) rounds these differently from one product with the zeros.
  S <- matrix(0, nrow(Z), J)
  for (j in seq_len(J)) {
    k <- which(spec$blk == j)
    S[, j] <- Z[, k, drop = FALSE] %*% w[k]
  }
  P <- S
  for (j in ord) {
    a <- which(spec$pre[, j] > 0)
    if (length(a)) P[, j] <- (if (ea) P else S)[, a, drop = FALSE] %*% B[a, j]
  }
  list(items = t(t(P %*% t(L)) * s + mu), scores = P)
}

# The linear-model benchmark's regressions: for each, the columns of the
# predictors and of the predicted indicators, and in by the positions in y
# of each endogenous construct's indicators, predicted one construct at a
# time. With earliest antecedents one fit serves every endogenous construct
# (.lm.fit() solves each column of y on its own, so the coefficients are
# those of one fit per construct); with direct antecedents the predictors
# are the indicators of each antecedent in turn, the antecedents in their
# order in the structural model (sm), which sets the rounding of the
# least-squares fit.
aw_pls_lm_sets <- function(spec, ea, sm) {
  endo <- which(colSums(spec$pre > 0) > 0)
  exo <- setdiff(seq_along(spec$cons), endo)
  if (ea) {
    y <- which(spec$blk %in% endo)
    return(list(list(x = which(spec$blk %in% exo), y = y, by = unname(split(seq_along(y), spec$blk[y])))))
  }
  lapply(endo, function(j) {
    ants <- match(sm[sm[, "target"] == spec$cons[j], "source"], spec$cons)
    y <- which(spec$blk == j)
    list(x = unlist(lapply(ants, function(a) which(spec$blk == a))), y = y, by = list(seq_along(y)))
  })
}

#' Cross-validated predictions of a PLS-SEM model, compiled
#'
#' A compiled replacement for `seminr::predict_pls()` (PLSpredict; Shmueli et
#' al., 2016, 2019). It takes a model estimated by `seminr::estimate_pls()`,
#' re-estimates it on the training sets of a k-fold cross-validation with
#' the compiled PLS-SEM engine of [aw_boot_pls()], and predicts the
#' indicators of the endogenous constructs of each held-out case, both with
#' the PLS model and with the linear-model benchmark. The folds are drawn
#' exactly as `predict_pls()` draws them, from the same random-number
#' stream: each repetition shuffles the cases with `sample(n)` and cuts them
#' into `noFolds` consecutive folds. The result is the object `predict_pls()`
#' returns, bit for bit, so `summary()` reports the same RMSE and MAE.
#'
#' Each training set is estimated with the model's own `maxIt` and stop
#' criterion, as `predict_pls()` does (a bootstrap uses seminr's defaults).
#' In-sample predictions average, for each case, the predictions of the
#' `noFolds - 1` models whose training set contained it; the item predictions
#' of the repetitions are averaged, and the composite predictions are those
#' of the last repetition, as in seminr. With `cores` above 1, the training
#' sets are estimated in that many worker processes; the results do not
#' depend on it.
#'
#' Models of composites (Mode A and Mode B) are native. Models with
#' `reflective()` constructs, which seminr 2.6.0 predicts with the
#' model-implied conditional expectation of de Rooij et al. (2023), models
#' with interaction terms, and the cases the PLS-SEM engine leaves to seminr
#' are passed to `seminr::predict_pls()` itself. So is every call under
#' seminr versions before 2.6.0, which repeated the same folds in each
#' repetition and could loop forever with earliest antecedents on some
#' structural models.
#'
#' @param model A model returned by `seminr::estimate_pls()`.
#' @param technique `seminr::predict_DA` (direct antecedents, the default) or
#'   `seminr::predict_EA` (earliest antecedents).
#' @param noFolds The number of folds; `NULL` for leave-one-out
#'   cross-validation.
#' @param reps The number of repetitions, each with new random folds; `NULL`
#'   for one.
#' @param cores The number of processes to estimate the training sets in;
#'   `NULL` for one.
#' @return An object of class `predict_pls_model`, laid out as
#'   `seminr::predict_pls()` returns it; `native` is `FALSE` when seminr
#'   produced it.
#' @export
#' @examplesIf requireNamespace("seminr", quietly = TRUE)
#' library(seminr)
#' mm <- constructs(composite("Image", multi_items("IMAG", 1:5)),
#'                  composite("Value", multi_items("PERV", 1:2)),
#'                  composite("Satisfaction", multi_items("CUSA", 1:3)))
#' sm <- relationships(paths(from = "Image", to = c("Value", "Satisfaction")),
#'                     paths(from = "Value", to = "Satisfaction"))
#' model <- estimate_pls(mobi, mm, sm)
#' set.seed(123)
#' pred <- aw_predict_pls(model, technique = predict_EA, noFolds = 10, reps = 10)
#' summary(pred)
aw_predict_pls <- function(model, technique = seminr::predict_DA, noFolds = NULL, reps = NULL, cores = NULL) {
  if (!requireNamespace("seminr", quietly = TRUE)) stop("aw_predict_pls() needs the seminr package")
  delegate <- function() {
    out <- seminr::predict_pls(model, technique = technique, noFolds = noFolds, reps = reps, cores = cores)
    out$native <- FALSE
    out
  }
  ea <- identical(technique, seminr::predict_EA)
  whole <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v == round(v)
  # The rules above are seminr 2.6.0's; earlier versions are left to seminr.
  if (!aw_seminr_260() || !inherits(model, "pls_model") || (!ea && !identical(technique, seminr::predict_DA))) return(delegate())
  spec <- aw_pls_spec(model)
  # Consistent PLS and interaction terms are left to seminr.
  if (is.null(spec) || any(spec$mode >= 3L) || !is.null(spec$stage1) || any(spec$ix > 0L)) return(delegate())
  n <- nrow(spec$X)
  k <- if (is.null(noFolds)) n else noFolds
  R <- if (is.null(reps)) 1L else reps
  ord <- aw_pls_order(spec$pre)
  if (!whole(k) || k < 2 || k > n || !whole(R) || R < 1 || is.null(ord)) return(delegate())
  k <- as.integer(k)
  R <- as.integer(R)
  cores <- if (is.null(cores)) 1L else cores
  X <- spec$X
  p <- ncol(X)
  J <- length(spec$cons)
  # The predicted indicators: those of each endogenous construct, the
  # constructs in their order as targets in the structural model.
  yi <- unlist(lapply(match(unique(model$smMatrix[, "target"]), spec$cons), function(j) which(spec$blk == j)))
  q <- length(yi)
  sets <- aw_pls_lm_sets(spec, ea, model$smMatrix)
  # The folds, from the caller's random-number stream; restored when the
  # work goes to seminr, which then draws the same folds.
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  back <- function() {
    if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
    delegate()
  }
  perms <- lapply(seq_len(R), function(r) sample(n))
  tps <- split(seq_len(n), cut(seq(1, n), breaks = k, labels = FALSE))
  if (length(tps) != k) return(back())
  # Every training set, estimated by the compiled kernel: one call per
  # training-set size (folds differ by one case at most).
  ntr <- n - lengths(tps)
  fits <- vector("list", R * k)
  for (sz in unique(ntr)) {
    fs <- which(ntr == sz)
    at <- as.vector(outer(fs, (seq_len(R) - 1L) * k, `+`))
    idx <- matrix(0L, sz, length(at))
    for (a in seq_along(at)) idx[, a] <- perms[[(at[a] - 1L) %/% k + 1L]][-tps[[(at[a] - 1L) %% k + 1L]]]
    res <- aw_pls_run_cores(spec, idx, cores)
    if (is.null(res) || any(res[, ncol(res)] != 0)) return(back())
    for (a in seq_along(at)) fits[[at[a]]] <- res[a, ]
  }
  reps_out <- vector("list", R)
  for (r in seq_len(R)) {
    oos <- matrix(NA_real_, n, q)
    oos_lm <- oos
    oos_sc <- matrix(NA_real_, n, J)
    ins <- array(NA_real_, c(n, q, k))
    ins_lm <- ins
    ins_sc <- array(NA_real_, c(n, J, k))
    for (f in seq_len(k)) {
      te <- perms[[r]][tps[[f]]]
      tr <- perms[[r]][-tps[[f]]]
      Xtr <- X[tr, , drop = FALSE]
      Xte <- X[te, , drop = FALSE]
      sc <- aw_pls_train_scale(Xtr)
      mu <- sc$mu
      s <- sc$s
      o <- fits[[(r - 1L) * k + f]]
      pin <- aw_pls_fold_predict(o, spec, Xtr, mu, s, ea, ord)
      pout <- aw_pls_fold_predict(o, spec, Xte, mu, s, ea, ord)
      oos[te, ] <- pout$items[, yi]
      oos_sc[te, ] <- pout$scores
      ins[tr, , f] <- pin$items[, yi]
      ins_sc[tr, , f] <- pin$scores
      for (st in sets) {
        A <- cbind(1, Xtr[, st$x, drop = FALSE])
        fit <- .lm.fit(A, Xtr[, st$y, drop = FALSE])
        if (fit$rank < ncol(A)) return(back())
        Ate <- cbind(1, Xte[, st$x, drop = FALSE])
        b <- as.matrix(fit$coefficients)
        for (g in st$by) {
          cols <- match(st$y[g], yi)
          oos_lm[te, cols] <- Ate %*% b[, g, drop = FALSE]
          ins_lm[tr, cols, f] <- A %*% b[, g, drop = FALSE]
        }
      }
    }
    reps_out[[r]] <- list(oos = oos, ins = rowSums(ins, dims = 2L, na.rm = TRUE) / (k - 1L), oos_lm = oos_lm,
                          ins_lm = rowSums(ins_lm, dims = 2L, na.rm = TRUE) / (k - 1L), oos_sc = oos_sc,
                          ins_sc = rowSums(ins_sc, dims = 2L, na.rm = TRUE) / (k - 1L))
  }
  avg <- function(el) aw_pls_rep_mean(lapply(reps_out, `[[`, el))
  items <- spec$items[yi]
  actual <- model$data[, items, drop = FALSE]
  # The data's row names label the predictions of a single run (reps =
  # NULL), and none once repetitions are asked for, even one.
  rn <- if (is.null(reps)) rownames(model$data)
  named <- function(m) {
    dimnames(m) <- list(rn, items)
    m
  }
  pls_oos <- named(avg("oos"))
  pls_ins <- named(avg("ins"))
  lm_oos <- named(avg("oos_lm"))
  lm_ins <- named(avg("ins_lm"))
  last <- reps_out[[R]]
  cs <- function(m) {
    dimnames(m) <- dimnames(model$construct_scores)
    m
  }
  # With a single predicted indicator seminr's matrices drop to vectors:
  # named by the data's row names for a single run; unnamed with
  # repetitions, where the benchmark's predictions stay n x 1 matrices. The
  # residuals then follow from the indicator's values minus those.
  if (q == 1L) {
    actual <- actual[[1L]]
    pls_oos <- pls_oos[, 1L]
    pls_ins <- pls_ins[, 1L]
    if (is.null(reps)) {
      lm_oos <- lm_oos[, 1L]
      lm_ins <- lm_ins[, 1L]
    }
  }
  out <- list(composites = list(composite_out_of_sample = cs(last$oos_sc), composite_in_sample = cs(last$ins_sc),
                                actuals_star = model$construct_scores),
              items = list(PLS_out_of_sample = pls_oos, PLS_in_sample = pls_ins, lm_out_of_sample = lm_oos, lm_in_sample = lm_ins,
                           item_actuals = model$data, PLS_out_of_sample_residuals = actual - pls_oos,
                           PLS_in_sample_residuals = actual - pls_ins, lm_out_of_sample_residuals = actual - lm_oos,
                           lm_in_sample_residuals = actual - lm_ins),
              model = model)
  class(out) <- "predict_pls_model"
  out$native <- TRUE
  out
}

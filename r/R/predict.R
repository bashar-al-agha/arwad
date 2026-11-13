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
#   differs on macOS arm64, which has no long double), and the composite
#   predictions are those of the last repetition.

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

# The average of the repetitions' predictions, element by element, as
# mean() computes it (x is a list of equal-sized matrices): with one
# repetition the predictions themselves plus 0 (mean() turns -0 into 0);
# with several, a compiled loop calling mean() once per element, which runs
# GNU R's algorithm, extended precision included, at C speed.
aw_pls_rep_mean <- function(x) {
  if (length(x) == 1L) return(x[[1L]] + 0)
  a <- simplify2array(x)
  d <- dim(a)
  dim(a) <- c(d[1L] * d[2L], d[3L])
  matrix(aw_pls_rep_mean_unit()(t(a)), d[1L], d[2L])
}

# The mean of each column of M; compiled on first use.
aw_pls_rep_mean_kernel <- function(M) {
  nc <- ncol(M)
  z <- numeric(nc)
  for (j in seq_len(nc)) z[j] <- mean(M[, j])
  z
}

aw_pls_rep_mean_unit <- function() {
  if (is.null(aw_state$rep_mean)) aw_state$rep_mean <- aw_compile(aw_pls_rep_mean_kernel, name = "aw_pls_rep_mean")
  aw_state$rep_mean
}

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

# Validates a model and cross-validation arguments for native PLSpredict/CVPAT,
# draws the R fold permutations from the caller's RNG stream, and estimates all
# R * k training sets in batches by training-set size via aw_pls_run_cores().
# Restores .Random.seed and returns NULL when any check or fold fails.
aw_pls_cv_fits <- function(model, technique, noFolds, reps, cores, need_lm = TRUE) {
  ea <- identical(technique, seminr::predict_EA)
  whole <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v == round(v)
  if (!aw_seminr_260() || !inherits(model, "pls_model") || isTRUE(model$hoc) ||
      (!ea && !identical(technique, seminr::predict_DA))) return(NULL)
  spec <- aw_pls_spec(model)
  if (is.null(spec) || any(spec$mode >= 3L) || !is.null(spec$stage1) || any(spec$ix > 0L)) return(NULL)
  n <- nrow(spec$X)
  k <- if (is.null(noFolds)) n else noFolds
  R <- if (is.null(reps)) 1L else reps
  ord <- aw_pls_order(spec$pre)
  if (!whole(k) || k < 2 || k > n || !whole(R) || R < 1 || is.null(ord)) return(NULL)
  k <- as.integer(k)
  R <- as.integer(R)
  cores <- if (is.null(cores)) 1L else cores
  if (!whole(cores) || cores < 1) return(NULL)
  endo <- unique(model$smMatrix[, "target"])
  itc <- lapply(endo, function(c) which(spec$blk == match(c, spec$cons)))
  yi <- unlist(itc, use.names = FALSE)
  sets <- if (need_lm) aw_pls_lm_sets(spec, ea, model$smMatrix)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  restore <- function() {
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    NULL
  }
  perms <- lapply(seq_len(R), function(r) sample(n))
  tps <- split(seq_len(n), cut(seq(1, n), breaks = k, labels = FALSE))
  if (length(tps) != k) return(restore())
  ntr <- n - lengths(tps)
  fits <- vector("list", R * k)
  for (sz in unique(ntr)) {
    fs <- which(ntr == sz)
    at <- as.vector(outer(fs, (seq_len(R) - 1L) * k, `+`))
    idx <- matrix(0L, sz, length(at))
    for (a in seq_along(at)) idx[, a] <- perms[[(at[a] - 1L) %/% k + 1L]][-tps[[(at[a] - 1L) %% k + 1L]]]
    res <- aw_pls_run_cores(spec, idx, cores)
    if (is.null(res) || any(res[, ncol(res)] != 0)) return(restore())
    for (a in seq_along(at)) fits[[at[a]]] <- res[a, ]
  }
  list(spec = spec, ea = ea, ord = ord, n = n, k = k, R = R, endo = endo, itc = itc, yi = yi,
       q = length(yi), sets = sets, perms = perms, tps = tps, fits = fits, restore = restore)
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
#'   produced it, and `q2_predict` holds the cross-validated \eqn{Q^2_{\text{predict}}}
#'   of each predicted indicator when `native` is `TRUE` (see [aw_q2_predict()]).
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
#' pred$q2_predict
aw_predict_pls <- function(model, technique = seminr::predict_DA, noFolds = NULL, reps = NULL, cores = NULL) {
  if (!requireNamespace("seminr", quietly = TRUE)) stop("aw_predict_pls() needs the seminr package")
  delegate <- function() {
    out <- seminr::predict_pls(model, technique = technique, noFolds = noFolds, reps = reps, cores = cores)
    # seminr returns NULL for models it does not predict (higher-order models).
    if (is.null(out)) return(out)
    out$native <- FALSE
    out
  }
  cv <- aw_pls_cv_fits(model, technique, noFolds, reps, cores, need_lm = TRUE)
  if (is.null(cv)) return(delegate())
  spec <- cv$spec
  ea <- cv$ea
  ord <- cv$ord
  n <- cv$n
  k <- cv$k
  R <- cv$R
  yi <- cv$yi
  q <- cv$q
  sets <- cv$sets
  perms <- cv$perms
  tps <- cv$tps
  fits <- cv$fits
  back <- function() {
    cv$restore()
    delegate()
  }
  X <- spec$X
  J <- length(spec$cons)
  reps_out <- vector("list", R)
  for (r in seq_len(R)) {
    oos <- matrix(NA_real_, n, q)
    oos_lm <- oos
    oos_ia <- oos
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
      oos_ia[te, ] <- rep(mu[yi], each = length(te))
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
                          ins_lm = rowSums(ins_lm, dims = 2L, na.rm = TRUE) / (k - 1L), oos_ia = oos_ia, oos_sc = oos_sc,
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
  ia_oos <- avg("oos_ia")
  Xyi <- X[, yi, drop = FALSE]
  q2_predict <- stats::setNames(1 - colSums((Xyi - pls_oos)^2) / colSums((Xyi - ia_oos)^2), items)
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
  out$q2_predict <- q2_predict
  out
}

#' Cross-validated \eqn{Q^2_{\text{predict}}} for PLSpredict
#'
#' Computes \eqn{Q^2_{\text{predict}}} for each predicted endogenous indicator
#' as defined by Shmueli et al. (2016, 2019) and Hair, Hult, Ringle and
#' Sarstedt (2022, pp. 201-209):
#' \deqn{Q^2_{\text{predict}} = 1 - \frac{\sum_i (y_i - \hat{y}_{i,\text{PLS}})^2}{\sum_i (y_i - \bar{y}_{i,\text{train}})^2}}
#' where \eqn{\hat{y}_{i,\text{PLS}}} is the out-of-sample PLS prediction
#' (averaged across repetitions when `reps > 1`) and \eqn{\bar{y}_{i,\text{train}}}
#' is the indicator's mean in the training fold that held case \eqn{i} out
#' (averaged across repetitions in the same way).
#'
#' When `x` is a native prediction object returned by [aw_predict_pls()], its
#' pre-computed `q2_predict` vector is returned directly. When `x` is a model
#' returned by `seminr::estimate_pls()`, cross-validated predictions are run
#' and \eqn{Q^2_{\text{predict}}} is returned.
#'
#' @param x A `predict_pls_model` returned by [aw_predict_pls()], or a model
#'   returned by `seminr::estimate_pls()`.
#' @param technique `seminr::predict_DA` (direct antecedents, the default) or
#'   `seminr::predict_EA` (earliest antecedents), used when `x` is a model.
#' @param noFolds The number of folds (`NULL` for leave-one-out), used when
#'   `x` is a model.
#' @param reps The number of repetitions (`NULL` for one), used when `x` is a
#'   model.
#' @param cores The number of worker processes (`NULL` for one), used when `x`
#'   is a model.
#' @return A named numeric vector of \eqn{Q^2_{\text{predict}}} values, one per
#'   predicted endogenous indicator; `NULL` for a higher-order model, which
#'   `seminr::predict_pls()` does not predict.
#' @references
#' Hair, J. F., Hult, G. T. M., Ringle, C. M., & Sarstedt, M. (2022). *A primer
#' on partial least squares structural equation modeling (PLS-SEM)* (3rd ed.).
#' SAGE.
#'
#' Shmueli, G., Ray, S., Velasquez Estrada, J. M., & Chatla, S. B. (2016). The
#' elephant in the room: Predictive performance of PLS models. *Journal of
#' Business Research, 69*(10), 4552-4564.
#'
#' Shmueli, G., Sarstedt, M., Hair, J. F., Cheah, J.-H., Ting, H., Vaithilingam,
#' S., & Ringle, C. M. (2019). Predictive model assessment in PLS-SEM:
#' Guidelines for using PLSpredict. *European Journal of Marketing, 53*(11),
#' 2322-2347.
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
#' aw_q2_predict(model, technique = predict_EA, noFolds = 10, reps = 10)
aw_q2_predict <- function(x, technique = seminr::predict_DA, noFolds = NULL, reps = NULL, cores = NULL) {
  if (inherits(x, "predict_pls_model")) {
    if (!is.null(x$q2_predict)) return(x$q2_predict)
    stop("this 'predict_pls_model' does not carry training-fold indicator means (it was produced by seminr::predict_pls() or delegated to it); pass the 'pls_model' to aw_q2_predict() directly",
         call. = FALSE)
  }
  if (!inherits(x, "pls_model")) {
    stop("'x' must be a 'predict_pls_model' from aw_predict_pls() or a 'pls_model' from seminr::estimate_pls()",
         call. = FALSE)
  }
  # The state the folds are drawn from, kept so that a delegated prediction's
  # folds can be drawn again below. A session without one is initialised
  # here, as its first draw would do; otherwise seminr's draw and the redraw
  # would each initialise it from a different clock time.
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) set.seed(NULL)
  old <- get(".Random.seed", envir = globalenv())
  pred <- aw_predict_pls(x, technique = technique, noFolds = noFolds, reps = reps, cores = cores)
  if (is.null(pred) || is.null(pred$items)) return(NULL)
  if (!is.null(pred$q2_predict)) return(pred$q2_predict)
  # When aw_predict_pls() delegated to seminr::predict_pls() (e.g. consistent
  # PLS or interaction terms under seminr >= 2.6.0, or any model under seminr
  # < 2.6.0), reconstruct the same folds from the starting seed to obtain the
  # training-fold indicator means, then restore the post-prediction .Random.seed.
  after <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(after)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", after, envir = globalenv()),
          add = TRUE)
  assign(".Random.seed", old, envir = globalenv())
  endo <- unique(x$smMatrix[, "target"])
  items <- unlist(lapply(endo, function(c) x$mmMatrix[x$mmMatrix[, "construct"] == c, "measurement"]), use.names = FALSE)
  Xyi <- as.matrix(x$data[, items, drop = FALSE])
  n <- nrow(Xyi)
  q <- ncol(Xyi)
  k <- if (is.null(noFolds)) n else as.integer(noFolds)
  R <- if (is.null(reps)) 1L else as.integer(reps)
  perms <- if (aw_seminr_260()) lapply(seq_len(R), function(r) sample(n)) else rep(list(sample(n)), R)
  tps <- split(seq_len(n), cut(seq(1, n), breaks = k, labels = FALSE))
  ia_reps <- lapply(seq_len(R), function(r) {
    oos_ia <- matrix(NA_real_, n, q)
    for (f in seq_len(k)) {
      te <- perms[[r]][tps[[f]]]
      tr <- perms[[r]][-tps[[f]]]
      # The training means, as scale() centres (colMeans()).
      oos_ia[te, ] <- rep(colMeans(Xyi[tr, , drop = FALSE]), each = length(te))
    }
    oos_ia
  })
  ia_oos <- aw_pls_rep_mean(ia_reps)
  res_pls <- as.matrix(pred$items$PLS_out_of_sample_residuals)
  stats::setNames(1 - colSums(res_pls^2) / colSums((Xyi - ia_oos)^2), items)
}

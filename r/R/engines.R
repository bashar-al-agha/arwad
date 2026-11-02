# Statistical engines (Phase 2): functions written in the compilable subset,
# shipped with the package and compiled on first use. Each returns exactly
# what GNU R returns when it runs the same algorithm with the same seed, and
# is validated on at least 50 datasets (tests/testthat/test-engines.R).

# The compiled kernel of aw_boot_ols(): case resampling and LINPACK QR.
aw_boot_ols_kernel <- function(X, y, B) {
  n <- nrow(X)
  k <- ncol(X)
  out <- matrix(0, B, k)
  Xb <- matrix(0, n, k)
  yb <- numeric(n)
  for (b in seq_len(B)) {
    idx <- sample.int(n, replace = TRUE)
    for (i in seq_len(n)) {
      yb[i] <- y[idx[i]]
      for (j in seq_len(k)) Xb[i, j] <- X[idx[i], j]
    }
    cf <- .lm.fit(Xb, yb)$coefficients
    for (j in seq_len(k)) out[b, j] <- cf[j]
  }
  out
}

#' Bootstrap OLS coefficients, compiled
#'
#' Case-resampling (nonparametric) bootstrap of least-squares coefficients.
#' Each replicate draws `nrow(X)` rows with replacement
#' (`sample.int(n, replace = TRUE)`) and fits them with [stats::.lm.fit()],
#' R's own LINPACK QR, which [stats::lm()] uses too. The loop runs as native
#' code when the Arwad toolchain is available ([aw_toolchain()]) and in GNU R
#' otherwise; for the same seed both give identical results. A replicate whose
#' resampled design is rank-deficient gets the coefficients `.lm.fit()` gives
#' (pivoted).
#'
#' @param X A numeric design matrix; include a column of ones for an intercept.
#' @param y A numeric response vector of length `nrow(X)`.
#' @param B The number of bootstrap replicates.
#' @return A `B` by `ncol(X)` matrix of coefficients, one row per replicate,
#'   with the column names of `X`.
#' @export
#' @examples
#' X <- cbind(intercept = 1, wt = mtcars$wt, hp = mtcars$hp)
#' set.seed(1)
#' b <- aw_boot_ols(X, mtcars$mpg, B = 200)
#' apply(b, 2, sd) # bootstrap standard errors
aw_boot_ols <- function(X, y, B = 1000L) {
  X <- as.matrix(X)
  cn <- colnames(X)
  Xd <- matrix(as.double(X), nrow(X), ncol(X))
  yd <- as.double(y)
  if (length(yd) != nrow(Xd)) stop("length(y) must equal nrow(X)")
  if (is.null(aw_state$boot_ols)) aw_state$boot_ols <- aw_compile(aw_boot_ols_kernel, name = "aw_boot_ols")
  out <- aw_state$boot_ols(Xd, yd, as.integer(B))
  colnames(out) <- cn
  out
}

# The compiled kernel of aw_boot_pls(): the PLS path algorithm (Wold 1982;
# Lohmoller 1989) with Mode A or Mode B outer weights and the path weighting
# scheme, for every column of `idx` (row indices of one sample each). Each
# numeric step was matched to seminr 2.5.0 and 2.6.0 by black-box tests on
# their outputs; see docs/adr/0006-pls-sem-engine.md. blk[k] is the construct
# of indicator k and pos[k] its place among the nb[blk[k]] indicators of that
# construct; mode[j] is 1 for Mode A, 2 for Mode B, 3 for a reflective
# construct (Mode A weights, then consistent PLS), 4 for a Mode A composite
# that seminr < 2.6.0 corrects too; pre[i, j] = m > 0 when construct i is
# the m-th antecedent of j. One row per sample: the weights (one per
# indicator), the loadings (p x J, by column; a cross-loading is the
# covariance times 0, so it keeps seminr's signed zero), path coefficients,
# total effects and HTMT (J x J each, by column), iterations, the last weight
# change, and a status: 1 when consistent PLS is undefined for the sample (a
# rho_A that is not positive and finite), 2 when its solution is
# inadmissible, as seminr 2.6.0 checks it.
# Two-stage interaction terms (Henseler and Chin, 2010): ix[k] = m > 0 when
# column k of X is the single indicator of the m-th term (a placeholder in X;
# the kernel fills it for each sample) and xj[j] = m when construct j is that
# term. The first stage is the model without its interaction terms; row r of
# W1 holds its weights for sample r (one per first-stage indicator, from a
# run of this kernel on it), k1[a] is the column of X of its a-th indicator
# and c1[a] that indicator's first-stage construct, and iv1[m] and md1[m]
# are the first-stage constructs whose scores the m-th term multiplies.
# wonly = 1 (the first-stage run) returns the weights alone, one row per
# sample, and skips the rest of the estimation, which the second stage
# does not use.
aw_pls_kernel <- function(X, blk, pos, nb, mode, pre, maxIt, crit, idx, ix, xj, k1, c1, iv1, md1, W1, wonly) {
  # The sample size is the number of rows of idx: n for a bootstrap
  # resample, n - 1 for a leave-one-out (jackknife) sample.
  n <- nrow(idx)
  p <- ncol(X)
  J <- nrow(pre)
  R <- ncol(idx)
  JJ <- J * J
  pJ <- p * J
  width <- p + pJ + 3L * JJ + 3L
  if (wonly == 1L) width <- p
  out <- matrix(0, R, width)
  Xb <- matrix(0, n, p)
  Z <- matrix(0, n, p)
  Ci <- matrix(0, p, p)
  W <- matrix(0, p, J)
  Wn <- matrix(0, p, J)
  Y <- matrix(0, n, J)
  Yt <- matrix(0, n, J)
  E <- matrix(0, J, J)
  Bm <- matrix(0, J, J)
  H <- matrix(0, J, J)
  Rc <- matrix(0, J, J)
  Sg <- matrix(0, J, J)
  Lc <- matrix(0, J, J)
  dn <- integer(J)
  rho <- numeric(J)
  wwv <- numeric(J)
  plsc <- 0L
  for (j in seq_len(J)) if (mode[j] == 3L) plsc <- 1L
  ni <- 0L
  for (k in seq_len(p)) if (ix[k] > ni) ni <- ix[k]
  p1 <- length(k1)
  J1 <- 0L
  for (a in seq_len(p1)) if (c1[a] > J1) J1 <- c1[a]
  xs <- numeric(ni)
  for (r in seq_len(R)) {
    for (k in seq_len(p)) if (ix[k] == 0L) for (i in seq_len(n)) Xb[i, k] <- X[idx[i, r], k]
    # The indicators, standardised as scale() does.
    cm <- colMeans(Xb)
    for (k in seq_len(p)) if (ix[k] == 0L) {
      v <- Xb[, k] - cm[k]
      s <- sqrt(sum(v^2) / (n - 1))
      for (i in seq_len(n)) Z[i, k] <- v[i] / s
    }
    # Each interaction term's indicator: the product of the first stage's
    # scores of its two constructs (the first-stage indicators, standardised,
    # times its weights), standardised in turn. xs[m] keeps the product's
    # standard deviation, sd() computed as sqrt(cov(x, x)).
    if (ni > 0L) {
      Z1 <- matrix(0, n, p1)
      W1m <- matrix(0, p1, J1)
      for (a in seq_len(p1)) {
        for (i in seq_len(n)) Z1[i, a] <- Z[i, k1[a]]
        W1m[a, c1[a]] <- W1[r, a]
      }
      S1 <- Z1 %*% W1m
      for (k in seq_len(p)) if (ix[k] > 0L) for (i in seq_len(n)) Xb[i, k] <- S1[i, iv1[ix[k]]] * S1[i, md1[ix[k]]]
      cm <- colMeans(Xb)
      for (k in seq_len(p)) if (ix[k] > 0L) {
        v <- Xb[, k] - cm[k]
        s <- sqrt(sum(v^2) / (n - 1))
        for (i in seq_len(n)) Z[i, k] <- v[i] / s
        xs[ix[k]] <- sqrt(cov(Xb[, k], Xb[, k]))
      }
    }
    # Mode B: the inverse of each such construct's indicator correlations,
    # which stays the same in every iteration of a sample (Ci[k, l] for two
    # indicators of one construct).
    for (j in seq_len(J)) if (mode[j] == 2L) {
      Zj <- matrix(0, n, nb[j])
      for (k in seq_len(p)) if (blk[k] == j) for (i in seq_len(n)) Zj[i, pos[k]] <- Z[i, k]
      Cj <- solve(cor(Zj))
      for (k in seq_len(p)) if (blk[k] == j) for (l in seq_len(p)) if (blk[l] == j) Ci[k, l] <- Cj[pos[k], pos[l]]
    }
    for (k in seq_len(p)) for (j in seq_len(J)) W[k, j] <- if (blk[k] == j) 1 else 0
    it <- 0L
    dw <- 0
    for (iter in 0:maxIt) {
      # Outer approximation: standardised construct scores.
      S <- Z %*% W
      cm <- colMeans(S)
      for (j in seq_len(J)) {
        v <- S[, j] - cm[j]
        s <- sqrt(sum(v^2) / (n - 1))
        for (i in seq_len(n)) Y[i, j] <- v[i] / s
      }
      # Inner weights, path weighting scheme: correlations with successors,
      # regression coefficients on antecedents.
      for (j in seq_len(J)) {
        na <- 0L
        for (i in seq_len(J)) {
          E[i, j] <- 0
          if (pre[i, j] > 0) na <- na + 1L
          if (pre[j, i] > 0) E[i, j] <- cor(Y[, i], Y[, j])
        }
        if (na > 0L) {
          Yp <- matrix(0, n, na)
          for (i in seq_len(J)) {
            q <- pre[i, j]
            if (q > 0) for (t in seq_len(n)) Yp[t, q] <- Y[t, i]
          }
          Ypt <- t(Yp)
          bj <- drop(solve(Ypt %*% Yp, Ypt %*% Y[, j]))
          for (i in seq_len(J)) if (pre[i, j] > 0) E[i, j] <- bj[pre[i, j]]
        }
      }
      # Inner proxies, standardised. Mode A weights are the indicators'
      # covariances with the proxy; Mode B weights the coefficients of the
      # proxy regressed on the indicators, solve(cor(Zj)) %*% cor(Zj, proxy).
      # Then unit-variance scores.
      S <- Y %*% E
      cm <- colMeans(S)
      for (j in seq_len(J)) {
        v <- S[, j] - cm[j]
        s <- sqrt(sum(v^2) / (n - 1))
        for (i in seq_len(n)) Yt[i, j] <- v[i] / s
      }
      for (k in seq_len(p)) for (j in seq_len(J)) Wn[k, j] <- if (blk[k] == j && mode[j] != 2L) cov(Z[, k], Yt[, j]) else 0
      for (j in seq_len(J)) if (mode[j] == 2L) {
        Cj <- matrix(0, nb[j], nb[j])
        cj <- numeric(nb[j])
        for (k in seq_len(p)) if (blk[k] == j) {
          cj[pos[k]] <- cor(Z[, k], Yt[, j])
          for (l in seq_len(p)) if (blk[l] == j) Cj[pos[k], pos[l]] <- Ci[k, l]
        }
        wj <- drop(Cj %*% cj)
        for (k in seq_len(p)) if (blk[k] == j) Wn[k, j] <- wj[pos[k]]
      }
      S <- Z %*% Wn
      for (j in seq_len(J)) {
        s <- sqrt(sum(S[, j]^2) / (n - 1))
        for (k in seq_len(p)) Wn[k, j] <- Wn[k, j] / s
      }
      dw <- sum(abs(Wn - W))
      W <- Wn
      it <- iter
      if (dw < crit) break
    }
    if (wonly == 1L) {
      for (k in seq_len(p)) out[r, k] <- W[k, blk[k]]
      next
    }
    # Final scores (not rescaled), path coefficients, total effects. seminr
    # gives an interaction term's score the standard deviation of its product
    # (the paths into and out of it, and its loading, follow from that score).
    Y <- Z %*% W
    for (j in seq_len(J)) if (xj[j] > 0L) for (i in seq_len(n)) Y[i, j] <- Y[i, j] * xs[xj[j]]
    # Consistent PLS (Dijkstra and Henseler, 2015) when the model has a
    # reflective construct: rho_A from the final weights w and the indicators'
    # covariances S, (w'w)^2 w'(S - diag S)w / w'(ww' - diag ww')w, for every
    # reflective construct with two or more indicators, and for the Mode A
    # composites that seminr < 2.6.0 corrects too (mode 4); 1 for other
    # composites and single indicators. The paths then come from the score
    # correlations divided by sqrt(rho_i rho_j).
    bad <- 0L
    if (plsc == 1L) {
      for (j in seq_len(J)) {
        rho[j] <- 1
        if (mode[j] == 3L || mode[j] == 4L) {
          wc <- matrix(0, nb[j], 1L)
          for (k in seq_len(p)) if (blk[k] == j) wc[pos[k], 1L] <- W[k, j]
          qm <- t(wc) %*% wc
          wwv[j] <- qm[1L, 1L]
          if (nb[j] > 1L) {
            Zc <- matrix(0, n, nb[j])
            for (k in seq_len(p)) if (blk[k] == j) for (i in seq_len(n)) Zc[i, pos[k]] <- Z[i, k]
            Sc <- cov(Zc)
            Wc <- wc %*% t(wc)
            for (a in seq_len(nb[j])) {
              Sc[a, a] <- 0
              Wc[a, a] <- 0
            }
            qn <- t(wc) %*% Sc %*% wc
            qd <- t(wc) %*% Wc %*% wc
            rho[j] <- wwv[j]^2 * (qn[1L, 1L] / qd[1L, 1L])
            if (is.na(rho[j]) || rho[j] <= 0 || rho[j] == 1 / 0) bad <- 1L
          }
        }
      }
      if (bad == 0L) {
        Rc <- cor(Y)
        for (j in seq_len(J)) for (i in seq_len(J)) Rc[i, j] <- if (i == j) 1 else Rc[i, j] / sqrt(rho[i] * rho[j])
      }
    }
    for (j in seq_len(J)) {
      na <- 0L
      for (i in seq_len(J)) {
        Bm[i, j] <- 0
        if (pre[i, j] > 0) na <- na + 1L
      }
      if (na > 0L && plsc == 1L) {
        if (bad == 0L) {
          Ap <- matrix(0, na, na)
          bp <- matrix(0, na, 1L)
          for (i in seq_len(J)) {
            q <- pre[i, j]
            if (q > 0) {
              bp[q, 1L] <- Rc[i, j]
              for (h in seq_len(J)) if (pre[h, j] > 0) Ap[q, pre[h, j]] <- Rc[i, h]
            }
          }
          bj <- drop(solve(Ap, bp))
          for (i in seq_len(J)) if (pre[i, j] > 0) Bm[i, j] <- bj[pre[i, j]]
        }
      } else if (na > 0L) {
        Yp <- matrix(0, n, na)
        for (i in seq_len(J)) {
          q <- pre[i, j]
          if (q > 0) for (t in seq_len(n)) Yp[t, q] <- Y[t, i]
        }
        Ypt <- t(Yp)
        bj <- drop(solve(Ypt %*% Yp, Ypt %*% Y[, j]))
        for (i in seq_len(J)) if (pre[i, j] > 0) Bm[i, j] <- bj[pre[i, j]]
      }
    }
    # Admissibility of a consistent-PLS solution, as seminr 2.6.0 checks it
    # (from 2.6.0 on, bootstrap_model() leaves inadmissible resamples out),
    # thresholds measured by bisection: a rho_A above one; a loading of a
    # reflective construct above 1 + 5e-9 in absolute value (a single
    # indicator's loading is one up to rounding); a corrected construct
    # correlation of one or more in absolute value; an R-squared outside
    # [0, 1]; or a construct correlation matrix implied by the structural
    # model whose smallest eigenvalue is below 1e-8. That matrix has the
    # corrected correlations of the exogenous constructs, unit variances, and
    # for each endogenous construct (after its antecedents) the correlations
    # its paths imply; its smallest eigenvalue is above 1e-8 exactly when the
    # matrix less 1e-8 on the diagonal has a Cholesky factor.
    adm <- 1L
    if (plsc == 1L && bad == 0L) {
      for (j in seq_len(J)) {
        if (mode[j] == 3L && rho[j] > 1) adm <- 0L
        if (mode[j] == 3L) for (k in seq_len(p)) if (blk[k] == j && abs(sqrt(rho[j]) / wwv[j] * W[k, j]) > 1 + 5e-9) adm <- 0L
        r2 <- 0
        na <- 0L
        for (i in seq_len(J)) {
          if (i != j && abs(Rc[i, j]) >= 1) adm <- 0L
          if (pre[i, j] > 0) {
            na <- na + 1L
            r2 <- r2 + Bm[i, j] * Rc[i, j]
          }
        }
        if (na > 0L && (r2 < 0 || r2 > 1)) adm <- 0L
      }
      if (adm == 1L) {
        for (j in seq_len(J)) {
          dn[j] <- 1L
          for (i in seq_len(J)) if (pre[i, j] > 0) dn[j] <- 0L
        }
        for (j in seq_len(J)) for (i in seq_len(J)) Sg[i, j] <- if (dn[i] == 1L && dn[j] == 1L) Rc[i, j] else 0
        for (pass in seq_len(J)) for (j in seq_len(J)) if (dn[j] == 0L) {
          rdy <- 1L
          for (i in seq_len(J)) if (pre[i, j] > 0 && dn[i] == 0L) rdy <- 0L
          if (rdy == 1L) {
            for (k in seq_len(J)) if (dn[k] == 1L) {
              s <- 0
              for (i in seq_len(J)) if (pre[i, j] > 0) s <- s + Bm[i, j] * Sg[i, k]
              Sg[j, k] <- s
              Sg[k, j] <- s
            }
            Sg[j, j] <- 1
            dn[j] <- 1L
          }
        }
        for (k in seq_len(J)) {
          d <- Sg[k, k] - 1e-8
          for (m in seq_len(k - 1L)) d <- d - Lc[k, m] * Lc[k, m]
          if (!(d > 0)) {
            adm <- 0L
            break
          }
          Lc[k, k] <- sqrt(d)
          for (i in seq_len(J)) if (i > k) {
            s <- Sg[i, k]
            for (m in seq_len(k - 1L)) s <- s - Lc[i, m] * Lc[k, m]
            Lc[i, k] <- s / Lc[k, k]
          }
        }
      }
    }
    Tm <- Bm
    Pw <- Bm
    if (J > 1L) for (m in 2:J) {
      Pw <- Pw %*% Bm
      Tm <- Tm + Pw
    }
    # HTMT (Henseler, Ringle and Sarstedt, 2015) from the indicators: the
    # mean absolute heterotrait correlation over the geometric mean of the
    # two mean absolute monotrait correlations (1 for a single indicator).
    Cx <- cor(Xb, Xb)
    mono <- numeric(J)
    for (j in seq_len(J)) {
      mv <- numeric(0)
      for (b in seq_len(p)) if (blk[b] == j) for (a in seq_len(p)) if (blk[a] == j && a < b) mv <- c(mv, abs(Cx[a, b]))
      mono[j] <- if (length(mv) > 0L) mean(mv) else 1
    }
    for (j in seq_len(J)) for (i in seq_len(J)) {
      H[i, j] <- NA_real_
      if (i < j) {
        hv <- numeric(0)
        for (b in seq_len(p)) if (blk[b] == j) for (a in seq_len(p)) if (blk[a] == i) hv <- c(hv, abs(Cx[a, b]))
        H[i, j] <- mean(hv) / sqrt(mono[i] * mono[j])
      }
    }
    for (k in seq_len(p)) out[r, k] <- W[k, blk[k]]
    # Loadings: covariances with the scores; for a reflective construct the
    # consistent loadings sqrt(rho_A) / (w'w) * w (rho_A = 1 for one indicator).
    # (A sample with status 1 is left out, so its loadings are not needed.)
    L <- cov(Z, Y)
    for (j in seq_len(J)) for (k in seq_len(p)) {
      out[r, p + (j - 1L) * p + k] <- if (blk[k] != j) L[k, j] * 0 else if (mode[j] != 3L) L[k, j] else if (bad == 0L) sqrt(rho[j]) / wwv[j] * W[k, j] else 0
    }
    for (j in seq_len(J)) for (i in seq_len(J)) {
      q <- (j - 1L) * J + i
      out[r, p + pJ + q] <- Bm[i, j]
      out[r, p + pJ + JJ + q] <- Tm[i, j]
      out[r, p + pJ + 2L * JJ + q] <- H[i, j]
    }
    out[r, p + pJ + 3L * JJ + 1L] <- it
    out[r, p + pJ + 3L * JJ + 2L] <- dw
    out[r, p + pJ + 3L * JJ + 3L] <- if (bad == 1L) 1 else if (adm == 0L) 2 else 0
  }
  out
}

# TRUE when the installed seminr is 2.6.0 or later. Its consistent PLS
# corrects only reflective constructs, and its bootstrap_model() also leaves
# out the resamples whose consistent-PLS solution is inadmissible (ADR-0006,
# amendment "seminr 2.6.0"). aw_boot_pls() follows the installed version,
# because its results must equal what bootstrap_model() would return.
aw_seminr_260 <- function() utils::packageVersion("seminr") >= "2.6.0"

# What the PLS kernel needs from a model estimated by seminr::estimate_pls(),
# or NULL when the model is outside the native subset: composites with Mode A
# or Mode B weights and reflective() constructs (type "C" in the mmMatrix,
# estimated by consistent PLS), two-stage interaction terms in a model of
# composites, the path weighting scheme, each indicator in one construct,
# complete numeric data (no product-indicator or orthogonal interaction
# terms, no higher-order constructs).
aw_pls_spec <- function(model) {
  if (!inherits(model, "pls_model")) stop("'seminr_model' must be a model estimated by seminr::estimate_pls()")
  mm <- model$mmMatrix
  sm <- model$smMatrix
  cons <- model$constructs
  # A matrix when every construct has the same number of indicators.
  items <- as.vector(model$mmVariables)
  raw <- model$rawdata
  if (!is.matrix(mm) || !is.matrix(sm) || !all(mm[, "type"] %in% c("A", "B", "C")) || anyDuplicated(mm[, "measurement"]) ||
      !setequal(mm[, "measurement"], items) || !all(mm[, "construct"] %in% cons) || !all(cons %in% mm[, "construct"]) ||
      !all(c(sm[, "source"], sm[, "target"]) %in% cons) || !identical(model$inner_weights, seminr::path_weighting) ||
      !is.data.frame(raw)) return(NULL)
  ints <- aw_pls_interactions(model, mm, cons, items, names(raw))
  if (is.null(ints)) return(NULL)
  obs <- items[!items %in% ints[, "item"]]
  X <- raw[, obs, drop = FALSE]
  if (!all(vapply(X, is.numeric, NA)) || anyNA(X) || nrow(X) < 3L) return(NULL)
  X <- as.matrix(X)
  dimnames(X) <- NULL
  storage.mode(X) <- "double"
  if (!all(is.finite(X))) return(NULL)
  # An interaction term's column is a placeholder: the kernel computes it
  # for each sample.
  Xf <- matrix(0, nrow(X), length(items))
  Xf[, match(obs, items)] <- X
  spec <- aw_pls_spec_parts(mm, sm, cons, items)
  if (is.null(spec)) return(NULL)
  spec <- c(list(X = Xf), spec, list(maxIt = as.integer(model$settings$maxIt), crit = 10^-model$settings$stopCriterion,
                                     cons = cons, items = items))
  spec$ix <- match(items, ints[, "item"], nomatch = 0L)
  spec$xj <- match(cons, ints[, "construct"], nomatch = 0L)
  if (!nrow(ints)) return(spec)
  # The first stage: the model without its interaction terms (and the paths
  # from them), estimated with seminr's defaults, maxIt = 300 and
  # stopCriterion = 7, whatever the model was estimated with, as seminr
  # does. Consistent PLS with interaction terms is left to seminr, and so is
  # a construct that the first stage would leave without a path.
  if (any(spec$mode >= 3L)) return(NULL)
  cons1 <- cons[!cons %in% ints[, "construct"]]
  sm1 <- sm[!sm[, "source"] %in% ints[, "construct"], , drop = FALSE]
  if (!all(cons1 %in% c(sm1[, "source"], sm1[, "target"]))) return(NULL)
  s1 <- aw_pls_spec_parts(mm[mm[, "construct"] %in% cons1, , drop = FALSE], sm1, cons1, obs)
  if (is.null(s1)) return(NULL)
  spec$stage1 <- c(list(X = X), s1, list(maxIt = 300L, crit = 10^-7, cons = cons1, items = obs, ix = integer(length(obs)),
                                         xj = integer(length(cons1)), wonly = 1L))
  spec$k1 <- match(obs, items)
  spec$c1 <- s1$blk
  spec$iv1 <- match(ints[, "iv"], cons1)
  spec$md1 <- match(ints[, "moderator"], cons1)
  spec
}

# The parts of a kernel spec that come from the measurement and structural
# model: each indicator's construct (blk) and place in it (pos), indicators
# per construct (nb), each construct's mode, and the antecedents (pre); NULL
# when a construct mixes types or is a Mode B construct with one indicator.
aw_pls_spec_parts <- function(mm, sm, cons, items) {
  J <- length(cons)
  pre <- matrix(0L, J, J)
  for (j in seq_len(J)) {
    src <- sm[sm[, "target"] == cons[j], "source"]
    pre[match(src, cons), j] <- seq_along(src)
  }
  blk <- match(mm[match(items, mm[, "measurement"]), "construct"], cons)
  types <- lapply(cons, function(cn) unique(mm[mm[, "construct"] == cn, "type"]))
  if (any(lengths(types) != 1L)) return(NULL)
  mode <- 1L + (unlist(types) == "B") + 2L * (unlist(types) == "C")
  # seminr < 2.6.0 corrects the Mode A composites of a model with a
  # reflective construct as well (mode 4); seminr 2.6.0 leaves every
  # composite as it is (rho_A = 1), as Dijkstra and Henseler (2015) do.
  if (any(mode == 3L) && !aw_seminr_260()) mode[mode == 1L] <- 4L
  nb <- tabulate(blk, J)
  # seminr rejects a Mode B construct with one indicator.
  if (any(mode == 2L & nb < 2L)) return(NULL)
  pos <- integer(length(items))
  for (j in seq_len(J)) pos[blk == j] <- seq_len(nb[j])
  list(blk = blk, pos = pos, nb = nb, mode = mode, pre = pre)
}

# The two-stage interaction terms of a model, one row each: the term's
# construct, its single indicator (the construct's name and "_intxn"; not in
# the data), and the two constructs it multiplies, from its name, which
# seminr documents as "iv*moderator". A matrix with no rows for a model
# without interaction terms; NULL for one outside the native subset
# (product-indicator, orthogonal or other interaction terms, or an indicator
# missing from the data that is no such term).
aw_pls_interactions <- function(model, mm, cons, items, vars) {
  out <- matrix(character(), 0L, 4L, dimnames = list(NULL, c("construct", "item", "iv", "moderator")))
  extra <- items[!items %in% vars]
  terms <- Filter(function(e) inherits(e, "interaction"), model$measurement_model)
  if (!length(extra) && !length(terms)) return(out)
  if (length(extra) != length(terms) || !all(vapply(terms, inherits, NA, "two_stage_interaction"))) return(NULL)
  tc <- mm[match(extra, mm[, "measurement"]), "construct"]
  if (anyDuplicated(tc) || !identical(extra, paste0(tc, "_intxn")) || any(tabulate(match(mm[, "construct"], tc), length(tc)) != 1L)) return(NULL)
  main <- cons[!cons %in% tc]
  for (cn in tc) {
    at <- gregexpr("*", cn, fixed = TRUE)[[1L]]
    pairs <- lapply(at[at > 0L], function(a) c(substr(cn, 1L, a - 1L), substr(cn, a + 1L, nchar(cn))))
    pairs <- Filter(function(pr) all(pr %in% main), pairs)
    if (length(pairs) != 1L) return(NULL)
    out <- rbind(out, c(cn, paste0(cn, "_intxn"), pairs[[1L]]))
  }
  out
}

# Calls `kernel` (aw_pls_kernel() or its compiled unit) on the samples in idx;
# W1 holds the first stage's weights for those samples when the model has
# interaction terms.
aw_pls_call <- function(kernel, spec, idx, W1 = NULL) {
  if (is.null(spec$ix)) spec$ix <- integer(length(spec$items))
  if (is.null(spec$xj)) spec$xj <- integer(length(spec$cons))
  if (is.null(W1)) W1 <- matrix(0, 1L, 1L)
  k1 <- if (is.null(spec$k1)) 0L else spec$k1
  c1 <- if (is.null(spec$c1)) 0L else spec$c1
  iv1 <- if (is.null(spec$iv1)) 0L else spec$iv1
  md1 <- if (is.null(spec$md1)) 0L else spec$md1
  wonly <- if (is.null(spec$wonly)) 0L else spec$wonly
  kernel(spec$X, spec$blk, spec$pos, spec$nb, spec$mode, spec$pre, spec$maxIt, spec$crit, idx, spec$ix, spec$xj, k1, c1, iv1,
         md1, W1, wonly)
}

aw_pls_unit <- function() {
  if (is.null(aw_state$pls)) aw_state$pls <- aw_compile(aw_pls_kernel, name = "aw_pls")
  aw_state$pls
}

# Runs the kernel on the samples in idx; NULL when the native run fails (the
# caller then hands the work to seminr). A model with interaction terms is
# estimated twice per sample: its first stage, whose weights then go into
# the second.
aw_pls_run <- function(spec, idx) {
  tryCatch({
    W1 <- NULL
    if (!is.null(spec$stage1)) W1 <- aw_pls_call(aw_pls_unit(), spec$stage1, idx)[, seq_along(spec$k1), drop = FALSE]
    aw_pls_call(aw_pls_unit(), spec, idx, W1)
  }, error = function(e) NULL)
}

# The model re-estimated on its n leave-one-out samples (the jackknife behind
# the acceleration of a BCa interval), with the model's own maxIt and stop
# criterion, as for the original sample. Row i leaves out case i; NULL when
# the native run fails.
aw_pls_jack <- function(spec) {
  n <- nrow(spec$X)
  idx <- matrix(0L, n - 1L, n)
  for (i in seq_len(n)) idx[, i] <- seq_len(n)[-i]
  aw_pls_run(spec, idx)
}

# A socket cluster of `cores` new R processes set up like this one: its
# library paths and every arwad option, with the cache and the Goose
# executable this process uses (however it was found), so the workers load
# the engines built here instead of running them in GNU R.
aw_cluster <- function(cores) {
  cl <- parallel::makeCluster(cores)
  opts <- options()[startsWith(names(options()), "arwad.")]
  opts$arwad.cache_dir <- aw_cache_dir()
  tc <- aw_toolchain()
  if (tc$ok) opts$arwad.goose <- tc$goose
  parallel::clusterCall(cl, function(lib, opts) {
    .libPaths(lib)
    options(opts)
    loadNamespace("arwad")
    NULL
  }, .libPaths(), opts)
  cl
}

# aw_pls_run() over `cores` processes, each taking a contiguous block of the
# resamples; the rows are bound back in order. The engine is built (or found
# in the on-disk cache) here first, so the workers only load it.
aw_pls_run_cores <- function(spec, idx, cores) {
  cores <- max(1L, min(as.integer(cores), ncol(idx)))
  if (cores == 1L || is.null(aw_pls_run(spec, idx[, 1L, drop = FALSE]))) return(aw_pls_run(spec, idx))
  parts <- parallel::splitIndices(ncol(idx), cores)
  if (.Platform$OS.type == "unix") {
    out <- parallel::mclapply(parts, function(cols) aw_pls_run(spec, idx[, cols, drop = FALSE]), mc.cores = cores)
  } else {
    cl <- aw_cluster(cores)
    on.exit(parallel::stopCluster(cl))
    out <- parallel::parLapply(cl, parts, function(cols, spec, idx) {
      get("aw_pls_run", envir = asNamespace("arwad"))(spec, idx[, cols, drop = FALSE])
    }, spec = spec, idx = idx)
  }
  if (any(vapply(out, function(o) !is.matrix(o), NA))) return(NULL)
  do.call(rbind, out)
}

# The descriptive table seminr attaches to a bootstrap: estimate, bootstrap
# mean and bootstrap standard deviation for each column, side by side. All
# rows keep the estimate's row names as they are (a matrix, for weights,
# when every construct has the same number of indicators).
aw_pls_describe <- function(est, boot, rows = NULL, cols = colnames(est)) {
  m <- apply(boot, c(1L, 2L), mean)
  s <- apply(boot, c(1L, 2L), stats::sd)
  rn <- if (is.null(rows)) dimnames(est)[[1L]] else rows
  if (is.null(rows)) rows <- seq_len(nrow(est))
  out <- cbind(est[rows, cols, drop = FALSE], m[rows, cols, drop = FALSE], s[rows, cols, drop = FALSE])
  dimnames(out) <- list(rn, c(paste(cols, "PLS Est."), paste(cols, "Boot Mean"), paste(cols, "Boot SD")))
  class(out) <- c("matrix", "array", "table_output")
  out
}

#' Bootstrap a PLS-SEM model, compiled
#'
#' A compiled replacement for `seminr::bootstrap_model()`. It takes a model
#' estimated by `seminr::estimate_pls()` and re-estimates it on `nboot` case
#' resamples with the PLS path algorithm (path weighting scheme, Mode A
#' weights for `composite()` constructs, Mode B weights for
#' `composite(..., weights = mode_B)` ones, and for `reflective()` constructs
#' Mode A weights followed by consistent PLS, as seminr applies it;
#' standardised data), drawn exactly as `bootstrap_model()` draws them:
#' resample `b` is `sample.int(n, replace = TRUE)` after `set.seed(seed + b)`.
#' Path coefficients, total effects, loadings and weights are bit-identical to
#' `bootstrap_model()`; HTMT ratios agree to within four units in the last
#' place. The caller's random-number state is left as it was. Models outside
#' the native subset (other weighting schemes, product-indicator or
#' orthogonal interaction terms, consistent PLS with interaction terms,
#' higher-order constructs, missing data) are passed to
#' `seminr::bootstrap_model()` itself.
#'
#' Two-stage interaction terms (moderation, seminr's default
#' `interaction_term(iv, moderator, method = two_stage)`, and
#' `quadratic_term()`) are native in models of composites. As in seminr, each
#' resample is estimated twice: first the model without its interaction
#' terms, with `maxIt = 300` and `stopCriterion = 7` (also for the original
#' sample, whatever the model was estimated with); then the full model, in
#' which each term's single indicator is the product of the first stage's
#' scores of its two constructs, and the term's score keeps that product's
#' standard deviation, so its paths are in the units seminr reports.
#'
#' The engine follows the installed version of seminr (2.5.0 or 2.6.0), so
#' that the numbers stay seminr's. Every resample is estimated with seminr's
#' default `maxIt` (300) and `stopCriterion` (7), whatever the model was
#' estimated with, as `bootstrap_model()` does. In a model with a reflective
#' construct, consistent PLS corrects the reflective constructs and treats
#' composites as fully reliable; seminr before 2.6.0 corrects every Mode A
#' construct with two or more indicators, `composite()` ones included.
#' Resamples for which consistent PLS is undefined (a negative rho_A) are
#' left out of the bootstrap and `boots` counts the others, with a message.
#' From seminr 2.6.0 on, so are the resamples whose consistent-PLS solution
#' is inadmissible (a rho_A above one, a loading of a reflective construct
#' above one, a corrected construct correlation of one or more, an R-squared
#' outside \[0, 1\], or a construct correlation matrix implied by the
#' structural model that is not positive definite), and the object also
#' records `boots_requested` and `boots_dropped`. When fewer than two
#' resamples are left, the model goes to `bootstrap_model()`, which then
#' stops.
#'
#' With `cores` above 1 the resamples are split over that many worker
#' processes (forked on Linux and macOS, a socket cluster on Windows). Each
#' resample is computed on its own, so the result is the same for any number
#' of cores; starting the workers takes time, so more cores pay off for
#' thousands of resamples, not for a few hundred.
#'
#' @param seminr_model A model returned by `seminr::estimate_pls()`.
#' @param nboot The number of bootstrap resamples.
#' @param seed The seed; when `NULL`, one is drawn with
#'   `sample.int(100000, 1)`, as `bootstrap_model()` does.
#' @param cores The number of processes to run the resamples in, also when
#'   the model is handed to `seminr::bootstrap_model()`.
#' @return An object of class `boot_seminr_model`, laid out as
#'   `bootstrap_model()` returns it, so `summary()` and seminr's other
#'   methods work on it; `native` is `FALSE` when seminr produced it.
#' @export
#' @examplesIf requireNamespace("seminr", quietly = TRUE)
#' library(seminr)
#' mm <- constructs(composite("Image", multi_items("IMAG", 1:5)),
#'                  composite("Value", multi_items("PERV", 1:2)),
#'                  composite("Satisfaction", multi_items("CUSA", 1:3)))
#' sm <- relationships(paths(from = "Image", to = c("Value", "Satisfaction")),
#'                     paths(from = "Value", to = "Satisfaction"))
#' model <- estimate_pls(mobi, mm, sm)
#' boot <- aw_boot_pls(model, nboot = 200, seed = 123)
#' summary(boot)$bootstrapped_paths
aw_boot_pls <- function(seminr_model, nboot = 500L, seed = NULL, cores = 1L) {
  if (!requireNamespace("seminr", quietly = TRUE)) stop("aw_boot_pls() needs the seminr package")
  if (is.null(seed)) seed <- sample.int(100000, 1)
  # The counts keep the type nboot was given, as in seminr.
  requested <- nboot
  boots <- nboot
  nboot <- as.integer(nboot)
  spec <- aw_pls_spec(seminr_model)
  delegate <- function() {
    out <- seminr::bootstrap_model(seminr_model, nboot = requested, cores = cores, seed = seed)
    out$native <- FALSE
    out
  }
  if (is.null(spec) || nboot < 1L) return(delegate())
  n <- nrow(spec$X)
  # The original sample, for HTMT and total effects; a sample for which
  # consistent PLS is undefined (status 1) goes to seminr, which warns. An
  # inadmissible estimate (status 2) is bootstrapped, as seminr does.
  orig <- aw_pls_run(spec, matrix(seq_len(n), ncol = 1L))
  if (is.null(orig) || orig[1L, ncol(orig)] == 1) return(delegate())
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  idx <- matrix(0L, n, nboot)
  for (b in seq_len(nboot)) {
    set.seed(seed + b)
    idx[, b] <- sample.int(n, replace = TRUE)
  }
  # bootstrap_model() estimates every resample with seminr's defaults,
  # maxIt = 300 and stopCriterion = 7, whatever the model was estimated with.
  bspec <- spec
  bspec$maxIt <- 300L
  bspec$crit <- 10^-7
  res <- aw_pls_run_cores(bspec, idx, cores)
  if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  if (is.null(res)) return(delegate())
  # Like bootstrap_model(), leave out the resamples for which consistent PLS
  # is undefined (seminr 2.5.0 reports them as failing to converge) and, from
  # seminr 2.6.0 on, those whose solution is inadmissible (seminr warns
  # there), and count only the others in `boots`; with fewer than two left,
  # seminr stops with an error and returns the model unbootstrapped, so that
  # case goes to seminr.
  v260 <- aw_seminr_260()
  st <- res[, ncol(res)]
  ok <- if (v260) st == 0 else st != 1
  if (sum(ok) < 2L) return(delegate())
  if (!all(ok)) {
    message(sprintf("aw_boot_pls(): %d of %d resamples left out of the bootstrap statistics: consistent PLS is %s for them, as in seminr::bootstrap_model().",
                    sum(!ok), nboot, if (v260) "undefined (rho_A not positive) or inadmissible" else "undefined (rho_A not positive)"))
    res <- res[ok, , drop = FALSE]
    boots <- boots - sum(!ok)
  }
  # seminr 2.6.0 also records the resamples requested and dropped (in the
  # type nboot was given), before the seed.
  counts <- if (v260) list(boots_requested = requested, boots_dropped = requested - boots)
  aw_pls_boot_object(seminr_model, spec, res, orig, nrow(res), boots, seed, counts)
}

# Unpacks kernel rows into the arrays and tables of a boot_seminr_model;
# orig is the kernel's row for the original sample.
aw_pls_boot_object <- function(model, spec, res, orig, nboot, boots, seed, counts = NULL) {
  cons <- spec$cons
  items <- spec$items
  p <- length(items)
  J <- length(cons)
  JJ <- J * J
  pJ <- p * J
  reps <- as.character(seq_len(nboot))
  sq <- function(off) array(t(res[, off + seq_len(JJ), drop = FALSE]), c(J, J, nboot), dimnames = list(cons, cons, reps))
  weights <- array(0, c(p, J, nboot), dimnames = list(items, cons, reps))
  for (k in seq_len(p)) weights[k, spec$blk[k], ] <- res[, k]
  out <- model
  out$boot_paths <- sq(p + pJ)
  out$boot_loadings <- array(t(res[, p + seq_len(pJ), drop = FALSE]), c(p, J, nboot), dimnames = list(items, cons, reps))
  out$boot_weights <- weights
  out$boot_HTMT <- sq(p + pJ + 2L * JJ)
  out$boot_total_paths <- sq(p + pJ + JJ)
  src <- cons[cons %in% model$smMatrix[, "source"]]
  tgt <- cons[cons %in% model$smMatrix[, "target"]]
  htmt <-matrix(orig[1L, p + pJ + 2L * JJ + seq_len(JJ)], J, J, dimnames = list(cons, cons))
  total <- matrix(orig[1L, p + pJ + JJ + seq_len(JJ)], J, J, dimnames = list(cons, cons))
  out$paths_descriptives <- aw_pls_describe(model$path_coef, out$boot_paths, src, tgt)
  out$loadings_descriptives <- aw_pls_describe(model$outer_loadings, out$boot_loadings)
  out$weights_descriptives <- aw_pls_describe(model$outer_weights, out$boot_weights)
  out$HTMT_descriptives <- aw_pls_describe(htmt, out$boot_HTMT)
  out$total_paths_descriptives <- aw_pls_describe(total, out$boot_total_paths, src, tgt)
  out$boots <- boots
  for (f in names(counts)) out[[f]] <- counts[[f]]
  out$seed <- seed
  out$native <- TRUE
  class(out) <- c("boot_seminr_model", "seminr_model")
  out
}

# The compiled kernel of aw_glm_fit() and aw_boot_glm(): iteratively
# reweighted least squares for three families with their canonical links,
# re-implemented from R's documentation (?glm, ?glm.control, ?family) and
# black-box tests against R's own family functions and glm.fit().
# fam: 1 binomial (logit), 2 poisson (log), 3 gaussian (identity).
# B = 0 fits the data once; B > 0 refits B case resamples. One row per fit:
# the coefficients, the deviance, the iterations, a status and the rank.
# Status 0: converged with nothing to warn about. Anything else (divergence,
# invalid fitted values, no convergence, fitted values at the boundary) is a
# case for which glm.fit() warns or errors; the R wrappers hand it to
# glm.fit() so that the result and its warnings are R's own.
aw_glm_kernel <- function(X, y, wt, fam, B, epsilon, maxit) {
  n <- nrow(X)
  p <- ncol(X)
  eps <- 2.220446049250313e-16
  tol <- min(1e-07, epsilon / 1000)
  reps <- if (B > 0L) B else 1L
  out <- matrix(0, reps, p + 4L)
  Xb <- matrix(0, n, p)
  yb <- numeric(n)
  wb <- numeric(n)
  off <- numeric(n)
  eta <- numeric(n)
  mu <- numeric(n)
  mev <- numeric(n)
  vm <- numeric(n)
  dr <- numeric(n)
  for (b in seq_len(reps)) {
    if (B > 0L) {
      idx <- sample.int(n, replace = TRUE)
    } else {
      idx <- seq_len(n)
    }
    for (i in seq_len(n)) {
      k <- idx[i]
      yb[i] <- y[k]
      wb[i] <- wt[k]
      for (j in seq_len(p)) Xb[i, j] <- X[k, j]
    }
    # Starting values: the family's mustart, through the link.
    for (i in seq_len(n)) {
      if (fam == 1L) {
        m <- (wb[i] * yb[i] + 0.5) / (wb[i] + 1)
        eta[i] <- log(m / (1 - m))
      } else if (fam == 2L) {
        eta[i] <- log(yb[i] + 0.1)
      } else {
        eta[i] <- yb[i]
      }
    }
    start <- numeric(p)
    piv <- seq_len(p)
    rank <- p
    status <- 0
    iter <- 0L
    conv <- FALSE
    devold <- 0
    dev <- 0
    for (it in 0:maxit) {
      # Fitted values, deviance residuals and their validity at eta.
      for (i in seq_len(n)) {
        e <- eta[i]
        if (is.na(e)) {
          status <- 2
          m <- e
        } else if (fam == 1L) {
          ex <- if (e < -30) eps else if (e > 30) 1 / eps else exp(e)
          m <- ex / (1 + ex)
          a1 <- if (yb[i] != 0) yb[i] * log(yb[i] / m) else 0
          y0 <- 1 - yb[i]
          a0 <- if (y0 != 0) y0 * log(y0 / (1 - m)) else 0
          dr[i] <- 2 * wb[i] * (a1 + a0)
          if (!(m > 0 && m < 1)) status <- 2
        } else if (fam == 2L) {
          m <- max(exp(e), eps)
          dr[i] <- if (yb[i] > 0) 2 * (wb[i] * (yb[i] * log(yb[i] / m) - (yb[i] - m))) else 2 * (m * wb[i])
          if (m == 1 / 0) status <- 2
        } else {
          m <- e
          dr[i] <- wb[i] * ((yb[i] - m)^2)
        }
        mu[i] <- m
      }
      if (status != 0) break
      dev <- sum(dr)
      if (it > 0L) {
        if (is.na(dev) || abs(dev) == 1 / 0) {
          status <- 1
          break
        }
        if (abs(dev - devold) / (abs(dev) + 0.1) < epsilon) {
          conv <- TRUE
          break
        }
      }
      devold <- dev
      if (it == maxit) break
      # One weighted least-squares step.
      for (i in seq_len(n)) {
        e <- eta[i]
        if (fam == 1L) {
          mev[i] <- if (e > 30 || e < -30) eps else exp(e) / ((1 + exp(e)) * (1 + exp(e)))
          vm[i] <- mu[i] * (1 - mu[i])
        } else if (fam == 2L) {
          mev[i] <- max(exp(e), eps)
          vm[i] <- mu[i]
        } else {
          mev[i] <- 1
          vm[i] <- 1
        }
        if (is.na(vm[i]) || vm[i] == 0 || is.na(mev[i]) || mev[i] == 0) status <- 3
      }
      if (status != 0) break
      z <- (eta - off) + (yb - mu) / mev
      w <- sqrt((wb * mev^2) / vm)
      fit <- .lm.fit(Xb * w, z * w, tol = tol)
      cf <- fit$coefficients
      for (j in seq_len(p)) if (is.na(cf[j]) || abs(cf[j]) == 1 / 0) status <- 4
      if (status != 0) break
      piv <- fit$pivot
      rank <- fit$rank
      start[piv] <- cf
      eta <- drop(Xb %*% start)
      eta <- eta + off
      iter <- it + 1L
    }
    if (status == 0 && !conv) status <- 5
    if (status == 0) {
      for (i in seq_len(n)) {
        if (fam == 1L && (mu[i] > 1 - 10 * eps || mu[i] < 10 * eps)) status <- 6
        if (fam == 2L && mu[i] < 10 * eps) status <- 6
      }
    }
    if (rank < p) for (j in (rank + 1L):p) start[piv[j]] <- NA_real_
    for (j in seq_len(p)) out[b, j] <- start[j]
    out[b, p + 1L] <- dev
    out[b, p + 2L] <- iter
    out[b, p + 3L] <- status
    out[b, p + 4L] <- rank
  }
  out
}

# The family code the kernel knows, or 0 for glm.fit().
aw_glm_family <- function(family) {
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) stop("'family' must be a family object, such as binomial()")
  code <- c("binomial logit" = 1L, "poisson log" = 2L, "gaussian identity" = 3L)[paste(family$family, family$link)]
  list(family = family, code = if (is.na(code)) 0L else unname(code))
}

# Data the kernel fits exactly as glm.fit() does: finite numbers, positive
# prior weights, and responses for which glm.fit() neither warns nor errors
# before it starts (binomial: 0 <= y <= 1 with whole numbers of successes;
# poisson: whole non-negative counts).
aw_glm_native_ok <- function(X, y, w, code) {
  if (code == 0L || !is.numeric(X) || !(is.numeric(y) || is.logical(y)) || nrow(X) == 0L || ncol(X) == 0L ||
      length(y) != nrow(X) || length(w) != nrow(X)) return(FALSE)
  if (!all(is.finite(X)) || !all(is.finite(y)) || !all(is.finite(w)) || any(w <= 0)) return(FALSE)
  switch(code,
    { m <- w * y; all(y >= 0 & y <= 1) && !any(abs(m - round(m)) > 0.001) },
    all(y >= 0 & y == round(y)),
    TRUE)
}

aw_glm_unit <- function() {
  if (is.null(aw_state$glm)) aw_state$glm <- aw_compile(aw_glm_kernel, name = "aw_glm")
  aw_state$glm
}

#' Fit a generalised linear model, compiled
#'
#' Iteratively reweighted least squares for binomial (logit), Poisson (log)
#' and Gaussian (identity) models, run as native code through Arwad and
#' identical to [stats::glm.fit()]: the same coefficients, deviance and
#' number of iterations, bit for bit. Each step solves its weighted least
#' squares problem with R's own LINPACK QR, as `glm.fit()` does. Anything the
#' native path does not cover is handed to `glm.fit()` itself: other families
#' and links, zero or missing weights, non-finite data, and every fit for
#' which `glm.fit()` warns (no convergence, fitted probabilities of 0 or 1,
#' step halving), so that the result and its warnings are R's.
#'
#' @param x A numeric design matrix; include a column of ones for an intercept.
#' @param y The response: for the binomial family, proportions between 0 and
#'   1 (with `weights` the numbers of trials) or 0/1 outcomes.
#' @param family A family object or function: `binomial()`, `poisson()` or
#'   `gaussian()` run natively with their default links.
#' @param weights Optional prior weights.
#' @param epsilon,maxit Convergence tolerance and iteration limit, as in
#'   [stats::glm.control()].
#' @return A list: `coefficients` (named by the columns of `x`; `NA` for
#'   columns dropped as linearly dependent), `deviance`, `iter`, `converged`,
#'   `rank`, and `native` (`FALSE` when `glm.fit()` produced the result).
#' @export
#' @examples
#' X <- cbind(intercept = 1, wt = mtcars$wt, hp = mtcars$hp)
#' fit <- aw_glm_fit(X, mtcars$am, family = binomial())
#' fit$coefficients
#' identical(fit$coefficients, glm.fit(X, mtcars$am, family = binomial())$coefficients)
aw_glm_fit <- function(x, y, family = stats::gaussian(), weights = NULL, epsilon = 1e-8, maxit = 25L) {
  fam <- aw_glm_family(family)
  ctl <- stats::glm.control(epsilon = epsilon, maxit = maxit)
  X <- as.matrix(x)
  w <- if (is.null(weights)) rep(1, NROW(y)) else weights
  r_fit <- function() {
    f <- stats::glm.fit(X, y, weights = weights, family = fam$family, control = ctl)
    list(coefficients = f$coefficients, deviance = f$deviance, iter = f$iter, converged = f$converged,
         rank = f$rank, native = FALSE)
  }
  if (!aw_glm_native_ok(X, y, w, fam$code)) return(r_fit())
  p <- ncol(X)
  res <- tryCatch(aw_glm_unit()(matrix(as.double(X), nrow(X), p), as.double(y), as.double(w), fam$code, 0L,
                                ctl$epsilon, as.integer(ctl$maxit)),
                  error = function(e) NULL)
  if (is.null(res) || res[1L, p + 3L] != 0) return(r_fit())
  cf <- res[1L, seq_len(p)]
  names(cf) <- colnames(X)
  list(coefficients = cf, deviance = res[1L, p + 1L], iter = as.integer(res[1L, p + 2L]), converged = TRUE,
       rank = as.integer(res[1L, p + 4L]), native = TRUE)
}

#' Bootstrap GLM coefficients, compiled
#'
#' Case-resampling (nonparametric) bootstrap of the coefficients of a
#' binomial (logit), Poisson (log) or Gaussian (identity) model. Each
#' replicate draws `nrow(x)` rows with replacement
#' (`sample.int(n, replace = TRUE)`) and fits them as [aw_glm_fit()] does. For
#' the same seed the result is identical to this loop in R, warnings
#' included: `for (b in 1:B) { i <- sample.int(n, replace = TRUE);
#' out[b, ] <- glm.fit(x[i, , drop = FALSE], y[i], weights = weights[i],
#' family = family)$coefficients }`. Replicates for which `glm.fit()` warns
#' are refitted by `glm.fit()` itself, from the same resample.
#'
#' @inheritParams aw_glm_fit
#' @param B The number of bootstrap replicates.
#' @return A `B` by `ncol(x)` matrix of coefficients, one row per replicate,
#'   with the column names of `x`.
#' @export
#' @examples
#' X <- cbind(intercept = 1, wt = mtcars$wt)
#' set.seed(1)
#' b <- aw_boot_glm(X, mtcars$vs, family = binomial(), B = 200)
#' apply(b, 2, sd) # bootstrap standard errors
aw_boot_glm <- function(x, y, family = stats::gaussian(), B = 1000L, weights = NULL, epsilon = 1e-8, maxit = 25L) {
  fam <- aw_glm_family(family)
  ctl <- stats::glm.control(epsilon = epsilon, maxit = maxit)
  X <- as.matrix(x)
  n <- nrow(X)
  p <- ncol(X)
  B <- as.integer(B)
  w <- if (is.null(weights)) rep(1, NROW(y)) else weights
  one <- function(i) {
    stats::glm.fit(X[i, , drop = FALSE], y[i], weights = if (!is.null(weights)) weights[i],
                   family = fam$family, control = ctl)$coefficients
  }
  r_loop <- function() {
    out <- matrix(NA_real_, B, p)
    for (b in seq_len(B)) out[b, ] <- one(sample.int(n, replace = TRUE))
    out
  }
  if (!aw_glm_native_ok(X, y, w, fam$code)) {
    out <- r_loop()
  } else {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(0L)
    seed0 <- get(".Random.seed", envir = globalenv())
    res <- tryCatch(aw_glm_unit()(matrix(as.double(X), n, p), as.double(y), as.double(w), fam$code, B,
                                  ctl$epsilon, as.integer(ctl$maxit)),
                    error = function(e) NULL)
    if (is.null(res)) {
      assign(".Random.seed", seed0, envir = globalenv())
      out <- r_loop()
    } else {
      out <- res[, seq_len(p), drop = FALSE]
      redo <- which(res[, p + 3L] != 0)
      if (length(redo)) {
        # Replay the resamples up to the last one to redo, then leave the
        # stream where the kernel left it.
        seed1 <- get(".Random.seed", envir = globalenv())
        assign(".Random.seed", seed0, envir = globalenv())
        for (b in seq_len(max(redo))) {
          i <- sample.int(n, replace = TRUE)
          if (b %in% redo) out[b, ] <- one(i)
        }
        assign(".Random.seed", seed1, envir = globalenv())
      }
    }
  }
  colnames(out) <- colnames(X)
  out
}

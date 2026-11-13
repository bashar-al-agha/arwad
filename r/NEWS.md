# arwad 0.1.0.9007

* New `aw_assess_cvpat()` and `aw_assess_cvpat_compare()`: the
  cross-validated predictive ability test (CVPAT; Liengaard et al., 2021;
  Sharma et al., 2023), compiled replacements for
  `seminrExtras::assess_cvpat()` and `seminrExtras::assess_cvpat_compare()`
  with their arguments and defaults. The cross-validated predictions come from
  the compiled PLSpredict of `aw_predict_pls()` and the bootstrap of the loss
  differences runs as native code. The tables, `Boot T value` and
  `Boot P Value` included, and the random-number state after the call are
  identical to `seminrExtras` 1.0.4 on the models of composites in the test
  set. As in `seminrExtras`, the indicator-average benchmark predicts with the
  item means of the whole sample. `testtype = "less"`, on which
  `seminrExtras` 1.0.4 fails, gives the lower-tail p-value. Models outside the
  native subset go to `seminrExtras` when it is installed (ADR-0008, amendment
  of 13 November 2026).
* New `aw_q2_predict()`, and `q2_predict` in the result of `aw_predict_pls()`:
  Q²predict of each predicted indicator, one minus the PLS model's sum of
  squared out-of-sample errors over that of the training-fold item means
  (Shmueli et al., 2016, 2019).

# arwad 0.1.0.9006

* New `aw_predict_pls()`: PLSpredict (k-fold cross-validated prediction), a
  compiled replacement for `seminr::predict_pls()` with its arguments and
  defaults. It draws the same folds from the same random-number stream,
  estimates every training set with the compiled PLS-SEM kernel (on several
  cores if asked), and returns `predict_pls()`'s object bit for bit, so
  `summary()` reports the same RMSE and MAE for the PLS model and the
  linear-model benchmark. Identical to `seminr` 2.6.0 on all 79 models of
  composites in the test set, with earliest and direct antecedents, 5 to 10
  folds, repetitions and leave-one-out, under R's reference BLAS and under
  OpenBLAS (as R uses it on Linux), whose rounding depends on the shape of
  each matrix product: Arwad forms the products in the shapes `seminr`
  uses. On the Primer's extended corporate reputation model, 10 folds and
  10 repetitions take 0.67 s on one core against 7.70 s for `seminr`
  (11.5 times faster), leave-one-out 2.73 s against 12.49 s on four cores
  (4.6 times; `bench/phase2/REPORT.md`, section 3e). Models of composites (Mode A and
  Mode B) are native with `seminr` 2.6.0 or later; models with reflective
  constructs (which `seminr` 2.6.0 predicts with a different algorithm),
  interaction terms, and every call under earlier `seminr` versions go to
  `seminr` (ADR-0008).

# arwad 0.1.0.9005

* Licence: from this version on, Arwad is released under the Apache License
  2.0 (ADR-0003, amendment of 2 November 2026). The licence adds a patent
  licence from every contributor, states that it grants no right to the
  name "Arwad", and asks redistributors to keep the `NOTICE` file. Every
  earlier version, the v0.1.0 release included, remains available under the
  MIT licence.

# arwad 0.1.0.9004

* `aw_boot_pls()`: the first stage of a model with two-stage interaction
  terms now stops after its weights, which are all the second stage uses.
  The loadings, paths, total effects and HTMT it no longer computes were
  never used, so the results are unchanged, bit for bit. On the Primer's
  moderation model the first stage takes 18% less time and a resample of
  both stages 9% less (one core, 1,000 resamples, the same compiled kernel
  with and without the shortcut; ADR-0006, amendment "two-stage interaction
  terms").
* Benchmark (bench/phase2/REPORT.md, §3d, now three rounds of 10,000
  resamples, seminr's fastest against Arwad's median): 57.6 s on four cores
  against 390.1 s for seminr 2.6.0 on the same cores, 6.8 times faster, and
  3.2 times from one core, with identical arrays. The 4.0 times given for
  0.1.0.9003 came from a single round in which the four workers were
  unusually slow (0.78 of the one-core time, against 0.47 in the three
  rounds); the shortcut above saves about 9% of a resample, and the rest of
  the difference is that round's variation.

# arwad 0.1.0.9003

* `aw_boot_pls()` runs models with two-stage interaction terms (moderation:
  `interaction_term(iv, moderator, method = two_stage)`, seminr's default,
  and `quadratic_term()`) natively, in models of Mode A and Mode B
  composites. As seminr does, every resample is estimated twice: the model
  without its interaction terms first (always with `maxIt = 300` and
  `stopCriterion = 7`), then the full model, in which each term's indicator
  is the product of the first stage's scores and the term's score keeps
  that product's standard deviation. On 34 test models (the Primer's
  corporate reputation model with switching costs moderating CUSA -> CUSL,
  ECSI variants, a quadratic term and 24 synthetic models) the estimates and
  bootstraps are bit-identical to seminr 2.6.0 and 2.5.0 (HTMT within four
  units in the last place). BCa intervals work for these models too.
  Product-indicator and orthogonal terms, and consistent PLS with
  interaction terms, are still handed to seminr (ADR-0006, amendment
  "two-stage interaction terms").
* Benchmark (bench/phase2/REPORT.md, §3d): 10,000 resamples of the Primer's
  corporate reputation model with switching costs moderating CUSA -> CUSL
  took 97.7 s on four cores against 386.9 s for seminr 2.6.0 on the same
  cores, 4.0 times faster (one round; 3.1 times from one core), with
  identical arrays. The gain is smaller than for models without an
  interaction term because every resample is estimated twice; the first
  stage still computes loadings, paths and HTMT that only the second needs,
  the next optimisation. Correction: the benchmark put one of its 45,000
  HTMT values four units in the last place from seminr's, so the entry
  above now says four units, the bound the tests check (it said two, from
  the 3-resample tests).

# arwad 0.1.0.9002

* `aw_pls_ci(type = "bca")`: bias-corrected and accelerated (BCa) intervals
  (Efron, 1987), which the PLS-SEM Primer names for bootstrap distributions
  with skewness beyond plus or minus 2. The acceleration comes from the
  jackknife, the model re-estimated on its n leave-one-out samples in one
  call of the compiled PLS kernel, which now takes the sample size from its
  index matrix (bootstrap results unchanged). The leave-one-out estimates
  equal `seminr::estimate_pls()` bit for bit (Mode A, Mode B, consistent
  PLS), and the limits equal Efron's formula; BCa intervals for the ECSI
  model took 0.43 s, against 4.7 s for the same jackknife with seminr.
  Models outside the native subset stop with the reason (ADR-0007,
  amendment).

# arwad 0.1.0.9001

* `aw_boot_pls()` follows the installed seminr, 2.5.0 or 2.6.0 (CRAN, 13
  October 2026), whose consistent PLS changed in two ways; both were matched
  by black-box tests (ADR-0006, amendment "seminr 2.6.0"). seminr 2.6.0
  treats composites as fully reliable (rho_A = 1), as Dijkstra and Henseler
  do, so path coefficients of models that mix `reflective()` and Mode A
  `composite()` constructs differ from 2.5.0's; and `bootstrap_model()`
  leaves out every resample whose solution is inadmissible (a rho_A above
  one, a reflective loading above one, a corrected construct correlation of
  one or more, an R-squared outside [0, 1], or a construct correlation
  matrix implied by the structural model whose smallest eigenvalue is below
  1e-8), and records `boots_requested` and `boots_dropped`. Under 2.6.0,
  `aw_boot_pls()` does the same, with the same counts. The CI now checks the
  package against both versions. Found when the CI of the consistent-PLS
  commit failed on all systems: CI installs the current CRAN release, this
  laptop still had 2.5.0.
* When `aw_boot_pls()` hands a model to seminr, it passes `nboot` as given,
  so the counts in the result keep its type, as when seminr is called
  directly (it passed an integer before), and it passes `cores`, as its
  documentation says (before, seminr chose the number of processes itself:
  all cores, which seminr 2.5.0 does even under `R CMD check`'s limit of
  two, and the new CI job with 2.5.0 failed there).
* Benchmark of consistent PLS (bench/phase2/REPORT.md, section 3c): the ECSI
  model with all seven constructs reflective, 10,000 resamples, took 20.0 s
  on four cores against 169.8 s for `seminr::bootstrap_model()` 2.6.0 on the
  same cores (8.5 times faster) and 26.1 s on one core (6.5 times faster);
  arrays identical, the same 3,679 resamples kept.

# arwad 0.1.0.9000

* `aw_boot_pls()` runs models with `reflective()` constructs as native code:
  Mode A weights followed by consistent PLS (PLSc; Dijkstra and Henseler,
  2015), the way seminr 2.5.0 applies it. Estimates and bootstrap arrays are
  bit-identical to seminr on the 46 test models with reflective constructs
  (all reflective, next to Mode A composites, next to Mode B constructs,
  single-indicator reflective constructs); HTMT stays within four units in
  the last place. Two seminr rules found by black-box tests are kept, so that
  the bootstrap belongs to the estimate seminr produced: in a model with a
  reflective construct the correction applies to every Mode A construct with
  two or more indicators, `composite()` ones included (the published method
  corrects only constructs modelled as common factors); and resamples for
  which consistent PLS is undefined (a negative rho_A) are left out, with
  `boots` counting the others and a message saying how many. Still handed to
  seminr: the factorial and centroid schemes, interaction terms,
  higher-order constructs and missing data.
* Fixed: the resamples are estimated with seminr's default `maxIt` (300) and
  `stopCriterion` (7), as `bootstrap_model()` does whatever the model was
  estimated with. Until now `aw_boot_pls()` used the model's own settings,
  which gave numbers different from seminr's for a model estimated with
  another `stopCriterion`, or for a resample needing more than 300
  iterations. The published comparisons and benchmarks met neither case
  (their arrays were identical to seminr's).
* With fewer than two usable resamples (for example `nboot = 1`),
  `aw_boot_pls()` hands the model to `bootstrap_model()`, which stops, as it
  does when called directly; before, it returned a one-resample bootstrap.

# arwad 0.1.0

First citable release (28 October 2026), archived on Zenodo:
<https://doi.org/10.5281/zenodo.23725078>.

* `aw_compile()` compiles R functions written in a statistical subset
  (scalars, vectors and matrices, loops, R's random-number streams, BLAS and
  LAPACK through R's own libraries) to native code through Goose; whatever
  the native code does not model is re-run in GNU R. Bit-identical to GNU R
  across a 201-function corpus on Linux, macOS and Windows, and clean under
  AddressSanitizer.
* Statistical engines, each validated on at least 50 datasets:
  `aw_boot_ols()`, `aw_glm_fit()` and `aw_boot_glm()` (identical to
  `.lm.fit()` and `glm.fit()`), `aw_boot_pls()` (a drop-in for
  `seminr::bootstrap_model()` with Mode A and Mode B composites) and
  `aw_pls_ci()` (bias-corrected percentile confidence intervals).
* The development entries below record how the release came about.

# arwad 0.0.0.9015

* New `aw_pls_ci()`: bias-corrected percentile confidence intervals, the
  intervals the PLS-SEM Primer asks for, for every coefficient of a
  bootstrapped model (paths, total and total indirect effects, loadings,
  weights, HTMT) and for specific indirect effects, from `aw_boot_pls()` or
  `seminr::bootstrap_model()`. With `type = "percentile"` the limits are
  seminr's `summary()` limits, bit for bit; the bias correction is
  `boot::boot.ci()`'s with zero acceleration. Each table adds the bias and
  the skewness of the bootstrap values. Unlike seminr's `summary()`, the
  total indirect effects follow `alpha` too.
* Fixed: with `cores` above 1 on Windows, `aw_boot_pls()`'s worker processes
  ran the engine in GNU R (same results, about 9 times slower) when Goose
  was set only by `options(arwad.goose = )`. Workers now receive every
  `arwad` option and the Goose executable this session uses. The published
  benchmarks set Goose by environment variable and were not affected.

# arwad 0.0.0.9014

* `aw_boot_pls()` runs models with formative constructs
  (`composite(..., weights = mode_B)`) as native code, including the Primer's
  extended corporate reputation model; results are bit-identical to seminr
  (HTMT within four units in the last place, see below). Only `reflective()`
  (PLSc) models, other weighting schemes, interactions, higher-order
  constructs and missing data still go to seminr.
* `cov(X)` and `cor(X)` of one matrix compile, matched to GNU R bit for bit
  (R's `cor(X)` rounds differently from `cor(X, X)`).
* Correction: HTMT agrees with seminr to within four units in the last
  place, not two as stated since 0.0.0.9010. Measured over 1,000 resamples:
  the ECSI model has values up to 3 units apart, the Primer's extended
  corporate reputation model one value 4 units apart; about three in four
  values are identical. The test now checks the distance in units.

# arwad 0.0.0.9013

* `cov(m[, j], k[, i])` and `cor(m[, j], k[, i])` of matrix columns read the
  columns in place and keep each column's mean (and standard deviation) until
  the matrix is written again, so repeated calls skip repeated work; values
  stay R's, bit for bit. The PLS-SEM engine uses it: 1,000 resamples of the
  ECSI model in 2.06 s instead of 2.86 s on the benchmark laptop, results
  unchanged.
* The PLS-SEM kernel builds each regression's transposed design once.

# arwad 0.0.0.9012

* `aw_boot_pls()` is 1.8 times faster (loadings and HTMT from one covariance
  and one correlation matrix, compiled matrix forms of `cov()` and `cor()`),
  and `cores =` runs resamples in parallel processes with identical results.

# arwad 0.0.0.9011

* `m[, j] <- v` and `m[i, ] <- v` compile. A function or signature that runs
  in GNU R instead of native code now says so once, with the reason
  (`options(arwad.quiet = TRUE)` turns this off).

# arwad 0.0.0.9010

* Third Phase 2 engine: `aw_boot_pls()`, a compiled drop-in for
  `seminr::bootstrap_model()` (Mode A composites, path weighting scheme).
  Path coefficients, total effects, loadings and weights are bit-identical to
  seminr on 54 models; HTMT agrees to within two units in the last place.
* `colMeans()`, `cov()` and `cor()` compile, matched to GNU R bit for bit; on
  64-bit ARM, `cov()` and `cor()` call R's own routine.

# arwad 0.0.0.9009

* Second Phase 2 engine: `aw_glm_fit()` and `aw_boot_glm()`, iteratively
  reweighted least squares for binomial (logit), Poisson (log) and Gaussian
  (identity) models, compiled. Coefficients, deviance and iteration counts
  are bit-identical to `glm.fit()` on 55 designs for each family; every fit
  for which `glm.fit()` warns or errors is handed to `glm.fit()`, so warnings
  are R's too. Re-implemented from R's documentation and black-box tests.
* `fit <- .lm.fit(x, y, tol)` compiles, with `fit$coefficients`,
  `$residuals`, `$rank` and `$pivot`; rank-deficient designs now run
  natively, pivoted exactly as R pivots them.
* `drop()` compiles; so does `x[i] <- v` with an index vector and a value
  vector (zero indices dropped, value recycled, as in R).
* Fix: `x ^ 2` is computed as `x * x`, as R's arithmetic does; R's C
  function `R_pow()` can round squares differently on 64-bit Windows.
* Fix: `sqrt()` is always the processor's correctly rounded square root
  (units are built with `-fno-math-errno`); GCC could otherwise call the C
  library's `sqrt()`, which on Windows rounds some results differently.
* Fix: the portable build path (R CMD SHLIB) works on macOS again; R's
  `R_ext/Applic.h` is no longer included.

# arwad 0.0.0.9008

* First Phase 2 engine: `aw_boot_ols()`, a compiled case-resampling
  bootstrap of OLS coefficients, bit-identical to R running the same
  bootstrap; 43.7x faster than an `lm()` loop and 1.9x faster than an
  `.lm.fit()` loop on one laptop (bench/phase2/REPORT.md).
* `.lm.fit(x, y)$coefficients` and `$residuals` compile, through R's own
  LINPACK `dqrls`; validated bit for bit on 55 designs.
* `stats` added to Imports.

# arwad 0.0.0.9007

* Corpus grown to 201 functions (Phase 1 target: 200): `slice4.R` and
  `slice5.R` add statistics, estimators, numerical methods, finance and
  time series, linear algebra (IRLS logistic regression, Jacobi, LU), rank
  statistics, integer algorithms and simulation (bootstrap, permutation
  tests, Gibbs, Metropolis). All run natively and match GNU R bit for bit.

# arwad 0.0.0.9006

* `sample()` and `sample.int()` with and without replacement, drawing
  indices through R's own `R_unif_index()`; bootstrap resampling such as
  `mean(x[sample.int(n, replace = TRUE)])` now compiles and reproduces R's
  stream. Corpus: 92 functions.

# arwad 0.0.0.9005

* Compiler slice 3, random numbers: `runif()` and `rnorm()` draw from GNU
  R's own generators through its C API, between `GetRNGstate()` and
  `PutRNGstate()`, so values and the final `.Random.seed` are identical to
  R's. A call that bails out restores `.Random.seed` before R re-runs it.
* `floor()` and `ceiling()`.
* Corpus: 86 functions; seeded entries compare the generator state too.

# arwad 0.0.0.9004

* Compiler slice 2b, part 2: `%*%`, `crossprod()` and `solve()` through R's
  own BLAS (`dgemv`, `dgemm`, `dsyrk`) and LAPACK (`dgesv`, `dgecon`), so
  OLS as `solve(crossprod(X), crossprod(X, y))` compiles and matches R bit
  for bit. Units now link `$(LAPACK_LIBS) $(BLAS_LIBS)`.
* Compiled `+` and `*` select the result NaN in C as each platform's
  hardware does for R's operand order (arm64 favours signalling NaNs); the
  earlier `x + 0.0` quieting was optimised away by clang.
* Corpus: 79 functions.

# arwad 0.0.0.9003

* Compiler slice 2b, part 1, matrices: matrix arguments and locals,
  `matrix()`, `m[i, j]` read and write, rows and columns, `t()`, `nrow()`,
  `ncol()`, `dim()`, and elementwise arithmetic that keeps dimensions.
  Results are returned with their `dim` attribute.
* `mean()` now follows GNU R's algorithm, including its fallback to the sum
  of `x[i] / n` when the plain sum is not finite (found by CI on macOS
  arm64, where `long double` is `double`).
* Compiled `+` and `*` return a quieted `NA` or `NaN`, as R's arithmetic
  does; arm64 hardware otherwise favours R's signalling `NA_real_`.
* Corpus: 72 functions.

# arwad 0.0.0.9002

* Compiler slice 2a, vectors: vectorised arithmetic, comparison, logic and
  math with R's recycling; `c()`, `mean()`, `min()`, `max()`; `a:b`,
  `seq_len()`, `seq_along()` as values; index vectors and logical masks for
  reading and for assigning a single value. Corpus: 59 functions, all native
  and bit-identical to GNU R.

# arwad 0.0.0.9001

* First compiler slice (Phase 1): `aw_compile()` for double, integer and
  logical scalars and vectors, `if`/`for`/`while`/`repeat`, indexing and
  sub-assignment, and `length`, `numeric`, `integer`, `logical`, `sum`, `abs`,
  `sqrt`, `exp`, `log`, `is.na`.
* Guarded native execution: cases the native code does not model re-run the
  call in GNU R, so values, errors and warnings are R's (ADR-0005).
* `aw_explain()`, `aw_source()`, `aw_toolchain()`, `aw_cache_clear()`.
* Equivalence corpus of 39 functions, all native and bit-identical to GNU R.
* Compiled `+` and `*` keep the first operand's `NA` or `NaN`, as GNU R does,
  whatever operand order the C compiler picks (found by CI on Linux).
* Direct build path: R's compile and link commands are read once from a dry
  run of `R CMD SHLIB` and the Goose runtime is compiled once per toolchain,
  so a new signature builds in about 1 s instead of 4-5 s.
  `options(arwad.build = "shlib")` keeps the old path, which is also the
  fallback.
* Dispatch: the native unit of the previous call checks the argument types
  in C, cutting the per-call cost from about 33 us to 4.5 us.

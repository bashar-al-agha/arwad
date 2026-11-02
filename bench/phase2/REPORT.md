# Phase 2 engines: compiled bootstraps of OLS, GLM and PLS-SEM models

**Date:** 27 October 2026 (PLS-SEM figures re-measured 28 October 2026, version 0.0.0.9013; formative constructs, §3b, version 0.0.0.9014; reflective constructs, §3c, 29 October 2026, version 0.1.0.9001 against seminr 2.6.0; moderation, §3d, 2 November 2026, version 0.1.0.9003 against seminr 2.6.0) · **Status:** preliminary, single laptop (the Phase 0 machine) · **Raw data:** [`ols_boot_timing.csv`](ols_boot_timing.csv), [`glm_boot_timing.csv`](glm_boot_timing.csv), [`pls_boot_timing.csv`](pls_boot_timing.csv), [`pls_modeB_timing.csv`](pls_modeB_timing.csv), [`pls_plsc_timing.csv`](pls_plsc_timing.csv), [`pls_moderation_timing.csv`](pls_moderation_timing.csv) · **Scripts:** [`ols_boot_bench.R`](ols_boot_bench.R), [`glm_boot_bench.R`](glm_boot_bench.R), [`pls_boot_bench.R`](pls_boot_bench.R), [`pls_modeB_bench.R`](pls_modeB_bench.R), [`pls_plsc_bench.R`](pls_plsc_bench.R), [`pls_moderation_bench.R`](pls_moderation_bench.R)

## Summary

`arwad::aw_boot_ols()` runs a case-resampling bootstrap of least-squares coefficients as native code. Each replicate draws rows with R's own `sample.int()` stream and fits them with R's own LINPACK QR (`dqrls`, the routine behind `.lm.fit()` and `lm()`), so the result is **bit-identical** to R running the same bootstrap with the same seed. For 2,000 replicates of a 200 by 5 design it took **60.5 ms**: **43.7 times faster** than the `lm()` loop most researchers write, and **1.9 times faster** than the fastest hand-optimised R loop.

`arwad::aw_glm_fit()` and `arwad::aw_boot_glm()` do the same for generalised linear models: iteratively reweighted least squares for binomial (logit), Poisson (log) and Gaussian (identity) models. For 1,000 replicates of a logistic regression on a 200 by 5 design, `aw_boot_glm()` took **382.4 ms**: **10.8 times faster** than a `glm()` loop and **3.25 times faster** than a `glm.fit()` loop, with **identical coefficients** to both.

`arwad::aw_boot_pls()` is a compiled replacement for `seminr::bootstrap_model()`, the reference implementation of the PLS-SEM bootstrap in R. It re-estimates a model on case resamples with the PLS path algorithm (Mode A and Mode B composites, and reflective constructs by consistent PLS; path weighting scheme) and returns the same `boot_seminr_model` object, so `summary()` works on it unchanged. Path coefficients, total effects, loadings and weights are **bit-identical** to seminr's; HTMT ratios agree to within four units in the last place (about three in four values identical). For **10,000 resamples** of the seven-construct ECSI model, the number the Primer recommends, it took **12.5 s on four cores** against **241.4 s** for seminr on the same four cores (seminr's fastest of three runs): **19.3 times faster**; on a single core it took **27.3 s** (**8.8 times faster** than seminr on four). For 1,000 resamples on one core it took **2.46 s**: 10.2 times faster than seminr on four cores and 18.7 times faster than seminr on one core. With every construct reflective (consistent PLS, §3c), 10,000 resamples of the same model took **20.0 s on four cores** against **169.8 s** for seminr 2.6.0: **8.5 times faster**, with identical arrays and the same resamples left out. With a two-stage interaction term (moderation, §3d), 10,000 resamples of the Primer's moderation model took **97.7 s on four cores** against **386.9 s** for seminr 2.6.0: **4.0 times faster**, with identical arrays; every resample is estimated twice, so the gain is smaller.

## 1. OLS bootstrap: what was measured

Design: n = 200 rows, k = 5 columns (intercept and four standard-normal predictors), B = 2,000 replicates. Every contender starts from `set.seed(1)`; times are medians of `bench::mark()` (3 to 10 iterations), garbage collections included. R code is byte-compiled (`compiler::cmpfun`).

| Contender | What it is | Median (ms) | `aw_boot_ols()` is |
|---|---|---|---|
| `aw_boot_ols()` | the compiled engine | 60.5 | — |
| `.lm.fit` loop | `.lm.fit(X[idx, , drop = FALSE], y[idx])$coefficients` per replicate: the fastest idiomatic R | 115.5 | 1.9× faster |
| same algorithm in R | the engine's own kernel run by GNU R | 449.9 | 7.4× faster |
| `lm()` loop | `coef(lm(y ~ ., data = d[idx, ]))` per replicate: what most users write | 2,645.0 | 43.7× faster |

Equivalence: the engine's 2,000 by 5 coefficient matrix is identical (`identical(..., num.eq = FALSE)`) to both the same algorithm in R and the `.lm.fit` loop. The first call, in a fresh cache, took 7.0 s (one-off setup plus building the unit); later calls and sessions reuse the native unit.

## 2. GLM bootstrap: what was measured

Design: n = 200 rows, k = 5 columns (intercept and four standard-normal predictors), a 0/1 outcome drawn from a logistic model (coefficients −0.3, 0.8, −0.5, 0.4, 0.2), B = 1,000 replicates. Same protocol as §1.

| Contender | What it is | Median (ms) | `aw_boot_glm()` is |
|---|---|---|---|
| `aw_boot_glm()` | the compiled engine | 382.4 | — |
| `glm.fit` loop | `glm.fit(X[idx, , drop = FALSE], y[idx], family = binomial())$coefficients` per replicate: the fastest idiomatic R | 1,241.0 | 3.25× faster |
| `glm()` loop | `coef(glm(y ~ ., data = d[idx, ], family = binomial))` per replicate: what most users write | 4,123.0 | 10.8× faster |
| same algorithm in R | the engine's own kernel run by GNU R | 4,695.0 | 12.3× faster |

Equivalence: the engine's 1,000 by 5 coefficient matrix is identical to all three R contenders; no replicate produced a warning. The first call, in a fresh cache, took 12.9 s.

## 3. PLS-SEM bootstrap: what was measured

Design: the European Customer Satisfaction Index model (Tenenhaus et al., 2005) on seminr's `mobi` data: 7 constructs, 24 indicators, 12 paths, n = 250; path weighting scheme, Mode A composites, seminr's defaults (`maxIt` = 300, stop criterion 1e-7); seed 123. seminr distributes resamples over a cluster of R processes, one per logical core by default (4 on this laptop: 2 cores, 4 threads); `aw_boot_pls()` runs in one process unless given `cores`. The 1,000-resample times are medians of `bench::mark()` (3 runs), the 10,000-resample times come from three rounds (below).

| Contender | Resamples | Time | `aw_boot_pls()` on 1 core is |
|---|---|---|---|
| `aw_boot_pls()` | 1,000 | 2,457 ms | — |
| `aw_boot_pls(cores = 4)` | 1,000 | 2,296 ms | 0.93 of its time |
| `bootstrap_model()`, all cores | 1,000 | 25,160 ms | 10.24× faster |
| `bootstrap_model(cores = 1)` | 1,000 | 45,870 ms | 18.67× faster |
| same algorithm in R | 1,000 | 37,210 ms | 15.15× faster |

Ten thousand resamples, three rounds of the same comparison on the same laptop (identical results in every round):

| Round | `bootstrap_model()`, 4 cores | `aw_boot_pls()`, 1 core | `aw_boot_pls(cores = 4)` |
|---|---|---|---|
| 1 | 309.7 s | 27.3 s | 12.5 s |
| 2 | 241.5 s | 27.7 s | 12.2 s |
| 3 | 241.4 s | 24.5 s | 13.3 s |
| median | 241.5 s | 27.3 s | 12.5 s |

seminr's time varies on this laptop (212 to 655 s over eight runs on 27 and 28 October; 241 to 310 s in this one), Arwad's much less. The headline ratios are therefore conservative: seminr's **fastest** run against Arwad's **median**: on the same four cores 241.4 / 12.5 = **19.3 times faster**, from a single core 241.4 / 27.3 = **8.8 times faster**. seminr ran evenly this time, so its median gives the same ratios.

Equivalence: the arrays of path coefficients, total effects, loadings and weights are identical (`identical(..., num.eq = FALSE)`) to seminr's with all cores and with one core, at 1,000 and at 10,000 resamples, and for any number of Arwad workers. Of the HTMT values, 75.1% are identical and the largest relative difference is 4.4e-16. The first call, in a fresh cache, took 27.7 s (one-off setup, building the engine, and the 1,000 resamples; 17.2 s in the previous run).

History: the first compiled version (27 October) took 4,584 ms for 1,000 resamples and 41.3 s for 10,000 on one core. Profiling showed half of each resample in correlations computed pair by pair (HTMT) and in cross-loadings, each recomputing column means; computing them as one covariance matrix and one correlation matrix (compiled `cov(X, Y)` and `cor(X, Y)`) cut the engine's time per resample from 5.15 to 2.85 ms, with identical results. That version (0.0.0.9012) measured 3,483 ms for 1,000 resamples (44,570 ms for seminr on one core: 12.80× faster) and, for 10,000, 15.3 s on four cores and 32.8 s on one core against seminr's fastest run of 239.9 s: 15.7 times faster and 7.3 times faster. A second profile then put 80% of a resample in the PLS iterations, where every `cov()` and `cor()` of two columns recomputed both means although the indicators do not change between iterations; the compiler now keeps each column's mean and standard deviation until the matrix is written again ([ADR-0006](../../docs/adr/0006-pls-sem-engine.md), amendment), which took 1,000 resamples from 2.86 to 2.06 s in a direct comparison of the two versions, with identical results, and gives the table above. Parallel workers add little at 1,000 resamples (starting them costs about as much as they save) and halve the time at 10,000 on this two-core laptop.

Other systems: every CI run times 100 resamples of the same model on each system ([`pls_quick.R`](pls_quick.R), posted to the run's summary). Version 0.0.0.9013, on GitHub's runners (the first results, with version 0.0.0.9011, were 0.36, 0.38 and 0.28 s):

| System | `aw_boot_pls()` | `bootstrap_model(cores = 1)` | Same algorithm in R |
|---|---|---|---|
| Linux x86-64 | 0.16 s | 1.39 s (8.9× slower) | 2.07 s (13.2× slower) |
| macOS arm64 | 0.21 s | 0.76 s (3.6× slower) | 1.39 s (6.6× slower) |
| Windows x86-64 | 0.21 s | 1.63 s (7.8× slower) | 2.10 s (10.0× slower) |

Every run was native (no fallback to GNU R) and identical to seminr. Three cautions: these machines run seminr faster than the laptop of the table above, so against seminr the gain there is 3.6 to 8.9 times; 100 resamples is a short run, in which fixed costs weigh more than at 1,000 or 10,000; and on 64-bit ARM the compiled `cov()` and `cor()` call R's own routine for every pair of columns (see the spec), so the kept column statistics save nothing there, which is why macOS gains least. Against the same algorithm in R, which measures what compilation itself adds, the gain is 6.6 to 13.2 times.

### 3b. Formative constructs: the Primer's extended corporate reputation model

Design: the extended corporate reputation model of Hair, Hult, Ringle and Sarstedt's *Primer* on seminr's `corp_rep_data`: 8 constructs, 31 indicators, 13 paths, n = 344; the four drivers QUAL, PERF, CSOR and ATTR are formative (Mode B), the other constructs Mode A composites (CUSA a single item); path weighting scheme, `maxIt` = 3,000 as the Primer recommends, stop criterion 1e-7, seed 123. Version 0.0.0.9014, the same laptop, 28 October 2026 ([`pls_modeB_bench.R`](pls_modeB_bench.R), [`pls_modeB_timing.csv`](pls_modeB_timing.csv)). The 1,000-resample times are medians of three runs; the 10,000-resample times come from one round, seminr first.

| Contender | Resamples | Time | `aw_boot_pls()` on 1 core is |
|---|---|---|---|
| `aw_boot_pls()` | 1,000 | 5,012 ms | — |
| `aw_boot_pls(cores = 4)` | 1,000 | 4,210 ms | 0.84 of its time |
| `bootstrap_model()`, all cores | 1,000 | 38,810 ms | 7.74× faster |
| `bootstrap_model(cores = 1)` | 1,000 | 77,240 ms | 15.41× faster |
| `aw_boot_pls()` | 10,000 | 52.2 s | — |
| `aw_boot_pls(cores = 4)` | 10,000 | 25.7 s | 0.49 of its time |
| `bootstrap_model()`, all cores | 10,000 | 370.1 s | 7.09× faster |

At 10,000 resamples on the same four cores, 370.1 / 25.7 = **14.4 times faster**; from one core, **7.1 times faster** than seminr on four. This is a single round, against three in §3, and seminr's time varies on this laptop, so these figures stand beside the headline of §3 rather than replacing it. A resample of this model takes the engine about 5 ms against about 2.5 ms for the ECSI model: the sample is larger (344 cases against 250), there are 31 indicators against 24, and each formative construct adds the inverse of its indicator correlations (once per resample) and a small matrix product (once per iteration).

Equivalence: the arrays of path coefficients, total effects, loadings and weights are identical to seminr's with all cores and with one core, at 1,000 and at 10,000 resamples. HTMT, counted in units in the last place over the 28,000 values of 1,000 resamples: 70.9% identical, 23.8% one unit apart, 5.2% two, 27 values three units and one value four units apart (largest relative difference 4.9e-16). The same count for the ECSI model of §3 gives 75.1% identical and at most three units apart; until this measurement the documentation said "within two units", which was too tight (corrected in [ADR-0006](../../docs/adr/0006-pls-sem-engine.md)). The first call, in a fresh cache, took 24.1 s.

### 3c. Reflective constructs: consistent PLS on the full ECSI model

Design: the ECSI model of §3 with all seven constructs `reflective()`, which seminr estimates by consistent PLS (PLSc): 7 constructs, 24 indicators (Complaints a single item), 12 paths, n = 250; path weighting scheme, `maxIt` = 3,000 as the Primer recommends (`bootstrap_model()` estimates the resamples with its default 300, and so does the engine; [ADR-0006](../../docs/adr/0006-pls-sem-engine.md)), stop criterion 1e-7, seed 123. Version 0.1.0.9001 against **seminr 2.6.0**, the current CRAN release (§1 to §3b used 2.5.0), the same laptop, 29 October 2026 ([`pls_plsc_bench.R`](pls_plsc_bench.R), [`pls_plsc_timing.csv`](pls_plsc_timing.csv)). The 1,000-resample times are medians of three runs; the 10,000-resample times come from one round, seminr first.

| Contender | Resamples | Time | `aw_boot_pls()` on 1 core is |
|---|---|---|---|
| `aw_boot_pls()` | 1,000 | 2,905 ms | — |
| `aw_boot_pls(cores = 4)` | 1,000 | 3,309 ms | 1.14 of its time |
| `bootstrap_model()`, all cores | 1,000 | 18,340 ms | 6.31× faster |
| `bootstrap_model(cores = 1)` | 1,000 | 44,520 ms | 15.33× faster |
| `aw_boot_pls()` | 10,000 | 26.1 s | — |
| `aw_boot_pls(cores = 4)` | 10,000 | 20.0 s | 0.77 of its time |
| `bootstrap_model()`, all cores | 10,000 | 169.8 s | 6.51× faster |

At 10,000 resamples on the same four cores, 169.8 / 20.0 = **8.5 times faster**; from one core, **6.5 times faster** than seminr on four. This is a single round, as in §3b, and seminr's time varies on this laptop, so these figures stand beside the headline of §3 rather than replacing it. At 1,000 resamples four worker processes cost more to start than they save (3.3 s against 2.9 s on one core); at 10,000 they pay off. Consistent PLS costs the engine about a fifth more than Mode A on the same model: 2,905 ms per 1,000 resamples against 2,457 ms in §3 (the reliability of each construct, path regressions on the corrected correlations, and the admissibility checks).

Resamples left out: seminr 2.6.0 kept 370 of 1,000 and 3,679 of 10,000 resamples. It leaves out every resample whose consistent-PLS solution is undefined or inadmissible (a rho_A or a reflective loading above one, among other checks; §6), and this model, with a single-item construct and seven corrected blocks, fails them often. Both programs estimate every resample, so both times include the ones left out; the engine leaves out the same resamples.

Equivalence: the arrays of path coefficients, total effects, loadings and weights, and the number of resamples kept, are identical to seminr's with all cores and with one core at 1,000 resamples, and at 10,000. HTMT, counted in units in the last place over the 7,770 values of the 370 kept resamples: 74.5% identical, 20.3% one unit apart, 5.1% two, 9 values three units apart (largest relative difference 4.2e-16). The first call, in a fresh cache, took 44.4 s: the kernel, with consistent PLS and its checks, is the largest function Arwad compiles so far.

### 3d. Moderation: a two-stage interaction term in the Primer's corporate reputation model

Design: the moderation model of Hair, Hult, Ringle and Sarstedt's *Primer* on seminr's `corp_rep_data`: the extended corporate reputation model of §3b (QUAL, PERF, CSOR and ATTR formative (Mode B); COMP, LIKE, CUSA (a single item) and CUSL Mode A composites) with switching costs (SC, four indicators, a Mode A composite) and the interaction term CUSA*SC, so that switching costs moderate the effect of CUSA on CUSL. The term is built by the two-stage approach, seminr's default and the one the Primer recommends: 10 constructs counting the term, 35 indicators and the term's product indicator, 15 paths, n = 344; path weighting scheme, `maxIt` = 3,000 as the Primer recommends (`bootstrap_model()` estimates the resamples with its default 300, and so does the engine), stop criterion 1e-7, seed 123. Both programs estimate every resample twice: first the model without the term and its path (always with `maxIt` = 300 and stop criterion 1e-7), whose scores of CUSA and SC give the term's indicator as their product, then the full model ([ADR-0006](../../docs/adr/0006-pls-sem-engine.md), amendment "two-stage interaction terms"). Version 0.1.0.9003 against seminr 2.6.0, the same laptop, 2 November 2026 ([`pls_moderation_bench.R`](pls_moderation_bench.R), [`pls_moderation_timing.csv`](pls_moderation_timing.csv)). The 1,000-resample times are medians of three runs; the 10,000-resample times come from one round, seminr first.

| Contender | Resamples | Time | `aw_boot_pls()` on 1 core is |
|---|---|---|---|
| `aw_boot_pls()` | 1,000 | 12,360 ms | — |
| `aw_boot_pls(cores = 4)` | 1,000 | 8,097 ms | 0.65 of its time |
| `bootstrap_model()`, all cores | 1,000 | 39,130 ms | 3.17× faster |
| `bootstrap_model(cores = 1)` | 1,000 | 101,600 ms | 8.22× faster |
| `aw_boot_pls()` | 10,000 | 126.0 s | — |
| `aw_boot_pls(cores = 4)` | 10,000 | 97.7 s | 0.78 of its time |
| `bootstrap_model()`, all cores | 10,000 | 386.9 s | 3.07× faster |

At 10,000 resamples on the same four cores, 386.9 / 97.7 = **4.0 times faster**; from one core, **3.1 times faster** than seminr on four. This is a single round, as in §3b and §3c, and seminr's time varies on this laptop, so these figures stand beside the headline of §3 rather than replacing it.

The gain is smaller than for the extended model alone (§3b, 14.4 times). One resample takes the engine about 12.4 ms on one core, against about 5.0 ms in §3b (same data, n = 344): every resample is estimated twice, and the model has SC's four indicators and the product besides. seminr's time per resample on one core grows less, from 77.2 ms to 101.6 ms, so core for core the engine is 8.2 times faster at 1,000 resamples, against 15.4 times in §3b. The engine's first stage still computes the loadings, path coefficients, total effects and HTMT of the model without the term, of which only the weights are used; leaving them out is the next optimisation of this path. Four worker processes took 0.65 of the one-core time at 1,000 resamples and 0.78 at 10,000 (0.84 and 0.49 in §3b).

Equivalence: the arrays of path coefficients, total effects, loadings and weights, and the number of resamples, are identical to seminr's with all cores and with one core at 1,000 resamples, and at 10,000. HTMT, counted in units in the last place over the 45,000 values of the 1,000 resamples (45 construct pairs each): 73.3% identical, 21.9% one unit apart, 4.7% two, 32 values three units and one value four units apart (largest relative difference 4.9e-16). The tests, with three resamples per model, had shown at most two units for models with interaction terms, and the first documentation of interaction terms said "within two units"; it now gives the bound the tests check, four units, as for every other model ([ADR-0006](../../docs/adr/0006-pls-sem-engine.md), amendment "two-stage interaction terms", and NEWS). The first call, in a fresh cache, took 55.3 s.

A first run, on 29 October 2026, overlapped other work on the laptop, and its one-core time for 10,000 resamples (220.8 s) was 1.8 times this run's, which scales with the number of resamples (126.0 s against ten times 12.36 s); the benchmark was therefore run again with no other work alongside, and this section reports the repeat. The first run's figures were more favourable to Arwad, because seminr took 521.0 s for 10,000 resamples in it (5.9 times faster on four cores); they are not used, and its log is kept outside the repository.

## 4. Validation (Phase 2 rule: at least 50 datasets)

[`r/tests/testthat/test-engines.R`](../../r/tests/testthat/test-engines.R), run in CI on Linux, macOS and Windows. The designs are 30 numeric data frames of R's `datasets` package (response = first numeric column, design = intercept plus the other numeric columns) and 25 seeded synthetic designs (n from 15 to 1,000, k from 2 to 6, predictor scales from 0.1 to 100).

- `.lm.fit()` compiled is bit-identical to R's on all 55 designs. Rank-deficient designs run natively and are pivoted exactly as R pivots them.
- `aw_boot_ols()` equals the same bootstrap run in R, bit for bit, on 50 designs.
- `aw_glm_fit()` equals `glm.fit()` on 55 designs for each family (binomial: the response plus seeded noise, split at the median; Poisson: the response shifted and rounded to counts; Gaussian: as is): the same coefficients, deviance, iteration count and warnings, bit for bit. In the local run, 163 of the 165 fits ran natively; the other two are fits for which `glm.fit()` warns, and they are handed to `glm.fit()` itself.
- `aw_boot_glm()` equals the same `glm.fit()` loop in R on 50 designs, with the same warnings and the same final random-number state. Replicates for which `glm.fit()` warns are refitted by `glm.fit()` from the same resample.
- `colMeans()`, `cov()` and `cor()`, compiled for the PLS-SEM engine, equal R's on 300 random cases per run, including data offset by 1e8, and so do `cov()` and `cor()` of matrix columns with kept column statistics when the matrices change between calls (40 cases).
- PLS-SEM: 125 models (54 until 28 October 2026, session 4, when Mode B was added; 80 until session 6, when consistent PLS was added): 14 on seminr's own datasets (the ECSI model in four variants, among them `maxIt` = 3,000 as the Primer recommends; the corporate reputation models of the Primer, simple and extended; a second corporate reputation sample; the influencer model; with formative (Mode B) constructs, the Primer's extended corporate reputation model with its four formative drivers and two ECSI variants; with `reflective()` constructs, which seminr estimates by consistent PLS, the ECSI model with four and with all seven constructs reflective (Complaints a single indicator), the full ECSI model with four reflective constructs next to three Mode A composites, and the Primer's extended model with reflective outcomes and formative drivers) and 111 seeded synthetic models (46 with Mode A composites: 2 to 7 constructs, 1 to 6 indicators, n from 40 to 800, Likert-type and rescaled data, shuffled construct and path orders; 23 of them again with their odd-numbered multi-indicator constructs in Mode B; and 42 with reflective constructs: 15 with every construct reflective, 15 with their odd-numbered constructs reflective next to Mode A composites, 12 Mode B models with their other constructs reflective). All 125 are in the native subset, 46 of them with reflective constructs and 39 with Mode B constructs. On every one the compiled estimate equals `estimate_pls()`: weights, the full loading matrix (signed zeros included), path coefficients, iteration count and final weight change. `aw_boot_pls()` equals `bootstrap_model()` on all 125: the bootstrap arrays and their descriptive tables bit for bit, HTMT within four units in the last place, and the resamples seminr leaves out are left out in the same way: those for which consistent PLS is undefined and, from seminr 2.6.0 on, those whose solution is inadmissible (a 30-resample test checks a model with 4 such resamples under seminr 2.5.0 and 27 under 2.6.0; a separate run of the 46 models with reflective constructs at 100 resamples each kept the same 2,922 of 4,600 resamples as seminr 2.6.0, [ADR-0006](../../docs/adr/0006-pls-sem-engine.md), amendment "seminr 2.6.0"). The tests run against seminr 2.6.0 on Linux, macOS and Windows and against 2.5.0 on Linux. Two further tests check that models estimated with other `maxIt` and `stopCriterion` values are bootstrapped with seminr's defaults, as `bootstrap_model()` does, and that a bootstrap with fewer than two usable resamples ends as seminr ends it. `summary()` gives identical tables of paths, weights, loadings and total effects, and the caller's random-number state is left as seminr leaves it. The test also checks that every call ran as native code (no fallback to GNU R).
- PLS-SEM with two-stage interaction terms (moderation; since 29 October 2026, version 0.1.0.9003): 34 more models, 159 PLS-SEM models in all. They are the Primer's corporate reputation model with switching costs moderating CUSA → CUSL (§3d); ECSI variants with the term first, between or after the constructs, with two terms, a Mode B iv, an endogenous iv, a term into a mediator, a quadratic term, and other `maxIt` and stop criteria; and 24 seeded synthetic models with an interaction effect in the data, some with Mode B constructs and some with Likert-type data. All 34 run natively. On every one the compiled estimate equals `estimate_pls()` (weights, loadings, paths, iteration count and final weight change), the bootstrap equals `bootstrap_model()` (the arrays and their descriptive tables bit for bit, HTMT within four units in the last place), and the leave-one-out estimates behind BCa intervals equal `estimate_pls()` on the data less one case, under seminr 2.6.0 and 2.5.0 ([ADR-0006](../../docs/adr/0006-pls-sem-engine.md), amendment "two-stage interaction terms"). Product-indicator and orthogonal terms, and interaction terms in models with `reflective()` constructs, are handed to seminr, and a test checks that they are.

## 5. Interpretation

1. **Against the code researchers actually write, the gain is large.** Most of the time of the `lm()` and `glm()` loops is spent building model frames and formula machinery for every replicate; the engines do only the arithmetic.
2. **Against the best hand-tuned R, the gain is modest but real:** 1.9 times for OLS and 3.25 times for GLM. Those loops already spend much of their time in the same Fortran QR; the engines save R's per-iteration allocation and copying.
3. **Same answers, by construction.** The engines call R's own random-number stream and R's own QR, and the GLM algorithm was matched to `glm.fit()` by black-box tests, so there is no "Arwad answer" to explain to a reviewer: it is R's answer, sooner.
4. **The PLS-SEM engine wins twice:** on arithmetic (one core of compiled code against seminr's R) and on overhead (seminr runs the whole estimation function for every resample and starts a cluster of R processes for every bootstrap). It does so on one core; with `cores = 4` it halves the time again for 10,000 resamples of the ECSI model on this two-core laptop (less with consistent PLS and with an interaction term, §3c and §3d). Where the engine's arithmetic doubles, as with the two estimations per resample of an interaction term, the gain shrinks (§3d: 4.0 times on four cores).
5. **The GLM engine has room to grow.** About 0.38 ms per logistic fit is well above what the arithmetic needs; profiling the kernel (allocation per iteration, repeated `exp()` calls) is the next optimisation.

## 6. Found on the way

Matching `glm.fit()` bit for bit exposed two rounding differences in the compiler, fixed for every compiled function: `x ^ 2` is now computed as `x * x`, as R's arithmetic computes it (R's C function `R_pow()` can round squares differently on 64-bit Windows), and `sqrt()` is always the processor's correctly rounded square root (units are built with `-fno-math-errno`; otherwise GCC may call the C library's `sqrt()`, which on Windows rounds some results to the neighbouring double).

Matching seminr bit for bit required numerical details that no paper specifies, all established from seminr's outputs (clean room; [ADR-0006](../../docs/adr/0006-pls-sem-engine.md)): weights are normalised by the uncentred root mean square of their scores, not by the standard deviation; scores and inner proxies are re-standardised at every iteration; resample `b` uses `set.seed(seed + b)`, which is why seminr's results do not depend on the number of cores; HTMT uses absolute correlations and treats a single indicator's within-construct correlation as 1; cross-loadings are covariances multiplied by zero, so a negative covariance leaves a negative zero; and seminr's list of indicators becomes a matrix when every construct has the same number of indicators. Three statistics had to be compiled to R's exact rounding: `colMeans()`, `cov()` and `cor()` (of two vectors and of two matrices; `cor()` of a single matrix rounds differently). Mode B added two more details: of ten numerical forms of the same regression, only the explicit inverse of the block's indicator correlations times their correlations with the inner proxy, `solve(cor(Zj)) %*% cor(Zj, proxy)`, reproduces seminr's weights (the other nine differ by 2e-16 to 3e-15); and that needs `cor()` of a single matrix, which R computes as the covariance matrix divided by the double-precision square roots of its diagonal (diagonal exactly 1, the rest clamped to [−1, 1]), now compiled to that rule.

Consistent PLS added three findings ([ADR-0006](../../docs/adr/0006-pls-sem-engine.md), amendment of session 6). First, in a model with a reflective construct, seminr 2.5.0 corrects the correlations of every Mode A construct with two or more indicators by its rho_A, `composite()` constructs included, whereas the published method corrects only constructs modelled as common factors; the engine keeps seminr's rule, so that the bootstrap belongs to the estimate seminr produced. Second, seminr leaves out the resamples for which consistent PLS is undefined (a negative rho_A) and counts only the others in `boots`; the engine does the same and says so in a message. Third, and not limited to consistent PLS: `bootstrap_model()` estimates every resample with its default `maxIt` (300) and stop criterion (1e-7), whatever the model was estimated with. Until 0.1.0 the engine used the model's own settings, which gives different numbers for a model estimated with another stop criterion, or for a resample that needs more than 300 iterations under `maxIt` = 3,000. None of the published comparisons met either case (their arrays were identical to seminr's), but the engine now follows seminr.

seminr 2.6.0 (CRAN, 13 October 2026) changed consistent PLS, and the CI found it: the CI installs the current release while this laptop still had 2.5.0, so the consistent-PLS commit passed here and failed on every CI system. Version 2.6.0 treats composites as fully reliable, the published rule (the first finding above now describes 2.5.0 only), and leaves out of a bootstrap every resample whose solution is inadmissible. The thresholds of its five admissibility checks had to be measured by bisection: a reflective loading counts as above one only beyond 1 + 5e-9, and the construct correlation matrix implied by the structural model fails when its smallest eigenvalue is below 1e-8 ([ADR-0006](../../docs/adr/0006-pls-sem-engine.md), amendment "seminr 2.6.0"). The engine now follows the installed version, and the CI tests both.

A mistake caught by the benchmark: the first version of the PLS-SEM kernel assigned whole matrix columns, which the compiler does not support yet, and `aw_compile()` then runs a function in GNU R without an error. Its first benchmark therefore timed interpreted R (36.8 s for 1,000 resamples, slower than seminr on four cores) and its first validation run compared interpreted R with seminr. The kernel now assigns element by element, every result above comes from native code, and the tests now fail if the engine's code is not compiled or any call falls back to GNU R.

## 7. Environment and next steps

Intel Core i7-7600U (2 cores, 4 logical, 16 GB of memory), Windows 11, R 4.6.0, seminr 2.5.0 (§1 to §3b) and 2.6.0 (§3c and §3d), Goose v2026.3, Rtools45 GCC 14.3. Mode B joined the native path on 28 October 2026, bias-corrected percentile intervals followed the same day (`aw_pls_ci()`, [ADR-0007](../../docs/adr/0007-pls-sem-confidence-intervals.md)), and consistent PLS in the evening (version 0.1.0.9000, matched to seminr 2.6.0 on 29 October, version 0.1.0.9001); BCa intervals and two-stage interaction terms followed on 29 October (versions 0.1.0.9002 and 0.1.0.9003). Next (ROADMAP Phase 2): product-indicator and orthogonal interaction terms, and higher-order constructs, in the native PLS-SEM path; then the remaining engines (DEA, panel models, ARIMA/GARCH, meta-analysis).

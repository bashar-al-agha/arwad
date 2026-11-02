# Roadmap

Durations are estimates for a solo project lead working with AI agent teams ([AGENT_TEAM.md](AGENT_TEAM.md)). They will be re-estimated once Phase 0 produces real measurements.

## Phase 0: Foundations (2–4 weeks)

| # | Task | Output | Status (24 Oct 2026) |
|---|---|---|---|
| 0.1 | Install Goose v2026.3 (Windows x64 locally, Linux x64 in CI); compile and run the bundled samples; record versions | `bench/phase0/REPORT.md` §2 | ✅ Windows, TinyCC back end |
| 0.2 | Hand-write a Goose OLS kernel: Householder QR without pivoting for full-rank designs | `bench/phase0/ols_bootstrap.goose` | ✅ spike |
| 0.3 | Hand-write a nonparametric case bootstrap: B = 10,000, n = 1,000, k = 10, with resampling indices identical to R's | `bench/phase0/ols_bootstrap.goose` | ✅ spike (shared LCG; R's own RNG in Phase 1) |
| 0.4 | Bind both to R through a minimal C shim and `.Call`; confirm inputs are read without copying | [`poc/`](../poc/) | ✅ 1.59× faster than R's best loop; 46× less R memory |
| 0.5 | Equivalence harness; report maximum absolute and relative error | `bench/phase0/compare.R`, `poc/run_poc.R` | ✅ max relative error 8.5e-15 |
| 0.6 | Benchmark harness ([BENCHMARK_METHOD.md](BENCHMARK_METHOD.md)): wall time and peak memory, including an optimising back end | `bench/phase0/` | ◐ TinyCC and GCC -O2 timed; R-side memory measured; parallel and whole-process memory pending |
| 0.7 | Report with raw data, whatever the result | [`bench/phase0/REPORT.md`](../bench/phase0/REPORT.md) | ✅ preliminary |

**Exit:** the report exists; equivalence is within 1e-9 relative error; a decision memo says *proceed*, *adjust* or *change back end*. **Met (24 Oct 2026):** the report's verdict is *proceed*, with the strategy adjustment in its §5 (call BLAS/LAPACK for dense work; target R's per-iteration overheads). Still open from 0.6: parallel workers and whole-process memory, now tracked under Phase 3 and Phase 1 respectively.

## Phase 1: Compiler vertical slice (2–4 months)

| # | Task | Status (26 Oct 2026) |
|---|---|---|
| 1.1 | `arwad` R package skeleton ([`r/`](../r/)) that passes `R CMD check` on Windows, Linux and macOS | ✅ Status OK on Linux, macOS (arm64) and Windows in CI, with Goose installed and the corpus compiled ([`.github/workflows/r-check.yml`](../.github/workflows/r-check.yml)) |
| 1.2 | Front end, type inference and Goose emitter; content-hash cache | ◐ slice 1: scalars and vectors, one type per variable, guarded native execution ([ADR-0005](adr/0005-guarded-native-execution.md)); typed IR and memory planner pending |
| 1.3 | `aw_compile(f)` for [spec/subset-v0.md](../spec/subset-v0.md), with fallback and a recorded reason | ◐ slice 1 of the subset (see the status column in the spec) |
| 1.4 | Build paths: an optimising C compiler through `R CMD SHLIB` (fallback); a direct-compiler path with the Goose runtime compiled once; dispatch checked in C; TinyCC for interactive use | ◐ direct path and C-side dispatch done: first call about 1 s per signature (was 3.8–5.3 s), 4.5 µs per call (was 33 µs); loops 11–33× faster than R's byte code ([bench/phase1/REPORT.md](../bench/phase1/REPORT.md)); TinyCC pending |
| 1.5 | Bit-exact R random-number generation in the static tier (Mersenne-Twister, inversion normals, rejection `sample()`) | ✅ `runif`, `rnorm`, `sample` and `sample.int` from R's own generators, bit-identical streams, seed restored on bailout (26 Oct 2026) |
| 1.6 | Test corpus of 200 functions (hand-written, agent-generated, and mined from CRAN where licences allow), each with an equivalence test | ✅ 201 functions ([`r/inst/corpus/`](../r/inst/corpus/)), all native and bit-identical to GNU R on Linux, macOS and Windows |
| 1.7 | Slice 2: vectorised arithmetic with recycling, `mean`/`min`/`max`, `c()`, masks and index vectors (2a); matrices and `m[i, j]`, `%*%`, `t`, `crossprod`, `solve` through LAPACK (2b) | ✅ 2a and 2b done (26 Oct 2026): vectors; matrices; `%*%`, `crossprod`, `solve` through R's BLAS and LAPACK |

**Exit:** every corpus function either compiles or falls back correctly; all compiled functions match GNU R within 1e-9; no memory errors under AddressSanitizer.

**Exit met (verified 27 Oct 2026).** All 201 corpus functions compile and run natively, and their results are bit-identical to GNU R (stricter than 1e-9) on Linux, macOS (arm64) and Windows in CI; the fallback and bailout tests show that unsupported constructs and R's error cases go to GNU R with a recorded reason; the `asan` CI job built every unit with `-fsanitize=address` and ran the whole suite, 1,554 corpus checks included, with no sanitizer error; `R CMD check` is OK on all three systems. Carried into Phase 2: the typed IR and memory planner of 1.2 (code generation is direct from R's language objects today) and the optional TinyCC path of 1.4.

## Phase 2: Statistical engines (4–8 months)

Each engine is validated against its reference package on 50 datasets.

| Engine | Reference | Status (28 Oct 2026) |
|---|---|---|
| OLS, WLS, GLM (IRLS) | `stats::lm`, `stats::glm` | ◐ OLS done: `.lm.fit()` compiled through R's own LINPACK, bit-identical on 55 designs; engine `aw_boot_ols()` (bootstrap) 43.7× faster than an `lm()` loop. GLM done (27 Oct 2026): `aw_glm_fit()` and `aw_boot_glm()` for binomial, Poisson and Gaussian models, bit-identical to `glm.fit()` on 55 designs per family, warnings included ([bench/phase2/REPORT.md](../bench/phase2/REPORT.md)). WLS: weighted Gaussian fits match `glm.fit()`; a check against `lm(weights = )` is pending |
| PLS-SEM path algorithm with 10,000-replicate bootstrap and bias-corrected percentile intervals | `seminr` | ✅ Done (28 Oct 2026): path algorithm and bootstrap (27 Oct 2026): `aw_boot_pls()`, a drop-in for `seminr::bootstrap_model()` (Mode A and, since 28 Oct 2026, Mode B composites; path weighting), bit-identical to `seminr` on every native model in the test set (HTMT within four units in the last place); 10,000 resamples of the ECSI model in 12.5 s on four cores instead of 241 s for seminr on the same cores (27 s on one core; 28 Oct 2026, after the compiler learnt to keep column statistics); `cores =` runs resamples in parallel with identical results ([bench/phase2/REPORT.md](../bench/phase2/REPORT.md), [ADR-0006](adr/0006-pls-sem-engine.md)). Bias-corrected percentile intervals for every coefficient and for specific indirect effects: `aw_pls_ci()`, percentile limits bit-identical to `seminr`'s `summary()` on all 80 test models, bias correction as `boot::boot.ci()` ([ADR-0007](adr/0007-pls-sem-confidence-intervals.md)). Consistent PLS for `reflective()` constructs native since 28 Oct 2026 (session 6), bit-identical to `seminr` on 46 test models, including the resamples it leaves out; the same change made the resamples use `seminr`'s default `maxIt` and stop criterion, as `bootstrap_model()` does ([ADR-0006](adr/0006-pls-sem-engine.md), amendment). Since 29 Oct 2026 the engine follows the installed `seminr`: 2.6.0 treats composites as fully reliable and leaves out inadmissible resamples, 2.5.0 corrects Mode A composites in mixed models; the CI tests both. 10,000 resamples of the ECSI model with every construct reflective: 20.0 s on four cores against 169.8 s for `seminr` 2.6.0 (8.5×; [REPORT §3c](../bench/phase2/REPORT.md)). BCa intervals since 29 Oct 2026: `aw_pls_ci(type = "bca")`, the jackknife in the compiled kernel ([ADR-0007](adr/0007-pls-sem-confidence-intervals.md), amendment). Two-stage interaction terms (moderation, `seminr`'s default, and quadratic terms) native since 29 Oct 2026, bit-identical to `seminr` on 34 test models ([ADR-0006](adr/0006-pls-sem-engine.md), amendment); 10,000 resamples of the Primer's moderation model: 57.6 s on four cores against 390.1 s for `seminr` 2.6.0 (6.8×, three rounds; [REPORT §3d](../bench/phase2/REPORT.md)); the first stage computes only the weights since version 0.1.0.9004. Next, beyond this row: the re-estimation analyses of the next row; higher-order constructs; product-indicator and orthogonal terms |
| PLS-SEM analyses that re-estimate the model hundreds or thousands of times, as the project lead's PLS-SEM protocol requires: PLSpredict and CVPAT (k-fold cross-validated prediction, 10 folds × 10 repetitions), MICOM and the permutation test of multigroup analysis, bootstrap multigroup analysis, FIMIX-PLS segmentation (many EM starts); later CTA-PLS and the Gaussian copula's bootstrap | `seminr` (`predict_pls()`, `estimate_pls_mga()`); `cSEM` (`testCVPAT()`, `testMICOM()`, `testMGD()`); FIMIX-PLS: SmartPLS (no R implementation on the reference machine) | planned (added 2 Nov 2026). Plan: reuse the compiled PLS kernel, which already estimates any set of samples given as an index matrix (bootstrap resamples, leave-one-out samples), the input that permutation groups and training folds also are; add the prediction of held-out cases; match each reference by black-box tests ([ADR-0003](adr/0003-mit-licence-dco-clean-room.md)), bit for bit where it is deterministic, on 50 models; FIMIX-PLS against its published algorithm and SmartPLS's output for the Primer's example |
| DEA (CCR and BCC; input and output orientation) | `Benchmarking` | planned |
| Panel fixed and random effects; Hausman test | `plm` | planned |
| ARIMA likelihood; GARCH(1,1) family | `stats::arima`, `rugarch` | planned |
| Meta-analysis: fixed effect and random effects (REML) | `metafor` | planned |

Data frames enter the static tier as column sets of known types.

**Exit:** every engine passes equivalence; benchmark report published internally; preprint drafted.

## Phase 3: Parallel and WebAssembly (2–4 months)

- Share-nothing parallel resampling through Goose worker programs.
- WebAssembly build of the static tier (Goose → C → wasm) and a browser demo: load a CSV and run a 5,000-replicate bootstrap locally.
- Methods paper submitted (*Journal of Statistical Software* or *The R Journal*).
- First public release, subject to the project lead's decision in [ADR-0003](adr/0003-mit-licence-dco-clean-room.md).

**Exit:** demo live; paper submitted; v0.1.0 tagged with a Zenodo DOI.

## Phase 4: Standalone runtime (12–24 months)

- ADR choosing the implementation language of the dynamic tier.
- Evaluator for full R semantics; base and recommended packages; package loader; R C-API layer prioritised by the symbol census ([COMPATIBILITY.md](COMPATIBILITY.md)).
- Nightly CRAN dashboard.

**Exit:** at least 50% of CRAN packages pass natively (L1 + L2).

## Phase 5: Full ecosystem (24 months and beyond)

- S4 and `methods` completeness, Rcpp coverage, graphics devices, Fortran toolchain integration, tcltk, and the long tail.

**Exit:** at least 95%, then 99%, of CRAN passes natively, and Arwad runs without GNU R installed.

## Decided now, decided later

| Decided now | Decided later (by ADR) |
|---|---|
| Two-tier architecture ([ADR-0001](adr/0001-two-tier-runtime.md)) | Implementation language of the dynamic tier (start of Phase 4) |
| Goose as the first native back end ([ADR-0002](adr/0002-goose-native-backend.md)) | |
| MIT, DCO, clean room, no patent, the name Arwad ([ADR-0003](adr/0003-mit-licence-dco-clean-room.md)) | |
| Numerical-equivalence gate ([ADR-0004](adr/0004-numerical-equivalence-gate.md)) | |
| Guarded native execution with GNU R as the fallback ([ADR-0005](adr/0005-guarded-native-execution.md)) | |

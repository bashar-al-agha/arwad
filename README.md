<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/arwad-banner-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="docs/assets/arwad-banner-light.png">
    <img alt="Arwad: R, compiled. A statistical programming language." src="docs/assets/arwad-banner-light.png" width="100%">
  </picture>
</h1>

<p align="center">
  <a href="LICENSE"><img alt="License: Apache 2.0" src="https://img.shields.io/badge/License-Apache%202.0-1E3A5F?style=flat-square"></a>
  <a href="https://doi.org/10.5281/zenodo.23725077"><img alt="DOI 10.5281/zenodo.23725077" src="https://zenodo.org/badge/DOI/10.5281/zenodo.23725077.svg"></a>
  <a href="docs/ROADMAP.md"><img alt="Status: Phase 2, three engines" src="https://img.shields.io/badge/Status-Phase%202%2C%20three%20engines-1E3A5F?style=flat-square"></a>
  <a href="https://github.com/bashar-al-agha/arwad/actions/workflows/r-check.yml"><img alt="R-CMD-check" src="https://github.com/bashar-al-agha/arwad/actions/workflows/r-check.yml/badge.svg"></a>
  <a href="bench/phase2/REPORT.md"><img alt="Benchmark report" src="https://img.shields.io/badge/Benchmarks-Phase%202%20report-1E3A5F?style=flat-square"></a>
  <a href="https://github.com/aardappel/goose"><img alt="Built on Goose" src="https://img.shields.io/badge/Built%20on-Goose%20v2026.3-1E3A5F?style=flat-square"></a>
  <a href="CITATION.cff"><img alt="Cite this repository" src="https://img.shields.io/badge/Cite-CITATION.cff-1E3A5F?style=flat-square"></a>
</p>

**A statistical programming language that keeps R's syntax and its package ecosystem, and compiles numerical work to native code with deterministic, garbage-collection-free memory.**

<p align="center">
  <a href="bench/phase2/REPORT.md"><img alt="Results: bit-identical to R and seminr" src="https://img.shields.io/badge/Results-bit--identical%20to%20R%20and%20seminr-A16207?style=for-the-badge"></a>
  <a href="bench/phase2/REPORT.md"><img alt="PLS-SEM bootstrap: 19.3 times faster than seminr on the same 4 cores" src="https://img.shields.io/badge/PLS--SEM%20bootstrap-19.3%C3%97%20faster%20than%20seminr-15803D?style=for-the-badge"></a>
  <a href="bench/phase2/REPORT.md"><img alt="Regression bootstrap: 43.7 times faster than an lm() loop" src="https://img.shields.io/badge/Regression%20bootstrap-43.7%C3%97%20faster%20than%20lm%28%29%20loops-15803D?style=for-the-badge"></a>
  <a href="https://github.com/bashar-al-agha/arwad/actions/workflows/r-check.yml"><img alt="Tested on Linux, macOS and Windows" src="https://img.shields.io/badge/Tested%20on-Linux%20%C2%B7%20macOS%20%C2%B7%20Windows-1E3A5F?style=for-the-badge"></a>
</p>

> [!TIP]
> **What this changes for R and statistics.** Researchers keep writing R, and reviewers get exactly the numbers R gives: bit for bit, checked on three operating systems in every test run. The difference is time. On one laptop, a PLS-SEM bootstrap with the 10,000 resamples the Primer recommends took 12.5 s instead of 241 s with `seminr` on the same four cores (seminr's fastest of three runs; Arwad needs 27 s on a single core), and a regression bootstrap ran 43.7 times faster than the usual `lm()` loop. Compilation alone makes the same algorithm 7 to 15 times faster than plain R.

<p align="center">
  <a href="bench/phase2/REPORT.md">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="docs/assets/results-dark.svg">
      <source media="(prefers-color-scheme: light)" srcset="docs/assets/results-light.svg">
      <img alt="Phase 2 results: 43.7 times faster than an lm() loop, 10.8 times faster than a glm() loop, 19.3 times faster than seminr on the same four cores for 10,000 PLS-SEM resamples, and 7 to 15 times faster than the same algorithm in plain R; all bit-identical to GNU R and seminr (HTMT within four units in the last place)." src="docs/assets/results-light.svg" width="100%">
    </picture>
  </a>
</p>

> [!NOTE]
> **Status:** first citable release, **v0.1.0** (28 Oct 2026) · Phase 1 complete (exit criteria met 27 Oct 2026) · Phase 2 under way: engine `aw_boot_ols()` bootstraps a regression 43.7× faster than an `lm()` loop, `aw_boot_glm()` a logistic regression 10.8× faster than a `glm()` loop, and `aw_boot_pls()` a PLS-SEM model 19.3× faster than `seminr::bootstrap_model()` on the same four cores (10,000 resamples: 12.5 s instead of 241 s; with formative (Mode B) constructs, the Primer's extended corporate reputation model 14.4× faster; with reflective constructs (consistent PLS), the ECSI model 8.5× faster than seminr 2.6.0; with a two-stage interaction term (moderation), the Primer's corporate reputation model 6.8× faster), all with bit-identical results (for PLS-SEM, HTMT within four units in the last place) ([bench/phase2/REPORT.md](bench/phase2/REPORT.md)), and `aw_pls_ci()` gives the bias-corrected percentile confidence intervals the PLS-SEM Primer asks for, and BCa intervals with a compiled jackknife ([ADR-0007](docs/adr/0007-pls-sem-confidence-intervals.md)) · the `arwad` R package compiles loop-heavy R functions through Goose ([bench/phase1/REPORT.md](bench/phase1/REPORT.md)): **11–33× faster** than R's byte code on loops on one laptop (11–40× across three runs), results **bit-identical** to GNU R across a 201-function corpus on Linux, macOS and Windows; about 1 s to build each new signature, and 4.5 µs of dispatch per call · Phase 0: [bench/phase0/REPORT.md](bench/phase0/REPORT.md) · MIT licence
>
> **Project lead:** Dr. Bashar Al Agha · [ORCID 0009-0004-5178-2345](https://orcid.org/0009-0004-5178-2345)

## Install

```r
# install.packages("remotes")
remotes::install_github("bashar-al-agha/arwad", subdir = "r")
```

Compiling needs Goose v2026.3 and the C compiler R uses for packages (Rtools on Windows); `aw_toolchain()` says whether both were found. Without them every function still runs, in GNU R. A step-by-step guide is at [langwi.com/get-started](https://langwi.com/get-started).

## Cite

> Al Agha, B. (2026). *Arwad: R, compiled* (Version 0.1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23725078

That DOI is version 0.1.0's; [10.5281/zenodo.23725077](https://doi.org/10.5281/zenodo.23725077) always points to the newest release. GitHub's *Cite this repository* button reads the same details from [CITATION.cff](CITATION.cff).

## The idea

R is the working language of statistics: CRAN alone hosts 25,504 packages (census of 24 October 2026, [research/](research/)). Its interpreter was designed for interactive analysis. Every value lives on a garbage-collected heap, modifying an object can copy it whole, and loops are slow unless rewritten in C or C++.

[Goose](https://github.com/aardappel/goose) (Apache-2.0, by Wouter van Oortmerssen) is a new systems language with **no heap at all**: values live inline on data stacks that the compiler assigns, leaving a scope is the free, and the compiler emits one plain C file. Goose functions marked `export fn` can be linked into a C program as a library, and slices cross the boundary as `{data, len}` without copying.

Arwad joins the two. People keep writing R. The parts of a program that do numerical work (resampling loops, iterative estimators, matrix kernels) are compiled through Goose into native code whose memory footprint is known before it runs.

## Why "Arwad"

Arwad is Syria's only inhabited island: the ancient Phoenician city of Arados, whose sailors carried the alphabet across the Mediterranean. Like Java and Kotlin, this language is named after an island. (Working title until 24 October 2026: R-Goose.)

## Architecture in brief

| Tier | Runs | Implementation |
|---|---|---|
| **Dynamic tier** | Everything R allows: closures, environments, lazy evaluation, S3/S4, `eval`, every CRAN package | Phases 1–3: GNU R itself, with Arwad shipped as an R package. Phase 4 onward: a new R runtime under Arwad's licence, with an R C-API compatibility layer |
| **Static tier** | Type- and shape-specialised numerical code | R function → typed IR → Goose → C → native code (TinyCC in-process, or clang/MSVC), called from R through `export fn` with zero-copy slices |

Full design: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## "All of CRAN" is a design requirement

Of the 25,504 CRAN packages, 19,923 (78.1%) are pure R and 5,581 (21.9%) contain compiled code; 3,378 link to Rcpp. Arwad is compatible with all of them from its first release because its dynamic tier *is* GNU R, and acceleration is added underneath. The standalone runtime (Phase 4) must re-earn that compatibility package by package, measured nightly and published. See [docs/COMPATIBILITY.md](docs/COMPATIBILITY.md).

## Roadmap

| Phase | Deliverable | Exit criterion (measured, never assumed) |
|---|---|---|
| 0 Foundations | Goose toolchain, hand-written OLS and bootstrap kernels, benchmark and equivalence harnesses | A reproducible benchmark report against GNU R, with whatever result it shows |
| 1 Compiler slice | `arwad::aw_compile()` for the subset in [spec/subset-v0.md](spec/subset-v0.md) | 200 test functions compile and match GNU R within 1e-9 relative error |
| 2 Statistical engines | OLS/GLM, PLS-SEM with bootstrap, DEA, panel FE/RE, ARIMA/GARCH likelihoods, meta-analysis | Each engine matches its reference package on 50 datasets |
| 3 Parallel and WebAssembly | Share-nothing parallel resampling; browser build | 5,000-replicate bootstrap in the browser; methods paper submitted |
| 4 Standalone runtime | Arwad executable with its own evaluator and R C-API layer | At least 50% of CRAN passes its own examples and tests |
| 5 Full ecosystem | S4, Rcpp, graphics devices, Fortran, the long tail | At least 95%, then 99%, of CRAN |

Details, durations and dependencies: [docs/ROADMAP.md](docs/ROADMAP.md).

## Ground rules

1. **No unmeasured claims.** A speed-up is stated only with a reproducible benchmark and its raw data ([docs/BENCHMARK_METHOD.md](docs/BENCHMARK_METHOD.md)).
2. **Numerical equivalence first.** Compiled results are checked against GNU R, by default within 1e-9 relative error, with bit-identical random-number streams ([ADR-0004](docs/adr/0004-numerical-equivalence-gate.md)).
3. **Clean room.** No GNU R source code (GPL) is copied into this repository; behaviour is re-implemented from documentation and black-box tests ([docs/IP_AND_LICENSING.md](docs/IP_AND_LICENSING.md)).
4. **Credit where due.** Goose is the work of Wouter van Oortmerssen; R is the work of the R Core Team. Arwad builds on both and is affiliated with neither.

## Repository map

| Path | Contents |
|---|---|
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | Two-tier design, compiler pipeline, risks, prior art |
| [docs/COMPATIBILITY.md](docs/COMPATIBILITY.md) | CRAN census and the compatibility strategy |
| [docs/ROADMAP.md](docs/ROADMAP.md) | Phases, tasks, exit criteria |
| [docs/IP_AND_LICENSING.md](docs/IP_AND_LICENSING.md) | Copyright, the Apache License 2.0, trademark, patent assessment, licence compatibility |
| [docs/AGENT_TEAM.md](docs/AGENT_TEAM.md) | The ten builder–verifier duets and how they are run |
| [docs/BENCHMARK_METHOD.md](docs/BENCHMARK_METHOD.md) | How performance is measured and reported |
| [docs/adr/](docs/adr/) | Architecture decision records |
| [spec/subset-v0.md](spec/subset-v0.md) | The R subset the Phase 1 compiler accepts, with implementation status |
| [r/](r/) | The `arwad` R package: `aw_compile()`, its tests, and the equivalence corpus |
| [research/](research/) | CRAN census data and the script that produces it |
| [bench/phase0/](bench/phase0/) | Phase 0 spike: R reference, Goose kernel, equivalence check, raw outputs and report |
| [bench/phase1/](bench/phase1/) | Compiler slice 1 against GNU R's byte-compiled code: script, raw timings and report |

## Licence

[Apache License 2.0](LICENSE) © 2026 Bashar Al Agha; redistributions keep the [NOTICE](NOTICE) file. Versions up to 0.1.0.9004, the v0.1.0 release included, were published under the MIT licence and remain available under it ([ADR-0003](docs/adr/0003-mit-licence-dco-clean-room.md), amendment). See [TRADEMARKS.md](TRADEMARKS.md) for use of the project name.

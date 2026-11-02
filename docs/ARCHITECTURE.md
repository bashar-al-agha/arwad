# Architecture

This document describes how Arwad is built and why. Decisions with lasting consequences are recorded separately as ADRs in [adr/](adr/).

## 1. Goals and non-goals

**Goals**

- Run R code unchanged, including any CRAN package.
- Compile the numerical parts of R programs to native code with memory use known before execution and no garbage-collection pauses.
- Produce results numerically equivalent to GNU R, including identical random-number streams for a given seed.
- Remain permissively licensed (the Apache License 2.0; MIT up to version 0.1.0.9004) and clean-room with respect to GNU R's GPL source.

**Non-goals (for now)**

- Changing the R language. Arwad accepts R as written; it does not introduce new syntax before Phase 4.
- Out-performing R on code that already spends its time inside optimised BLAS/LAPACK or C libraries. The gains are expected where R's interpreter and heap dominate: loops, resampling, element-wise work, many small allocations.

## 2. Two sets of constraints

Arwad sits between a very dynamic language and a deliberately restrictive one. The design follows from where they conflict.

| R requires | Goose provides or forbids |
|---|---|
| First-class closures that capture environments with unbounded lifetimes | No escaping closures, no function pointers, no dynamic dispatch beyond ADT tags |
| Environments as mutable, cyclic graphs; promises for lazy evaluation | No heap; lifetimes follow scopes (or `reusable` pools when they do not nest) |
| Recursion is ordinary | Recursive functions cannot keep growable data across a recursive call |
| Lists and character vectors with variable-size elements, freely indexed | Arrays of variable-size elements iterate but do not index (short of a side index of references) |
| Separate loading of thousands of packages | Whole-program compilation only; no shared libraries of Goose code |
| One interpreter thread calling into C | Exported functions are not thread-safe and must be called from the thread that ran `goose_init`; threads are share-nothing workers exchanging flat values over typed queues |

Goose is therefore an excellent target for **numerical kernels** (fixed-type vectors and matrices, loops, resampling, iterative estimation) and an unsuitable host for **R's object model**. Hence two tiers.

## 3. Two tiers

```
            R source (unchanged)
                   │
        ┌──────────┴───────────┐
        │     Dynamic tier     │  closures, environments, promises, S3/S4,
        │  (GNU R, phases 1–3) │  eval, every CRAN package
        └──────────┬───────────┘
                   │ aw_compile(f) / automatic hot-spot selection
                   ▼
        ┌──────────────────────┐
        │     Static tier      │  typed IR → memory plan → Goose → C → native
        │   (compiled kernels) │  zero-copy slices in, R-allocated buffers out
        └──────────────────────┘
```

The boundary is deliberately narrow: data crosses as contiguous numeric buffers, which is how R already stores numeric, integer and logical vectors and matrices.

## 4. Static tier pipeline

### 4.1 Front end

- **Phases 1–3:** R's own parser (`parse()` with `getParseData()`) produces the syntax tree, so Arwad parses exactly what R parses at no cost.
- **Phase 4:** the standalone runtime needs its own parser. Candidate: [tree-sitter-r](https://github.com/r-lib/tree-sitter-r) (MIT), validated against R's parser on all CRAN sources.

### 4.2 Intermediate representation (RG-IR)

A typed SSA representation of one R function at a time.

- **Types:** logical, integer, double, complex; character only in later versions. Each value carries a static *shape*: scalar, vector of length *n*, or matrix *n × k*, with `dim` and `names` attributes tracked where known.
- **R semantics preserved explicitly:** 1-based indexing; recycling rules and their warnings; `NA` propagation (integer `NA` is `INT_MIN`; `NA_real_` is a NaN with payload 1954, distinct from other NaNs); integer overflow producing `NA` with a warning; `/` always returning double; R's definitions of `%%` and `%/%`; `1:0` yielding `c(1, 0)`; `for` over zero-length sequences.

### 4.3 Specialisation and type inference

Functions are compiled for the concrete argument types and shapes seen at the call (as in Julia or Numba). The compiled entry checks its assumptions; if a later call violates them, the call falls back to the R version and may trigger a new specialisation. Inference is abstract interpretation over RG-IR.

### 4.4 Memory planner

The planner maps every IR value to Goose storage:

- fixed-size arrays when the shape is known;
- grow-only arrays (`T[>..]`) for results built incrementally;
- loop-body scopes for per-iteration temporaries, so a resampling iteration's scratch memory is released by a single store when the iteration ends. This is the main expected win for bootstrap and Monte Carlo work.

The plan yields the static tier's peak memory before execution; results report *planned bytes* next to *measured bytes*.

### 4.5 Code generation

The emitter writes one Goose module per compilation unit, with `export fn` entry points that take `f64[:]`, `i32[:]` and similar slices, and compiles it with `goose --standalone --header`. Linear algebra calls BLAS and LAPACK through `extern fn`, using the **same** BLAS that the host R uses (reference BLAS, OpenBLAS or MKL), so differences reflect algorithms rather than libraries. Benchmarks always state which BLAS was used.

### 4.6 Native build and cache

- Interactive use: Goose's bundled TinyCC backend compiles in-process with no external toolchain.
- Production: clang, GCC or MSVC at `-O2`.
- Each compiled unit is cached by a hash of (RG-IR, argument types, Goose version, C compiler and flags).
- Phase 1 builds with R's own compile and link commands, read once from a dry run of `R CMD SHLIB`, and links one Goose runtime object, compiled once per toolchain, into every unit (about 1 s per unit on the Phase 0 laptop); `R CMD SHLIB` with a standalone unit is the fallback.
- Because Goose compiles whole programs, each unit becomes its own small native library with its own `goose_init`; the binding layer keeps their symbols separate.

### 4.7 Binding to R

The R package `arwad` loads compiled units and calls them through `.Call`:

- inputs: `REAL(x)` / `INTEGER(x)` / `LOGICAL(x)` with `XLENGTH(x)` passed as `{data, len}` slices, read without copying;
- outputs: allocated by R (`allocVector`) and passed in as writable slices, so nothing is copied back;
- calls happen on R's main thread, which is also the thread that ran `goose_init`, satisfying Goose's threading rule.

### 4.8 Parallelism

Parallel resampling runs inside Goose: worker programs (`thread_fn`) receive flat jobs (seeds or index vectors) over typed queues and return flat results. No R process is forked and no dataset is duplicated per worker beyond what each job needs.

## 5. Dynamic tier

### 5.1 Phases 1–3: GNU R

Arwad ships as an R package (`library(arwad)`). All of R and CRAN keeps working because nothing is replaced, and `aw_compile()` returns the original function, with the reason recorded, whenever it cannot compile it.

### 5.2 Phase 4 onward: standalone runtime

A new R runtime, written from documentation, with:

- an evaluator for full R semantics (closures, promises, environments, conditions and restarts, S3, S4, reference classes). These are graph-shaped objects with unbounded lifetimes, so this tier needs a garbage collector and is written in a GC-capable implementation language (Rust, C++ or Zig; decided by ADR at the start of Phase 4);
- an **R C-API compatibility layer** that offers the documented API of *Writing R Extensions* so compiled CRAN packages load unchanged; it is prioritised by measured symbol usage across CRAN ([COMPATIBILITY.md](COMPATIBILITY.md));
- a package loader for `DESCRIPTION`, `NAMESPACE` and lazy-load databases;
- re-implementations of the base and recommended packages, which is the largest single body of work.

The static tier is unchanged; it simply gains a second host.

## 6. Determinism and numerical equivalence

- **Random numbers:** R's default generators (Mersenne-Twister; inversion for normals; rejection sampling for `sample()` since R 3.6.0) are re-implemented bit-for-bit, so `set.seed(n)` yields identical streams in R and Arwad.
- **Floating point:** no fast-math. Reductions use a documented order. R accumulates `sum()` and `mean()` in `long double` (80-bit on x86-64 GCC builds, including R for Windows); the static tier matches that per back end, and back ends that cannot are flagged.
- **Reference data:** equivalence suites include the NIST Statistical Reference Datasets for linear regression (Longley, Norris, Pontius, Filip, Wampler), whose certified values test ill-conditioned problems.

## 7. Risks and open questions

| Risk | Mitigation |
|---|---|
| Goose is very young (repository created 30 Aug 2026; release v2026.3 on 21 Oct 2026; single author; README states it was produced almost entirely with AI from a human design) | Pin Goose versions; contribute fixes upstream; keep RG-IR back-end-agnostic so a direct C or LLVM back end can replace Goose without touching the front end |
| Whole-program compilation | One native library per compilation unit, cached by content hash |
| No closures or function pointers in Goose | Higher-order calls (`vapply(x, f, 0)`, `optim(par, fn)`) compile only when `f` is known at compile time and can be specialised; otherwise they stay in the dynamic tier |
| Recursion limits | Recursive R functions stay in the dynamic tier until a pool-based lowering is designed |
| Lists, character vectors and data frames | Phase 1 handles numeric data only; data frames arrive as column sets of known types in Phase 2 |
| Compatibility long tail in Phase 4 | Measured nightly; GNU R remains the host until the standalone runtime passes the agreed thresholds |

## 8. Prior art and lessons

- **R's byte-code compiler** (the `compiler` package, since R 2.13.0; JIT by default since R 3.4.0) shows how far R can go without changing its object model.
- **FastR** (Oracle Labs, GraalVM) reached broad package compatibility through a C-API emulation layer but did not displace GNU R.
- **Renjin** (JVM) compiled R, and even packages' C and Fortran code, to JVM bytecode; its development has since stalled.
- **pqR** and **Riposte** showed large speed-ups on specific workloads through interpreter and vector-trace optimisations.
- **Numba** and **JAX** in Python show the accelerator pattern Arwad starts with: compile a well-defined subset, leave the rest to the host language.

The lesson Arwad takes from all of them: compatibility is the long tail, so deliver value through acceleration first and earn independence from GNU R measurably.

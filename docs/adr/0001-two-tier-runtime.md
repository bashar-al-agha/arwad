# ADR-0001: Two-tier runtime

- **Status:** Accepted
- **Date:** 2026-10-24

## Context

R's semantics (first-class closures, mutable cyclic environments, promises, S3/S4 dispatch, `eval`) require objects with unbounded, graph-shaped lifetimes. Goose forbids a heap, escaping closures and function pointers, and compiles whole programs only. CRAN holds 25,504 packages, 21.9% of them with compiled code written against the R C API.

## Decision

Split Arwad into a **dynamic tier** that implements R's object model and hosts packages, and a **static tier** that compiles type- and shape-specialised numerical code through Goose. In Phases 1–3 the dynamic tier is GNU R and Arwad is an R package; from Phase 4 a standalone dynamic tier under Arwad's own permissive licence ([ADR-0003](0003-mit-licence-dco-clean-room.md)) is built with an R C-API compatibility layer.

## Consequences

- Full CRAN compatibility from the first release, because GNU R hosts everything.
- Value is delivered early, through acceleration of the code that benefits most.
- The standalone runtime becomes a separate, measurable programme of work rather than a precondition.
- The boundary between tiers must stay narrow: contiguous numeric buffers in, R-allocated buffers out.

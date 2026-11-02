# Intellectual property and licensing

*This document records the project's position and the reasoning behind it. It is not legal advice; obtain professional advice before filing any trademark or patent application.*

## 1. Copyright

Copyright in the code arises automatically when it is written (Berne Convention; in India, computer programs are literary works under the Copyright Act, 1957). No registration is needed. Optional registration with the Indian Copyright Office adds evidential weight in a dispute. The copyright holder for code written by or for the project lead is **Bashar Al Agha**.

## 2. The Apache License 2.0

Arwad is licensed under the Apache License 2.0 from version 0.1.0.9005 (2 November 2026). Earlier versions, the v0.1.0 release included, were published under the MIT licence ([ADR-0003](adr/0003-mit-licence-dco-clean-room.md), amendment).

- **It is self-applied.** There is no application, approval or fee: placing the licence text in `LICENSE` *is* the licence grant. GitHub then labels the repository "Apache-2.0".
- **What it allows:** anyone may use, copy, modify, distribute, sublicense and sell the software, provided they keep the copyright notices, give recipients a copy of the licence, mark the files they changed, and pass on the `NOTICE` file (Section 4).
- **What it adds to MIT:**
  - every contributor grants a licence to the patents their contribution practises, and that licence ends for anyone who sues claiming the software infringes a patent (Section 3);
  - it grants no permission to use the project's name or trademarks (Section 6), which supports the trademark registration in §3;
  - the `NOTICE` file, which carries the attribution to the copyright holder into redistributions and derivative works (Section 4(d)).
- **What it does not do:** like every permissive licence, it does not stop anyone from copying, selling or making a closed-source fork of the code; only a copyleft licence does that, at the cost of adoption. The *name* is protected by the trademark, not by the licence.
- **What the author keeps:** the copyright itself, and the requirement that every copy carry the author's notices.
- **What cannot be undone:** a version released under a licence stays available under it to whoever has it. Every version published before 0.1.0.9005 therefore remains available under MIT; later versions may again be released under different terms, because the copyright holder decides the licence of each release.

## 3. How an open-source project is "owned"

| Lever | What it protects | Action |
|---|---|---|
| Copyright | The code as written | Automatic; keep the notices in `LICENSE` and `NOTICE` |
| **Trademark** | The *name* and *logo*: only the owner can call a release "Arwad" | Clearance search, then register "Arwad" (India: classes 9 and 42), then extend abroad through the Madrid Protocol |
| Control of the official repository and releases | Which build is the authentic one | Project lead holds owner rights; releases are signed and tagged |
| Priority and credit | Who did it first | Dated commits; Zenodo DOI per release; peer-reviewed methods paper |

Trademark is the strongest lever: Linux, Python and Rust are open source, yet their names are controlled.

## 4. The project name

**Arwad** (confirmed 24 October 2026; the working title was "R-Goose"). Arwad is Syria's only inhabited island, the ancient Phoenician city of Arados; the Phoenicians carried the alphabet across the Mediterranean. Like Java and Kotlin, the language is named after an island.

| Candidate | Why it was set aside |
|---|---|
| R-Goose | "Goose" is the name of Wouter van Oortmerssen's language and of Block's open-source AI agent *goose* (Linux Foundation, Agentic AI Foundation): weak as a trademark and suggests an affiliation that does not exist |
| R++ | Already in use: commercial statistics software by Christophe Genolini that generates and runs R code, and earlier a Bell Labs rule-based extension of C++; `+` is not allowed in GitHub repository names, CRAN package names or domain names, and search engines ignore it |
| Anser, Ebla, Statera | Squatted GitHub and npm names; search results dominated by other meanings |

Availability on 24 October 2026: `arwad` was free on CRAN, PyPI and npm; `arwad.dev`, `arwadlang.org` and `arwad-lang.org` showed no registration; a quick search found no software trademark. A professional clearance search must precede the trademark application.

Arwad does not use the R logo and does not imply endorsement by the R Foundation or the R Core Team. Decision recorded in [ADR-0003](adr/0003-mit-licence-dco-clean-room.md).

## 5. Patent assessment

| Question | Assessment |
|---|---|
| Does publishing destroy novelty? | Yes in India and most countries. India's grace provisions are narrow (for example, certain exhibitions or papers read before learned societies, within 12 months). The United States allows one year after the inventor's own disclosure; Europe has essentially no grace period. Pushing code to a public repository is publication. |
| Is software patentable in India? | "A computer programme per se" is excluded by section 3(k) of the Patents Act; an application must claim a technical effect (for example, a measured reduction in memory use by a specific method), following the CRI Guidelines. |
| What prior art exists? | Goose's memory model is public (since 30 August 2026); compiling R or R subsets to native code has a long history (Section 8 of [ARCHITECTURE.md](ARCHITECTURE.md)). A claim would have to rest on a specific, non-obvious technical method, such as the memory planner for statistical resampling. |
| Cost | Professional drafting and prosecution, plus foreign filings for international protection. |

**Decision (24 October 2026): no patent.** The plan is open source worldwide (MIT, the Apache License 2.0 since 2 November 2026), plus a registered trademark (India first, then the Madrid Protocol for other countries), plus a peer-reviewed methods paper. Recorded in [ADR-0003](adr/0003-mit-licence-dco-clean-room.md).

## 6. Licence compatibility

| Component | Licence | Role in Arwad | Obligation |
|---|---|---|---|
| Goose compiler and runtime | Apache-2.0 | Build tool; its runtime C is linked into compiled kernels | Keep Goose's `LICENSE` and `NOTICE` when redistributing it or its runtime; the same licence as Arwad |
| GNU R | GPL-2 or GPL-3 | Host runtime in Phases 1–3 | Never copy GPL source into this repository (clean room) |
| R API headers (`R.h`, `Rinternals.h`, …) | LGPL-2.1 or later (per R's `COPYRIGHTS` file) | The `arwad` package compiles against them | Permitted for packages under any licence CRAN accepts, Apache 2.0 included, as thousands of CRAN packages do |
| tree-sitter-r (Phase 4 candidate) | MIT | Parser for the standalone runtime | Keep its notice |
| Reference BLAS/LAPACK, OpenBLAS | BSD-style | Linear algebra | Keep notices |
| CRAN packages | Each package's own licence (64.8% mention GPL, 28.7% MIT) | Run by Arwad; not bundled in its distribution | Each package keeps its own licence |

The Free Software Foundation considers the Apache License 2.0 compatible with GPL-3 but not with GPL-2 alone: a package licensed "GPL-2" only cannot copy Arwad's code into itself, while packages under GPL (>= 2) or GPL-3 can. Calling Arwad from a package, or listing it in `Imports`, copies none of its code.

## 7. Contributions

- **Developer Certificate of Origin** (`git commit -s`), inbound = outbound: contributions are licensed under the Apache License 2.0, as its Section 5 provides, patent grant included. Because the licence permits sublicensing, the project lead can still ship contributed code in other editions, provided the notices are kept, so a Contributor License Agreement is not required.
- **Clean room:** contributors and agents must not paste or translate GNU R source code. Behaviour is re-implemented from documentation (*R Language Definition*, *Writing R Extensions*, *R Internals*) and verified by black-box comparison with GNU R. Pull requests that contain GPL-derived code are rejected.

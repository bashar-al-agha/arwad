# ADR-0003: Licence (MIT, then the Apache License 2.0), DCO, clean room, no patent; name Arwad

- **Status:** Accepted; amended 2026-11-02 (licence: the Apache License 2.0 from version 0.1.0.9005)
- **Date:** 2026-10-24
- **Decided by:** the project lead (the patent question was delegated by him to the engineering team on 2026-10-24 and decided as below)

## Context

The project lead wants an open, globally adopted language that he owns and that carries his name, as a decisive step in his career. Ownership in open source comes from copyright, trademark, control of releases and recognised priority. Public disclosure ends the possibility of patenting in most jurisdictions. GNU R is GPL-licensed; Goose is Apache-2.0.

## Decision

1. **Licence:** MIT, copyright Bashar Al Agha, worldwide. (A dual "MIT OR Apache-2.0" licence, as used by Rust, remains an option that adds an explicit patent grant for users.)
2. **Contributions:** Developer Certificate of Origin sign-off; inbound = outbound.
3. **Clean room:** no GNU R source is copied or translated; behaviour is re-implemented from documentation and black-box tests.
4. **Patent: none.** Reasons: India's section 3(k) excludes computer programs per se; Goose's memory model is public prior art (since 30 August 2026); compiling R subsets to native code has a long history; costs are high and enforcement unrealistic for an open-source language. Global protection of the project's identity comes instead from a **registered trademark**: first in India (classes 9 and 42), then internationally through the Madrid Protocol (WIPO).
5. **Priority and credit:** a public repository with dated commits, a Zenodo DOI per release, and a peer-reviewed methods paper.

## Name

**Arwad**, confirmed by the project lead on 2026-10-24 (his first preference, R++, is already in use and cannot be a repository, package or domain name; see [IP_AND_LICENSING.md](../IP_AND_LICENSING.md) §4). The working title "R-Goose" is retired. The trademark application follows a professional clearance search.

## Consequences

- The repository became public on 2026-10-24, under the MIT licence.
- Trademark registration, not the licence, protects the project's identity worldwide.

## Amendment (2 November 2026): the Apache License 2.0

**Decided by** the project lead on 2026-11-02, choosing among the Apache License 2.0 alone (recommended to him), a dual "Apache 2.0 OR MIT" licence, and keeping MIT.

**Context.** The project lead asked whether the Apache License 2.0 could be adopted to protect Arwad. No permissive licence stops copying, selling or closed-source forks; only a copyleft licence does, at the cost of adoption by companies and package authors. What the Apache License 2.0 adds to MIT:

1. every contributor grants a licence to the patents their contribution practises, and that licence ends for anyone who sues claiming the software infringes a patent (Section 3);
2. it grants no permission to use the project's name or trademarks (Section 6), which supports the trademark registration of decision 4;
3. a `NOTICE` file that redistributions and derivative works must carry (Section 4(d)), so the attribution to the copyright holder travels with the code;
4. modified files must say that they were changed (Section 4(b)).

A dual licence was not chosen because anyone could pick MIT and skip points 1 to 4.

**Decision.** From version 0.1.0.9005 Arwad is licensed under the Apache License 2.0, copyright Bashar Al Agha, with a `NOTICE` file at the repository root and in the R package (`inst/NOTICE`). Contributions stay inbound = outbound (Section 5 of the licence) with the DCO sign-off. Decisions 3 to 5 are unchanged.

**Consequences.**

- Everything published before this change, the v0.1.0 release (DOI 10.5281/zenodo.23725078) included, remains available under the MIT licence to whoever has it. The change applies from version 0.1.0.9005 on.
- The relicensing needed no one else's consent: all 59 commits before it are the project lead's own.
- The R package declares `License: Apache License (== 2.0)`, a licence CRAN accepts, and no longer carries the MIT `LICENSE` stub.
- Goose, whose runtime is linked into compiled kernels, is under the same licence.
- The Free Software Foundation considers the Apache License 2.0 compatible with GPL-3 but not with GPL-2 alone, so a package licensed "GPL-2" only cannot copy Arwad's code into itself ([IP_AND_LICENSING.md](../IP_AND_LICENSING.md) §6).
- The file name of this record keeps "mit-licence" so that existing links keep working.

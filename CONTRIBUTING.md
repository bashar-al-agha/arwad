# Contributing to Arwad

Arwad is in Phase 0 ([docs/ROADMAP.md](docs/ROADMAP.md)). These rules apply equally to the project lead's agent team and to outside contributors. Questions and ideas are welcome as issues.

## Non-negotiable rules

1. **Clean room.** Never paste, translate or closely paraphrase GNU R source code (GPL). Work from the R manuals and from black-box comparison with GNU R. A pull request containing GPL-derived code is rejected.
2. **Equivalence gate.** Any change that can affect a computed value includes an equivalence test against GNU R ([ADR-0004](docs/adr/0004-numerical-equivalence-gate.md)).
3. **No unmeasured claims.** Performance statements follow [docs/BENCHMARK_METHOD.md](docs/BENCHMARK_METHOD.md) and include raw data.
4. **Developer Certificate of Origin.** Sign off every commit with `git commit -s`. Contributions are licensed under the Apache License 2.0 (inbound = outbound, as its Section 5 provides).

## Workflow

1. Pick or open an issue; for architectural changes, write an ADR first ([GOVERNANCE.md](GOVERNANCE.md)).
2. Branch from `main` (`feat/...`, `fix/...`), keep commits focused, and use Conventional Commits.
3. Open a pull request; a verifier from a different model family, or a different person, must approve it.
4. The project lead merges.

General conduct, security reporting and support follow the account-wide files in [bashar-al-agha/.github](https://github.com/bashar-al-agha/.github).

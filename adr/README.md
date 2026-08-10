# Architecture Decision Records

One file per architectural decision, dated and immutable. An ADR records *why a
choice was made at a point in time*, including what was rejected and what the
choice costs. It is history, not documentation of current behaviour.

Same conventions as GiottoClass's and GiottoDisk's `adr/`, so a decision that
spans packages reads the same way in each repo.

## Why these exist alongside the other docs

Keep the boundary sharp — each of these answers a question the others do not:

| Doc | Answers | Tense |
|---|---|---|
| roxygen blocks in `R/` | What a function does and what its arguments mean | present, per-function |
| `vignettes/` | How a workflow hangs together end to end | present, narrative |
| `NEWS.md` | What changed, per release | past, per-version |
| `adr/` (here) | Why we chose this over the alternatives, and when | past, immutable |

Practical test for where something belongs:

- "`filterData()` returns `list(feats_keep, cells_keep)`" → roxygen on the
  generic.
- "The subcellular workflow aggregates before overlapping" → a vignette.
- "Scaled expression became a `ScaledMatrix` in this release" → NEWS.md.
- "We reverted the forced log-norm write because snapshot adoption already
  covers it" → ADR.

Giotto has no `CLAUDE.md` / `AGENTS.md` today. If one is added it takes the
"invariants that hold right now, and where the code is" row, as in GiottoClass —
and this table should be updated rather than left to drift.

The overlap is intentional and one-directional: the roxygen docs and vignettes
state the *outcome* of an ADR without rehearsing the argument; the ADR is where
the argument and the discarded options live. When they disagree, the code and
its docs win for current behaviour and the ADR wins for intent — and the
disagreement is itself a signal that a superseding ADR is owed.

## Branch note

This directory was started on `docs/adr`, branched from `gsource` — the
long-running integration branch where the disk-backed and verb-split arc is
staged before it flows into `suite_dev`. The records are not gsource-only
history; they describe decisions the merged package carries.

## Writing one

1. Copy `0000-template.md` to `NNNN-short-kebab-title.md`, taking the next free
   number. Numbers are record order, not decision order — an ADR backfilled
   today for a 2025 decision still takes the next number and carries the older
   date.
2. Fill it in. Keep it to a page; if it needs more, the extra belongs in a
   vignette and the ADR should link to it.
3. Add a row to the index below.

## Amending one

An accepted ADR is not edited, with two exceptions: its `Status` line, and links
added to it (`Superseded by`). Everything else changes by writing a new ADR that
supersedes it. A wrong ADR that got reversed is more useful than no record of
the reversal.

Statuses: **Proposed** · **Accepted** · **Superseded by NNNN** · **Reversed**
(tried, undone, nothing replaced it) · **Deprecated** (still true, no longer
load-bearing).

## Scope — what earns an ADR

Something an outsider (or you in six months) would otherwise change by accident.
Roughly: a decision that constrains future code, was contested or non-obvious,
or has a cost worth remembering.

Not: bug fixes, refactors that preserve behaviour, or performance tuning whose
only artifact is a changed default — those belong in NEWS.md.

## Finding these

The index below is the fallback. The primary path is a pointer **from the code
the decision constrains** — `adr/NNNN` in a comment at the site someone would
edit. That is where an ADR gets read: at the moment you are about to change the
thing, not while browsing.

**None of the records below carry code pointers yet** — they were written as a
batch ahead of the pointer pass. Adding them is owed; the sites are named in
each record's *References*.

## Index

| # | Title | Status | Date |
|---|---|---|---|
| [0001](0001-wrappers-keep-names-dispatch-moves-to-verbs.md) | Wrappers keep their names; the numeric core moves to the verb generics | Accepted | 2026-05-19 |
| [0002](0002-lazy-matrices-stay-lazy-until-snapshot.md) | The normalization chain stays lazy; the vault commits at snapshot | Accepted | 2026-05-03 |
| [0003](0003-readers-detect-paths-once-into-paths-slot.md) | Readers detect paths once into `@paths` and inject them into the init frame | Accepted | 2026-05-20 |
| [0004](0004-autopcaparam-sentinel.md) | `"auto"` is a sentinel param class, resolved by the substrate at dispatch | Accepted | 2026-07-09 |

## Backfill candidates

Decisions already argued out elsewhere in the repo or in commit messages, not
yet written up. Not a queue — write one when it next comes up in conversation,
so the ADR captures the argument while it is fresh.

- **`runUMAP` moved to `uwot::umap2`, with `n_epochs` and `init` scaled to
  dataset size.** Benchmarked retune rather than a taste change, and it moved
  `RcppHNSW` / `rnndescent` into Suggests. Worth a record if the defaults are
  ever questioned; otherwise NEWS.md carries it. (`406e08257`, `3c81d8945`,
  `c4f678be7`, `8972443ca`). The commit trail also notes that a `umapParam` via
  `reduceData` is still owed — that is roadmap, not ADR.
- **`use_parallel` is skipped for `IterableMatrix` in the `covGroups` /
  `covLoess` HVF paths.** Parallelizing over a BPCells-backed matrix costs more
  than it saves. Constrains anyone adding a parallel branch to an HVF method.
  (`40983e461`)
- **`var_p_resid` uses a sparse-aware row variance**, and the Giotto and
  GiottoDisk implementations of it currently disagree on what they return —
  calibrated Pearson-residual variance from raw counts vs `rowVars()` of
  log-norm. The divergence is unresolved; a record is owed once it is, and it
  should say which definition won. (`ebc069bae`, `7de812377`)
- **`processData(allMatrix, BiocSingularParam)` is a compatibility shim** left
  from the single-surface first pass (0001). Whether it stays is a decision no
  one has made explicitly. (`59c1ee43e`)
- **Leiden `n_iterations` default raised to 20.** Same category as the UMAP
  retune — NEWS.md unless the default is contested. (`8972443ca`)
- **`instructions()` replaced the `*GiottoInstructions` family.** A deprecation
  with a migration, not an architectural choice; listed only so the next reader
  does not go looking for an ADR. (`9934b4d42`)

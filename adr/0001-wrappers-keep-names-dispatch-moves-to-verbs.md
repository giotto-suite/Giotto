# 0001. Wrappers keep their names; the numeric core moves to the verb generics

- **Status:** Accepted
- **Date:** 2026-05-19
- **Supersedes:** —
- **Superseded by:** —

## Context

Giotto's public surface is a set of workflow functions users call by name —
`filterGiotto()`, `normalizeGiotto()`, `runPCA()`, `calculateHVF()`,
`addCellStatistics()` / `addFeatStatistics()`, `filterCombinations()`. Every
tutorial, script and downstream package calls them, so their names and
signatures are effectively frozen.

The expression matrix inside a `giotto` object, however, is no longer one thing:
it can be a `dgCMatrix`, a `DelayedArray`, a BPCells `IterableMatrix`, or a
GiottoDisk `parquetExprStore` that must be streamed rather than loaded. The
numeric core of each workflow function needs a different implementation per
representation, and the streaming implementations live in GiottoDisk — which
depends on Giotto, not the other way round. Giotto cannot branch on a class it
does not import.

## Decision

The wrapper functions keep their names and remain the user entry point. They own
everything gobject-shaped — resolving `spat_unit` / `feat_type` defaults,
provenance, parameter history, setting results back into the object — and
delegate the numeric core to a verb generic with a `Param` object:

| wrapper | generic | param |
|---|---|---|
| `filterGiotto()`, `filterCombinations()` | `filterData()` | `filterParam` |
| `normalizeGiotto()` (library, log, scale) | `processData()` | `normParam`, `scaleParam` |
| `runPCA()` | `reduceData()` | `pcaParam` family |
| `calculateHVF()` | `analyzeData()` | HVF `analyzeParam`s |
| `addCellStatistics()` / `addFeatStatistics()` | `analyzeData()` | qc `analyzeParam`s |

Giotto defines the `Param` families and the default methods on `allMatrix`;
GiottoDisk registers methods for its store classes against the same generics.

This shape arrived in two passes. The first (2026-05-03) routed everything
through `processData()` — one generic, with the `Param` class carrying the
distinction. On 2026-05-19 the routings were moved onto the verb whose return
contract they actually match, following the split recorded in GiottoClass ADR
0003. One compatibility method remains from the first pass:
`processData(allMatrix, BiocSingularParam)`.

## Consequences

- A new expression backend is added by registering methods, with no edit to
  Giotto. That is the property the whole arrangement exists to buy.
- The numeric core becomes testable on a bare matrix, without a `giotto` object.
- The wrappers are thin but not trivial, and the boundary is a rule rather than a
  suggestion: gobject semantics in the wrapper, arithmetic in the method. A
  wrapper that computes anything is a bug in this design.
- Each `Param` family duplicates its wrapper's arguments and must be kept in
  sync by hand; adding a knob means touching both.
- The intermediate contracts are now API. `filterData()` returning
  `list(feats_keep, cells_keep)` and `reduceData()` returning
  `list(u, d, v, sdev, ...)` are relied on by methods in another package.
- Debugging costs one hop: the stack goes wrapper → generic → method, and with
  `autoPcaParam` (0004) two.
- The 05-03 → 05-19 move is the measured cost of choosing the wrong verb: wide
  and mechanical, touching every routing site and its tests. Route new work by
  return contract the first time.

## Alternatives considered

- **Dispatch on `giotto` itself** — the container class is identical regardless
  of what matrix it holds, so there is nothing to dispatch on. The varying thing
  is inside the object.
- **`if (inherits(x, "parquetExprStore"))` branches in the wrappers** — would
  make Giotto import GiottoDisk, inverting the dependency.
- **Keep the single `processData()` surface** — tried for two weeks and
  abandoned: the return type became `Param`-dependent, so no consumer could be
  written against the generic. See GiottoClass ADR 0003.
- **Have users call the generics directly** — breaks every existing script and
  tutorial, and moves gobject bookkeeping into user code.

## References

- `R/filter.R` — `filterGiotto()`, `filterCombinations()` → `filterData()`
- `R/normalize.R` — `normalizeGiotto()` → `processData()`; `normParam` /
  `scaleParam` / `adjustParam` / `threshParam` families
- `R/dimension_reduction.R` — `runPCA()` → `reduceData()`
- `R/variable_genes.R` — `calculateHVF()` → `analyzeData()` + its methods
- `R/auxiliary_giotto.R` — `addCellStatistics()` / `addFeatStatistics()`
- First pass: `a5e36a7c4`, `4b2cca8db`, `58605fe2e`, `b3907b631`, `d8494ee8a`
- Verb split: `3b2c3c655`, `8b65ca669`, `ffaef3155`, `df326d613`, `5f0a05114`,
  `75a4c7e1b`, `59c1ee43e`

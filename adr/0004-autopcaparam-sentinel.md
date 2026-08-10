# 0004. `"auto"` is a sentinel param class, resolved by the substrate at dispatch

- **Status:** Accepted
- **Date:** 2026-07-09
- **Supersedes:** —
- **Superseded by:** —

## Context

The right PCA algorithm depends on where the matrix lives. In-memory matrices
want IRLBA. A streaming `parquetExprStore` wants Halko-style randomized SVD or a
gram-eigen decomposition, chosen by GiottoDisk's own heuristics — and it must
ignore `scale`, since scaling densifies and breaks the O(N·k) streaming budget.

`runPCA(method = )` is a user-facing argument, so `"auto"` had to mean something.
But the knowledge needed to resolve it lives in the substrate's package, which is
downstream of Giotto (0001). Resolving `"auto"` inside `pcaParam()` is also too
early: the param is constructed before the generic sees `x`, so at construction
time there is nothing to inspect.

## Decision

`pcaParam("auto", ...)` returns an `autoPcaParam` — a `pcaParam` subclass that
names no algorithm. The routing is a method, not a switch:

```r
setMethod("reduceData", signature("allMatrix", "autoPcaParam"), ...)  # → irlba
setMethod("reduceData", signature("ANY", "autoPcaParam"), ...)        # → warn + irlba
# GiottoDisk registers reduceData(parquetExprStore, autoPcaParam) → random / gram-eigen
```

Each method strips `method` and `dry_run` from the param's knob list, rebuilds a
concrete flavor with `do.call(pcaParam, knobs)`, and re-dispatches on it.
`dry_run = TRUE` returns the resolved concrete param instead of running the
decomposition, and is set **only** on the auto flavor — a concrete param always
executes, so carrying the field there would make a resolved param ambiguous to
inspect.

## Consequences

- A substrate encodes its routing by registering one method. Adding a backend
  never edits a resolution table in Giotto.
- Dispatch happens twice: `reduceData(x, auto)` → `reduceData(x, concrete)`. Stack
  traces and profiles carry the extra frame, on top of the wrapper hop from 0001.
- The `ANY` fallback warns and uses IRLBA rather than failing. This is deliberate
  — a cryptic "unable to find an inherited method" is worse than a working
  fallback that names the missing method. The cost is that a substrate shipping
  without its routing method looks fine to anyone who suppresses warnings, and
  will quietly attempt an in-memory algorithm on data that may not fit.
- Knobs are forwarded by `do.call()`, so a new `pcaParam()` argument is carried
  through auto resolution automatically — but a new argument named `method` or
  `dry_run` would be silently dropped by the strip step.
- Revisit if a substrate needs to resolve on something other than `class(x)`
  (matrix dimensions, available memory); the current shape can express that
  inside a method, but only for classes that have one.

## Alternatives considered

- **`switch(class(x), ...)` inside `runPCA()`** — requires Giotto to know
  GiottoDisk's classes, inverting the dependency.
- **Resolve `"auto"` when the param is constructed** — `pcaParam()` never sees
  `x`, and params are routinely built before the data is in hand.
- **Resolve inside each concrete `reduceData()` method** — every method would
  need the full routing table, and adding a backend would mean editing all of
  them.
- **Drop `"auto"` and require an explicit method** — makes the user responsible
  for knowing which substrate their expression matrix is on, which is what the
  store classes exist to hide.

## References

- `R/pca-param.R` — `autoPcaParam` class, `pcaParam()` factory, the two
  `reduceData(*, autoPcaParam)` methods, `dry_run` handling
- `R/dimension_reduction.R` — `runPCA()` building the param
- Commits `438b94726` (sentinel), `703030dae` (`ANY` fallback + warning),
  `518a82a8e` (`do.call(pcaParam, knobs)` refactor)

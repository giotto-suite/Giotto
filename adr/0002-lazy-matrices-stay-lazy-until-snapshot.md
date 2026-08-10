# 0002. The normalization chain stays lazy; the vault commits at snapshot

- **Status:** Accepted
- **Date:** 2026-05-03
- **Supersedes:** —
- **Superseded by:** —

## Context

On a vault-backed project the normalization chain runs over a lazy
representation: a BPCells `IterableMatrix` reading from disk, or a
`parquetExprStore`. The chain is library normalization → `log1p` → centering and
scaling, and the steps differ in one way that matters:

- `log1p` **preserves sparsity**. Materializing its result is affordable.
- centering **densifies**. Materializing it is not.

`standardise_flex()` already encodes that asymmetry: it returns a
`ScaledMatrix` (delayed, centering applied at read) for `matrix` / `Matrix` /
`DelayedMatrix`, and a BPCells `add_rows()` / `multiply_rows()` chain for an
`IterableMatrix`. Neither ever builds a dense centered matrix.

That left the log-norm step. On 2026-05-01, as part of wiring GiottoDisk's
artifact lifecycle, `normalizeGiotto()` was changed to
`setGiotto(..., norm_expr, write = TRUE)` — forcing the lazy chain to
materialize into the vault at normalize time, on the reasoning that `log1p` is
sparsity-safe so the write is cheap and gives the project a real artifact.

## Decision

Reverted, on 2026-05-03. No step in `normalizeGiotto()` forces a write; both
`norm_expr` and `norm_scaled_expr` are set with the default `write = FALSE`.

Materialization has exactly two points, and normalize time is not one of them:

1. **Set time** — format conversion only. A setter converts an *in-memory*
   payload to a vault artifact (GiottoClass ADR 0002); anything already
   disk-backed passes through untouched.
2. **Snapshot time** — vault commitment. `snapshotSave()`'s
   `.ss_gdsrc_register_external_expr()` walks the expression slot and
   `sourceAdopt()`s any `IterableMatrix` not already in the vault.

The forced write was redundant with (2): the chain the write would have
materialized is adopted at snapshot anyway.

**The revert commit records no rationale** (`90dd884c0`, and a later
side-branch removal mislabelled with the original message, `7f3cd8b81`). The
above is what the surrounding code expresses, not a motive recovered from the
history. Treat the reasoning as reconstructed.

## Consequences

- A session that normalizes and never snapshots holds a lazy chain over vault
  inputs. Nothing is lost — but every read re-runs the chain, so repeated access
  to `"normalized"` pays for the arithmetic each time.
- `"scaled"` is a delayed matrix in the object. Any consumer must accept that:
  it exposes `dim()` but not a `dgeMatrix`-style `@Dim` slot, which is what the
  merFISH test assertions had to be changed to reflect.
- **Any new densifying step must not materialize into the gobject.** That is the
  standing constraint this record exists to protect.
- Revisit if chains get deep enough that re-read cost dominates. The fix then is
  an opt-in checkpoint at set time — `write = TRUE` at a chosen step — not an
  unconditional write, since the argument already exists for exactly this.

## Alternatives considered

- **Force the write after log-norm** — tried (2026-05-01) and reverted
  (2026-05-03). Redundant with snapshot adoption, and it puts I/O inside
  `normalizeGiotto()` where the caller cannot see or skip its cost.
- **Materialize the scaled matrix too** — densifies. This is the entire reason
  `ScaledMatrix` and the BPCells row-op chain are used instead of a plain scale.
- **Materialize nothing, ever, including at snapshot** — a snapshot that merely
  *references* matrices living in `tempdir()` loses them at session end. That
  was a real bug, and snapshot-time adoption is its fix.
- **Materialize at set time for every disk-backed payload** — re-writes
  artifacts already in the vault and duplicates them on disk.

## References

- `R/normalize.R` — `normalizeGiotto()`, the two `setGiotto()` calls at the end
- `GiottoClass::standardise_flex()` (`R/flex_functions.R`) — `ScaledMatrix` /
  BPCells lazy branches
- `GiottoDisk/R/methods-snapshotSave.R` —
  `.ss_gdsrc_register_external_expr()`, `sourceAdopt()`
- GiottoClass `adr/0002` — setter write-through and the `write =` argument
- Commits `3945c5da5` (added), `90dd884c0` (reverted), `7f3cd8b81` (later
  removal, misleading message), `1a738b419` (test assertions follow)

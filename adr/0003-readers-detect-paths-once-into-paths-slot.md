# 0003. Readers detect paths once into `@paths` and inject them into the init frame

- **Status:** Accepted
- **Date:** 2026-05-20
- **Supersedes:** —
- **Superseded by:** —

## Context

The platform readers — `CosmxReader`, `XeniumReader`, `VisiumHDReader`,
`StereoSeqReader` — exist to find files by convention inside a vendor output
directory so the user does not have to name each one. Discovery happened inline
in `initialize()`, and the closures the reader builds there
(`@calls$load_image`, `@calls$load_mask`, the `gobject` function) captured the
discovered paths as default arguments.

GiottoDisk subclasses these readers (`*DiskReader`). A subclass `initialize()`
runs in its own frame and needs the same discovered paths to build its own
closures. With detection inline, the subclass frame had two options: re-run
discovery, or reach back into variables it does not have.

Re-running is what happened, and it was broken in a way inline detection makes
easy: `StereoSeqReader`'s init called `.stereoseq_find_mask(p)` with a `p` that
did not exist in that frame. It never errored only because the mask path was
normally supplied explicitly — the auto-detect branch was dead code that would
have failed the moment it ran.

## Decision

Each reader class gains `paths = "list"` (prototype `list()`). `initialize()`
calls a single `.<platform>_detect_paths()` helper, assigns the result to
`.Object@paths`, and then injects the names into the init frame:

```r
.Object@paths <- .cosmx_detect_paths(.Object@cosmx_dir)
list2env(.Object@paths, envir = environment())
```

Both the base init body and any subclass init frame reference the paths as plain
variable names (`gef_path`, `bin1_gef_path`, `image_dir`, `mask_path`, …) via the
same convention. Detection runs once, in one place, per reader.

## Consequences

- Discovery happens once per reader construction and its result is inspectable on
  the object rather than trapped in a closure's default argument.
- **The names returned by `.<platform>_detect_paths()` are now API for
  subclasses.** Renaming a key silently breaks a `*DiskReader` init frame in
  another package, and the failure surfaces as object-not-found at
  construction, not at the rename.
- `list2env()` makes the bindings implicit: reading the init body, nothing at the
  use site says where `gef_path` came from. That is the readability cost paid to
  keep the closure default arguments short.
- A reader that needs a path not in the detect helper's return list must add it
  there rather than calling a finder inline — inline finders are exactly what
  produced the stereoseq bug.
- Revisit if the injected set grows large enough that explicit
  `.Object@paths$gef_path` at each use site becomes the clearer read.

## Alternatives considered

- **Re-run detection in each frame that needs it** — duplicated directory I/O,
  and the stereoseq `.stereoseq_find_mask(p)` bug is the concrete failure mode.
- **Explicit `.Object@paths$name` at every use site** — clearer provenance, but
  the paths are consumed mainly as closure default arguments, so every signature
  would carry a slot expression. Considered and traded away for readable
  closures.
- **Pass paths in as constructor arguments** — moves discovery back to the
  caller, which is the problem the readers exist to solve.
- **An environment instead of a list slot** — reference semantics would break
  the value semantics that let a reader be copied and serialized.

## References

- `R/convenience_cosmx.R` — `paths` slot, `.cosmx_detect_paths()`
- `R/convenience_xenium.R` — `.xenium_detect_paths()`
- `R/convenience_visiumHD.R` — `.visiumhd_detect_paths()`
- `R/convenience_stereoseq.r` — `.stereoseq_detect_paths()`, which absorbed
  `.stereoseq_find_mask()`
- Commits `8957d3975` (slot + helpers), `e34508a1e` (stereoseq mask path fix)

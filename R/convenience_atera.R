# Atera reader.
#
# Atera output is currently byte-for-byte the same layout as Xenium output, so
# `AteraReader` subclasses `XeniumReader` and overrides nothing. The point of
# the subclass is that it gives Atera its own public entry points and its own
# place to hang divergent behaviour later: when the formats split, the change is
# an override on `AteraReader`, with no effect on Xenium callers.
#
# `@include` is required, not cosmetic -- Collate is alphabetical, so without it
# `convenience_atera.R` is sourced before `convenience_xenium.R` and
# `contains = "XeniumReader"` fails at build time.

#' @include convenience_xenium.R
NULL

#' @title Atera reader
#' @name AteraReader-class
#' @description
#' Reader for Atera output directories. Inherits the `XeniumReader`
#' implementation, since the two output layouts are presently identical.
#' @keywords internal
setClass("AteraReader",
    contains = "XeniumReader",
    # only override: path-detection messages should name Atera
    prototype = list(platform = "Atera")
)

#' @export
`.DollarNames.AteraReader` <- `.DollarNames.XeniumReader`

#' @name importAtera
#' @title Import an Atera dataset
#' @description
#' Creates an `AteraReader` instance with reader functions for converting
#' individual pieces of Atera data into Giotto-compatible representations.
#'
#' The Atera output layout currently matches Xenium's, so this reader inherits
#' the Xenium implementation. It exists as a separate entry point so that the
#' two can diverge without affecting Xenium users.
#' @param atera_dir Atera output directory
#' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included
#' as a subcellular transcript detection (default = 20). Only applies when
#' transcript-level data is present.
#' @param backend (optional) a `gsource`-inheriting project backend. When
#' provided, creates the `giotto` object as a managed on-disk project.
#' @returns `AteraReader` object, or `AteraDiskReader` when `backend` is set
#' @examples
#' \dontrun{
#' x <- importAtera("path/to/atera_outs")
#' g <- x$create_gobject(load_transcripts = FALSE)
#' }
#' @export
importAtera <- function(atera_dir = NULL, qv_threshold = 20, backend = NULL) {
    if (!is.null(backend)) {
        package_check(
            "GiottoDisk",
            repository = "github:giotto-suite/GiottoDisk"
        )
        return(GiottoDisk::importAteraDisk(
            atera_dir = atera_dir,
            backend = backend,
            qv_threshold = qv_threshold
        ))
    }

    a <- list(Class = "AteraReader")
    # the inherited slot is named `xenium_dir`; only the public argument is
    # renamed, so nothing about the parent implementation has to change
    if (!is.null(atera_dir)) {
        a$xenium_dir <- atera_dir
    }
    a$qv <- qv_threshold

    do.call(new, args = a)
}

#' @name createGiottoAteraObject
#' @title Create a Giotto object from Atera output
#' @description
#' Convenience wrapper around [importAtera()] that builds a `giotto` object in
#' one call. Parameters and defaults match [createGiottoXeniumObject()], since
#' the output layouts are presently identical.
#'
#' Note that `load_transcripts` defaults to `TRUE` to match the Xenium
#' behaviour. Atera exports that ship only the aggregated
#' `cell_feature_matrix` need `load_transcripts = FALSE`.
#' @inheritParams createGiottoXeniumObject
#' @param atera_dir Atera output directory
#' @returns a `giotto` object
#' @examples
#' \dontrun{
#' # aggregated-matrix export (no transcripts.parquet)
#' g <- createGiottoAteraObject(
#'     atera_dir = "path/to/atera_outs",
#'     load_transcripts = FALSE
#' )
#' }
#' @export
createGiottoAteraObject <- function(
        atera_dir,
        transcript_path = NULL,
        bounds_path = list(
            cell = "cell",
            nucleus = "nucleus"
        ),
        gene_panel_json_path = NULL,
        expression_path = NULL,
        cell_metadata_path = NULL,
        feat_type = c(
            "rna",
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword",
            "GenomicControl"
        ),
        split_keyword = list(
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword",
            "GenomicControl"
        ),
        qv_threshold = 20,
        load_images = "focus",
        load_aligned_images = NULL,
        load_transcripts = TRUE,
        load_expression = TRUE,
        load_cellmeta = FALSE,
        backend = NULL,
        instructions = NULL,
        verbose = NULL) {
    x <- importAtera(atera_dir, backend = backend)
    x$qv <- qv_threshold

    a <- list(
        load_bounds = bounds_path,
        feat_type = feat_type,
        split_keyword = split_keyword,
        load_images = load_images,
        load_aligned_images = load_aligned_images,
        load_transcripts = load_transcripts,
        load_expression = load_expression,
        load_cellmeta = load_cellmeta,
        instructions = instructions,
        verbose = verbose
    )

    if (!is.null(transcript_path)) a$transcript_path <- transcript_path
    if (!is.null(gene_panel_json_path)) {
        a$gene_panel_json_path <- gene_panel_json_path
    }
    if (!is.null(expression_path)) a$expression_path <- expression_path
    if (!is.null(cell_metadata_path)) a$metadata_path <- cell_metadata_path

    g <- do.call(x$create_gobject, args = a)
    return(g)
}

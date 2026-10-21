# Atera reader.
#
# Atera v1 output ships zarr archives only and is read on disk, through
# GiottoDisk's own Atera reader (`backend =` -> GiottoDisk::importAteraDisk()).
# There is no in-memory Atera reader for now: without a backend, importAtera()
# stops. Atera preview data in the Xenium layout is read with the Xenium
# readers.

#' @name importAtera
#' @md
#' @title Import an Atera dataset
#' @description
#' Creates a reader with functions for converting individual pieces of Atera
#' data into Giotto-compatible representations.
#'
#' Atera v1 output ships zarr archives only and is read on disk: `backend` is
#' required and returns a `GiottoDisk::importAteraDisk()` reader, which writes
#' each data type to the on-disk project. There is no in-memory Atera reader;
#' without a `backend` this stops. Atera preview data in the Xenium layout
#' (parquet / h5 files) is read with [importXenium()] and
#' [createGiottoXeniumObject()].
#' @param atera_dir Atera output directory
#' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included
#' as a subcellular transcript detection (default = 20).
#' @param backend a `gsource`-inheriting project backend (typically
#' `GiottoDisk::sourceCreate()`), or a path to one. Required.
#' @returns an `AteraDiskReader` object
#' @examples
#' \dontrun{
#' x <- importAtera("path/to/atera_outs", backend = "path/to/project")
#' g <- x$create_gobject()
#' }
#' @export
importAtera <- function(atera_dir = NULL, qv_threshold = 20, backend = NULL) {
    if (is.null(backend)) {
        stop(wrap_txt(
            "Atera output is read on disk: pass `backend =` (a GiottoDisk
            project, e.g. from GiottoDisk::sourceCreate(), or a path for one).
            Atera preview data in the Xenium layout reads with importXenium() /
            createGiottoXeniumObject()."
        ), call. = FALSE)
    }
    package_check(
        "GiottoDisk",
        repository = "github:giotto-suite/GiottoDisk"
    )
    GiottoDisk::importAteraDisk(
        atera_dir = atera_dir,
        backend = backend,
        qv_threshold = qv_threshold
    )
}

#' @name createGiottoAteraObject
#' @md
#' @title Create a Giotto object from Atera output
#' @description
#' Convenience wrapper around [importAtera()] that builds a `giotto` object in
#' one call. Atera output is read on disk, so `backend` is required (see
#' [importAtera()]): transcripts, polygons
#' and expression are written to the on-disk project, and `load_images =
#' "focus"` loads the `morphology_2d` images. Parameters match
#' [createGiottoXeniumObject()]; `gene_panel_json_path` does not apply, since
#' Atera ships no gene panel. Atera preview data in the Xenium layout is read
#' with [createGiottoXeniumObject()].
#' @inheritParams createGiottoXeniumObject
#' @param atera_dir Atera output directory
#' @returns a `giotto` object
#' @examples
#' \dontrun{
#' g <- createGiottoAteraObject(
#'     atera_dir = "path/to/atera_outs",
#'     backend = "path/to/project"
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
    x <- importAtera(atera_dir, qv_threshold = qv_threshold, backend = backend)

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

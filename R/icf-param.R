# icfParam ####
#
# Interaction-changed features as an `analyzeData` verb, following the motif and
# deconvolution layers (spatial_motif.R, deconv-param.R): a VIRTUAL parent under
# `analyzeParam`, one concrete class per engine, a factory, and `analyzeData`
# methods that dispatch on the *expression object*. A disk-backed store therefore
# arrives intact and can route to its own streaming method (GiottoDisk) with no
# change here.
#
# The algorithm is not in this package: it is the edgeicg package
# (github: drieslab/edgeicg), which Giotto suggests, as the motif layer uses
# smotif. Giotto owns the param class, the verb and the gobject plumbing; the
# carrier method below is the only code here that calls the engine. The classic
# findInteractionChangedFeats() algorithm is unchanged and stays the default.

#' @title Interaction-changed feature params
#' @name icf_param
#' @description
#' Parameter objects for interaction-changed feature (ICF) analysis through
#' [analyzeData()]. Use them through [findInteractionChangedFeats()] with
#' `method = "edgeicg"`, or pass one to [analyzeData()] with a `giotto` object.
#'
#' `"edgeicg"` (package \pkg{edgeicg}) tests, for every receiver cell type,
#' sender cell type and feature, whether the receiver's expression changes with
#' its exposure to the sender at three scales: touching cells (`contact`),
#' short-range neighbours (`paracrine`) and the regional neighbourhood
#' (`niche`). Transcripts leaked from touching neighbours are modelled and
#' removed per cell before testing; niche hits are calibrated with a spatial
#' (toroidal-shift) null and attributed to senders.
#' @param method character. Only `"edgeicg"` for now.
#' @param ... method-specific parameters; see Details.
#' @details
#' edgeicg:
#' \describe{
#'   \item{`expression_values`}{raw counts to use (default `"raw"`). The model
#'     needs counts, not normalized values.}
#'   \item{`use_polygons`}{logical. Build the spillover graph (touching cells)
#'     from the cell polygons of the spatial unit when they exist (default
#'     `TRUE`); otherwise k nearest neighbours.}
#'   \item{`block_genes`}{integer. features per block; bounds memory
#'     (default 500).}
#'   \item{`workers`}{integer. `mirai` daemons (default 0, sequential).
#'     Results are identical.}
#'   \item{`window_cells`}{cells per window. Above this, cells are visited in
#'     spatially ordered windows and every statistic is summed over windows,
#'     so memory does not grow with the dataset. Default: one window for
#'     in-memory counts; for a GiottoDisk store, one window per worker (at
#'     most 50,000 cells), with the store baked once in window order and read
#'     in contiguous cell ranges.}
#'   \item{`backend`}{`"auto"` (default: the Rust kernels of \pkg{edgeicgrs}
#'     when installed), `"rust"` or `"r"`.}
#'   \item{`calibrate`}{logical. spatial calibration of the niche scale
#'     (default `TRUE`).}
#'   \item{`min_shift`}{numeric. minimum calibration shift in microns
#'     (default 450).}
#'   \item{`report_q`}{numeric. entries with q below this are returned in
#'     `ICFscores` (default 0.1); the full fit is kept in `$edgeicg`.}
#'   \item{`min_lfc`}{numeric. entries must also change by at least this
#'     log2 fold change between receivers at the 10th and 90th percentile of
#'     their exposure (default 0.25; 0 = significance only). At whole-section
#'     scale significance alone also picks up small confounded associations
#'     (zonation, sub-states within a cell type).}
#'   \item{`spatial_covariates`}{logical. adjust the niche screen for smooth
#'     spatial trends (`edgeicg::spatial_basis()`), which removes most
#'     large-scale zonation from niche hits; slower (default `FALSE`).}
#'   \item{`spatial_df`, `covariate_scope`}{spline degrees of freedom per axis
#'     of the spatial basis (default 6) and where it enters: `"niche"`
#'     (default) or `"all"` tests.}
#'   \item{other}{passed to `edgeicg::edgeicg()` (e.g. `min_cells`,
#'     `min_det`, `min_exposed`, `calib_cond`).}
#' }
#' @md
#' @returns an [icfParam-class]-inheriting object
#' @examples
#' p <- icfParam("edgeicg", workers = 4)
#' p$workers
#' @seealso [findInteractionChangedFeats()], [analyze_param]
NULL

#' @rdname icf_param
#' @exportClass icfParam
setClass("icfParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname icf_param
#' @exportClass edgeIcfParam
setClass("edgeIcfParam", contains = "icfParam")

# param factory ####
#' @rdname icf_param
#' @export
icfParam <- function(method = "edgeicg", ...) {
    method <- match.arg(tolower(method), "edgeicg")
    switch(method,
        "edgeicg" = .icf_param_edgeicg(...)
    )
}

#' @keywords internal
#' @noRd
.icf_param_edgeicg <- function(...) {
    p <- new("edgeIcfParam", param = list(...))
    p$expression_values <- p$expression_values %null% "raw"
    p$use_polygons <- !isFALSE(p$use_polygons)
    p$block_genes <- as.integer(p$block_genes %null% 500L)
    p$workers <- as.integer(p$workers %null% 0L)
    p$calibrate <- !isFALSE(p$calibrate)
    p$min_shift <- as.numeric(p$min_shift %null% 450)
    p$report_q <- as.numeric(p$report_q %null% 0.1)
    p$min_lfc <- as.numeric(p$min_lfc %null% 0.25)
    p$spatial_covariates <- isTRUE(p$spatial_covariates)
    p$spatial_df <- as.integer(p$spatial_df %null% 6L)
    p$covariate_scope <- p$covariate_scope %null% "niche"
    p$verbose <- isTRUE(p$verbose)
    if (p$block_genes < 1L) {
        stop("[icfParam] `block_genes` must be a positive integer.",
            call. = FALSE)
    }
    p
}

# analyzeData(<expression>, edgeIcfParam) ####
#' @title edgeicg interaction-changed features
#' @name icf_edgeicg
#' @description
#' The [analyzeData()] method for an [edgeIcfParam-class] on an expression
#' matrix. The computation is done by the \pkg{edgeicg} package, which is
#' suggested rather than imported: install it with
#' `remotes::install_github("drieslab/edgeicg")`.
#' @param x raw counts, features by cells: a `matrix` or `Matrix`, or a
#'   GiottoDisk expression store, which edgeicg reads in windows without
#'   loading it.
#' @param param an [edgeIcfParam-class].
#' @param coords numeric matrix of cell centroids (x, y), row names = cell IDs.
#' @param cell_type cell type labels named by cell ID.
#' @param lib library size per cell over all features, named by cell ID.
#' @param spill_edges optional touching-cell edges (`recv`, `send`, `w`) with
#'   `recv` / `send` indexing `colnames(x)`.
#' @param feats optional subset of features to test.
#' @param ... unused. Declared because [analyzeData()] carries it.
#' @md
#' @returns an `EdgeICG` fit (see `edgeicg::edgeicg()`).
#' @aliases analyzeData,ANY,edgeIcfParam-method
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "edgeIcfParam"),
    function(x, param, ..., coords, cell_type, lib, spill_edges = NULL,
    feats = NULL) {
        package_check("edgeicg", repository = "github:drieslab/edgeicg")
        ids <- colnames(x)
        if (is.null(ids)) .gstop("expression needs cell IDs as column names")
        missing <- setdiff(ids, names(cell_type))
        if (length(missing)) {
            .gstop(length(missing), "cell(s) have no cell type label")
        }
        if (!is.null(feats)) {
            keep <- which(rownames(x) %in% feats)
            # a matrix is subset in memory; a store keeps a lazy view
            x <- if (is.matrix(x) || inherits(x, "Matrix")) {
                x[keep, , drop = FALSE]
            } else {
                x[keep, ]
            }
        }
        a <- param@param
        extra <- a[setdiff(names(a), c("expression_values", "use_polygons",
            "report_q", "min_lfc", "spatial_covariates", "spatial_df",
            "covariate_scope", "block_genes", "workers", "calibrate",
            "min_shift", "verbose"))]
        if (isTRUE(param$spatial_covariates)) {
            extra$covariates <- edgeicg::spatial_basis(
                coords[ids, , drop = FALSE], df = param$spatial_df
            )
            extra$covariate_scope <- param$covariate_scope
        }
        do.call(edgeicg::edgeicg, c(list(
            counts = x,
            coords = coords[ids, , drop = FALSE],
            types = as.character(cell_type[ids]),
            lib = as.numeric(lib[ids]),
            spill_edges = spill_edges,
            block_genes = param$block_genes,
            workers = param$workers,
            calibrate = param$calibrate,
            min_shift = param$min_shift,
            verbose = param$verbose
        ), extra))
    }
)

# analyzeData(<giotto>, icfParam) ####
# Dispatches on the VIRTUAL class, so it serves every engine. Owns everything
# gobject-shaped: defaults, metadata, alignment by cell ID, polygons.
#' @rdname icf_edgeicg
#' @param spat_unit,feat_type spatial unit and feature type
#' @param cluster_column cell metadata column holding cell type labels
#' @param selected_feats optional subset of features to test
setMethod("analyzeData",
    signature(x = "giotto", param = "icfParam"),
    function(x, param, spat_unit = NULL, feat_type = NULL,
    cluster_column = NULL, selected_feats = NULL, ...) {
        spat_unit <- set_default_spat_unit(gobject = x, spat_unit = spat_unit)
        feat_type <- set_default_feat_type(
            gobject = x, spat_unit = spat_unit, feat_type = feat_type
        )
        if (is.null(cluster_column)) .gstop("cluster_column is required")
        expr <- getExpression(x,
            spat_unit = spat_unit, feat_type = feat_type,
            values = param$expression_values, output = "exprObj"
        )
        ids <- colnames(expr[])
        meta <- getCellMetadata(x,
            spat_unit = spat_unit, feat_type = feat_type,
            output = "data.table", copy_obj = TRUE
        )
        if (!cluster_column %in% colnames(meta)) {
            .gstop(sprintf(
                "cluster_column '%s' is not a cell metadata column",
                cluster_column
            ))
        }
        cell_type <- stats::setNames(
            as.character(meta[[cluster_column]]), meta[["cell_ID"]]
        )
        sl <- getSpatialLocations(x,
            spat_unit = spat_unit, output = "data.table", copy_obj = TRUE
        )
        coords <- as.matrix(sl[, c("sdimx", "sdimy")])
        rownames(coords) <- sl[["cell_ID"]]
        if (anyNA(match(ids, rownames(coords)))) {
            .gstop("some expression cells have no spatial location")
        }
        # library size over all features, through the same dispatch as
        # addCellStatistics(), so streaming backends supply their own
        cs <- analyzeData(expr[], analyzeParam("cell_stats"))
        lib <- stats::setNames(cs[["total_expr"]], cs[["cells"]])
        spill <- if (isTRUE(param$use_polygons)) {
            .icf_touching_edges(x, spat_unit, ids)
        }
        fit <- analyzeData(expr[], param,
            coords = coords, cell_type = cell_type, lib = lib,
            spill_edges = spill, feats = selected_feats, ...
        )
        .icf_from_edgeicg(fit,
            values = param$expression_values, cluster = cluster_column,
            spill = if (is.null(spill)) "knn" else "polygons",
            report_q = param$report_q, min_lfc = param$min_lfc
        )
    }
)

# helpers ####
# Touching-cell edges from the spatial unit's polygons (NULL when there are
# none). Vertices of all parts, outer rings only.
#' @keywords internal
#' @noRd
.icf_touching_edges <- function(gobject, spat_unit, ids) {
    gp <- tryCatch(
        getPolygonInfo(gobject,
            name = spat_unit,
            return_giottoPolygon = TRUE, verbose = FALSE
        ),
        error = function(e) NULL
    )
    if (is.null(gp)) return(NULL)
    sv <- gp[]
    # GiottoDisk-backed polygons: read the outlines once (vertices are needed)
    if (!inherits(sv, "SpatVector")) {
        sv <- GiottoDisk::storeRead(sv, output = "terra")
    }
    geo <- terra::geom(sv)
    geo <- geo[geo[, "hole"] == 0, , drop = FALSE]
    pid <- terra::values(sv)[["poly_ID"]][geo[, "geom"]]
    keep <- pid %in% ids
    if (!any(keep)) return(NULL)
    poly <- data.frame(
        cell = match(pid[keep], ids), x = geo[keep, "x"], y = geo[keep, "y"]
    )
    edgeicg::contact_from_polygons(poly, length(ids))
}

# EdgeICG fit -> icfObject. ICFscores keeps the classic columns where they
# have a meaning (feats, cell_type = receiver, int_cell_type = sender, diff,
# log2fc, p.value, p.adj) and adds scale, se, z, class and, for niche hits,
# the attributed sender.
#' @keywords internal
#' @noRd
.icf_from_edgeicg <- function(fit, values, cluster, spill, report_q,
    min_lfc = 0.25) {
    r <- edgeicg::edgeicg_results(fit, q_thr = report_q, min_lfc = min_lfc)
    if (nrow(r)) r <- edgeicg::edgeicg_classify(fit, r)
    # NSE variables
    receiver <- sender <- gene <- effect <- z <- q <- NULL
    scores <- data.table::data.table(
        feats = r$gene,
        cell_type = r$receiver,
        int_cell_type = r$sender,
        scale = r$scale,
        diff = r$effect,
        log2fc = r$lfc,
        se = r$se,
        z = r$z,
        p.value = 2 * stats::pnorm(-abs(r$z)),
        p.adj = r$q,
        class = r$class,
        regional = r$regional,
        top_sender = if ("top_sender" %in% names(r)) r$top_sender else NA,
        alpha = if ("alpha" %in% names(r)) r$alpha else NA_real_
    )
    if (nrow(scores)) {
        scores[, `:=`(
            unif_int = paste0(cell_type, "--", int_cell_type),
            spec_int = paste0(cell_type, "--", int_cell_type),
            type_int = "hetero"
        )]
    }
    structure(
        .Data = list(
            ICFscores = scores,
            Giotto_info = list(
                "values" = values,
                "cluster" = cluster,
                "spatial network" = paste0("edgeicg (spillover: ", spill, ")")
            ),
            test_info = list(
                "test" = "edgeicg",
                "p.adj" = "BH within each scale",
                "report_q" = report_q,
                "min_lfc" = min_lfc,
                "calibrated niche" = !is.null(fit$calibration)
            ),
            edgeicg = fit
        ),
        class = "icfObject"
    )
}

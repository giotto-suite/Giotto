# cccParam ####
#
# Cell-cell communication (ligand-receptor co-location) as an `analyzeData`
# verb, following the ICF layer (icf-param.R): a VIRTUAL parent under
# `analyzeParam`, one concrete class per engine, a factory, and `analyzeData`
# methods that dispatch on the *expression object*, so a GiottoDisk store
# arrives intact and is read per feature block.
#
# The algorithm is not in this package: it is the lrmap package
# (github: drieslab/lrmap), which Giotto suggests. Giotto owns the param class,
# the verb and the gobject plumbing. The classic spatCellCellcom() test
# (Delaunay neighbours, mean ligand + receptor against random cell positions)
# is unchanged and kept as a baseline (method = "classic"); lrmap is the
# default.

#' @title Cell-cell communication params
#' @name ccc_param
#' @description
#' Parameter objects for ligand-receptor co-location through [analyzeData()].
#' Use them through [spatCellCellcom()] with `method = "lrmap"`, or pass one to
#' [analyzeData()] with a `giotto` object.
#'
#' `"lrmap"` (package \pkg{lrmap}) scores, for every sender cell type,
#' receiver cell type and ligand-receptor pair, how much ligand from sender
#' cells surrounds receiver cells that express the receptor, at three scales:
#' touching cells (`contact`), short-range neighbours (`paracrine`) and the
#' regional neighbourhood (`niche`). Significance comes from toroidal shifts
#' of the sender layer. The default test asks whether ligand-high senders sit
#' nearer receptor-high receivers than the two cell types do on average, so
#' a call is specific to the pair rather than to the cell-type layout.
#' Transcripts leaked between touching cells are estimated and removed first.
#' @param method character. Only `"lrmap"` for now.
#' @param ... method-specific parameters; see Details.
#' @details
#' lrmap:
#' \describe{
#'   \item{`expression_values`}{raw counts to use (default `"raw"`).}
#'   \item{`scales`}{any of `"contact"`, `"paracrine"`, `"niche"` and
#'     `"expression"` (the expression-only reference: mean ligand x mean
#'     receptor against shuffled cell labels). Default: the three spatial
#'     scales.}
#'   \item{`leak`}{`TRUE` (default): estimate the leak between touching cells
#'     on all features with `lrmap::lrmap_leak()` and remove it; `FALSE`: no
#'     leak removal; or leak coefficients (an `lrmap_leak()` table or an
#'     edgeicg fit).}
#'   \item{`leak_pool`}{how pairs of cell types with little contact get their
#'     leak: `"median"` (default), `"q1"` or `"none"` (see
#'     `lrmap::lrmap_leak()`).}
#'   \item{`test`}{which null decides `p.adj`, `log2fc` and `PI` in the
#'     result: `"type"` (default; beyond the co-location of the two cell
#'     types), `"shuffle"` (within-type shuffle) or `"colocation"` (includes
#'     the co-location of the two cell types). All are kept as columns.}
#'   \item{`threads`}{threads for the Rust kernels of \pkg{lrmaprs} when
#'     installed (default 1).}
#'   \item{other}{passed to `lrmap::lrmap()`, e.g. `k`, `sigma_p`, `sigma_n`
#'     (scale definitions), `B`, `min_shift`, `min_det`, `min_cells`,
#'     `same_type`, `seed`.}
#' }
#' @md
#' @returns a [cccParam-class]-inheriting object
#' @examples
#' p <- cccParam("lrmap", sigma_p = 40, threads = 2)
#' p$sigma_p
#' @seealso [spatCellCellcom()], [analyze_param]
NULL

#' @rdname ccc_param
#' @exportClass cccParam
setClass("cccParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname ccc_param
#' @exportClass lrmapCccParam
setClass("lrmapCccParam", contains = "cccParam")

# param factory ####
#' @rdname ccc_param
#' @export
cccParam <- function(method = "lrmap", ...) {
    method <- match.arg(tolower(method), "lrmap")
    switch(method,
        "lrmap" = .ccc_param_lrmap(...)
    )
}

#' @keywords internal
#' @noRd
.ccc_param_lrmap <- function(...) {
    p <- new("lrmapCccParam", param = list(...))
    p$expression_values <- p$expression_values %null% "raw"
    p$scales <- p$scales %null% c("contact", "paracrine", "niche")
    p$leak <- p$leak %null% TRUE
    p$leak_pool <- p$leak_pool %null% "median"
    p$test <- p$test %null% "type"
    p$threads <- as.integer(p$threads %null% 1L)
    p$verbose <- isTRUE(p$verbose)
    bad <- setdiff(p$scales, c("contact", "paracrine", "niche", "expression"))
    if (length(bad)) {
        stop("[cccParam] unknown scales: ", paste(bad, collapse = ", "),
            call. = FALSE)
    }
    if (!p$test %in% c("type", "shuffle", "colocation")) {
        stop("[cccParam] `test` must be \"type\", \"shuffle\" or ",
            "\"colocation\".", call. = FALSE)
    }
    p
}

# analyzeData(<expression>, lrmapCccParam) ####
#' @title lrmap ligand-receptor co-location
#' @name ccc_lrmap
#' @description
#' The [analyzeData()] method for an [lrmapCccParam-class] on an expression
#' matrix. The computation is done by the \pkg{lrmap} package, which is
#' suggested rather than imported: install it with
#' `remotes::install_github("drieslab/lrmap")` (and optionally
#' `drieslab/lrmaprs` for its Rust kernels).
#' @param x raw counts, features by cells: a `matrix` or `Matrix`, or a
#'   GiottoDisk expression store, of which only the ligand and receptor
#'   features (and, for the leak estimate, feature blocks) are read.
#' @param param an [lrmapCccParam-class].
#' @param coords numeric matrix of cell centroids (x, y), row names = cell IDs.
#' @param cell_type cell type labels named by cell ID.
#' @param lib library size per cell over all features, named by cell ID.
#' @param lr ligand-receptor table: columns `ligand`, `receptor` and
#'   optionally `ligand_subunits`, `receptor_subunits` (genes joined by `_`).
#' @param ... unused. Declared because [analyzeData()] carries it.
#' @md
#' @returns an `LRmap` fit (see `lrmap::lrmap()`).
#' @aliases analyzeData,ANY,lrmapCccParam-method
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "lrmapCccParam"),
    function(x, param, ..., coords, cell_type, lib, lr) {
        package_check("lrmap", repository = "github:drieslab/lrmap")
        ids <- colnames(x)
        if (is.null(ids)) .gstop("expression needs cell IDs as column names")
        missing <- setdiff(ids, names(cell_type))
        if (length(missing)) {
            .gstop(length(missing), "cell(s) have no cell type label")
        }
        counts <- .ccc_counts(x, ids)
        types <- as.character(cell_type[ids])
        xy <- coords[ids, , drop = FALSE]
        lb <- as.numeric(lib[ids])
        leak <- param$leak
        if (isTRUE(leak)) {
            leak <- lrmap::lrmap_leak(counts, xy, types, lib = lb,
                genes = rownames(x), pool = param$leak_pool,
                verbose = param$verbose)
        } else if (isFALSE(leak)) {
            leak <- NULL
        }
        a <- param@param
        extra <- a[setdiff(names(a), c("expression_values", "scales", "leak",
            "leak_pool", "test", "threads", "verbose"))]
        do.call(lrmap::lrmap, c(list(
            counts = counts, coords = xy, types = types, lr = lr, lib = lb,
            scales = param$scales, leak = leak, leak_pool = param$leak_pool,
            threads = param$threads, verbose = param$verbose
        ), extra))
    }
)

# analyzeData(<giotto>, cccParam) ####
# Dispatches on the VIRTUAL class, so it serves every engine. Owns everything
# gobject-shaped: defaults, metadata, alignment by cell ID.
#' @rdname ccc_lrmap
#' @param spat_unit,feat_type spatial unit and feature type
#' @param cluster_column cell metadata column holding cell type labels
setMethod("analyzeData",
    signature(x = "giotto", param = "cccParam"),
    function(x, param, spat_unit = NULL, feat_type = NULL,
    cluster_column = NULL, lr = NULL, ...) {
        spat_unit <- set_default_spat_unit(gobject = x, spat_unit = spat_unit)
        feat_type <- set_default_feat_type(
            gobject = x, spat_unit = spat_unit, feat_type = feat_type
        )
        if (is.null(cluster_column)) .gstop("cluster_column is required")
        if (is.null(lr) || !nrow(lr)) .gstop("a ligand-receptor table is required")
        inp <- .ccc_inputs(x, spat_unit, feat_type, cluster_column,
            param$expression_values)
        fit <- analyzeData(inp$expr, param,
            coords = inp$coords, cell_type = inp$cell_type, lib = inp$lib,
            lr = lr, ...
        )
        .ccc_from_lrmap(fit, lr = lr, values = param$expression_values,
            cluster = cluster_column, test = param$test)
    }
)

# helpers ####
# Expression (matrix or store), cell types, coordinates and library size of a
# spatial unit, aligned by cell ID.
#' @keywords internal
#' @noRd
.ccc_inputs <- function(gobject, spat_unit, feat_type, cluster_column,
    values) {
    expr <- getExpression(gobject,
        spat_unit = spat_unit, feat_type = feat_type,
        values = values, output = "exprObj"
    )[]
    ids <- colnames(expr)
    meta <- getCellMetadata(gobject,
        spat_unit = spat_unit, feat_type = feat_type,
        output = "data.table", copy_obj = TRUE
    )
    if (!cluster_column %in% colnames(meta)) {
        .gstop(sprintf(
            "cluster_column '%s' is not a cell metadata column", cluster_column
        ))
    }
    cell_type <- stats::setNames(
        as.character(meta[[cluster_column]]), meta[["cell_ID"]]
    )
    sl <- getSpatialLocations(gobject,
        spat_unit = spat_unit, output = "data.table", copy_obj = TRUE
    )
    coords <- as.matrix(sl[, c("sdimx", "sdimy")])
    rownames(coords) <- sl[["cell_ID"]]
    if (anyNA(match(ids, rownames(coords)))) {
        .gstop("some expression cells have no spatial location")
    }
    # library size over all features, through the same dispatch as
    # addCellStatistics(), so streaming backends supply their own
    cs <- analyzeData(expr, analyzeParam("cell_stats"))
    lib <- stats::setNames(cs[["total_expr"]], cs[["cells"]])
    list(expr = expr, ids = ids, cell_type = cell_type, coords = coords,
        lib = lib)
}

# lrmap's counts input: the matrix itself, or for a GiottoDisk store a reader
# function(genes) returning those features for the cells `ids` (sparse)
#' @keywords internal
#' @noRd
.ccc_counts <- function(x, ids) {
    if (is.matrix(x) || inherits(x, "Matrix")) return(x)
    rn <- rownames(x)
    function(genes) {
        g <- intersect(genes, rn)
        m <- GiottoDisk::storeRead(x[g, ], output = "dgcmatrix",
            max_rows = Inf, max_cols = Inf)
        m[g, ids, drop = FALSE]
    }
}

# LRmap fit -> data.table with the spatCellCellcom() columns where they have a
# meaning (LR_comb, lig_cell_type, lig_expr, ligand, rec_cell_type, rec_expr,
# receptor, LR_expr, rand_expr, av_diff, log2fc, pvalue, p.adj, LR_cell_comb,
# PI, lig_nr, rec_nr) plus the lrmap columns. One row per pair, cell-type pair and scale.
# The fit, the LR table and the call settings are kept as attributes.
#' @keywords internal
#' @noRd
.ccc_from_lrmap <- function(fit, lr, values, cluster, test = "type") {
    tb <- data.table::as.data.table(fit$table)
    ex <- tb$scale == "expression"
    qq <- switch(test, type = tb$q_type, shuffle = tb$q_shuffle,
        colocation = tb$q)
    pp <- switch(test, type = tb$p_type, shuffle = tb$p_shuffle,
        colocation = tb$p)
    zz <- switch(test, type = tb$z_type, shuffle = tb$z_shuffle,
        colocation = tb$z)
    ff <- if (test == "colocation") tb$fold else tb$fold_type
    qq[ex] <- tb$q[ex]; pp[ex] <- tb$p[ex]; zz[ex] <- tb$z[ex]
    ff[ex] <- tb$fold[ex]
    lfc <- log2(pmax(ff, 1e-12))
    lowest <- suppressWarnings(min(qq[!is.na(qq) & qq > 0]))
    q_eff <- ifelse(!is.na(qq) & qq == 0, lowest, qq)
    res <- data.table::data.table(
        LR_comb = paste0(tb$ligand, "-", tb$receptor),
        lig_cell_type = tb$sender, lig_expr = tb$mean_L, ligand = tb$ligand,
        rec_cell_type = tb$receiver, rec_expr = tb$mean_R,
        receptor = tb$receptor, scale = tb$scale,
        LR_expr = tb$score, rand_expr = tb$score / ff,
        av_diff = tb$score - tb$score / ff,
        log2fc = lfc, z = zz, pvalue = pp, p.adj = qq,
        LR_cell_comb = paste0(tb$sender, "--", tb$receiver),
        PI = lfc * (-log10(q_eff)),
        lig_nr = tb$n_sender, rec_nr = tb$n_receiver,
        lig_det = tb$det_L, rec_det = tb$det_R,
        fold_coloc = tb$fold, p.adj_coloc = tb$q,
        fold_type = tb$fold_type, p.adj_type = tb$q_type,
        p.adj_shuffle = tb$q_shuffle
    )
    attr(res, "lrmap") <- fit
    attr(res, "lr") <- lr
    attr(res, "ccc_info") <- list(method = "lrmap", test = test,
        values = values, cluster = cluster)
    res[]
}

# The classic consumers of spatCellCellcom() output (combCCcom(), the CCcom
# heatmap and dotplot) key rows by LR_comb and LR_cell_comb; an lrmap result
# has one row per scale, so it must be subset to one scale first.
#' @keywords internal
#' @noRd
.ccc_one_scale <- function(x, what) {
    if ("scale" %in% names(x) && length(unique(x[["scale"]])) > 1L) {
        .gstop(what, "takes one scale of an lrmap result; subset it first,",
            "e.g. res[scale == \"contact\"]")
    }
    invisible(x)
}

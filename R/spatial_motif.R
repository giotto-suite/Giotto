# Multicellular motif enrichment on a spatial network.
#
# Shape follows GiottoClass ADR 0003 (five verb generics) and Giotto ADR 0001
# (wrappers keep their names, dispatch moves to verbs): `cellProximityMotifs()`
# owns everything gobject-shaped and computes nothing, while the numeric work
# is an `analyzeData()` method keyed on a param class.
#
# The engine is a pluggable backend, not a hard-wired dependency. `motifParam`
# is VIRTUAL and defines the contract; `smotifParam` is one implementation of
# it. A second engine is a new concrete subclass plus its own `analyzeData`
# method, contributable from another package with no change here -- the same
# arrangement `pcaParam` has, where GiottoDisk supplies `gramEigenPcaParam`
# from outside this package. Only the `smotifParam` method below knows that
# smotif exists; the router, the result contract and the plots do not.


# classes ####

#' @name motif_param
#' @title Motif enrichment parameters
#' @description
#' Parameter classes for [analyzeData()] motif enrichment. `motifParam` is the
#' VIRTUAL contract; concrete subclasses select an engine. `motifParam()` is
#' the factory.
#'
#' A backend supplies a concrete subclass and an `analyzeData` method on it.
#' The method receives an `igraph` whose vertex names are cell IDs, plus the
#' cell type labels, and must return the result contract described in
#' [cellProximityMotifs()].
#' @param method engine to use. `"auto"` picks the best available:
#'   \pkg{smotif} when installed.
#' @param size motif size: 2, 3 or 4.
#' @param null null model. `"label"` permutes cell type labels over all cells;
#'   `"stratified"` permutes within strata; `"conditional"` holds the observed
#'   pairwise composition fixed, so a 3- or 4-cell motif that is still enriched
#'   is enriched beyond what its constituent pairs already explain.
#' @param n_perm number of null draws.
#' @param set_seed,seed_number seed control.
#' @param ... further engine-specific arguments.
#' @returns a `motifParam` subclass object
NULL

#' @rdname motif_param
#' @exportClass motifParam
setClass("motifParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname motif_param
#' @exportClass smotifParam
setClass("smotifParam", contains = "motifParam")

#' @rdname motif_param
#' @exportClass autoMotifParam
setClass("autoMotifParam", contains = "motifParam")


# param factory ####

#' @rdname motif_param
#' @export
motifParam <- function(method = "auto",
    size = 3L,
    null = c("label", "stratified", "conditional"),
    n_perm = 1000L,
    set_seed = TRUE,
    seed_number = 1234,
    ...) {
    method <- match.arg(tolower(method), c("auto", "smotif"))
    null <- match.arg(null)
    size <- as.integer(size)
    checkmate::assert_int(size, lower = 2L, upper = 4L)
    checkmate::assert_int(as.integer(n_perm), lower = 1L)

    cls <- switch(method,
        "auto" = "autoMotifParam",
        "smotif" = "smotifParam"
    )
    p <- methods::new(cls, param = list(...))
    p$size <- size
    p$null <- null
    p$n_perm <- as.integer(n_perm)
    p$set_seed <- isTRUE(set_seed)
    p$seed_number <- as.numeric(seed_number)
    p
}


# analyzeData(<igraph>, autoMotifParam) ####

# The "auto" sentinel resolves to a concrete engine from what is installed,
# then dispatches. Same shape as pcaParam("auto").
#' @rdname motif_param
setMethod(
    "analyzeData", signature(x = "igraph", param = "autoMotifParam"),
    function(x, param, ...) {
        if (!requireNamespace("smotif", quietly = TRUE)) {
            .gstop(
                "no motif engine is available.",
                "Install smotif:",
                '  package_check("smotif", repository = "github:drieslab/smotif")',
                "or supply another engine as a motifParam subclass.",
                sep = "\n"
            )
        }
        p <- methods::new("smotifParam", param = param@param)
        analyzeData(x, p, ...)
    }
)


# analyzeData(<igraph>, smotifParam) ####

# The only smotif-aware code in this package.
#' @rdname motif_param
setMethod(
    "analyzeData", signature(x = "igraph", param = "smotifParam"),
    function(x, param, cell_type = NULL, strata = NULL, anchored_on = NULL,
    ...) {
        package_check("smotif", repository = "github:drieslab/smotif")
        if (is.null(cell_type)) {
            .gstop("cell_type labels are required")
        }
        n <- igraph::vcount(x)
        if (length(cell_type) != n) {
            .gstop("cell_type must have one entry per graph vertex")
        }
        el <- igraph::as_edgelist(x, names = FALSE)
        sg <- smotif::SpatialGraph(
            nodes = data.table::data.table(
                cell_id = as.character(igraph::V(x)$name %null%
                    seq_len(n)),
                x = 0, y = 0,
                cell_type = as.character(cell_type),
                sample_id = "s",
                region_id = if (is.null(strata)) {
                    NA_character_
                } else {
                    as.character(strata)
                }
            ),
            edges = data.table::data.table(
                source = as.character(igraph::V(x)$name %null%
                    seq_len(n))[el[, 1L]],
                target = as.character(igraph::V(x)$name %null%
                    seq_len(n))[el[, 2L]],
                sample_id = "s"
            )
        )
        if (isTRUE(param$set_seed)) set.seed(param$seed_number)
        smotif::motif_enrichment(
            sg,
            size = param$size,
            n_perm = param$n_perm,
            seed = as.integer(param$seed_number),
            null = param$null,
            strata_column = if (is.null(strata)) NULL else "region_id",
            anchored_on = anchored_on,
            ...
        )
    }
)


# analyzeData(<giotto>, motifParam) ####

# Dispatches on the VIRTUAL class, so it serves every engine.
#' @rdname motif_param
#' @param spat_unit,feat_type spatial unit and feature type
#' @param spatial_network_name name of the spatial network to use
#' @param cluster_column cell metadata column holding cell type labels
#' @param strata_column optional cell metadata column to stratify the null by
#' @param anchored_on optional character vector of cell IDs
#' @param x a `giotto` object or an `igraph`
#' @param cell_type cell type labels, one per vertex (igraph method only)
#' @param strata optional strata, one per vertex (igraph method only)
setMethod(
    "analyzeData", signature(x = "giotto", param = "motifParam"),
    function(x, param,
    spat_unit = NULL, feat_type = NULL,
    spatial_network_name = "Delaunay_network",
    cluster_column = NULL,
    strata_column = NULL,
    anchored_on = NULL,
    ...) {
        spat_unit <- set_default_spat_unit(gobject = x, spat_unit = spat_unit)
        feat_type <- set_default_feat_type(
            gobject = x, spat_unit = spat_unit, feat_type = feat_type
        )
        if (is.null(cluster_column)) {
            .gstop("cluster_column is required")
        }

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

        # Fast path: a disk-backed network with nothing pending can be read
        # straight from parquet in the backend, skipping the igraph
        # materialization entirely. Anything else -- a subset store, another
        # engine, an in-memory network -- goes the ordinary way.
        store <- .motif_edge_store_paths(x, spat_unit, spatial_network_name)
        if (!is.null(store) && is.null(strata_column) &&
            is.null(anchored_on) && is(param, "smotifParam") &&
            param$null %in% c("label", "conditional")) {
            order_dt <- smotif::store_node_order(store$nodes)
            lab <- .motif_align(meta, order_dt$node_id, cluster_column)
            res <- smotif::motif_enrichment_store(
                nodes_path = store$nodes, edges_path = store$edges,
                cell_type = lab, size = param$size, n_perm = param$n_perm,
                seed = as.integer(param$seed_number), null = param$null, ...
            )
            .motif_check_contract(res)
            return(res)
        }

        ig <- .motif_network_as_igraph(x, spat_unit, spatial_network_name)
        vids <- igraph::V(ig)$name
        if (is.null(vids)) {
            .gstop(
                "the spatial network has no vertex names, so its nodes",
                "cannot be matched to cell metadata"
            )
        }
        lab <- .motif_align(meta, vids, cluster_column)
        strata <- NULL
        if (!is.null(strata_column)) {
            if (!strata_column %in% colnames(meta)) {
                .gstop(sprintf(
                    "strata_column '%s' is not a cell metadata column",
                    strata_column
                ))
            }
            strata <- .motif_align(meta, vids, strata_column)
        }

        res <- analyzeData(ig, param,
            cell_type = lab, strata = strata,
            anchored_on = anchored_on, ...
        )
        .motif_check_contract(res)
        res
    }
)


# helpers ####

#' @keywords internal
#' @noRd
.motif_align <- function(meta, vids, column) {
    v <- meta[[column]]
    names(v) <- meta[["cell_ID"]]
    out <- as.character(v[vids])
    if (anyNA(out)) {
        .gstop(sprintf(
            "%d network node(s) have no '%s' value in cell metadata",
            sum(is.na(out)), column
        ))
    }
    out
}

# The network slot is polymorphic as of GiottoClass 0.6.0: an igraph, or a
# GiottoDisk dataStore on a backed project. ADR 0004 is explicit that every new
# consumer of @network owes a branch, so this is it.
#' @keywords internal
#' @noRd
.motif_network_as_igraph <- function(gobject, spat_unit, name) {
    sn <- getSpatialNetwork(gobject,
        spat_unit = spat_unit, name = name,
        output = "spatialNetworkObj", verbose = FALSE
    )
    net <- sn[]
    if (inherits(net, "dataStore")) {
        package_check("GiottoDisk",
            repository = "github:giotto-suite/GiottoDisk"
        )
        return(GiottoDisk::storeRead(net, output = "igraph"))
    }
    if (!inherits(net, "igraph")) {
        .gstop(sprintf(
            "unsupported spatial network storage: %s",
            paste(class(net), collapse = "/")
        ))
    }
    net
}


# Paths to a parquetEdgeStore's parquet files, or NULL when the network is not
# one, when ops are pending on it (the files on disk do not reflect a pending
# subset), or when the backend that reads them is absent.
#' @keywords internal
#' @noRd
.motif_edge_store_paths <- function(gobject, spat_unit, name) {
    if (!requireNamespace("smotif", quietly = TRUE) ||
        !requireNamespace("smotifrs", quietly = TRUE)) {
        return(NULL)
    }
    sn <- try(
        getSpatialNetwork(gobject,
            spat_unit = spat_unit, name = name,
            output = "spatialNetworkObj", verbose = FALSE
        ),
        silent = TRUE
    )
    if (inherits(sn, "try-error")) return(NULL)
    net <- sn[]
    if (!inherits(net, "parquetEdgeStore")) return(NULL)
    if (length(methods::slot(net, "ops")) > 0L) return(NULL)
    root <- methods::slot(net, "path")
    nodes <- list.files(file.path(root, "nodes"), "[.]parquet$",
        full.names = TRUE
    )
    edges <- list.files(file.path(root, "edges"), "[.]parquet$",
        full.names = TRUE
    )
    # one file per subdir today; hive-partitioned writes are a future change
    # on the GiottoDisk side, and this path must not silently read only part
    if (length(nodes) != 1L || length(edges) != 1L) return(NULL)
    list(nodes = nodes, edges = edges)
}

# The contract every engine must satisfy. Checked here rather than trusted, so
# a backend that drifts fails at the boundary instead of inside a plot.
#' @keywords internal
#' @noRd
.motif_check_contract <- function(res) {
    required <- c(
        "motif_id", "topology", "size", "color_tuple", "observed",
        "expected", "sd_null", "z", "fold", "p_enrich", "p_deplete", "p_adj"
    )
    if (!inherits(res, "data.frame")) {
        .gstop("a motif engine must return a data.frame/data.table")
    }
    missing <- setdiff(required, colnames(res))
    if (length(missing)) {
        .gstop(
            "the motif engine's result is missing required column(s):",
            paste(missing, collapse = ", ")
        )
    }
    if (!is.list(res[["color_tuple"]])) {
        .gstop("`color_tuple` must be a list column, one entry per motif")
    }
    invisible(TRUE)
}


# user-facing wrapper ####

#' @title cellProximityMotifs
#' @name cellProximityMotifs
#' @description
#' Enrichment of multicellular motifs -- recurrent arrangements of 2, 3 or 4
#' neighbouring cells -- in a spatial network.
#'
#' Where [cellProximityEnrichment()] asks which *pairs* of cell types sit
#' together more than chance, this asks the same of larger neighbourhoods, and
#' distinguishes their shape: a triangle of three cell types is a different
#' motif from a chain of the same three.
#' @inheritParams data_access_params
#' @param spatial_network_name name of the spatial network to use
#' @param cluster_column cell metadata column holding cell type labels
#' @param size motif size: 2, 3 or 4
#' @param null null model. `"label"` permutes labels over all cells;
#'   `"stratified"` permutes within `strata_column`; `"conditional"` holds the
#'   observed pairwise composition fixed.
#' @param n_perm number of null draws
#' @param strata_column optional cell metadata column to stratify the null by
#' @param anchored_on optional character vector of cell IDs; only motifs
#'   touching one of them are counted
#' @param method motif engine. `"auto"` uses whatever is installed.
#' @param set_seed,seed_number seed control
#' @param ... passed to the engine
#' @returns a `data.table`, one row per motif class, with columns
#'   `motif_id`, `topology`, `size`, `color_tuple` (list column of cell type
#'   labels in canonical orbit order), `observed`, `expected`, `sd_null`, `z`,
#'   `fold`, `p_enrich`, `p_deplete` and `p_adj`.
#' @details
#' Under the default `"label"` null, a 3- or 4-cell motif built out of an
#' attracting pair will look enriched simply because the pair is. Use
#' `null = "conditional"` to ask the higher-order question: it holds the
#' observed pairwise cell-type composition fixed, so what remains is enrichment
#' beyond what the pairs already explain.
#'
#' `topology` names the arrangement: `edge`; `open` (a chain) and `closed` (a
#' triangle) at size 3; and `claw`, `path`, `cycle`, `paw`, `diamond` and `K4`
#' at size 4. Positions in `color_tuple` are structural roles, not sorted, so
#' `A-B-B` and `B-A-A` are different motifs.
#'
#' Requires an engine. \pkg{smotif} is the reference one, and its
#' \pkg{smotifrs} backend is what makes size 4 tractable at scale.
#' @section Building the network for a large section:
#' On a large section the network, not the motif counting, is the slow step --
#' at 170,000 cells the enrichment takes seconds while the default Delaunay
#' construction takes minutes. Two arguments to
#' [GiottoClass::createSpatialNetwork()] matter:
#'
#' ```r
#' g <- createSpatialNetwork(g,
#'     method = "Delaunay",
#'     delaunay_method = "delaunayn_geometry",   # same graph, far faster
#'     maximum_distance_delaunay = <explicit>    # trim gap-spanning edges
#' )
#' ```
#'
#' `"delaunayn_geometry"` returns the *identical* triangulation to the default
#' `"deldir"` -- the same edges, not an approximation -- in 0.94 s against
#' ~363 s at 200,000 points. It needs the \pkg{geometry} package.
#'
#' Set `maximum_distance_delaunay` explicitly rather than leaving it at
#' `"auto"`. A Delaunay triangulation joins cells across empty space, so a
#' section with gaps gets long edges that are not adjacency in any biological
#' sense, and those edges become motifs. `"auto"` does trim them -- it is the
#' boxplot upper whisker -- but it is one global threshold on a section whose
#' density varies, and on a clustered 20,000-cell set it stranded 124 cells with
#' no edges at all. A cell with no edges is not in the network and so is in no
#' motif; `createSpatialNetwork()` now reports how many. Picking the cutoff from
#' a high quantile of the unfiltered edge lengths keeps that trade-off visible.
#'
#' kNN is faster still but is a different neighbourhood definition, not a speed
#' setting: measured against Delaunay on the same points it shares only ~59% of
#' edges, and `k = 6` gives mean degree 7.1 against Delaunay's 5.9, so it is
#' denser than it looks. [GiottoClass::radiusNetworkParam()] is the option
#' closest in spirit to "cells within this distance".
#' @seealso [cellProximityEnrichment()] for the pairwise case,
#'   [motifParam()] for the parameter object.
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' m <- cellProximityMotifs(g, cluster_column = "leiden_clus", size = 3L)
#' head(m)
#' }
#' @export
cellProximityMotifs <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    spatial_network_name = "Delaunay_network",
    cluster_column,
    size = 3L,
    null = c("label", "stratified", "conditional"),
    n_perm = 1000L,
    strata_column = NULL,
    anchored_on = NULL,
    method = "auto",
    set_seed = TRUE,
    seed_number = 1234,
    ...) {
    p <- motifParam(
        method = method, size = size, null = match.arg(null),
        n_perm = n_perm, set_seed = set_seed, seed_number = seed_number
    )
    analyzeData(gobject, p,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column,
        strata_column = strata_column,
        anchored_on = anchored_on,
        ...
    )
}

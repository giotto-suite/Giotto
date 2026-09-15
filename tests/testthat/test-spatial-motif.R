# The motif layer is meant to be engine-agnostic: `motifParam` is the VIRTUAL
# contract and `smotifParam` is one implementation, the same arrangement
# `pcaParam` has where GiottoDisk contributes `gramEigenPcaParam` from outside
# this package. That claim is only worth making if a second engine can actually
# be attached without touching any file here, so the first test below does
# exactly that with a throwaway engine backed by igraph::motifs(). If it ever
# needs to reach into smotif, the abstraction has leaked.

# --- a second engine, defined entirely in this test file -------------------

# testthat runs tests in a child of the package namespace, which is locked, so
# both registrations name an assignable environment explicitly.
setClass("fakeMotifParam", contains = "motifParam", where = globalenv())

setMethod(
    "analyzeData", signature(x = "igraph", param = "fakeMotifParam"),
    where = globalenv(),
    definition = function(x, param, cell_type = NULL, strata = NULL,
                          anchored_on = NULL, ...) {
        # counts connected triples, ignores colour; enough to fill the contract
        n <- sum(igraph::count_triangles(x)) / 3
        data.table::data.table(
            motif_id = "size3_closed_any",
            topology = "closed",
            size = 3L,
            color_tuple = list(c("any", "any", "any")),
            observed = as.numeric(n),
            expected = as.numeric(n),
            sd_null = 1,
            z = 0,
            fold = 1,
            p_enrich = 1,
            p_deplete = 1,
            p_adj = 1
        )
    }
)

test_that("a second engine attaches without any change to the motif layer", {
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    p <- methods::new("fakeMotifParam", param = list(
        size = 3L, null = "label", n_perm = 10L,
        set_seed = TRUE, seed_number = 1
    ))
    # routed through the same analyzeData(giotto, motifParam) method that the
    # real engine uses -- dispatch is on the VIRTUAL class
    res <- analyzeData(g, p, cluster_column = "leiden_clus")
    expect_s3_class(res, "data.table")
    expect_identical(res$topology, "closed")
    expect_gt(res$observed, 0)
})

test_that("motifParam is virtual and its concrete subclasses are not", {
    expect_true(isVirtualClass("motifParam"))
    expect_false(isVirtualClass("smotifParam"))
    expect_true(extends("smotifParam", "motifParam"))
    expect_true(extends("smotifParam", "analyzeParam"))
    # the seam other packages extend through must be exported
    expect_true("motifParam" %in% getNamespaceExports("Giotto"))
})

test_that("the factory validates and records its settings", {
    p <- motifParam(size = 4L, null = "conditional", n_perm = 50L)
    expect_s4_class(p, "autoMotifParam")
    expect_identical(p$size, 4L)
    expect_identical(p$null, "conditional")
    expect_identical(p$n_perm, 50L)
    expect_s4_class(motifParam(method = "smotif"), "smotifParam")

    expect_error(motifParam(size = 5L))
    expect_error(motifParam(size = 1L))
    expect_error(motifParam(n_perm = 0L))
    expect_error(motifParam(null = "nonsense"))
    expect_error(motifParam(method = "nonsense"))
})

test_that("the result contract is checked, not trusted", {
    ok <- data.table::data.table(
        motif_id = "a", topology = "closed", size = 3L,
        color_tuple = list(c("A", "B", "C")), observed = 1, expected = 1,
        sd_null = 1, z = 0, fold = 1, p_enrich = 1, p_deplete = 1, p_adj = 1
    )
    expect_true(.motif_check_contract(ok))

    expect_error(.motif_check_contract(data.table::data.table(motif_id = "a")),
        "missing required column"
    )
    expect_error(.motif_check_contract(list(a = 1)), "data.frame")

    flat <- data.table::copy(ok)
    flat$color_tuple <- "A-B-C"
    expect_error(.motif_check_contract(flat), "list column")
})

test_that("cluster_column problems are reported clearly", {
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    p <- motifParam(size = 3L, n_perm = 5L)
    expect_error(analyzeData(g, p), "cluster_column is required")
    expect_error(
        analyzeData(g, p, cluster_column = "not_a_column"),
        "not a cell metadata column"
    )
    expect_error(
        analyzeData(g, p,
            cluster_column = "leiden_clus",
            strata_column = "not_a_column"
        ),
        "not a cell metadata column"
    )
})

test_that("cellProximityMotifs returns the documented contract", {
    skip_if_not_installed("GiottoData")
    skip_if_not_installed("smotif")
    g <- GiottoData::loadGiottoMini("visium")
    m <- cellProximityMotifs(g,
        cluster_column = "leiden_clus",
        size = 3L, n_perm = 49L
    )
    expect_identical(names(m), c(
        "motif_id", "topology", "size", "color_tuple", "observed", "expected",
        "sd_null", "z", "fold", "p_enrich", "p_deplete", "p_adj"
    ))
    expect_true(all(m$topology %in% c("open", "closed")))
    expect_true(all(m$size == 3L))
    expect_equal(sum(m$observed), attr(m, "n_instances"))
    expect_true(all(m$p_enrich > 0 & m$p_enrich <= 1))
    # positions are structural roles, so a wedge centre is not sorted away
    expect_true(any(vapply(m$color_tuple, function(x) x[1] != x[2], NA)))
})

test_that("size 2 motifs agree with the pairwise proximity counts", {
    skip_if_not_installed("GiottoData")
    skip_if_not_installed("smotif")
    g <- GiottoData::loadGiottoMini("visium")
    m2 <- cellProximityMotifs(g,
        cluster_column = "leiden_clus", size = 2L, n_perm = 19L
    )
    cp <- cellProximityEnrichment(g,
        cluster_column = "leiden_clus", number_of_simulations = 19L
    )
    # the two entry points count the same edges
    expect_equal(sum(m2$observed), sum(cp$enrichm_res$original))
})

# --- disk-backed networks --------------------------------------------------
#
# As of GiottoClass 0.6.0 the @network slot is polymorphic: an igraph, or a
# GiottoDisk dataStore on a backed project. ADR 0004 asks every new consumer to
# owe that branch. There is also a fast path that reads the store's parquet
# directly, and the trap it has to avoid is a *subsetted* store -- its pending
# ops are not reflected in the files on disk, so reading them raw would
# silently analyse the whole network instead of the subset.

.backed_gobject <- function(dir, n = 200L, seed = 5L) {
    set.seed(seed)
    locs <- data.table::data.table(
        cell_ID = sprintf("c%03d", seq_len(n)),
        sdimx = runif(n, 0, 100), sdimy = runif(n, 0, 100)
    )
    m <- matrix(rpois(n * 5L, 3), nrow = 5L,
        dimnames = list(paste0("g", 1:5), locs$cell_ID)
    )
    g <- GiottoClass::createGiottoObject(
        expression = m, spatial_locs = locs, backend = dir
    )
    g <- GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(
            cell_ID = locs$cell_ID, ct = sample(c("A", "B", "C"), n, TRUE)
        ),
        by_column = TRUE, column_cell_ID = "cell_ID"
    )
    GiottoClass::createSpatialNetwork(g,
        method = "Delaunay", name = "Delaunay_network"
    )
}

test_that("a disk-backed network gives the same answer as an in-memory one", {
    skip_if_not_installed("GiottoDisk")
    skip_if_not_installed("smotif")
    withr::local_options(giotto.check_valid = FALSE)
    dir <- file.path(withr::local_tempdir(), "proj")
    g <- .backed_gobject(dir)

    net <- GiottoClass::getSpatialNetwork(g,
        name = "Delaunay_network", output = "spatialNetworkObj"
    )[]
    expect_true(inherits(net, "dataStore"))

    disk <- cellProximityMotifs(g,
        cluster_column = "ct", size = 3L, n_perm = 49L
    )

    # the same data, unbacked
    meta <- GiottoClass::pDataDT(g)
    locs <- GiottoClass::getSpatialLocations(g, output = "data.table")
    m <- matrix(1, nrow = 2L, ncol = nrow(locs),
        dimnames = list(c("g1", "g2"), locs$cell_ID)
    )
    g2 <- GiottoClass::createGiottoObject(expression = m, spatial_locs = locs)
    g2 <- GiottoClass::addCellMetadata(g2,
        new_metadata = data.frame(cell_ID = meta$cell_ID, ct = meta$ct),
        by_column = TRUE, column_cell_ID = "cell_ID"
    )
    g2 <- GiottoClass::createSpatialNetwork(g2,
        method = "Delaunay", name = "Delaunay_network"
    )
    mem <- cellProximityMotifs(g2,
        cluster_column = "ct", size = 3L, n_perm = 49L
    )

    expect_identical(attr(disk, "n_instances"), attr(mem, "n_instances"))
    expect_setequal(disk$motif_id, mem$motif_id)
    j <- merge(
        disk[, list(motif_id, a = observed)],
        mem[, list(motif_id, b = observed)],
        by = "motif_id"
    )
    expect_equal(j$a, j$b)
})

test_that("the direct-parquet fast path declines a subsetted store", {
    skip_if_not_installed("GiottoDisk")
    skip_if_not_installed("smotifrs")
    withr::local_options(giotto.check_valid = FALSE)
    dir <- file.path(withr::local_tempdir(), "proj")
    g <- .backed_gobject(dir)

    # a clean store can be read straight from disk
    expect_false(is.null(
        .motif_edge_store_paths(g, "cell", "Delaunay_network")
    ))

    # a subsetted one carries pending ops the parquet files do not reflect,
    # so the fast path must refuse it rather than analyse the whole network
    sn <- GiottoClass::getSpatialNetwork(g,
        name = "Delaunay_network", output = "spatialNetworkObj"
    )
    ids <- GiottoClass::pDataDT(g)$cell_ID[1:40]
    sn[] <- sn[][ids]
    g2 <- GiottoClass::setSpatialNetwork(g, sn, verbose = FALSE)
    expect_gt(length(methods::slot(GiottoClass::getSpatialNetwork(g2,
        name = "Delaunay_network", output = "spatialNetworkObj"
    )[], "ops")), 0L)
    expect_null(.motif_edge_store_paths(g2, "cell", "Delaunay_network"))

    # and the subset is still analysable, through the ordinary path
    skip_if_not_installed("smotif")
    r <- cellProximityMotifs(g2, cluster_column = "ct", size = 3L, n_perm = 19L)
    expect_gt(nrow(r), 0L)
    expect_lt(attr(r, "n_instances"), 1e5)
})

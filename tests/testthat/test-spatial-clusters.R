# spatialSplitCluster() and identifyTMAcores() had no coverage. Both were
# unusable on any in-memory spatial network for the whole 0.6.0 cycle --
# `spatIDs()` and the old GiottoClass helper both read `@network` as the
# from/to table it stopped holding -- and nothing caught it. These pin the
# shape both functions need from `.spat_net_to_igraph()`.

rlang::local_options(lifecycle_verbosity = "quiet")

.net_gobject <- function(method = "Delaunay") {
    g <- test_data$vis
    GiottoClass::createSpatialNetwork(
        g,
        method = method,
        name = paste0(method, "_network"),
        verbose = FALSE
    )
}


# .spat_net_to_igraph ####

test_that(".spat_net_to_igraph returns an undirected graph with named nodes", {
    g <- .net_gobject("Delaunay")
    sn <- GiottoClass::getSpatialNetwork(
        g, name = "Delaunay_network", output = "spatialNetworkObj"
    )

    net <- .spat_net_to_igraph(sn)
    expect_s3_class(net, "igraph")
    expect_false(igraph::is_directed(net))
    # vertex names are the cell_IDs both callers match metadata on
    expect_true(all(nzchar(names(igraph::V(net)))))
    expect_gt(igraph::vcount(net), 0L)
    expect_gt(igraph::ecount(net), 0L)
})

test_that(".spat_net_to_igraph drops edge attributes unless asked", {
    g <- .net_gobject("Delaunay")
    sn <- GiottoClass::getSpatialNetwork(
        g, name = "Delaunay_network", output = "spatialNetworkObj"
    )
    stored <- igraph::as.igraph(sn)
    # the stored graph carries them; the shaped one should not
    expect_true(length(igraph::edge_attr_names(stored)) > 0L)

    expect_length(igraph::edge_attr_names(.spat_net_to_igraph(sn)), 0L)
    expect_identical(
        igraph::edge_attr_names(.spat_net_to_igraph(sn, attr = "distance")),
        "distance"
    )
})

test_that(".spat_net_to_igraph undirects kNN without collapsing edges", {
    g <- .net_gobject("kNN")
    sn <- GiottoClass::getSpatialNetwork(
        g, name = "kNN_network", output = "spatialNetworkObj"
    )
    stored <- igraph::as.igraph(sn)
    skip_if_not(igraph::is_directed(stored), "kNN network not stored directed")

    net <- .spat_net_to_igraph(sn)
    expect_false(igraph::is_directed(net))
    # mode = "each": reciprocal pairs stay two edges rather than collapsing
    expect_equal(igraph::ecount(net), igraph::ecount(stored))
    expect_equal(igraph::vcount(net), igraph::vcount(stored))
})


# the two callers ####

test_that("spatialSplitCluster splits a cluster column by spatial adjacency", {
    g <- .net_gobject("Delaunay")
    cm <- GiottoClass::pDataDT(g)
    skip_if_not("leiden_clus" %in% colnames(cm), "no leiden_clus in mini")

    out <- spatialSplitCluster(
        g,
        spatial_network_name = "Delaunay_network",
        cluster_col = "leiden_clus",
        return_gobject = TRUE
    )
    res <- GiottoClass::pDataDT(out)

    expect_true("leiden_clus_split" %in% colnames(res))
    # the regression: an empty vertex set produced no assignments at all
    expect_false(all(is.na(res$leiden_clus_split)))
    expect_gte(
        data.table::uniqueN(res$leiden_clus_split),
        data.table::uniqueN(res$leiden_clus)
    )
})

test_that("identifyTMAcores assigns core ids", {
    g <- .net_gobject("Delaunay")

    out <- identifyTMAcores(
        g,
        spatial_network_name = "Delaunay_network",
        return_gobject = TRUE
    )
    res <- GiottoClass::pDataDT(out)

    expect_true("core_id" %in% colnames(res))
    expect_false(all(is.na(res$core_id)))
    expect_gt(data.table::uniqueN(res$core_id), 0L)
})

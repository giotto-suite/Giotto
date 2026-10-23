# Analysis verbs on a giottoMulti: results are written to the multi's joint
# cell metadata, keyed by `sample::cell_ID`.

.mini_multi <- function() {
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    mg <- createGiottoMulti(list(a = g, b = g))
    suppressMessages(normalizeGiotto(mg, verbose = FALSE))
}

test_that("addStatistics computes area per sample on a multi", {
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    mg <- .mini_multi()
    keep <- head(spatIDs(mg), 700)  # all of a, part of b
    mg <- subset(mg, cells = keep)

    mg <- suppressMessages(addStatistics(mg))
    cm <- getCellMetadata(mg, output = "data.table")
    expect_setequal(cm$cell_ID, keep)
    expect_false(anyNA(cm$area))

    ref <- expanse(getPolygonInfo(g, return_giottoPolygon = TRUE),
        output = "data.table")
    a <- cm[list_ID == "a"]
    expect_equal(a$area,
        ref$area[match(sub("^a::", "", a$cell_ID), ref$cell_ID)])
})

test_that("Leiden clustering on a multi writes joint cell metadata", {
    mg <- .mini_multi()
    mg <- suppressMessages(runPCA(mg, feats_to_use = NULL))
    mg <- suppressMessages(
        createNearestNetwork(mg, dimensions_to_use = 1:10, k = 10))

    mg <- doLeidenCluster(mg, resolution = 0.5, n_iterations = 10,
        name = "leiden_multi")
    cm <- getCellMetadata(mg, output = "data.table")
    expect_true("leiden_multi" %in% names(cm))
    expect_false(anyNA(cm$leiden_multi))
    expect_setequal(unique(cm$list_ID), c("a", "b"))
    expect_false("leiden_multi" %in%
        names(getCellMetadata(mg[["a"]], output = "data.table")))
})

test_that("clusterData with a BlusterParam runs on a multi", {
    mg <- .mini_multi()
    mg <- suppressMessages(runPCA(mg, feats_to_use = NULL))

    mg <- clusterData(mg, param = clusterParam("kmeans", centers = 3),
        what = "dimension_reduction", name = "km")
    km <- getCellMetadata(mg, output = "data.table")$km
    expect_length(km, length(spatIDs(mg)))
    expect_length(unique(km), 3L)
})

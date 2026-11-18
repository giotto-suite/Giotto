# runBlockPCA(): block-scaled concatenation of per-feature-type embeddings.
# The two blocks are on very different scales (sd 10 vs 0.1), so an
# unscaled concatenation would be all RNA.

.block_fixture <- function(n = 300, seed = 3, backend = NULL) {
    set.seed(seed)
    cells <- sprintf("c_%04d", seq_len(n))
    xr <- matrix(stats::rnorm(n * 12, sd = 10), n, 12,
        dimnames = list(cells, paste0("Dim.", 1:12)))
    xp <- matrix(stats::rnorm(n * 4, sd = 0.1), n, 4,
        dimnames = list(cells, paste0("Dim.", 1:4)))
    mr <- Matrix::Matrix(matrix(stats::rpois(5 * n, 5), 5, n,
        dimnames = list(paste0("g", 1:5), cells)), sparse = TRUE)
    mp <- Matrix::Matrix(matrix(stats::rpois(3 * n, 5), 3, n,
        dimnames = list(paste0("p", 1:3), cells)), sparse = TRUE)
    g <- createGiottoObject(
        expression = list(cell = list(rna = list(raw = mr),
            protein = list(raw = mp))),
        backend = backend, verbose = FALSE
    )
    g <- setGiotto(g, createDimObj(coordinates = xr, name = "pca",
        spat_unit = "cell", feat_type = "rna", method = "pca"),
        verbose = FALSE)
    g <- setGiotto(g, createDimObj(coordinates = xp, name = "protein.pca",
        spat_unit = "cell", feat_type = "protein", method = "pca"),
        verbose = FALSE)
    list(g = g, xr = xr, xp = xp)
}

.block_ref <- function(blocks, w) {
    do.call(cbind, Map(function(x, wi) {
        x * sqrt(wi / sum(apply(x, 2L, stats::var)))
    }, blocks, w))
}

.block_get <- function(g, name = "blockpca", feat_type = "rna",
    output = "matrix") {
    getDimReduction(g, feat_type = feat_type, reduction_method = "blockpca",
        name = name, output = output)
}

test_that("blocks are scaled to their weights and stored per feature type", {
    fx <- .block_fixture()
    g <- runBlockPCA(fx$g, rotate = FALSE)
    ref <- .block_ref(list(fx$xr, fx$xp), c(0.5, 0.5))
    for (ft in c("rna", "protein")) {
        x <- .block_get(g, feat_type = ft)
        expect_equal(unname(x), unname(ref))
        expect_identical(rownames(x), rownames(fx$xr))
    }
    v <- apply(x, 2L, stats::var)
    expect_equal(c(sum(v[1:12]), sum(v[13:16])), c(0.5, 0.5))

    bi <- .block_get(g, output = "dimObj")@misc$block_input
    expect_identical(bi$feat_types, c("rna", "protein"))
    expect_identical(bi$reduction_names, c("pca", "protein.pca"))
    expect_identical(bi$dimensions_to_use, list(rna = 1:12, protein = 1:4))
    expect_equal(unname(bi$block_share), cbind(rep(1:0, c(12, 4)),
        rep(0:1, c(12, 4))))

    g2 <- runBlockPCA(fx$g, block_weights = c(protein = 1, rna = 3),
        rotate = FALSE, name = "w")
    v2 <- apply(.block_get(g2, name = "w"), 2L, stats::var)
    expect_equal(c(sum(v2[1:12]), sum(v2[13:16])), c(0.75, 0.25))
    expect_equal(.block_get(g2, name = "w", output = "dimObj")@misc$
        block_input$block_weights, c(rna = 0.75, protein = 0.25))
})

test_that("the rotation keeps distances and orders dimensions by variance", {
    fx <- .block_fixture()
    a <- .block_get(runBlockPCA(fx$g, rotate = FALSE))
    g <- runBlockPCA(fx$g)
    b <- .block_get(g)
    expect_equal(as.matrix(stats::dist(b)), as.matrix(stats::dist(a)))
    vb <- apply(b, 2L, stats::var)
    expect_true(all(diff(vb) <= 1e-12))
    obj <- .block_get(g, output = "dimObj")
    expect_equal(obj@misc$eigenvalues, unname(vb))
    expect_equal(sum(obj@misc$eigenvalues), 1)
    sh <- obj@misc$block_input$block_share
    expect_identical(colnames(sh), c("rna", "protein"))
    expect_equal(unname(rowSums(sh)), rep(1, 16))
    expect_equal(unname(colSums(sh * vb)), c(0.5, 0.5))
    expect_identical(b, .block_get(runBlockPCA(fx$g)))
})

test_that("dimensions, weights and feature types are resolved and checked", {
    fx <- .block_fixture()
    g <- runBlockPCA(fx$g, rotate = FALSE,
        dimensions_to_use = list(protein = 1:2, rna = c(1:3, 50)))
    expect_equal(unname(.block_get(g)),
        unname(.block_ref(list(fx$xr[, 1:3], fx$xp[, 1:2]), c(0.5, 0.5))))
    g <- runBlockPCA(fx$g, rotate = FALSE, dimensions_to_use = 1:2)
    expect_identical(ncol(.block_get(g)), 4L)

    expect_error(runBlockPCA(fx$g, dimensions_to_use = list(1:2)),
        "one entry per feature type")
    expect_error(runBlockPCA(fx$g,
        dimensions_to_use = list(rna = 1:2, missing = 1:2)), "feature types")
    expect_error(runBlockPCA(fx$g, dimensions_to_use = 30:40),
        "no stored dimensions")
    expect_error(runBlockPCA(fx$g, block_weights = c(1, 0)), "positive")
    expect_error(runBlockPCA(fx$g, block_weights = c(rna = 1, missing = 1)),
        "feature types")
    expect_error(runBlockPCA(fx$g, block_weights = 1), "Must have length")
    expect_error(runBlockPCA(fx$g, reduction_names = c("a", "b", "c")),
        "\\[runBlockPCA\\] `reduction_names` needs one value")
    expect_error(runBlockPCA(fx$g, feat_types = "rna"))
    expect_error(runBlockPCA(fx$g, feat_types = c("rna", "missing")),
        "\\[runBlockPCA\\] no dimension reduction for feature type\\(s\\)")
})

.block_edges <- function(g, nn_type, name) {
    e <- getNearestNetwork(g, feat_type = "rna", nn_type = nn_type,
        name = name, output = "data.table")
    data.table::setorder(e[, .(from, to, weight)], from, to)[]
}

test_that("the integrated kNN and sNN match createNearestNetwork() on all joint dimensions", {
    fx <- .block_fixture()
    g <- runBlockPCA(fx$g, name = "joint", integrated_k = 15,
        integrated_snn = "giotto")
    ref <- createNearestNetwork(g, feat_type = "rna", type = "sNN", k = 15,
        dim_reduction_to_use = "blockpca", dim_reduction_name = "joint",
        dimensions_to_use = 1:16, keep_knn = TRUE, name = "ref",
        verbose = FALSE)
    expect_true(all(c("joint_kNN", "joint_sNN") %in%
        names(c(g@nn_network$cell$rna$kNN, g@nn_network$cell$rna$sNN))))
    expect_equal(.block_edges(g, "kNN", "joint_kNN"),
        .block_edges(ref, "kNN", "kNN.blockpca"))
    expect_equal(.block_edges(g, "sNN", "joint_sNN"),
        .block_edges(ref, "sNN", "ref"))

    a <- doLeidenCluster(g, feat_type = "rna", nn_network_to_use = "sNN",
        network_name = "joint_sNN", resolution = 0.5, name = "leiden")
    b <- doLeidenCluster(ref, feat_type = "rna", nn_network_to_use = "sNN",
        network_name = "ref", resolution = 0.5, name = "leiden")
    expect_identical(pDataDT(a)$leiden, pDataDT(b)$leiden)
    expect_identical(getNearestNetwork(g, feat_type = "rna", nn_type = "kNN",
        name = "joint_kNN")@misc$nn_input$dimensions_to_use, 1:16)
})

test_that("the integrated kNN feeds runIntegratedUMAP(graph = \"knn\")", {
    fx <- .block_fixture()
    g <- runBlockPCA(fx$g, integrated_k = 15)
    g <- runIntegratedUMAP(g, graph = "knn", network_name = "blockpca_kNN",
        k = 15, name = "joint.umap")
    for (ft in c("rna", "protein")) {
        expect_identical(dim(getDimReduction(g, feat_type = ft,
            reduction_method = "umap", name = "joint.umap",
            output = "matrix")), c(300L, 2L))
    }
    expect_error(runIntegratedUMAP(g, graph = "knn",
        network_name = "blockpca_kNN", k = 16), "larger than the integrated kNN")
})

test_that("integrated network names, sNN choices and checks", {
    fx <- .block_fixture()
    g <- runBlockPCA(fx$g)
    expect_identical(names(g@nn_network$cell$rna$kNN), "blockpca_kNN")
    expect_null(g@nn_network$cell$rna$sNN)
    expect_null(g@nn_network$cell$protein)

    g <- runBlockPCA(fx$g, integrated_k = 10, integrated_snn = "jaccard",
        integrated_network_name = "block", integrated_snn_name = "block_snn")
    expect_identical(names(g@nn_network$cell$rna$kNN), "block")
    expect_identical(names(g@nn_network$cell$rna$sNN), "block_snn")
    e <- .block_edges(g, "sNN", "block_snn")
    expect_true(all(e$weight >= 1 / 15 & e$weight <= 1))

    expect_error(runBlockPCA(fx$g, integrated_k = 300),
        "integrated_k \\(300\\) must be smaller than the number of cells")
    expect_error(runBlockPCA(fx$g, integrated_snn = "louvain"))
})

test_that("runBlockPCA works on a GiottoDisk-backed object", {
    skip_if_not_installed("GiottoDisk")
    vault <- withr::local_tempdir()
    a <- .block_get(runBlockPCA(.block_fixture()$g))
    g <- runBlockPCA(.block_fixture(backend = file.path(vault, "v"))$g)
    expect_false(is.null(g@source))
    expect_equal(.block_get(g), a)
})

test_that("cells are aligned by ID and bad embeddings are refused", {
    fx <- .block_fixture()
    ref <- .block_get(runBlockPCA(fx$g, rotate = FALSE))
    o <- rev(seq_len(nrow(fx$xp)))
    g <- setGiotto(fx$g, createDimObj(coordinates = fx$xp[o, ],
        name = "protein.pca", spat_unit = "cell", feat_type = "protein",
        method = "pca"), verbose = FALSE)
    expect_equal(.block_get(runBlockPCA(g, rotate = FALSE)), ref)

    g <- setGiotto(fx$g, createDimObj(coordinates = fx$xp[-(1:4), ],
        name = "protein.pca", spat_unit = "cell", feat_type = "protein",
        method = "pca"), verbose = FALSE)
    expect_error(runBlockPCA(g), "4 missing from protein")

    xp <- fx$xp
    xp[] <- 1
    g <- setGiotto(fx$g, createDimObj(coordinates = xp,
        name = "protein.pca", spat_unit = "cell", feat_type = "protein",
        method = "pca"), verbose = FALSE)
    expect_error(runBlockPCA(g), "protein has no variance")

    xp <- fx$xp
    xp[3, 2] <- NA
    g <- setGiotto(fx$g, createDimObj(coordinates = xp,
        name = "protein.pca", spat_unit = "cell", feat_type = "protein",
        method = "pca"), verbose = FALSE)
    expect_error(runBlockPCA(g), "protein has NA or infinite values for 1")
})

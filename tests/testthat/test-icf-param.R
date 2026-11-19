# The edgeicg route of findInteractionChangedFeats(): the wrapper owns the
# gobject plumbing (defaults, labels, locations, library size, polygons) and the
# numeric core is analyzeData(<expression>, edgeIcfParam), which calls the
# edgeicg engine. The tests check the class wiring, that the wrapper reproduces
# a direct edgeicg::edgeicg() call on the same inputs, and that the classic
# route is unchanged.

test_that("icfParam is virtual and edgeIcfParam is its concrete subclass", {
    expect_true(isVirtualClass("icfParam"))
    expect_false(isVirtualClass("edgeIcfParam"))
    expect_true(extends("edgeIcfParam", "analyzeParam"))
    expect_true("icfParam" %in% getNamespaceExports("Giotto"))
    p <- icfParam("edgeicg", workers = 2, block_genes = 50)
    expect_s4_class(p, "edgeIcfParam")
    expect_identical(p$workers, 2L)
    expect_identical(p$expression_values, "raw")
    expect_error(icfParam("edgeicg", block_genes = 0), "block_genes")
})

.icf_sim_giotto <- function() {
    sim <- edgeicg::simulate_edgeicg(n = 2500, G = 100, seed = 7, gamma = 1.2)
    ids <- colnames(sim$counts)
    g <- createGiottoObject(
        expression = sim$counts,
        spatial_locs = data.frame(
            sdimx = sim$coords[, 1], sdimy = sim$coords[, 2], cell_ID = ids
        ),
        verbose = FALSE
    )
    g <- addCellMetadata(g,
        new_metadata = data.frame(cell_ID = ids, ct = sim$types),
        by_column = TRUE, column_cell_ID = "cell_ID"
    )
    list(g = g, sim = sim)
}

test_that("method = 'edgeicg' reproduces a direct edgeicg() call", {
    skip_if_not_installed("edgeicg")
    s <- .icf_sim_giotto()
    res <- findInteractionChangedFeats(s$g,
        cluster_column = "ct", method = "edgeicg",
        min_shift = 150, B = 8, report_q = 1
    )
    expect_s3_class(res, "icfObject")
    sc <- res$ICFscores
    expect_true(all(c("feats", "cell_type", "int_cell_type", "scale", "z",
        "p.adj", "class") %in% colnames(sc)))
    direct <- edgeicg::edgeicg(s$sim$counts, s$sim$coords, s$sim$types,
        lib = Matrix::colSums(s$sim$counts), min_shift = 150, B = 8,
        verbose = FALSE
    )
    d <- edgeicg::edgeicg_results(direct, "contact", q_thr = 1)
    m <- merge(sc[scale == "contact"], d,
        by.x = c("cell_type", "int_cell_type", "feats"),
        by.y = c("receiver", "sender", "gene")
    )
    expect_equal(nrow(m), nrow(d))
    expect_equal(m$z.x, m$z.y, tolerance = 1e-8)
    # a true effect of the simulator is found
    hit <- sc[scale == "contact" & cell_type == "t1" & int_cell_type == "t2" &
        feats %in% s$sim$truth[receiver == "t1", gene] & p.adj < 0.05]
    expect_gt(nrow(hit), 5)
})

test_that("the classic route ignores the new arguments", {
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    a <- findInteractionChangedFeats(g,
        cluster_column = "leiden_clus",
        selected_feats = c("Gna12", "Ccnd2", "Btbd17"),
        nr_permutations = 10, do_parallel = FALSE
    )
    b <- findInteractionChangedFeats(g,
        cluster_column = "leiden_clus",
        selected_feats = c("Gna12", "Ccnd2", "Btbd17"),
        nr_permutations = 10, do_parallel = FALSE, method = "classic"
    )
    expect_identical(a$ICFscores, b$ICFscores)
    expect_identical(a$test_info$test, "permutation")
})

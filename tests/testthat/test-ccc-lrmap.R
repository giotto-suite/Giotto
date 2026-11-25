# The lrmap route of spatCellCellcom(): the wrapper owns the gobject plumbing
# (defaults, labels, locations, library size) and the numeric core is
# analyzeData(<expression>, lrmapCccParam), which calls the lrmap engine. The
# tests check the class wiring, that the wrapper reproduces a direct
# lrmap::lrmap() call on the same inputs (with and without leak removal), the
# plots, and that the classic route still runs.

test_that("cccParam is virtual and lrmapCccParam is its concrete subclass", {
    expect_true(isVirtualClass("cccParam"))
    expect_false(isVirtualClass("lrmapCccParam"))
    expect_true(extends("lrmapCccParam", "analyzeParam"))
    p <- cccParam("lrmap", threads = 2, sigma_p = 40)
    expect_s4_class(p, "lrmapCccParam")
    expect_identical(p$threads, 2L)
    expect_identical(p$expression_values, "raw")
    expect_true(isTRUE(p$leak))
    expect_error(cccParam("lrmap", scales = "far"), "unknown scales")
    expect_error(cccParam("lrmap", test = "nope"), "test")
})

.ccc_sim_giotto <- function(...) {
    sim <- lrmap::lrmap_simulate("touch", size = 700, n_clusters = 6,
        induce = 2, seed = 4, ...)
    ids <- paste0("c", seq_len(ncol(sim$counts)))
    colnames(sim$counts) <- ids
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

.ccc_run <- function(g, lr, ...) {
    suppressWarnings(spatCellCellcom(g,
        cluster_column = "ct", feat_set_1 = lr$ligand,
        feat_set_2 = lr$receptor, B = 15, min_shift = 150, ...
    ))
}

test_that("method = 'lrmap' reproduces a direct lrmap() call (no leak)", {
    skip_if_not_installed("lrmap")
    s <- .ccc_sim_giotto()
    lr <- s$sim$lr[1:6, c("ligand", "receptor")]
    res <- .ccc_run(s$g, lr, leak = FALSE)
    expect_s3_class(res, "data.table")
    expect_true(all(c("LR_comb", "lig_cell_type", "rec_cell_type", "scale",
        "LR_expr", "log2fc", "p.adj", "PI", "fold_type", "p.adj_coloc") %in%
        colnames(res)))
    direct <- suppressWarnings(lrmap::lrmap(s$sim$counts, s$sim$coords,
        s$sim$types, lr, lib = Matrix::colSums(s$sim$counts), B = 15,
        min_shift = 150, verbose = FALSE))
    d <- data.table::as.data.table(direct$table)
    m <- merge(res, d,
        by.x = c("lig_cell_type", "rec_cell_type", "ligand", "receptor", "scale"),
        by.y = c("sender", "receiver", "ligand", "receptor", "scale")
    )
    expect_equal(nrow(m), nrow(d))
    expect_equal(m$LR_expr, m$score, tolerance = 1e-10)
    expect_equal(m$p.adj, m$q_type, tolerance = 1e-10)
    expect_equal(m$p.adj_coloc, m$q, tolerance = 1e-10)
    # the planted interface signal (L in A -> R in B) is called at contact
    hit <- res[scale == "contact" & lig_cell_type == "A" &
        rec_cell_type == "B" & ligand == "L"]
    expect_gt(hit$z, 2)
    expect_false(is.null(attr(res, "lrmap")))
})

test_that("leak = TRUE estimates the leak on all features, as lrmap_leak()", {
    skip_if_not_installed("lrmap")
    s <- .ccc_sim_giotto(leak = 0.05)
    lr <- s$sim$lr[1:4, c("ligand", "receptor")]
    res <- .ccc_run(s$g, lr, scales = "contact")
    lib <- Matrix::colSums(s$sim$counts)
    lk <- lrmap::lrmap_leak(s$sim$counts, s$sim$coords, s$sim$types,
        lib = lib, verbose = FALSE)
    direct <- suppressWarnings(lrmap::lrmap(s$sim$counts, s$sim$coords,
        s$sim$types, lr, lib = lib, leak = lk, B = 15, min_shift = 150,
        scales = "contact", verbose = FALSE))
    expect_equal(sort(res$LR_expr), sort(direct$table$score), tolerance = 1e-10)
})

test_that("lr_table with complexes, ignored classic arguments, one-scale guard", {
    skip_if_not_installed("lrmap")
    s <- .ccc_sim_giotto()
    lrt <- data.frame(ligand = c("L", "L"), receptor = c("R", "RG1"),
        receptor_subunits = c("R", "R_G1"))
    expect_message(
        res <- suppressWarnings(spatCellCellcom(s$g, cluster_column = "ct",
            lr_table = lrt, random_iter = 5, B = 10, min_shift = 150,
            leak = FALSE, scales = c("contact", "paracrine"))),
        "Ignoring `random_iter`"
    )
    expect_true("L-RG1" %in% res$LR_comb)
    expect_error(combCCcom(res, res), "one scale")
    expect_error(spatCellCellcom(s$g, cluster_column = "ct"),
        "feat_set_1")
})

test_that("lrmap plots return ggplots", {
    skip_if_not_installed("lrmap")
    s <- .ccc_sim_giotto()
    lr <- s$sim$lr[1:6, c("ligand", "receptor")]
    res <- .ccc_run(s$g, lr, leak = FALSE,
        scales = c("contact", "paracrine", "niche", "expression"))
    args <- list(gobject = s$g, cccResult = res, show_plot = FALSE,
        return_plot = TRUE, save_plot = FALSE)
    expect_s3_class(do.call(plotLRmapOverview, args), "ggplot")
    for (v in c("tissue", "zoom")) {
        expect_s3_class(do.call(plotLRmapPair, c(args, list(sender = "A",
            receiver = "B", ligand = "L", receptor = "R", view = v,
            zoom_size = 200))), "ggplot")
    }
    expect_s3_class(do.call(plotLRmapCompare, args), "ggplot")
    # a subset keeps the fit, so the plots still work on it
    expect_s3_class(plotLRmapOverview(s$g, res[scale == "contact"],
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE), "ggplot")
})

test_that("method = 'classic' still runs the original test", {
    skip_if_not_installed("lrmap")
    s <- .ccc_sim_giotto()
    g <- normalizeGiotto(s$g, verbose = FALSE)
    g <- createSpatialNetwork(g, method = "Delaunay", verbose = FALSE)
    res <- spatCellCellcom(g, cluster_column = "ct", feat_set_1 = "L",
        feat_set_2 = "R", method = "classic", random_iter = 5,
        do_parallel = FALSE, verbose = "none")
    expect_true(all(c("LR_comb", "LR_cell_comb", "PI") %in% colnames(res)))
    expect_false("scale" %in% colnames(res))
})

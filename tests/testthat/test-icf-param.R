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

test_that("edgeicg is the default method; classic stays available", {
    expect_identical(eval(formals(findInteractionChangedFeats)$method)[1], "edgeicg")
    expect_identical(findICF, findInteractionChangedFeats)
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    b <- findInteractionChangedFeats(g,
        cluster_column = "leiden_clus",
        selected_feats = c("Gna12", "Ccnd2", "Btbd17"),
        nr_permutations = 10, do_parallel = FALSE, method = "classic"
    )
    expect_identical(b$test_info$test, "permutation")
})

test_that("classic-test arguments are reported as ignored by the default route", {
    skip_if_not_installed("edgeicg")
    s <- .icf_sim_giotto()
    expect_message(
        findInteractionChangedFeats(s$g,
            cluster_column = "ct", nr_permutations = 10, min_shift = 150, B = 4,
            selected_feats = rownames(s$sim$counts)[1:10]
        ),
        "Ignoring `nr_permutations`"
    )
    expect_no_message(
        findInteractionChangedFeats(s$g,
            cluster_column = "ct", min_shift = 150, B = 4,
            selected_feats = rownames(s$sim$counts)[1:10]
        ),
        message = "Ignoring"
    )
})

test_that("selected_feats: the leak is fitted on all features unless leak_genes = 'tested'", {
    skip_if_not_installed("edgeicg", minimum_version = "0.3.0")
    s <- .icf_sim_giotto()
    sel <- rownames(s$sim$counts)[1:10]
    run <- function(...) suppressMessages(findInteractionChangedFeats(s$g,
        cluster_column = "ct", min_shift = 150, B = 4, selected_feats = sel,
        verbose = FALSE, ...))
    full <- suppressMessages(findInteractionChangedFeats(s$g,
        cluster_column = "ct", min_shift = 150, B = 4, verbose = FALSE))
    a <- run(); tt <- run(leak_genes = "tested")
    expect_equal(sort(a$edgeicg$genes), sort(intersect(sel, full$edgeicg$genes)))
    expect_match(a$edgeicg$settings$leak_genes, "^all")
    expect_match(tt$edgeicg$settings$leak_genes, "^tested")
    expect_equal(a$edgeicg$pairs$lambda, full$edgeicg$pairs$lambda)
})

test_that("a GiottoDisk store is read by edgeicg without loading and equals the matrix route", {
    skip_if_not_installed("edgeicg")
    skip_if_not_installed("GiottoDisk")
    sim <- edgeicg::simulate_edgeicg(n = 2500, G = 100, seed = 7, gamma = 1.2)
    path <- tempfile(fileext = ".parquet")
    on.exit(unlink(path, recursive = TRUE))
    st <- GiottoDisk::storeWrite(GiottoDisk::parquetExprStore(path = path), sim$counts)
    ids <- colnames(sim$counts)
    coords <- sim$coords
    rownames(coords) <- ids
    p <- icfParam("edgeicg", min_shift = 150, B = 8, window_cells = 700)
    args <- list(param = p, coords = coords,
        cell_type = stats::setNames(sim$types, ids),
        lib = stats::setNames(Matrix::colSums(sim$counts), ids))
    f_mem <- do.call(analyzeData, c(list(x = sim$counts), args))
    f_st <- do.call(analyzeData, c(list(x = st), args))
    expect_gt(f_st$n_windows, 3L)
    expect_equal(f_st$scale_fits$contact$z, f_mem$scale_fits$contact$z,
        tolerance = 1e-8)
    expect_equal(f_st$screen$z, f_mem$screen$z, tolerance = 1e-8)
    # feature subsetting keeps the store lazy
    f_sub <- do.call(analyzeData, c(list(x = st,
        feats = rownames(sim$counts)[1:60]), args))
    expect_identical(f_sub$genes, rownames(sim$counts)[1:60])
})

test_that("touching edges are the same from in-memory and GiottoDisk polygons", {
    skip_if_not_installed("edgeicg")
    skip_if_not_installed("GiottoDisk")
    # a row of 4 touching squares plus one apart
    sq <- function(i, x, r = 5) cbind(id = i, part = 1,
        x = x + c(-r, r, r, -r), y = c(-r, -r, r, r), hole = 0)
    m <- do.call(rbind, lapply(1:5, function(i) sq(i, c(10, 20, 30, 40, 90)[i])))
    sv <- terra::vect(m, type = "polygons")
    sv$poly_ID <- paste0("c", 1:5)
    gp <- createGiottoPolygon(sv, name = "cell", verbose = FALSE)
    g <- createGiottoObject(
        expression = Matrix::Matrix(1, 2, 5, sparse = TRUE,
            dimnames = list(c("a", "b"), paste0("c", 1:5))),
        spatial_locs = data.frame(sdimx = c(10, 20, 30, 40, 90), sdimy = 0,
            cell_ID = paste0("c", 1:5)),
        verbose = FALSE
    )
    g <- setPolygonInfo(g, gp, verbose = FALSE, initialize = FALSE)
    ids <- paste0("c", 1:5)
    e_mem <- .icf_touching_edges(g, "cell", ids)
    gp_disk <- gp
    gp_disk[] <- GiottoDisk::storeWrite(GiottoDisk::parquetGeomStore(), sv)
    g <- setPolygonInfo(g, gp_disk, verbose = FALSE, initialize = FALSE)
    e_disk <- .icf_touching_edges(g, "cell", ids)
    expect_equal(e_disk, e_mem)
    expect_setequal(paste(e_mem$recv, e_mem$send),
        c("1 2", "2 1", "2 3", "3 2", "3 4", "4 3"))
})

test_that("plotICF draws an edgeicg icfObject", {
    skip_if_not_installed("edgeicg")
    s <- .icf_sim_giotto()
    res <- findInteractionChangedFeats(s$g,
        cluster_column = "ct", method = "edgeicg",
        min_shift = 150, B = 8, report_q = 1
    )
    sc <- res$ICFscores[scale == "contact" & cell_type == "t1"]
    feats <- head(unique(sc$feats), 3)
    pl <- plotICF(s$g, res,
        source_type = "t1", source_markers = feats[1],
        ICF_feats = stats::setNames(feats[-1], rep("t2", 2)),
        return_plot = TRUE, show_plot = FALSE, save_plot = FALSE
    )
    expect_s3_class(pl, "ggplot")
})

test_that("the edgeicg plots draw from an icfObject", {
    skip_if_not_installed("edgeicg")
    s <- .icf_sim_giotto()
    res <- findInteractionChangedFeats(s$g,
        cluster_column = "ct", method = "edgeicg",
        min_shift = 150, B = 8, report_q = 1
    )
    out <- list(return_plot = TRUE, show_plot = FALSE, save_plot = FALSE)
    tr <- s$sim$truth[receiver == "t1"]
    plots <- list(
        do.call(plotICFOverview, c(list(s$g, res), out)),
        do.call(plotICFLeak, c(list(s$g, res), out)),
        do.call(plotICFPair, c(list(s$g, res, receiver = "t1",
            sender = "t2"), out)),
        do.call(plotICFExposure, c(list(s$g, res, receiver = "t1",
            sender = "t2", feats = head(tr$gene, 2), scale = "contact"), out))
    )
    for (p in plots) expect_s3_class(p, "ggplot")
    at <- data.table::as.data.table(res$edgeicg$attribution)
    if (nrow(at)) {
        r1 <- res$edgeicg$types[at$receiver[1]]
        p <- do.call(plotICFAttribution, c(list(s$g, res, receiver = r1,
            feats = at$gene[1]), out))
        expect_s3_class(p, "ggplot")
    }
    expect_error(plotICFPair(s$g, res, receiver = "t1", sender = "nope"),
        "no tests")
})

test_that("min_lfc floors ICFscores and spatial_covariates reach the engine", {
    skip_if_not_installed("edgeicg")
    s <- .icf_sim_giotto()
    run <- function(...) findInteractionChangedFeats(s$g,
        cluster_column = "ct", method = "edgeicg", min_shift = 150, B = 8,
        report_q = 0.05, ...)
    a <- run()
    expect_true(all(abs(a$ICFscores$log2fc) >= 0.25))
    expect_identical(a$test_info$min_lfc, 0.25)
    b <- run(min_lfc = 0)
    expect_gte(nrow(b$ICFscores), nrow(a$ICFscores))
    c <- run(spatial_covariates = TRUE, spatial_df = 3)
    expect_gt(c$edgeicg$settings$covariates, 0)
    expect_identical(c$edgeicg$settings$covariate_scope, "niche")
})

test_that("classic downstream functions refuse an edgeicg result clearly", {
    skip_if_not_installed("edgeicg")
    s <- .icf_sim_giotto()
    res <- findInteractionChangedFeats(s$g,
        cluster_column = "ct", min_shift = 150, B = 4,
        selected_feats = rownames(s$sim$counts)[1:10]
    )
    expect_error(filterICF(res), "classic test")
    expect_error(combineICF(res), "classic test")
    expect_error(plotCellProximityFeats(s$g, res), "classic test")
})

test_that("region_column gives region effects, $ICFregion and plotICFRegion()", {
    skip_if_not_installed("edgeicg", "0.2.0")
    sim <- edgeicg::simulate_edgeicg(n = 3000, G = 160, seed = 3, region_gamma = 1.2)
    ids <- colnames(sim$counts)
    g <- createGiottoObject(
        expression = sim$counts,
        spatial_locs = data.frame(
            sdimx = sim$coords[, 1], sdimy = sim$coords[, 2], cell_ID = ids
        ),
        verbose = FALSE
    )
    g <- addCellMetadata(g,
        new_metadata = data.frame(cell_ID = ids, ct = sim$types, half = sim$region),
        by_column = TRUE, column_cell_ID = "cell_ID"
    )
    res <- findInteractionChangedFeats(g, cluster_column = "ct",
        region_column = "half", min_shift = 150, B = 4)
    expect_true(data.table::is.data.table(res$ICFregion))
    expect_true(all(c("feats", "cell_type", "int_cell_type", "region",
        "effect", "q", "q_het", "pattern") %in% names(res$ICFregion)))
    expect_setequal(unique(res$ICFregion$region), c("left", "right"))
    f <- head(sim$truth[region == "right", gene], 2)
    pl <- plotICFRegion(g, res, receiver = "t1", sender = "t2", feats = f,
        scale = "contact", return_plot = TRUE, show_plot = FALSE,
        save_plot = FALSE)
    expect_s3_class(pl, "ggplot")
    expect_error(
        findInteractionChangedFeats(g, cluster_column = "ct",
            region_column = "nope", min_shift = 150, B = 4),
        "not a cell metadata column"
    )
    res0 <- findInteractionChangedFeats(g, cluster_column = "ct",
        min_shift = 150, B = 4)
    expect_null(res0$ICFregion)
    expect_error(plotICFRegion(g, res0, "t1", "t2", f), "no region effects")
})

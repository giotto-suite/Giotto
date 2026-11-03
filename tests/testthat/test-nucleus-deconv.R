# Tests for runNucleusDeconv() and hierarchicalCellTypes(): the gobject
# plumbing around the bag model (test-bag-model.R) and spatialDWLS's
# expression loss (tested in spatialDWLS). Synthetic and small: 100 spots on a
# grid, 3-8 cells per spot plus some between spots, features that separate
# the three cell types.

.nd_fixture <- function(seed = 1L) {
    set.seed(seed)
    types <- c("A", "B", "C")
    grid <- expand.grid(i = 0:9, j = 0:9)
    spots <- data.table::data.table(cell_ID = sprintf("s%02d", seq_len(nrow(grid))),
        sdimx = grid$i * 100, sdimy = grid$j * 100)
    n <- sample(3:8, nrow(spots), replace = TRUE)
    inside <- data.table::data.table(
        spot = rep(spots$cell_ID, n),
        sdimx = rep(spots$sdimx, n) + stats::runif(sum(n), -15, 15),
        sdimy = rep(spots$sdimy, n) + stats::runif(sum(n), -15, 15))
    between <- data.table::data.table(spot = NA_character_,
        sdimx = stats::runif(30, 0, 900), sdimy = stats::runif(30, 0, 900))
    between <- between[(sdimx %% 100 > 40 & sdimx %% 100 < 60)]
    nuc <- rbind(inside, between)
    nuc[, cell_ID := sprintf("n%04d", .I)]
    nuc[, type := sample(types, .N, replace = TRUE, prob = c(0.5, 0.3, 0.2))]

    genes <- paste0("g", 1:30)
    sig <- matrix(1, 30, 3, dimnames = list(genes, types))
    sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
    tab <- table(factor(nuc$spot, spots$cell_ID), factor(nuc$type, types))
    counts <- sig %*% t(matrix(tab, nrow(tab))) * 2
    counts[] <- stats::rpois(length(counts), counts)
    dimnames(counts) <- list(genes, spots$cell_ID)
    prop <- tab / rowSums(tab)

    centre <- matrix(stats::rnorm(3 * 10, sd = 1.5), 3, 10,
        dimnames = list(types, NULL))
    feats <- t(centre[nuc$type, ] +
        matrix(stats::rnorm(nrow(nuc) * 10), nrow(nuc), 10))
    dimnames(feats) <- list(paste0("f", 1:10), nuc$cell_ID)

    g <- gwith_options(list(giotto.use_conda = FALSE, giotto.verbose = FALSE,
            giotto.no_python_warn = TRUE, giotto.has_conda = FALSE), {
        g <- createGiottoObject(expression = counts,
            spatial_locs = spots[, .(sdimx, sdimy, cell_ID)], verbose = FALSE)
        g <- setGiotto(g, createExprObj(log2(counts + 1), name = "normalized",
            spat_unit = "cell", feat_type = "rna"), verbose = FALSE)
        g <- addCellMetadata(g, new_metadata = data.frame(
            cell_ID = spots$cell_ID, leiden_clus = rep(1:2, 50)),
            by_column = TRUE, column_cell_ID = "cell_ID")
        g <- setGiotto(g, createGiottoPolygon(polyStamp(
            circleVertices(27.5, 30), spots, verbose = FALSE), name = "cell",
            verbose = FALSE), verbose = FALSE)
        # starting proportions: the truth blurred towards uniform
        start <- 0.6 * prop + 0.4 / 3
        g <- setGiotto(g, createSpatEnrObj(data.table::data.table(
            cell_ID = spots$cell_ID, matrix(start, nrow(start),
                dimnames = list(NULL, types))),
            name = "DWLS", spat_unit = "cell", feat_type = "rna",
            method = "DWLS", verbose = FALSE), verbose = FALSE)
        g <- setGiotto(g, createGiottoPolygon(polyStamp(circleVertices(4, 12),
            nuc[, .(sdimx, sdimy, cell_ID)], verbose = FALSE),
            name = "nucleus", verbose = FALSE), verbose = FALSE)
        g <- addSpatialCentroidLocations(g, poly_info = "nucleus",
            verbose = FALSE)
        setGiotto(g, createExprObj(feats, name = "emb", spat_unit = "nucleus",
            feat_type = "img"), verbose = FALSE)
    })
    list(g = g, sig = sig, nuc = nuc, types = types)
}

skip_if_no_dwls_loss <- function() {
    skip_if_not_installed("spatialDWLS")
    skip_if_not(exists("dwls_design", asNamespace("spatialDWLS")),
        "spatialDWLS >= 0.1.2 needed")
}

test_that("runNucleusDeconv() types the cells and writes results to their spatial unit", {
    skip_if_no_dwls_loss()
    f <- .nd_fixture()
    g <- runNucleusDeconv(f$g, spat_unit = "nucleus", spot_unit = "cell",
        sign_matrix = f$sig, epochs = 60L, lr = 1e-2, batch_size = 16L,
        verbose = FALSE)
    enr <- getSpatialEnrichment(g, spat_unit = "nucleus", feat_type = "img",
        name = "nucleusDeconv", output = "spatEnrObj")
    expect_identical(sort(enr[]$cell_ID), sort(f$nuc$cell_ID))
    expect_setequal(colnames(enr[]), c("cell_ID", f$types))
    expect_equal(unname(rowSums(as.matrix(enr[][, f$types, with = FALSE]))),
        rep(1, nrow(f$nuc)))
    expect_type(enr@misc$model$theta, "list")

    md <- pDataDT(g, spat_unit = "nucleus", feat_type = "img")
    md <- md[match(f$nuc$cell_ID, cell_ID)]
    expect_gt(mean(md$nucleusDeconv_type == f$nuc$type), 0.85)
    # between-spot cells are typed too, and have no spot
    expect_identical(md$nucleusDeconv_spot, f$nuc$spot)
    expect_gt(mean(md[is.na(nucleusDeconv_spot)]$nucleusDeconv_type ==
        f$nuc[is.na(spot)]$type), 0.7)
    expect_true(all(md$nucleusDeconv_confidence > 0 &
        md$nucleusDeconv_confidence <= 1))
    expect_match(tail(names(g@parameters), 1), "_nucleus_deconvolution$")
    # the spot unit is untouched
    expect_identical(pDataDT(g)$leiden_clus, pDataDT(f$g)$leiden_clus)
})

test_that("runNucleusDeconv() returns the enrichment object and checks its inputs", {
    skip_if_no_dwls_loss()
    f <- .nd_fixture()
    enr <- runNucleusDeconv(f$g, spat_unit = "nucleus", spot_unit = "cell",
        sign_matrix = f$sig, epochs = 2L, return_gobject = FALSE,
        name = "nd2", verbose = FALSE)
    expect_s4_class(enr, "spatEnrObj")
    expect_identical(objName(enr), "nd2")
    expect_identical(spatUnit(enr), "nucleus")
    expect_error(runNucleusDeconv(f$g, spat_unit = "nucleus",
        spot_unit = "cell", sign_matrix = f$sig, cluster_column = "nope",
        epochs = 1L), "cluster column")
    bad <- f$sig
    colnames(bad)[3] <- "D"
    expect_error(runNucleusDeconv(f$g, spat_unit = "nucleus",
        spot_unit = "cell", sign_matrix = bad, epochs = 1L), "no column")
    expect_error(runNucleusDeconv(f$g, spat_unit = "cell", spot_unit = "cell",
        sign_matrix = f$sig), "same spatial unit")
    expect_error(runNucleusDeconv(f$g, spat_unit = "nucleus",
        spot_unit = "cell", sign_matrix = f$sig, deconv_name = NULL,
        loss_weights = c(dwls = 0)), "no loss")
})

test_that("hierarchicalCellTypes() reports the lineage and a confident subtype", {
    P <- rbind(c(0.8, 0.1, 0.1), c(0.45, 0.40, 0.15), c(0.1, 0.2, 0.7))
    colnames(P) <- c("A1", "A2", "B")
    r <- .hier_read_out(P, c(A1 = "A", A2 = "A"), tau = 0.7)
    expect_identical(r$lineage, c("A", "A", "B"))
    expect_identical(r$subtype, c("A1", NA, "B"))
    expect_identical(r$label, c("A1", "A", "B"))
    expect_equal(r$share, c(0.8 / 0.9, 0.45 / 0.85, 1))

    skip_if_no_dwls_loss()
    f <- .nd_fixture()
    g <- runNucleusDeconv(f$g, spat_unit = "nucleus", spot_unit = "cell",
        sign_matrix = f$sig, epochs = 2L, verbose = FALSE)
    g <- hierarchicalCellTypes(g, spat_unit = "nucleus",
        lineage_map = c(A = "AB", B = "AB"), tau = 0.9)
    md <- pDataDT(g, spat_unit = "nucleus", feat_type = "img")
    expect_true(all(c("cell_type_hier", "cell_type_hier_lineage",
        "cell_type_hier_subtype", "cell_type_hier_share") %in% colnames(md)))
    expect_true(all(md$cell_type_hier_lineage %in% c("AB", "C")))
})

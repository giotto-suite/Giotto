# Tests for Giotto's DWLS wrapper.
#
# The algorithm is the spatialDWLS package (drieslab/spatialDWLS); its own
# tests cover the solver on synthetic data (single-type clusters, NA spots,
# RNG state, parallel = serial, sparse = dense, windows). What is tested here
# is Giotto's side: the param, the verb, the gobject plumbing and the
# functions that read the stored result.
#
# The pinned values were re-pinned when the solver was rewritten. The rewrite
# reproduces the previous solver to ~1e-14 on this fixture and on the whole
# mini object. The pins changed only because the earlier ones predate
# f8fa23b5, which keyed the clusters to the expression columns: on this subset
# 118 of 120 cells had been solved against another cell's cluster.
#
# Deliberately small: 120 cells and 100 signature genes, about half a second.

skip_if_no_mini <- function() {
    skip_if_not_installed("GiottoData")
    skip_if_not_installed("spatialDWLS")
}

.deconv_fixture <- function(n_cells = 120L, n_genes = 100L, seed = 1L) {
    gwith_options(list(giotto.use_conda = FALSE, giotto.verbose = FALSE,
                       giotto.no_python_warn = TRUE,
                       giotto.has_conda = FALSE), {
        g0 <- suppressMessages(
            GiottoData::loadGiottoMini("visium", python_path = NA)
        )
        g <- GiottoClass::subsetGiotto(
            g0, cell_ids = pDataDT(g0)$cell_ID[seq_len(n_cells)]
        )
        md <- pDataDT(g)
        set.seed(seed)
        genes <- rownames(GiottoClass::getExpression(g, values = "normalized",
                                                     output = "matrix"))
        sm <- makeSignMatrixDWLS(g,
            sign_gene = sample(genes, n_genes),
            cell_type_vector = as.character(md$leiden_clus)
        )
        list(g = g, sm = sm, md = md)
    })
}


test_that("makeSignMatrixDWLS builds a mean-expression reference", {
    skip_if_no_mini()
    f <- .deconv_fixture()

    # NOT a binary membership matrix -- that is makeSignMatrixPAGE(). This one
    # is mean expression per cell type, which is why the two families document
    # `sign_matrix` separately rather than sharing one description.
    expect_true(is.matrix(f$sm))
    expect_identical(dim(f$sm), c(100L, 3L))
    expect_setequal(colnames(f$sm), c("1", "2", "5"))
    expect_false(all(f$sm %in% c(0, 1)))
    expect_gt(max(f$sm), 1)
    expect_equal(head(f$sm[, 1], 3),
                 c(17.037257, 12.750889, 15.206172),
                 tolerance = 1e-5, ignore_attr = TRUE)
})


test_that("runDWLSDeconv output is pinned", {
    skip_if_no_mini()
    f <- .deconv_fixture()

    res <- runDWLSDeconv(f$g, sign_matrix = f$sm,
                         cluster_column = "leiden_clus",
                         return_gobject = FALSE)
    expect_s4_class(res, "spatEnrObj")
    expect_identical(res@method, "DWLS")
    expect_identical(res@name, "DWLS")
    expect_setequal(names(res@misc), c(
        "expr_values_used", "logbase", "cluster_column_used",
        "number_of_cells_per_spot", "used_cut_off"
    ))

    dt <- res[]
    expect_identical(nrow(dt), 120L)
    expect_setequal(names(dt), c("cell_ID", "1", "2", "5"))
    expect_equal(head(dt[["2"]], 3), c(0.128027, 0.082551, 1.000000),
                 tolerance = 1e-5)
    expect_equal(head(dt[["5"]], 3), c(0.047666, 0.615664, 0.000000),
                 tolerance = 1e-5)
    expect_equal(head(dt[["1"]], 3), c(0.824307, 0.301786, 0.000000),
                 tolerance = 1e-5)
})


test_that("DWLS returns proportions, not scores", {
    skip_if_no_mini()
    f <- .deconv_fixture()

    # The contract that separates this family from enrichment: every cell's
    # values are a composition. Enrichment scores are unbounded and signed.
    dt <- runDWLSDeconv(f$g, sign_matrix = f$sm,
                        cluster_column = "leiden_clus",
                        return_gobject = FALSE)[]
    prop <- as.matrix(dt[, setdiff(names(dt), "cell_ID"), with = FALSE])

    # Bounded to within floating point, not exactly. The constrained solve
    # leaves a few values a hair below zero -- measured 3 of 360 here, the
    # largest -6e-17. A validator written as a hard `>= 0` would reject real
    # DWLS output, so the tolerance is part of the contract, not a concession.
    expect_gt(min(prop), -1e-12)
    expect_lte(max(prop), 1)
    expect_equal(unname(rowSums(prop)), rep(1, nrow(prop)), tolerance = 1e-8)
})


test_that("the deconvolution result lands in the gobject", {
    skip_if_no_mini()
    f <- .deconv_fixture()

    g2 <- runDWLSDeconv(f$g, sign_matrix = f$sm,
                        cluster_column = "leiden_clus")
    expect_s4_class(g2, "giotto")
    expect_true("DWLS" %in% GiottoClass::list_spatial_enrichments_names(
        g2, spat_unit = "cell", feat_type = "rna"))

    stored <- GiottoClass::getSpatialEnrichment(g2, name = "DWLS",
                                                output = "spatEnrObj")
    direct <- runDWLSDeconv(f$g, sign_matrix = f$sm,
                            cluster_column = "leiden_clus",
                            return_gobject = FALSE)
    expect_equal(stored[], direct[])

    # parameter history is appended under its own key, not the enrichment one
    expect_gt(length(g2@parameters), length(f$g@parameters))
    last_key <- names(g2@parameters)[length(g2@parameters)]
    expect_match(last_key, "_spatial_deconvolution$")
    expect_identical(
        unname(g2@parameters[[length(g2@parameters)]][["method used"]]), "DWLS")
})


test_that("runSpatialDeconv routes to the same result as the direct call", {
    skip_if_no_mini()
    f <- .deconv_fixture()

    via <- runSpatialDeconv(f$g, deconv_method = "DWLS", sign_matrix = f$sm,
                            cluster_column = "leiden_clus",
                            return_gobject = FALSE)
    direct <- runDWLSDeconv(f$g, sign_matrix = f$sm,
                            cluster_column = "leiden_clus",
                            return_gobject = FALSE)
    expect_equal(via[], direct[])
})


test_that("analyzeData() on a matrix matches the giotto wrapper", {
    skip_if_no_mini()
    f <- .deconv_fixture()

    expr <- GiottoClass::getExpression(f$g, values = "normalized",
                                       output = "matrix")
    clus <- stats::setNames(f$md$leiden_clus, f$md$cell_ID)
    verb <- analyzeData(expr, deconvParam("DWLS"),
        sign_matrix = f$sm, cell_clusters = clus)
    wrapped <- runDWLSDeconv(f$g, sign_matrix = f$sm,
                             cluster_column = "leiden_clus",
                             return_gobject = FALSE)[]
    data.table::setcolorder(verb, names(wrapped))
    expect_equal(verb, wrapped, ignore_attr = TRUE)
})


test_that("the param reaches the engine, and giotto.bpparam is honoured", {
    skip_if_not_installed("spatialDWLS")
    skip_if_not_installed("BiocParallel")

    seen <- NULL
    local_mocked_bindings(
        dwls = function(expr, sign_matrix, clusters, ...) {
            seen <<- list(...)
            prop <- matrix(1 / 2, 2, ncol(expr),
                dimnames = list(c("A", "B"), colnames(expr)))
            prop
        },
        .package = "spatialDWLS"
    )
    expr <- matrix(1, 3, 4, dimnames = list(paste0("g", 1:3),
        paste0("s", 1:4)))
    sig <- matrix(1, 3, 2, dimnames = list(paste0("g", 1:3), c("A", "B")))

    bp <- BiocParallel::SnowParam(2)
    old <- options(giotto.bpparam = bp)
    on.exit(options(old), add = TRUE)
    res <- analyzeData(expr, deconvParam("DWLS", n_cell = 7, cutoff = 1.5,
        logbase = 10, na_as_zero = TRUE, chunk_size = 9L),
        sign_matrix = sig, cell_clusters = rep("k", 4))

    expect_identical(seen$n_cell, 7)
    expect_identical(seen$cutoff, 1.5)
    expect_identical(seen$logbase, 10)
    expect_true(seen$na_as_zero)
    expect_identical(seen$chunk_size, 9L)
    expect_identical(seen$BPPARAM, bp)

    # an explicit BPPARAM wins over the session option
    serial <- BiocParallel::SerialParam()
    analyzeData(expr, deconvParam("DWLS", BPPARAM = serial),
        sign_matrix = sig, cell_clusters = rep("k", 4))
    expect_identical(seen$BPPARAM, serial)

    # the engine's matrix comes back as the enrichment contract
    expect_s3_class(res, "data.table")
    expect_identical(names(res), c("cell_ID", "A", "B"))
    expect_identical(res$cell_ID, paste0("s", 1:4))
})


test_that("deconvParam validates its inputs", {
    expect_s4_class(deconvParam(), "dwlsDeconvParam")
    expect_error(deconvParam("DWLS", BPPARAM = "two"), "BiocParallelParam")
    expect_identical(deconvParam("dwls", n_cell = 10)$n_cell, 10)
    expect_error(deconvParam("DWLS", n_cell = 0), "n_cell")
    expect_error(deconvParam("DWLS", n_cell = c(5, 10)), "n_cell")
    expect_error(deconvParam("RCTD"))
})


test_that("the spot interaction functions read the stored result", {
    skip_if_no_mini()
    f <- .deconv_fixture()
    g <- runDWLSDeconv(f$g, sign_matrix = f$sm, cluster_column = "leiden_clus")

    # cellProximityEnrichmentSpots() used to drop the first column, assuming
    # it was cell_ID; the stored enrichment keeps cell_ID last.
    res <- cellProximityEnrichmentSpots(g, cells_in_spot = 5,
        number_of_simulations = 5)
    expect_s3_class(res$enrichm_res, "data.table")
    expect_gt(nrow(res$enrichm_res), 0L)

    pred <- featExpDWLS(g, ave_celltype_exp = f$sm)
    expect_identical(dim(pred), c(nrow(f$sm), 120L))
})

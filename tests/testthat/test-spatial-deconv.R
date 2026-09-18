# Characterization tests for DWLS deconvolution.
#
# These pin what runDWLSDeconv() and makeSignMatrixDWLS() produce *today*,
# before deconvolution is split into its own param family. Nothing here asserts
# the numbers are right -- only that they do not change. There was no test file
# for DWLS at all, and the enrichment round showed what that costs: four live
# bugs sat in the sibling functions unnoticed, one of which meant
# runRankEnrich() could not be called.
#
# Deliberately small: 120 cells and 100 signature genes runs in ~8 s. DWLS
# optimizes per cell, so the full mini object takes ~20 s for no extra coverage.

skip_if_no_mini <- function() skip_if_not_installed("GiottoData")

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


test_that("runDWLSDeconv output is unchanged", {
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
    expect_equal(head(dt[["2"]], 3), c(0.976896, 0.000000, 1.000000),
                 tolerance = 1e-5)
    expect_equal(head(dt[["5"]], 3), c(0.023104, 0.854972, 0.000000),
                 tolerance = 1e-5)
    expect_equal(head(dt[["1"]], 3), c(0.000000, 0.145028, 0.000000),
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

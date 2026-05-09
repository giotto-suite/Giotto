# =============================================================================
# Tests for createGiottoVisiumHDObject{,Bin,Cell}(expression_backend = "parquet")
# =============================================================================
# Verifies:
#   1. default backend = "matrix" preserves dgCMatrix expression bit-for-bit
#   2. backend = "parquet" returns a parquetExprStore-backed exprObj
#   3. dim, feat set, and cell set match across backends — bin + cell
# Skips automatically when GiottoDisk or the local fixture isn't present.
# =============================================================================

VHD_DIR <- "/Users/rubendries/Documents/Datasets/visiumHD/Visium_HD_Tiny_3prime_Dataset_outs"

skip_if_no_fixture <- function() {
    if (!dir.exists(VHD_DIR) ||
        !dir.exists(file.path(VHD_DIR, "binned_outputs")) ||
        !dir.exists(file.path(VHD_DIR, "segmented_outputs"))) {
        skip("VisiumHD fixture not available")
    }
}

skip_if_no_giottodisk <- function() {
    if (!requireNamespace("GiottoDisk", quietly = TRUE)) {
        skip("GiottoDisk not installed (Suggests)")
    }
}


# ---- BIN -------------------------------------------------------------------

bin_args <- function(backend) {
    list(
        binned_outputs_dir       = file.path(VHD_DIR, "binned_outputs"),
        bin                      = 8,
        expression_source        = "filtered",
        load_image               = FALSE,
        load_metadata            = FALSE,
        load_spatlocs            = FALSE,
        create_tessellated_polys = FALSE,
        expression_backend       = backend,
        verbose                  = FALSE
    )
}

test_that("bin reader: default backend = matrix preserves dgCMatrix expr", {
    skip_if_no_fixture()
    g <- suppressWarnings(do.call(createGiottoVisiumHDObjectBin,
                                    bin_args("matrix")))
    em <- GiottoClass::getExpression(g, output = "matrix")
    expect_true(inherits(em, "Matrix"))
    expect_false(inherits(em, "parquetExprStore"))
})

test_that("bin reader: backend = parquet returns parquetExprStore exprObj", {
    skip_if_no_fixture(); skip_if_no_giottodisk()
    g <- suppressWarnings(do.call(createGiottoVisiumHDObjectBin,
                                    bin_args("parquet")))
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    expect_s4_class(eo, "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
})

test_that("bin reader: matrix and parquet backends produce identical dim/feat/cell sets", {
    skip_if_no_fixture(); skip_if_no_giottodisk()
    g_mem <- suppressWarnings(do.call(createGiottoVisiumHDObjectBin,
                                        bin_args("matrix")))
    g_pq  <- suppressWarnings(do.call(createGiottoVisiumHDObjectBin,
                                        bin_args("parquet")))
    em <- GiottoClass::getExpression(g_mem, output = "matrix")
    pe <- slot(GiottoClass::getExpression(g_pq, output = "exprObj"),
                "exprMat")
    expect_equal(nrow(em), nrow(pe))
    expect_equal(ncol(em), ncol(pe))
    expect_setequal(rownames(em), pe@feat_ids)
    expect_setequal(colnames(em), pe@cell_ids)
})


# ---- CELL ------------------------------------------------------------------

cell_args <- function(backend) {
    list(
        segmented_outputs_dir = file.path(VHD_DIR, "segmented_outputs"),
        expression_source     = "filtered",
        load_image            = FALSE,
        load_polygons         = NULL,
        expression_backend    = backend,
        verbose               = FALSE
    )
}

test_that("cell reader: default backend = matrix preserves dgCMatrix expr", {
    skip_if_no_fixture()
    g <- suppressWarnings(do.call(createGiottoVisiumHDObjectCell,
                                    cell_args("matrix")))
    em <- GiottoClass::getExpression(g, output = "matrix")
    expect_true(inherits(em, "Matrix"))
    expect_false(inherits(em, "parquetExprStore"))
})

test_that("cell reader: backend = parquet returns parquetExprStore exprObj", {
    skip_if_no_fixture(); skip_if_no_giottodisk()
    g <- suppressWarnings(do.call(createGiottoVisiumHDObjectCell,
                                    cell_args("parquet")))
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    expect_s4_class(eo, "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
})

test_that("cell reader: matrix and parquet backends produce identical dim/feat/cell sets", {
    skip_if_no_fixture(); skip_if_no_giottodisk()
    g_mem <- suppressWarnings(do.call(createGiottoVisiumHDObjectCell,
                                        cell_args("matrix")))
    g_pq  <- suppressWarnings(do.call(createGiottoVisiumHDObjectCell,
                                        cell_args("parquet")))
    em <- GiottoClass::getExpression(g_mem, output = "matrix")
    pe <- slot(GiottoClass::getExpression(g_pq, output = "exprObj"),
                "exprMat")
    expect_equal(nrow(em), nrow(pe))
    expect_equal(ncol(em), ncol(pe))
    expect_setequal(rownames(em), pe@feat_ids)
    expect_setequal(colnames(em), pe@cell_ids)
})

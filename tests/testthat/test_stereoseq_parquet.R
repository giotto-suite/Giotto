# =============================================================================
# Tests for createGiottoStereoSeqObject{Bin,Cell}(expression_backend = "parquet")
# =============================================================================
# Mirrors test_xenium_parquet.R: verifies that the new expression_backend
# argument:
#   1. preserves bit-for-bit current behavior when backend = "matrix"
#   2. produces a parquetExprStore-backed exprObj when backend = "parquet"
#   3. yields identical dim, feat set, and cell set across backends
#       — for both bin and cellbin readers
# Skips automatically when GiottoDisk, rhdf5, or the local fixture isn't
# present (so CI doesn't fail).
# =============================================================================

STEREOSEQ_DIR <- "/Users/rubendries/Documents/Datasets/stereoseq/C04687E314_backup/outs"

skip_if_no_fixture <- function() {
    if (!dir.exists(STEREOSEQ_DIR) ||
        !dir.exists(file.path(STEREOSEQ_DIR, "feature_expression"))) {
        skip("Stereo-seq fixture not available")
    }
}

skip_if_no_giottodisk <- function() {
    if (!requireNamespace("GiottoDisk", quietly = TRUE)) {
        skip("GiottoDisk not installed (Suggests)")
    }
}

skip_if_no_rhdf5 <- function() {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        skip("rhdf5 not installed (Suggests)")
    }
}


# ---- BIN -------------------------------------------------------------------

bin_args <- function(backend) {
    list(
        stereoseq_dir      = STEREOSEQ_DIR,
        bin_size           = "bin100",
        load_image         = FALSE,
        load_mask          = FALSE,
        load_spatlocs      = FALSE,
        expression_backend = backend,
        verbose            = FALSE
    )
}

test_that("bin reader: default backend = matrix preserves dgCMatrix expr", {
    skip_if_no_fixture(); skip_if_no_rhdf5()
    g <- suppressWarnings(do.call(createGiottoStereoSeqObjectBin,
                                    bin_args("matrix")))
    em <- GiottoClass::getExpression(g, output = "matrix")
    expect_true(inherits(em, "Matrix"))
    expect_false(inherits(em, "parquetExprStore"))
})

test_that("bin reader: backend = parquet returns parquetExprStore exprObj", {
    skip_if_no_fixture(); skip_if_no_rhdf5(); skip_if_no_giottodisk()
    g <- suppressWarnings(do.call(createGiottoStereoSeqObjectBin,
                                    bin_args("parquet")))
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    expect_s4_class(eo, "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
})

test_that("bin reader: matrix and parquet backends produce identical dim/feat/cell sets", {
    skip_if_no_fixture(); skip_if_no_rhdf5(); skip_if_no_giottodisk()
    g_mem <- suppressWarnings(do.call(createGiottoStereoSeqObjectBin,
                                        bin_args("matrix")))
    g_pq  <- suppressWarnings(do.call(createGiottoStereoSeqObjectBin,
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
        stereoseq_dir      = STEREOSEQ_DIR,
        load_image         = FALSE,
        load_polygons      = FALSE,
        load_mask          = FALSE,
        load_spatlocs      = FALSE,
        expression_backend = backend,
        verbose            = FALSE
    )
}

test_that("cell reader: default backend = matrix preserves dgCMatrix expr", {
    skip_if_no_fixture(); skip_if_no_rhdf5()
    g <- suppressWarnings(do.call(createGiottoStereoSeqObjectCell,
                                    cell_args("matrix")))
    em <- GiottoClass::getExpression(g, output = "matrix")
    expect_true(inherits(em, "Matrix"))
    expect_false(inherits(em, "parquetExprStore"))
})

test_that("cell reader: backend = parquet returns parquetExprStore exprObj", {
    skip_if_no_fixture(); skip_if_no_rhdf5(); skip_if_no_giottodisk()
    g <- suppressWarnings(do.call(createGiottoStereoSeqObjectCell,
                                    cell_args("parquet")))
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    expect_s4_class(eo, "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
})

test_that("cell reader: matrix and parquet backends produce identical dim/feat/cell sets", {
    skip_if_no_fixture(); skip_if_no_rhdf5(); skip_if_no_giottodisk()
    g_mem <- suppressWarnings(do.call(createGiottoStereoSeqObjectCell,
                                        cell_args("matrix")))
    g_pq  <- suppressWarnings(do.call(createGiottoStereoSeqObjectCell,
                                        cell_args("parquet")))
    em <- GiottoClass::getExpression(g_mem, output = "matrix")
    pe <- slot(GiottoClass::getExpression(g_pq, output = "exprObj"),
                "exprMat")
    expect_equal(nrow(em), nrow(pe))
    expect_equal(ncol(em), ncol(pe))
    expect_setequal(rownames(em), pe@feat_ids)
    expect_setequal(colnames(em), pe@cell_ids)
})

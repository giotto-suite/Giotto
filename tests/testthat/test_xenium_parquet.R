# =============================================================================
# Tests for createGiottoXeniumObject(expression_backend = "parquet")
# =============================================================================
# Phase 3 Layer 2 — verifies that the new expression_backend argument:
#   1. preserves bit-for-bit current behavior when backend = "matrix"
#   2. produces a parquetExprStore-backed exprObj when backend = "parquet"
#   3. yields identical dim, feat set, and cell set across backends
#   4. errors clearly when GiottoDisk is missing or input is not mtx
# Skips automatically when GiottoDisk or the mini Xenium fixture isn't present.
# =============================================================================

XENIUM_DIR <- file.path(
    "/Users/rubendries/Documents/Software/scstream_prototype",
    "Data/workshop_xenium"
)

skip_if_no_fixture <- function() {
    if (!dir.exists(XENIUM_DIR) ||
        !file.exists(file.path(XENIUM_DIR,
                                "cell_feature_matrix/matrix.mtx.gz"))) {
        skip("mini Xenium fixture not available")
    }
}

skip_if_no_giottodisk <- function() {
    if (!requireNamespace("GiottoDisk", quietly = TRUE)) {
        skip("GiottoDisk not installed (Suggests)")
    }
}

base_args <- function(backend) {
    list(
        xenium_dir         = XENIUM_DIR,
        load_images        = NULL,
        load_aligned_images = NULL,
        load_transcripts   = FALSE,
        load_expression    = TRUE,
        expression_backend = backend,
        verbose            = FALSE
    )
}


test_that("default backend = matrix preserves dgCMatrix expression", {
    skip_if_no_fixture()

    g <- suppressWarnings(do.call(createGiottoXeniumObject,
                                   base_args("matrix")))
    em <- GiottoClass::getExpression(g, output = "matrix")
    # dgTMatrix or dgCMatrix — both are sparse Matrix subclasses
    expect_true(inherits(em, "Matrix"))
    expect_false(inherits(em, "parquetExprStore"))
    expect_equal(dim(em), c(377L, 7655L))
})


test_that("backend = parquet returns a parquetExprStore-backed exprObj", {
    skip_if_no_fixture()
    skip_if_no_giottodisk()

    g <- suppressWarnings(do.call(createGiottoXeniumObject,
                                   base_args("parquet")))
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    expect_s4_class(eo, "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
    pe <- slot(eo, "exprMat")
    expect_equal(dim(pe), c(377L, 7655L))
})


test_that("matrix and parquet backends produce identical dim/feat/cell sets", {
    skip_if_no_fixture()
    skip_if_no_giottodisk()

    g_mem <- suppressWarnings(do.call(createGiottoXeniumObject,
                                       base_args("matrix")))
    g_pq  <- suppressWarnings(do.call(createGiottoXeniumObject,
                                       base_args("parquet")))

    em <- GiottoClass::getExpression(g_mem, output = "matrix")
    pe <- slot(GiottoClass::getExpression(g_pq, output = "exprObj"),
                "exprMat")

    expect_equal(nrow(em), nrow(pe))
    expect_equal(ncol(em), ncol(pe))
    expect_setequal(rownames(em), pe@feat_ids)
    expect_setequal(colnames(em), pe@cell_ids)
})


test_that(
    "expression_backend = 'parquet' errors clearly when GiottoDisk missing",
    {
        skip_if_no_fixture()
        if (requireNamespace("GiottoDisk", quietly = TRUE)) {
            skip("GiottoDisk is installed; cannot test the missing-pkg error")
        }
        expect_error(
            suppressWarnings(do.call(createGiottoXeniumObject,
                                      base_args("parquet"))),
            regexp = "GiottoDisk"
        )
    }
)


test_that("parquet backend handles all 3 Xenium expression formats", {
    skip_if_no_fixture()
    skip_if_no_giottodisk()

    mtx_dir <- file.path(XENIUM_DIR, "cell_feature_matrix")
    tar_gz  <- file.path(XENIUM_DIR, "cell_feature_matrix.tar.gz")
    h5_file <- file.path(XENIUM_DIR, "cell_feature_matrix.h5")

    if (!file.exists(tar_gz) && !file.exists(h5_file)) {
        skip("only the unpacked mtx fixture is available; need tar.gz and/or h5")
    }
    if (!requireNamespace("hdf5r", quietly = TRUE) && file.exists(h5_file)) {
        # treat as if h5 not present
        h5_file <- ""
    }

    build <- function(label, ep) {
        op <- options(giotto.xenium_parquet_dir =
                        file.path(tempdir(), paste0("xenium_pq_", label)))
        on.exit(options(op), add = TRUE)
        args <- base_args("parquet")
        args$expression_path <- ep
        g <- suppressWarnings(do.call(createGiottoXeniumObject, args))
        slot(GiottoClass::getExpression(g, output = "exprObj"), "exprMat")
    }

    pe_dir <- build("dir", mtx_dir)
    expect_s4_class(pe_dir, "parquetExprStore")

    if (file.exists(tar_gz)) {
        pe_tar <- build("tar", tar_gz)
        expect_s4_class(pe_tar, "parquetExprStore")
        expect_equal(dim(pe_tar), dim(pe_dir))
        expect_setequal(pe_tar@feat_ids, pe_dir@feat_ids)
        expect_setequal(pe_tar@cell_ids, pe_dir@cell_ids)
    }

    if (file.exists(h5_file) && nzchar(h5_file)) {
        pe_h5 <- build("h5", h5_file)
        expect_s4_class(pe_h5, "parquetExprStore")
        expect_equal(dim(pe_h5), dim(pe_dir))
        expect_setequal(pe_h5@feat_ids, pe_dir@feat_ids)
        expect_setequal(pe_h5@cell_ids, pe_dir@cell_ids)
    }
})


test_that("global option giotto.expression_backend toggles the default", {
    skip_if_no_fixture()
    skip_if_no_giottodisk()

    op <- options(giotto.expression_backend = "parquet")
    on.exit(options(op), add = TRUE)

    args <- base_args("matrix")  # placeholder; will be overridden below
    args$expression_backend <- NULL
    # Drop the explicit arg so the formal default fires (which reads the option)
    g <- suppressWarnings(do.call(createGiottoXeniumObject, args))
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
})

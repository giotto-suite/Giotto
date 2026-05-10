# =============================================================================
# Tests for createGiottoCosMxObject(expression_backend = "parquet")
# =============================================================================
# Verifies (via the .cosmx_expression internal directly, since the full
# creator requires a polygon-fixture path that isn't always available):
#   1. default backend = "matrix" preserves dgCMatrix expression
#   2. backend = "parquet" returns parquetExprStore-backed exprObjs
#   3. dim, feat set, and cell set match across backends — per feat_type
#      after split_keyword
# Skips automatically when GiottoDisk or the local fixture is absent.
# =============================================================================

CSV_PATH <- file.path(
    "/Users/rubendries/Documents/Datasets/CosMx/Breast_cancer",
    "Flatfiles_RNA/flatFiles/BreastCancer",
    "BreastCancer_exprMat_file.csv.gz"
)

skip_if_no_fixture <- function() {
    if (!file.exists(CSV_PATH)) skip("CosMx exprMat fixture not available")
}

skip_if_no_giottodisk <- function() {
    if (!requireNamespace("GiottoDisk", quietly = TRUE)) {
        skip("GiottoDisk not installed (Suggests)")
    }
}

base_args <- function(backend) {
    list(
        path = CSV_PATH, slide = 1, fovs = NULL,
        feat_type     = c("rna", "negprobes"),
        split_keyword = list("SystemControl"),
        expression_backend = backend,
        verbose = FALSE
    )
}


test_that("default backend = matrix preserves dgCMatrix expression", {
    skip_if_no_fixture()
    exlist <- suppressWarnings(do.call(Giotto:::.cosmx_expression,
                                        base_args("matrix")))
    expect_true(length(exlist) >= 1L)
    em <- exlist[[1L]][]
    expect_true(inherits(em, "Matrix"))
    expect_false(inherits(em, "parquetExprStore"))
})

test_that("backend = parquet returns parquetExprStore-backed exprObjs", {
    skip_if_no_fixture(); skip_if_no_giottodisk()
    exlist <- suppressWarnings(do.call(Giotto:::.cosmx_expression,
                                        base_args("parquet")))
    expect_true(length(exlist) >= 1L)
    for (eo in exlist) {
        expect_s4_class(eo, "exprObj")
        expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
    }
})

test_that("matrix and parquet backends produce identical dim/feat/cell sets per feat_type", {
    skip_if_no_fixture(); skip_if_no_giottodisk()
    mat_list <- suppressWarnings(do.call(Giotto:::.cosmx_expression,
                                            base_args("matrix")))
    pq_list  <- suppressWarnings(do.call(Giotto:::.cosmx_expression,
                                            base_args("parquet")))

    expect_equal(length(mat_list), length(pq_list))

    by_ft <- function(lst) setNames(lst, vapply(lst,
        function(eo) slot(eo, "feat_type"), character(1)))
    mat_by_ft <- by_ft(mat_list)
    pq_by_ft  <- by_ft(pq_list)
    expect_setequal(names(mat_by_ft), names(pq_by_ft))

    for (ft in names(mat_by_ft)) {
        em <- mat_by_ft[[ft]][]
        pe <- slot(pq_by_ft[[ft]], "exprMat")
        expect_equal(nrow(em), nrow(pe), label = sprintf("[%s] nrow", ft))
        expect_equal(ncol(em), ncol(pe), label = sprintf("[%s] ncol", ft))
        expect_setequal(rownames(em), pe@feat_ids)
        expect_setequal(colnames(em), pe@cell_ids)
    }
})


# ---- csv_to_parquetExprStore unit-style round-trip ------------------------

test_that("csv_to_parquetExprStore round-trips a synthetic wide CSV", {
    skip_if_no_giottodisk()
    set.seed(1)
    n_cells <- 100L; n_genes <- 50L
    mat <- matrix(rpois(n_cells * n_genes, 0.5), nrow = n_cells)
    mat[rbinom(length(mat), 1, 0.9) == 1L] <- 0L
    df <- data.table::data.table(
        fov = sample(1:5, n_cells, replace = TRUE),
        cell_ID = seq_len(n_cells),
        mat
    )
    gene_names <- paste0("gene_", sprintf("%03d", seq_len(n_genes)))
    data.table::setnames(df, c("fov", "cell_ID", gene_names))
    csv_path <- tempfile(fileext = ".csv")
    data.table::fwrite(df, csv_path)

    out_dir <- tempfile(pattern = "csv_pq_")
    pe <- GiottoDisk::csv_to_parquetExprStore(
        csv_path = csv_path, output_path = out_dir,
        cell_id_col = "cell_ID", skip_cols = "fov",
        batch_rows = 30L, overwrite = TRUE
    )
    expect_s4_class(pe, "parquetExprStore")
    expect_equal(nrow(pe), n_genes)
    expect_equal(ncol(pe), n_cells)
    expect_equal(pe@feat_ids, gene_names)
    expect_equal(pe@cell_ids, as.character(df$cell_ID))
})

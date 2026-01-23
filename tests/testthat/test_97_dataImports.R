# python_path = NULL
# if(is.null(python_path)) {
#   installGiottoEnvironment()
# }



### TESTS FOR DATA IMPORT FUNCTIONS
# ------------------------------------

test_that("Expression matrix is read correctly", {
    # getSpatialDataset
    expect_no_error(suppressWarnings(GiottoData::getSpatialDataset(
        dataset = "scRNA_mouse_brain",
        directory = paste0(getwd(), "/testdata/")
    )))

    expect_true(file.exists(paste0(getwd(), "/testdata/brain_sc_expression_matrix.txt.gz")))

    # readExprMatrix
    expr_mat <- readExprMatrix(
        paste0(getwd(), "/testdata/brain_sc_expression_matrix.txt.gz"))

    expect_s4_class(expr_mat, "dgCMatrix")
    expect_equal(expr_mat@Dim, c(27998, 8039))

    # check a few genes
    expect_equal(expr_mat@Dimnames[[1]][20], "Sgcz")
    expect_equal(expr_mat@Dimnames[[1]][50], "Zfp804a")
})

# get10Xmatrix_h5
# TODO

# stitchFieldCoordinates
# TODO

# stitchTileCoordinates
# TODO

# -----------------------------
# remove files after testing
if (file.exists("./testdata/brain_sc_expression_matrix.txt.gz")) {
    unlink("./testdata/brain_sc_expression_matrix.txt.gz")
}

if (file.exists("./testdata/brain_sc_metadata.csv")) {
    unlink("./testdata/brain_sc_metadata.csv")
}

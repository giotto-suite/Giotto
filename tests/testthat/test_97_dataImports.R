# python_path = NULL
# if(is.null(python_path)) {
#   installGiottoEnvironment()
# }



### TESTS FOR DATA IMPORT FUNCTIONS
# ------------------------------------

test_that("Expression matrix is read correctly", {
    # getSpatialDataset returns the directory it wrote to, which is a
    # per-dataset subdirectory of `directory`
    expect_no_error(
        data_dir <- GiottoData::getSpatialDataset(
            dataset = "scRNA_mouse_brain",
            directory = file.path(getwd(), "testdata")
        )
    )

    expr_file <- file.path(data_dir, "brain_sc_expression_matrix.txt.gz")
    expect_true(file.exists(expr_file))

    # readExprMatrix
    expr_mat <- readExprMatrix(expr_file)

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
if (dir.exists("./testdata/scRNA_mouse_brain")) {
    unlink("./testdata/scRNA_mouse_brain", recursive = TRUE, force = TRUE)
}

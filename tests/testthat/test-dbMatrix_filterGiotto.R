# silence deprecated internal functions
rlang::local_options(lifecycle_verbosity = "quiet")

# ---------------------------------------------------------------------------- #
# Setup data
visium <- GiottoData::loadGiottoMini(dataset = "visium")
dgc <- getExpression(visium, output = "matrix")

con <- DBI::dbConnect(duckdb::duckdb(), ":memory:")

dbsm <- dbMatrix::dbMatrix(
    value = dgc,
    con = con,
    name = "dgc",
    class = "dbSparseMatrix",
    overwrite = TRUE
)

# Create exprObj with dbsm
expObj_db <- createExprObj(
    expression_data = dbsm,
    expression_matrix_class = "dbSparseMatrix",
    name = "raw"
)

# Create giotto object
gobject_db <- suppressWarnings(createGiottoObject(expression = expObj_db))

# ---------------------------------------------------------------------------- #
# Perform filtering
visium_filtered <- filterGiotto(visium,
    spat_unit = "cell",
    feat_type = "rna",
    expression_values = "raw"
)

gobject_db_filtered <- filterGiotto(gobject_db,
    spat_unit = "cell",
    feat_type = "rna",
    expression_values = "raw"
)

# Get filtered matrix
dgc_visium <- getExpression(visium_filtered, output = "matrix")
mat_db <- getExpression(gobject_db_filtered, output = "matrix")
dgc_db <- dbMatrix:::as_matrix(mat_db)

# ---------------------------------------------------------------------------- #
# Test filterGiotto() equivalence between dbMatrix and dgCMatrix

test_that("dbMatrix equivalent to dgCMatrix after filterGiotto()", {
    expect_equal(dgc_visium, dgc_db)
})

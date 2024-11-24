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

# ---------------------------------------------------------------------------- #
# Perform library normalization and scaling
visium_filtered <- normalizeGiotto(
    gobject = visium_filtered,
    spat_unit = "cell",
    feat_type = "rna",
    expression_values = "raw",
    library_size_norm = TRUE,
    log_norm = FALSE,
    scale_feats = FALSE,
    scale_cells = FALSE
)


gobject_db_filtered <- normalizeGiotto(
    gobject = gobject_db_filtered,
    spat_unit = "cell",
    feat_type = "rna",
    expression_values = "raw",
    library_size_norm = TRUE,
    log_norm = FALSE,
    scale_feats = FALSE,
    scale_cells = FALSE
)
# Get normalized matrix
dgc_visium <- getExpression(visium_filtered, output = "matrix", values = "normalized")
mat_db <- getExpression(gobject_db_filtered, output = "matrix", values = "normalized")
dgc_db <- dbMatrix:::as_matrix(mat_db)

# ---------------------------------------------------------------------------- #
# Test normalizeGiotto() equivalence between dbMatrix and dgCMatrix
test_that("dbMatrix equivalent to dgCMatrix after normalizeGiotto(library_size_norm = TRUE)", {
    expect_equal(dgc_visium, dgc_db)
})

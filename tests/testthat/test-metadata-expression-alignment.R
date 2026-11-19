# Verbs that pair a cell metadata column with expression columns have to key
# on cell_ID rather than assume the two arrive in the same order.
#
# Expression and cell metadata are fetched through independent accessors, and
# the object model neither guarantees a shared cell order nor restores one: a
# fresh object has them aligned, the setters preserve whatever order is there,
# and initialize() does not put a drifted metadata table back onto the
# expression axis. The shipped mini visium object is already in that state --
# its metadata agrees with the expression columns at zero of 624 positions.
#
# The property asserted below is the weakest one that catches positional
# pairing, and it needs no reference implementation: permuting the metadata
# row order must not change the result.

skip_if_no_mini <- function() skip_if_not_installed("GiottoData")

.mini <- function() {
    gwith_options(list(giotto.use_conda = FALSE, giotto.verbose = FALSE,
                       giotto.no_python_warn = TRUE,
                       giotto.has_conda = FALSE), {
        suppressMessages(GiottoData::loadGiottoMini("visium", python_path = NA))
    })
}

# Reorder @cell_metadata rows without touching anything else.
.permute_meta <- function(g, seed = 42L) {
    cm <- GiottoClass::getCellMetadata(g, output = "cellMetaObj",
        copy_obj = TRUE)
    set.seed(seed)
    cm[] <- cm[][sample.int(nrow(cm[])), ]
    GiottoClass::setCellMetadata(g, cm, verbose = FALSE, initialize = FALSE)
}


test_that("the fixture really does disagree on cell order", {
    skip_if_no_mini()
    # The premise of the rest of the file. If this ever fails, the fixture
    # changed and these tests stop covering anything.
    g <- .mini()
    ex <- GiottoClass::getExpression(g, values = "normalized",
        output = "matrix")
    ids <- GiottoClass::getCellMetadata(g, output = "data.table")[["cell_ID"]]
    expect_setequal(colnames(ex), ids)
    expect_false(identical(colnames(ex), ids))
})

test_that("adjustGiottoMatrix is invariant to metadata row order", {
    skip_if_no_mini()
    skip_if_not_installed("limma")
    g <- .mini()

    a <- suppressMessages(adjustGiottoMatrix(g, expression_values = "normalized",
        batch_columns = "leiden_clus", return_gobject = FALSE))
    b <- suppressMessages(adjustGiottoMatrix(.permute_meta(g),
        expression_values = "normalized",
        batch_columns = "leiden_clus", return_gobject = FALSE))

    expect_equal(as.matrix(a), as.matrix(b))
})

test_that("adjustGiottoMatrix corrects each cell against its own batch", {
    skip_if_no_mini()
    skip_if_not_installed("limma")
    g <- .mini()

    # Independent ground truth: call limma directly, with the batch vector
    # put into the matrix's column order by name.
    ex <- GiottoClass::getExpression(g, values = "normalized",
        output = "matrix")
    md <- GiottoClass::getCellMetadata(g, output = "data.table")
    batch <- md[["leiden_clus"]][match(colnames(ex), md[["cell_ID"]])]
    ref <- limma::removeBatchEffect(x = ex, batch = batch)

    got <- suppressMessages(adjustGiottoMatrix(g,
        expression_values = "normalized",
        batch_columns = "leiden_clus", return_gobject = FALSE))

    expect_equal(as.matrix(got), as.matrix(ref))
})

test_that("a metadata table missing cells is refused, not silently recycled", {
    skip_if_no_mini()
    skip_if_not_installed("limma")
    g <- .mini()
    cm <- GiottoClass::getCellMetadata(g, output = "cellMetaObj",
        copy_obj = TRUE)
    cm[] <- cm[][-1L, ]
    g2 <- GiottoClass::setCellMetadata(g, cm, verbose = FALSE,
        initialize = FALSE)

    expect_error(
        suppressMessages(adjustGiottoMatrix(g2,
            expression_values = "normalized",
            batch_columns = "leiden_clus", return_gobject = FALSE)),
        "do not all match"
    )
})

test_that("runGiottoHarmony is invariant to metadata row order", {
    skip_if_no_mini()
    skip_if_not_installed("harmony")
    g <- .mini()

    # returns a dimObj; `[` gets the embedding
    a <- suppressMessages(runGiottoHarmony(g, vars_use = "leiden_clus",
        return_gobject = FALSE))
    b <- suppressMessages(runGiottoHarmony(.permute_meta(g),
        vars_use = "leiden_clus", return_gobject = FALSE))

    expect_equal(a[], b[])
})

test_that("runDWLSDeconv is invariant to metadata row order", {
    skip_if_no_mini()
    # Deliberately small, for the same reason test-spatial-deconv.R is: DWLS
    # optimizes per cell, and this runs it twice.
    g0 <- .mini()
    g <- GiottoClass::subsetGiotto(g0,
        cell_ids = pDataDT(g0)$cell_ID[seq_len(60L)])
    set.seed(1L)
    genes <- rownames(GiottoClass::getExpression(g, values = "normalized",
        output = "matrix"))
    sm <- suppressMessages(makeSignMatrixDWLS(g,
        sign_gene = sample(genes, 60L),
        cell_type_vector = pDataDT(g)$leiden_clus))

    a <- suppressMessages(runDWLSDeconv(g, sign_matrix = sm,
        cluster_column = "leiden_clus", return_gobject = FALSE))
    b <- suppressMessages(runDWLSDeconv(.permute_meta(g), sign_matrix = sm,
        cluster_column = "leiden_clus", return_gobject = FALSE))

    expect_equal(a, b)
})

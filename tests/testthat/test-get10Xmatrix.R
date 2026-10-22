# get10Xmatrix() on a small MEX folder: duplicated gene symbols get "--n"
# suffixes, and split_by_type separates the control probes.

.write_mex <- function(dir) {
    m <- Matrix::sparseMatrix(i = c(1, 2, 3, 4, 1), j = c(1, 1, 2, 3, 3),
                              x = c(5, 2, 7, 1, 3), dims = c(4, 3))
    Matrix::writeMM(m, file.path(dir, "matrix.mtx"))
    R.utils::gzip(file.path(dir, "matrix.mtx"))
    writeLines(c("g1\tGENE1\tGene Expression",
                 "g2\tGENE2\tGene Expression",
                 "g3\tGENE1\tGene Expression",
                 "g4\tNEG\tNegative Control Probe"),
               gzfile(file.path(dir, "features.tsv.gz")))
    writeLines(paste0("c", 1:3), gzfile(file.path(dir, "barcodes.tsv.gz")))
}

test_that("duplicated gene symbols are made unique", {
    skip_if_not_installed("R.utils")
    d <- withr::local_tempdir()
    .write_mex(d)
    res <- get10Xmatrix(d, gene_column_index = 2, split_by_type = TRUE)
    expect_setequal(names(res), c("Gene Expression", "Negative Control Probe"))
    rna <- res[["Gene Expression"]]
    expect_setequal(rownames(rna), c("GENE1--1", "GENE2", "GENE1--2"))
    expect_equal(as.character(colnames(rna)), paste0("c", 1:3))
    expect_equal(sum(rna), 17)
    expect_equal(as.character(rownames(res[["Negative Control Probe"]])), "NEG")
})

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

# A gene repeated within its type (G1) and a gene that shares its name with
# an antibody (CD4), written as an h5 and an mtx folder
# (helper-visiumhd-fixture.R).
.write_tenx_types <- function(dir) {
    m <- Matrix::sparseMatrix(i = c(1, 2, 3, 4, 5, 1), j = c(1, 1, 2, 3, 3, 2),
                              x = c(5, 2, 7, 1, 3, 4), dims = c(5, 3),
                              dimnames = list(c("G1", "G1", "CD4", "CD4", "G2"),
                                              paste0("c", 1:3)))
    .vhd_fx_write_matrix(m, dir, "raw_feature_bc_matrix", types = c(
        rep("Gene Expression", 3), "Antibody Capture", "Gene Expression"))
}

test_that("mtx and h5 give the same names, counted within each type", {
    skip_if_not_installed("hdf5r")
    d <- withr::local_tempdir()
    .write_tenx_types(d)
    mtx <- get10Xmatrix(file.path(d, "raw_feature_bc_matrix"),
                        gene_column_index = 2)
    h5 <- get10Xmatrix_h5(file.path(d, "raw_feature_bc_matrix.h5"))
    for (r in list(mtx, h5)) {
        expect_identical(rownames(r[["Gene Expression"]]),
                         c("G1--1", "G1--2", "CD4", "G2"))
        expect_identical(rownames(r[["Antibody Capture"]]), "CD4")
    }
    for (ft in names(mtx)) {
        expect_s4_class(mtx[[ft]], "dgCMatrix")
        expect_null(names(rownames(mtx[[ft]])))
        expect_null(names(colnames(mtx[[ft]])))
        expect_equal(mtx[[ft]], h5[[ft]])
    }
})

test_that("unsplit matrices count repeats across the whole list", {
    skip_if_not_installed("hdf5r")
    d <- withr::local_tempdir()
    .write_tenx_types(d)
    mtx <- get10Xmatrix(file.path(d, "raw_feature_bc_matrix"),
                        gene_column_index = 2, split_by_type = FALSE)
    h5 <- get10Xmatrix_h5(file.path(d, "raw_feature_bc_matrix.h5"),
                          split_by_type = FALSE)[[1L]]
    nms <- c("G1--1", "G1--2", "CD4--1", "CD4--2", "G2")
    expect_identical(rownames(mtx), nms)
    expect_identical(rownames(h5), nms)
})

test_that("a one-column features file falls back to column 1", {
    skip_if_not_installed("R.utils")
    d <- withr::local_tempdir()
    .write_mex(d)
    con <- gzfile(file.path(d, "features.tsv.gz"), "w")
    writeLines(c("GENE1", "GENE2", "GENE3", "GENE4"), con)
    close(con)
    expect_warning(
        res <- get10Xmatrix(d, gene_column_index = 2),
        "using column 1"
    )
    expect_identical(rownames(res), paste0("GENE", 1:4))
})

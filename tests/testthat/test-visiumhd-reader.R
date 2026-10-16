# In-memory VisiumHDReader on a synthetic Space Ranger outs directory
# (helper-visiumhd-fixture.R).

skip_if_not_installed("arrow")
skip_if_not_installed("hdf5r")

.fx <- .write_visiumhd_fixture()
.bo <- file.path(.fx$root, "binned_outputs")
.so <- file.path(.fx$root, "segmented_outputs")
.q <- function(x) suppressWarnings(suppressMessages(x))

test_that("loading transcripts keeps the bin's own cell metadata", {
    g <- .q(importVisiumHD(.bo, bin = 8)$create_gobject(
        load_transcripts = TRUE, load_image = FALSE, verbose = FALSE))
    cm <- pDataDT(g, spat_unit = "bin008")
    expect_identical(cm$cell_ID, spatIDs(g, spat_unit = "bin008"))
    expect_false(anyNA(cm$array_row))
})

test_that("cell metadata rows follow the requested barcodes", {
    bc <- rev(.fx$grid[, unique(bc8)])
    cm <- .q(Giotto:::.visiumhd_cellmeta(file.path(.bo, "square_008um"),
        bin = 8, barcodes = bc, verbose = FALSE))
    expect_identical(cm[]$cell_ID, bc)
})

test_that("expression loads from a direct matrix path", {
    m8 <- file.path(.bo, "square_008um")
    for (p in c(file.path(m8, "raw_feature_bc_matrix"),
                file.path(m8, "raw_feature_bc_matrix.h5"))) {
        ex <- .q(Giotto:::.visiumhd_expression(p, bin = 8, verbose = FALSE))
        expect_identical(spatUnit(ex[[1L]]), "bin008")
        expect_identical(featType(ex[[1L]]), "rna")
    }
    ex <- .q(Giotto:::.visiumhd_expression(
        file.path(.so, "raw_feature_cell_matrix.h5"), verbose = FALSE))
    expect_identical(spatUnit(ex[[1L]]), "cell")
})

test_that("an h5 source gives the rna feat_type, as an mtx one does", {
    g <- .q(importVisiumHD(.bo, bin = 8)$create_gobject(
        expression_path = file.path(.bo, "square_008um",
                                    "raw_feature_bc_matrix.h5"),
        load_image = FALSE, verbose = FALSE))
    expect_identical(names(g@expression$bin008), "rna")
})

test_that("get10Xmatrix_h5 gives symbols without splitting by type", {
    m <- get10Xmatrix_h5(file.path(.bo, "square_008um",
                                   "raw_feature_bc_matrix.h5"),
                         gene_ids = "symbols", split_by_type = FALSE,
                         remove_zero_rows = FALSE)
    expect_identical(rownames(m[[1L]]), rownames(.fx$m8))
})

test_that("debug messages name the paths they load", {
    expect_no_error(.q(Giotto:::.visiumhd_tessellate(
        file.path(.bo, "square_008um"), shape_size = 4, bin = 8,
        verbose = "debug")))
})

test_that("load_expression extracts a bin tar that is not extracted yet", {
    bo <- withr::local_tempdir()
    withr::with_dir(.bo, utils::tar(file.path(bo, "008um_outputs.tar"),
                                    "square_008um"))
    rd <- .q(importVisiumHD(bo, bin = 8))
    ex <- .q(rd$load_expression(verbose = FALSE))
    expect_identical(spatIDs(ex[[1L]]), .fx$grid[, unique(bc8)])
})

test_that("micron = TRUE works with the detected (directory) paths", {
    g <- .q(importVisiumHD(.bo, bin = 8)$create_gobject(micron = TRUE,
        load_image = FALSE, verbose = FALSE))
    mpp <- jsonlite::read_json(file.path(.bo, "square_008um", "spatial",
        "scalefactors_json.json"))$microns_per_pixel
    tp <- arrow::read_parquet(file.path(.bo, "square_008um", "spatial",
        "tissue_positions.parquet"))
    sl <- getSpatialLocations(g, spat_unit = "bin008", output = "data.table")
    sl <- sl[match(tp$barcode, cell_ID)]
    expect_equal(sl$sdimx, tp$pxl_col_in_fullres * mpp)
    expect_equal(sl$sdimy, -tp$pxl_row_in_fullres * mpp)
})


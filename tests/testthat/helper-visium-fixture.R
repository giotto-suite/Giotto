# Copy of GiottoDisk's tests/testthat/helper-visium-fixture.R; keep the two
# in step.

# A synthetic Space Ranger `outs` directory for standard Visium, built with
# the Visium HD fixture writers (helper-visiumhd-fixture.R): a 4 x 4 spot grid
# with raw and filtered matrices (mtx directory and .h5), spot positions,
# scalefactors and images. The last array row is out of tissue, so the
# filtered matrix keeps 12 of 16 spots. Gene G5 has no counts.
#
# `positions`: "csv" writes the Space Ranger >= 2 `tissue_positions.csv`
# (with header), "list" the Space Ranger 1 `tissue_positions_list.csv`.
.write_visium_fixture <- function(root = tempfile("visium_outs_"),
                                  positions = "csv") {
    set.seed(7)
    grid <- data.table::CJ(array_row = 0:3, array_col = 0:3)
    grid[, `:=`(
        barcode = sprintf("AAAC%04d-1", .I),
        in_tissue = as.integer(array_row < 3L),
        pxl_row_in_fullres = 100 + 50 * array_row,
        pxl_col_in_fullres = 100 + 50 * array_col
    )]
    genes <- sprintf("G%d", 1:5)
    m <- Matrix::rsparsematrix(length(genes), nrow(grid), 0.4,
        rand.x = function(k) sample(1:4, k, replace = TRUE))
    m[5L, ] <- 0
    m <- Matrix::drop0(m)
    dimnames(m) <- list(genes, grid$barcode)

    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    .vhd_fx_write_matrix(m, root, "raw_feature_bc_matrix")
    .vhd_fx_write_matrix(m[, grid$in_tissue == 1L], root,
                         "filtered_feature_bc_matrix")
    .vhd_fx_write_spatial(root, scalefactors = list(
        spot_diameter_fullres = 40, tissue_hires_scalef = 1,
        fiducial_diameter_fullres = 60, tissue_lowres_scalef = 0.5))

    pos <- grid[, .(barcode, in_tissue, array_row, array_col,
                    pxl_row_in_fullres, pxl_col_in_fullres)]
    sdir <- file.path(root, "spatial")
    if ("csv" %in% positions) {
        data.table::fwrite(pos, file.path(sdir, "tissue_positions.csv"))
    }
    if ("list" %in% positions) {
        data.table::fwrite(pos, file.path(sdir, "tissue_positions_list.csv"),
                           col.names = FALSE)
    }
    list(root = root, grid = grid, counts = m)
}

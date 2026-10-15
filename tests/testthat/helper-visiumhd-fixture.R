# Copy of GiottoDisk's tests/testthat/helper-visiumhd-fixture.R; keep the two
# in step.

# A synthetic Space Ranger v4 Visium HD `outs` directory, small enough to
# write in every test run. It carries the files the readers touch -- binned
# 2 and 8 um outputs (mtx directory and .h5), `tissue_positions.parquet`,
# `scalefactors_json.json`, images, segmented outputs (cell matrix, cell and
# nucleus geojson) and `barcode_mappings.parquet` -- built from one 2 um
# count matrix, so every coarser matrix is its exact aggregate.
#
# Geometry, in fullres pixels at 1 um / px: a 16 x 16 grid of 2 um bins,
# bin (r, c) centred at (x = 2c + 1, y = 2r + 1); the 8 um bins group 4 x 4 of
# them. Cell 1 covers [0, 6]^2, cell 2 covers [6, 12]^2, so 8 um bin (0, 0)
# holds 2 um bins of both cells. The last 8 um row (2 um rows 12-15) is out
# of tissue. Gene G5 has no
# counts, so zero-feature removal has something to drop.

.vhd_fx_bc <- function(um, r, c) sprintf("s_%03dum_%05d_%05d-1", um, r, c)

.vhd_fx_gzip <- function(p) {
    con <- gzfile(paste0(p, ".gz"), "w")
    writeLines(readLines(p), con)
    close(con)
    unlink(p)
}

# 10x feature-barcode h5 (CSC by barcode) + the mtx directory beside it
.vhd_fx_write_matrix <- function(m, dir, stem) {
    genes <- data.frame(id = sprintf("ENSG%05d", seq_len(nrow(m))),
                        name = rownames(m), type = "Gene Expression")
    h5 <- hdf5r::H5File$new(file.path(dir, paste0(stem, ".h5")), mode = "w")
    g <- h5$create_group("matrix")
    m <- methods::as(m, "CsparseMatrix")
    g[["data"]] <- as.integer(m@x)
    g[["indices"]] <- as.integer(m@i)
    g[["indptr"]] <- as.integer(m@p)
    g[["shape"]] <- as.integer(dim(m))
    g[["barcodes"]] <- colnames(m)
    f <- g$create_group("features")
    f[["id"]] <- genes$id
    f[["name"]] <- genes$name
    f[["feature_type"]] <- genes$type
    f[["genome"]] <- rep("GRCh38", nrow(genes))
    f[["_all_tag_keys"]] <- "genome"
    h5$close_all()

    mdir <- file.path(dir, stem)
    dir.create(mdir, showWarnings = FALSE)
    Matrix::writeMM(m, file.path(mdir, "matrix.mtx"))
    .vhd_fx_gzip(file.path(mdir, "matrix.mtx"))
    data.table::fwrite(genes, file.path(mdir, "features.tsv.gz"),
                       sep = "\t", col.names = FALSE)
    data.table::fwrite(data.table::data.table(colnames(m)),
                       file.path(mdir, "barcodes.tsv.gz"), col.names = FALSE)
}

.vhd_fx_write_spatial <- function(dir, positions = NULL, bin_um = NULL) {
    sdir <- file.path(dir, "spatial")
    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    sf <- list(microns_per_pixel = 1, tissue_lowres_scalef = 0.5,
               fiducial_diameter_fullres = 10, tissue_hires_scalef = 1,
               regist_target_img_scalef = 1)
    if (!is.null(bin_um)) {
        sf <- c(list(spot_diameter_fullres = bin_um, bin_size_um = bin_um), sf)
    }
    jsonlite::write_json(sf, file.path(sdir, "scalefactors_json.json"),
                         auto_unbox = TRUE, digits = NA)
    for (nm in c("tissue_hires_image.png", "tissue_lowres_image.png")) {
        r <- terra::rast(nrows = 32, ncols = 32, nlyrs = 3, vals = 128,
                         extent = terra::ext(0, 32, 0, 32))
        suppressWarnings(terra::writeRaster(r, file.path(sdir, nm),
            datatype = "INT1U", overwrite = TRUE))
    }
    if (!is.null(positions)) {
        arrow::write_parquet(positions,
                             file.path(sdir, "tissue_positions.parquet"))
    }
}

.vhd_fx_square <- function(id, lo, hi) {
    list(type = "Feature", properties = list(cell_id = id),
         geometry = list(type = "Polygon", coordinates = list(list(
             c(lo, lo), c(hi, lo), c(hi, hi), c(lo, hi), c(lo, lo)))))
}

.write_visiumhd_fixture <- function(root = tempfile("vhd_outs_")) {
    set.seed(11)
    n <- 16L
    grid <- data.table::CJ(r = 0:(n - 1L), c = 0:(n - 1L))
    grid[, `:=`(
        barcode = .vhd_fx_bc(2L, r, c),
        x = 2 * c + 1, y = 2 * r + 1,
        bc8 = .vhd_fx_bc(8L, r %/% 4L, c %/% 4L)
    )]
    genes <- sprintf("G%d", 1:5)
    m2 <- Matrix::rsparsematrix(length(genes), nrow(grid), 0.3,
        rand.x = function(k) sample(1:4, k, replace = TRUE))
    m2[5L, ] <- 0
    m2 <- Matrix::drop0(m2)
    dimnames(m2) <- list(genes, grid$barcode)

    in_sq <- function(x, y, lo, hi) x > lo & x < hi & y > lo & y < hi
    grid[, cell := data.table::fifelse(in_sq(x, y, 0, 6), 1L,
                   data.table::fifelse(in_sq(x, y, 6, 12), 2L, NA_integer_))]
    grid[, nuc := in_sq(x, y, 1, 5) | in_sq(x, y, 7, 11)]
    cell_ids <- sprintf("cellid_%09d-1", 1:2)

    # coarser matrices as exact aggregates of the 2 um one
    agg <- function(groups, levels) {
        g <- Matrix::sparseMatrix(
            i = seq_along(groups)[!is.na(groups)],
            j = match(groups[!is.na(groups)], levels),
            x = 1, dims = c(length(groups), length(levels)))
        out <- m2 %*% g
        dimnames(out) <- list(genes, levels)
        out
    }
    bc8 <- unique(grid$bc8)
    m8 <- agg(grid$bc8, bc8)
    mc <- agg(cell_ids[grid$cell], cell_ids)

    bo <- file.path(root, "binned_outputs")
    for (um in c(2L, 8L)) {
        d <- file.path(bo, sprintf("square_%03dum", um))
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        m <- if (um == 2L) m2 else m8
        .vhd_fx_write_matrix(m, d, "raw_feature_bc_matrix")
        rows <- if (um == 2L) grid else unique(grid[, .(
            barcode = bc8, r = r %/% 4L, c = c %/% 4L)], by = "barcode")
        pos <- data.frame(
            barcode = rows$barcode,
            in_tissue = as.integer(if (um == 2L) rows$r < n - 4L
                                   else rows$r < (n %/% 4L) - 1L),
            array_row = as.integer(rows$r), array_col = as.integer(rows$c),
            pxl_row_in_fullres = um * rows$r + um / 2,
            pxl_col_in_fullres = um * rows$c + um / 2
        )
        .vhd_fx_write_spatial(d, positions = pos, bin_um = um)
    }

    so <- file.path(root, "segmented_outputs")
    dir.create(so, recursive = TRUE, showWarnings = FALSE)
    .vhd_fx_write_matrix(mc, so, "raw_feature_cell_matrix")
    .vhd_fx_write_spatial(so)
    geo <- function(fs) list(type = "FeatureCollection", features = fs)
    jsonlite::write_json(
        geo(list(.vhd_fx_square(1L, 0, 6), .vhd_fx_square(2L, 6, 12))),
        file.path(so, "cell_segmentations.geojson"), auto_unbox = TRUE)
    jsonlite::write_json(
        geo(list(.vhd_fx_square(1L, 1, 5), .vhd_fx_square(2L, 7, 11))),
        file.path(so, "nucleus_segmentations.geojson"), auto_unbox = TRUE)

    arrow::write_parquet(data.frame(
        square_002um = grid$barcode,
        square_008um = grid$bc8,
        cell_id = cell_ids[grid$cell],
        in_nucleus = grid$nuc & !is.na(grid$cell),
        in_cell = !is.na(grid$cell)
    ), file.path(root, "barcode_mappings.parquet"))

    list(root = root, m2 = m2, m8 = m8, mc = mc, grid = grid)
}

# In-memory VisiumReader on a synthetic Space Ranger outs directory
# (helper-visium-fixture.R).

skip_if_not_installed("hdf5r")

.fx <- .write_visium_fixture(positions = c("csv", "list"))
.q <- function(x) suppressWarnings(suppressMessages(x))
# get10Xmatrix() returns a dgTMatrix with named dimnames, get10Xmatrix_h5() a
# dgCMatrix: compare values
.counts <- function(g) {
    m <- as.matrix(getExpression(g, output = "matrix"))
    dimnames(m) <- lapply(dimnames(m), unname)
    m
}

test_that("an empty reader has its loaders and no paths", {
    x <- importVisium()
    expect_s4_class(x, "VisiumReader")
    expect_true(all(c("load_expression", "load_tissue_position",
        "load_cellmeta", "load_scalefactors", "load_image", "load_polygon",
        "create_gobject") %in% names(x@calls)))
    expect_true(all(vapply(x@paths, is.null, logical(1L))))
})

test_that("paths prefer tissue_positions.csv and the matrix directory", {
    p <- importVisium(.fx$root)@paths
    expect_identical(basename(p$positions_path), "tissue_positions.csv")
    expect_identical(basename(p$expr_path), "raw_feature_bc_matrix")
    expect_identical(basename(p$img_path), "tissue_hires_image.png")

    d <- .write_visium_fixture(positions = "list")$root
    unlink(file.path(d, "raw_feature_bc_matrix"), recursive = TRUE)
    p <- importVisium(d)@paths
    expect_identical(basename(p$positions_path), "tissue_positions_list.csv")
    expect_identical(basename(p$expr_path), "raw_feature_bc_matrix.h5")
})

test_that("both position layouts give the same spots", {
    x <- importVisium(.fx$root)
    sl_csv <- x$load_tissue_position(output = "data.table")
    sl_list <- x$load_tissue_position(
        path = file.path(.fx$root, "spatial", "tissue_positions_list.csv"),
        output = "data.table")
    expect_equal(sl_csv, sl_list)
    expect_identical(sl_csv$cell_ID, .fx$grid$barcode)
    expect_equal(sl_csv$sdimx, .fx$grid$pxl_col_in_fullres)
    expect_equal(sl_csv$sdimy, -.fx$grid$pxl_row_in_fullres)
})

test_that("create_gobject builds spots, metadata, image and polygons", {
    g <- .q(importVisium(.fx$root)$create_gobject(verbose = FALSE))
    expect_identical(spatIDs(g), .fx$grid$barcode)
    expect_setequal(featIDs(g), sprintf("G%d", 1:4))
    cm <- pDataDT(g)
    expect_identical(cm$cell_ID, spatIDs(g))
    expect_true(all(c("in_tissue", "array_row", "array_col") %in% names(cm)))
    expect_identical(list_spatial_info_names(g), "cell")
    expect_length(getPolygonInfo(g, name = "cell", return_giottoPolygon = FALSE),
        16L)
    expect_identical(list_images_names(g), "image")
    # spots two array columns apart are 100 px apart in the fixture
    expect_equal(instructions(g, "micron_scale"), 1)
    pe <- terra::ext(getPolygonInfo(g, name = "cell",
                                    return_giottoPolygon = FALSE)[1])
    expect_equal(pe$xmax - pe$xmin, 55, tolerance = 0.01, ignore_attr = TRUE)
})

test_that("filtered expression keeps the in-tissue spots", {
    x <- importVisium(.fx$root, expression_source = "filtered")
    g <- .q(x$create_gobject(load_image = FALSE, verbose = FALSE))
    expect_identical(spatIDs(g), .fx$grid[in_tissue == 1L, barcode])
    expect_identical(getSpatialLocations(g, output = "data.table")$cell_ID,
        spatIDs(g))
})

test_that("an h5 matrix gives the same object as the mtx directory", {
    x <- importVisium(.fx$root)
    g_mtx <- .q(x$create_gobject(load_image = FALSE, verbose = FALSE))
    g_h5 <- .q(x$create_gobject(
        expression_path = file.path(.fx$root, "raw_feature_bc_matrix.h5"),
        load_image = FALSE, verbose = FALSE))
    expect_equal(.counts(g_h5), .counts(g_mtx))
    expect_equal(getSpatialLocations(g_h5, output = "data.table"),
        getSpatialLocations(g_mtx, output = "data.table"))
})

test_that("setting a reader slot rebuilds the paths", {
    x <- importVisium(.fx$root)
    x$expression_source <- "filtered"
    expect_identical(basename(x@paths$expr_path), "filtered_feature_bc_matrix")
    expect_error(x$expression_source <- "nope", "expression_source")
    expect_error(x$create_gobject <- NULL, "can be set")
})

test_that("createGiottoVisiumObject matches the reader in both modes", {
    g_rd <- .q(importVisium(.fx$root, feature_id_type = "ensembl")$create_gobject(
        verbose = FALSE))
    g_dir <- .q(createGiottoVisiumObject(.fx$root, gene_column_index = 1,
        verbose = FALSE))
    expect_equal(.counts(g_dir), .counts(g_rd))

    sdir <- file.path(.fx$root, "spatial")
    g_h5 <- .q(createGiottoVisiumObject(
        h5_visium_path = file.path(.fx$root, "raw_feature_bc_matrix.h5"),
        h5_gene_ids = "ensembl",
        h5_tissue_positions_path = file.path(sdir, "tissue_positions.csv"),
        h5_json_scalefactors_path = file.path(sdir, "scalefactors_json.json"),
        verbose = FALSE))
    expect_equal(.counts(g_h5), .counts(g_rd))
    expect_length(list_images_names(g_h5), 0L)
    expect_identical(list_spatial_info_names(g_h5), "cell")
})

test_that("backend = returns the GiottoDisk reader", {
    skip_if_not_installed("GiottoDisk")
    skip_if_not("importVisiumDisk" %in% getNamespaceExports("GiottoDisk"))
    x <- importVisium(.fx$root, backend = withr::local_tempdir())
    expect_s4_class(x, "VisiumDiskReader")
})

test_that("the micron scale prefers microns_per_pixel, then the spot grid", {
    pos <- .visium_tissue_positions(
        file.path(.fx$root, "spatial", "tissue_positions.csv"))
    sf <- list(spot_diameter_fullres = 40)
    expect_equal(.visium_micron_scale(sf, pos), 1)
    expect_equal(.visium_micron_scale(c(sf, microns_per_pixel = 0.5), pos), 0.5)
    expect_warning(ms <- .visium_micron_scale(sf, pos[1:2]), "approximate")
    expect_equal(ms, 65 / 40)
    expect_error(.visium_micron_scale(NULL, NULL), "positions or scalefactors")
})

test_that("scalefactors with extra entries and oddly named images load", {
    d <- .write_visium_fixture()$root
    sdir <- file.path(d, "spatial")
    sf <- jsonlite::read_json(file.path(sdir, "scalefactors_json.json"))
    sf$some_new_entry <- 1
    jsonlite::write_json(sf, file.path(sdir, "scalefactors_json.json"),
                         auto_unbox = TRUE, digits = NA)
    expect_no_warning(.visium_read_scalefactors(
        file.path(sdir, "scalefactors_json.json")))
    file.copy(file.path(sdir, "tissue_hires_image.png"),
              file.path(sdir, "aligned_image.png"))
    img <- .q(importVisium(d)$load_image(path = file.path(sdir, "aligned_image.png")))
    expect_s4_class(img[[1L]], "giottoLargeImage")
})

test_that("addVisiumPolygons uses the object's micron scale", {
    g <- .q(importVisium(.fx$root)$create_gobject(load_polygons = FALSE,
        load_image = FALSE, verbose = FALSE))
    expect_length(list_spatial_info_names(g), 0L)
    g <- .q(addVisiumPolygons(g))
    p <- getPolygonInfo(g, name = "cell", return_giottoPolygon = FALSE)
    e <- terra::ext(p[1])
    expect_equal(e$xmax - e$xmin, 55, tolerance = 0.01, ignore_attr = TRUE)
})

# Tests for calculateImageFeatures(): tiles and masks per polygon, for plain
# and affine images, the morphology provider and where the result is stored.
# The UNI2-h provider needs Python and a gated model; its parity with the
# reference embeddings is checked outside the package tests.

# A 200 x 200 px RGB image (1 unit per px) with four discs of radius 8 in
# different colours, and the discs as polygons.
.img_fixture <- function() {
    centres <- data.frame(x = c(40, 150, 60, 160), y = c(50, 40, 150, 160))
    cols <- rbind(c(200, 30, 30), c(30, 200, 30), c(30, 30, 200),
        c(120, 60, 180))
    r <- terra::rast(nrows = 200, ncols = 200, nlyrs = 3,
        extent = terra::ext(0, 200, 0, 200))
    xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
    v <- matrix(255, terra::ncell(r), 3)
    for (i in seq_len(nrow(centres))) {
        inside <- (xy[, 1] - centres$x[i])^2 + (xy[, 2] - centres$y[i])^2 <= 64
        v[inside, ] <- rep(cols[i, ], each = sum(inside))
    }
    terra::values(r) <- v
    disc <- data.table::data.table(circleVertices(8, 60))
    polys <- polyStamp(disc, data.table::data.table(cell_ID = paste0("n", 1:4),
        sdimx = centres$x, sdimy = centres$y), verbose = FALSE)
    g <- gwith_options(list(giotto.verbose = FALSE, giotto.use_conda = FALSE,
            giotto.no_python_warn = TRUE), {
        g <- createGiottoObject(expression = matrix(1, 1, 1,
            dimnames = list("x", "dummy")), verbose = FALSE)
        g <- setGiotto(g, createGiottoPolygon(polys, name = "nucleus",
            verbose = FALSE), verbose = FALSE)
        img <- createGiottoLargeImage(r, name = "he", use_rast_ext = TRUE,
            verbose = FALSE)
        addGiottoLargeImage(g, list(img))
    })
    list(g = g, cols = cols, centres = centres)
}

# custom provider: the centre pixel's colour, the mask at the centre, the
# mask's size and the tile side
.centre_provider <- function(tiles, masks, um_per_px, polys) {
    h <- (dim(tiles)[2] + 1) / 2
    cbind(r = tiles[, h, h, 1], g = tiles[, h, h, 2], b = tiles[, h, h, 3],
        in_mask = masks[, h, h], mask_px = apply(masks, 1, sum),
        side = dim(tiles)[2], um_per_px = um_per_px)
}

.feats <- function(g, values) {
    t(as.matrix(getExpression(g, spat_unit = "nucleus", feat_type = "img",
        values = values, output = "matrix")))
}

test_that("tiles are centred on each polygon and its mask", {
    f <- .img_fixture()
    g <- calculateImageFeatures(f$g, spat_unit = "nucleus", image_name = "he",
        method = .centre_provider, tile_um = 40, block_px = 64L,
        batch_size = 3L)
    m <- .feats(g, "custom")
    expect_identical(rownames(m), paste0("n", 1:4))
    expect_equal(unname(m[, c("r", "g", "b")]), f$cols)
    expect_true(all(m[, "in_mask"] == 1))
    expect_true(all(abs(m[, "mask_px"] - pi * 64) < 0.15 * pi * 64))
    expect_true(all(m[, "side"] == 41))
    expect_true(all(m[, "um_per_px"] == 1))
})

test_that("an affine image is read through its transform", {
    f <- .img_fixture()
    img <- getGiottoImage(f$g, name = "he")
    # rotate and scale image and polygons alike: the tiles must not change
    rot <- matrix(c(cos(0.3), sin(0.3), -sin(0.3), cos(0.3)), 2) * 2
    aimg <- affine(img, rot)
    expect_s4_class(aimg, "giottoAffineImage")
    objName(aimg) <- "he_affine"
    g <- addGiottoLargeImage(f$g, list(aimg))
    poly <- getPolygonInfo(g, name = "nucleus", return_giottoPolygon = TRUE)
    g <- setGiotto(g, affine(poly, rot), verbose = FALSE)
    g <- calculateImageFeatures(g, spat_unit = "nucleus",
        image_name = "he_affine", method = .centre_provider, tile_um = 80)
    m <- .feats(g, "custom")
    expect_equal(unname(m[, c("r", "g", "b")]), f$cols)
    expect_true(all(m[, "in_mask"] == 1))
    expect_equal(unname(m[1, "um_per_px"]), 2, tolerance = 1e-6)
    expect_true(all(m[, "side"] == 41))
})

test_that("the morphology provider measures shape and stain", {
    f <- .img_fixture()
    g <- calculateImageFeatures(f$g, spat_unit = "nucleus", image_name = "he",
        method = "morphology", tile_um = 40)
    m <- .feats(g, "morphology")
    expect_identical(colnames(m), c("area", "perimeter", "circularity",
        "elongation", "mean_r", "mean_g", "mean_b", "sd_r", "sd_g", "sd_b",
        "mean_od"))
    expect_equal(unname(m[, "area"]), rep(pi * 64, 4), tolerance = 0.01)
    expect_true(all(m[, "circularity"] > 0.98))
    expect_true(all(abs(m[, "elongation"] - 1) < 0.01))
    # the rasterised disc is the disc's colour, up to its rim
    expect_equal(unname(m[, c("mean_r", "mean_g", "mean_b")]), f$cols,
        tolerance = 0.05)
    expect_true(all(m[, "mean_od"] > 0))
})

test_that("a polygon outside the image gets a white tile, not an error", {
    f <- .img_fixture()
    far <- polyStamp(data.table::data.table(circleVertices(8, 30)),
        data.table::data.table(cell_ID = "n5", sdimx = 400, sdimy = 400),
        verbose = FALSE)
    poly <- getPolygonInfo(f$g, name = "nucleus", return_giottoPolygon = TRUE)
    all_p <- createGiottoPolygon(rbind(poly[],
        createGiottoPolygon(far, verbose = FALSE)[]), name = "nucleus",
        verbose = FALSE)
    g <- setGiotto(f$g, all_p, verbose = FALSE)
    g <- calculateImageFeatures(g, spat_unit = "nucleus", image_name = "he",
        method = .centre_provider, tile_um = 40)
    m <- .feats(g, "custom")
    expect_equal(unname(m["n5", c("r", "g", "b")]), c(255, 255, 255))
    expect_equal(unname(m[paste0("n", 1:4), c("r", "g", "b")]), f$cols)
})

test_that("features are stored as expression of the spatial unit", {
    f <- .img_fixture()
    g <- calculateImageFeatures(f$g, spat_unit = "nucleus", image_name = "he",
        method = .centre_provider, tile_um = 40, feat_type = "img2",
        name = "centre")
    ex <- getExpression(g, spat_unit = "nucleus", feat_type = "img2",
        values = "centre", output = "exprObj")
    expect_identical(spatUnit(ex), "nucleus")
    expect_identical(colnames(ex[]), paste0("n", 1:4))
    expect_identical(ex@misc$image, "he")
    expect_error(calculateImageFeatures(f$g, spat_unit = "nucleus",
        image_name = "he", method = "nope"), "should be one of")
})

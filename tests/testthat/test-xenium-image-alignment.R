# 10x alignment matrices act on pixel coordinates with y pointing down;
# Giotto images sit at negative y
r <- terra::rast(matrix(1:12, nrow = 3), extent = terra::ext(0, 4, -3, 0))
img <- GiottoClass::createGiottoLargeImage(r, name = "img", verbose = FALSE)
micron <- 0.5

# where the image's pixel-space corner (px, py) should land in Giotto space
place <- function(aff, px, py) {
    p <- aff[1:2, 1:2] %*% c(px, py) + aff[1:2, 3]
    c(p[[1]], -p[[2]]) * micron
}
corners <- function(x) {
    terra::crds(GiottoClass::affine(terra::as.polygons(terra::ext(r)),
        x@affine))[1:4, ]
}

test_that("a sheared 10x alignment places the image corners", {
    aff <- matrix(c(
        1.2, 0.3, 5,
        -0.1, 0.9, 2,
        0, 0, 1
    ), nrow = 3, byrow = TRUE)
    out <- .tenx_img_affine(img, aff, micron)
    got <- corners(out)
    # the raster's corner (x, y) in Giotto space is pixel (x, -y)
    want <- t(apply(terra::crds(terra::as.polygons(terra::ext(r)))[1:4, ], 1,
        function(p) place(aff, p[[1]], -p[[2]])))
    expect_equal(unname(got), unname(want), tolerance = 1e-9)
})

# Stereo-seq bin1 points carry the DNB geometry: a DNB covers [x, x + 1) in
# the GEF, and cells take DNBs by area.

.bin1_gef <- function() {
    list(
        bin1DT = data.table::data.table(x = c(0L, 0L, 1L, 5L),
            y = c(0L, 0L, 2L, 3L), gene_idx = c(1L, 2L, 1L, 2L),
            count = c(1L, 3L, 2L, 1L)),
        bin1GeneDT = data.table::data.table(geneName = c("g1", "g2"))
    )
}

test_that("bin1 points carry the DNB geometry and the area rule", {
    for (neg in c(TRUE, FALSE)) {
        gbp <- .stereoseq_build_binpoints(.bin1_gef(), "geneName",
            negative_y = neg, spat_unit = "bin1", verbose = FALSE)
        expect_identical(gbp@bin_width, 1)
        expect_identical(gbp@bin_offset, c(0.5, if (neg) -0.5 else 0.5))
        expect_identical(gbp@bin_overlap, "area")
    }
})

test_that("cells take the DNBs their outline covers", {
    gbp <- .stereoseq_build_binpoints(.bin1_gef(), "geneName",
        negative_y = FALSE, spat_unit = "bin1", verbose = FALSE)
    # DNB (1, 2) covers [1, 2) x [2, 3): the outline misses its point (the
    # corner) but covers 0.9 of it
    gp <- createGiottoPolygon(data.frame(poly_ID = "c",
        x = c(1.1, 3, 3, 1.1), y = c(1.9, 1.9, 3.5, 3.5)))
    by_area <- binOverlaps(calculateOverlap(gp, gbp), gbp)
    by_point <- binOverlaps(calculateOverlap(gp, gbp, bin_overlap = "inside"),
        gbp)
    expect_identical(by_area$bin_ID, "bin_2")
    expect_identical(nrow(by_point), 0L)
})

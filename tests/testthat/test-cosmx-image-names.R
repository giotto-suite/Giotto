# CosMx exports name their per-FOV images with the FOV number padded to three
# digits (CellComposite_F001.jpg) or, in AtoMx exports, five
# (CellComposite_F00033.jpg). The file name format is read from the files.

.cosmx_img_dir <- function(names) {
    d <- tempfile("cosmx_img_")
    dir.create(d)
    file.create(file.path(d, names))
    d
}

test_that("three-digit FOV image names keep the F%03d format", {
    d <- .cosmx_img_dir(c("CellComposite_F001.jpg", "CellComposite_F002.jpg"))
    expect_identical(.cosmx_file_fov_fmt(d), "F%03d")
})

test_that("five-digit FOV image names are read as F%05d", {
    d <- .cosmx_img_dir(sprintf("CellComposite_F%05d.jpg", 33:34))
    expect_identical(.cosmx_file_fov_fmt(d), "F%05d")
    # a vector of file paths, as the mask loader may receive
    expect_identical(.cosmx_file_fov_fmt(file.path(d, list.files(d))), "F%05d")
    # and the format finds the files the FOV numbers point to
    expect_length(Sys.glob(file.path(d, sprintf("*%s*",
        sprintf(.cosmx_file_fov_fmt(d), 33L)))), 1L)
})

test_that("mixed or unreadable names fall back to F%03d", {
    expect_identical(.cosmx_file_fov_fmt(
        .cosmx_img_dir(c("a_F001.jpg", "a_F00002.jpg"))), "F%03d")
    expect_identical(.cosmx_file_fov_fmt(.cosmx_img_dir(character())), "F%03d")
})

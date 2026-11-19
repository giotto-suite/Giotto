# Stereo-seq file finding: a SAW `outs/` directory and a flat STOmics download
# must resolve to the same files.

.stereoseq_layout <- function(flat) {
    root <- tempfile("stereoseq_")
    gef_dir <- if (flat) root else file.path(root, "feature_expression")
    img_dir <- if (flat) root else file.path(root, "image")
    dir.create(gef_dir, recursive = TRUE)
    dir.create(img_dir, showWarnings = FALSE)
    for (f in c("S.tissue.gef", "S.cellbin.gef", "S.adjusted.cellbin.gef",
        "S.raw.gef")) {
        file.create(file.path(gef_dir, f))
    }
    file.create(file.path(img_dir, "S_DAPI_regist.tif"))
    file.create(file.path(img_dir, "S_HE_mask.tif"))
    root
}

test_that("a flat folder and an outs/ directory find the same files", {
    for (flat in c(FALSE, TRUE)) {
        root <- .stereoseq_layout(flat)
        gef_dir <- if (flat) root else file.path(root, "feature_expression")
        img_dir <- if (flat) root else file.path(root, "image")

        bin <- .stereoseq_detect_paths(root, type = "bin", gef_type = "tissue")
        expect_identical(bin$gef_path, file.path(gef_dir, "S.tissue.gef"))
        expect_identical(bin$image_dir, img_dir)
        expect_identical(bin$mask_path, file.path(img_dir, "S_HE_mask.tif"))

        cell <- .stereoseq_detect_paths(root, type = "cell",
            gef_type = "adjusted_cellbin")
        expect_identical(cell$gef_path,
            file.path(gef_dir, "S.adjusted.cellbin.gef"))
        expect_identical(cell$bin1_gef_path, file.path(gef_dir, "S.tissue.gef"))
        expect_identical(
            .stereoseq_find_gef(root, "cellbin"),
            file.path(gef_dir, "S.cellbin.gef")
        )
        unlink(root, recursive = TRUE)
    }
})

test_that("a missing directory or GEF type is reported", {
    expect_error(.stereoseq_find_gef(tempfile("nope_"), "tissue"),
        "directory not found")
    root <- tempfile("stereoseq_empty_")
    dir.create(root)
    expect_error(.stereoseq_find_gef(root, "tissue"),
        "No gef_type=\"tissue\" file found")
    unlink(root, recursive = TRUE)
})

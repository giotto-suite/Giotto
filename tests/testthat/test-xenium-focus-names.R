# Xenium Onboard Analysis 2.x names its focus images morphology_focus_000k;
# the channel of each comes from the OME header.

test_that("XOA 2.x focus images are named by their OME channel", {
    skip_on_cran()
    skip_if_offline()
    d <- GiottoData::getSpatialDataset("xenium_mini_lung_ffpe",
        directory = tempdir())
    f <- list.files(file.path(d, "morphology_focus"), full.names = TRUE)
    nm <- .xenium_focus_channel_names(f)
    expect_identical(unname(nm[basename(f)]),
        c("dapi", "18s", "atp1a1_cd45_e-cadherin", "alphasma_vimentin"))
})

test_that("without an OME header there are no channel names", {
    f <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::rast(matrix(1:4, 2)), f)
    expect_null(.xenium_focus_channel_names(f))
    g <- sub("\\.tif$", ".ome.tif", f)
    file.copy(f, g)
    expect_null(.xenium_focus_channel_names(g))
})

test_that("a Xenium folder without morphology_focus loads without images", {
    skip_on_cran()
    skip_if_offline()
    d <- GiottoData::getSpatialDataset("xenium_mini_lung_ffpe",
        directory = tempdir())
    nofocus <- file.path(tempdir(), "xenium_nofocus")
    unlink(nofocus, recursive = TRUE)
    dir.create(nofocus)
    keep <- setdiff(list.files(d), "morphology_focus")
    file.copy(file.path(d, keep), nofocus, recursive = TRUE)

    expect_warning(
        g <- createGiottoXeniumObject(nofocus, load_transcripts = FALSE,
            verbose = FALSE),
        "No morphology_focus image found"
    )
    expect_s4_class(g, "giotto")
    expect_length(GiottoClass::list_images_names(g, img_type = "largeImage"), 0L)
})

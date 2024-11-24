# download files for Adult Mouse Olfactory Bulb (Space Ranger v1.3.0) ####
# https://www.10xgenomics.com/datasets/adult-mouse-olfactory-bulb-1-standard
# chosen for the small size of the dataset

urls <- c(
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_raw_feature_bc_matrix.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_raw_feature_bc_matrix.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz"
)

datadir <- paste0(getwd(), "/testdata/vis_1_3_0")
if (!dir.exists(datadir)) dir.create(datadir)

# download files
lapply(
    urls,
    function(url) {
        myfilename <- basename(url)
        mydestfile <- file.path(datadir, myfilename)
        download.file(url = url, destfile = mydestfile, quiet = TRUE)
    }
)

manifest <- file.path(datadir, basename(urls))
names(manifest) <- gsub("Visium_Mouse_Olfactory_Bulb_", "", basename(manifest))

lapply(
    manifest[c(
        "filtered_feature_bc_matrix.tar.gz",
        "raw_feature_bc_matrix.tar.gz",
        "spatial.tar.gz"
    )],
    untar,
    exdir = datadir
)

options("giotto.use_conda" = FALSE)
# warnings related to this option must be suppressed



# create: directory ####

ext_to_rounded_num <- function(e) {
    num <- round(ext(e)[])
    names(num) <- NULL
    return(num)
}

test_that("visium create dir raw is working", {
    # unfiltered (all spots)
    g_nofil <- suppressWarnings(createGiottoVisiumObject(
        visium_dir = datadir,
        expr_data = "raw",
        verbose = FALSE
    ))
    expect_true(validObject(g_nofil))
    # meta
    expect_true(nrow(pDataDT(g_nofil)) == 4992)
    expect_true(nrow(fDataDT(g_nofil)) == 19500)
    expect_true("ENSMUSG00000051951" %in% fDataDT(g_nofil)$feat_ID)
    # expression
    expect_true(identical(dim(getExpression(g_nofil)), c(19500L, 4992L)))
    # spatlocs
    expect_equal(nrow(getSpatialLocations(g_nofil)), 4992L)
    # image
    i_nofil <- getGiottoImage(g_nofil, image_type = "largeImage")
    expect_s4_class(i_nofil, "giottoLargeImage")
    expect_identical(dim(i_nofil), c(2000, 2000, 3))
    expect_equal(ext_to_rounded_num(i_nofil), c(0, 1e4, -1e4, 0))


    # polys
    p_nofil <- getPolygonInfo(g_nofil, return_giottoPolygon = TRUE)
    expect_s4_class(p_nofil, "giottoPolygon")
    expect_true(nrow(p_nofil) == 4992)
    expect_equal(ext_to_rounded_num(p_nofil), c(1057, 8726, -8705, -1431))
})

test_that("visium create dir filtered is working", {
    # filtered (in_tissue spots only)
    g_fil <- suppressWarnings(createGiottoVisiumObject(
        visium_dir = datadir,
        expr_data = "filter",
        verbose = FALSE
    ))
    expect_true(validObject(g_fil))
    # meta
    expect_true(nrow(pDataDT(g_fil)) == 1185)
    expect_true(nrow(fDataDT(g_fil)) == 19332)
    expect_true("ENSMUSG00000051951" %in% fDataDT(g_fil)$feat_ID)
    # expression
    expect_true(identical(dim(getExpression(g_fil)), c(19332L, 1185L)))
    # spatlocs
    expect_equal(nrow(getSpatialLocations(g_fil)), 1185L)
    # image
    i_fil <- getGiottoImage(g_fil, image_type = "largeImage")
    expect_s4_class(i_fil, "giottoLargeImage")
    expect_identical(dim(i_fil), c(2000, 2000, 3))
    expect_equal(ext_to_rounded_num(i_fil), c(0, 1e4, -1e4, 0))

    # polys
    p_fil <- getPolygonInfo(g_fil, return_giottoPolygon = TRUE)
    expect_s4_class(p_fil, "giottoPolygon")
    expect_true(nrow(p_fil) == 1185)
    expect_equal(ext_to_rounded_num(p_fil), c(3917, 7244, -7625, -2737))
})





# create: H5 ####
test_that("visium create H5 raw is working", {
    g_nofil <- suppressWarnings(createGiottoVisiumObject(
        h5_visium_path = manifest["raw_feature_bc_matrix.h5"],
        h5_gene_ids = "symbols",
        h5_json_scalefactors_path = file.path(
            datadir, "spatial", "scalefactors_json.json"
        ),
        h5_image_png_path = file.path(
            datadir, "spatial", "tissue_lowres_image.png"
        ),
        h5_tissue_positions_path = file.path(
            datadir, "spatial", "tissue_positions_list.csv"
        )
    ))

    expect_true(validObject(g_nofil))
    # meta
    expect_true(nrow(pDataDT(g_nofil)) == 4992)
    expect_true(nrow(fDataDT(g_nofil)) == 19500)
    expect_true("Xkr4" %in% fDataDT(g_nofil)$feat_ID)
    # expression
    expect_true(identical(dim(getExpression(g_nofil)), c(19500L, 4992L)))
    # spatlocs
    expect_equal(nrow(getSpatialLocations(g_nofil)), 4992L)
    # image
    i_nofil <- getGiottoImage(g_nofil, image_type = "largeImage")
    expect_s4_class(i_nofil, "giottoLargeImage")
    expect_identical(dim(i_nofil), c(600, 600, 3))
    expect_equal(ext_to_rounded_num(i_nofil), c(0, 1e4, -1e4, 0))


    # polys
    p_nofil <- getPolygonInfo(g_nofil, return_giottoPolygon = TRUE)
    expect_s4_class(p_nofil, "giottoPolygon")
    expect_true(nrow(p_nofil) == 4992)
    expect_equal(ext_to_rounded_num(p_nofil), c(1057, 8726, -8705, -1431))
})


test_that("visium create H5 filtered is working", {
    g_fil <- suppressWarnings(createGiottoVisiumObject(
        h5_visium_path = manifest["filtered_feature_bc_matrix.h5"],
        h5_gene_ids = "ensembl",
        h5_json_scalefactors_path = file.path(
            datadir, "spatial", "scalefactors_json.json"
        ),
        h5_image_png_path = file.path(
            datadir, "spatial", "tissue_lowres_image.png"
        ),
        h5_tissue_positions_path = file.path(
            datadir, "spatial", "tissue_positions_list.csv"
        )
    ))

    expect_true(validObject(g_fil))
    # meta
    expect_true(nrow(pDataDT(g_fil)) == 1185)
    expect_true(nrow(fDataDT(g_fil)) == 19332)
    expect_true("ENSMUSG00000051951" %in% fDataDT(g_fil)$feat_ID)
    # expression
    expect_true(identical(dim(getExpression(g_fil)), c(19332L, 1185L)))
    # spatlocs
    expect_equal(nrow(getSpatialLocations(g_fil)), 1185L)
    # image
    i_fil <- getGiottoImage(g_fil, image_type = "largeImage")
    expect_s4_class(i_fil, "giottoLargeImage")
    expect_identical(dim(i_fil), c(600, 600, 3))
    expect_equal(ext_to_rounded_num(i_fil), c(0, 1e4, -1e4, 0))

    # polys
    p_fil <- getPolygonInfo(g_fil, return_giottoPolygon = TRUE)
    expect_s4_class(p_fil, "giottoPolygon")
    expect_true(nrow(p_fil) == 1185)
    expect_equal(ext_to_rounded_num(p_fil), c(3917, 7244, -7625, -2737))
})

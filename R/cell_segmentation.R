#' @title doCellSegmentation
#' @name doCellSegmentation
#' @description segment cells in Dapi image
#' @param raster_img raster image; Dapi image of cells
#' @param folder_path character; where to save the file
#' @param reduce_resolution numeric; the original Dapi image from Vizgen works
#' better in the Mesmer algorithm if its resolution is reduced 4 times.
#' @param overlapping_pixels numeric; the number of pixels to overlap when
#' calculating the rolling window
#' @param python_path specify specific path to python if required
#' @returns segmentation results
#' @details
#' This function is a wrapper for the Mesmer algorithm implemented in python.
#' It segments the cells in tissue by applying a rolling window on the complete
#' image. It saves all the files in the specified folder with the coordinates
#' of the tile: sx (start x), ex (end x), sy, and ey.
#'
#' @export
doCellSegmentation <- function(
        raster_img,
        folder_path,
        reduce_resolution = 4,
        overlapping_pixels = 50,
        python_path = NULL) {
    package_check("deepcell", repository = "pip")
    package_check("PIL", repository = "pip")

    # prepare python path and segmentation script
    reticulate::use_python(required = TRUE, python = python_path)
    python_segmentation_function <- system.file("python",
        "python_segmentation.py",
        package = "Giotto"
    )
    reticulate::source_python(file = python_segmentation_function)

    # create mesmer app
    mesmer_app <- python_create_mesmer_app()

    # get params for rolling window
    img_dim <- dim(raster_img)
    xdim <- img_dim[2]
    ydim <- img_dim[1]

    mesmer_tile_dim <- 512
    tile_dim <- mesmer_tile_dim * reduce_resolution
    margin <- overlapping_pixels
    window_size <- tile_dim - margin

    nxwindow <- xdim %/% window_size
    nywindow <- ydim %/% window_size

    # sliding window
    start_x <- 0
    end_x <- start_x + tile_dim
    for (i in seq_len(nxwindow)) {
        start_y <- 0
        end_y <- start_y + tile_dim
        for (j in seq_len(nywindow)) {
            ext_crop <- terra::ext(c(start_x, end_x, start_y, end_y))
            img_crop <- terra::crop(raster_img, ext_crop, snap = "in")
            img_crop_rescaled <- terra::aggregate(img_crop, reduce_resolution)
            img_coordinates <- as.character(paste("sx", start_x,
                "ex", end_x,
                "sy", start_y,
                "ey", end_y,
                sep = "_"
            ))
            file_name <- file.path(folder_path, paste0(img_coordinates, ".png"))
            segmentation_result <- python_segment_image(
                mesmer_app,
                Matrix::as.matrix(img_crop_rescaled,
                    wide = TRUE
                ),
                file_name
            )

            start_y <- end_y - margin
            end_y <- start_y + tile_dim
        }
        start_x <- end_x - margin
        end_x <- start_x + tile_dim
    }


    print(segmentation_result)
}





#'
#' @title perform cellpose segmentation
#' @description
#'
#' perform the Giotto Wrapper of cellpose segmentation. This is for a model
#' inference to generate segmentation mask file from input image.
#' main parameters needed
#' @name doCellposeSegmentation
#' @param image_dir character, required. Provide a path to a gray scale or a
#' three channel image.
#' @param python_path python environment with cellpose installed.
#' default = "giotto_cellpose".
#' @param mask_output required. Provide a path to the output mask file.
#' @param channel_1 channel number for cytoplasm, default to 0(gray scale)
#' @param channel_2 channel number for Nuclei, default to 0(gray scale)
#' @param model_name Name of the model to run inference. Default to 'cyto3',
#' if you want to run cutomized trained model, place your model file in
#' ~/.cellpose/models and specify your model name.
#' @param batch_size Cellpose Parameter, Number of 224x224 patches to run
#' simultaneously on the GPU. Can make smaller or bigger depending on GPU
#' memory usage. Defaults to 8.
#' @param resample Cellpose Parameter
#' @param channel_axis Cellpose Parameter
#' @param z_axis Cellpose Parameter
#' @param normalize Cellpose Parameter
#' @param invert Cellpose Parameter
#' @param rescale Cellpose Parameter
#' @param diameter Cellpose Parameter
#' @param flow_threshold Cellpose Parameter
#' @param cellprob_threshold Cellpose Parameter
#' @param do_3D Cellpose Parameter
#' @param anisotropy Cellpose Parameter
#' @param stitch_threshold Cellpose Parameter
#' @param min_size Cellpose Parameter
#' @param niter Cellpose Parameter
#' @param augment Cellpose Parameter
#' @param tile Cellpose Parameter
#' @param tile_overlap Cellpose Parameter
#' @param bsize Cellpose Parameter
#' @param interp Cellpose Parameter
#' @param compute_masks Cellpose Parameter
#' @param progress Cellpose Parameter
#' @returns No return variable, as this will write directly to output path
#' provided.
#' @examples
#' # example code
#' doCellposeSegmentation(
#'     image_dir = input_image,
#'     mask_output = output, channel_1 = 2,
#'     channel_2 = 1, model_name = "cyto3", batch_size = 4
#' )
#' @export
doCellposeSegmentation <- function(
        python_env = "giotto_cellpose",
        image_dir,
        mask_output,
        channel_1 = 0,
        channel_2 = 0,
        model_name = "cyto3",
        batch_size = 8,
        resample = TRUE,
        channel_axis = NULL,
        z_axis = NULL,
        normalize = TRUE,
        invert = FALSE,
        rescale = NULL,
        diameter = NULL,
        flow_threshold = 0.4,
        cellprob_threshold = 0.0,
        do_3D = FALSE,
        anisotropy = NULL,
        stitch_threshold = 0.0,
        min_size = 15,
        niter = NULL,
        augment = FALSE,
        tile = TRUE,
        tile_overlap = 0.1,
        bsize = 224,
        interp = TRUE,
        compute_masks = TRUE,
        progress = NULL,
        verbose = TRUE, ...) {
    # Check Input arguments
    model_name <- match.arg(
        model_name, unique(c("cyto3", "cyto2", "cyto", "nuclei", model_name))
    )
    ## Load required python libraries
    GiottoClass::set_giotto_python_path(python_env)
    GiottoUtils::package_check("cellpose", repository = "pip")

    cellpose <- reticulate::import("cellpose")
    np <- reticulate::import("numpy")
    cv2 <- reticulate::import("cv2")
    torch <- reticulate::import("torch")
    message("successfully loaded giotto environment with cellpose.")

    if (!(torch$cuda$is_available())) {
        warning("GPU is not available for this session, inference may be slow.")
    }

    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE, "Loading Image from ",
        image_dir
    )

    img <- cellpose$io$imread(image_dir)
    GiottoUtils::vmsg(.v = verbose, .is_debug = FALSE, "Loading Model...")

    model_to_seg <- cellpose$models$Cellpose(
        model_type = model_name,
        gpu = torch$cuda$is_available()
    )
    channel_to_seg <- as.integer(c(channel_1, channel_2))

    GiottoUtils::vmsg(.v = verbose, .is_debug = FALSE, "Segmenting Image...")
    segmentation <- model_to_seg$eval

    result <- segmentation(img,
        diameter = diameter,
        channels = channel_to_seg,
        batch_size = batch_size,
        resample = resample,
        channel_axis = channel_axis,
        z_axis = z_axis,
        normalize = normalize,
        invert = invert,
        rescale = rescale,
        flow_threshold = flow_threshold,
        cellprob_threshold = cellprob_threshold,
        do_3D = do_3D,
        anisotropy = anisotropy,
        stitch_threshold = stitch_threshold,
        min_size = min_size,
        niter = niter,
        augment = augment,
        tile = tile,
        tile_overlap = tile_overlap,
        bsize = bsize,
        interp = interp,
        compute_masks = compute_masks,
        progress = progress
    )
    masks <- result[[1]]
    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE,
        "Segmentation finished... Saving mask file..."
    )
    GiottoUtils::package_check("terra")
    rast <- terra::rast(masks)
    terra::writeRaster(rast, mask_output, overwrite = TRUE)
}

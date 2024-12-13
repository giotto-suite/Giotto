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
doCellSegmentation <- function(raster_img,
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
#' Perform the Giotto Wrapper of cellpose segmentation. This is for a model
#' inference to generate segmentation mask file from input image.
#' main parameters needed
#' @name doCellposeSegmentation
#' @param input character, required. Provide a path to a gray scale or a
#' three channel image.
#' @param python_path python environment with cellpose installed.
#' default = "giotto_segmentation".
#' @param mask_output required. Provide a path to the output mask file.
#' @param channel_1 channel number for cytoplasm, default to 0(gray scale)
#' @param channel_2 channel number for Nuclei, default to 0(gray scale)
#' @param model_name Name of the model to run inference. Default to 'cyto3',
#' if you want to run cutomized trained model, place your model file in
#' ~/.cellpose/models and specify your model name.
#' @param batch_size Cellpose Parameter, Number of 224x224 patches to run
#' simultaneously on the GPU. Can make smaller or bigger depending on GPU
#' memory usage. Defaults to 8.
#' @param resample (logical, optional) – run dynamics at original image size
#' (will be slower but create more accurate boundaries). Defaults to True.
#' @param channel_axis (int, optional) – channel axis in element of list x, or
#' of np.ndarray x. if NULL, channels dimension is attempted to be
#' automatically determined. Defaults to NULL.
#' @param z_axis (int, optional) – z axis in element of list x, or of
#' np.ndarray x. if NULL, z dimension is attempted to be automatically
#' determined. Defaults to NULL.
#' @param normalize Logical or list. Controls image normalization:
#'        If TRUE, normalizes to 1st-99th percentile
#'        Can be a list with parameters:
#'        * lowhigh: Numeric vector c(low, high) for manual normalization values
#'        * sharpen: Numeric, image sharpening factor (1/4-1/8 of cell diameter in pixels)
#'        * normalize: Logical, whether to run normalization
#'        * percentile: Numeric vector c(low_perc, high_perc)
#'        * tile_norm: Integer, window size in pixels for tile normalization
#'        * norm3D: Logical, normalize across full z-stack
#'        Default: TRUE
#' @param invert Logical. Inverts pixel intensity before processing.
#' Default: FALSE
#' @param rescale Numeric. Resize factor for images.
#'        Default: NULL (sets to 1.0)
#' @param diameter Numeric. Cell diameter for each image.
#'        Default: NULL (uses diam_mean or diam_train if available)
#' @param flow_threshold Numeric. Flow error threshold for cell retention (2D only).
#'        Default: 0.4
#' @param cellprob_threshold Numeric. Threshold for mask pixel retention.
#'        Default: 0.0
#' @param do_3D Logical. Enable 3D segmentation for 3D/4D inputs.
#'        Default: FALSE
#' @param anisotropy Numeric. Z-axis rescaling factor for 3D segmentation.
#'        Default: NULL
#' @param stitch_threshold Numeric. Threshold for 3D mask stitching (2D mode only).
#'        Default: 0.0
#' @param min_size Integer. Minimum ROI size in pixels.
#'        Default: 15
#' @param max_size_fraction Numeric. Maximum mask size as fraction of image.
#'        Default: 0.4
#' @param niter Integer. Number of dynamics computation iterations.
#'        Default: NULL (set based on diameter)
#' @param augment Logical. Enable tile augmentation with overlapping.
#'        Default: FALSE
#' @param tile_overlap Numeric. Tile overlap fraction for flow computation.
#'        Default: 0.1
#' @param bsize Integer. Block size for tiles (recommended: 224).
#'        Default: 224
#' @param dP_smooth Integer. Gaussian smoothing standard deviation for 3D flows.
#' @param interp Logical. Enable interpolation for 2D dynamics.
#'        Default: TRUE
#' @param compute_masks Logical. Compute dynamics and return masks.
#'        Default: TRUE
#' @param progress progress bar. Defaults to NULL
#' @returns No return variable, as this will write directly to output path
#' provided.
#' @examples
#' doCellposeSegmentation(
#'     input = input_image,
#'     mask_output = output, channel_1 = 2,
#'     channel_2 = 1, model_name = "cyto3", batch_size = 4
#' )
#' @export
doCellposeSegmentation <- function(input,
    mask_output,
    python_env = "giotto_segmentation",
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
    max_size_fraction = 0.4,
    niter = NULL,
    augment = FALSE,
    tile_overlap = 0.1,
    bsize = 224,
    dP_smooth,
    interp = TRUE,
    compute_masks = TRUE,
    progress = NULL,
    verbose = NULL,
    ...) {
    ## Load required python libraries
    GiottoClass::set_giotto_python_path(python_env)
    GiottoUtils::package_check("cellpose>=3.1.0", repository = "pip")

    # Check Input arguments
    model_name <- match.arg(
        model_name, unique(c("cyto3", "cyto2", "cyto", "nuclei", model_name))
    )
    if (!is.null(channel_axis)) channel_axis <- as.integer(channel_axis)

    cellpose <- reticulate::import("cellpose")
    np <- reticulate::import("numpy")
    cv2 <- reticulate::import("cv2")
    torch <- reticulate::import("torch")
    message("successfully loaded giotto environment with cellpose.")

    if (!(torch$cuda$is_available())) {
        warning(
            "GPU is not available for this session, inference may be slow."
        )
    }

    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE, "Loading Image from ", input
    )

    img <- cellpose$io$imread(input)
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
        batch_size = as.integer(batch_size),
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
        min_size = as.integer(min_size),
        max_size_fraction = max_size_fraction,
        niter = niter,
        augment = augment,
        tile_overlap = tile_overlap,
        bsize = as.integer(bsize),
        dP_smooth = as.integer(dP_smooth),
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




#'
#' @title perform Mesmer(Deepcell) segmentation
#' @description
#' Perform the Giotto Wrapper of mesmer segmentation. This is for a model
#' inference to generate segmentation mask file from input image.
#' main parameters needed
#' @name doMesmerSegmentation
#' @param input character, required. Provide a path to a IF image.
#' @param python_env python environment with deepcell installed.
#' default = "giotto_segmentation". See deepcell official website for more details.
#' @param mask_output required. Provide a path to the output mask file.
#' @param nucleus_channel channel number for Nuclei, default to 1
#' @param membrane_channel channel number for cell boundary, default to 2
#' @param micron_scale numeric. Multiplicative scalefactor to convert pixel
#' dimensions to physical microns.
#' @returns No return variable, as this will write directly to output path
#' provided.
#' @examples
#' doMesmerSegmentation(
#'     input = input_image,
#'     mask_output = output,
#'     nucleus_channel = 1,
#'     membrane_channel = 2,
#'     micron_scale = 0.5
#' )
#' @export
doMesmerSegmentation <- function(
        input,
        mask_output,
        python_env = "giotto_segmentation",
        nucleus_channel = 1,
        membrane_channel = 2,
        micron_scale = 0.25,
        verbose = NULL,
        ...) {
    ## Load required python libraries
    GiottoClass::set_giotto_python_path(python_env)
    GiottoUtils::package_check("deepcell", repository = "pip")
    tiff <- reticulate::import("tifffile")
    np <- reticulate::import("numpy")
    deepcell <- reticulate::import("deepcell.applications")
    message("successfully loaded giotto environment with deepcell.")

    # Initialize the Mesmer application from DeepCell
    mesmer <- deepcell$Mesmer()

    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE, "Loading Image... ",
    )
    GiottoUtils::package_check("terra")
    rast <- terra::rast(input)
    # Convert the R matrix to a NumPy array explicitly
    nucleus_channel_np <- np$array(drop(terra::as.array(rast[[as.numeric(nucleus_channel)]])))
    membrane_channel_np <- np$array(drop(terra::as.array(rast[[as.numeric(membrane_channel)]])))
    stacked_array <- np$stack(list(nucleus_channel_np, membrane_channel_np), axis = as.integer(-1))
    # Add a new axis to the stacked array to fit Mesmer input
    stacked_array <- np$expand_dims(stacked_array, axis = as.integer(0))

    GiottoUtils::vmsg(.v = verbose, .is_debug = FALSE, "Segmenting Image...")

    segmentation_predictions <- mesmer$predict(
        stacked_array,
        image_mpp = micron_scale
    )
    mask <- segmentation_predictions[1, , , 1]
    mask_r <- reticulate::py_to_r(mask)

    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE,
        "Segmentation finished... Saving mask file..."
    )

    rast <- terra::rast(mask_r)
    terra::writeRaster(rast, mask_output, overwrite = TRUE)
}



#'
#' @title perform Stardist segmentation
#' @description
#' Perform the Giotto wrapper of Stardist 2D segmentation. This is for a model
#' inference to generate segmentation mask file from input image.
#' main parameters needed
#' @name doStardistSegmentation
#' @param input character, required. Provide a path to an image.
#' @param python_env python environment with Stardist installed.
#' default = "giotto_segmentation". See Stardist official website for more details.
#' @param mask_output required. Provide a path to the output mask file.
#' @param model_name Name of the model to run inference. Default to '2D_versatile_fluo'.
#' If using HE model, input image must be RGB, else the nuclei_channel must be given
#' @param nuclei_channel Required using IF based nuclei segmentation, channel number of the nuclei staining.
#' @param prob_thresh prob_thresh for model (if not given use model default)
#' @param nms_thresh nms_thresh for model (if not given use model default)
#' @returns No return variable, as this will write directly to output path provided.
#' @examples
#' # example code
#' doStardistSegmentation(
#'     input = input_image,
#'     mask_output = output,
#'     model_name = "2D_versatile_fluo",
#'     nuclei_channel = 3
#' )
#'
#' @export
doStardistSegmentation <- function(
        input,
        mask_output,
        python_env = "giotto_segmentation",
        model_name = "2D_versatile_fluo",
        nuclei_channel = NULL,
        prob_thresh = NULL,
        nms_thresh = NULL,
        verbose = NULL,
        ...) {
    # Import the necessary Python modules
    ## Load required python libraries
    GiottoClass::set_giotto_python_path(python_env)
    GiottoUtils::package_check("stardist", repository = "pip")
    stardist <- reticulate::import("stardist.models")
    csbdeep <- reticulate::import("csbdeep.utils")
    np <- reticulate::import("numpy")

    # Load the StarDist2D model
    model_name <- match.arg(
        model_name, unique(c("2D_versatile_fluo", "2D_versatile_he", "2D_paper_dsb2018", "2D_demo", model_name))
    )
    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE, "Loading model ",
        model_name
    )

    model <- stardist$StarDist2D$from_pretrained(model_name)

    # Load the image
    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE, "Loading Image from ",
        input
    )
    GiottoUtils::package_check("terra")
    rast <- terra::rast(input)
    if (model_name != "2D_versatile_he" && is.null(nuclei_channel)) {
        stop("using IF based nuclei segmentation, please specify nuclei channel")
    } else if (model_name == "2D_versatile_he") {
        img <- np$array(terra::as.array(rast))
    } else {
        img <- np$array(drop(terra::as.array(rast[[as.numeric(nuclei_channel)]])))
    }

    # Normalize the image
    normalized_img <- csbdeep$normalize(img)

    # Perform prediction with StarDist2D model
    results <- model$predict_instances(normalized_img,
        prob_thresh = prob_thresh,
        nms_thresh = nms_thresh
    )

    # Extract the labels (first output from predict_instances)
    mask <- results[[1]]
    GiottoUtils::vmsg(
        .v = verbose, .is_debug = FALSE,
        "Segmentation finished... Saving mask file..."
    )
    rast_m <- terra::rast(mask)
    terra::writeRaster(rast_m, mask_output, overwrite = TRUE)
}

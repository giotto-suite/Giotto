# Spatial Method-Specific Convenience Functions for Giotto Object Creation #



# Common Utility Functions ####

#' @title Read a structured folder of exported data
#' @name read_data_folder
#' @description Framework function for reading the exported folder of a spatial
#' method and detecting the presence of needed files. NULL values denote missing
#' items.\cr
#' `.read_data_folder()` should not be called directly. Instead, specific
#' reader functions should be built using it as a base.
#' @param spat_method spatial method for which the data is being read
#' @param data_dir exported data directory to read from
#' @param dir_items named list of directory items to expect and keywords to
#' match
#' @param data_to_use character. Which type(s) of expression data to build the
#' gobject with. Values should match with a *workflow* item in require_data_DT
#' (see details)
#' @param require_data_DT data.table detailing if expected data items are
#' required or optional for each \code{data_to_use} *workflow*
#' @param cores cores to use
#' @param verbose be verbose
#' @param toplevel stackframes back where the user-facing function was called.
#' default is one stackframe above `.read_data_folder`.
#' @returns data.table
#' @details
#' **Steps performed:**
#' \itemize{
#'   \item{1. detection of items within \code{data_dir} by looking for keywords
#'   assigned through \code{dir_items}}
#'   \item{2. check of detected items to see if everything needed has been
#'   found. Dictionary of necessary vs optional items for each
#'   \code{data_to_use} *workflow* is provided through \code{require_data_DT}}
#'   \item{3. if multiple filepaths are found to be matching then select the
#'   first one. This function is only intended to find the first level
#'   subdirectories and files.}
#' }
#'
#' **Example reader implementation:**
#' \preformatted{
#'   foo <- function(x_dir,
#'                   data_to_use,
#'                   cores = NA,
#'                   verbose = NULL) {
#'     dir_items <- list(
#'       data1 = "regex_pattern1",
#'       data2 = "regex_pattern2",
#'       data3 = "regex_pattern3"
#'     )
#'
#'     # DT of info to check directory for. Has 3 cols
#'     require_data_DT <- data.table::data.table(
#'       workflow = "a", # data_to_use is matched against this
#'       item = c(
#'         "data1",
#'         "data2",
#'         "data3"
#'       ),
#'       needed = c(
#'         FALSE, # data1 optional for this workflow (if missing: warn)
#'         TRUE,  # data2 vital for this workflow (if missing: error)
#'         TRUE   # data3 vital for this workflow (if missing: error)
#'       )
#'     )
#'
#'     .read_data_folder(
#'       spat_method = "x_method",
#'       data_dir = x_dir,
#'       dir_items = dir_items,
#'       data_to_use = data_to_use,
#'       require_data_DT = require_data_DT,
#'       cores = cores,
#'       verbose = verbose
#'     )
#'   }
#' }
#'
#' @md
NULL

#' @describeIn read_data_folder Should not be used directly
#' @keywords internal
.read_data_folder <- function(spat_method = NULL,
    data_dir = NULL,
    dir_items,
    data_to_use,
    load_format = NULL,
    require_data_DT,
    cores = NA,
    verbose = NULL,
    toplevel = 2L) {
    ch <- box_chars()

    # 0. check params
    if (is.null(data_dir) ||
        !dir.exists(data_dir)) {
        .gstop(
            .n = toplevel, "The full path to a", spat_method,
            "directory must be given."
        )
    }
    vmsg(.v = verbose, "A structured", spat_method, "directory will be used")
    if (!data_to_use %in% require_data_DT$workflow) {
        .gstop(
            .n = toplevel,
            "Data requirements for data_to_use not found in require_data_DT"
        )
    }

    # 1. detect items
    dir_items <- lapply_flex(dir_items, function(x) {
        Sys.glob(paths = file.path(data_dir, x))
    }, cores = cores)
    # (length = 1 if present, length = 0 if missing)
    dir_items_lengths <- lengths(dir_items)

    # 2. check directory contents
    vmsg(.v = verbose, "Checking directory contents...")

    for (item in names(dir_items)) {
        # IF ITEM FOUND

        if (dir_items_lengths[[item]] > 0) {
            # print found items if verbose = "debug"
            if (isTRUE(verbose)) {
                vmsg(
                    .v = verbose, .is_debug = TRUE,
                    .initial = paste0(ch$s, "> "),
                    item, " found"
                )
                for (item_i in seq_along(dir_items[[item]])) {
                    # print found item names
                    subItem <- gsub(
                        pattern = ".*/", replacement = "",
                        x = dir_items[[item]][[item_i]]
                    )
                    vmsg(
                        .v = verbose, .is_debug = TRUE,
                        .initial = paste0(ch$s, ch$s, ch$l, ch$h, ch$h),
                        subItem
                    )
                }
            }
        } else {
            # IF ITEM MISSING
            # necessary (error)
            # optional (warning)

            # data.table variables
            workflow <- needed <- filetype <- NULL


            require_data_DT <- require_data_DT[workflow == data_to_use, ]
            if (!is.null(load_format)) {
                require_data_DT <- require_data_DT[filetype == load_format, ]
            }

            if (item %in% require_data_DT[needed == TRUE, item]) {
                stop(item, " is missing")
            }
            if (item %in% require_data_DT[needed == FALSE, item]) {
                warning(item, "is missing (optional)")
            }
        }
    }

    # 3. select first path in list if multiple are detected
    if (any(dir_items_lengths > 1)) {
        warning(wrap_txt("Multiple matches for expected directory item(s).
                        First matching item selected"))

        multiples <- which(dir_items_lengths > 1)
        for (mult_i in multiples) {
            message(names(dir_items)[[mult_i]], "multiple matches found:")
            print(dir_items[[mult_i]])
            dir_items[[mult_i]] <- dir_items[[mult_i]][[1]]
        }
    }
    vmsg(.v = verbose, "Directory check done")

    return(dir_items)
}






.reader_fun_prints <- function(x, pre) {
    nfun <- length(x@calls)
    funs <- names(x@calls)
    if (nfun > 0L) {
        pre_funs <- format(c(pre, rep("", nfun - 1L)))
        for (i in seq_len(nfun)) {
            cat(pre_funs[i], " ", funs[i], "()\n", sep = "")
        }
    }
}

.filetype_prints <- function(x, pre) {
    nftype <- length(x@filetype)
    datatype <- format(names(x@filetype))
    pre_ftypes <- format(c(pre, rep("", nftype - 1L)))
    cat(
        sprintf(
            "%s %s -- %s\n",
            pre_ftypes,
            datatype,
            x@filetype
        ),
        sep = ""
    )
}

# pattern - list.files pattern to use to search for specific files/dirs
# warn - whether to warn when a pattern does not find any files
# first - whether to only return the first match
.detect_in_dir <- function(path, pattern, recursive = FALSE,
    platform, warn = TRUE, first = TRUE) {
    f <- list.files(path,
        pattern = pattern, recursive = recursive,
        full.names = TRUE
    )
    lenf <- length(f)
    if (lenf == 1L) {
        return(f)
    } # one match
    else if (lenf == 0L) { # no matches
        if (warn) {
            warning(
                sprintf(
                    "%s not detected in %s directory",
                    pattern,
                    platform
                ),
                call. = FALSE
            )
        }
        return(NULL)
    }

    # more than one match
    if (first) {
        return(f[[1L]])
    } else {
        return(f)
    }
}



# *---- object creation ----* ####






## Visium ####

#' @title Create a giotto object from 10x visium data
#' @name createGiottoVisiumObject
#' @description Create Giotto object directly from a 10X visium folder. Also
#' accepts visium H5 outputs.
#'
#' @param visium_dir path to the 10X visium directory [required]
#' @param expr_data raw or filtered data (see details)
#' @param gene_column_index which column index to select (see details)
#' @param h5_visium_path path to visium 10X .h5 file
#' @param h5_gene_ids gene names as symbols (default) or ensemble gene ids
#' @param h5_tissue_positions_path path to tissue locations (.csv file)
#' @param h5_image_png_path path to tissue .png file (optional). Image
#' autoscaling looks for matches in the filename for either "hires" or "lowres"
#' @param h5_json_scalefactors_path path to .json scalefactors (optional)
#' @param png_name select name of png to use (see details)
#' @param do_manual_adj deprecated
#' @param xmax_adj deprecated
#' @param xmin_adj deprecated
#' @param ymax_adj deprecated
#' @param ymin_adj deprecated
#' @param instructions list of instructions or output result from
#' \code{\link[GiottoClass]{createGiottoInstructions}}
#' @param cores how many cores or threads to use to read data if paths are
#' provided
#' @param expression_matrix_class class of expression matrix to use
#' (e.g. "dgCMatrix", "DelayedArray")
#' @param h5_file optional path to create an on-disk h5 file
#' @param verbose be verbose
#'
#' @returns giotto object
#' @details
#' If starting from a Visium 10X directory:
#' \itemize{
#'   \item{expr_data: raw will take expression data from
#'   raw_feature_bc_matrix and filter from filtered_feature_bc_matrix}
#'   \item{gene_column_index: which gene identifiers (names) to use if there
#'   are multiple columns (e.g. ensemble and gene symbol)}
#'   \item{png_name: by default the first png will be selected, provide the png
#'   name to override this (e.g. myimage.png)}
#'   \item{the file scalefactors_json.json will be detected automatically and
#'   used to attempt to align the data}
#' }
#'
#' If starting from a Visium 10X .h5 file
#' \itemize{
#'   \item{h5_visium_path: full path to .h5 file: /your/path/to/visium_file.h5}
#'   \item{h5_tissue_positions_path: full path to spatial locations file:
#'   /you/path/to/tissue_positions_list.csv}
#'   \item{h5_image_png_path: full path to png:
#'   /your/path/to/images/tissue_lowres_image.png}
#'   \item{h5_json_scalefactors_path: full path to .json file:
#'   /your/path/to/scalefactors_json.json}
#' }
#'
#' @export
createGiottoVisiumObject <- function(visium_dir = NULL,
    expr_data = c("raw", "filter"),
    gene_column_index = 1,
    h5_visium_path = NULL,
    h5_gene_ids = c("symbols", "ensembl"),
    h5_tissue_positions_path = NULL,
    h5_image_png_path = NULL,
    h5_json_scalefactors_path = NULL,
    png_name = NULL,
    do_manual_adj = FALSE, # deprecated
    xmax_adj = 0, # deprecated
    xmin_adj = 0, # deprecated
    ymax_adj = 0, # deprecated
    ymin_adj = 0, # deprecated
    instructions = NULL,
    expression_matrix_class = c("dgCMatrix", "DelayedArray"),
    h5_file = NULL,
    cores = NA,
    verbose = NULL) {
    # NSE vars
    barcode <- row_pxl <- col_pxl <- in_tissue <- array_row <- array_col <- NULL

    # handle deprecations
    img_dep_msg <- "The params 'do_manual_adj', 'xmax_adj', 'xmin_adj',
    'ymax_adj', 'ymin_adj' are no longer used.
    Please use the automated workflow."
    if (!isFALSE(do_manual_adj) ||
        xmax_adj != 0 ||
        xmin_adj != 0 ||
        ymax_adj != 0 ||
        ymin_adj != 0) {
        stop(wrap_txt(img_dep_msg))
    }

    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)


    # get arguments list for object creation
    if (!is.null(h5_visium_path)) {
        argslist <- .visium_read_h5(
            h5_visium_path = h5_visium_path, # expression matrix file
            h5_gene_ids = h5_gene_ids, # symbol or ensembl
            h5_tissue_positions_path = h5_tissue_positions_path,
            h5_image_png_path = h5_image_png_path,
            h5_json_scalefactors_path = h5_json_scalefactors_path,
            verbose = verbose
        )
    } else {
        argslist <- .visium_read_folder(
            visium_dir = visium_dir,
            expr_data = expr_data, # type of expression matrix to load
            gene_column_index = gene_column_index, # symbol or ensembl
            png_name = png_name,
            verbose = verbose
        )
    }

    # additional args to pass to object creation
    argslist$verbose <- verbose
    argslist$expression_matrix_class <- expression_matrix_class
    argslist$h5_file <- h5_file
    argslist$instructions <- instructions

    giotto_object <- do.call(.visium_create, args = argslist)

    return(giotto_object)
}








.visium_create <- function(
        expr_counts_path,
        h5_gene_ids = NULL, # h5
        gene_column_index = NULL, # folder
        tissue_positions_path,
        image_path = NULL,
        scale_json_path = NULL,
        png_name = NULL,
        instructions = NULL,
        expression_matrix_class = c("dgCMatrix", "DelayedArray"),
        h5_file = NULL,
        verbose = NULL) {
    # NSE vars
    barcode <- cell_ID <- row_pxl <- col_pxl <- in_tissue <- array_row <-
        array_col <- NULL

    # Assume path checking has been done

    # 1. expression
    if (!is.null(h5_gene_ids)) {
        expr_results <- get10Xmatrix_h5(
            path_to_data = expr_counts_path,
            gene_ids = h5_gene_ids
        )
    } else {
        expr_results <- get10Xmatrix(
            path_to_data = expr_counts_path,
            gene_column_index = gene_column_index
        )
    }

    # if expr_results is not a list, make it a list compatible with downstream
    if (!is.list(expr_results)) {
        expr_results <- list(
            "Gene Expression" = expr_results
        )
    }

    # format expected data into list to be used with readExprData()
    raw_matrix_list <- list("cell" = list("rna" = list(
        "raw" = expr_results[["Gene Expression"]]
    )))

    # add protein expression data to list if it exists
    if ("Antibody Capture" %in% names(expr_results)) {
        raw_matrix_list$cell$protein$raw <- expr_results[["Antibody Capture"]]
    }


    # 2. spatial locations
    spatial_results <- data.table::fread(tissue_positions_path)
    colnames(spatial_results) <- c(
        "barcode", "in_tissue", "array_row",
        "array_col", "col_pxl", "row_pxl"
    )
    spatial_results <- spatial_results[match(colnames(
        raw_matrix_list$cell[[1]]$raw
    ), barcode)]
    data.table::setnames(spatial_results, old = "barcode", new = "cell_ID")
    spatial_locs <- spatial_results[, .(cell_ID, row_pxl, -col_pxl)]
    # flip x and y
    colnames(spatial_locs) <- c("cell_ID", "sdimx", "sdimy")


    # 3. scalefactors (optional)
    json_info <- .visium_read_scalefactors(scale_json_path)


    # 4. image (optional)
    if (!is.null(image_path)) {
        visium_png_list <- .visium_image(
            image_path = image_path,
            json_info = json_info,
            verbose = verbose
        )
    }

    # 5. metadata
    meta_results <- spatial_results[
        , .(cell_ID, in_tissue, array_row, array_col)
    ]
    expr_types <- names(raw_matrix_list$cell)
    meta_list <- list()
    for (etype in expr_types) {
        meta_list[[etype]] <- meta_results
    }


    # 6. giotto object
    giotto_object <- createGiottoObject(
        expression = raw_matrix_list,
        spatial_locs = spatial_locs,
        instructions = instructions,
        cell_metadata = meta_list,
        images = visium_png_list
    )


    # 7. polygon information
    if (!is.null(json_info)) {
        visium_polygons <- .visium_spot_poly(
            spatlocs = spatial_locs,
            json_scalefactors = json_info
        )
        giotto_object <- setPolygonInfo(
            gobject = giotto_object,
            x = visium_polygons,
            centroids_to_spatlocs = FALSE,
            verbose = FALSE,
            initialize = TRUE
        )

        ms <- 65 / json_info$spot_diameter_fullres
        instructions(giotto_object, "micron_scale") <- ms
    }

    return(giotto_object)
}



# Find and check the filepaths within a structured visium directory
.visium_read_folder <- function(
        visium_dir = NULL,
        expr_data = c("raw", "filter"),
        gene_column_index = 1,
        png_name = NULL,
        verbose = NULL) {
    vmsg(.v = verbose, "A structured visium directory will be used")

    ## check arguments
    if (is.null(visium_dir)) {
        .gstop("visium_dir needs to be a path to a visium directory")
    }
    visium_dir <- path.expand(visium_dir)
    if (!dir.exists(visium_dir)) .gstop(visium_dir, " does not exist!")
    expr_data <- match.arg(expr_data, choices = c("raw", "filter"))


    ## 1. check expression
    expr_counts_path <- switch(expr_data,
        "raw" = file.path(visium_dir, "raw_feature_bc_matrix"),
        "filter" = file.path(visium_dir, "filtered_feature_bc_matrix")
    )
    if (!file.exists(expr_counts_path)) {
        .gstop(expr_counts_path, "does not exist!")
    }


    ## 2. check spatial locations
    spatial_dir <- file.path(visium_dir, "spatial")
    
    spatial_files <- list.files(spatial_dir)
    
    tissue_positions_path <- file.path(
        spatial_dir,
        spatial_files[grep("tissue_positions*", spatial_files)])

    ## 3. check spatial image
    if (is.null(png_name)) {
        png_list <- list.files(spatial_dir, pattern = "*.png")
        png_name <- png_list[1]
    }
    png_path <- file.path(spatial_dir, png_name)
    if (!file.exists(png_path)) .gstop(png_path, " does not exist!")


    ## 4. check scalefactors
    scalefactors_path <- file.path(
        spatial_dir,
        spatial_files[grep("scalefactors_json.json", spatial_files)])
    
    if (!file.exists(scalefactors_path)) {
        .gstop(scalefactors_path, "does not exist!")
    }


    list(
        expr_counts_path = expr_counts_path,
        gene_column_index = gene_column_index,
        tissue_positions_path = tissue_positions_path,
        image_path = png_path,
        scale_json_path = scalefactors_path
    )
}



.visium_read_h5 <- function(
        h5_visium_path = h5_visium_path, # expression matrix
        h5_gene_ids = h5_gene_ids,
        h5_tissue_positions_path = h5_tissue_positions_path,
        h5_image_png_path = h5_image_png_path,
        h5_json_scalefactors_path = h5_json_scalefactors_path,
        verbose = NULL) {
    # 1. filepaths
    vmsg(
        .v = verbose,
        "A path to an .h5 10X file was provided and will be used"
    )
    if (!file.exists(h5_visium_path)) {
        .gstop("The provided path ", h5_visium_path, " does not exist")
    }
    if (is.null(h5_tissue_positions_path)) {
        .gstop("A path to the tissue positions (.csv) needs to be provided to
                h5_tissue_positions_path")
    }
    if (!file.exists(h5_tissue_positions_path)) {
        .gstop(
            "The provided path ", h5_tissue_positions_path,
            " does not exist"
        )
    }
    if (!is.null(h5_image_png_path)) {
        if (!file.exists(h5_image_png_path)) {
            .gstop(
                "The provided h5 image path ", h5_image_png_path,
                "does not exist.
            Set to NULL to exclude or provide the correct path."
            )
        }
    }
    if (!is.null(h5_json_scalefactors_path)) {
        if (!file.exists(h5_json_scalefactors_path)) {
            warning(wrap_txt(
                "No file found at h5_json_scalefactors_path.
                Scalefactors are needed for proper image alignment and
                polygon generation"
            ))
        }
    }

    list(
        expr_counts_path = h5_visium_path,
        h5_gene_ids = h5_gene_ids,
        tissue_positions_path = h5_tissue_positions_path,
        image_path = h5_image_png_path,
        scale_json_path = h5_json_scalefactors_path
    )
}









# Visium Polygon Creation

#' @title Add Visium Polygons to Giotto Object
#' @name addVisiumPolygons
#' @param gobject Giotto Object created with visium data, containing spatial
#' locations corresponding to spots
#' @param scalefactor_path path to scalefactors_json.json Visium output
#' @returns Giotto Object with to-scale circular polygons added at each spatial
#' location
#' @details
#' Adds circular giottoPolygons to the spatial_info slot of a Giotto Object
#' for the "cell" spatial unit.
#' @export
addVisiumPolygons <- function(gobject,
    scalefactor_path = NULL) {
    assert_giotto(gobject)

    visium_spat_locs <- getSpatialLocations(
        gobject = gobject,
        spat_unit = "cell"
    )

    scalefactors_list <- .visium_read_scalefactors(
        json_path = scalefactor_path
    )

    visium_polygons <- .visium_spot_poly(
        spatlocs = visium_spat_locs,
        json_scalefactors = scalefactors_list
    )

    gobject <- addGiottoPolygons(
        gobject = gobject,
        gpolygons = list(visium_polygons)
    )

    return(gobject)
}





#' @title Read Visium ScaleFactors
#' @name .visium_read_scalefactors
#' @param json_path path to scalefactors_json.json for Visium experimental data
#' @returns scalefactors within the provided json file as a named list,
#' or NULL if not discovered
#' @details asserts the existence of and reads in a .json file
#' containing scalefactors for Visium data in the expected format.
#' Returns NULL if no path is provided or if the file does not exist.
#' @keywords internal
.visium_read_scalefactors <- function(json_path = NULL) {
    if (!checkmate::test_file_exists(json_path)) {
        if (!is.null(json_path)) {
            warning("scalefactors not discovered at: \n",
                json_path,
                call. = FALSE
            )
        }
        return(NULL)
    }

    json_scalefactors <- read_json(json_path)

    # Intial assertion that json dimensions are appropriate
    checkmate::assert_list(
        x = json_scalefactors,
        types = "numeric",
        min.len = 4L,
        max.len = 5L
    )

    expected_json_names <- c(
        "regist_target_img_scalef", # NEW as of 2023
        "spot_diameter_fullres",
        "tissue_hires_scalef",
        "fiducial_diameter_fullres",
        "tissue_lowres_scalef"
    )

    # Visium assay with chemistry v2 contains an additional
    # keyword in the json file
    new_format_2023 <- checkmate::test_list(
        x = json_scalefactors,
        types = "numeric",
        len = 5L
    )

    # If the scalefactors are of size 4 (older assay), clip the new keyword
    if (!new_format_2023) expected_json_names <- expected_json_names[2:5]

    if (!setequal(names(json_scalefactors), expected_json_names)) {
        warning(GiottoUtils::wrap_txt(
            "h5 scalefactors json names differ from expected.
            [Expected]:", expected_json_names, "\n",
            "[Actual]:", names(json_scalefactors)
        ))
    }

    return(json_scalefactors)
}


#' @title Calculate Pixel to Micron Scalefactor
#' @name visium_micron_scalefactor
#' @param json_scalefactors list of scalefactors from
#' .visium_read_scalefactors()
#' @returns scale factor for converting pixel to micron
#' @details
#' Calculates pixel to micron scalefactor.
#' Visium xy coordinates are based on the fullres image
#' The values provided are directly usable for generating polygon information
#' or calculating the micron size relative to spatial coordinates for this set
#' of spatial information.
#' @keywords internal
.visium_micron_scale <- function(json_scalefactors) {
    # visium spots diameter                       : 55 micron
    # diameter of a spot at this spatial scaling  : scalefactor_list$spot_diameter_fullres
    px_to_micron <- 55 / json_scalefactors$spot_diameter_fullres
    return(px_to_micron)
}


#' @title Create Polygons for Visium Data
#' @name .visium_spot_poly
#' @param spatlocs spatial locations data.table or `spatLocsObj` containing
#' centroid locations of visium spots
#' @param json_scalefactors list of scalefactors from
#' .visium_read_scalefactors()
#' @returns giottoPolygon object
#' @details
#' Creates circular polygons for spatial representation of
#' Visium spots.
#' @keywords internal
#' @md
.visium_spot_poly <- function(spatlocs = NULL,
    json_scalefactors) {
    if (inherits(spatlocs, "spatLocsObj")) {
        spatlocs <- spatlocs[]
    }

    spot_adj <- 55 / 65

    vis_spot_poly <- GiottoClass::circleVertices(
        radius = json_scalefactors$spot_diameter_fullres / 2 * spot_adj,
        npoints = 100
    )

    GiottoClass::polyStamp(
        stamp_dt = vis_spot_poly,
        spatlocs = spatlocs,
        verbose = FALSE
    ) %>%
        createGiottoPolygonsFromDfr(
            calc_centroids = TRUE,
            verbose = FALSE
        )
}






# json_info expects the list read output from .visium_read_scalefactors
# image_path should be expected to be full filepath
# should only be used when do_manual_adj (deprecated) is FALSE
.visium_image <- function(
        image_path,
        json_info = NULL,
        micron_scale = FALSE,
        verbose = NULL) {
    # assume image already checked
    vmsg(.v = verbose, .initial = " - ", "found image")

    # 1. determine image scalefactor to use ---------------------------------- #
    if (!is.null(json_info)) checkmate::assert_list(json_info)
    png_name <- basename(image_path) # used for name pattern matching only

    if (is.null(json_info)) { # if none provided
        warning(wrap_txt(
            "No scalefactors json info provided.
            Visium image scale_factor defaulting to 1"
        ))
        scale_factor <- 1
    } else { # if provided

        scale_factor <- NULL # initial value

        # determine type of visium image
        visium_img_type <- NULL
        possible_types <- c("lowres", "hires")
        for (img_type in possible_types) {
            if (grepl(img_type, png_name)) visium_img_type <- img_type
        }

        if (is.null(visium_img_type)) { # if not recognized visium image type
            .gstop(
                "\'image_path\' filename did not partial match either
                \'lowres\' or \'hires\'. Ensure specified image is either the
                Visium lowres or hires image and rename it accordingly"
            )
        }

        vmsg(
            .v = verbose, .initial = " - ",
            "found scalefactors. attempting automatic alignment for the",
            str_quote(visium_img_type), "image\n\n"
        )

        scale_factor <- switch(visium_img_type,
            "lowres" = json_info[["tissue_lowres_scalef"]],
            "hires" = json_info[["tissue_hires_scalef"]]
        )
    }

    if (isTRUE(micron_scale)) {
        scale_factor <- scale_factor * .visium_micron_scale(json_info)
    }

    # 2. create image -------------------------------------------------------- #
    visium_img <- createGiottoLargeImage(
        raster_object = image_path,
        name = "image",
        negative_y = TRUE,
        scale_factor = (1 / scale_factor)
    )

    visium_img_list <- list(visium_img)
    names(visium_img_list) <- c("image")

    return(visium_img_list)
}











## MERSCOPE ####


#' @title Create Vizgen MERSCOPE largeImage
#' @name createMerscopeLargeImage
#' @description
#' Read MERSCOPE stitched images as giottoLargeImage. Images will also be
#' transformed to match the spatial coordinate reference system of the paired
#' points and polygon data.
#' @param image_file character. Path to one or more MERSCOPE images to load
#' @param transforms_file character. Path to MERSCOPE transforms file. Usually
#' in the same folder as the images and named
#' 'micron_to_mosaic_pixel_transform.csv'
#' @param name character. name to assign the image. Multiple should be provided
#' if image_file is a list.
#' @param flip_axis character. Axis along which to reflect the image after the
#' spatial shift: \code{"none"} (default), \code{"y"}, \code{"x"}, or
#' \code{"both"}. The reflection origin is the midpoint of the image extent in
#' micron space.
#' @returns giottoLargeImage
#' @export
createMerscopeLargeImage <- function(image_file,
    transforms_file,
    name = "image",
    flip_axis = c("none", "y", "x", "both")) {
    flip_axis <- match.arg(flip_axis)

    if (flip_axis != "none") {
        message(
            "The image is flipped around its own spatial midpoint (microns).",
            " If you also used flip_axis (!= 'none) in createGiottoMerscopeObject,",
            " verify that the image aligns with the polygons and transcripts --",
            " the flip centers may differ if the image does not cover the same",
            " area as the transcript coordinates."
        )
    }

    checkmate::assert_character(transforms_file)
    tfsDT <- data.table::fread(transforms_file)
    if (inherits(image_file, "character")) {
        image_file <- as.list(image_file)
    }
    checkmate::assert_list(image_file)

    scalef  <- c(1 / tfsDT[[1, 1]], 1 / tfsDT[[2, 2]])
    x_shift <- -tfsDT[[1, 3]] / tfsDT[[1, 1]]
    y_shift <- -tfsDT[[2, 3]] / tfsDT[[2, 2]]

    out <- lapply(seq_along(image_file), function(i) {
        gimg <- createGiottoLargeImage(
            raster_object = image_file[[i]],
            name          = name[[i]],
            scale_factor  = scalef,
            negative_y    = FALSE
        )
        gimg <- spatShift(gimg, dx = x_shift, dy = y_shift)

        if (flip_axis != "none") {
            if (flip_axis %in% c("y", "both")) {
                gimg  <- flip(gimg, direction = "vertical")
            }
            if (flip_axis %in% c("x", "both")) {
                gimg  <- flip(gimg, direction = "horizontal")
            }
        }

        return(gimg)
    })

    if (length(out) == 1L) {
        out <- out[[1L]]
    }

    return(out)
}







#' @title Create Vizgen MERSCOPE Giotto Object
#' @name createGiottoMerscopeObject
#' @description Given the path to a MERSCOPE experiment directory, creates a
#' Giotto object. Supports both Parquet and HDF5 cell boundary formats (v2).
#' Features automatic 2D vs 3D architecture detection with vertex validation.
#' @param merscope_dir full path to the exported merscope directory
#' @param version integer. Version of logic to use (1 = Legacy/Original,
#' 2 = New/Optimized)
#' @param data_to_use which of either the 'subcellular' or 'aggregate'
#' information to use for object creation
#' @param FOVs which FOVs to use when building the subcellular object.
#' (default is NULL). NULL loads all FOVs (very slow)
#' @param polygon_format format of the boundary files, either 'parquet' or
#' 'hdf5' (v2 only)
#' @param poly_z_indices which z-indices to use for the polygons (Default 0:6
#' for v2, 1:7 for v1)
#' @param flip_axis character. \strong{Version 2 only.} Axis along which to
#' reflect polygon vertices and transcript coordinates after loading:
#' \code{"none"} (default, data read as-is), \code{"y"}, \code{"x"}, or
#' \code{"both"}. The reflection origin is the midpoint of the transcript
#' coordinate range so both layers stay registered. Silently forced to
#' \code{"none"} for \code{version = 1}, which handles orientation internally
#' via \code{polygon_mask_list_params} (\code{flip_vertical = TRUE}).
#' @param spat_unit character. Name assigned to the cell spatial unit (default
#' \code{"cell"}). Controls the polygon layer name when a single Z-plane is
#' present and serves as the default \code{new_spat_unit} in
#' \code{aggregate_stack_param} when that field is not already supplied.
#' @param calculate_overlap whether to run \code{\link{calculateOverlapRaster}}
#' @param overlap_to_matrix whether to run \code{\link{overlapToMatrix}}
#' @param aggregate_stack whether to run \code{\link{aggregateStacks}}
#' @param aggregate_stack_param params to pass to \code{\link{aggregateStacks}}
#' @param split_keyword character vector of keywords to separate feature
#' detections (e.g. "Blank")
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @returns a giotto object
#' @export
createGiottoMerscopeObject <- function(merscope_dir,
    version = 1,
    data_to_use = c("subcellular", "aggregate"),
    FOVs = NULL,
    polygon_format = c("parquet", "hdf5"),
    poly_z_indices = NULL,
    flip_axis = c("none", "y", "x", "both"),
    spat_unit = "cell",
    calculate_overlap = TRUE,
    overlap_to_matrix = TRUE,
    aggregate_stack = TRUE,
    aggregate_stack_param = list(
        summarize_expression = "sum",
        summarize_locations = "mean",
        new_spat_unit = "cell"
    ),
    split_keyword = c("Blank"),
    instructions = NULL,
    cores = NA,
    verbose = TRUE) {

    # 1. version-specific defaults and validation
    merscope_dir <- path.expand(merscope_dir)
    data_to_use  <- match.arg(arg = data_to_use, choices = c("subcellular", "aggregate"))
    flip_axis    <- match.arg(flip_axis)

    if (version == 1 && flip_axis != "none") {
        message(
            "'flip_axis' is only available for version = 2. ",
            "Version 1 handles orientation internally via ",
            "polygon_mask_list_params (flip_vertical = TRUE). ",
            "Setting flip_axis = 'none'."
        )
        flip_axis <- "none"
    }

    if (isTRUE(verbose)) {
        if (inherits(future::plan(), "sequential")) {
            message(
                "Running sequentially. For better performance consider:\n",
                "  future::plan(future::multisession, workers = <n>)"
            )
        }
    }
    # suppress Giotto's internal sequential warning (surfaced above)
    options("giotto.warn_sequential" = FALSE)

    if (version == 1) {
        if (is.null(poly_z_indices)) poly_z_indices <- seq(from = 1, to = 7)
        poly_z_indices <- as.integer(poly_z_indices)
        if (any(poly_z_indices < 1)) {
            stop(wrap_txt(
                "Version 1 requires 1-based indexing (poly_z_indices >= 1).",
                errWidth = TRUE
            ))
        }
        polygon_format <- "hdf5"
        fovs <- NULL
    } else {
        if (is.null(poly_z_indices)) poly_z_indices <- 0:6
        poly_z_indices <- as.integer(poly_z_indices)
        if (any(poly_z_indices < 0)) {
            stop(wrap_txt(
                "Version 2 requires 0-based indexing (poly_z_indices >= 0).",
                errWidth = TRUE
            ))
        }
        polygon_format <- match.arg(arg = polygon_format, choices = c("parquet", "hdf5"))
        fovs <- FOVs
    }

    # 2. directory scanning
    dir_items <- .read_merscope_folder(
        merscope_dir   = merscope_dir,
        data_to_use    = data_to_use,
        polygon_format = polygon_format,
        version        = version,
        cores          = cores,
        verbose        = verbose
    )

    # 3. data loading
    data_list <- .load_merscope_folder(
        dir_items      = dir_items,
        data_to_use    = data_to_use,
        polygon_format = polygon_format,
        poly_z_indices = poly_z_indices,
        fovs           = fovs,
        version        = version,
        cores          = cores,
        verbose        = verbose
    )

    # 4. version 2: optimized 2D vs 3D auto-detection
    if (version == 2 && data_to_use == "subcellular" && !is.null(data_list$poly_info)) {

        # bypass aggregation if only 1 Z-plane loaded
        if (length(data_list$poly_info) == 1 && isTRUE(aggregate_stack)) {
            if (isTRUE(verbose)) {
                message(
                    "--> [Auto-Correction] Only 1 Z-plane loaded.",
                    " Setting 'aggregate_stack = FALSE'."
                )
            }
            aggregate_stack <- FALSE
        }

        if (length(data_list$poly_info) > 1) {
            ids_zA <- data_list$poly_info[[1]]@unique_ID_cache
            ids_zB <- data_list$poly_info[[2]]@unique_ID_cache

            if (isTRUE(setequal(ids_zA, ids_zB))) {
                sv_A   <- slot(data_list$poly_info[[1]], "spatVector")
                sv_B   <- slot(data_list$poly_info[[2]], "spatVector")
                crds_A <- head(terra::crds(sv_A), 100)
                crds_B <- head(terra::crds(sv_B), 100)

                if (isTRUE(all.equal(crds_A, crds_B))) {
                    if (isTRUE(verbose)) {
                        message(
                            "--> [Auto-Detection] Identical IDs AND vertices",
                            " detected (Replicated 2D Data)."
                        )
                        message(
                            "--> [Auto-Correction] Dropping redundant Z-planes.",
                            " Setting 'aggregate_stack = FALSE'."
                        )
                    }
                    data_list$poly_info <- data_list$poly_info[1]
                    poly_z_indices      <- poly_z_indices[1]
                    aggregate_stack     <- FALSE
                }
            }
        }
    }

    # 5. object creation
    if (data_to_use == "subcellular") {
        merscope_gobject <- .createGiottoMerscopeObject_subcellular(
            data_list             = data_list,
            version               = version,
            flip_axis             = flip_axis,
            spat_unit             = spat_unit,
            calculate_overlap     = calculate_overlap,
            overlap_to_matrix     = overlap_to_matrix,
            aggregate_stack       = aggregate_stack,
            aggregate_stack_param = aggregate_stack_param,
            split_keyword         = split_keyword,
            instructions          = instructions,
            cores                 = cores,
            verbose               = verbose
        )
    } else {
        merscope_gobject <- .createGiottoMerscopeObject_aggregate(
            data_list = data_list,
            cores     = cores,
            verbose   = verbose
        )
    }

    return(merscope_gobject)
}




#' @describeIn createGiottoMerscopeObject Create giotto object with
#' 'subcellular' workflow
#' @param data_list list of loaded data from \code{\link{load_merscope_folder}}
#' @keywords internal
.createGiottoMerscopeObject_subcellular <- function(data_list,
    version               = 1,
    flip_axis             = c("none", "y", "x", "both"),
    spat_unit             = "cell",
    calculate_overlap     = TRUE,
    overlap_to_matrix     = TRUE,
    aggregate_stack       = TRUE,
    aggregate_stack_param = list(
        summarize_expression = "sum",
        summarize_locations  = "mean"
    ),
    split_keyword         = c("Blank"),
    instructions          = NULL,
    cores                 = NA,
    verbose               = TRUE) {

    tx_dt     <- data_list$tx_dt
    poly_info <- data_list$poly_info
    flip_axis <- match.arg(flip_axis)

    # optional spatial flip of transcripts and polygons (version 2 only)
    if (version == 2 && flip_axis != "none") {
        x_mid <- mean(range(tx_dt$global_x))
        y_mid <- mean(range(tx_dt$global_y))

        if (flip_axis %in% c("y", "both")) tx_dt[, global_y := y_mid * 2 - global_y]
        if (flip_axis %in% c("x", "both")) tx_dt[, global_x := x_mid * 2 - global_x]

        poly_info <- lapply(poly_info, function(gpoly) {
            if (flip_axis %in% c("y", "both")) gpoly <- GiottoClass::flip(gpoly, y0 = y_mid)
            if (flip_axis %in% c("x", "both")) gpoly <- GiottoClass::flip(gpoly, x0 = x_mid)
            gpoly
        })
    }

    split_list <- if (!is.null(split_keyword)) as.list(split_keyword) else NULL

    gpoints <- GiottoClass::createGiottoPoints(
        tx_dt,
        x_colname       = "global_x",
        y_colname       = "global_y",
        feat_ID_colname = "gene",
        split_keyword   = split_list
    )

    gpoints_list <- list("rna" = gpoints[["rna"]])
    if (!is.null(split_list)) {
        for (kw in unlist(split_list)) {
            if (kw %in% names(gpoints)) {
                list_name                 <- ifelse(kw == "Blank", "neg_probe", kw)
                gpoints_list[[list_name]] <- gpoints[[kw]]
            }
        }
    }

    gpolygons <- if (length(poly_info) == 1) {
        setNames(list(poly_info[[1]]), spat_unit)
    } else {
        poly_info
    }

    if (version == 1) {
        z_sub <- GiottoClass::createGiottoObjectSubcellular(
            gpoints      = gpoints_list,
            gpolygons    = gpolygons,
            instructions = instructions,
            cores        = cores,
            polygon_mask_list_params = list(
                mask_method           = "guess",
                flip_vertical         = TRUE,
                flip_horizontal       = FALSE,
                shift_horizontal_step = FALSE
            )
        )
    } else {
        z_sub <- GiottoClass::createGiottoObjectSubcellular(
            gpoints      = gpoints_list,
            gpolygons    = gpolygons,
            instructions = instructions,
            cores        = cores,
            verbose      = verbose
        )
    }

    # STEP 1: Aggregate stacks FIRST (before overlap)
    if (isTRUE(aggregate_stack)) {
        if (isTRUE(verbose)) message("Aggregating 3D stacks...")
        agg_params         <- aggregate_stack_param
        agg_params$gobject <- z_sub
        if (is.null(agg_params$spat_units))    agg_params$spat_units    <- names(gpolygons)
        if (is.null(agg_params$feat_type))     agg_params$feat_type     <- "rna"
        if (is.null(agg_params$new_spat_unit)) agg_params$new_spat_unit <- spat_unit
        z_sub <- do.call(GiottoClass::aggregateStacks, agg_params)
    }

    # STEP 2: Calculate overlap on final polygon layer(s)
    # If aggregate_stack ran, overlap is computed on the new aggregated unit.
    # If not, it falls back to the original gpolygons names.
    if (isTRUE(calculate_overlap)) {

        final_poly_names <- if (isTRUE(aggregate_stack)) {
            agg_params$new_spat_unit
        } else {
            names(gpolygons)
        }

        for (poly_name in final_poly_names) {
            if (isTRUE(verbose)) message("Calculating overlap for polygon layer: ", poly_name)

            z_sub <- GiottoClass::calculateOverlap(
                gobject = z_sub, spat_info = poly_name, feat_info = "rna"
            )
            if (isTRUE(overlap_to_matrix)) {

                z_sub <- GiottoClass::overlapToMatrix(
                    gobject = z_sub, spat_info = poly_name, feat_info = "rna", name = "raw"
                )
            }
        }
    }

    return(z_sub)
}




#' @describeIn createGiottoMerscopeObject Create giotto object with 'aggregate'
#' workflow
#' @param data_list list of loaded data from \code{\link{load_merscope_folder}}
#' @keywords internal
.createGiottoMerscopeObject_aggregate <- function(data_list,
    cores = NA,
    verbose = TRUE) {
    # unpack data_list
    micronToPixelScale <- data_list$micronToPixelScale
    expr_dt <- data_list$expr_dt
    cell_meta <- data_list$expr_mat
    image_list <- data_list$images

    # split expr_dt by expression and blank

    # feat_id_all =
}




## Spatial Genomics ####

#' @title Create Spatial Genomics Giotto Object
#' @name createSpatialGenomicsObject
#' @param sg_dir full path to the exported Spatial Genomics directory
#' @param instructions new instructions
#' (e.g. result from createGiottoInstructions)
#' @returns giotto object
#' @description Given the path to a Spatial Genomics data directory, creates a
#' Giotto object.
#' @export
createSpatialGenomicsObject <- function(sg_dir = NULL,
    instructions = NULL) {
    # Find files in Spatial Genomics directory
    dapi <- list.files(sg_dir, full.names = TRUE, pattern = "DAPI")
    mask <- list.files(sg_dir, full.names = TRUE, pattern = "mask")
    tx <- list.files(sg_dir, full.names = TRUE, pattern = "transcript")
    # Create Polygons
    gpoly <- createGiottoPolygonsFromMask(
        mask,
        shift_vertical_step = FALSE,
        shift_horizontal_step = FALSE,
        flip_horizontal = FALSE,
        flip_vertical = FALSE
    )
    # Create Points
    tx <- data.table::fread(tx)
    gpoints <- createGiottoPoints(tx)
    dim(tx)
    # Create object and add image
    gimg <- createGiottoLargeImage(dapi, use_rast_ext = TRUE)
    sg <- createGiottoObjectSubcellular(
        gpoints = list("rna" = gpoints),
        gpolygons = list("cell" = gpoly),
        instructions = instructions
    )
    sg <- addGiottoImage(sg, images = list(image = gimg))
    # Return SG object
    return(sg)
}















# *---- folder reading and detection ----* ####


#' @describeIn read_data_folder Read a structured MERSCOPE folder
#' @keywords internal
.read_merscope_folder <- function(merscope_dir,
    data_to_use,
    polygon_format = "hdf5",
    version = 1,
    cores = NA,
    verbose = NULL) {
    # boundary pattern: v2 supports Parquet in addition to HDF5
    boundary_pattern <- if (version == 2 && polygon_format == "parquet") {
        "*cell_boundaries.parquet*"
    } else {
        "*cell_boundaries*"
    }

    dir_items <- list(
        `boundary info`       = boundary_pattern,
        `image info`          = "*images*",
        `cell feature matrix` = "*cell_by_gene*",
        `cell metadata`       = "*cell_metadata*",
        `raw transcript info` = "*transcripts*"
    )

    sub_reqs <- data.table::data.table(
        workflow = c("subcellular"),
        item     = c(
            "boundary info",
            "raw transcript info",
            "image info",
            "cell feature matrix",
            "cell metadata"
        ),
        needed = c(TRUE, TRUE, FALSE, FALSE, FALSE)
    )

    agg_reqs <- data.table::data.table(
        workflow = c("aggregate"),
        item     = c(
            "boundary info",
            "raw transcript info",
            "image info",
            "cell feature matrix",
            "cell metadata"
        ),
        needed = c(FALSE, FALSE, FALSE, TRUE, TRUE)
    )

    return(.read_data_folder(
        spat_method     = "MERSCOPE",
        data_dir        = merscope_dir,
        dir_items       = dir_items,
        data_to_use     = data_to_use,
        require_data_DT = rbind(sub_reqs, agg_reqs),
        cores           = cores,
        verbose         = verbose
    ))
}










# * ---- folder loading ---- * ####



## MERSCOPE ####

#' @title Load MERSCOPE data from folder
#' @name load_merscope_folder
#' @param dir_items list of full filepaths from
#' \code{\link{.read_merscope_folder}}
#' @inheritParams createGiottoMerscopeObject
#' @returns list of loaded-in MERSCOPE data
NULL

#' @rdname load_merscope_folder
#' @keywords internal
.load_merscope_folder <- function(dir_items,
    data_to_use,
    polygon_format,
    fovs           = NULL,
    poly_z_indices = 0:6,
    version        = 1,
    cores          = NA,
    verbose        = TRUE) {
    # 1. load workflow-specific data
    if (data_to_use == "subcellular") {
        data_list <- .load_merscope_folder_subcellular(
            dir_items      = dir_items,
            data_to_use    = data_to_use,
            polygon_format = polygon_format,
            fovs           = fovs,
            poly_z_indices = poly_z_indices,
            version        = version,
            cores          = cores,
            verbose        = verbose
        )
    } else {
        data_list <- .load_merscope_folder_aggregate(
            dir_items   = dir_items,
            data_to_use = data_to_use,
            cores       = cores,
            verbose     = verbose
        )
    }

    # 2. load images if available
    if (!is.null(dir_items$`image info`)) {
        ## micron to px scaling factor
        micronToPixelScale <- Sys.glob(
            paths = file.path(dir_items$`image info`, "*micron_to_mosaic_pixel_transform*")
        )[[1]]
        data_list$micronToPixelScale <- data.table::fread(micronToPixelScale, nThread = cores)

        ## staining images
        images_filenames   <- list.files(dir_items$`image info`, pattern = ".tif")
        bound_stains_types <- unique(sapply(strsplit(images_filenames, "_"), `[`, 2))

        data_list$images <- GiottoUtils:::lapply_flex(bound_stains_types, function(stype) {
            img_paths <- Sys.glob(paths = file.path(
                dir_items$`image info`, paste0("*", stype, "*")
            ))
            GiottoUtils:::lapply_flex(img_paths, createGiottoLargeImage, cores = cores)
        }, cores = cores)
    }

    return(data_list)
}



#' @describeIn load_merscope_folder Load items for 'subcellular' workflow
#' @keywords internal
.load_merscope_folder_subcellular <- function(dir_items,
    data_to_use,
    polygon_format,
    cores          = NA,
    poly_z_indices = 0L:6L,
    version        = 1,
    verbose        = TRUE,
    fovs           = NULL) {
    if (isTRUE(verbose)) message("Loading transcript level info...")
    if (is.null(fovs)) {
        tx_dt <- data.table::fread(dir_items$`raw transcript info`, nThread = cores)
    } else {
        tx_dt <- GiottoUtils::read_colmatch(
            file            = dir_items$`raw transcript info`,
            col             = "fov",
            values_to_match = fovs,
            nThread         = cores
        )
    }
    tx_dt[, c("x", "y") := NULL]
    data.table::setcolorder(tx_dt, c("gene", "global_x", "global_y", "global_z"))

    if (isTRUE(verbose)) message("Loading polygon info...")
    if (version == 2 && polygon_format == "parquet") {
        poly_info <- readPolygonVizgenParquet(
            file           = dir_items$`boundary info`,
            z_index        = poly_z_indices,
            calc_centroids = TRUE,
            verbose        = verbose
        )
    } else {
        poly_info <- readPolygonFilesVizgenHDF5(
            boundaries_path = dir_items$`boundary info`,
            z_indices       = poly_z_indices,
            fovs            = fovs
        )
    }

    return(list(
        "poly_info"          = poly_info,
        "tx_dt"              = tx_dt,
        "micronToPixelScale" = NULL,
        "expr_dt"            = NULL,
        "cell_meta"          = NULL,
        "images"             = NULL
    ))
}



#' @describeIn load_merscope_folder Load items for 'aggregate' workflow
#' @keywords internal
.load_merscope_folder_aggregate <- function(dir_items,
    data_to_use,
    cores = NA,
    verbose = TRUE) {
    # metadata is polygon-related measurements
    vmsg("Loading cell metadata...", .v = verbose)
    cell_metadata_file <- data.table::fread(
        dir_items$`cell metadata`,
        nThread = cores
    )

    vmsg("Loading expression matrix", .v = verbose)
    expr_dt <- data.table::fread(
        dir_items$`cell feature matrix`,
        nThread = cores
    )


    data_list <- list(
        "poly_info" = NULL,
        "tx_dt" = NULL,
        "micronToPixelScale" = NULL,
        "expr_dt" = expr_dt,
        "cell_meta" = cell_metadata_file,
        "images" = NULL
    )
}

















## ArchR ####

#' Create an ArchR project and run LSI dimension reduction
#'
#' @param fragmentsPath A character vector containing the paths to the input
#' files to use to generate the ArrowFiles.
#' These files can be in one of the following formats: (i) scATAC tabix files,
#' (ii) fragment files, or (iii) bam files.
#' @param genome A string indicating the default genome to be used for all ArchR
#' functions. Currently supported values include "hg19","hg38","mm9", and
#' "mm10".
#' This value is stored as a global environment variable, not part of the
#' ArchRProject.
#' This can be overwritten on a per-function basis using the given function's
#' geneAnnotationand genomeAnnotation parameter. For something other than one of
#' the currently supported, see createGeneAnnnotation() and
#' createGenomeAnnnotation()
#' @param createArrowFiles_params list of parameters passed to
#' `ArchR::createArrowFiles`
#' @param ArchRProject_params list of parameters passed to `ArchR::ArchRProject`
#' @param addIterativeLSI_params list of parameters passed to
#' `ArchR::addIterativeLSI`
#' @param threads number of threads to use. Default = `ArchR::getArchRThreads()`
#' @param force Default = FALSE
#' @param verbose Default = TRUE
#'
#' @returns An ArchR project with GeneScoreMatrix, TileMatrix, and
#' TileMatrix-based LSI
#' @export
createArchRProj <- function(fragmentsPath,
    genome = c("hg19", "hg38", "mm9", "mm10"),
    createArrowFiles_params = list(
        sampleNames = "sample1",
        minTSS = 0,
        minFrags = 0,
        maxFrags = 1e+07,
        minFragSize = 10,
        maxFragSize = 2000,
        offsetPlus = 0,
        offsetMinus = 0,
        TileMatParams = list(tileSize = 5000)
    ),
    ArchRProject_params = list(
        outputDirectory = getwd(),
        copyArrows = FALSE
    ),
    addIterativeLSI_params = list(),
    threads = ArchR::getArchRThreads(),
    force = FALSE,
    verbose = TRUE) {
    if (!requireNamespace("ArchR")) {
        message('ArchR is needed. Install the package using
                remotes::install_github("GreenleafLab/ArchR")')
    }

    ## Add reference genome
    message("Loading reference genome")
    ArchR::addArchRGenome(genome)

    # Creating Arrow Files
    message("Creating Arrow files")
    ArrowFiles <- do.call(
        ArchR::createArrowFiles,
        c(
            inputFiles = fragmentsPath,
            verbose = verbose,
            force = force,
            createArrowFiles_params
        )
    )

    # Creating an ArchRProject
    message("Creating ArchRProject")
    proj <- do.call(
        ArchR::ArchRProject,
        c(list(ArrowFiles = ArrowFiles),
            threads = threads,
            ArchRProject_params
        )
    )

    # Data normalization and dimensionality reduction
    message("Running dimension reduction")
    proj <- do.call(
        ArchR::addIterativeLSI,
        c(
            ArchRProj = proj,
            verbose = verbose,
            name = "IterativeLSI",
            threads = threads,
            force = force,
            addIterativeLSI_params
        )
    )
}

#' Create a Giotto object from an ArchR project
#'
#' @param archRproj ArchR project
#' @param expression expression information
#' @param expression_feat Giotto object available features (e.g. atac, rna, ...)
#' @param spatial_locs data.table or data.frame with coordinates for cell
#' centroids
#' @param sampleNames A character vector containing the ArchR project sample
#' name
#' @param ... additional arguments passed to `createGiottoObject`
#'
#' @returns A Giotto object with at least an atac or epigenetic modality
#'
#' @export
createGiottoObjectfromArchR <- function(archRproj,
    expression = NULL,
    expression_feat = "atac",
    spatial_locs = NULL,
    sampleNames = "sample1",
    ...) {
    # extract GeneScoreMatrix
    GeneScoreMatrix_summarizedExperiment <- ArchR::getMatrixFromProject(
        archRproj
    )
    GeneScoreMatrix <- slot(
        slot(
            GeneScoreMatrix_summarizedExperiment, "assays"
        ),
        "data"
    )[["GeneScoreMatrix"]]

    ## get cell names
    cell_names <- colnames(GeneScoreMatrix)
    cell_names <- gsub(paste0(sampleNames, "#"), "", cell_names)
    cell_names <- gsub("-1", "", cell_names)

    ## get gene names
    gene_names <- slot(
        GeneScoreMatrix_summarizedExperiment,
        "elementMetadata"
    )[["name"]]

    ## replace colnames with cell names
    colnames(GeneScoreMatrix) <- cell_names

    ## replace rownames with gene names
    rownames(GeneScoreMatrix) <- gene_names


    if (!is.null(expression)) {
        expression_matrix <- data.table::fread(expression)

        expression_cell_names <- colnames(expression_matrix)
        cell_names <- intersect(cell_names, expression_cell_names)

        expression_matrix <- Matrix::Matrix(as.matrix(expression_matrix[, -1]),
            dimnames = list(
                expression_matrix[[1]],
                colnames(expression_matrix[, -1])
            ),
            sparse = TRUE
        )

        expression <- expression_matrix[, cell_names]

        GeneScoreMatrix <- GeneScoreMatrix[, cell_names]
    }


    ## filter spatial locations
    if (!is.null(spatial_locs)) {
        x <- read.csv(spatial_locs)
        x <- x[x$cell_ID %in% cell_names, ]
        spatial_locs <- x
    }

    # Creating GiottoObject
    message("Creating GiottoObject")

    if (!is.null(expression)) {
        gobject <- createGiottoObject(
            expression = list(
                GeneScoreMatrix = GeneScoreMatrix,
                raw = expression
            ),
            expression_feat = expression_feat,
            spatial_locs = spatial_locs,
            ...
        )
    } else {
        gobject <- createGiottoObject(
            expression = list(GeneScoreMatrix = GeneScoreMatrix),
            expression_feat = expression_feat,
            spatial_locs = spatial_locs,
            ...
        )
    }

    # add LSI dimension reduction
    coordinates <- slot(archRproj, "reducedDims")[["IterativeLSI"]][["matSVD"]]

    ## clean cell names
    lsi_cell_names <- rownames(coordinates)
    lsi_cell_names <- gsub(paste0(sampleNames, "#"), "", lsi_cell_names)
    lsi_cell_names <- gsub("-1", "", lsi_cell_names)

    rownames(coordinates) <- lsi_cell_names

    coordinates <- coordinates[cell_names, ]

    dimension_reduction <- Giotto::createDimObj(
        coordinates = coordinates,
        name = "lsi",
        spat_unit = "cell",
        feat_type = expression_feat[1],
        method = "lsi"
    )
    gobject <- setDimReduction(gobject,
        dimension_reduction,
        spat_unit = "cell",
        feat_type = expression_feat[1],
        name = "lsi",
        reduction_method = "lsi"
    )

    return(gobject)
}

## CLASS ####
# ------- ###


setClass(
    "VisiumHDReader",
    slots = list(
        visiumHD_dir = "character",
        expression_source = "character",
        gene_column_index = "numeric",
        barcodes = "character",
        array_subset_row = "numeric",
        array_subset_col = "numeric",
        pxl_subset_row = "numeric",
        pxl_subset_col = "numeric",
        calls = "list"
    ),
    prototype = list(
        expression_source = "raw",
        gene_column_index = 2,
        barcodes = NULL,
        array_subset_row = NULL,
        array_subset_col = NULL,
        pxl_subset_row = NULL,
        pxl_subset_col = NULL,
        calls = list()
    )
)



# * show ####
setMethod("show", signature("VisiumHDReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "VisiumHDReader"))
    print_slots <- c(
        "dir", "expression_source", "gene_column_index",
        "barcodes", "array_subset_row", "array_subset_col",
        "pxl_subset_row", "pxl_subset_col",
        "funs"
    )
    pre <- sprintf("%s :", format(print_slots))
    names(pre) <- print_slots

    # dir
    d <- object@visiumHD_dir
    if (length(d) > 0L) {
        nch <- nchar(d)
        d <- abbrev_path(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }

    # expression_source
    expression_source <- object@expression_source
    cat(pre["expression_source"], expression_source, "\n")

    # gene_column_index
    gene_column_index <- object@gene_column_index
    cat(pre["gene_column_index"], gene_column_index, "\n")

    # barcodes
    barcodes <- ifelse(!is.null(object@barcodes), "found", "none")
    cat(pre["barcodes"], barcodes, "\n")

    # array_subset_row
    array_subset_row <- ifelse(!is.null(object@array_subset_row),
        "found", "none"
    )
    cat(pre["array_subset_row"], array_subset_row, "\n")

    # array_subset_col
    array_subset_col <- ifelse(!is.null(object@array_subset_col),
        "found", "none"
    )
    cat(pre["array_subset_col"], array_subset_col, "\n")

    # pxl_subset_row
    pxl_subset_row <- ifelse(!is.null(object@pxl_subset_row), "found", "none")
    cat(pre["pxl_subset_row"], pxl_subset_row, "\n")

    # pxl_subset_col
    pxl_subset_col <- ifelse(!is.null(object@pxl_subset_col), "found", "none")
    cat(pre["pxl_subset_col"], pxl_subset_col, "\n")

    # funs
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("VisiumHDReader"), function(x, ...) show(x))



#' @title Import a Visium HD assay
#' @name importVisiumHD
#' @description
#' Giotto import functionalities for Visium HD datasets. This function generates
#' a `VisiumHDReader` instance that has convenient reader functions for
#' converting individual pieces of Visium HD data into Giotto-compatible
#' representations when the param `visiumHD_dir` is provided.
#' A function that creates the full `giotto` object is also available.
#' These functions should have all param values provided as defaults, but
#' can be flexibly modified to do things such as look in alternative
#' directories or paths.
#' @param visiumHD_dir Visium HD output directory (e.g. square_016um)
#' @param expression_source character. Raw or filter expression data. Defaults
#' to "raw"
#' @param gene_column_index numeric. Expression column to use for gene names
#' 1 = Ensembl and 2 = gene symbols
#' @param barcodes character vector. (optional) Use if you only want to load
#' a subset of the pixel barcodes
#' @param array_subset_row numeric vector. (optional) Vector with min and max
#' values to subset based on array rows
#' @param array_subset_col numeric vector. (optional) Vector with min and max
#' values to subset based on array columns
#' @param pxl_subset_row numeric vector. (optional) Vector with min and max
#' values to subset based on row pixels
#' @param pxl_subset_col numeric vector. (optional) Vector with min and max
#' values to subset based on column pixels
#' @details
#' Loading functions are generated after the `visiumHD_dir` is added.
#' @returns VisiumHDReader object
#' @examples
#' \dontrun{
#' # Create a `VisiumHDReader` object
#' reader <- importVisiumHD()
#'
#' # Set the visiumHD_dir
#' reader$visiumHD_dir <- "path to visium HD dir"
#' readerHD$visiumHD_dir <- visiumHD_dir
#'
#' # Load tissue positions or create cell metadata
#' tissue_pos <- readerHD$load_tissue_position()
#' metadata <- readerHD$load_metadata()
#'
#' # Load matrix or create expression object
#' matrix <- readerHD$load_matrix()
#' expression_obj <- readerHD$load_expression()
#'
#' # Load transcript data (cell metadata, expression object, and transcripts
#' # per pixel)
#' my_transcripts <- readerHD$load_transcripts(
#'     array_subset_row = c(500, 1000),
#'     array_subset_col = c(500, 1000)
#' )
#'
#' # Create a `giotto` object and add the loaded data
#' # TODO
#' }
#' @export
importVisiumHD <- function(visiumHD_dir = NULL,
    expression_source = "raw",
    gene_column_index = 2,
    barcodes = NULL,
    array_subset_row = NULL,
    array_subset_col = NULL,
    pxl_subset_row = NULL,
    pxl_subset_col = NULL,
    shape = "hexagon",
    shape_size = 400) {
    # get params
    a <- list(Class = "VisiumHDReader")

    if (!is.null(visiumHD_dir)) {
        a$visiumHD_dir <- visiumHD_dir
    }

    a$expression_source <- expression_source
    a$gene_column_index <- gene_column_index

    if (!is.null(barcodes)) {
        a$barcodes <- barcodes
    }

    if (!is.null(array_subset_row)) {
        a$array_subset_row <- array_subset_row
    }

    if (!is.null(array_subset_col)) {
        a$array_subset_col <- array_subset_col
    }

    if (!is.null(pxl_subset_row)) {
        a$pxl_subset_row <- pxl_subset_row
    }

    if (!is.null(pxl_subset_col)) {
        a$pxl_subset_col <- pxl_subset_col
    }

    do.call(new, args = a)
}


# * init ####
setMethod("initialize", signature("VisiumHDReader"), function(.Object, visiumHD_dir,
    expression_source,
    gene_column_index,
    barcodes,
    array_subset_row,
    array_subset_col,
    pxl_subset_row,
    pxl_subset_col) {
    # provided params (if any)
    if (!missing(visiumHD_dir)) {
        checkmate::assert_directory_exists(visiumHD_dir)
        .Object@visiumHD_dir <- visiumHD_dir
    }

    if (!missing(expression_source)) {
        .Object@expression_source <- expression_source
    }

    if (!missing(gene_column_index)) {
        .Object@gene_column_index <- gene_column_index
    }

    if (!missing(barcodes)) {
        .Object@barcodes <- barcodes
    }

    if (!missing(array_subset_row)) {
        .Object@array_subset_row <- array_subset_row
    }

    if (!missing(array_subset_col)) {
        .Object@array_subset_col <- array_subset_col
    }

    if (!missing(pxl_subset_row)) {
        .Object@pxl_subset_row <- pxl_subset_row
    }

    if (!missing(pxl_subset_col)) {
        .Object@pxl_subset_col <- pxl_subset_col
    }

    # NULL case
    if (length(.Object@visiumHD_dir) == 0) {
        return(.Object) # return early if no path given
    }


    # detect paths and subdirs
    p <- .Object@visiumHD_dir


    .visiumHD_detect <- function(pattern, path = p, recursive = FALSE) {
        .detect_in_dir(
            pattern = pattern, path = path,
            recursive = recursive, platform = "visiumHD"
        )
    }


    filter_expr_dir <- .visiumHD_detect(
        pattern = "filtered_feature_bc_matrix",
        path = p
    )
    raw_expr_dir <- .visiumHD_detect(
        pattern = "raw_feature_bc_matrix",
        path = p
    )

    s <- .Object@expression_source
    if (s == "raw") {
        expr_dir <- raw_expr_dir
    } else if (s == "filter") {
        expr_dir <- filter_expr_dir
    } else {
        stop("expression source for visiumHD can only be raw or filter")
    }

    spatial_dir <- .visiumHD_detect(pattern = "spatial", path = p)


    c_index <- .Object@gene_column_index
    if (!c_index %in% c(1, 2)) {
        stop("gene column index can only be 1 (Ensembl) or 2 (gene symbols)")
    }

    read_folder_fun <- function(path = spatial_dir,
    gene_column_index = c_index,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    verbose = NULL) {
        .visiumHD_read_folder(
            path = path,
            expr_data = c("raw", "filter"),
            gene_column_index = 1,
            png_name = NULL,
            verbose = verbose
        )
    }
    .Object@calls$read_folder <- read_folder_fun

    ## matrix load call
    matrix_fun <- function(path = expr_dir,
    gene_column_index = c_index,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    verbose = NULL) {
        .visiumHD_matrix(
            path = path,
            gene_column_index = gene_column_index,
            remove_zero_rows = remove_zero_rows,
            split_by_type = split_by_type,
            verbose = verbose
        )
    }
    .Object@calls$load_matrix <- matrix_fun



    ## expression load call
    expression_fun <- function(path = expr_dir,
    gene_column_index = c_index,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    verbose = NULL) {
        .visiumHD_expression(
            path = path,
            gene_column_index = gene_column_index,
            remove_zero_rows = remove_zero_rows,
            split_by_type = split_by_type,
            verbose = verbose
        )
    }
    .Object@calls$load_expression <- expression_fun



    ## tissue position load call
    tissue_position_fun <- function(path = spatial_dir,
    verbose = NULL) {
        .visiumHD_tissue_positions(
            path = path,
            verbose = verbose
        )
    }
    .Object@calls$load_tissue_position <- tissue_position_fun

    ## scale factor load call
    read_scalefactors <- function(path = spatial_dir,
    verbose = NULL) {
        .visiumHD_read_scalefactors(
            path = path,
            verbose = verbose
        )
    }
    .Object@calls$load_scalefactor <- read_scalefactors

    load_image_fun <- function(path = spatial_dir,
    image_name = c("hires", "lowres"),
    scale_factor_name = c("tissue_hires_scalef", "tissue_lowres_scalef"),
    verbose = NULL) {
        .visiumHD_image(
            image_path = path,
            json_info = json_info,
            micron_scale = micron_scale,
            verbose = verbose
        )
    }
    .Object@calls$load_image <- load_image_fun

    load_poly_fun <- function(path = expr_dir,
    shape = "hexagon",
    shape_size = 400,
    name = "hex400",
    verbose = NULL) {
        .visiumHD_poly(
            path = path,
            shape = shape,
            shape_size = shape_size,
            name = name,
            verbose = verbose
        )
    }
    .Object@calls$load_polygon <- load_poly_fun

    ## metadata load call
    meta_fun <- function(path = spatial_dir,
    verbose = NULL) {
        .visiumHD_meta(
            path = path,
            verbose = verbose
        )
    }
    .Object@calls$load_metadata <- meta_fun



    ## transcript load call
    transcript_fun <- function(expr_path = expr_dir,
    tissue_positions_path = spatial_dir,
    barcodes = .Object@barcodes,
    array_subset_row = .Object@array_subset_row,
    array_subset_col = .Object@array_subset_col,
    pxl_subset_row = .Object@pxl_subset_row,
    pxl_subset_col = .Object@pxl_subset_col) {
        .visiumHD_transcript(
            expr_path = expr_path,
            tissue_positions_path = tissue_positions_path,
            barcodes = barcodes,
            array_subset_row = array_subset_row,
            array_subset_col = array_subset_col,
            pxl_subset_row = pxl_subset_row,
            pxl_subset_col = pxl_subset_col,
            verbose = TRUE
        )
    }
    .Object@calls$load_transcripts <- transcript_fun

    giotto_object_fun <- function(visiumHD_dir = visiumHD_dir,
    expr_data = c("raw", "filter"),
    gene_column_index = 2,
    tissue_positions_path = spatial_dir,
    expression_path = expr_path,
    metadata_path = spatial_dir,
    shape = "hexagon",
    shape_size = 400,
    load_expression = TRUE,
    load_metadata = TRUE,
    instructions = NULL,
    png_name = NULL,
    verbose = NULL) {
        load_expression <- as.logical(load_expression)
        load_metadata <- as.logical(load_metadata)

        funs <- .Object@calls

        # init gobject
        g <- giotto()
        if (!is.null(instructions)) {
            instructions(g) <- instructions
        }

        # transcripts
        tx_list <- funs$load_transcripts(
            expr_path = expr_dir,
            tissue_positions_path = spatial_dir,
            barcodes = .Object@barcodes,
            array_subset_row = .Object@array_subset_row,
            array_subset_col = .Object@array_subset_col,
            pxl_subset_row = .Object@pxl_subset_row,
            pxl_subset_col = .Object@pxl_subset_col
        )

        g <- setGiotto(g, tx_list$gpoints[["rna"]])

        polys <- funs$load_polygon(
            path = expr_dir,
            shape = shape,
            shape_size = shape_size,
            name = paste0(shape, shape_size),
            verbose = NULL
        )
        g <- setGiotto(g, polys)
        g <- addSpatialCentroidLocations(
            gobject = g,
            poly_info = paste0(shape, shape_size)
        )
        g <- calculateOverlap(g,
            spatial_info = paste0(shape, shape_size),
            feat_info = "rna"
        )
        g <- overlapToMatrix(g,
            poly_info = paste0(shape, shape_size),
            feat_info = "rna",
            name = "raw"
        )
        # images

        images <- funs$load_image(
            path = spatial_dir,
            image_name = c("hires", "lowres"),
            scale_factor_name = c("tissue_hires_scalef", "tissue_lowres_scalef"),
            verbose = NULL
        )

        g <- setGiotto(g, images)

        # expression & meta
        # Need to check that names agree for poly/expr/meta
        allowed_ids <- spatIDs(polys)

        if (load_expression) {
            exlist <- funs$load_expression(
                path = expr_dir,
                gene_column_index = c_index,
                remove_zero_rows = TRUE,
                split_by_type = TRUE
            )
            # have to run this for loop its one expression matrice rn

            # only keep allowed cells and set into gobject
            # for (ex in exlist) {
            #     bool <- colnames(ex[]) %in% allowed_ids
            #     ex[] <- ex[][, bool]
            #   g <- setGiotto(g, ex)
            # }
            g <- setGiotto(g, exlist[[1]])
        }

        if (load_metadata) {
            cx <- funs$load_metadata(
                path = metadata_path
            )
            # check this later causing to appear empty cellmetadata
            # cx[] <- cx[][cell_ID %in% allowed_ids,]
            g <- setGiotto(g, cx)
        }

        return(g)
    }
    .Object@calls$create_gobject <- giotto_object_fun

    return(.Object)
})


# * access ####

#' @export
setMethod("$", signature("VisiumHDReader"), function(x, name) {
    basic_info <- c(
        "visiumHD_dir", "expression_source", "gene_column_index",
        "barcodes", "array_subset_row", "array_subset_col",
        "pxl_subset_row", "pxl_subset_col"
    )
    if (name %in% basic_info) {
        return(methods::slot(x, name))
    }

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("VisiumHDReader"), function(x, name, value) {
    basic_info <- c(
        "visiumHD_dir", "expression_source", "gene_column_index",
        "barcodes", "array_subset_row", "array_subset_col",
        "pxl_subset_row", "pxl_subset_col"
    )
    if (name %in% basic_info) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }

    stop(sprintf(
        "Only items in '%s' can be set",
        paste0(basic_info, collapse = "', '")
    ))
})

#' @export
`.DollarNames.VisiumHDReader` <- function(x, pattern) {
    dn <- c(
        "visiumHD_dir", "expression_source", "gene_column_index",
        "barcodes", "array_subset_row", "array_subset_col",
        "pxl_subset_row", "pxl_subset_col"
    )
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}



.visiumHD_read_folder <- function(path,
    expr_data = c("raw", "filter"),
    gene_column_index = 1,
    png_name = NULL,
    verbose = NULL) {
    vmsg(.v = verbose, "A structured visium directory will be used")

    if (is.null(path)) {
        .gstop("path needs to be a path to a visium directory")
    }
    path <- path.expand(path)
    path <- dirname(path)
    if (!dir.exists(path)) .gstop(path, " does not exist!")
    expr_data <- match.arg(expr_data, choices = c("raw", "filter"))

    ## 1. check expression
    expr_counts_path <- switch(expr_data,
        "raw" = paste0(path, "/", "raw_feature_bc_matrix/"),
        "filter" = paste0(path, "/", "filtered_feature_bc_matrix/")
    )
    if (!file.exists(expr_counts_path)) {
        .gstop(expr_counts_path, "does not exist!")
    }

    ## 2. check spatial locations
    spatial_dir <- paste0(path, "/", "spatial")
    tissue_positions_path <- Sys.glob(paths = file.path(
        spatial_dir,
        "tissue_positions*"
    ))

    ## 3. check spatial image
    if (is.null(png_name)) {
        png_list <- list.files(spatial_dir, pattern = "*.png")
        png_name <- png_list[1]
    }
    png_path <- paste0(spatial_dir, "/", png_name)
    if (!file.exists(png_path)) .gstop(png_path, " does not exist!")

    ## 4. check scalefactors
    scalefactors_path <- paste0(spatial_dir, "/", "scalefactors_json.json")
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

.visiumHD_matrix <- function(path = path,
    gene_column_index = 2,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    verbose = TRUE) {
    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to matrix file provided or auto-detected"
        ), call. = FALSE)
    }

    # check existence and access rights of files
    checkmate::assert_directory_exists(path)

    vmsg(.v = verbose, "loading expression matrix ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # load expression results with the 10X default matrix function
    matrix_results <- get10Xmatrix(
        path_to_data = path,
        gene_column_index = gene_column_index,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )

    return(matrix_results)
}





.visiumHD_expression <- function(path,
    gene_column_index = 2,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    verbose = TRUE) {
    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to matrix file provided or auto-detected"
        ), call. = FALSE)
    }

    # check existence and access rights of files
    checkmate::assert_directory_exists(path)

    vmsg(.v = verbose, "loading expression matrix ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # load expression results with the 10X default matrix function
    matrix_results <- get10Xmatrix(
        path_to_data = path,
        gene_column_index = gene_column_index,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )


    exprObj <- createExprObj(
        expression_data = matrix_results,
        spat_unit = "pixel",
        feat_type = "rna",
        name = "raw",
        provenance = "pixel"
    )

    return(list("rna" = exprObj))
}

.visiumHD_tissue_positions <- function(path,
    verbose = TRUE) {
    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to tissue positions file provided or auto-detected"
        ), call. = FALSE)
    }

    # check existence and access rights of files
    checkmate::assert_directory_exists(path)

    vmsg(.v = verbose, "loading tissue positions file ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # check existence and access rights of files
    tissue_positions_path <- file.path(path, "tissue_positions.parquet")
    checkmate::assert_file_exists(tissue_positions_path)

    # read with parquet and data.table
    # tissue_positions = data.table::as.data.table(x = arrow::read_parquet(
    # tissue_positions_path))
    tissue_positions <- arrow::read_parquet(tissue_positions_path)

    return(tissue_positions)
}

.check_new_format <- function(json_scalefactors) {
    return(checkmate::test_list(
        x = json_scalefactors,
        types = "numeric",
        len = 5L
    ))
}

.adjust_expected_names <- function(new_format_2023, expected_json_names) {
    if (!new_format_2023) {
        expected_json_names <- expected_json_names[2:5]
    }
    return(expected_json_names)
}

.validate_json_names <- function(json_scalefactors, expected_names) {
    if (!setequal(names(json_scalefactors)[1:5], expected_names)) {
        warning(GiottoUtils::wrap_txt(
            "h5 scalefactors json names differ from expected.
       [Expected]:", expected_names, "\n",
            "[Actual]:", names(json_scalefactors)
        ))
    }
}

.visiumHD_read_scalefactors <- function(path, verbose = TRUE) {
    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to scale factors file provided or auto-detected"
        ), call. = FALSE)
    }
    # check existence and access rights of files
    checkmate::assert_directory_exists(path)

    vmsg(.v = verbose, "loading scale factors file ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # check existence and access rights of files
    json_path <- file.path(path, "scalefactors_json.json")
    checkmate::assert_file_exists(json_path)

    json_scalefactors <- read_json(json_path)

    expected_json_names <- c(
        "regist_target_img_scalef", # NEW as of 2023
        "spot_diameter_fullres",
        "tissue_hires_scalef",
        "fiducial_diameter_fullres",
        "tissue_lowres_scalef"
    )
    new_format_2023 <- .check_new_format(json_scalefactors)
    expected_json_names <- .adjust_expected_names(
        new_format_2023,
        expected_json_names
    )
    .validate_json_names(json_scalefactors, expected_json_names)

    return(json_scalefactors)
}

.visiumHD_micron_scale <- function(json_scalefactors) {
    # Check if json_scalefactors is a list and contains the required field
    if (!is.list(json_scalefactors)) {
        stop("json_scalefactors must be a list")
    }

    if (!"microns_per_pixel" %in% names(json_scalefactors)) {
        stop("microns_per_pixel field is missing from json_scalefactors")
    }

    # Extract the microns_per_pixel value
    px_to_micron <- json_scalefactors$microns_per_pixel

    return(px_to_micron)
}
.get_image_type <- function(png_name) {
    possible_types <- c("lowres", "hires")
    for (img_type in possible_types) {
        if (grepl(img_type, png_name)) {
            return(img_type)
        }
    }
    stop("image_path filename did not match either 'lowres' or 'hires'.
        Ensure the image is named accordingly.")
}

.get_scale_factor <- function(visiumHD_img_type, json_info) {
    if (is.null(json_info)) {
        warning("No scalefactors json info provided. VisiumHD image
                scale_factor defaulting to 1.")
        return(1)
    }

    checkmate::assert_list(json_info)

    scale_factor <- switch(visiumHD_img_type,
        "lowres" = json_info[["tissue_lowres_scalef"]],
        "hires" = json_info[["tissue_hires_scalef"]],
        stop("Unexpected image type: ", visiumHD_img_type)
    )

    if (is.null(scale_factor)) {
        stop(
            "Scale factor for ", visiumHD_img_type,
            " image not found in json_info."
        )
    }

    return(scale_factor)
}

.apply_micron_scale <- function(scale_factor, json_info, px_to_micron) {
    if (isTRUE(px_to_micron)) {
        scale_factor <- scale_factor * .visiumHD_micron_scale(json_info)
    }
    return(scale_factor)
}

.visiumHD_image <- function(image_path,
    json_info = NULL,
    micron_scale = FALSE,
    verbose = NULL) {
    # Assume image already checked
    vmsg(.v = verbose, .initial = " - ", "found image")

    # Determine image scalefactor to use
    if (missing(image_path)) {
        stop(wrap_txt(
            "No path to image file provided or auto-detected"
        ), call. = FALSE)
    }

    # 1. determine image scalefactor to use ---------------------------------- #
    image_path <- .visiumHD_read_folder(path = image_path)
    image_path <- image_path[[4]]
    json_info <- .visiumHD_read_scalefactors(dirname(image_path))
    if (!is.null(json_info)) checkmate::assert_list(json_info)
    png_name <- basename(image_path) # used for name pattern matching only

    if (is.null(json_info)) { # if none provided
        warning(wrap_txt(
            "No scalefactors json info provided.
            VisiumHD image scale_factor defaulting to 1"
        ))
        scale_factor <- 1
    } else { # if provided

        scale_factor <- NULL # initial value
    }
    visiumHD_img_type <- .get_image_type(png_name)
    scale_factor <- .get_scale_factor(visiumHD_img_type, json_info)
    # scale_factor <- .apply_micron_scale(scale_factor, json_info, px_to_micron)

    # 2. create image -------------------------------------------------------- #
    visiumHD_img <- createGiottoLargeImage(
        raster_object = image_path,
        name = "image",
        negative_y = TRUE,
        scale_factor = (1 / scale_factor)
    )

    visiumHD_img_list <- list(visiumHD_img)
    names(visiumHD_img_list) <- c("image")

    return(visiumHD_img_list)
}

.visiumHD_poly <- function(path,
    shape = "hexagon",
    shape_size = 400,
    name = "hex400",
    verbose = TRUE) {
    if (!shape %in% c("hexagon", "square", "circle")) {
        stop("Invalid shape. Please choose either
            'hexagon', 'square', or 'circle'.")
    }

    if (shape_size <= 0) {
        stop("Size must be a positive number.")
    }

    tp <- arrow::read_parquet(
        file = .visiumHD_read_folder(
            path,
            verbose = FALSE
        )[3]$tissue_positions_path,
        as_data_frame = FALSE
    )


    original_feat_ext <- tp %>%
        dplyr::summarise(
            xmin = min(pxl_row_in_fullres, na.rm = TRUE),
            xmax = max(pxl_row_in_fullres, na.rm = TRUE),
            ymin = min(pxl_col_in_fullres, na.rm = TRUE),
            ymax = max(pxl_col_in_fullres, na.rm = TRUE)
        ) %>%
        dplyr::collect() %>%
        {
            ext(.$xmin, .$xmax, .$ymin, .$ymax)
        }

    # ext(gpoints$rna@spatVector)
    message(paste("Creating a", shape, "polygon with size", shape_size))
    polygons <- tessellate(
        extent = original_feat_ext,
        shape = shape,
        shape_size = shape_size,
        name = name
    )

    return(polygons)
}

.visiumHD_meta <- function(path,
    verbose = TRUE) {
    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to tissue positions file provided or auto-detected"
        ), call. = FALSE)
    }

    # check existence and access rights of files
    checkmate::assert_directory_exists(path)

    vmsg(.v = verbose, "loading tissue positions file ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # check existence and access rights of files
    tissue_positions_path <- file.path(path, "tissue_positions.parquet")
    checkmate::assert_file_exists(tissue_positions_path)

    # read with parquet and data.table
    tissue_positions <- data.table::as.data.table(
        x = arrow::read_parquet(tissue_positions_path)
    )

    vmsg(.v = verbose, "creating metadata ...")

    data.table::setnames(tissue_positions, "barcode", "cell_ID")

    cx <- createCellMetaObj(
        metadata = tissue_positions,
        spat_unit = "pixel",
        feat_type = "rna",
        provenance = "pixel",
        verbose = verbose
    )
    return(cx)
}



.visiumHD_transcript <- function(expr_path,
    gene_column_index = 2,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    tissue_positions_path,
    barcodes = NULL,
    array_subset_row = c(500, 1000),
    array_subset_col = c(500, 1000),
    pxl_subset_row = NULL,
    pxl_subset_col = NULL,
    verbose = TRUE) {
    # function to create expression matrix
    matrix <- .visiumHD_matrix(
        path = expr_path,
        gene_column_index = gene_column_index,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type,
        verbose = verbose
    )


    # function to create tissue position data.table
    tissue_positions <- data.table::as.data.table(.visiumHD_tissue_positions(
        path = tissue_positions_path,
        verbose = verbose
    ))



    vmsg(
        .v = verbose,
        "creating visiumHD tissue position x expression data file ..."
    )

    # subset data
    if (!is.null(barcodes)) {
        vmsg(.v = verbose, "subsetting visiumHD on barcodes")
        tissue_positions <- tissue_positions[barcode %in% barcodes]
    }

    if (!is.null(array_subset_row)) {
        if (is.vector(array_subset_row) & length(array_subset_row) == 2) {
            vmsg(.v = verbose, "subsetting visiumHD on array rows")
            tissue_positions <- tissue_positions[
                array_row > array_subset_row[1] &
                    array_row < array_subset_row[2]
            ]
        } else {
            stop("array_subset_row was provided but is not a vector with
                length 2")
        }
    }

    if (!is.null(array_subset_col)) {
        if (is.vector(array_subset_col) & length(array_subset_col) == 2) {
            vmsg(.v = verbose, "subsetting visiumHD on array columns")
            tissue_positions <- tissue_positions[
                array_col > array_subset_col[1] &
                    array_col < array_subset_col[2]
            ]
        } else {
            stop("array_subset_col was provided but is not a vector with
                length 2")
        }
    }

    if (!is.null(pxl_subset_row)) {
        if (is.vector(pxl_subset_row) & length(pxl_subset_row) == 2) {
            vmsg(.v = verbose, "subsetting visiumHD on row pixels")
            tissue_positions <- tissue_positions[
                pxl_row_in_fullres > pxl_subset_row[1] &
                    pxl_row_in_fullres < pxl_subset_row[2]
            ]
        } else {
            cat("pxl_subset_row is ", pxl_subset_row)
            stop("pxl_subset_row was provided but is not a vector with
                length 2")
        }
    }

    if (!is.null(pxl_subset_col)) {
        if (is.vector(pxl_subset_col) & length(pxl_subset_col) == 2) {
            vmsg(.v = verbose, "subsetting visiumHD on column pixels")
            tissue_positions <- tissue_positions[
                pxl_col_in_fullres > pxl_subset_col[1] &
                    pxl_col_in_fullres < pxl_subset_col[2]
            ]
        } else {
            cat(pxl_subset_col)
            stop("pxl_subset_col was provided but is not a vector with
                length 2")
        }
    }

    # also subset matrix if needed
    if (any(!is.null(c(
        barcodes,
        array_subset_row, array_subset_col,
        pxl_subset_row, pxl_subset_col
    )))) {
        vmsg(.v = verbose, "subsetting visiumHD on expression matrix")
        matrix <- matrix[, colnames(matrix) %in% tissue_positions$barcode]
    }

    # convert expression matrix to minimal data.table object
    matrix_tile_dt <- data.table::as.data.table(Matrix::summary(matrix))
    genes <- matrix@Dimnames[[1]]
    samples <- matrix@Dimnames[[2]]
    matrix_tile_dt[, gene := genes[i]]
    matrix_tile_dt[, pixel := samples[j]]


    # merge data.table matrix and spatial coordinates to create input for
    # Giotto Polygons
    gpoints <- data.table::merge.data.table(
        matrix_tile_dt, tissue_positions,
        by.x = "pixel", by.y = "barcode"
    )
    gpoints <- gpoints[
        , .(pixel, pxl_row_in_fullres, pxl_col_in_fullres, gene, x)
    ]
    colnames(gpoints) <- c("pixel", "x", "y", "gene", "counts")

    gpoints <- createGiottoPoints(x = gpoints[, .(x, y, gene, pixel, counts)])

    # ensure output is always a list
    if (!is.list(gpoints)) {
        gpoints <- list(gpoints)
        names(gpoints) <- objName(gpoints[[1L]])
    }

    return(list(
        "matrix" = matrix,
        "tissue_positions" = tissue_positions,
        "gpoints" = gpoints
    ))
}

#' @title Create 10x VisiumHD Giotto Object
#' @name createGiottoVisiumHDObject
#' @description Given the path to a VisiumHD output folder, creates a
#' Giotto object
#' @param VisiumHD_dir full path to the exported visiumHD directory
#' @param expr_data raw or filtered data (see details)
#' @param gene_column_index which column index to select (see details)
#' @param expression_matrix_class class of expression matrix to use
#' (e.g. "dgCMatrix", "DelayedArray")
#' @param cores how many cores or threads to use to read data if paths are
#' provided
#' @param verbose be verbose
#' @inheritParams get10Xmatrix
#' @returns giotto object
#' @details
#' \itemize{
#'   \item{expr_data: raw will take expression data from raw_feature_bc_matrix
#'   and filter from filtered_feature_bc_matrix}
#'   \item{gene_column_index: which gene identifiers (names) to use if there
#'   are multiple columns (e.g. ensemble and gene symbol)}
#' }
#' @export
createGiottoVisiumHDObject <- function(visiumHD_dir = NULL,
    expr_data = c("raw", "filter"),
    gene_column_index = 1,
    instructions = NULL,
    expression_matrix_class = c("dgCMatrix", "DelayedArray"),
    cores = NA,
    verbose = FALSE) {
    # NSE vars
    barcode <- row_pxl <- col_pxl <- in_tissue <- array_row <- array_col <- NULL

    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    readerHD <- importVisiumHD()
    readerHD$visiumHD_dir <- visiumHD_dir

    argslist <- readerHD$read_folder()
    argslist$verbose <- verbose
    argslist$expression_matrix_class <- expression_matrix_class
    argslist$instructions <- instructions

    giotto_object <- do.call(.visiumHD_create, c(argslist, readerHD = readerHD))

    return(giotto_object)
}

.visiumHD_create <- function(expr_counts_path,
    gene_column_index = NULL, # folder
    tissue_positions_path,
    image_path = NULL,
    scale_json_path = NULL,
    png_name = NULL,
    instructions = NULL,
    expression_matrix_class = c("dgCMatrix", "DelayedArray"),
    readerHD = readerHD,
    verbose = NULL) {
    # NSE vars
    barcode <- cell_ID <- row_pxl <- col_pxl <- in_tissue <- array_row <-
        array_col <- NULL

    if (is.null(readerHD)) {
        stop("readerHD is not provided")
    }

    expr_counts_path <- readerHD$read_folder()[[1]]
    # 1. expression
    expr_results <- get10Xmatrix(
        path_to_data = expr_counts_path,
        gene_column_index = gene_column_index
    )

    # if expr_results is not a list, make it a list compatible with downstream
    if (!is.list(expr_results)) {
        expr_results <- list("Gene Expression" = expr_results)
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
    tissue_positions_path <- readerHD$read_folder()[[2]]
    spatial_results <- readerHD$load_tissue_position()
    data.table::setnames(spatial_results, old = "barcode", new = "cell_ID")
    spatial_locs <- spatial_results[
        , .(cell_ID, pxl_row_in_fullres, -pxl_col_in_fullres)
    ] # flip x and y
    colnames(spatial_locs) <- c("cell_ID", "sdimx", "sdimy")

    # 3. scalefactors (optional)
    json_info <- readerHD$load_scalefactor()

    # 4. image (optional)
    visium_png_list <- readerHD$load_image()

    # 5. metadata
    meta_results <- spatial_results[, .(cell_ID, in_tissue, array_row, array_col)]
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
    visium_polygons <- readerHD$load_polygon()
    giotto_object <- setPolygonInfo(
        gobject = giotto_object,
        x = visium_polygons,
        centroids_to_spatlocs = FALSE,
        verbose = FALSE,
        initialize = TRUE
    )

    return(giotto_object)
}

## CLASS ####
# ------- ###

setClassUnion("nullOrCharacter", c("NULL", "character"))
setClassUnion("nullOrNumeric", c("NULL", "numeric"))


setClass(
    "VisiumHDReader",
    slots = list(
        visiumhd_dir = "character",
        bin = "integer",
        micron = "logical",
        outdir = "nullOrCharacter",
        expression_source = "character",
        feature_id_type = "character", # gene column index
        tissue_only = "logical",
        barcodes = "nullOrCharacter",
        array_subset_row = "nullOrNumeric",
        array_subset_col = "nullOrNumeric",
        pxl_subset_row = "nullOrNumeric",
        pxl_subset_col = "nullOrNumeric",
        filter = "ANY",
        filter_coverage = "numeric",
        calls = "list"
    ),
    prototype = list(
        bin = 8L,
        expression_source = "raw",
        feature_id_type = "symbols",
        tissue_only = FALSE,
        barcodes = NULL,
        array_subset_row = NULL,
        array_subset_col = NULL,
        pxl_subset_row = NULL,
        pxl_subset_col = NULL,
        filter = NULL,
        filter_coverage = 0.5,
        calls = list()
    )
)



# * show ####
setMethod("show", signature("VisiumHDReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", class(object)))
    print_slots <- c(
        "dir", "bin", "micron", "outdir", "expression_source", 
        "feature_id_type", "tissue_only", "barcodes", 
        "array_subset_row", "array_subset_col",
        "pxl_subset_row", "pxl_subset_col", "filter_poly", 
        "filter_coverage_cutoff", "funs"
    )
    pre <- sprintf("%s :", format(print_slots))
    names(pre) <- print_slots

    # dir
    d <- object@visiumhd_dir
    if (length(d) > 0L) {
        d <- GiottoUtils::str_abbreviate(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }
    
    is_seg <- .visiumhd_is_segmentation_output(object@visiumhd_dir)
    
    # bin
    bin <- object@bin
    if (is_seg) {
        cat(pre["bin"], "cell segmentation", "\n")
    } else {
        cat(pre["bin"], bin, "\n")
    }
    
    
    # micron
    micron <- object@micron
    cat(pre["micron"], micron, "\n")
    
    # outdir
    if (!is_seg) {
        od <- object@outdir
        if (!is.null(od)) {
            od <- GiottoUtils::str_abbreviate(od)
            cat(pre["outdir"], od, "\n")
        } else {
            cat(pre["outdir"], "default\n")
        }
    }
    
    # expression_source
    expression_source <- object@expression_source
    cat(pre["expression_source"], expression_source, "\n")

    # feature_id_type
    feature_id_type <- object@feature_id_type
    cat(pre["feature_id_type"], feature_id_type, "\n")
    
    # tissue_only
    if (!is_seg) {
        tissue_only <- object@tissue_only
        cat(pre["tissue_only"], tissue_only, "\n")
    }

    # barcodes
    barcodes <- ifelse(!is.null(object@barcodes), "found", "none")
    cat(pre["barcodes"], barcodes, "\n")

    if (!is_seg ) {
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
        
        # filter
        filter <- ifelse(!is.null(object@filter), "found", "none")
        cat(pre["filter_poly"], filter, "\n")
        
        # filter_coverage
        filter_coverage <- object@filter_coverage
        cat(pre["filter_coverage_cutoff"], filter_coverage, "\n")
    }

    # funs
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("VisiumHDReader"), function(x, ...) show(x))



#' @title Import a Visium HD assay
#' @name importVisiumHD
#' @description Giotto import functionalities for Visium HD datasets. This
#' function generates a `VisiumHDReader` instance that has convenient reader
#' functions for converting individual pieces of Visium HD data into
#' Giotto-compatible representations when the param `visiumhd_dir` is provided.
#' A function that creates the full `giotto` object is also available. These
#' functions should have all param values provided as defaults, but can be
#' flexibly modified to do things such as look in alternative directories or
#' paths.
#' @param visiumhd_dir Visium HD output directory containing tar outputs or
#' their extracted `square_?um` directories.
#' @param bin numeric. Bin level to load from. Currently one of 2, 8, or 16.
#' @param micron logical. Set `TRUE` to load in micron scale instead of fullres
#'   image mapping.
#' @param outdir (optional) directory to unpack bin output tar contents into.
#'   (Default is the same directory as the tarfile.)
#' @param expression_source character. Raw or filter expression data. Defaults
#'   to "raw"
#' @param feature_id_type character. One of `"symbol"` or `"ensembl"`.
#'   Determines which to use as the feature identifiers.
#' @param tissue_only logical. Whether to only load information tagged as
#'   `in_tissue` (in `tissue_positions.parquet`). This is ignored by transcript
#'   loading.
#' @param barcodes (optional) character. Specific pixel barcodes to keep.
#' @param array_subset_row (optional) numeric vector, length = 2. Min/max of
#'   array rows to keep
#' @param array_subset_col (optional) numeric vector, length = 2. Min/max of
#'   array cols to keep
#' @param pxl_subset_row (optional) numeric vector, length = 2. Min/max of
#'   fullres image mapped rows to keep. Note that values are inverted into the
#'   negatives.
#' @param pxl_subset_col (optional) numeric vector, length = 2. Min/max of
#'   fullres image mapped cols to keep.
#' @param filter (optional) a `SpatVector`, `sf`, or `giottoPolygon` to
#'   spatially filter the data by.
#' @param filter_coverage_cutoff numeric between 0 and 1. Minimal fraction of
#'   pixel coverage by `filter` in order to be selected.
#' @details Loading functions are generated after the `visiumhd_dir` is added.
#' @returns `VisiumHDReader` object
#' @examples
#' \dontrun{
#' # Create a `VisiumHDReader` object
#' reader <- importVisiumHD(bin = 16, tissue_only = TRUE)
#' 
#' # visualization
#' plot(reader)
#' plot(reader, image = TRUE)
#' reader$micron <- TRUE
#' plot(reader, image = TRUE)
#' reader$micron <- FALSE # reset
#'
#' # Set the visiumhd_dir
#' reader$visiumhd_dir <- "path to visium HD dir"
#' readerHD$visiumhd_dir <- visiumhd_dir
#'
#' # Load tissue positions or create cell metadata
#' tissue_pos <- readerHD$load_tissue_position()
#' metadata <- readerHD$load_cellmeta()
#'
#' # Load matrix or create expression object
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
#' g <- reader$create_gobject()
#' }
#' @export
importVisiumHD <- function(
    visiumhd_dir = NULL,
    bin = 8,
    micron = FALSE,
    outdir = NULL,
    expression_source = "raw",
    feature_id_type = c("symbol", "ensembl"),
    tissue_only = FALSE,
    barcodes = NULL,
    array_subset_row = NULL,
    array_subset_col = NULL,
    pxl_subset_row = NULL,
    pxl_subset_col = NULL,
    filter = NULL,
    filter_coverage_cutoff = 0.5) {
    # get params
    a <- list(Class = "VisiumHDReader")

    if (!is.null(visiumhd_dir)) {
        a$visiumhd_dir <- visiumhd_dir
    }

    a$bin <- as.integer(bin)
    a$micron <- as.logical(micron)
    if (!is.null(outdir)) {
        a$outdir <- outdir
    }
    a$expression_source <- match.arg(expression_source, c("raw", "filtered"))
    a$feature_id_type <- match.arg(feature_id_type, c("symbol", "ensembl"))
    a$tissue_only <- as.logical(tissue_only)

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
    
    if (!is.null(filter)) {
        a$filter <- filter
    }
    
    if (!is.null(filter_coverage_cutoff) && !is.na(filter_coverage_cutoff)) {
        a$filter_coverage <- filter_coverage_cutoff
    }

    do.call(new, args = a)
}

#' @rdname importVisiumHD
#' @param x `VisiumHDReader`
#' @param image logical. Whether to plot the image
#' @param ... additional params to pass
#' @export
setMethod("plot", signature("VisiumHDReader", "missing"), 
    function(x, image = FALSE, ...) {
    checkmate::assert_logical(image)
    sl <- x$load_tissue_position(output = "spatLocsObj")
    if (image) {
        img <- x$load_image()
        plot(img, ...)
        plot(hull(sl), add = TRUE, border = "red")
    } else {
        plot(hull(sl), ...)
    }

})

# * init ####
setMethod("initialize", signature("VisiumHDReader"), function(.Object,
    visiumhd_dir,
    bin,
    micron,
    outdir,
    expression_source,
    feature_id_type,
    tissue_only,
    barcodes,
    array_subset_row,
    array_subset_col,
    pxl_subset_row,
    pxl_subset_col,
    filter,
    filter_coverage) {
    # provided params (if any) -------------------------------------------- #
    if (!missing(visiumhd_dir)) {
        checkmate::assert_directory_exists(visiumhd_dir)
        .Object@visiumhd_dir <- visiumhd_dir
    }
    
    if (!missing(bin)) {
        .Object@bin <- as.integer(bin)
    }
    
    if (!missing(micron)) {
        .Object@micron <- micron
    }
    
    if (!missing(outdir)) {
        .Object@outdir <- outdir
    }

    if (!missing(expression_source)) {
        .Object@expression_source <- expression_source
    }

    if (!missing(feature_id_type)) {
        .Object@feature_id_type <- feature_id_type
    }
    
    if (!missing(tissue_only)) {
        .Object@tissue_only <- tissue_only
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
    
    if (!missing(filter)) {
        .Object@filter <- filter
    }
    
    if (!missing(filter_coverage)) {
        .Object@filter_coverage <- filter_coverage
    }

    # NULL case
    if (length(.Object@visiumhd_dir) == 0) {
        return(.Object) # return early if no path given
    }

    # validate settings --------------------------------------------------- #
    ex_src <- .Object@expression_source 
    if (!ex_src %in% c("raw", "filtered")) {
        stop("[VisiumHD] `$expression_source` can only be \"raw\" or \"filtered\"\n",
             call. = FALSE)
    }
    fidtype <- .Object@feature_id_type
    if (!fidtype %in% c("symbol", "ensembl")) {
        stop("[VisiumHD] `$feature_id_type` can only be \"symbol\" or \"ensembl\"\n",
             call. = FALSE)
    }
    bin <- .Object@bin
    if (!bin %in% c(2, 8, 16)) {
        stop("[VisiumHD] `$bin` can only be one of 2, 8, or 16\n",
             call. = FALSE)
    }
    poly_filter <- .Object@filter
    if (!is.null(poly_filter)) {
        if (!inherits(poly_filter, c("SpatVector", "giottoPolygon", "sf"))) {
            stop("[VisiumHD] `$filter` must be giottoPolygon, sf, or SpatVector\n",
                 call. = FALSE)
        }
    }

    # detect paths and subdirs -------------------------------------------- #
    p <- .Object@visiumhd_dir
    
    # test if provided directory has expected output structure
    is_output <- .visiumhd_validate_output_dir(p,
        bin = bin, expression_source = ex_src, warn = FALSE
    )
    if (is_output) {
        binpath <- p # directly use if output dir structure match
    } else {
        # otherwise, check expected .tar or subdir naming
        binpath <- .visiumhd_detect_output_dir(p, 
            bin = bin, expression_source = ex_src
        )
    }
    
    # if it is a segmentation output directory
    is_seg <- .visiumhd_is_segmentation_output(binpath)
    
    if (!is_seg) {
        binpath2 <- .visiumhd_detect_output_dir(p, 
            bin = 2L, expression_source = ex_src
        )
    }

    # setup closures ------------------------------------------------------ #
    ## expression load call
    expression_fun <- function(
        path = binpath,
        feature_id_type = .Object@feature_id_type,
        remove_zero_rows = TRUE,
        split_by_type = TRUE,
        outdir = .Object@outdir,
        bin = .Object@bin,
        expression_source = .Object@expression_source,
        barcodes = .Object@barcodes,
        force_untar = FALSE,
        untar_params = list(),
        verbose = NULL) {
        a <- get_args_list()
        do.call(.visiumhd_expression, a)
    }
    .Object@calls$load_expression <- expression_fun

    ## tissue position load call
    if (!is_seg) {
        tissue_position_fun <- function(
        path = binpath,
        outdir = .Object@outdir,
        bin = .Object@bin,
        micron = .Object@micron,
        scalefactors_path = binpath,
        barcodes = .Object@barcodes,
        tissue_only = .Object@tissue_only,
        array_subset_row = .Object@array_subset_row,
        array_subset_col = .Object@array_subset_col,
        pxl_subset_row = .Object@pxl_subset_row,
        pxl_subset_col = .Object@pxl_subset_col,
        filter = .Object@filter,
        filter_coverage_cutoff = .Object@filter_coverage,
        force_untar = FALSE,
        output = c("spatLocsObj", "arrow", "data.frame", "barcodes"),
        verbose = NULL,
        ...) {
            a <- get_args_list(...)
            do.call(.visiumhd_tissue_positions, a)
        }
        .Object@calls$load_tissue_position <- tissue_position_fun
    }

    ## scale factor load call
    read_scalefactors <- function(
        path = binpath,
        outdir = .Object@outdir,
        bin = .Object@bin,
        force_untar = FALSE,
        verbose = NULL,
        ...) {
        a <- get_args_list(...)
        do.call(.visiumhd_scalefactors, a)
    }
    .Object@calls$load_scalefactors <- read_scalefactors

    ## image load call
    load_image_fun <- function(
        path = binpath,
        outdir = .Object@outdir,
        bin = .Object@bin,
        image_type = NULL,
        micron = .Object@micron,
        scalefactors_path = binpath,
        force_untar = FALSE,
        untar_params = list(),
        verbose = NULL) {
        a <- get_args_list()
        do.call(.visiumhd_image, a)
    }
    .Object@calls$load_image <- load_image_fun

    ## tessellated poly create call
    if (!is_seg) {
        tess_poly_fun <- function(
        tissue_positions_path = binpath,
        shape = "hexagon",
        shape_size = 400,
        name = sprintf("%s%d", shape, as.integer(shape_size)),
        bin = .Object@bin,
        micron = .Object@micron,
        scalefactors_path = binpath,
        outdir = .Object@outdir,
        force_untar = FALSE,
        untar_params = list(),
        verbose = NULL) {
            a <- get_args_list()
            do.call(.visiumhd_tessellate, a)
        }
        .Object@calls$tessellate_polygon <- tess_poly_fun
    }
    
    if (is_seg) {
        poly_fun <- function(
        path = binpath,
        type = c("cell", "nucleus"),
        name = NULL,
        graphclust_annotated = FALSE,
        scalefactors_path = NULL,
        micron = FALSE,
        id_fmt = "cellid_%09d-1",
        flip_vertical = TRUE,
        verbose = NULL
        ) {
            a <- get_args_list()
            do.call(.visiumhd_poly, a)
        }
        .Object@calls$load_polygon <- poly_fun 
    }

    ## metadata load call
    if (!is_seg) {
        meta_fun <- function(
        tissue_positions_path = binpath,
        outdir = .Object@outdir,
        bin = .Object@bin,
        barcodes = .Object@barcodes,
        tissue_only = .Object@tissue_only,
        force_untar = FALSE,
        untar_params = list(),
        verbose = NULL) {
            a <- get_args_list()
            do.call(.visiumhd_cellmeta, a)
        }
        .Object@calls$load_cellmeta <- meta_fun
    }

    ## transcript load call
    if (!is_seg) {
        transcript_fun <- function(
        expr_path = binpath2,
        feature_id_type = .Object@feature_id_type,
        remove_zero_rows = TRUE,
        split_by_type = TRUE,
        expression_source = .Object@expression_source,
        tissue_positions_path = binpath2,
        scalefactors_path = binpath2,
        micron = .Object@micron,
        outdir = .Object@outdir,
        bin = 2,
        tissue_only = .Object@tissue_only,
        barcodes = .Object@barcodes,
        array_subset_row = .Object@array_subset_row,
        array_subset_col = .Object@array_subset_col,
        pxl_subset_row = .Object@pxl_subset_row,
        pxl_subset_col = .Object@pxl_subset_col,
        filter = .Object@filter,
        filter_coverage_cutoff = .Object@filter_coverage,
        force_untar = FALSE,
        untar_params = list(),
        output = c("giottoPoints", "full"),
        verbose = NULL,
        ...) {
            a <- get_args_list(...) # expr_list and spatlocs can also be passed
            do.call(.visiumhd_transcript, a)
        }
        .Object@calls$load_transcripts <- transcript_fun
    }

    if (!is_seg) {
        giotto_object_fun <- function(
        # what to add to gobject
            load_expression = TRUE,
            load_spatlocs = TRUE,
            load_metadata = TRUE,
            load_transcripts = FALSE,
            load_image = TRUE,
            create_tessellated_polys = FALSE,
            # loading modifiers & subsetting
            bin = .Object@bin,
            micron = .Object@micron,
            tissue_only = .Object@tissue_only,
            barcodes = .Object@barcodes,
            array_subset_row = .Object@array_subset_row,
            array_subset_col = .Object@array_subset_col,
            pxl_subset_row = .Object@pxl_subset_row,
            pxl_subset_col = .Object@pxl_subset_col,
            filter = .Object@filter,
            filter_coverage_cutoff = .Object@filter_coverage,
            # item-specific params
            # [expression]
            expression_source = .Object@expression_source,
            feature_id_type = .Object@feature_id_type,
            expression_remove_zero_rows = TRUE,
            expression_split_by_type = TRUE,
            # [image]
            image_type = NULL, # not actually used if full image_path is provided.
            # [tessellate]
            tessellate_shape = "hexagon",
            tessellate_shape_size = 400,
            tessellate_name = sprintf("%s%d", 
                                      tessellate_shape, as.integer(tessellate_shape_size)),
            # filepaths to data
            tissue_positions_path = binpath,
            scalefactors_path = binpath,
            expression_path = binpath,
            image_path = binpath,
            # untar and filesys params
            outdir = .Object@outdir,
            force_untar = FALSE, # should only apply once
            untar_params = list(), # additional params for `untar()`
            # general params
            instructions = NULL,
            verbose = NULL) {
            
            checkmate::assert_list(untar_params)
            load_expression <- as.logical(load_expression)
            load_spatlocs <- as.logical(load_spatlocs)
            load_metadata <- as.logical(load_metadata)
            load_transcripts <- as.logical(load_transcripts)
            load_image <- as.logical(load_image)
            create_tessellated_polys <- as.logical(create_tessellated_polys)
            
            # validate settings
            pre <- "[VisiumHD]"
            if (load_spatlocs && !load_expression) {
                stop(wrap_txtf("%s load_spatlocs = TRUE & load_expression = FALSE:
                spatlocs can only be loaded with expression info.", pre))
            }
            if (load_metadata && !load_expression) {
                stop(wrap_txtf("%s load_metadata = TRUE & load_expression = FALSE:
                metadata can only be loaded with expression info.", pre))
            }
            
            # setup param lists
            basic_params <- untar_params
            basic_params$bin <- bin
            basic_params$force <- FALSE
            basic_params$outdir <- outdir
            basic_params$verbose <- verbose
            
            spat_filter_params <- list(
                array_subset_row = array_subset_row,
                array_subset_col = array_subset_col,
                pxl_subset_row = pxl_subset_row,
                pxl_subset_col = pxl_subset_col,
                filter = filter,
                filter_coverage_cutoff = filter_coverage_cutoff
            )
            
            mat_params <- list(
                feature_id_type = feature_id_type,
                remove_zero_rows = expression_remove_zero_rows,
                split_by_type = expression_split_by_type,
                expression_source = expression_source
            )
            
            scalef_params <- list(
                micron = micron,
                scalefactors_path = scalefactors_path
            )
            
            # directly attempt untar if `force_untar` = TRUE
            if (force_untar) {
                force_params <- basic_params
                force_params$force <- TRUE
                force_params$spatial_only <- FALSE
                do.call(.visiumhd_untar_if_not, force_params)
            }
            
            # data loading
            funs <- .Object@calls
            
            sl_barcodes <- sl <- NULL
            if (load_spatlocs) {
                tp_params <- basic_params
                tp_params$path <- tissue_positions_path
                tp_params$barcodes <- barcodes
                tp_params$tissue_only <- tissue_only
                tp_params$output <- "spatLocsObj"
                sl <- do.call(.visiumhd_tissue_positions, 
                              args = c(tp_params, spat_filter_params, scalef_params))
                # update barcodes with any sl results
                sl_barcodes <- spatIDs(sl)
                barcodes <- sl_barcodes %null% barcodes
            }
            
            ex_barcodes <- expr_list <- NULL
            if (load_expression) {
                expr_params <- basic_params
                expr_params$path <- expression_path
                expr_params$barcodes <- barcodes
                expr_list <- do.call(.visiumhd_expression, 
                                     args = c(expr_params, mat_params)) # this is a list output
                # update barcodes with any expr results
                ex_barcodes <- spatIDs(expr_list[[1L]])
                barcodes <- ex_barcodes %null% barcodes
                # match expression barcodes back onto spatlocs just in case
                if (!is.null(sl)) sl <- sl[barcodes]
            }
            
            tx_list <- NULL
            if (load_transcripts) {
                tx_params <- basic_params
                tx_params$bin <- 2L # override basic_param
                if (bin == 2) {
                    # reuse if available and of bin 2
                    tx_params$expr_list <- expr_list 
                    tx_params$spatlocs <- sl
                }
                tx_params$output <- "full"
                # only used if `sl` is still NULL or not of bin 2
                tx_params$expr_path <- binpath2
                # only used if `expr_list` is still NULL or not of bin 2
                tx_params$tissue_positions_path <- binpath2
                tx_params$tissue_only <- FALSE
                tx_params$barcodes <- NULL # ignore
                # custom spat filter params
                tx_params$pxl_subset_row <- pxl_subset_row
                tx_params$pxl_subset_col = pxl_subset_col
                tx_params$filter = filter
                tx_params$filter_coverage_cutoff = filter_coverage_cutoff
                if (!is.null(array_subset_row)) {
                    tx_params$array_subset_row = array_subset_row * (bin / 2)
                }
                if (!is.null(array_subset_col)) {
                    tx_params$array_subset_col = array_subset_col * (bin / 2)
                }
                res <- do.call(funs$load_transcripts, 
                               args = c(tx_params, scalef_params, mat_params))
                tx_list <- res$gpoints
                # extract and update barcodes
                barcodes <- spatIDs(res$tissue_positions)
                rm(res) # cleanup
            }
            
            cmeta <- NULL
            if (load_metadata) {
                cm_params <- basic_params
                cm_params$tissue_positions_path <- tissue_positions_path
                cm_params$barcodes <- barcodes
                cm_params$tissue_only <- tissue_only
                cmeta <- do.call(.visiumhd_cellmeta, cm_params)
            }
            
            gimg <- NULL
            if (load_image) {
                img_params <- basic_params
                img_params$path <- image_path
                img_params$image_type <- image_type
                gimg <- do.call(.visiumhd_image, 
                                c(img_params, scalef_params))
            }
            
            tess_poly <- NULL
            if (create_tessellated_polys) {
                tess_params <- basic_params
                tess_params$tissue_positions_path <- tissue_positions_path
                tess_params$shape <- tessellate_shape
                tess_params$shape_size <- tessellate_shape_size
                tess_params$name <- tessellate_name
                tess_poly <- do.call(.visiumhd_tessellate, 
                                     args = c(tess_params, scalef_params))
            }
            
            # assemble giotto
            vmsg(.v = verbose, "assembling object ...")
            g <- giotto(instructions = instructions) # init
            if (!is.null(expr_list)) {
                g <- setGiotto(g, expr_list, verbose = verbose)
            }
            if (!is.null(sl)) {
                g <- setGiotto(g, sl, verbose = verbose)
            }
            if (!is.null(cmeta)) { # cmeta added here since there is a check in setCellMetadata
                g <- setGiotto(g, cmeta, verbose = verbose)
            }
            if (!is.null(gimg)) {
                g <- setGiotto(g, gimg, verbose = verbose)
            }
            if (!is.null(tess_poly)) {
                g <- setGiotto(g, tess_poly, verbose = verbose)
            }
            if (!is.null(tx_list)) {
                g <- setGiotto(g, tx_list, verbose = verbose)
            }
            gc(verbose = FALSE)
            g
        }
    } else {
        giotto_object_fun <- function( # what to add to gobject
            load_expression = TRUE,
            load_polygons = c("cell", "nucleus"),
            graphclust_annotated = FALSE,
            load_image = TRUE,
            # loading modifiers & subsetting
            micron = .Object@micron,
            barcodes = .Object@barcodes,
            # item-specific params
            # [expression]
            expression_source = .Object@expression_source,
            feature_id_type = .Object@feature_id_type,
            expression_remove_zero_rows = TRUE,
            expression_split_by_type = TRUE,
            # [image]
            image_type = NULL, # not actually used if full image_path is provided.
            # filepaths to data
            scalefactors_path = binpath,
            expression_path = binpath,
            image_path = binpath,
            geojson_path = binpath,
            # general params
            instructions = NULL,
            verbose = NULL,
            ...) {
            
            not_used <- names(list(...))
            if (length(not_used) > 0L) {
                vmsg(.v = verbose, 
                     "[visiumHD] params:", toString(not_used), 
                     "not used with segmentation outputs")
            }
            
            load_expression <- as.logical(load_expression)
            load_image <- as.logical(load_image)
            
            # setup param lists
            basic_params <- list()
            basic_params$bin <- NA_integer_ # dummy value
            basic_params$force <- FALSE # not needed for v4 outputs
            basic_params$outdir <- NULL # not needed for v4 outputs
            basic_params$verbose <- verbose
            
            mat_params <- list(
                feature_id_type = feature_id_type,
                remove_zero_rows = expression_remove_zero_rows,
                split_by_type = expression_split_by_type,
                expression_source = expression_source
            )
            
            scalef_params <- list(
                micron = micron,
                scalefactors_path = scalefactors_path
            )
            
            # data loading
            funs <- .Object@calls
            
            poly_list <- NULL
            if (length(load_polygons) > 0L) {
                load_polygons <- match.arg(load_polygons, 
                    choices = c("cell", "nucleus"), 
                    several.ok = TRUE
                )
                poly_params <- list(verbose = verbose)
                poly_params$path <- geojson_path
                poly_params$graphclust_annotated <- graphclust_annotated
                poly_params$scalefactors_path <- scalefactors_path
                poly_params$micron <- micron
                poly_params$verbose <- verbose
                poly_params$barcodes <- barcodes
                poly_list <- lapply(load_polygons, function(ptype) {
                    do.call(.visiumhd_poly, c(list(type = ptype), poly_params))
                })
            }
            
            ex_barcodes <- expr_list <- NULL
            if (load_expression) {
                expr_params <- basic_params
                expr_params$path <- expression_path
                expr_params$barcodes <- barcodes
                expr_list <- do.call(.visiumhd_expression, 
                    args = c(expr_params, mat_params)) # this is a list output
            }
            
            gimg <- NULL
            if (load_image) {
                img_params <- basic_params
                img_params$path <- image_path
                img_params$image_type <- image_type
                gimg <- do.call(.visiumhd_image, 
                    c(img_params, scalef_params))
            }
            
            # assemble giotto
            vmsg(.v = verbose, "assembling object ...")
            g <- giotto(instructions = instructions) # init
            if (!is.null(poly_list)) {
                g <- setGiotto(g, poly_list, 
                    centroids_to_spatlocs = TRUE, verbose = verbose
                )
            }
            if (!is.null(expr_list)) {
                g <- setGiotto(g, expr_list, verbose = verbose)
            }
            if (!is.null(gimg)) {
                g <- setGiotto(g, gimg, verbose = verbose)
            }
            gc(verbose = FALSE)
            g
        }
    }
    .Object@calls$create_gobject <- giotto_object_fun

    return(.Object)
})


# * access ####

.visiumhd_editable_slots <- c(
    "visiumhd_dir", "bin", "micron", "outdir", "expression_source", 
    "feature_id_type", "tissue_only", "barcodes",
    "array_subset_row", "array_subset_col",
    "pxl_subset_row", "pxl_subset_col",
    "filter", "filter_coverage"
)

#' @export
setMethod("$", signature("VisiumHDReader"), function(x, name) {
    if (name %in% .visiumhd_editable_slots) {
        return(methods::slot(x, name))
    }

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("VisiumHDReader"), function(x, name, value) {
    if (!name %in% .visiumhd_editable_slots) {
        stop(sprintf(
            "Only items in '%s' can be set",
            toString(.visiumhd_editable_slots)
        ))
    }
    
    if (name == "bin") value <- as.integer(value)
    
    methods::slot(x, name) <- value
    return(initialize(x))
})

#' @export
.DollarNames.VisiumHDReader <- function(x, pattern) {
    dn <- .visiumhd_editable_slots
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}


# MODULAR ####

.visiumhd_is_segmentation_output <- function(path) {
    p <- path[[1L]]
    dir.exists(p) && grepl("segmented_output", basename(p))
}

# detect default expected locations for output directory
.visiumhd_detect_output_dir <- function(p, bin, 
    expression_source = c("raw", "filtered"),
    verbose = NULL) {
    p <- normalizePath(p)
    bin <- as.integer(bin)
    expression_source <- match.arg(expression_source, c("raw", "filtered"))
    
    # template detection within dir function
    .detect <- function(pattern, path = p, recursive = FALSE) {
        .detect_in_dir(
            pattern = pattern, path = path,
            recursive = recursive, platform = "visiumHD"
        )
    }
    
    binpath <- NULL
    # detect if input is already the extracted bin folder (based on naming)
    is_output <- grepl(pattern = sprintf(".*square_%03dum$", bin), p)
    if (is_output) binpath <- p
    
    # if not already extracted...
    if (is.null(binpath)) {
        # detect if tar file is present (preferred)
        binpath <- .detect(sprintf("%03dum_outputs.tar", bin))
        # if tar missing, try extracted default name as subdirectory
        if (is.null(binpath)) {
            binpath <- .detect(sprintf("square_%03dum$", bin))
            # validation of bin subdir contents (if missing, send warning)
            if (!is.null(binpath)) {
                # no returns. This is just for checking
                .visiumhd_validate_output_dir(binpath, bin, expression_source, 
                    warn = TRUE
                )
            }
        }
    }
    # if still missing, throw warning
    if (is.null(binpath)) {
        warning(wrap_txtf(
            "Expected path to %d bin .tar or extracted 'square_%03dum' subdirectory not discovered.
            `Use importVisiumHD() for more specific filepaths`", bin, bin),
            call. = FALSE
        )
    }
    binpath
}

# * validate output dir structure ####

# check if a provided directory path is likely a bin output directory
# contains:
# 
# <binned data>
# - a set of expression data (.h5 or matrix market)
# - spatial data
#   - tissue positions
#   - scalefactors
#   - hi/lowres images
# 
# <segmented data>
# - cell_segmentations (any)
# - nucleus_segmentations (any)
# - graphclust_annotated_cell_segmentations.geojson
# - a set of expression data (.h5 or matrix market)
# - spatial data
#   - scalefactors
#   - hi/lowres images

.visiumhd_validate_output_dir <- function(binpath, bin, expression_source, 
    warn = TRUE) {
    bin <- as.integer(bin)
    expression_source <- match.arg(expression_source, c("raw", "filtered"))
    p <- binpath # shorten naming
    
    if (is.null(p)) return(FALSE)
    
    # determine if segemented or bin outputs
    type <- ifelse(grepl("segmented_outputs", basename(p)), "seg", "bin")
    
    output_base <- switch(type,
        "bin" = sprintf("visiumHD bin %d output", bin),
        "seg" = "visiumHD segmented output"
    )
    
    .detect <- function(pattern, path = p, recursive = FALSE) {
        .detect_in_dir(
            pattern = pattern, path = path, recursive = recursive,
            platform = output_base, warn = warn
        )
    }

    # check that at least one of the matrix market or .h5 are present
    expr_pattern <- switch(expression_source,
        "raw" = "raw_feature",
        "filtered" = "filtered_feature"
    )
    expr <- .detect(expr_pattern)
    if (is.null(expr)) return(FALSE)
    spatial <- .detect("spatial")
    if (is.null(spatial)) return(FALSE)
    
    # proceed with spatial subdir checks if spatial found
    .detect_spat <- function(pattern, 
        path = spatial, recursive = FALSE) {
        .detect_in_dir(
            pattern = pattern, path = path, recursive = recursive,
            platform = paste(output_base, "spatial"), warn = warn
        )
    }
    
    # no return of these values. This is just for checking/warning throws
    # when these items are expected from a global structured output
    # 
    # actual detection is handled per modular .visiumhd_* load function
    if (type == "bin") .detect_spat("tissue_positions.parquet")
    # tissue positions is not found in segmented outputs
    contents <- c(
        .detect_spat("scalefactors_json.json"),
        .detect_spat("hires_image"),
        .detect_spat("lowres_image")
    )
    if (length(contents) == 0L) return(FALSE)
    TRUE
}

.visiumhd_expression <- function(path,
    feature_id_type = c("symbols", "ensembl"),
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    outdir = NULL,
    bin = 8,
    expression_source = NULL,
    barcodes = NULL,
    force_untar = FALSE,
    untar_params = list(),
    verbose = NULL,
    verbose2 = TRUE) {
    checkmate::assert_list(untar_params)
    checkmate::assert_character(outdir, null.ok = TRUE)
    checkmate::assert_integerish(bin)
    checkmate::assert_character(barcodes, null.ok = TRUE)
    checkmate::assert_character(expression_source, null.ok = TRUE)
    checkmate::assert_logical(split_by_type)
    checkmate::assert_logical(remove_zero_rows)
    checkmate::assert_logical(force_untar)
    feature_id_type <- match.arg(feature_id_type, c("symbols", "ensembl"))

    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to matrix file provided or auto-detected"
        ), call. = FALSE)
    }
    
    if (!file.exists(path)) { # check if file OR dir
        stop("[VisiumHD] filepath does not exist\n", call. = FALSE)
    }

    if (verbose2 && !isTRUE(verbose %in% c("debug", "debug_log"))) {
        vmsg(.v = verbose, "loading expression matrix ...")
    }
    vmsg(.v = verbose, .is_debug = TRUE, path)
    
    # determine expression source
    expression_path <- NULL
    if ((.is_tar(path) || checkmate::test_directory_exists(path)) &&
        !.visiumhd_is_mm_dir(path)) {
        # if tar or the tar extracted dir...
        # - `expression_source` param not needed in other cases
        if (is.null(expression_source)) expression_source <- "raw" # default
        
        agg_type <- ifelse(grepl("segmented", basename(path)), "seg", "bin")
        expression_source <- paste(expression_source, agg_type, sep = "_")
        
        expression_path <- switch(expression_source,
            "raw_bin" = "raw_feature_bc_matrix",
            "filtered_bin" = "filtered_feature_bc_matrix",
            "raw_seg" = "raw_feature_cell_matrix",
            "filtered_seg" = "filtered_feature_cell_matrix",
            stop("[VisiumHD] expression_source:", expression_source,
                 "unrecognized\n", call. = FALSE)
        )
        
        # point to extracted matrix market format
        untar_params$path <- path
        untar_params$outdir <- outdir
        untar_params$outdir_subpath <- expression_path
        untar_params$bin <- bin
        untar_params$force <- force_untar
        untar_params$verbose <- verbose
        untar_params$spatial_only <- FALSE
        # update path
        path <- do.call(.visiumhd_untar_if_not, untar_params)
    }

    expr_params <- list(
        path = path,
        feature_id_type = feature_id_type,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )
    if (.visiumhd_is_mm_dir(path)) { # if fullpath to matrix market dir
        m <- do.call(.visiumhd_expression_mm, expr_params)
    } else if (.visiumhd_is_h5(path)) {
        m <- do.call(.visiumhd_expression_h5, expr_params)
    } else {
        stop("[VisiumHD] unrecognized expression matrix format", call. = FALSE)
    }
    if (!is.list(m)) m <- list(rna = m)
    
    # `m` may be a list depending on split_by_type and the data contained.
    lapply(names(m), function(feat_type) {
        x <- m[[feat_type]]
        if (!is.null(barcodes)) {
            bool <- unname(colnames(x)) %in% barcodes
            x <- x[, bool, drop = FALSE]
        }
        
        if (agg_type == "bin") {
            su <- sprintf("bin%03d", as.integer(bin))
        } else {
            su <- "cell"
        }
        createExprObj(x,
            name = "raw",
            spat_unit = su,
            feat_type = feat_type,
            provenance = su
        )
    })
}

.visiumhd_expression_mm <- function(path, 
    feature_id_type = c("symbols", "ensembl"),
    remove_zero_rows = TRUE,
    split_by_type = TRUE) {
    if (feature_id_type == "symbols") gene_column_index <- 2
    else gene_column_index <- 1
    
    # load expression results with the 10X default matrix function
    get10Xmatrix(path,
        gene_column_index = gene_column_index,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )
}

.visiumhd_expression_h5 <- function(path,
    feature_id_type = c("symbols", "ensembl"),
    remove_zero_rows = TRUE,
    split_by_type = TRUE) {
    feature_id_type <- match.arg(feature_id_type, c("symbols", "ensembl"))
    
    get10Xmatrix_h5(path,
        gene_ids = feature_id_type,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )
}

.visiumhd_is_mm_dir <- function(path) {
    if (!file.exists(path)) { # works for files and dirs
        stop("[VisiumHD] filepath does not exist\n", call. = FALSE)
    }
    if (!checkmate::test_directory_exists(path)) return(FALSE)
    fnames <- list.files(path)
    required_files <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
    all(required_files %in% fnames)
}

.visiumhd_is_h5 <- function(path) {
    if (!file.exists(path)) {
        stop("[VisiumHD] filepath does not exist\n", call. = FALSE)
    }
    "h5" %in% GiottoUtils::file_extension(path)
}

# get tissue positions information and possibly perform several filters.
# spatial filtering all happens through this function
# all outputs have filtering and micron rescale (if needed) applied.
# arrow and data.frame both contain all types of info from tissue_positions.parquet in convenient format
# barcode just has barcodes
# spatLocsObj output contains only spatlocs info (barcodes and xy in px or micron)
.visiumhd_tissue_positions <- function(path,
    outdir = NULL,
    bin = 8,
    micron = FALSE,
    scalefactors_path = NULL,
    barcodes = NULL,
    tissue_only = FALSE,
    array_subset_row = NULL,
    array_subset_col = NULL,
    pxl_subset_row = NULL,
    pxl_subset_col = NULL,
    filter = NULL,
    filter_coverage_cutoff = 0.5,
    force_untar = FALSE,
    verbose = NULL,
    verbose2 = TRUE, # escape hatch when called by other functions
    output = c("spatLocsObj", "arrow", "data.frame", "barcodes"),
    ...) {
    package_check("dplyr")
    package_check("arrow", custom_msg = sprintf(
        "package 'arrow' is not yet installed\n\n To install:\n%s\n%s%s",
        "Sys.setenv(ARROW_WITH_ZSTD = \"ON\") ",
        "install.packages(\"arrow\", ",
        "repos = c(\"https://apache.r-universe.dev\"))"
    ))
    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to tissue positions file provided or auto-detected"
        ), call. = FALSE)
    }
    .assert_numeric_range <- function(x, what) {
        if (is.null(x)) return(invisible())
        checkmate::assert_numeric(x, len = 2L)
        if (x[[1L]] > x[[2L]]) {
            stop(wrap_txtf(
                "'%s' range element 2 must be greater than element 1\n", what
            ), call. = FALSE)
        }
    }
    
    checkmate::assert_logical(tissue_only)
    checkmate::assert_logical(micron)
    checkmate::assert_character(barcodes, null.ok = TRUE)
    .assert_numeric_range(array_subset_row, what = "array_subset_row")
    .assert_numeric_range(array_subset_col, what = "array_subset_col")
    .assert_numeric_range(pxl_subset_row, what = "pxl_subset_row")
    .assert_numeric_range(pxl_subset_col, what = "pxl_subset_col")
    checkmate::assert_logical(force_untar)
    output <- match.arg(output, c("spatLocsObj", "arrow", "data.frame", "barcodes"))

    if (verbose2 && !isTRUE(verbose %in% c("debug", "debug_log"))) {
        vmsg(.v = verbose, "loading tissue positions file ...")
    }
    vmsg(.v = verbose, .is_debug = TRUE, path)

    spatlocs_path <- .visiumhd_untar_if_not(path,
        outdir = outdir, 
        outdir_subpath = file.path("spatial", "tissue_positions.parquet"),
        bin = bin,
        force = force_untar,
        verbose = verbose,
        spatial_only = TRUE,
        ...
    )
    data <- arrow::open_dataset(spatlocs_path)
    
    # spatial processing steps ------------------------------------------ #
    # invert y values
    data <- dplyr::mutate(data, pxl_row_in_fullres = pxl_row_in_fullres * -1L)
    
    # get micron info if needed
    if (micron) {
        if (!checkmate::test_file_exists(scalefactors_path)) {
            stop("[VisiumHD] if micron = TRUE, scalefactors_path must be filepath to scalefactors_json.json or the bin output h5 or its exported data.")
        }
        json_info <- .visiumhd_scalefactors(scalefactors_path, 
            outdir = outdir, 
            bin = bin, 
            force_untar = force_untar, 
            verbose = verbose,
            verbose2 = FALSE,
            ...
        )
        px2um <- .visiumhd_micron_scale(json_info)
        # apply micron scaling
        data <- dplyr::mutate(data, pxl_row_in_fullres = pxl_row_in_fullres * px2um)
        data <- dplyr::mutate(data, pxl_col_in_fullres = pxl_col_in_fullres * px2um)
    }
    
    # Spatial filters --------------------------------------------------- #
    # [XY] filtering
    if (!is.null(pxl_subset_row)) {
        data <- dplyr::filter(data, 
            pxl_row_in_fullres >= pxl_subset_row[[1L]] & pxl_row_in_fullres <= pxl_subset_row[[2L]]
        )
    }
    if (!is.null(pxl_subset_col)) {
        data <- dplyr::filter(data, 
            pxl_col_in_fullres >= pxl_subset_col[[1L]] & pxl_col_in_fullres <= pxl_subset_col[[2L]])
    }
    
    # [POLY] filtering
    if (!is.null(filter)) { # filter by poly
        data <- .visiumhd_array_xy_filter(
            tissue_positions_arrow = data,
            bin = bin,
            filter = filter,
            coverage_cutoff = filter_coverage_cutoff,
            output = "arrow",
            verbose = verbose
        )
    }
    
    # Non spatial filters ----------------------------------------------- #
    # [ARRAY] filtering
    if (!is.null(array_subset_row)) {
        data <- dplyr::filter(data, 
            array_row >= array_subset_row[[1L]] & array_row <= array_subset_row[[2L]]
        )
    }
    if (!is.null(array_subset_col)) {
        data <- dplyr::filter(data, 
            array_col >= array_subset_col[[1L]] & array_col <= array_subset_col[[2L]]
        )
    }
    # [ID] filtering
    if (!is.null(barcodes)) {
        data <- dplyr::filter(data, barcode %in% barcodes)
    }
    # [IN TISSUE] attribute filtering
    if (tissue_only) {
        data <- dplyr::filter(data, in_tissue == 1L)
    }
    
    # col renaming
    lookup <- c(sdimx = "pxl_col_in_fullres", sdimy = "pxl_row_in_fullres")
    data <- dplyr::rename(data, dplyr::all_of(lookup))
    
    switch(output,
        "arrow" = return(data),
        "data.frame" = return(dplyr::collect(data)), 
        "barcodes" = return(dplyr::pull(data, "barcode", as_vector = TRUE))
    )
    
    spat_df <- dplyr::collect(dplyr::select(data, 
        c("cell_ID" = "barcode", "sdimx", "sdimy")
    ))
    
    # output spatLocsObj
    su <- sprintf("bin%03d", as.integer(bin))
    sl <- createSpatLocsObj(spat_df,
        name = "raw",
        spat_unit = su,
        provenance = su,
        verbose = FALSE
    )
    sl
}

.visiumhd_validate_json_names <- function(json_scalefactors, expected_names) {
    if (!all(expected_names %in% names(json_scalefactors))) {
        warning(GiottoUtils::wrap_txt(
            "scalefactors json names differ from expected. [Expected]:", 
            expected_names, "\n",
            "[Actual]:", names(json_scalefactors)
        ), call. = FALSE)
    }
}

.visiumhd_scalefactors <- function(path, 
    outdir = NULL, 
    bin = 8,
    force_untar = FALSE, 
    verbose = TRUE,
    verbose2 = TRUE,
    ...) {
    # check if path is provided
    if (missing(path)) {
        stop(wrap_txt(
            "No path to scale factors file provided or auto-detected"
        ), call. = FALSE)
    }

    if (verbose2 && !isTRUE(verbose %in% c("debug", "debug_log"))) {
        vmsg(.v = verbose, "loading scale factors file ...")
    }
    vmsg(.v = verbose, .is_debug = TRUE, path)
    
    json_path <- .visiumhd_untar_if_not(path,
        outdir = outdir, 
        outdir_subpath = file.path("spatial", "scalefactors_json.json"),
        bin = bin,
        force = force_untar,
        verbose = verbose,
        spatial_only = TRUE,
        ...
    )

    json_scalefactors <- read_json(json_path)

    expected_json_names <- c(
        # "spot_diameter_fullres",
        # "bin_size_um",
        "microns_per_pixel",
        "regist_target_img_scalef",
        "tissue_lowres_scalef",
        "fiducial_diameter_fullres",
        "tissue_hires_scalef"
    )
    .visiumhd_validate_json_names(json_scalefactors, expected_json_names)
    json_scalefactors
}

.visiumhd_micron_scale <- function(json_scalefactors) {
    # Check if json_scalefactors is a list and contains the required field
    if (!is.list(json_scalefactors)) {
        stop("json_scalefactors must be a list")
    }

    if (!"microns_per_pixel" %in% names(json_scalefactors)) {
        stop("microns_per_pixel field is missing from json_scalefactors")
    }

    # Extract the microns_per_pixel value
    json_scalefactors$microns_per_pixel
}

.visiumhd_get_image_type <- function(path) {
    img_name <- basename(path)
    possible_types <- c("lowres", "hires")
    for (image_type in possible_types) {
        if (grepl(image_type, img_name)) {
            return(image_type)
        }
    }
    return("unknown")
}

.visiumhd_get_scale_factor <- function(image_type, json_info) {
    if (is.null(json_info)) {
        warning("No scalefactors json info provided. VisiumHD image
                scale_factor defaulting to 1.")
        return(1)
    }

    checkmate::assert_list(json_info)

    scale_factor <- switch(image_type,
        "lowres" = json_info[["tissue_lowres_scalef"]],
        "hires" = json_info[["tissue_hires_scalef"]],
        "unknown" = {
            vmsg("[image] not hires or lowres, assuming image scale = 1")
            1
        },
        stop("Unexpected image type: ", image_type)
    )

    if (is.null(scale_factor)) {
        stop(
            "Scale factor for ", image_type,
            " image not found in json_info."
        )
    }

    return(scale_factor)
}

.visiumhd_image <- function(path,
    outdir = NULL,
    bin = 8,
    image_type = NULL,
    micron = FALSE,
    scalefactors_path = NULL,
    force_untar = FALSE,
    untar_params = list(),
    verbose = NULL) {
    checkmate::assert_list(untar_params)

    if (missing(path)) {
        stop(wrap_txt(
            "[VisiumHD] No path to image file provided or auto-detected"
        ), call. = FALSE)
    }
    
    if (!file.exists(path)) { # works for files and dirs
        stop("[VisiumHD] filepath does not exist\n", call. = FALSE)
    }
    
    vmsg(.v = verbose, "loading image ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # determine image type
    if (.is_tar(path) || checkmate::test_directory_exists(path)) {
        # image_type only needed if path is to tarfile or the extracted data
        if (is.null(image_type)) image_type <- "hires" # default
        image_type <- match.arg(image_type, c("hires", "lowres"))
        # get image filepath
        img_path <- switch(image_type,
            "hires" = "tissue_hires_image.png",
            "lowres" = "tissue_lowres_image.png"
        )
        img_params <- untar_params
        img_params$path <- path
        img_params$outdir <- outdir
        img_params$outdir_subpath <- file.path("spatial", img_path)
        img_params$bin <- bin
        img_params$force <- force_untar
        img_params$verbose = verbose
        img_params$spatial_only <- TRUE
        path <- do.call(.visiumhd_untar_if_not, img_params)
    } else {
        # otherwise should be filepath -> pattern match from basename
        image_type <- .visiumhd_get_image_type(path)
    }
    
    scalef_params <- untar_params
    scalef_params$path <- scalefactors_path
    scalef_params$outdir <- outdir
    scalef_params$bin <- bin
    scalef_params$force_untar <- force_untar
    scalef_params$verbose <- verbose
    scalef_params$verbose2 <- FALSE
    json_info <- do.call(.visiumhd_scalefactors, scalef_params)
    
    scalef <- .visiumhd_get_scale_factor(image_type, json_info)
    
    # load image
    gimg <- createGiottoLargeImage(path,
        name = "image",
        negative_y = TRUE,
        scale_factor = 1 / scalef
    )
    if (micron) {
        px2um <- .visiumhd_micron_scale(json_info)
        gimg <- rescale(gimg, px2um, x0 = 0, y0 = 0)
    }
    gimg
}

# new to spaceranger V4 outputs
# path can either be 'segmented_outputs' dir or direct to .geojson.
# if direct .geojson, scalefactor_path is needed separately if micron = TRUE
.visiumhd_poly <- function(path, 
    type = c("cell", "nucleus"),
    name = NULL,
    graphclust_annotated = FALSE,
    scalefactors_path = NULL,
    barcodes = NULL,
    micron = FALSE,
    id_fmt = "cellid_%09d-1",
    flip_vertical = TRUE,
    verbose = NULL) {
    
    type <- match.arg(type, c("cell", "nucleus"))
    if (is.null(name)) name <- type
    checkmate::assert_character(id_fmt, len = 1L)
    checkmate::assert_character(name, len = 1L)
    checkmate::assert_character(barcodes, null.ok = TRUE)
    flip_vertical <- as.logical(flip_vertical)
    micron <- as.logical(micron)
    
    if (missing(path)) {
        stop(wrap_txt(
            "No path to polygons provided or auto-detected"
        ), call. = FALSE)
    }
    
    vmsg(.v = verbose, sprintf("loading %s%s polygons ...", 
        type, ifelse(graphclust_annotated, " graphclust", "")))
    vmsg(.v = verbose, .is_debug = TRUE, path)

    if (".geojson" %in% GiottoUtils::file_extension(path)) { # direct filepath
        poly_path <- path
        if (micron && is.null(scalefactors_path)) {
            stop("[VisiumHD] <polygons> `scalefactors_path` needed for micron = TRUE\n",
                 call. = FALSE)
        } 
    } else if (checkmate::test_directory_exists(path)) { # output dir
        pattern <- c(type, "segmentations.geojson")
        if (graphclust_annotated) pattern <- c("graphclust_annotated", pattern)
        pattern <- paste(pattern, collapse = "_")
        poly_path <- .detect_in_dir(path,
            pattern = pattern, 
            recursive = FALSE, 
            platform = "VisiumHD"
        )
        scalefactors_path <- scalefactors_path %null%
            .detect_in_dir(path,
                pattern = "scalefactors_json", 
                recursive = TRUE,
                platform = "VisiumHD"
            )
    } else { # unrecognized
        stop("[VisiumHD] <polygons> unrecognized filepath\n", call. = FALSE)
    }
    
    # if still no path found, abort
    if (is.null(poly_path)) {
        stop(sprintf("[VisiumHD] <polygons> no %s%s polygons file found",
            type, ifelse(graphclust_annotated, " graphclust", "")),
            call. = FALSE)
    }

    # createGiottoPolygon() seems to have some difficulty with the default ids
    sv <- terra::vect(poly_path, crs = "local")
    
    if (micron) {
        px2um <- GiottoUtils::read_json(scalefactors_path)$microns_per_pixel
        sv <- rescale(sv, px2um, x0 = 0, y0 = 0)
    }
    
    # .geojson comes with an attribute col called cell_id with integers
    sv$cell_id <- sprintf(id_fmt, as.integer(sv$cell_id))
    
    # assemble giottopolygon
    p <- createGiottoPolygon(sv,
        name = name,
        calc_centroids = TRUE,
        verbose = FALSE
    )
    
    if (flip_vertical) p <- flip(p, y0 = 0)
    
    if (!is.null(barcodes)) {
        p <- p[barcodes]
    }
    p
}

# run tessellate() based on extent of VisiumHD data
# path should be a path for tissue_positions
.visiumhd_tessellate <- function(tissue_positions_path,
    shape = "hexagon",
    shape_size = 400,
    name = sprintf("%s%d", shape, as.integer(shape_size)),
    bin = 2,
    # only needed for micron
    micron = FALSE,
    scalefactors_path = NULL,
    # optional tar params
    outdir = NULL,
    force_untar = FALSE,
    untar_params = list(),
    verbose = NULL) {
    shape <- match.arg(shape, c("hexagon", "square"))
    checkmate::assert_list(untar_params)
    if (shape_size <= 0) {
        stop("Size must be a positive number.\n", call. = FALSE)
    }
    
    if (missing(tissue_positions_path)) {
        stop(wrap_txt(
            "No path to tissue positions file provided or auto-detected"
        ), call. = FALSE)
    }
    
    vmsg(.v = verbose, "building tessellated polygons ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    tp_params <- untar_params
    tp_params$path <- tissue_positions_path
    tp_params$outdir <- outdir
    tp_params$bin <- bin
    tp_params$micron <- micron
    tp_params$scalefactors_path <- scalefactors_path
    tp_params$force_untar <- force_untar
    tp_params$verbose <- verbose
    tp_params$verbose2 <- FALSE
    tp_params$output = "arrow"
    tp_arrow <- do.call(.visiumhd_tissue_positions, tp_params)
    e <- .visiumhd_px_extent(tp_arrow)

    vmsg(.v = verbose, paste("Creating a", shape, "tessellation with size", shape_size))
    tessellate(
        extent = e,
        shape = shape,
        shape_size = shape_size,
        name = name
    )
}

.visiumhd_cellmeta <- function(tissue_positions_path,
    outdir = NULL,
    bin = 8,
    barcodes = NULL,
    tissue_only = FALSE,
    force_untar = FALSE,
    untar_params = list(),
    verbose = NULL) {
    checkmate::assert_list(untar_params)
    # check if path is provided
    if (missing(tissue_positions_path)) {
        stop(wrap_txt(
            "No path to tissue positions file provided or auto-detected"
        ), call. = FALSE)
    }
    
    vmsg(.v = verbose, "loading cell metadata ...")
    vmsg(.v = verbose, .is_debug = TRUE, tissue_positions_path)

    # only cellmeta info is in the tissue positions.
    # Keep the array row/col and in_tissue info
    tp_params <- untar_params
    tp_params$path <- tissue_positions_path
    tp_params$outdir <- outdir
    tp_params$bin <- bin
    tp_params$barcodes <- barcodes
    tp_params$tissue_only <- tissue_only
    tp_params$force_untar <- force_untar
    tp_params$verbose <- verbose
    tp_params$verbose2 <- FALSE
    tp_params$output <- "arrow"
    
    tp_arrow <- do.call(.visiumhd_tissue_positions, tp_params)
    cx_dt <- tp_arrow |>
        dplyr::select(c("cell_ID" = "barcode", "in_tissue", "array_col", "array_row")) |>
        dplyr::collect() |>
        data.table::as.data.table()

    su <- sprintf("bin%03d", as.integer(bin))
    createCellMetaObj(cx_dt,
        spat_unit = su,
        feat_type = "rna",
        provenance = su,
        verbose = verbose
    )
}

.visiumhd_transcript <- function(expr_path,
    feature_id_type = c("symbols", "ensembl"),
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    expression_source = NULL,
    tissue_positions_path,
    scalefactors_path = NULL,
    micron = FALSE,
    outdir = NULL,
    bin = 2,
    tissue_only = FALSE,
    barcodes = NULL,
    array_subset_row = c(500, 1000),
    array_subset_col = c(500, 1000),
    pxl_subset_row = NULL,
    pxl_subset_col = NULL,
    filter = NULL,
    filter_coverage_cutoff = 0.5,
    force_untar = FALSE,
    untar_params = list(),
    output = c("giottoPoints", "full"),
    expr_list = NULL, # alternate input (list of exprObj)
    spatlocs = NULL, # alternate input (spatLocsObj)
    verbose = NULL) {
    checkmate::assert_list(untar_params)
    output <- match.arg(output, c("giottoPoints", "full"))
    # check if path is provided
    if (missing(expr_path)) {
        stop(wrap_txt(
            "No path to expression value file provided or auto-detected"
        ), call. = FALSE)
    }
    if (missing(tissue_positions_path)) {
        stop(wrap_txt(
            "No path to tissue positions file provided or auto-detected"
        ), call. = FALSE)
    }
    
    vmsg(.v = verbose, "creating point detections ...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # common params
    untar_params$outdir <- outdir
    untar_params$bin <- bin
    untar_params$force_untar <- force_untar
    untar_params$verbose <- verbose
    untar_params$verbose2 <- FALSE
    
    # get spatlocs
    if (!is.null(spatlocs)) {
        checkmate::assert_class(spatlocs, "spatLocsObj")
        sl <- spatlocs
    } else {
        tp_params <- untar_params
        tp_params$path <- tissue_positions_path
        tp_params$micron <- micron
        tp_params$scalefactors_path <- scalefactors_path
        tp_params$barcodes <- barcodes
        tp_params$tissue_only <- tissue_only
        tp_params$array_subset_row <- array_subset_row
        tp_params$array_subset_col <- array_subset_col
        tp_params$pxl_subset_row <- pxl_subset_row
        tp_params$pxl_subset_col <- pxl_subset_col
        tp_params$filter <- filter
        tp_params$filter_coverage_cutoff <- filter_coverage_cutoff
        tp_params$output <- "spatLocsObj"
        
        sl <- do.call(.visiumhd_tissue_positions, tp_params)
    }
    
    sl_barcodes <- spatIDs(sl)
    
    # get expression info
    if (!is.null(expr_list)) {
        checkmate::assert_list(expr_list)
        list_classes <- vapply(expr_list, class, FUN.VALUE = character(1L))
        if (!all(list_classes == "exprObj")) {
            stop("'expr_list' input must be a list of `exrpObj`\n", 
                 call. = FALSE)
        }
    } else {
        expr_params <- untar_params
        expr_params$path <- expr_path
        expr_params$feature_id_type <- feature_id_type
        expr_params$remove_zero_rows <- remove_zero_rows
        expr_params$split_by_type <- split_by_type
        expr_params$expression_source <- expression_source
        expr_params$barcodes <- sl_barcodes
        
        expr_list <- do.call(.visiumhd_expression, expr_params) # this is a list
        # possibly more than 1 exprObj with same barcodes, but maybe multiple
        # feat_types
    }
    
    # match expression barcodes back onto spatlocs just in case
    ex_barcodes <- spatIDs(expr_list[[1L]])
    sl <- sl[ex_barcodes]

    # create pts
    res_list <- lapply(expr_list, function(x) {
        # convert expression matrix to minimal data.table object
        mat <- x[]
        matrix_tile_dt <- data.table::as.data.table(Matrix::summary(mat))
        matrix_tile_dt[, feat_ID := mat@Dimnames[[1L]][i]]
        matrix_tile_dt[, pixel := mat@Dimnames[[2L]][j]]
        
        # get data.table from spatlocs
        spatlocs_dt <- sl[]
        
        # merge data.table matrix and spatial coordinates to create input for
        # Giotto Polygons
        gpoints <- data.table::merge.data.table(matrix_tile_dt, spatlocs_dt,
            by.x = "pixel", by.y = "cell_ID"
        )
        gpoints <- gpoints[,
            c("sdimx", "sdimy", "feat_ID", "x", "pixel"), with = FALSE
        ]
        data.table::setnames(gpoints, old = "x", new = "count")
        createGiottoPoints(gpoints, feat_type = featType(x), verbose = FALSE)
    })
    
    if (output == "full") {
        res <- list(
            "expr" = expr_list,
            "tissue_positions" = sl,
            "gpoints" = res_list
        )
    } else {
        res <- res_list
    }
    res
}

# ** untarring ####

.is_tar <- function(path) {
    fext <- GiottoUtils::file_extension(path)
    any(fext %in% "tar")
}

# work with one of:
# 1. fullpath to data of interest
# 2. fullpath to dir of extracted data
# 3. fullpath to tar file -> expected location of data in extracted tar data
# Get the filepath of an extracted element from binlevel's tar file. If not
# extracted yet, perform the extraction then return the extracted filepath
#   params:
# path - filepath to binlevel tar OR extracted file
# ... - additional params for untar()
.visiumhd_untar_if_not <- function(path,
    # untar only
    outdir = NULL,
    outdir_subpath,
    bin = 8, 
    force = FALSE, 
    verbose = NULL, 
    spatial_only = FALSE,
    ...) {
    if (!file.exists(path)) { # checks if path is file OR dir
        stop("[VisiumHD] filepath does not exist\n", call. = FALSE)
    }

    .complete_path <- function(dir) { # for tar and dir cases
        fullpath <- file.path(dir, outdir_subpath)
        if (!file.exists(fullpath)) {
            stop("[VisiumHD] <untar_if_not> filepath does not exist\n",
                 call. = FALSE)
        }
        fullpath
    }
    
    # paths for extracted data
    if (!.is_tar(path)) {
        if (checkmate::test_file_exists(path)) return(path)
        if (checkmate::test_directory_exists(path)) {
            return(.complete_path(path))
        }
    }
    
    # paths based on .tar file
    if (spatial_only) {
        subdir <- .visiumhd_untar_bin_spatial(path,
            outdir = outdir,
            bin = bin,
            force = force,
            verbose = verbose,
            ...
        )
    } else {
        subdir <- .visiumhd_untar_bin(path,
            outdir = outdir,
            bin = bin,
            force = force,
            verbose = verbose,
            ...
        )
    }
    
    .complete_path(subdir)
}

# Get the subdirectory path of a bin level. If the extracted version of the bin
# data does not exist yet, extract it at `outdir` then return the finalized
# extracted subdirectory location.
#   params:
# f - filepath to tar file, NOT an expected file output of the tar
# outdir - output directory (defaults to the same directory as the .tar)
# bin - what bin to pull from (needed for square_* subdirectory naming)
# pattern - a specific grep pattern to search for inside the tar
# force - whether to force a re-extract.
#   without this, untarring only happens if data subdir has less than 2 items
#   at the top level.
# verbose - be verbose
# ... additional params for untar()
.visiumhd_untar_bin <- function(f, 
    outdir = NULL, bin = 8, pattern = NULL, force = FALSE, verbose = NULL, ...) {
    checkmate::assert_file_exists(f, extension = c("tar", "tar.gz")) # f must be tar
    checkmate::assert_character(pattern, null.ok = TRUE)
    bin <- as.integer(bin)

    if (is.null(outdir)) outdir <- dirname(f)
    # toplevel of tar is expected to be a single folder with format
    # - "square_%03dum"
    # More than one folder is expected underneath, unless only the spatial
    # info was extracted
    
    # expected path
    data_subdir <- file.path(outdir, sprintf("square_%03dum", bin))
    
    # check if data already has been extracted
    # spatial may already be extracted separately from everything else
    # set 2 as the default min
    data_info <- list.files(data_subdir, include.dirs = TRUE, recursive = FALSE)
    if (setequal(data_info, "spatial")) force <- TRUE
    
    if (!force) return(data_subdir) # return early if no need to untar
    
    files <- NULL # set a default
    if (!is.null(pattern)) {
        tar_manifest <- untar(f, list = TRUE, ...)
        files <- tar_manifest[grepl(pattern = pattern, x = tar_manifest)]
    }
    
    vmsg(.v = verbose, "untarring...\n", f)
    untar(f, exdir = outdir, files = files, ...)
    data_subdir
}

# untar only the spatial subdirectory for a bin level.
.visiumhd_untar_bin_spatial <- function(f,
    outdir = NULL, bin = 8, force = FALSE, verbose = NULL, ...) {
    
    bin <- as.integer(bin)
    if (is.null(outdir)) outdir <- dirname(f)
    
    # expected path
    data_subdir <- file.path(outdir, sprintf("square_%03dum", bin))
    data_info <- list.files(data_subdir, include.dirs = TRUE, recursive = FALSE)
    # force if spatial data is not present.
    if (!"spatial" %in% data_info) force <- TRUE
    if (!force) return(data_subdir)
    
    .visiumhd_untar_bin(f,
        outdir = outdir,
        bin = bin,
        pattern = "spatial",
        force = force,
        verbose = verbose,
        ...
    )
}

# ** poly filtering ####

# create a representative 2 layer raster encoding row and col array values
# layer 1 named row, layer 2 named col
.visiumhd_array_raster <- function(tissue_positions_arrow) {
    dim <- .visiumhd_array_dim(tissue_positions_arrow = tissue_positions_arrow)
    corners <- .visiumhd_array_corners(tissue_positions_arrow = tissue_positions_arrow)
    checkmate::assert_integerish(dim, len = 2L)
    checkmate::assert_data_frame(corners, 
        max.rows = 4L, max.cols = 2L, types = "numeric"
    )

    # setup xy raster with default extent
    r <- terra::rast(
        nrows = dim[[1L]], 
        ncols = dim[[2L]],
        extent = ext(0, dim[[2L]], -dim[[1L]], 0),
        nlyr = 2L
    )
    terra::values(r) <- cbind( # values are zero indexed
        rep(seq_len(dim[[1L]]), each = dim[[2L]]) - 1L,
        rep(seq_len(dim[[2L]]), times = dim[[1L]]) - 1L
    )
    names(r) <- c("row", "col")
    gimg <- createGiottoLargeImage(r)
    
    # affine transform xy raster to match corners
    # get centroids of corners
    rast_corners <- data.frame(
        x = c(0.5, dim[[2L]] - 0.5, 0.5, dim[[2L]] - 0.5), # L,R,L,R
        y = c(-0.5, -0.5, -dim[[1L]] + 0.5, -dim[[1L]] + 0.5) # T,T,B,B
    )
    a <- calculateAffineMatrixFromLandmarks(
        source = rast_corners,
        target = corners
    )
    affine(gimg, a, pre_multiply = TRUE)
}

.visiumhd_array_corners <- function(tissue_positions_arrow) {
    # Get array dimensions to know the corner positions
    dims <- .visiumhd_array_dim(tissue_positions_arrow)
    # Convert back to zero index
    max_row <- dims[1] - 1L
    max_col <- dims[2] - 1L
    
    # Query for the 4 corner positions
    tissue_positions_arrow |>
        dplyr::filter(
            (array_row == 0L & array_col == 0L) |           # top-left
            (array_row == 0L & array_col == max_col) |      # top-right
            (array_row == max_row & array_col == 0L) |      # bottom-left
            (array_row == max_row & array_col == max_col)   # bottom-right
        ) |>
        dplyr::arrange(array_row, array_col) |>  # Consistent ordering
        dplyr::select("pxl_col_in_fullres", "pxl_row_in_fullres") |> # xy
        dplyr::collect() |>
        as.data.frame()
}

.visiumhd_array_dim <- function(tissue_positions_arrow) {
    result <- tissue_positions_arrow |>
        dplyr::summarize(
            nr = max(array_row, na.rm = TRUE),
            nc = max(array_col, na.rm = TRUE)
        ) |>
        dplyr::collect()
    
    c(result$nr, result$nc) + 1L # zero indexed
}

.visiumhd_px_extent <- function(tissue_positions_arrow, ext = TRUE) {
    if (all(c("sdimx", "sdimy") %in% names(tissue_positions_arrow))) {
        result <- tissue_positions_arrow |>
            dplyr::summarize(
                xmin = min(sdimx, na.rm = TRUE),
                xmax = max(sdimx, na.rm = TRUE),
                ymin = min(sdimy, na.rm = TRUE),
                ymax = max(sdimy, na.rm = TRUE)
            ) |>
            dplyr::collect()
    } else {
        result <- tissue_positions_arrow |>
            dplyr::summarize(
                xmin = min(pxl_col_in_fullres, na.rm = TRUE),
                xmax = max(pxl_col_in_fullres, na.rm = TRUE),
                ymin = min(pxl_row_in_fullres, na.rm = TRUE),
                ymax = max(pxl_row_in_fullres, na.rm = TRUE)
            ) |>
            dplyr::collect()
    }

    e <- c(result$xmin, result$xmax, result$ymin, result$ymax)
    if (ext) e <- ext(e)
    e
}

# filter xy indices by extracting xy values from a representative raster
# * micron scaled filters are not addressed here
#   see `.visiumhd_tissue_positions()` instead
# 
# params:
# - tissue_positions_arrow: (lazy arrow query to tissue_positions.parquet)
# - filter: polygon annotation
# - coverage_cutoff: coverage percentage of raster cell by annotation to be selected
# - output: type of output
# - verbose: verbosity
# 
# returns: 
# - data.table with integer cols c("row", "col)
# - visiumhd barcodes
# - lazy arrow query
.visiumhd_array_xy_filter <- function(tissue_positions_arrow,
    bin = 8L,
    filter = NULL,
    coverage_cutoff = 0.5,
    output = c("arrow", "barcode", "data.table"),
    verbose = NULL) {
    if (is.null(filter)) return(NULL) # skip if null
    output <- match.arg(output, c("arrow", "barcode", "data.table"))
    vmsg(.v = verbose, "filtering based on provided polygon...")
    GiottoUtils::package_check("exactextractr")

    # filter checks
    # ensure filter is SpatVector geomtype == polygons of length 1
    if (!inherits(filter, c("SpatVector", "giottoPolygon", "sf"))) {
        stop("[VisiumHD] filter must be giottoPolygon, sf, or SpatVector\n",
             call. = FALSE)
    }
    if (inherits(filter, "sf")) {
        if (!sf::st_geometry_type(filter) == "POLYGONS") {
            stop("[VisiumHD] if filter is sf, it must be st_geometry_type = \"POLYGONS\"\n",
                 call. = FALSE)
        }
        filter <- as.terra(filter)
    }
    if (inherits(filter, "SpatVector")) {
        if (!terra::geomtype(filter) == "polygons") {
            stop("[VisiumHD] if filter is SpatVector, it must be geomtype = \"polygons\"\n",
                 call. = FALSE)
        }
    }
    if (nrow(filter) > 1L) {
        warning("[VisiumHD] only first poly in 'filter' will be used.", call. = FALSE)
        filter <- filter[1]
    }
    if (inherits(filter, "giottoPolygon")) {
        filter <- filter[] # drop to SpatVector
    }
    
    coverage_fraction <- row <- col <- NULL # NSE vars

    # array xy representation
    r <- .visiumhd_array_raster(tissue_positions_arrow)
    # xy array index extraction - returns as 1 item list
    aff <- r@affine # get affine info
    filter_aff <- affine(filter, aff, inv = TRUE) # apply reverse affine to polys
    filter_aff <- as.sf(filter_aff) # coerce to sf after affine for exactextractr
    res <- exactextractr::exact_extract(r@raster_object, filter_aff)[[1L]]
    
    data.table::setDT(res)
    res[, row := as.integer(row)]
    res[, col := as.integer(col)]

    # coverage handling
    if (coverage_cutoff != 0) {
        res <- res[coverage_fraction >= coverage_cutoff]
    }
    # keep only 'row' and 'col' info
    res <- res[, c("row", "col"), with = FALSE]
    if (output == "data.table") return(res)
    
    filtered_arrow <- tissue_positions_arrow |>
        dplyr::filter(barcode %in% .visiumhd_barcode(res$row, res$col, bin = bin))
    if (output == "arrow") return(filtered_arrow)
    
    # barcode output
    dplyr::pull(filtered_arrow, "barcode", as_vector = TRUE)
}

.visiumhd_barcode <- function(row, col, bin = 2L) {
    sprintf("s_%03dum_%05d_%05d-%d", 
        as.integer(bin), 
        as.integer(row), 
        as.integer(col), 
        1L
    )
}

# wrapper ####

#' @title Create 10x VisiumHD Giotto Object
#' @name createGiottoVisiumHDObject
#' @description Given the path to a VisiumHD output folder, creates a
#' Giotto object. For lower-level and independent loading of specific pieces
#' of data, see [importVisiumHD()].
#' 
#' VisiumHD is a sequencing array-based spatial transcriptomic
#' assay. The array is composed of regular 2x2 micron "pixels". In addition
#' to the 2 micron resolution, lower resolution 8 and 16 binned outputs are
#' also provided. Individual pixels are referred to by barcode identifiers.
#' @param visiumhd_dir filepath to the exported visiumHD directory
#' @param bin numeric. One of 2, 8, 16. Which binning resolution to load
#' expression and spatial locations from.
#' @param micron logical. Set `TRUE` to load in micron scale instead of fullres 
#' image mapping.
#' @param load_expression logical. Whether to load in expression matrix
#' @param load_spatlocs logical. Whether to load in spatial locations
#' @param load_metadata logical. Whether to include array row/col and in_tissue
#' metadata information.
#' @param load_image logical. Whether to load in paired image information.
#' @param load_transcripts logical. Whether to load in bin 2 micron data as
#' transcripts information. Very memory expensive. Using `filter` or subsets
#' on array or pixel row/col is recommended.
#' @param create_tessellated_polys logical. Whether to generate tessellated
#' polys across the dataset.
#' @param tissue_only logical. Whether to only load information tagged as
#' `in_tissue` (in `tissue_positions.parquet`). This is ignored by transcript
#' loading.
#' @param barcodes character. Specific pixel barcodes to keep.
#' @param array_subset_row numeric vector, length = 2. Min/max of array rows to
#' keep
#' @param array_subset_col numeric vector, length = 2. Min/max of array cols to
#' keep
#' @param pxl_subset_row numeric vector, length = 2. Min/max of fullres image
#' mapped rows to keep. Note that values are inverted into the negatives.
#' @param pxl_subset_col numeric vector, length = 2. Min/max of fullres image
#' mapped cols to keep.
#' @param filter a `SpatVector`, `sf`, or `giottoPolygon` to spatially filter
#' the data by.
#' @param filter_coverage_cutoff numeric between 0 and 1. Minimal fraction of 
#' pixel coverage by `filter` in order to be selected.
#' @param expression_source character. One of `"raw"` or `"filtered"`.
#' Designates whether to pull expression value from raw or filtered matrix
#' outputs. These refer to whether they only include `in_tissue` pixels.
#' @param feature_id_type character. One of `"symbol"` or `"ensembl"`.
#' Determines which to use as the feature identifiers.
#' @param expression_remove_zero_rows logical (default = `TRUE`). Whether to
#' remove features with no detections.
#' @param expression_split_by_type logical (default = `TRUE`). Whether to
#' split expression information (and generated transcripts) by feature types
#' in the dataset (if multi modalities present).
#' @param image_type character. One of `"hires"` (default) or `"lowres"`.
#' Determines which image output to load. Ignored if `image_path` is provided.
#' Fullres image should be created separately with `createGiottoLargeImage()`
#' and attached with `setGiotto()`.
#' @param tessellate_shape character, One of `"hexagon"` or `"square"`. Shape
#' of poly to tessellate if `create_tessellated_polys = TRUE`
#' @param tessellate_shape_size numeric. Size of shape to tessellate 
#' (see [GiottoClass::tessellate()]).
#' @param tessellate_name name of tessellated polygons to create.
#' @param tissue_positions_path (optional) filepath to `tissue_positions.parquet`.
#' @param scalefactors_path (optional) filepath to `scalefactors_json.json`
#' @param expression_path (optional) filepath to .h5 or matrix market directory
#' @param image_path (optional) filepath to image to use
#' @param outdir (optional) directory to unpack bin output tar contents into.
#' (Default is the same directory as the tarfile.)
#' @param force_untar logical. Whether to force a untarring operation. Useful
#' when files have changed or a previous untar operation was interrupted.
#' @param untar_params list. Additional named params to pass to [untar()].
#' @param instructions giotto instructions to apply.
#' @param verbose verbosity
#' @returns giotto object
#' @examples
#' if (FALSE) {
#' data_dir <- "path/to/visiumhd/dir/containing/tar/outputs"
#' g <- createGiottoVisiumHDObject(data_dir, tissue_only = TRUE, bin = 16)
#' 
#' # example spatial filter (this would have to be relative to your fullres)
#' sv <- createSpatLocsObj(c(1e4, -1.5e4)) |> 
#'     as.points() |> 
#'     buffer(5000)
#' g <- createGiottoVisiumHDObject(data_dir, 
#'     tissue_only = TRUE, 
#'     load_transcripts = TRUE, 
#'     bin = 16, 
#'     filter = sv
#' )
#' }
#' @seealso [importVisiumHD()]
#' @export
createGiottoVisiumHDObject <- function(visiumhd_dir,
    bin = 8,
    micron = FALSE,
    load_expression = TRUE, 
    load_spatlocs = TRUE, 
    load_metadata = TRUE, 
    load_image = TRUE,
    load_transcripts = FALSE, 
    create_tessellated_polys = FALSE,
    tissue_only = FALSE, # ignored by tx
    barcodes = NULL, # ignored by tx
    array_subset_row = NULL, 
    array_subset_col = NULL, 
    pxl_subset_row = NULL, 
    pxl_subset_col = NULL,
    filter = NULL, 
    filter_coverage_cutoff = 0.5,
    expression_source = "raw", 
    feature_id_type = c("symbol", "ensembl"),
    expression_remove_zero_rows = TRUE, 
    expression_split_by_type = TRUE, 
    image_type = "hires", 
    tessellate_shape = "hexagon", 
    tessellate_shape_size = 400,
    tessellate_name = sprintf("%s%d", 
        tessellate_shape, as.integer(tessellate_shape_size)),
    tissue_positions_path = NULL,
    scalefactors_path = NULL,
    expression_path = NULL,
    image_path = NULL,
    outdir = NULL,
    force_untar = FALSE,
    untar_params = list(),
    instructions = NULL,
    verbose = NULL) {
    
    reader <- importVisiumHD(
        visiumhd_dir = visiumhd_dir,
        bin = bin,
        micron = micron,
        outdir = outdir,
        expression_source = expression_source,
        feature_id_type = feature_id_type,
        tissue_only = tissue_only,
        barcodes = barcodes,
        array_subset_row = array_subset_row,
        array_subset_col = array_subset_col,
        pxl_subset_row = pxl_subset_row,
        pxl_subset_col = pxl_subset_col,
        filter = filter,
        filter_coverage_cutoff = filter_coverage_cutoff
    )
    
    read_args <- list(
        load_expression = load_expression,
        load_spatlocs = load_spatlocs,
        load_metadata = load_metadata,
        load_transcripts = load_transcripts,
        load_image = load_image,
        create_tessellated_polys = create_tessellated_polys,
        expression_remove_zero_rows = expression_remove_zero_rows,
        expression_split_by_type = expression_split_by_type,
        image_type = image_type,
        tessellate_shape = tessellate_shape,
        tessellate_shape_size = tessellate_shape_size,
        tessellate_name = tessellate_name,
        force_untar = force_untar,
        untar_params = untar_params,
        instructions = instructions,
        verbose = verbose
    )
    
    if (!is.null(tissue_positions_path)) {
        read_args$tissue_positions_path <- tissue_positions_path
    }
    if (!is.null(scalefactors_path)) {
        read_args$scalefactors_path <- scalefactors_path
    }
    if (!is.null(expression_path)) {
        read_args$expression_path <- expression_path
    }
    if (!is.null(image_path)) {
        read_args$image_path <- image_path
    }
    
    do.call(reader$create_gobject, read_args)
}

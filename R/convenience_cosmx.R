# CLASS ####


setClass(
    "CosmxReader",
    slots = list(
        cosmx_dir = "character",
        version = "character",
        slide = "numeric",
        fovs = "numeric",
        micron = "logical",
        px2um = "numeric",
        offsets = "ANY",
        calls = "list"
    ),
    prototype = list(
        version = "default",
        slide = 1,
        micron = FALSE,
        px2um = 0.12028, # from cosmx output help files
        offsets = NULL,
        calls = list()
    )
)

# * show ####
setMethod("show", signature("CosmxReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "CosmxReader"))
    print_slots <- c("version", "dir", "slide", "fovs",
                     "micron", "offsets", "funs")
    pre <- sprintf(
        "%s :", format(print_slots)
    )
    names(pre) <- print_slots

    # dir
    d <- object@cosmx_dir
    if (length(d) > 0L) {
        nch <- nchar(d)
        d <- GiottoUtils::str_abbreviate(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }
    
    # version
    v <- object@version
    cat(pre["version"], v, "\n")

    # slide
    slide <- object@slide
    cat(pre["slide"], slide, "\n")

    # fovs
    fovs <- object@fovs %none% "all"
    cat(pre["fovs"], paste(fovs, collapse = ", "), "\n")

    # micron scaling
    micron <- ifelse(object@micron, object@px2um, FALSE)
    cat(pre["micron"], micron, "\n")

    # offsets
    offs_status <- ifelse(nrow(object@offsets) > 0L, "found", "none")
    cat(pre["offsets"], offs_status, "\n")

    # funs
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("CosmxReader"), function(x, ...) show(x))

# * plot ####
setMethod(
    "plot", signature(x = "CosmxReader", y = "missing"),
    function(x, cex = 0.8, ...) {
        a <- list(...)
        dat <- x@offsets

        if (x@micron) {
            px2um <- x@px2um
            dat$x <- dat$x * px2um
            dat$y <- dat$y * px2um
        }
        
        if (is.null(dat)) { # don't run if no offsets
            cat("no offsets to plot\n")
            return(invisible(NULL))
        }

        plot(y ~ x, data = dat, asp = 1L, type = "n", ...)
        text(y ~ x, data = dat, labels = dat$fov, cex = cex, ...)
    }
)




#' @title Import a Nanostring CosMx Assay
#' @name importCosMx
#' @description
#' Giotto import functionalities for CosMx datasets. This function generates
#' a `CosmxReader` instance that has convenient reader functions for converting
#' individual pieces of CosMx data into Giotto-compatible representations when
#' the params `cosmx_dir` and `fovs` (if only a subset is desired) are provided.
#' A function that creates the full `giotto` object is also available.
#' These functions should have all param values provided as defaults, but
#' can be flexibly modified to do things such as look in alternative
#' directories or paths.
#' @param cosmx_dir CosMx output directory
#' @param slide numeric. Slide number. Defaults to 1
#' @param fovs numeric. (optional) If provided, will load specific fovs.
#' Otherwise, all FOVs will be loaded
#' @param micron logical. Whether to scale spatial information as micron
#' instead of the default pixels
#' @param px2um numeric. Scalefactor from pixels to micron. Defaults to 0.12028
#' based on `CosMx-ReadMe.html` info
#' @details
#' Loading functions are generated after the `cosmx_dir` is added.
#' Transcripts, expression, and metadata loading are all expected to be done
#' from the top level of the directory. Loading of polys, and any image sets
#' are expected to be from specific subdirectories containing only those
#' images for the set of FOVs.
#' @returns CosmxReader object
#' @examples
#' # Create a `CosmxReader` object
#' reader <- importCosMx()
#'
#' \dontrun{
#' # Set the cosmx_dir and fov parameters
#' reader$cosmx_dir <- "path to cosmx dir"
#' reader$fov <- c(1, 4)
#'
#' plot(reader) # displays FOVs (top left corner) in px scale.
#'
#' # Load polygons, transcripts, and images
#' polys <- reader$load_polys()
#' tx <- reader$load_transcripts()
#' imgs <- reader$load_images()
#'
#' # Create a `giotto` object and add the loaded data
#' g <- giotto()
#' g <- setGiotto(g, tx[["rna"]])
#' g <- setGiotto(g, polys)
#' g <- addGiottoLargeImage(g, largeImages = imgs)
#' force(g)
#' }
#' @export
importCosMx <- function(cosmx_dir = NULL, slide = 1, fovs = NULL,
    micron = FALSE, px2um = 0.12028) {
    # get params
    a <- list(Class = "CosmxReader")
    if (!is.null(cosmx_dir)) {
        a$cosmx_dir <- cosmx_dir
    }
    if (!is.null(fovs)) {
        a$fovs <- fovs
    }
    a$slide <- slide
    a$micron <- micron
    a$px2um <- px2um

    do.call(new, args = a)
}

# * init ####
setMethod(
    "initialize", signature("CosmxReader"),
    function(.Object, cosmx_dir, version, slide, fovs, micron, px2um) {
        # provided params (if any)
        if (!missing(cosmx_dir)) {
            checkmate::assert_directory_exists(cosmx_dir)
            .Object@cosmx_dir <- cosmx_dir
        }
        if (!missing(version)) {
            version <- match.arg(tolower(version), 
                c("default", "v6", "legacy")
            )
            .Object@version <- version
        }
        if (!missing(slide)) {
            .Object@slide <- slide
        }
        if (!missing(fovs)) {
            .Object@fovs <- fovs
        }
        if (!missing(micron)) {
            .Object@micron <- micron
        }
        if (!missing(px2um)) {
            .Object@px2um <- px2um
        }

        # NULL case
        if (length(.Object@cosmx_dir) == 0) {
            return(.Object) # return early if no path given
        }


        # detect paths and subdirs
        p <- .Object@cosmx_dir
        .cosmx_detect <- function(pattern) {
            .detect_in_dir(pattern = pattern, path = p, platform = "CosMx")
        }

        v <- .Object@version
        if (v == "default") v <- "v6"

        shifts_path <- .cosmx_detect("fov_positions_file")
        meta_path <- .cosmx_detect("metadata_file")
        tx_path <- .cosmx_detect("tx_file")
        mask_dir <- .cosmx_detect("CellLabels")
        expr_path <- .cosmx_detect("exprMat_file")
        composite_img_dir <- .cosmx_detect("CellComposite")
        overlay_img_dir <- .cosmx_detect("CellOverlay")
        compart_img_dir <- .cosmx_detect("CompartmentLabels")


        # load fov offsets through one of several methods
        if (is.null(.Object@offsets)) { # only run if not already existing
            pos <- NULL

            if (!is.null(shifts_path)) {
                pos <- try(.cosmx_fov_shift(shifts_path), silent = TRUE)
                if (inherits(pos, "try-error")) pos <- NULL
            }

            # proceed with other possible methods of inferring shifts if present
            if (!is.null(meta_path) && is.null(pos)) {
                pos <- .cosmx_infer_fov_shifts(
                    meta_dt = data.table::fread(meta_path),
                    flip_loc_y = TRUE
                )
            } else if (!is.null(tx_path) && is.null(pos)) {
                warning(wrap_txt(
                    "metadata_file not found:
                Detecting fov shifts from tx_file. (This is slower)"
                ), call. = FALSE)
                pos <- .cosmx_infer_fov_shifts(
                    tx_dt = data.table::fread(tx_path),
                    flip_loc_y = TRUE
                )
            } else if (is.null(pos)) {
                pos <- data.table::data.table()
                warning(wrap_txt(
                    "NO FOV SHIFTS.
                fov_positions_file, tx_file,
                and metadata_file not auto detected.
                One of these must be provided to infer FOV shifts.\n
                Alternatively, directly supply a data.table with:
                fov(int), x(numeric), y(numeric) in px scaling to `$offsets`"
                ), call. = FALSE)
            }

            .Object@offsets <- pos
        }



        # transcripts load call
        dcols <- switch(v,
            "v6" = c("x_local_px", "y_local_px", "cell_ID", "cell"),
            "legacy" = c("x_local_px", "y_local_px", "cell_ID")
        )
        tx_fun <- function(
        path = tx_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        dropcols = dcols,
        cores = determine_cores(),
        verbose = NULL) {
            .cosmx_transcript(
                path = path,
                fovs = .Object@fovs %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                dropcols = dropcols,
                micron = .Object@micron,
                px2um = .Object@px2um,
                cores = cores,
                verbose = verbose
            )
        }
        .Object@calls$load_transcripts <- tx_fun



        # mask load call
        vstep <- switch(v,
            "v6" = FALSE,
            "legacy" = 1
        )
        
        mask_fun <- function(
        path = mask_dir,
        # VERTICAL FLIP + NO VERTICAL SHIFT
        flip_vertical = TRUE,
        flip_horizontal = FALSE,
        shift_vertical_step = vstep,
        shift_horizontal_step = FALSE,
        remove_background_polygon = TRUE,
        verbose = NULL) {
            .cosmx_poly(
                path = path,
                fovs = .Object@fovs %none% NULL,
                flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal,
                shift_vertical_step = shift_vertical_step,
                shift_horizontal_step = shift_horizontal_step,
                remove_background_polygon = remove_background_polygon,
                micron = .Object@micron,
                px2um = .Object@px2um,
                offsets = .Object@offsets,
                verbose = verbose
            )
        }
        .Object@calls$load_polys <- mask_fun


        # expression load call
        expr_fun <- function(
        path = expr_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb")) {
            .cosmx_expression(
                path = path,
                fovs = .Object@fovs %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword
            )
        }
        .Object@calls$load_expression <- expr_fun


        # images load call
        negy <- switch(v,
            "v6" = TRUE,
            "legacy" = FALSE,
        )
        
        img_fun <- function(
        path = composite_img_dir,
        img_type = "composite",
        img_name_fmt = paste0(img_type, "_fov%03d"),
        negative_y = negy,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        verbose = NULL) {
            .cosmx_image(
                path = path,
                fovs = .Object@fovs %none% NULL,
                img_type = img_type,
                img_name_fmt = img_name_fmt,
                negative_y = negative_y,
                flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal,
                micron = .Object@micron,
                px2um = .Object@px2um,
                offsets = .Object@offsets,
                verbose = verbose
            )
        }
        .Object@calls$load_images <- img_fun


        # meta load call
        meta_fun <- function(
        path = meta_path,
        cores = determine_cores(),
        dropcols = c(
            "CenterX_local_px",
            "CenterY_local_px",
            "CenterX_global_px",
            "CenterY_global_px",
            "cell_id"
        ),
        verbose = NULL) {
            .cosmx_cellmeta(
                path = path,
                fovs = .Object@fovs %none% NULL,
                dropcols = dropcols,
                cores = cores,
                verbose = verbose
            )
        }
        .Object@calls$load_cellmeta <- meta_fun


        # build gobject call
        gobject_fun <- function(
        transcript_path = tx_path,
        cell_labels_dir = mask_dir,
        expression_path = expr_path,
        metadata_path = meta_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list(
            "NegPrb"
        ),
        load_images = list(
            composite = "composite",
            overlay = "overlay"
        ),
        load_expression = FALSE,
        load_cellmeta = TRUE,
        instructions = NULL,
        cores = determine_cores(),
        verbose = NULL) {
            load_expression <- as.logical(load_expression)
            load_cellmeta <- as.logical(load_cellmeta)

            if (!is.null(load_images)) {
                checkmate::assert_list(load_images)
                if (is.null(names(load_images))) {
                    stop("Images directories provided to
                    'load_images' must be named")
                }
            }

            funs <- .Object@calls

            # init gobject
            g <- giotto()
            if (!is.null(instructions)) {
                instructions(g) <- instructions
            }

            # transcripts
            tx_list <- funs$load_transcripts(
                path = transcript_path,
                feat_type = feat_type,
                split_keyword = split_keyword,
                cores = cores,
                verbose = verbose
            )
            for (tx in tx_list) {
                g <- setGiotto(g, tx)
            }

            # polys
            polys <- funs$load_polys(
                path = cell_labels_dir,
                verbose = FALSE
            )
            g <- setGiotto(g, polys, verbose = verbose)

            # images
            if (!is.null(load_images)) {
                # replace convenient shortnames
                load_images[load_images == "composite"] <- composite_img_dir
                load_images[load_images == "overlay"] <- overlay_img_dir

                imglist <- list()
                dirnames <- names(load_images)
                for (imdir_i in seq_along(load_images)) {
                    dir_imgs <- funs$load_images(
                        path = load_images[[imdir_i]],
                        img_type = dirnames[[imdir_i]],
                        verbose = verbose
                    )
                    imglist <- c(imglist, dir_imgs)
                }
                g <- addGiottoLargeImage(g, largeImages = imglist, 
                                         verbose = FALSE)
            }

            # expression & meta
            # Need to check that names agree for poly/expr/meta
            allowed_ids <- spatIDs(polys)

            if (load_expression) {
                exlist <- funs$load_expression(
                    path = expression_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword
                )

                # only keep allowed cells and set into gobject
                for (ex in exlist) {
                    bool <- colnames(ex[]) %in% allowed_ids
                    ex[] <- ex[][, bool]
                    g <- setGiotto(g, ex, verbose = verbose)
                }
            }

            if (load_cellmeta) {
                cx <- funs$load_cellmeta(
                    path = metadata_path,
                    cores = cores,
                    verbose = verbose
                )

                cx[] <- cx[][cell_ID %in% allowed_ids, ]
                g <- setGiotto(g, cx, verbose = verbose)
            }
            
            # create spatlocs
            g <- addSpatialCentroidLocations(g, verbose = FALSE)
            
            # add fovs
            g$fov <- gsub("^c_\\d+_(\\d+)_\\d+$", "\\1", pDataDT(g)$cell_ID)

            return(g)
        }
        .Object@calls$create_gobject <- gobject_fun

        return(.Object)
    }
)





# * access ####

#' @export
setMethod("$", signature("CosmxReader"), function(x, name) {
    basic_info <- c("cosmx_dir", "version", "slide", "fovs", 
                    "micron", "px2um", "offsets")
    if (name %in% basic_info) {
        return(methods::slot(x, name))
    }

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("CosmxReader"), function(x, name, value) {
    basic_info <- c("cosmx_dir", "version", "slide", "fovs", "micron", "px2um")
    if (name %in% basic_info) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }

    if (name == "offsets") {
        value <- data.table::setDT(value)
        data.table::setnames(value, new = c("fov", "x", "y"))
        methods::slot(x, name) <- value
        return(initialize(x))
    }

    stop(sprintf(
        "Only items in '%s' can be set",
        paste0(basic_info, collapse = "', '")
    ))
})

#' @export
`.DollarNames.CosmxReader` <- function(x, pattern) {
    dn <- c("cosmx_dir", "version", "slide", "fovs",
            "micron", "px2um", "offsets")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}





# MODULAR ####

.cosmx_fov_shift <- function(path) {
    if (is.null(path)) return(NULL) # return early (empty case)
    fov_shifts <- data.table::fread(path)
    fs_colnames <- colnames(fov_shifts)
    if (!"X_mm" %in% fs_colnames) return(NULL) # return early (WTX)
    
    offset_colnames <- c("fov", "x", "y")
    
    if (identical(fs_colnames, 
        c("fov", "x_global_px", "y_global_px")
    )) {
        data.table::setnames(fov_shifts, new = offset_colnames)
    } else if (identical(fs_colnames,
        c("FOV", "x_global_px", "y_global_px", "x_global_mm", "y_global_mm")
    )) {
        fov_shifts <- fov_shifts[, c("FOV", "x_global_px", "y_global_px")]
        data.table::setnames(fov_shifts, new = offset_colnames)
    } else {
        return(NULL) # unknown read method
    }
    fov_shifts
}

.cosmx_transcript <- function(
        path,
        fovs = NULL,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        dropcols = c(
            "x_local_px",
            "y_local_px",
            "cell_ID",
            "cell"
        ),
        micron = FALSE,
        px2um = 0.12028,
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to tx file provided or auto-detected"
        ), call. = FALSE)
    }

    checkmate::assert_file_exists(path)

    vmsg(.v = verbose, "loading feature detections...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    tx <- data.table::fread(input = path, nThread = cores, drop = dropcols)
    if (!is.null(fovs)) {
        # subset to only needed FOVs
        tx <- tx[fov %in% as.numeric(fovs), ]
    }

    # micron scaling if desired
    if (micron) {
        tx[, x_global_px := x_global_px * px2um]
        tx[, y_global_px := y_global_px * px2um]
    }

    # giottoPoints ----------------------------------------------------- #

    # static gpoints params
    gpoints_params <- list()
    gpoints_params$feat_type <- feat_type
    gpoints_params$split_keyword <- split_keyword
    gpoints_params$x_colname <- "x_global_px"
    gpoints_params$y_colname <- "y_global_px"
    gpoints_params$feat_ID_colname <- "target"
    gpoints_params$verbose = FALSE

    gpoints <- do.call(createGiottoPoints, c(list(x = tx), gpoints_params))
    # ensure output is always a list
    if (!is.list(gpoints)) {
        gpoints <- list(gpoints)
        names(gpoints) <- objName(gpoints[[1L]])
    }

    return(gpoints)
}

#' @name .cosmx_infer_fov_shifts
#' @title Infer CosMx local to global shifts
#' @description
#' From NanoString CosMx spatial info, infer the FOV shifts needed. These
#' values are needed for anything that requires the use of images, since those
#' do not come with spatial extent information embedded.
#' @param tx_dt transcript data.table input to use
#' (Only one of tx_dt or meta_dt should be used)
#' @param meta_dt cell metadata data.table input to use
#' (Only one of tx_dt or meta_dt should be used)
#' @param navg max n values to check per FOV to find average shift
#' @param flip_loc_y whether a y flip needs to be performed on the local y
#' values before comparing with global y values. See details
#' @returns data.table with three columns. 1. FOV (integer), xshift (numeric),
#' yshift (numeric). Values should always be in pixels
#' @details
#' Shifts are found by looking at the average of differences between xy global
#' and local coordinates in either the metadata or transcripts file. The number
#' of shift value to average across is determined with `navg`. The average is
#' in place to get rid of small differences in shifts, likely due to rounding
#' errors. Across the different versions of the CosMx exports, whether the
#' local y values are flipped compared to the global values has differed, so
#' there is also a step that checks the variance of y values per sampled set
#' per fov. In cases where the shift is calculated with the correct (inverted
#' or non-inverted) y local values, the variance is expected to be very low.
#' When the variance is higher than 0.001, the function is re-run with the
#' opposite `flip_loc_y` value.
#' @keywords internal
.cosmx_infer_fov_shifts <- function(tx_dt, meta_dt,
    flip_loc_y = TRUE, navg = 100L) {
    fov <- NULL # NSE vars
    if (!missing(tx_dt)) {
        tx_head <- tx_dt[, head(.SD, navg), by = fov]
        x <- tx_head[, mean(x_global_px - x_local_px), by = fov]
        if (flip_loc_y) {
            # test if flip is needed
            # Usual yshift variance / fov expected when correct is 0 to 1e-22
            # if var is too high for any fov, swap `flip_loc_y` value
            y <- tx_head[, var(y_global_px + y_local_px), by = fov]
            if (y[, any(V1 > 0.001)]) {
                return(.cosmx_infer_fov_shifts(
                    tx_dt = tx_dt, flip_loc_y = FALSE, navg = navg
                ))
            }

            # use +y if local y values are flipped
            y <- tx_head[, mean(y_global_px + y_local_px), by = fov]
        } else {
            y <- tx_head[, mean(y_global_px - y_local_px), by = fov]
        }
    } else if (!missing(meta_dt)) {
        meta_head <- meta_dt[, head(.SD, navg), by = fov]
        x <- meta_head[, mean(CenterX_global_px - CenterX_local_px), by = fov]
        if (flip_loc_y) {
            # test if flip is needed
            # Usual yshift variance / fov expected when correct is 0 to 1e-22
            # if var is too high for any fov, swap `flip_loc_y` value
            y <- meta_head[
                , var(CenterY_global_px + CenterY_local_px),
                by = fov
            ]
            if (y[, any(V1 > 0.001)]) {
                return(.cosmx_infer_fov_shifts(
                    meta_dt = meta_dt, flip_loc_y = FALSE, navg = navg
                ))
            }

            # use +y if local y values are flipped
            y <- meta_head[, mean(CenterY_global_px + CenterY_local_px),
                by = fov
            ]
        } else {
            y <- meta_head[, mean(CenterY_global_px - CenterY_local_px),
                by = fov
            ]
        }
    } else {
        stop("One of tx_dt or meta_dt must be provided\n")
    }

    res <- merge(x, y, by = "fov")
    data.table::setnames(res, new = c("fov", "x", "y"))

    return(res)
}

.cosmx_imgname_fovparser <- function(path) {
    im_names <- list.files(path)
    fovs <- as.numeric(sub(".*F(\\d+)\\..*", "\\1", im_names))
    if (any(is.na(fovs))) {
        warning(wrap_txt(
            "Images to load should be sets of images/fov in subdirectories.
            No other files should be present."
        ))
    }
    return(fovs)
}

.cosmx_poly <- function(
        path,
        slide = 1,
        fovs = NULL,
        name = "cell",
        # VERTICAL FLIP + NO SHIFTS
        flip_vertical = TRUE,
        flip_horizontal = FALSE,
        shift_vertical_step = FALSE,
        shift_horizontal_step = FALSE,
        remove_background_polygon = TRUE,
        micron = FALSE,
        px2um = 0.12028,
        offsets,
        verbose = NULL) {
    # NSE params
    f <- x <- y <- NULL

    if (missing(path)) {
        stop(wrap_txt(
            "No path to polys subdirectory provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, "loading segmentation masks...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    mask_params <- list(
        # static params
        mask_method = "multiple",
        # A background poly for nanostring masks sometimes shows up.
        # removal works by looking for any polys with size more than 90% of the
        # total FOV along either x or y axis
        remove_background_polygon = remove_background_polygon,
        fill_holes = TRUE,
        calc_centroids = TRUE,
        remove_unvalid_polygons = TRUE,
        # input params
        name = name,
        flip_vertical = flip_vertical,
        flip_horizontal = flip_horizontal,
        shift_vertical_step = shift_vertical_step,
        shift_horizontal_step = shift_horizontal_step,
        verbose = FALSE
    )

    fovs <- fovs %null% .cosmx_imgname_fovparser(path) # ALL if NULL
    with_pbar({
        p <- pbar(along = fovs)

        gpolys <- lapply(fovs, function(f) {
            segfile <- Sys.glob(paths = sprintf("%s/*F%03d*", path, f))
            # naming format: c_SLIDENUMBER_FOVNUMBER_CELLID
            mask_params$ID_fmt <- paste0(
                sprintf("c_%d_%d_", slide, f), "%d"
            )

            gpoly <- do.call(
                createGiottoPolygonsFromMask,
                args = c(list(maskfile = segfile), mask_params)
            )

            xshift <- offsets[fov == f, x]
            yshift <- offsets[fov == f, y]

            # if micron scale
            if (micron) {
                gpoly <- rescale(
                    gpoly,
                    fx = px2um, fy = px2um, x0 = 0, y0 = 0
                )
                xshift <- xshift * px2um
                yshift <- yshift * px2um
            }

            gpoly <- spatShift(x = gpoly, dx = xshift, dy = yshift)
            p(message = sprintf("F%03d", f))
            return(gpoly)
        })
    })

    if (length(gpolys) > 1L) {
        gpolys <- do.call(rbind, args = gpolys)
    }

    # never return lists. Only the single merged gpoly
    return(gpolys)
}

.cosmx_cellmeta <- function(
        path,
        slide = 1,
        fovs = NULL,
        dropcols = c(
            "CenterX_local_px",
            "CenterY_local_px",
            "CenterX_global_px",
            "CenterY_global_px",
            "cell_id"
        ),
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to metadata file provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, "loading cell metadata...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    verbose <- verbose %null% TRUE

    meta_dt <- data.table::fread(input = path, nThread = cores)

    # remove unneeded cols
    dropcols <- dropcols[dropcols %in% colnames(meta_dt)]
    meta_dt[, (dropcols) := NULL] # remove dropcols

    # subset to needed fovs
    if (!is.null(fovs)) {
        fovs <- as.integer(fovs)
        meta_dt <- meta_dt[fov %in% fovs, ]
    }

    # create cell ID as `c_SLIDENUMBER_FOVNUMBER_CELLID`
    if ("cell" %in% colnames(meta_dt)) {
        # assume already formatted (current datasets Mar-27-2024)
        meta_dt[, c("fov", "cell_ID") := NULL]
        data.table::setnames(meta_dt, old = "cell", "cell_ID")
    } else {
        # older datasets
        meta_dt[, cell_ID := sprintf("c_%d_%d_%d", slide, fov, cell_ID)]
        # remove fov
        meta_dt[, fov := NULL]
    }


    # TODO figure out what to do about protein expression here.
    cx <- createCellMetaObj(
        metadata = meta_dt,
        spat_unit = "cell",
        feat_type = "rna",
        provenance = "cell",
        verbose = verbose
    )
    return(cx)
}

.cosmx_expression <- function(
        path,
        slide = 1,
        fovs = NULL,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to exprMat file provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, "loading expression matrix...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    expr_dt <- data.table::fread(input = path, nThread = cores)

    # subset to needed fovs
    if (!is.null(fovs)) {
        fovs <- as.integer(fovs)
        expr_dt <- expr_dt[fov %in% fovs, ]
    }

    # remove background values (cell 0)
    expr_dt <- expr_dt[cell_ID != 0L, ]

    # create cell ID as `c_SLIDENUMBER_FOVNUMBER_CELLID`
    expr_dt[, cell_ID := sprintf("c_%d_%d_%d", slide, fov, cell_ID)]
    # remove fov
    expr_dt[, fov := NULL]

    # convert to Matrix
    expr_mat <- dt_to_matrix(expr_dt)
    expr_mat <- t_flex(expr_mat)

    # split expression for rna / negprb if any split keywords provided.
    # Output of this chunk should always be a named list of 1 or more matrices
    if (length(split_keyword) > 0) {
        expr_list <- vector(mode = "list", length = length(feat_type))
        names(expr_list) <- feat_type
        # iterate through other expr types
        for (key_i in seq_along(split_keyword)) {
            feat_ids <- rownames(expr_mat)
            bool <- grepl(pattern = split_keyword[[key_i]], x = feat_ids)
            # subset and store split matrix
            sub_mat <- expr_mat[bool, ]
            expr_list[[key_i + 1L]] <- sub_mat
            # remaining matrix
            expr_mat <- expr_mat[!bool, ]
        }
        # assign the main expr
        expr_list[[1L]] <- expr_mat
    } else {
        expr_list <- list(expr_mat)
        names(expr_list) <- feat_type[[1L]]
    }

    expr_list <- lapply(seq_along(expr_list), function(expr_i) {
        createExprObj(
            expression_data = expr_list[[expr_i]],
            spat_unit = "cell",
            feat_type = names(expr_list)[[expr_i]],
            name = "raw",
            provenance = "cell"
        )
    })

    return(expr_list)
}

.cosmx_image <- function(
        path,
        fovs = NULL,
        img_type = "composite",
        img_name_fmt = paste(img_type, "_fov%03d"),
        negative_y = FALSE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        micron = FALSE,
        px2um = 0.12028,
        offsets,
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to image subdirectory to load provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, sprintf("loading %s images...", img_type))
    vmsg(.v = verbose, .is_debug = TRUE, path)

    fovs <- fovs %null% .cosmx_imgname_fovparser(path) # ALL if NULL
    verbose <- verbose %null% TRUE

    with_pbar({
        p <- pbar(along = fovs)

        gimg_list <- lapply(fovs, function(f) {
            imgfile <- Sys.glob(paths = sprintf("%s/*F%03d*", path, f))
            img_name <- sprintf(img_name_fmt, f)

            gimg <- createGiottoLargeImage(
                raster_object = imgfile,
                name = img_name,
                negative_y = negative_y,
                flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal,
                verbose = verbose
            )

            xshift <- offsets[fov == f, x]
            yshift <- offsets[fov == f, y]

            if (micron) {
                gimg <- rescale(
                    gimg,
                    fx = px2um, fy = px2um, x0 = 0, y0 = 0
                )
                xshift <- xshift * px2um
                yshift <- yshift * px2um
            }

            gimg <- spatShift(x = gimg, dx = xshift, dy = yshift)
            p(message = sprintf("F%03d", f))
            return(gimg)
        })
    })


    return(gimg_list)
}









# wrapper ####


#' @title Create Nanostring CosMx Giotto Object
#' @name createGiottoCosMxObject
#' @description Given the path to a CosMx experiment directory, creates a Giotto
#' object.
#' @param cosmx_dir full path to the exported cosmx directory
#' @param version character. Version of CosMx output. Current selections are
#' either "default", "v6", and "legacy" (for the NSCLC dataset).
#' @param FOVs field of views to load (only affects subcellular data and images)
#' @param feat_type character. feature type. Provide more than one value if
#' using the `split_keyword` param. For each set of keywords to split by, an
#' additional feat_type should be provided in the same order. Affects how
#' the transcripts information is loaded. Helpful for separating out the
#' QC probes. See details.
#' @param split_keyword list of character vectors of keywords to split the
#' transcripts based on their feat_ID. Keywords will be `grepl()`
#' matched against the feature IDs information. See details.
#' @param load_images named list of filepaths to directories. Loads the 
#' composite and overlay images by default.
#' @param load_expression logical. (Default = FALSE) whether to load provided
#' expression matrix
#' @param load_cellmeta logical. (Default = TRUE) whether to load provided
#' cell metadata
#' @param fov_shifts_path Optional. Filepath to fov_positions_file
#' @param transcript_path Optional. Filepath to desired transcripts file to
#' load.
#' @param cell_labels_dir Optional. Path to directory containing CellLabels
#' images to load as polygons.
#' @param expression_path Optional. Filepath to cell feature matrix to load.
#' @param metadata_path Optional. Filepath to metadata file to load.
#' @param cores nthreads to use when loading in cell metadata and transcripts
#' @param data_to_use deprecated. Not used
#' @param remove_background_polygon deprecated. Now always done
#' @param background_algo deprecated. Not used
#' @param remove_unvalid_polygons deprecated. Now always done
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @inheritDotParams importCosMx -cosmx_dir -fovs
#' @returns a giotto object
#' @details
#' \[**Expected Directory**\] This function generates a giotto object when
#' given a link to a cosmx output directory. It expects the following items
#' within the directory where the \strong{bolded} portions are what this
#' function matches against:
#' \itemize{
#'   \item{\strong{CellComposite} (folder of images)}
#'   \item{\strong{CellLabels} (folder of images)}
#'   \item{\strong{CellOverlay} (folder of images)}
#'   \item{\strong{CompartmentLabels} (folder of images)}
#'   \item{experimentname_\strong{exprMat_file}.csv (file)}
#'   \item{experimentname_\strong{fov_positions_file}.csv (file)}
#'   \item{experimentname_\strong{metadata_file}.csv (file)}
#'   \item{experimentname_\strong{tx_file}.csv (file)}
#' }
#'
#' \[**feat_type and split_keyword**\]
#' Additional QC probe information is in the subcellular feature detections
#' information and must be separated from the gene expression information
#' during processing.
#' The QC probes have prefixes that allow them to be selected from the rest of
#' the feature IDs.
#' Giotto uses `feat_type` and `split_keyword` params to select these QC
#' probes out as separate feature types. See examples in
#' `[GiottoClass::createGiottoPoints]` for how this works.
#'
#' The Gene expression subset labeled as `rna` is accepted as the subset of
#' feat_IDs that do not get matched to any of the `split_keywords`.
#'
#' \[**Images**\] Images in the expected CellComposite and CellOverlay
#' folders will be loaded as giotto largeImage objects by default.
#' \code{\link{showGiottoImageNames}} can be used to see the available images.
#' @md
#' @export
createGiottoCosMxObject <- function(
        cosmx_dir,
        version = "default",
        FOVs = NULL,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        load_images = list(composite = "composite", 
                           overlay = "overlay"),
        load_expression = FALSE,
        load_cellmeta = FALSE,
        
        # optional filepaths
        fov_shifts_path = NULL,
        transcript_path = NULL,
        cell_labels_dir = NULL,
        expression_path = NULL,
        metadata_path = NULL,
        cores = determine_cores(),
        verbose = NULL,
        instructions = NULL,
        
        # deprecated params
        remove_unvalid_polygons = deprecated(), # not used
        data_to_use = deprecated(), # not used
        remove_background_polygon = deprecated(), # not used
        background_algo = deprecated(), # not used
        ...) {

    # handle deprecations
    if (!missing(remove_unvalid_polygons)) {
        deprecate_warn(
            "4.2.0", "createGiottoCosMxObject(remove_unvalid_polygons)")
    }
    if (!missing(remove_background_polygon)) {
        deprecate_warn(
            "4.2.0", "createGiottoCosMxObject(remove_background_polygon)")
    }
    if (!missing(data_to_use)) {
        deprecate_warn(
            "4.2.0", "createGiottoCosMxObject(data_to_use)")
    }
    if (!missing(background_algo)) {
        deprecate_warn(
            "4.2.0", "createGiottoCosMxObject(background_algo)")
    }

    
    # setup importer
    x <- importCosMx(cosmx_dir = cosmx_dir, fovs = FOVs, ...)
    x$version <- version
    if (!is.null(fov_shifts_path)) {
        checkmate::assert_file_exists(fov_shifts_path)
        x$offsets <- data.table::fread(fov_shifts_path)
    }
    
    # gobject creation args setup
    load_args <- list(
        feat_type = feat_type,
        split_keyword = split_keyword,
        load_images = load_images,
        instructions = instructions,
        cores = cores,
        verbose = verbose
    )
    # pass optional paths if provided
    if (!is.null(transcript_path)) {
        load_args$transcript_path <- transcript_path
    }
    if (!is.null(cell_labels_dir)) {
        load_args$cell_labels_dir <- cell_labels_dir
    }
    if (!is.null(expression_path)) {
        load_args$expression_path <- expression_path
    }
    if (!is.null(metadata_path)) {
        load_args$metadata_path <- metadata_path
    }
    
    do.call(x$create_gobject, load_args)
}




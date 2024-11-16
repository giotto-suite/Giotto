# CLASS ####


setClass(
    "CosmxReader",
    slots = list(
        cosmx_dir = "character",
        slide = "numeric",
        fovs = "numeric",
        micron = "logical",
        px2mm = "numeric",
        offsets = "ANY",
        calls = "list"
    ),
    prototype = list(
        slide = 1,
        micron = FALSE,
        px2mm = 0.12028, # from cosmx output help files
        offsets = NULL,
        calls = list()
    )
)

# * show ####
setMethod("show", signature("CosmxReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "CosmxReader"))
    print_slots <- c("dir", "slide", "fovs", "micron", "offsets", "funs")
    pre <- sprintf(
        "%s :", format(print_slots)
    )
    names(pre) <- print_slots

    # dir
    d <- object@cosmx_dir
    if (length(d) > 0L) {
        nch <- nchar(d)
        d <- abbrev_path(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }

    # slide
    slide <- object@slide
    cat(pre["slide"], slide, "\n")

    # fovs
    fovs <- object@fovs %none% "all"
    cat(pre["fovs"], paste(fovs, collapse = ", "), "\n")

    # micron scaling
    micron <- ifelse(object@micron, object@px2mm / 1000, FALSE)
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
#' @param px2mm numeric. Scalefactor from pixels to mm. Defaults to 0.12028
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
    micron = FALSE, px2mm = 0.12028) {
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
    a$px2mm <- px2mm

    do.call(new, args = a)
}

# * init ####
setMethod(
    "initialize", signature("CosmxReader"),
    function(.Object, cosmx_dir, slide, fovs, micron, px2mm) {
        # provided params (if any)
        if (!missing(cosmx_dir)) {
            checkmate::assert_directory_exists(cosmx_dir)
            .Object@cosmx_dir <- cosmx_dir
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
        if (!missing(px2mm)) {
            .Object@px2mm <- px2mm
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
                fov_shifts <- data.table::fread(shifts_path)
                if (!"X_mm" %in% colnames(fov_shifts)) {
                    # older version has fov, x, y (all numeric) in px shifts
                    data.table::setnames(fov_shifts, new = c("fov", "x", "y"))
                    pos <- fov_shifts
                }
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
            } else {
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
        tx_fun <- function(
        path = tx_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        dropcols = c(
            "x_local_px",
            "y_local_px",
            "cell_ID",
            "cell"
        ),
        verbose = NULL) {
            .cosmx_transcript(
                path = path,
                fovs = .Object@fovs %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                dropcols = dropcols,
                micron = .Object@micron,
                px2mm = .Object@px2mm,
                cores = determine_cores(),
                verbose = verbose
            )
        }
        .Object@calls$load_transcripts <- tx_fun



        # mask load call
        mask_fun <- function(
        path = mask_dir,
        # VERTICAL FLIP + NO VERTICAL SHIFT
        flip_vertical = TRUE,
        flip_horizontal = FALSE,
        shift_vertical_step = FALSE,
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
                px2mm = .Object@px2mm,
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
        img_fun <- function(
        path = composite_img_dir,
        img_type = "composite",
        img_name_fmt = paste0(img_type, "_fov%03d"),
        negative_y = TRUE,
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
                px2mm = .Object@px2mm,
                offsets = .Object@offsets,
                verbose = verbose
            )
        }
        .Object@calls$load_images <- img_fun


        # meta load call
        meta_fun <- function(
        path = meta_path,
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
                cores = determine_cores(),
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
        load_cellmeta = FALSE,
        instructions = NULL) {
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
                split_keyword = split_keyword
            )
            for (tx in tx_list) {
                g <- setGiotto(g, tx)
            }

            # polys
            polys <- funs$load_polys(
                path = cell_labels_dir,
                verbose = FALSE
            )
            g <- setGiotto(g, polys)

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
                    )
                    imglist <- c(imglist, dir_imgs)
                }
                g <- addGiottoLargeImage(g, largeImages = imglist)
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
                    g <- setGiotto(g, ex)
                }
            }

            if (load_cellmeta) {
                cx <- funs$load_cellmeta(
                    path = metadata_path
                )

                cx[] <- cx[][cell_ID %in% allowed_ids, ]
                g <- setGiotto(g, cx)
            }

            return(g)
        }
        .Object@calls$create_gobject <- gobject_fun

        return(.Object)
    }
)





# * access ####

#' @export
setMethod("$", signature("CosmxReader"), function(x, name) {
    basic_info <- c("cosmx_dir", "slide", "fovs", "micron", "px2mm", "offsets")
    if (name %in% basic_info) {
        return(methods::slot(x, name))
    }

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("CosmxReader"), function(x, name, value) {
    basic_info <- c("cosmx_dir", "slide", "fovs", "micron", "px2mm")
    if (name %in% basic_info) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }

    if (name == "offsets") {
        methods::slot(x, name) <- data.table::setDT(value)
        return(initialize(x))
    }

    stop(sprintf(
        "Only items in '%s' can be set",
        paste0(basic_info, collapse = "', '")
    ))
})

#' @export
`.DollarNames.CosmxReader` <- function(x, pattern) {
    dn <- c("cosmx_dir", "slide", "fovs", "micron", "px2mm", "offsets")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}





# MODULAR ####

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
        px2mm = 0.12028,
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
        px2micron <- px2mm / 1000
        tx[, x_global_px := x_global_px * px2micron]
        tx[, y_global_px := y_global_px * px2micron]
    }

    # giottoPoints ----------------------------------------------------- #

    # static gpoints params
    gpoints_params <- list()
    gpoints_params$feat_type <- feat_type
    gpoints_params$split_keyword <- split_keyword
    gpoints_params$x_colname <- "x_global_px"
    gpoints_params$y_colname <- "y_global_px"
    gpoints_params$feat_ID_colname <- "target"

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
        px2mm = 0.12028,
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
                px2micron <- px2mm / 1000
                gpoly <- rescale(
                    gpoly,
                    fx = px2micron, fy = px2micron, x0 = 0, y0 = 0
                )
                xshift <- xshift * px2micron
                yshift <- yshift * px2micron
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
        negative_y = TRUE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        micron = FALSE,
        px2mm = 0.12028,
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
                px2micron <- px2mm / 1000
                gimg <- rescale(
                    gimg,
                    fx = px2micron, fy = px2micron, x0 = 0, y0 = 0
                )
                xshift <- xshift * px2micron
                yshift <- yshift * px2micron
            }

            gimg <- spatShift(x = gimg, dx = xshift, dy = yshift)
            p(message = sprintf("F%03d", f))
            return(gimg)
        })
    })


    return(gimg_list)
}



#' @title Load CosMx folder subcellular info
#' @name .load_cosmx_folder_subcellular
#' @description loads in the feature detections information. Note that the mask
#' images are still required for a working subcellular object, and those are
#' loaded in \code{\link{.createGiottoCosMxObject_subcellular}}
#' @inheritParams createGiottoCosMxObject
#' @returns list
#' @keywords internal
.load_cosmx_folder_subcellular <- function(
        dir_items,
        FOVs = NULL,
        cores,
        verbose = TRUE) {
    vmsg(.v = verbose, "Loading subcellular information...")

    # subcellular checks
    if (!file.exists(dir_items$`transcript locations file`)) {
        stop(wrap_txt("No transcript locations file (.csv) detected"))
    }
    if (!file.exists(dir_items$`fov positions file`)) {
        stop(wrap_txt("No fov positions file (.csv) detected"))
    }

    # FOVs to load
    vmsg(.v = verbose, "Loading FOV offsets...")
    fov_offset_file <- fread(
        input = dir_items$`fov positions file`, nThread = cores
    )
    if (is.null(FOVs)) FOVs <- fov_offset_file$fov # default to ALL FOVs
    FOV_ID <- as.list(sprintf("%03d", FOVs))

    # TODO Load only relevant portions of file?

    vmsg(.v = verbose, "Loading transcript level info...")
    tx_coord_all <- fread(
        input = dir_items$`transcript locations file`, nThread = cores
    )
    vmsg(.v = verbose, "Subcellular load done")

    data_list <- list(
        "FOV_ID" = FOV_ID,
        "fov_offset_file" = fov_offset_file,
        "tx_coord_all" = tx_coord_all
    )

    return(data_list)
}



#' @title Load CosMx folder aggregate info
#' @name .load_cosmx_folder_aggregate
#' @inheritParams createGiottoCosMxObject
#' @returns list
#' @keywords internal
.load_cosmx_folder_aggregate <- function(
        dir_items,
        cores,
        verbose = TRUE) {
    # data.table vars
    fov <- cell_ID <- fov_cell_ID <- CenterX_global_px <-
        CenterY_global_px <- CenterX_local_px <-
        CenterY_local_px <- x_shift <- y_shift <- NULL

    # load aggregate information
    vmsg(.v = verbose, "Loading provided aggregated information...")

    # aggregate checks
    if (!file.exists(dir_items$`expression matrix file`)) {
        stop(wrap_txt("No expression matrix file (.csv) detected"))
    }
    if (!file.exists(dir_items$`metadata file`)) {
        stop(wrap_txt("No metadata file (.csv) detected. Needed for cell
                      spatial locations."))
    }

    # read in aggregate data
    expr_mat <- fread(
        input = dir_items$`expression matrix file`, nThread = cores
    )
    metadata <- fread(input = dir_items$`metadata file`, nThread = cores)

    # setorder expression and spatlocs
    data.table::setorder(metadata, fov, cell_ID)
    data.table::setorder(expr_mat, fov, cell_ID)


    # generate unique cell IDs
    expr_mat[, cell_ID := paste0(
        "fov", sprintf("%03d", fov), "-", "cell_", cell_ID
    )]
    expr_mat <- expr_mat[, fov := NULL]

    metadata[, fov_cell_ID := cell_ID]
    metadata[, cell_ID := paste0(
        "fov", sprintf("%03d", fov), "-", "cell_", cell_ID
    )]
    # reorder
    data.table::setcolorder(x = metadata, c("cell_ID", "fov", "fov_cell_ID"))


    # extract spatial locations
    spatlocs <- metadata[, .(CenterX_global_px, CenterY_global_px, cell_ID)]
    spatlocs_fov <- metadata[, .(CenterX_local_px, CenterY_local_px, cell_ID)]
    # regenerate FOV shifts
    metadata[, x_shift := CenterX_global_px - CenterX_local_px]
    metadata[, y_shift := CenterY_global_px - CenterY_local_px]
    fov_shifts <- metadata[, .(mean(x_shift), mean(y_shift)), fov]
    colnames(fov_shifts) <- c("fov", "x_shift", "y_shift")


    # rename spatloc column names
    spatloc_oldnames <- c("CenterX_global_px", "CenterY_global_px", "cell_ID")
    spatloc_oldnames_fov <- c("CenterX_local_px", "CenterY_local_px", "cell_ID")
    spatloc_newnames <- c("sdimx", "sdimy", "cell_ID")
    data.table::setnames(
        spatlocs,
        old = spatloc_oldnames, new = spatloc_newnames
    )
    data.table::setnames(
        spatlocs_fov,
        old = spatloc_oldnames_fov, new = spatloc_newnames
    )

    # cleanup metadata and spatlocs
    metadata <- metadata[, c(
        "CenterX_global_px", "CenterY_global_px",
        "CenterX_local_px", "CenterY_local_px"
    ) := NULL]
    # find unique cell_IDs present in both expression and metadata
    giotto_cell_ID <- unique(intersect(expr_mat$cell_ID, metadata$cell_ID))

    # subset to only unique cell_IDs
    expr_mat <- expr_mat[cell_ID %in% giotto_cell_ID, ]
    metadata <- metadata[cell_ID %in% giotto_cell_ID, ]


    # convert protein metadata to expr mat
    # take all mean intensity protein information except for
    # MembraneStain and DAPI
    protein_meta_cols <- colnames(metadata)
    protein_meta_cols <- protein_meta_cols[
        grepl(pattern = "Mean.*", x = protein_meta_cols)
    ]
    protein_meta_cols <- protein_meta_cols[
        !protein_meta_cols %in% c("Mean.MembraneStain", "Mean.DAPI")
    ]
    protein_meta_cols <- c("cell_ID", protein_meta_cols)

    prot_expr <- metadata[, protein_meta_cols, with = FALSE]
    prot_cell_ID <- metadata[, cell_ID]
    protM <- Matrix::Matrix(as.matrix(prot_expr[, -1]),
        dimnames = list(
            prot_expr[[1]],
            colnames(prot_expr[, -1])
        ),
        sparse = FALSE
    )
    protM <- t_flex(protM)

    # convert expression to sparse matrix
    spM <- Matrix::Matrix(as.matrix(expr_mat[, -1]),
        dimnames = list(
            expr_mat[[1]],
            colnames(expr_mat[, -1])
        ),
        sparse = TRUE
    )
    spM <- t_flex(spM)

    ## Ready for downstream aggregate gobject creation or appending into
    # existing subcellular Giotto object ##

    data_list <- list(
        "spatlocs" = spatlocs,
        "spatlocs_fov" = spatlocs_fov,
        "metadata" = metadata,
        "protM" = protM,
        "spM" = spM,
        "fov_shifts" = fov_shifts
    )

    return(data_list)
}









# OLD ####


#' @title Create Nanostring CosMx Giotto Object
#' @name createGiottoCosMxObject
#' @description Given the path to a CosMx experiment directory, creates a Giotto
#' object.
#' @param cosmx_dir full path to the exported cosmx directory
#' @param data_to_use which type(s) of expression data to build the gobject with
#' Default is \code{'all'} information available. \code{'subcellular'} loads
#' the transcript coordinates only. \code{'aggregate'} loads the provided
#' aggregated expression matrix.
#' @param FOVs field of views to load (only affects subcellular data and images)
#' @param remove_background_polygon try to remove background polygon
#' (default: FALSE)
#' @param background_algo algorithm to remove background polygon
#' @param remove_unvalid_polygons remove unvalid polygons (default: TRUE)
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @returns a giotto object
#' @details
#' [\strong{Expected Directory}] This function generates a giotto object when
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
#' [\strong{Workflows}] Workflow to use is accessed through the data_to_use
#' param
#' \itemize{
#'   \item{'all' - loads and requires subcellular information from tx_file and
#'   fov_positions_file
#'   and also the existing aggregated information
#'   (expression, spatial locations, and metadata)
#'   from exprMat_file and metadata_file.}
#'   \item{'subcellular' - loads and requires subcellular information from
#'   tx_file and
#'   fov_positions_file only.}
#'   \item{'aggregate' - loads and requires the existing aggregate information
#'   (expression, spatial locations, and metadata) from exprMat_file and
#'   metadata_file.}
#' }
#'
#' [\strong{Images}] Images in the default CellComposite, CellLabels,
#' CompartmentLabels, and CellOverlay
#' folders will be loaded as giotto largeImage objects in all workflows as
#' long as they are available. Additionally, CellComposite images will be
#' converted to giotto image objects, making plotting with
#' these image objects more responsive when accessing them from a server.
#' \code{\link{showGiottoImageNames}} can be used to see the available images.
#' @export
createGiottoCosMxObject <- function(
        cosmx_dir = NULL,
        data_to_use = c("all", "subcellular", "aggregate"),
        remove_background_polygon = TRUE,
        background_algo = c("range"),
        remove_unvalid_polygons = TRUE,
        FOVs = NULL,
        instructions = NULL,
        cores = determine_cores(),
        verbose = TRUE) {
    # 0. setup
    cosmx_dir <- path.expand(cosmx_dir)

    # determine data to use
    data_to_use <- match.arg(
        arg = data_to_use, choices = c("all", "subcellular", "aggregate")
    )
    if (data_to_use %in% c("all", "aggregate")) {
        stop(wrap_txt('Convenience workflows "all" and "aggregate" are not
                      available yet'))
    }

    # Define for data.table
    fov <- target <- x_local_px <- y_local_px <- z <- cell_ID <-
        CenterX_global_px <- CenterY_global_px <-
        CenterX_local_px <- CenterY_local_px <- NULL


    # 1. test if folder structure exists and is as expected
    dir_items <- .read_cosmx_folder(
        cosmx_dir = cosmx_dir,
        verbose = verbose
    )


    # 2. load and create giotto object
    cosmx_gobject <- switch(data_to_use,
        "subcellular" = .createGiottoCosMxObject_subcellular(
            dir_items,
            FOVs = FOVs,
            remove_background_polygon = remove_background_polygon,
            background_algo = background_algo,
            remove_unvalid_polygons = remove_unvalid_polygons,
            cores = cores,
            verbose = verbose,
            instructions = instructions
        ),
        "aggregate" = .createGiottoCosMxObject_aggregate(
            dir_items,
            cores = cores,
            verbose = verbose,
            instructions = instructions
        ),
        "all" = .createGiottoCosMxObject_all(
            dir_items,
            FOVs = FOVs,
            remove_background_polygon = remove_background_polygon,
            background_algo = background_algo,
            remove_unvalid_polygons = remove_unvalid_polygons,
            cores = cores,
            verbose = verbose,
            instructions = instructions
        )
    )


    # load in subcellular information, subcellular FOV objects, then join


    # load in pre-generated aggregated expression matrix
    if (data_to_use == "aggregate" | data_to_use == "all") {

    }

    vmsg(.v = verbose, "done")
    return(cosmx_gobject)
}



#' @title Load and create a CosMx Giotto object from subcellular info
#' @name .createGiottoCosMxObject_subcellular
#' @inheritParams createGiottoCosMxObject
#' @returns giotto object
#' @keywords internal
.createGiottoCosMxObject_subcellular <- function(
        dir_items,
        FOVs = NULL,
        remove_background_polygon = TRUE,
        background_algo = c("range"),
        remove_unvalid_polygons = TRUE,
        cores,
        verbose = TRUE,
        instructions = NULL) {
    target <- fov <- NULL

    # load tx detections and FOV offsets ------------------------------------- #
    data_list <- .load_cosmx_folder_subcellular(
        dir_items = dir_items,
        FOVs = FOVs,
        cores = cores,
        verbose = verbose
    )

    # unpack data_list
    FOV_ID <- data_list$FOV_ID
    fov_offset_file <- data_list$fov_offset_file
    tx_coord_all <- data_list$tx_coord_all

    # remove global xy values and cell_ID
    tx_coord_all[, c("x_global_px", "y_global_px", "cell_ID") := NULL]

    data.table::setcolorder(
        tx_coord_all, c("target", "x_local_px", "y_local_px", "z", "fov")
    )

    # feature detection type splitting --------------------------------------- #

    if (isTRUE(verbose)) message("Splitting detections by feature vs neg probe")
    all_IDs <- tx_coord_all[, unique(target)]
    neg_IDs <- all_IDs[grepl(pattern = "NegPrb", all_IDs)]
    feat_IDs <- all_IDs[!all_IDs %in% neg_IDs]

    # split detections DT
    feat_coords_all <- tx_coord_all[target %in% feat_IDs]
    neg_coords_all <- tx_coord_all[target %in% neg_IDs]

    if (isTRUE(verbose)) {
        message("  > Features: ", feat_coords_all[, .N])
        message("  > NegProbes: ", neg_coords_all[, .N])
    }

    # FOV-based processing --------------------------------------------------- #

    fov_gobjects_list <- lapply(FOV_ID, function(x) {
        # images --------------------------------------------------- #
        # build image paths
        if (isTRUE(verbose)) message("Loading image information...")

        composite_dir <- Sys.glob(paths = file.path(
            dir_items$`CellComposite folder`, paste0("*", x, "*")
        ))
        cellLabel_dir <- Sys.glob(paths = file.path(
            dir_items$`CellLabels folder`, paste0("*", x, "*")
        ))
        compartmentLabel_dir <- Sys.glob(paths = file.path(
            dir_items$`CompartmentLabels folder`, paste0("*", x, "*")
        ))
        cellOverlay_dir <- Sys.glob(paths = file.path(
            dir_items$`CellOverlay folder`, paste0("*", x, "*")
        ))

        # Missing warnings
        if (length(composite_dir) == 0) {
            warning("[ FOV ", x, " ] No composite images found")
            composite_dir <- NULL
        }
        if (length(cellLabel_dir) == 0) {
            stop("[ FOV ", x, " ] No cell mask images found")
        } # cell masks are necessary
        if (length(compartmentLabel_dir) == 0) {
            warning("[ FOV ", x, " ] No compartment label images found")
            compartmentLabel_dir <- NULL
        }
        if (length(cellOverlay_dir) == 0) {
            warning("[ FOV ", x, " ] No cell polygon overlay images found")
            cellOverlay_dir <- NULL
        }

        if (isTRUE(verbose)) message("Image load done")

        if (isTRUE(verbose)) wrap_msg("[ FOV ", x, "]")


        # transcripts ---------------------------------------------- #
        # get FOV specific tx locations
        if (isTRUE(verbose)) message("Assigning FOV feature detections...")


        # feature info
        coord_oldnames <- c("target", "x_local_px", "y_local_px")
        coord_newnames <- c("feat_ID", "x", "y")

        feat_coord <- feat_coords_all[fov == as.numeric(x)]
        data.table::setnames(
            feat_coord,
            old = coord_oldnames, new = coord_newnames
        )
        # neg probe info
        neg_coord <- neg_coords_all[fov == as.numeric(x)]
        data.table::setnames(
            neg_coord,
            old = coord_oldnames, new = coord_newnames
        )


        # build giotto object -------------------------------------- #
        if (isTRUE(verbose)) message("Building subcellular giotto object...")
        fov_subset <- createGiottoObjectSubcellular(
            gpoints = list(
                "rna" = feat_coord,
                "neg_probe" = neg_coord
            ),
            gpolygons = list("cell" = cellLabel_dir),
            polygon_mask_list_params = list(
                mask_method = "guess",
                flip_vertical = TRUE,
                flip_horizontal = FALSE,
                shift_horizontal_step = FALSE,
                remove_background_polygon = remove_background_polygon,
                background_algo = background_algo,
                remove_unvalid_polygons = remove_unvalid_polygons
            ),
            instructions = instructions,
            cores = cores
        )


        # find centroids as spatial locations ---------------------- #
        if (isTRUE(verbose)) {
            message("Finding polygon centroids as cell spatial locations...")
        }
        fov_subset <- addSpatialCentroidLocations(
            fov_subset,
            poly_info = "cell",
            spat_loc_name = "raw"
        )


        # create and add giotto image objects ---------------------- #
        if (isTRUE(verbose)) {
            message("Attaching image files...")
            print(composite_dir)
            print(cellOverlay_dir)
            print(compartmentLabel_dir)
        }

        gImage_list <- list()

        # load image if files are found
        if (!is.null(composite_dir)) {
            gImage_list$composite <- createGiottoLargeImage(
                raster_object = composite_dir,
                negative_y = FALSE,
                name = "composite"
            )
        }
        if (!is.null(cellOverlay_dir)) {
            gImage_list$overlay <- createGiottoLargeImage(
                raster_object = cellOverlay_dir,
                negative_y = FALSE,
                name = "overlay"
            )
        }
        if (!is.null(compartmentLabel_dir)) {
            gImage_list$compartment <- createGiottoLargeImage(
                raster_object = compartmentLabel_dir,
                negative_y = FALSE,
                name = "compartment"
            )
        } # TODO



        if (length(gImage_list) > 0) {
            fov_subset <- addGiottoImage(
                gobject = fov_subset,
                images = gImage_list
            )

            # convert to MG for faster loading (particularly relevant for
            # pulling from server)
            # TODO remove this
            fov_subset <- convertGiottoLargeImageToMG(
                giottoLargeImage = gImage_list$composite,
                gobject = fov_subset,
                return_gobject = TRUE,
                verbose = FALSE
            )
        } else {
            message("No images found for fov")
        }
    }) # lapply end

    # returning -------------------------------------------------------------- #

    if (length(FOVs) == 1) {
        return(fov_gobjects_list[[1]])
    } else {
        # join giotto objects according to FOV positions file
        if (isTRUE(verbose)) message("Joining FOV gobjects...")
        new_gobj_names <- paste0("fov", FOV_ID)
        id_match <- match(as.numeric(FOV_ID), fov_offset_file$fov)
        x_shifts <- fov_offset_file[id_match]$x_global_px
        y_shifts <- fov_offset_file[id_match]$y_global_px

        # Join giotto objects
        cosmx_gobject <- joinGiottoObjects(
            gobject_list = fov_gobjects_list,
            gobject_names = new_gobj_names,
            join_method = "shift",
            x_shift = x_shifts,
            y_shift = y_shifts
        )
        return(cosmx_gobject)
    }
}



#' @title Load and create a CosMx Giotto object from aggregate info
#' @name .createGiottoCosMxObject_aggregate
#' @inheritParams createGiottoCosMxObject
#' @returns giotto object
#' @keywords internal
.createGiottoCosMxObject_aggregate <- function(
        dir_items,
        cores,
        verbose = TRUE,
        instructions = NULL) {
    data_to_use <- fov <- NULL

    data_list <- .load_cosmx_folder_aggregate(
        dir_items = dir_items,
        cores = cores,
        verbose = verbose
    )

    # unpack data_list
    spatlocs <- data_list$spatlocs
    spatlocs_fov <- data_list$spatlocs_fov
    metadata <- data_list$metadata
    protM <- data_list$protM
    spM <- data_list$spM
    fov_shifts <- data_list$fov_shifts


    # create standard gobject from aggregate matrix
    if (data_to_use == "aggregate") {
        # Create aggregate gobject
        if (isTRUE(verbose)) message("Building giotto object...")
        cosmx_gobject <- createGiottoObject(
            expression = list("raw" = spM, "protein" = protM),
            cell_metadata = list("cell" = list(
                "rna" = metadata,
                "protein" = metadata
            )),
            spatial_locs = spatlocs,
            instructions = instructions,
            cores = cores
        )


        # load in images
        img_ID <- data.table::data.table(
            fov = fov_shifts[, fov],
            img_name = paste0(
                "fov",
                sprintf("%03d", fov_shifts[, fov]), "-image"
            )
        )

        if (isTRUE(verbose)) message("Attaching image files...")
        composite_dir <- Sys.glob(paths = file.path(
            dir_items$`CellComposite folder`, paste0("/*")
        ))
        cellLabel_dir <- Sys.glob(paths = file.path(
            dir_items$`CellLabels folder`, paste0("/*")
        ))
        compartmentLabel_dir <- Sys.glob(paths = file.path(
            dir_items$`CompartmentLabels folder`, paste0("/*")
        ))
        overlay_dir <- Sys.glob(paths = file.path(
            dir_items$`CellOverlay folder`, paste0("/*")
        ))

        if (length(cellLabel_imgList) > 0) {
            cellLabel_imgList <- lapply(cellLabel_dir, function(x) {
                createGiottoLargeImage(x, name = "cellLabel", negative_y = TRUE)
            })
        }
        if (length(composite_imgList) > 0) {
            composite_imgList <- lapply(composite_dir, function(x) {
                createGiottoLargeImage(x, name = "composite", negative_y = TRUE)
            })
        }
        if (length(compartmentLabel_dir) > 0) {
            compartmentLabel_imgList <- lapply(
                compartmentLabel_dir, function(x) {
                    createGiottoLargeImage(x,
                        name = "composite",
                        negative_y = TRUE
                    )
                }
            )
        }
        if (length(overlay_dir) > 0) {
            overlay_imgList <- lapply(overlay_dir, function(x) {
                createGiottoLargeImage(x,
                    name = "composite",
                    negative_y = TRUE
                )
            })
        }
    }
}




#' @title Load and create a CosMx Giotto object from subcellular and aggregate
#' info
#' @name .createGiottoCosMxObject_all
#' @param dir_items list of full directory paths from \code{.read_cosmx_folder}
#' @inheritParams createGiottoCosMxObject
#' @returns giotto object
#' @details Both \emph{subcellular}
#' (subellular transcript detection information) and
#' \emph{aggregate} (aggregated detection count matrices by cell polygon from
#' NanoString)
#' data will be loaded in. The two will be separated into "cell" and "cell_agg"
#' spatial units in order to denote the difference in origin of the two.
#' @seealso createGiottoCosMxObject .createGiottoCosMxObject_aggregate
#' .createGiottoCosMxObject_subcellular
#' @keywords internal
.createGiottoCosMxObject_all <- function(
        dir_items,
        FOVs,
        remove_background_polygon = TRUE,
        background_algo = "range",
        remove_unvalid_polygons = TRUE,
        cores,
        verbose = TRUE,
        instructions = NULL,
        ...) {
    # 1. create subcellular giotto as spat_unit 'cell'
    cosmx_gobject <- .createGiottoCosMxObject_subcellular(
        dir_items = dir_items,
        FOVs = FOVs,
        remove_background_polygon = remove_background_polygon,
        background_algo = background_algo,
        remove_unvalid_polygons = remove_unvalid_polygons,
        cores = cores,
        verbose = verbose,
        instructions = instructions
    )

    # 2. load and append aggregated information in spat_unit 'cell_agg'
    agg_data <- .load_cosmx_folder_aggregate(
        dir_items = dir_items,
        cores = cores,
        verbose = verbose
    )

    # unpack data_list
    spatlocs <- agg_data$spatlocs
    spatlocs_fov <- agg_data$spatlocs_fov
    metadata <- agg_data$metadata
    protM <- agg_data$protM
    spM <- agg_data$spM

    # add in pre-generated aggregated expression matrix information for 'all'
    # workflow

    # Add aggregate expression information
    if (isTRUE(verbose)) {
        wrap_msg(
            'Appending provided aggregate expression data as...
                               spat_unit: "cell_agg"
                               feat_type: "rna"
                               name: "raw"'
        )
    }
    # add expression data to expression slot
    s4_expr <- createExprObj(
        name = "raw",
        expression_data = spM,
        spat_unit = "cell_agg",
        feat_type = "rna",
        provenance = "cell_agg"
    )

    cosmx_gobject <- setGiotto(cosmx_gobject, s4_expr)

    # Add spatial locations
    if (isTRUE(verbose)) {
        wrap_msg(
            'Appending metadata provided spatial locations data as...
                               --> spat_unit: "cell_agg" name: "raw"
                               --> spat_unit: "cell" name: "raw_fov"'
        )
    }
    if (isTRUE(verbose)) {
        wrap_msg(
            'Polygon centroid derived spatial locations assigned as...
                               --> spat_unit: "cell" name: "raw" (default)'
        )
    }

    locsObj <- create_spat_locs_obj(
        name = "raw",
        coordinates = spatlocs,
        spat_unit = "cell_agg",
        provenance = "cell_agg"
    )
    locsObj_fov <- create_spat_locs_obj(
        name = "raw_fov",
        coordinates = spatlocs_fov,
        spat_unit = "cell_agg",
        provenance = "cell_agg"
    )

    cosmx_gobject <- setGiotto(cosmx_gobject, locsObj, initialize = FALSE)
    cosmx_gobject <- setGiotto(cosmx_gobject, locsObj_fov, initialize = FALSE)

    # initialize cell and feat IDs and metadata slots for 'cell_agg' spat_unit
    agg_cell_ID <- colnames(s4_expr[])
    agg_feat_ID <- rownames(s4_expr[])

    sub_feat_ID <- featIDs(cosmx_gobject, feat_type = "rna")
    feat_ID_new <- unique(c(agg_feat_ID, sub_feat_ID))

    # cell metadata

    # Add metadata to both the given and the poly spat_units
    if (isTRUE(verbose)) message("Appending provided cell metadata...")
    cosmx_gobject <- addCellMetadata(cosmx_gobject,
        spat_unit = "cell",
        feat_type = "rna",
        new_metadata = metadata,
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )
    cosmx_gobject <- addCellMetadata(cosmx_gobject,
        spat_unit = "cell_agg",
        feat_type = "rna",
        new_metadata = metadata,
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )

    initialize(cosmx_gobject)
}



#' @title Read a structured CosMx folder
#' @name .read_cosmx_folder
#' @inheritParams createGiottoCosMxObject
#' @seealso createGiottoCosMxObject load_cosmx_folder
#' @returns path_list a list of cosmx files discovered and their filepaths. NULL
#' values denote missing items
#' @keywords internal
.read_cosmx_folder <- function(
        cosmx_dir,
        verbose = TRUE) {
    ch <- box_chars()

    if (is.null(cosmx_dir) | !dir.exists(cosmx_dir)) {
        stop("The full path to a cosmx directory must be given.")
    }
    vmsg("A structured CosMx directory will be used\n", .v = verbose)

    # find directories (length = 1 if present, length = 0 if missing)
    dir_items <- list(
        `CellLabels folder` = "*CellLabels",
        `CompartmentLabels folder` = "*CompartmentLabels",
        `CellComposite folder` = "*CellComposite",
        `CellOverlay folder` = "*CellOverlay",
        `transcript locations file` = "*tx_file*",
        `fov positions file` = "*fov_positions_file*",
        `expression matrix file` = "*exprMat_file*",
        `metadata file` = "*metadata_file*"
    )
    dir_items <- lapply(
        dir_items, function(x) Sys.glob(paths = file.path(cosmx_dir, x))
    )
    dir_items_lengths <- lengths(dir_items)

    if (isTRUE(verbose)) {
        message("Checking directory contents...")
        for (item in names(dir_items)) {
            if (dir_items_lengths[[item]] > 0) {
                message(ch$s, "> ", item, " found")
            } else {
                warning(item, " is missing\n")
            }
        }
    }

    # select first directory in list if multiple are detected
    if (any(dir_items_lengths > 1)) {
        warning("Multiple matches for expected subdirectory item(s).\n
                First matching item selected")

        multiples <- which(dir_items_lengths > 1)
        for (mult_i in multiples) {
            message(names(dir_items)[[mult_i]], "multiple matches found:")
            print(dir_items[[mult_i]])
            dir_items[[mult_i]] <- dir_items[[mult_i]][[1]]
        }
    }
    vmsg("Directory check done", .v = verbose)

    return(dir_items)
}

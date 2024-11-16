# modular reader functions layout #
# <class & methods> #
# - initialize method for reader object
#   - filepath detection based on directory path
#   - register modular load functions as tx_fun(), poly_fun, etc. in the object.
#     include expected defaults that update based on filepath info
#   - `create_gobject()` single function that utilizes other functions
#     registered to the reader object in previous step. Returns gobject with
#     desired data contents
#     For params that `create_gobject()` is in charge of, the registered funs
#     should use the params passed to `create_gobject()` instead of the baked
#     in defaults
#
# <importXYZ> #
# - exported function to create a `XYZReader` class object
#
# <modular load functions> #
# - from `path` and other minimal args, create a giotto subobject with access
#   to specific ways to load and manipulate data





# CLASS ####



setClass(
    "XeniumReader",
    slots = list(
        xenium_dir = "character",
        filetype = "list",
        qv = "ANY",
        micron = "numeric",
        calls = "list"
    ),
    prototype = list(
        filetype = list(
            transcripts = "parquet",
            boundaries = "parquet",
            expression = "h5",
            cell_meta = "parquet"
        ),
        qv = 20,
        calls = list()
    )
)

# * show ####
setMethod("show", signature("XeniumReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "XeniumReader"))
    print_slots <- c("dir", "filetype", "qv_cutoff", "funs")
    pre <- sprintf(
        "%s :", format(print_slots)
    )
    names(pre) <- print_slots

    # dir
    d <- object@xenium_dir
    if (length(d) > 0L) {
        d <- abbrev_path(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }

    # qv
    qv <- object@qv
    cat(pre["qv_cutoff"], paste(qv, collapse = ", "), "\n")

    # filetype
    .filetype_prints(x = object, pre = pre["filetype"])

    # funs
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("XeniumReader"), function(x, ...) show(x))

# * init ####
setMethod(
    "initialize", signature("XeniumReader"),
    function(
        .Object,
        xenium_dir,
        filetype,
        qv_cutoff,
        micron) {
        obj <- callNextMethod(.Object)

        # provided params (if any)
        if (!missing(xenium_dir)) {
            checkmate::assert_directory_exists(xenium_dir)
            obj@xenium_dir <- xenium_dir
        }
        if (!missing(filetype)) {
            obj@filetype <- filetype
        }
        if (!missing(qv_cutoff)) {
            obj@qv <- qv_cutoff
        }
        if (!missing(micron)) {
            obj@micron <- micron
        }


        # check filetype
        ftype_data <- c("transcripts", "boundaries", "expression", "cell_meta")
        if (!all(ftype_data %in% names(obj@filetype))) {
            stop(wrap_txt(
                "`$filetype` must have entries for each of:\n",
                paste(ftype_data, collapse = ", ")
            ))
        }

        ftype <- obj@filetype
        ft_tab <- c("csv", "parquet")
        ft_exp <- c("h5", "mtx") # zarr not yet supported
        if (!ftype$transcripts %in% ft_tab) {
            stop(
                wrap_txt(
                    "`$filetype$transcripts` must be one of",
                    paste(ft_tab, collapse = ", ")
                ),
                call. = FALSE
            )
        }
        if (!ftype$boundaries %in% ft_tab) {
            stop(
                wrap_txt(
                    "`$filetype$boundaries` must be one of",
                    paste(ft_tab, collapse = ", ")
                ),
                call. = FALSE
            )
        }
        if (!ftype$cell_meta %in% ft_tab) {
            stop(
                wrap_txt(
                    "`$filetype$cell_meta` must be one of",
                    paste(ft_tab, collapse = ", ")
                ),
                call. = FALSE
            )
        }
        if (!ftype$expression %in% ft_exp) {
            stop(
                wrap_txt(
                    "`$filetype$expression` must be one of",
                    paste(tf_exp, collapse = ", ")
                ),
                call. = FALSE
            )
        }


        # detect paths and subdirs
        p <- obj@xenium_dir
        .xenium_detect <- function(pattern, ...) {
            .detect_in_dir(
                pattern = pattern, ...,
                path = p, platform = "Xenium",
            )
        }

        cell_meta_path <- .xenium_detect("cells", first = FALSE)
        panel_meta_path <- .xenium_detect("panel") # json
        experiment_info_path <- .xenium_detect(".xenium") # json

        # 3D stack - DAPI
        img_path <- .xenium_detect("morphology.", warn = FALSE)
        # 2D fusion images
        # - DAPI
        # - stainings for multimodal segmentation
        img_focus_path <- .xenium_detect("morphology_focus", warn = FALSE)
        # Maximum intensity projection (MIP) of the morphology image.
        # (Xenium Outputs v1.0 - 1.9. only)
        img_mip_path <- .xenium_detect("morphology_mip", warn = FALSE)

        tx_path <- .xenium_detect("transcripts", first = FALSE)
        cell_bound_path <- .xenium_detect("cell_bound", first = FALSE)
        nuc_bound_path <- .xenium_detect("nucleus_bound", first = FALSE)

        expr_path <- .xenium_detect("cell_feature_matrix", first = FALSE)

        .xenium_ftype <- function(paths, ftype) {
            paths[grepl(pattern = paste0(".", ftype), x = paths)]
        }


        # select file formats based on reader settings
        tx_path <- .xenium_ftype(tx_path, ftype$transcripts)
        cell_bound_path <- .xenium_ftype(cell_bound_path, ftype$boundaries)
        nuc_bound_path <- .xenium_ftype(nuc_bound_path, ftype$boundaries)
        cell_meta_path <- .xenium_ftype(cell_meta_path, ftype$cell_meta)

        # for mtx, check if directory instead
        if (ftype$expression == "mtx") {
            is_dir <- vapply(
                expr_path, checkmate::test_directory,
                FUN.VALUE = logical(1L)
            )
            expr_path <- expr_path[is_dir]
        } else {
            expr_path <- .xenium_ftype(expr_path, ftype$expression)
        }

        # decide micron scaling
        if (length(obj@micron) == 0) { # if no value already set
            if (!is.null(experiment_info_path)) {
                obj@micron <- fromJSON(
                    experiment_info_path
                )$pixel_size
            } else {
                warning(wrap_txt("No .xenium file found.
                        Guessing 0.2125 as micron scaling"))
                obj@micron <- 0.2125 # default
            }
        }

        # transcripts load call
        tx_fun <- function(
        path = tx_path,
        feat_type = c(
            "rna",
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        split_keyword = list(
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        flip_vertical = TRUE,
        dropcols = c(),
        qv_threshold = obj@qv,
        cores = determine_cores(),
        verbose = NULL) {
            .xenium_transcript(
                path = path,
                feat_type = feat_type,
                split_keyword = split_keyword,
                flip_vertical = flip_vertical,
                dropcols = dropcols,
                qv_threshold = qv_threshold,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_transcripts <- tx_fun

        # load polys call
        poly_fun <- function(
        path = cell_bound_path,
        name = "cell",
        flip_vertical = TRUE,
        calc_centroids = TRUE,
        cores = determine_cores(),
        verbose = NULL) {
            .xenium_poly(
                path = path,
                name = name,
                flip_vertical = flip_vertical,
                calc_centroids = calc_centroids,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_polys <- poly_fun

        # load cellmeta
        cmeta_fun <- function(
        path = cell_meta_path,
        dropcols = c("x_centroid", "y_centroid"),
        cores = determine_cores(),
        verbose = NULL) {
            .xenium_cellmeta(
                path = path,
                dropcols = dropcols,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_cellmeta <- cmeta_fun

        # load featmeta
        fmeta_fun <- function(
        path = panel_meta_path,
        gene_ids = "symbols",
        dropcols = c(),
        cores = determine_cores(),
        verbose = NULL) {
            .xenium_featmeta(
                path = path,
                gene_ids = gene_ids,
                dropcols = dropcols,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_featmeta <- fmeta_fun

        # load expression call
        expr_fun <- function(
        path = expr_path,
        gene_ids = "symbols",
        remove_zero_rows = TRUE,
        split_by_type = TRUE,
        verbose = NULL) {
            .xenium_expression(
                path = path,
                gene_ids = gene_ids,
                remove_zero_rows = remove_zero_rows,
                split_by_type = split_by_type,
                verbose = verbose
            )
        }
        obj@calls$load_expression <- expr_fun

        # load image call
        img_fun <- function(
        path = img_focus_path,
        name = "image",
        output_dir,
        micron = obj@micron,
        negative_y = TRUE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        verbose = NULL,
        ...) {
            .xenium_image(
                path = path,
                name = name,
                output_dir = output_dir,
                micron = micron,
                negative_y = negative_y,
                flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_image <- img_fun

        # load aligned image call
        img_aff_fun <- function(
        path,
        imagealignment_path,
        name = "aligned_image",
        micron = obj@micron,
        verbose = NULL) {
            read10xAffineImage(
                file = path,
                imagealignment_path = imagealignment_path,
                name = name,
                micron = micron,
                verbose = verbose
            )
        }
        obj@calls$load_aligned_image <- img_aff_fun


        # create giotto object call
        gobject_fun <- function(
        transcript_path = tx_path,
        load_bounds = list(
            cell = "cell",
            nucleus = "nucleus"
        ),
        gene_panel_json_path = panel_meta_path,
        expression_path = expr_path,
        metadata_path = cell_meta_path,
        feat_type = c(
            "rna",
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        split_keyword = list(
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        load_images = "focus",
        load_aligned_images = NULL,
        load_transcripts = TRUE,
        load_expression = FALSE,
        load_cellmeta = FALSE,
        instructions = NULL,
        verbose = NULL) {
            load_transcripts <- as.logical(load_transcripts)
            load_expression <- as.logical(load_expression)
            load_cellmeta <- as.logical(load_cellmeta)

            if (!load_transcripts && !load_expression) {
                warning(wrap_txt(
                    "One of either transcripts or expression info should be loaded for a fully functioning object"
                ))
            }

            if (!is.null(load_aligned_images)) {
                checkmate::assert_list(load_aligned_images)
                if (is.null(names(load_aligned_images))) {
                    stop(wrap_txt(
                        "'load_aligned_images' must be a named list"
                    ))
                }
                if (any(lengths(load_aligned_images) != 2L) ||
                    any(!vapply(load_aligned_images, is.character,
                        FUN.VALUE = logical(1L)
                    ))) {
                    stop(wrap_txt(
                        "'load_aligned_images' must be character with length 2:
                        1. image path
                        2. alignment matrix path"
                    ))
                }
            }
            if (!is.null(load_bounds)) {
                checkmate::assert_list(load_bounds)
                if (is.null(names(load_bounds))) {
                    stop("'load_bounds' must be named list of filepaths\n")
                }
            }


            # place calls in new variable for easier access
            funs <- obj@calls

            # init gobject
            g <- giotto(instructions = instructions)


            # transcripts
            if (load_transcripts) {
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    verbose = verbose
                )
                g <- setGiotto(g, tx_list, verbose = FALSE) # lists are fine
            }


            # polys
            if (!is.null(load_bounds)) {
                # replace convenient shortnames
                load_bounds[load_bounds == "cell"] <- cell_bound_path
                load_bounds[load_bounds == "nucleus"] <- nuc_bound_path

                blist <- list()
                bnames <- names(load_bounds)
                for (b_i in seq_along(load_bounds)) {
                    b <- funs$load_polys(
                        path = load_bounds[[b_i]],
                        name = bnames[[b_i]],
                        verbose = verbose
                    )
                    blist <- c(blist, b)
                }
                g <- setGiotto(g, blist, verbose = FALSE)
            }


            # expression
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    gene_ids = "symbols",
                    remove_zero_rows = TRUE,
                    split_by_type = TRUE,
                    verbose = verbose
                )
                g <- setGiotto(g, ex)
            }


            # feat metadata
            fx <- funs$load_featmeta(
                path = gene_panel_json_path,
                # ID = symbols makes sense with the subcellular feat_IDs
                gene_ids = "symbols",
                # no dropcols
                verbose = verbose
            )
            g <- setGiotto(g, fx, verbose = FALSE)


            # cell metadata
            if (load_cellmeta) {
                cx <- funs$load_cellmeta(
                    path = metadata_path,
                    verbose = verbose
                )
                # verbose = FALSE to silence warning that meta is already generated
                g <- setGiotto(g, cx, verbose = FALSE)
            }


            # images
            if (!is.null(load_images)) {
                load_images <- lapply(load_images, normalizePath, mustWork = FALSE)
                img_focus_path <- normalizePath(img_focus_path, mustWork = FALSE)

                # [exception] handle focus image dir
                is_focus <- load_images == "focus" | load_images == img_focus_path
                # split the focus image dir away from other entries
                load_images <- load_images[!is_focus]

                if (any(is_focus)) {
                    focus_dir <- img_focus_path
                    focus_files <- list.files(focus_dir, full.names = TRUE)
                    focus_files <- focus_files[!dir.exists(focus_files)] # ignore matches to export dir
                    nbound <- length(focus_files) - 1L
                    focus_names <- c("dapi", sprintf("bound%d", seq_len(nbound)))
                    names(focus_files) <- focus_names

                    # append to rest of entries
                    load_images <- c(load_images, focus_files)
                }

                # ensure that input is list
                checkmate::assert_list(load_images)
                if (is.null(names(load_images))) {
                    stop("'load_images' must be a named list of filepaths\n")
                }

                imglist <- list()
                imnames <- names(load_images)
                for (impath_i in seq_along(load_images)) {
                    im <- funs$load_image(
                        path = load_images[[impath_i]],
                        name = imnames[[impath_i]]
                    )
                    imglist <- c(imglist, im)
                }
                g <- setGiotto(g, imglist)
            }

            # aligned images can be placed in random places and do not have
            # a standardized naming scheme. Cannot load with expected default.

            if (!is.null(load_aligned_images)) {
                aimglist <- list()
                aimnames <- names(load_aligned_images)
                for (aim_i in seq_along(load_aligned_images)) {
                    vmsg(.v = verbose, sprintf(
                        "loading aligned image as '%s'",
                        aimnames[[aim_i]]
                    ))
                    aim <- funs$load_aligned_image(
                        path = load_aligned_images[[aim_i]][1],
                        imagealignment_path = load_aligned_images[[aim_i]][2],
                        name = aimnames[[aim_i]]
                    )
                    aimglist <- c(aimglist, aim)
                }
                g <- setGiotto(g, aimglist)
            }

            # centroids
            vmsg(.v = verbose, "calculating centroids")
            spat_units_to_calc <- list_spatial_info_names(g)
            g <- addSpatialCentroidLocations(g,
                poly_info = spat_units_to_calc,
                provenance = as.list(spat_units_to_calc),
                verbose = FALSE
            )

            vmsg(.v = verbose, "done")

            return(g)
        }
        obj@calls$create_gobject <- gobject_fun


        return(obj)
    }
)




# access ####

#' @export
setMethod("$", signature("XeniumReader"), function(x, name) {
    basic_info <- c("xenium_dir", "filetype", "qv")
    if (name %in% basic_info) {
        return(methods::slot(x, name))
    }

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("XeniumReader"), function(x, name, value) {
    basic_info <- c("xenium_dir", "filetype", "qv")
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
`.DollarNames.XeniumReader` <- function(x, pattern) {
    dn <- c("xenium_dir", "filetype", "qv")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}




# CREATE READER ####

#' @title Import a 10X Xenium Assay
#' @name importXenium
#' @description
#' Giotto import functionalities for Xenium datasets. This function creates a
#' `XeniumReader` instance that has convenient reader functions for converting
#' individual pieces of Xenium data into Giotto-compatible representations.
#'
#' These functions should have all param values provided as defaults, but
#' can be flexibly modified to do things such as look in alternative
#' directories or paths
#' @param xenium_dir Xenium output directory
#' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included
#' as a subcellular transcript detection (default = 20)
#' @returns `XeniumReader` object
#' @export
importXenium <- function(xenium_dir = NULL, qv_threshold = 20) {
    a <- list(Class = "XeniumReader")
    if (!is.null(xenium_dir)) {
        a$xenium_dir <- xenium_dir
    }
    a$qv <- qv_threshold

    do.call(new, args = a)
}





# MODULAR ####


## transcript ####

.xenium_transcript <- function(
        path,
        feat_type = c(
            "rna",
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        split_keyword = list(
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        flip_vertical = TRUE,
        dropcols = c(),
        qv_threshold = 20,
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to tx file provided or auto-detected"
        ), call. = FALSE)
    }

    checkmate::assert_file_exists(path)
    e <- file_extension(path) %>%
        head(1L) %>%
        tolower()
    vmsg(.v = verbose, .is_debug = TRUE, "[TX_READ] FMT =", e)
    vmsg(.v = verbose, .is_debug = TRUE, path)

    # read in as data.table
    a <- list(
        path = path,
        dropcols = dropcols,
        qv_threshold = qv_threshold,
        verbose = verbose
    )
    vmsg("Loading transcript level info...", .v = verbose)
    # pass to specific reader fun based on filetype
    # return as data.table with colnames `feat_ID`, `x`, `y`
    tx <- switch(e,
        "csv" = do.call(.xenium_transcript_csv,
            args = c(a, list(cores = cores))
        ),
        "parquet" = do.call(.xenium_transcript_parquet, args = a),
        "zarr" = stop("zarr not yet supported")
    )

    # flip values vertically
    y <- NULL # NSE var
    if (flip_vertical) tx[, y := -y]

    # create gpoints
    gpointslist <- createGiottoPoints(
        x = tx,
        feat_type = feat_type,
        split_keyword = split_keyword,
        verbose = FALSE
    )

    if (inherits(gpointslist, "list")) {
        gpointslist <- list(gpointslist)
    }

    return(gpointslist)
}


.xenium_transcript_csv <- function(
        path,
        dropcols = c(),
        qv_threshold = 20,
        cores = determine_cores(),
        verbose = NULL) {
    tx_dt <- data.table::fread(
        path,
        nThread = cores,
        colClasses = c(transcript_id = "character"),
        drop = dropcols
    )
    data.table::setnames(
        x = tx_dt,
        old = c("feature_name", "x_location", "y_location"),
        new = c("feat_ID", "x", "y")
    )

    # qv filtering
    if (!is.null(qv_threshold)) {
        n_before <- tx_dt[, .N]
        tx_dt <- tx_dt[qv >= qv_threshold]
        n_after <- tx_dt[, .N]

        vmsg(
            .v = verbose,
            sprintf(
                "QV cutoff: %d\n Feature points removed: %d, out of %d",
                qv_threshold,
                n_before - n_after,
                n_before
            )
        )
    }

    return(tx_dt)
}

.xenium_transcript_parquet <- function(
        path,
        dropcols = c(),
        qv_threshold = 20,
        verbose = NULL) {
    package_check("dplyr")
    package_check("arrow", custom_msg = sprintf(
        "package 'arrow' is not yet installed\n\n To install:\n%s\n%s%s",
        "Sys.setenv(ARROW_WITH_ZSTD = \"ON\") ",
        "install.packages(\"arrow\", ",
        "repos = c(\"https://apache.r-universe.dev\"))"
    ))

    # setup tx parquet query
    tx_arrow <- arrow::read_parquet(file = path, as_data_frame = FALSE) %>%
        dplyr::mutate(transcript_id = cast(transcript_id, arrow::string())) %>%
        dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
        dplyr::mutate(feature_name = cast(feature_name, arrow::string())) %>%
        dplyr::select(-dplyr::any_of(dropcols))

    # qv filtering
    if (!is.null(qv_threshold)) {
        .nr <- function(x) {
            dplyr::tally(x) %>%
                dplyr::collect() %>%
                as.numeric()
        }
        n_before <- .nr(tx_arrow)
        tx_arrow <- dplyr::filter(tx_arrow, qv > qv_threshold)
        n_after <- .nr(tx_arrow)

        vmsg(.v = verbose, sprintf(
            "QV cutoff: %f\n Feature points removed: %d, out of %d",
            qv_threshold, n_before - n_after, n_before
        ))
    }

    # pull into memory as data.table
    tx_dt <- as.data.frame(tx_arrow) %>% data.table::setDT()
    data.table::setnames(
        x = tx_dt,
        old = c("feature_name", "x_location", "y_location"),
        new = c("feat_ID", "x", "y")
    )
    return(tx_dt)
}


## polygon ####

.xenium_poly <- function(
        path,
        name = "cell",
        flip_vertical = TRUE,
        calc_centroids = TRUE,
        cores = determine_cores(),
        verbose = NULL) {
    checkmate::assert_file_exists(path)
    checkmate::assert_character(name, len = 1L)

    e <- file_extension(path) %>%
        head(1L) %>%
        tolower()

    a <- list(path = path)
    vmsg(sprintf("Loading boundary info '%s'", name), .v = verbose)
    vmsg(.v = verbose, .is_debug = TRUE, "[POLY_READ] FMT =", e)
    vmsg(.v = verbose, .is_debug = TRUE, path)
    # pass to specific load function based on file extension
    # returns as data.table with colnames `cell_id`, `vertex_x`, `vertex_y`
    polys <- switch(e,
        "csv" = do.call(.xenium_poly_csv, args = c(a, list(cores = cores))),
        "parquet" = do.call(.xenium_poly_parquet, args = a),
        "zarr" = stop("zarr not yet supported")
    )

    vertex_y <- NULL # NSE var
    if (flip_vertical) polys[, vertex_y := -vertex_y]

    # create gpolys
    verbose <- verbose %null% FALSE
    gpolys <- createGiottoPolygon(
        x = polys,
        name = name,
        calc_centroids = calc_centroids,
        verbose = verbose
    )
    return(gpolys)
}

.xenium_poly_csv <- function(path, cores = determine_cores()) {
    data.table::fread(
        path,
        nThread = cores,
        colClasses = c(cell_id = "character")
    )
}

.xenium_poly_parquet <- function(path) {
    package_check(
        pkg_name = c("arrow", "dplyr"),
        repository = c("CRAN:arrow", "CRAN:dplyr")
    )
    # read & convert to DT
    arrow::read_parquet(file = path, as_data_frame = FALSE) %>%
        dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
        as.data.frame() %>%
        data.table::setDT()
}


## cellmeta ####

.xenium_cellmeta <- function(
        path,
        dropcols = c(),
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to metadata file provided or auto-detected"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)

    e <- file_extension(path) %>%
        head(1L) %>%
        tolower()
    a <- list(path = path, dropcols = dropcols)
    vmsg("Loading 10X cell metadata...", .v = verbose)
    vmsg(.v = verbose, .is_debug = TRUE, "[CMETA_READ] FMT =", e)
    vmsg(.v = verbose, .is_debug = TRUE, path)
    verbose <- verbose %null% TRUE
    cx <- switch(e,
        "csv" = do.call(
            .xenium_cellmeta_csv,
            args = c(a, list(cores = cores))
        ),
        "parquet" = do.call(.xenium_cellmeta_parquet, args = a)
    )
    data.table::setnames(cx, "cell_id", "cell_ID")

    cx <- createCellMetaObj(
        metadata = cx,
        spat_unit = "cell",
        feat_type = "rna",
        provenance = "cell",
        verbose = verbose
    )
    return(cx)
}

.xenium_cellmeta_csv <- function(path, dropcols = c(),
    cores = determine_cores()) {
    data.table::fread(path, nThread = cores, drop = dropcols)
}

.xenium_cellmeta_parquet <- function(path, dropcols = c()) {
    arrow::read_parquet(file = path, as_data_frame = FALSE) %>%
        dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
        dplyr::select(-dplyr::any_of(dropcols)) %>%
        data.table::as.data.table()
}


## featmeta ####

.xenium_featmeta <- function(
        path,
        gene_ids = "symbols",
        dropcols = c(),
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to panel metadata file provided or auto-detected"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    vmsg("Loading feature metadata...", .v = verbose)
    vmsg(.v = verbose, .is_debug = TRUE, path)
    # updated for pipeline v1.6 json format
    fdata_ext <- GiottoUtils::file_extension(path)
    if ("json" %in% fdata_ext) {
        feat_meta <- .load_xenium_panel_json(
            path = path, gene_ids = gene_ids
        )
    } else {
        feat_meta <- data.table::fread(path, nThread = cores)
        colnames(feat_meta)[[1]] <- "feat_ID"
    }

    dropcols <- dropcols[dropcols %in% colnames(feat_meta)]
    # remove dropcols
    if (length(dropcols) > 0L) feat_meta[, (dropcols) := NULL]

    fx <- createFeatMetaObj(
        metadata = feat_meta,
        spat_unit = "cell",
        feat_type = "rna",
        provenance = "cell",
        verbose = verbose
    )

    return(fx)
}


.load_xenium_panel_json <- function(path, gene_ids = "symbols") {
    gene_ids <- match.arg(gene_ids, c("symbols", "ensembl"))

    # tested on v1.6
    j <- fromJSON(path)
    # j$metadata # dataset meta
    # j$payload # main content
    # j$payload$chemistry # panel chemistry used
    # j$payload$customer # panel customer
    # j$payload$designer # panel designer
    # j$payload$spec_version # versioning
    # j$payload$panel # dataset panel stats

    panel_info <- j$payload$targets$type %>%
        data.table::as.data.table()

    switch(gene_ids,
        "symbols" = data.table::setnames(
            panel_info,
            old = c("data.id", "data.name", "descriptor"),
            new = c("ensembl", "feat_ID", "type")
        ),
        "ensembl" = data.table::setnames(
            panel_info,
            old = c("data.id", "data.name", "descriptor"),
            new = c("feat_ID", "symbol", "type")
        )
    )
    return(panel_info)
}



## expression ####

.xenium_expression <- function(
        path,
        gene_ids = "symbols",
        remove_zero_rows = TRUE,
        split_by_type = TRUE,
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to expression dir (mtx) or file (h5) provided or
            auto-detected"
        ), call. = FALSE)
    }
    if (!file.exists(path)) stop("filepath or directory does not exist.\n")
    a <- list(
        path = path,
        gene_ids = gene_ids,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )

    if (checkmate::test_directory_exists(path)) {
        e <- "mtx" # assume mtx dir
        # zarr can also be unzipped into a dir, but zarr implementation with
        # 32bit UINT support is not available in R yet (needed for cell_IDs).
    } else {
        e <- file_extension(path) %>%
            head(1L) %>%
            tolower()
    }

    vmsg("Loading 10x pre-aggregated expression...", .v = verbose)
    vmsg(.v = verbose, .is_debug = TRUE, "[EXPR_READ] FMT =", e)
    vmsg(.v = verbose, .is_debug = TRUE, path)
    verbose <- verbose %null% TRUE
    ex_list <- switch(e,
        "mtx" = do.call(.xenium_expression_mtx, args = a),
        "h5" = do.call(.xenium_expression_h5, args = a),
        "zarr" = stop("zarr reading not yet implemented")
    )

    # ensure list
    if (!inherits(ex_list, "list")) ex_list <- list(ex_list)
    # set correct feature name
    fname <- "rna"
    if (length(names(ex_list)) > 1L) fname <- names(ex_list)

    # specific feat type naming updates
    fname[fname == "Gene Expression"] <- "rna"
    fname[fname == "Negative Control Codeword"] <- "NegControlCodeword"
    fname[fname == "Negative Control Probe"] <- "NegControlProbe"
    fname[fname == "Blank Codeword"] <- "UnassignedCodeword" # from legacy Xenium pipeline
    fname[fname == "Genomic Control"] <- "GenomicControl"
    fname[fname == "Unassigned Codeword"] <- "UnassignedCodeword"
    fname[fname == "Deprecated Codeword"] <- "DeprecatedCodeword"

    # catch for " " characters in feat type
    # (no major reason for doing this. spaces just make it harder to read)
    fname <- gsub(" ", "_", fname)

    # lapply to process more than one if present
    eo_list <- lapply(seq_along(ex_list), function(ex_i) {
        createExprObj(
            expression_data = ex_list[[ex_i]],
            name = "raw",
            spat_unit = "cell",
            feat_type = fname[[ex_i]],
            provenance = "cell"
        )
    })

    return(eo_list)
}

.xenium_expression_h5 <- function(
        path,
        gene_ids = "symbols",
        remove_zero_rows = TRUE,
        split_by_type = TRUE) {
    get10Xmatrix_h5(
        path_to_data = path,
        gene_ids = gene_ids,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )
}

.xenium_expression_mtx <- function(
        path,
        gene_ids = "symbols",
        remove_zero_rows = TRUE,
        split_by_type = TRUE) {
    gene_ids <- switch(gene_ids,
        "ensembl" = 1,
        "symbols" = 2
    )
    get10Xmatrix(
        path_to_data = path,
        gene_column_index = gene_ids,
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type
    )
}



## image ####


.xenium_image <- function(
        path,
        name,
        output_dir,
        micron,
        negative_y = TRUE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        verbose = NULL,
        ...) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to image file provided or auto-detected"
        ), call. = FALSE)
    }

    # *** whether .ome or not does not matter for this function *** #

    # set default if still missing
    if (missing(name)) name <- "image"
    # default = new tif_exports folder one layer down
    if (missing(output_dir)) output_dir <- "default"

    # [names]
    if (length(name) != length(path) &&
        length(name) != 1) {
        stop("length of `name` should be same as length of `path`")
    }
    if (length(name) == 1 &&
        length(path) > 1) {
        name <- sprintf("%s_%d", name, seq_along(path))
    }

    # [directory input] -> load as individual image paths
    # these need to be expanded then appended to running named list of images
    is_dir <- vapply(path, dir.exists, FUN.VALUE = logical(1L))
    dir_path <- path[is_dir]
    dir_name <- name[is_dir]
    path <- path[!is_dir]
    name <- name[!is_dir]

    # expand directory inputs
    if (length(dir_path) > 0L) {
        for (dir_i in seq_along(dir_path)) {
            dp_i <- dir_path[[dir_i]] # dir path
            dn_i <- dir_name[[dir_i]] # dir name
            vmsg(.is_debug = TRUE, "img dir input:", dp_i)

            # expand and update to per-image
            dfp_i <- list.files(dp_i, full.names = TRUE) # dir file paths
            dfp_i <- dfp_i[!dir.exists(dfp_i)] # ignore dir matches
            # (such as the export directory)
            dfn_i <- sprintf("%s_%d", dn_i, seq_along(dfp_i)) # dir file names
            vmsg(.is_debug = TRUE, "* [img paths]:\n", paste(dfp_i, collapse = "\n"))
            vmsg(.is_debug = TRUE, "* [img names]:\n", paste(dfn_i, collapse = "\n"))

            # append to single file lists
            path <- c(path, dfp_i)
            name <- c(name, dfn_i)
        }
    }

    # [paths]
    # check files exist
    vapply(path, checkmate::assert_file_exists, FUN.VALUE = character(1L))

    # micron
    checkmate::assert_numeric(micron)

    with_pbar({
        p <- pbar(along = path)

        gimg_list <- lapply(seq_along(path), function(img_i) {
            # handle .ome conversion and image subobject creation
            gimg <- .xenium_image_single(
                path = path[[img_i]],
                name = name[[img_i]],
                output_dir = output_dir,
                micron = micron,
                negative_y = negative_y,
                flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal,
                verbose = verbose,
                ...
            )
            p()
            return(gimg)
        })
    })
    return(gimg_list)
}

# per image...
# if .ome : check that converted output path file exists.
#   if exists && if overwrite : remove converted image
#   if still not exist : create converted image
# use image
.xenium_image_single <- function(
        path,
        name = "image",
        output_dir,
        micron,
        negative_y = TRUE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        overwrite = FALSE,
        verbose = NULL) {
    vmsg(.v = verbose, sprintf("loading image as '%s'", name))
    vmsg(.v = verbose, .is_debug = TRUE, path)
    vmsg(
        .v = verbose, .is_debug = TRUE,
        sprintf(
            "negative_y: %s\nflip_vertical: %s\nflip_horizontal: %s",
            negative_y, flip_vertical, flip_horizontal
        ),
        .prefix = ""
    )

    # terra::rast() and gdal still have difficulties with 10x single channel
    # .ome.tif images. May be related to JP2OpenJPEG driver but even loading
    # this does not seem to fix it.
    if ("ome" %in% file_extension(path)) {
        if (output_dir == "default") {
            # default output dir is a new folder under the same directory
            output_dir <- file.path(dirname(path), "tif_exports")
        }

        # check for existence of converted tiff file in output dir
        # fullpath of tiff to write
        tiff_path <- file.path(output_dir, basename(path))
        tiff_path <- gsub(".ome.tif", ".tif", tiff_path)
        if (checkmate::test_file_exists(tiff_path)) {
            vmsg(.is_debug = TRUE, sprintf(
                "converted tiff already present\n%s", tiff_path
            ))
            # if found AND overwrite, remove it to be regenerated downstream
            if (isTRUE(overwrite)) {
                unlink(tiff_path, force = TRUE)
            }
            # the convenience fun can be run multiple times on the dataset
            # So, we allow directly using already converted imgs
        }

        # check the fullpath again
        if (!checkmate::test_file_exists(tiff_path)) {
            vmsg(.is_debug = TRUE, sprintf(
                "converting ome to tif\n%s", tiff_path
            ))
            # if missing, do conversion
            # output is expected at `tiff_path`
            ometif_to_tif(
                input_file = path,
                output_dir = output_dir,
                overwrite = overwrite
            )
        }

        path <- tiff_path
    }

    img <- createGiottoLargeImage(path,
        name = name,
        flip_vertical = flip_vertical,
        flip_horizontal = flip_horizontal,
        negative_y = negative_y,
        verbose = verbose
    )
    img <- rescale(img, micron, x0 = 0, y0 = 0)
    return(img)
}

# for aligned_image (affine), see the `XeniumReader` init method




# wrapper ####


#' @title Create 10x Xenium Giotto Object
#' @name createGiottoXeniumObject
#' @description Create a Giotto object from a Xenium experiment output folder.
#' Only the `xenium_dir`, `load_images`, and `load_aligned_images` params
#' need to be supplied when defaults are sufficient. All other params have
#' defaults set and are there in case of non-standard directory layouts or
#' alternative preference in file format to load from.\cr
#' When possible, `.parquet` files are loaded. This requires the additional
#' installation of \pkg{arrow} with zstd support. See details. `h5` is also
#' used by default if the 10x provided expression matrix is loaded.\cr
#' The 10X provided aggregated expression matrix and cell metdata are not
#' loaded by default since the results may be slightly different from those
#' that Giotto spatially aggregates.
#' @param xenium_dir Full path to the exported xenium directory
#' @param transcript_path Optional. Filepath to desired transcripts file to
#' load. Either the `.parquet` or `.csv` files can be used.
#' @param bounds_path Optional. Named list of filepaths to desired Xenium
#' bounds/polygon files to load. Either the `.parquet` or `.csv` files can be
#' used. The default is to load the `.parquets` of both cell and nucleus.
#' @param gene_panel_json_path Optional. Filepath to panel json. This json
#' contains feature metadata information and ENSG names.
#' @param expression_path Optional. Filepath to cell feature matrix. Accepts
#' either the `.h5` or the unzipped directory containing `.mtx` files.
#' @param cell_metadata_path Optional. Filepath to `cells.csv.gz` or
#' `cells.parquet`
#' which contain cell metadata information.
#' @param feat_type character. feature type. Provide more than one value if
#' using the `split_keyword` param. For each set of keywords to split by, an
#' additional feat_type should be provided in the same order. Affects how
#' the transcripts information is loaded. Helpful for separating out the
#' QC probes. See details.
#' @param split_keyword list of character vectors of keywords to split the
#' transcripts based on their feat_ID. Keywords will be `grepl()`
#' matched against the feature IDs information. See details.
#' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included
#' as a subcellular transcript detection (default = 20)
#' @param load_images Named list of filepaths to `.tif` images, usually the
#' ones in the `morphology_focus` directory. These `ome.tif` images are not
#' compatible and must be converted to `tif` using
#' `[GiottoClass::ometif_to_tif()]`.
#' @param load_aligned_images Named list of filepaths. The list names are used
#' as the image names when loaded. Two filepaths are expected per entry. The
#' first one should be to the `.tif` image. The second path is to the `.csv`
#' alignment matrix file. `ome.tif` images will work, but they are currently
#' slower in our imaging pipeline.
#' @param load_expression logical. Default = FALSE. Whether to load in 10X
#' provided expression matrix.
#' @param load_cellmeta logical. Default = FALSE. Whether to load in 10X
#' provided cell metadata information
#' @param instructions list of instructions or output result from
#' [createGiottoInstructions()]
#' @param verbose logical or NULL. NULL uses the `giotto.verbose` option
#' setting and defaults to TRUE.
#' @returns `giotto` object
#' @details
#'
#' \[\strong{arrow zstd support}\]
#' Xenium parquets have zstd compression. \pkg{arrow} is used to access
#' parquets, however it may not install on all systems with zstd by default.
#' You can check whether zstd support is installed by running:
#' `arrow::arrow_info()$capabilities[["zstd"]]`. If `FALSE`, it needs to be
#' reinstalled with the following:
#' \preformatted{
#'  Sys.setenv(ARROW_WITH_ZSTD = "ON")
#'  install.packages("arrow", repos = c("https://apache.r-universe.dev"))
#' }
#'
#' \[\strong{QC feature types}\]
#' Xenium provides info on feature detections that include more than only the
#' Gene Expression specific probes. Additional probes for QC are included:
#' \emph{blank codeword}, \emph{negative control codeword}, and
#' \emph{negative control probe}. These additional QC probes each occupy and
#' are treated as their own feature types so that they can largely remain
#' independent of the gene expression information.
#'
#' \[\strong{feat_type and split_keyword}\]
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
#' @md
#' @export
createGiottoXeniumObject <- function(
        xenium_dir,
        transcript_path = NULL, # optional
        bounds_path = list( # looks for parquets by default
            cell = "cell",
            nucleus = "nucleus"
        ),
        gene_panel_json_path = NULL, # optional
        expression_path = NULL, # optional
        cell_metadata_path = NULL, # optional
        feat_type = c(
            "rna",
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        split_keyword = list(
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        qv_threshold = 20,
        load_images = "focus",
        load_aligned_images = NULL,
        load_transcripts = TRUE,
        load_expression = FALSE,
        load_cellmeta = FALSE,
        instructions = NULL,
        verbose = NULL) {
    x <- importXenium(xenium_dir)
    # apply reader params
    x$qv <- qv_threshold

    # directly passed
    a <- list(
        load_bounds = bounds_path,
        feat_type = feat_type,
        split_keyword = split_keyword,
        load_images = load_images,
        load_aligned_images = load_aligned_images,
        load_transcripts = load_transcripts,
        load_expression = load_expression,
        load_cellmeta = load_cellmeta,
        instructions = instructions,
        verbose = verbose
    )

    # only passed if not null
    if (!is.null(transcript_path)) a$transcript_path <- transcript_path
    if (!is.null(gene_panel_json_path)) {
        a$gene_panel_json_path <- gene_panel_json_path
    }
    if (!is.null(expression_path)) a$expression_path <- expression_path
    if (!is.null(cell_metadata_path)) a$metadata_path <- cell_metadata_path

    g <- do.call(x$create_gobject, args = a)
    return(g)
}




#' #' @title Create 10x Xenium Giotto Object
#' #' @name createGiottoXeniumObject
#' #' @description Given the path to a Xenium experiment output folder, creates a
#' #' Giotto object
#' #' @param xenium_dir full path to the exported xenium directory
#' #' @param data_to_use which type(s) of expression data to build the gobject with
#' #' (e.g. default: \strong{'subcellular'}, 'aggregate', or 'all')
#' #' @param load_format files formats from which to load the data. Either `csv` or
#' #' `parquet` currently supported.
#' #' @param h5_expression (boolean) whether to load cell_feature_matrix from .h5
#' #' file. Default is \code{TRUE}
#' #' @param h5_gene_ids use gene symbols (default) or ensembl ids for the .h5 gene
#' #' expression matrix
#' #' @param bounds_to_load vector of boundary information to load
#' #' (e.g. \code{'cell'}
#' #' or \code{'nucleus'} by themselves or \code{c('cell', 'nucleus')} to load both
#' #' at the same time.)
#' #' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included
#' #' as a subcellular transcript detection (default = 20)
#' #' @param key_list (advanced) list of grep-based keywords to split the
#' #' subcellular feature detections by feature type. See details
#' #' @inheritParams get10Xmatrix
#' #' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' #' @returns giotto object
#' #' @details
#' #'
#' #' \[\strong{QC feature types}\]
#' #' Xenium provides info on feature detections that include more than only the
#' #' Gene Expression specific probes. Additional probes for QC are included:
#' #' \emph{blank codeword}, \emph{negative control codeword}, and
#' #' \emph{negative control probe}. These additional QC probes each occupy and
#' #' are treated as their own feature types so that they can largely remain
#' #' independent of the gene expression information.
#' #'
#' #' \[\strong{key_list}\]
#' #' Related to \code{data_to_use = 'subcellular'} workflow only:
#' #' Additional QC probe information is in the subcellular feature detections
#' #' information and must be separated from the gene expression information
#' #' during processing.
#' #' The QC probes have prefixes that allow them to be selected from the rest of
#' #' the feature IDs.
#' #' Giotto uses a named list of keywords (\code{key_list}) to select these QC
#' #' probes, with the list names being the names that will be assigned as the
#' #' feature type of these feature detections. The default list is used when
#' #' \code{key_list} = NULL.
#' #'
#' #' Default list:
#' #' \preformatted{
#' #'  list(blank_code = 'BLANK_',
#' #'       neg_code = 'NegControlCodeword_',
#' #'       neg_probe = c('NegControlProbe_|antisense_'))
#' #' }
#' #'
#' #' The Gene expression subset is accepted as the subset of feat_IDs that do not
#' #' map to any of the keys.
#' #'
#' #' @export
#' createGiottoXeniumObject <- function(
#'         xenium_dir,
#'         data_to_use = c("subcellular", "aggregate"),
#'         load_format = "csv",
#'         h5_expression = TRUE,
#'         h5_gene_ids = c("symbols", "ensembl"),
#'         gene_column_index = 1,
#'         bounds_to_load = c("cell"),
#'         qv_threshold = 20,
#'         key_list = NULL,
#'         instructions = NULL,
#'         cores = NA,
#'         verbose = TRUE
#' ) {
#'     # 0. setup
#'     xenium_dir <- path.expand(xenium_dir)
#'
#'     # Determine data to load
#'     data_to_use <- match.arg(
#'         arg = data_to_use, choices = c("subcellular", "aggregate"))
#'
#'     # Determine load formats
#'     load_format <- "csv" # TODO Remove this and add as param once other options
#'     # are available
#'     load_format <- match.arg(
#'         arg = load_format, choices = c("csv", "parquet", "zarr"))
#'
#'     # set number of cores automatically, but with limit of 10
#'     cores <- determine_cores(cores)
#'     data.table::setDTthreads(threads = cores)
#'
#'     # 1. detect xenium folder and find filepaths to load
#'
#'     # path_list contents:
#'     # tx_path
#'     # bound_paths
#'     # cell_meta_path
#'     # agg_expr_path
#'     # panel_meta_path
#'     path_list <- .read_xenium_folder(
#'         xenium_dir = xenium_dir,
#'         data_to_use = data_to_use,
#'         bounds_to_load = bounds_to_load,
#'         load_format = load_format,
#'         h5_expression = h5_expression,
#'         verbose = verbose
#'     )
#'
#'
#'     # 2. load in data
#'
#'     # data_list contents:
#'     # feat_meta
#'     # tx_dt
#'     # bound_dt_list
#'     # cell_meta
#'     # agg_expr
#'     data_list <- .load_xenium_folder(
#'         path_list = path_list,
#'         load_format = load_format,
#'         data_to_use = data_to_use,
#'         h5_expression = h5_expression,
#'         h5_gene_ids = h5_gene_ids,
#'         gene_column_index = gene_column_index,
#'         cores = cores,
#'         verbose = verbose
#'     )
#'
#'
#'     # TODO load images
#'
#'
#'     # 3. Create giotto objects
#'
#'     if (data_to_use == "subcellular") {
#'         # ** feat type search keys **
#'         if (is.null(key_list)) {
#'             key_list <- list(
#'                 blank_code = "BLANK_",
#'                 neg_code = "NegControlCodeword_",
#'                 neg_probe = c("NegControlProbe_|antisense_")
#'             )
#'         }
#'
#'         # needed:
#'         # feat_meta
#'         # tx_dt
#'         # bound_dt_list
#'         xenium_gobject <- .createGiottoXeniumObject_subcellular(
#'             data_list = data_list,
#'             qv_threshold = qv_threshold,
#'             key_list = key_list,
#'             instructions = instructions,
#'             cores = cores,
#'             verbose = verbose
#'         )
#'     }
#'
#'     if (data_to_use == "aggregate") {
#'         # needed:
#'         # feat_meta
#'         # cell_meta
#'         # agg_expr
#'         # optional?
#'         # tx_dt
#'         # bound_dt_list
#'         xenium_gobject <- .createGiottoXeniumObject_aggregate(
#'             data_list = data_list,
#'             instructions = instructions,
#'             cores = cores,
#'             verbose = verbose
#'         )
#'     }
#'
#'     return(xenium_gobject)
#' }
#'
#'
#'
#'
#' #' @title Create a Xenium Giotto object from subcellular info
#' #' @name .createGiottoXeniumObject_subcellular
#' #' @description Subcellular workflow for createGiottoXeniumObject
#' #' @param data_list list of data loaded by \code{\link{.load_xenium_folder}}
#' #' @param key_list regex-based search keys for feature IDs to allow separation
#' #' into separate giottoPoints objects by feat_type
#' #' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included
#' #' as a subcellular transcript detection (default = 20)
#' #' @inheritParams get10Xmatrix
#' #' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' #' @returns giotto object
#' #' @seealso createGiottoXeniumObject .createGiottoXeniumObject_aggregate
#' #' @keywords internal
#' .createGiottoXeniumObject_subcellular <- function(
#'         data_list,
#'         key_list = NULL,
#'         qv_threshold = 20,
#'         instructions = NULL,
#'         cores = NA,
#'         verbose = TRUE
#' ) {
#'     # data.table vars
#'     qv <- NULL
#'
#'     # Unpack data_list info
#'     feat_meta <- data_list$feat_meta
#'     tx_dt <- data_list$tx_dt
#'     bound_dt_list <- data_list$bound_dt_list
#'
#'     # define for data.table
#'     cell_id <- feat_ID <- feature_name <- NULL
#'
#'     vmsg("Building subcellular giotto object...", .v = verbose)
#'     # Giotto points object
#'     vmsg("> points data prep...", .v = verbose)
#'
#'     # filter by qv_threshold
#'     vmsg("> filtering feature detections for Phred score >= ",
#'          qv_threshold, .v = verbose)
#'     n_before <- tx_dt[, .N]
#'     tx_dt_filtered <- tx_dt[qv >= qv_threshold]
#'     n_after <- tx_dt_filtered[, .N]
#'
#'     if (verbose) {
#'         cat(
#'             "Number of feature points removed: ",
#'             n_before - n_after,
#'             " out of ", n_before, "\n"
#'         )
#'     }
#'
#'     vmsg("> splitting detections by feat_type", .v = verbose)
#'     # discover feat_IDs for each feat_type
#'     all_IDs <- tx_dt_filtered[, unique(feat_ID)]
#'     feat_types_IDs <- lapply(
#'         key_list, function(x) all_IDs[grepl(pattern = x, all_IDs)])
#'     rna <- list("rna" = all_IDs[!all_IDs %in% unlist(feat_types_IDs)])
#'     feat_types_IDs <- append(rna, feat_types_IDs)
#'
#'     # separate detections by feature type
#'     points_list <- lapply(
#'         feat_types_IDs,
#'         function(types) {
#'             tx_dt_filtered[feat_ID %in% types]
#'         }
#'     )
#'
#'     # Giotto polygons object
#'     vmsg("> polygons data prep...", .v = verbose)
#'     polys_list <- lapply(
#'         bound_dt_list,
#'         function(bound_type) {
#'             bound_type[, cell_id := as.character(cell_id)]
#'         }
#'     )
#'
#'     xenium_gobject <- createGiottoObjectSubcellular(
#'         gpoints = points_list,
#'         gpolygons = polys_list,
#'         instructions = instructions,
#'         cores = cores,
#'         verbose = verbose
#'     )
#'
#'     # generate centroids
#'     vmsg("Calculating polygon centroids...", .v = verbose)
#'     xenium_gobject <- addSpatialCentroidLocations(
#'         xenium_gobject,
#'         poly_info = c(names(bound_dt_list)),
#'         provenance = as.list(names(bound_dt_list))
#'     )
#'
#'     return(xenium_gobject)
#' }
#'
#'
#'
#'
#'
#' #' @title Create a Xenium Giotto object from aggregate info
#' #' @name .createGiottoXeniumObject_aggregate
#' #' @description Aggregate workflow for createGiottoXeniumObject
#' #' @param data_list list of data loaded by \code{.load_xenium_folder}
#' #' @inheritParams get10Xmatrix
#' #' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' #' @returns giotto object
#' #' @seealso createGiottoXeniumObject .createGiottoXeniumObject_subcellular
#' #' @keywords internal
#' .createGiottoXeniumObject_aggregate <- function(
#'         data_list,
#'         # include_analysis = FALSE,
#'         instructions = NULL,
#'         cores = NA,
#'         verbose = TRUE
#' ) {
#'     # Unpack data_list info
#'     feat_meta <- data_list$feat_meta
#'     cell_meta <- data_list$cell_meta
#'     agg_expr <- data_list$agg_expr
#'
#'     # define for data.table
#'     cell_ID <- x_centroid <- y_centroid <- NULL
#'
#'     # clean up names for aggregate matrices
#'     names(agg_expr) <- gsub(pattern = " ", replacement = "_", names(agg_expr))
#'     geneExpMat <- which(names(agg_expr) == "Gene_Expression")
#'     names(agg_expr)[[geneExpMat]] <- "raw"
#'
#'     # set cell_id as character
#'     cell_meta <- cell_meta[, data.table::setnames(.SD, "cell_id", "cell_ID")]
#'     cell_meta <- cell_meta[, cell_ID := as.character(cell_ID)]
#'
#'     # set up spatial locations
#'     agg_spatlocs <- cell_meta[, .(x_centroid, y_centroid, cell_ID)]
#'
#'     # set up metadata
#'     agg_meta <- cell_meta[, !c("x_centroid", "y_centroid")]
#'
#'     vmsg("Building aggregate giotto object...", .v = verbose)
#'     xenium_gobject <- createGiottoObject(
#'         expression = agg_expr,
#'         spatial_locs = agg_spatlocs,
#'         instructions = instructions,
#'         cores = cores,
#'         verbose = verbose
#'     )
#'
#'     # append aggregate metadata
#'     xenium_gobject <- addCellMetadata(
#'         gobject = xenium_gobject,
#'         new_metadata = agg_meta,
#'         by_column = TRUE,
#'         column_cell_ID = "cell_ID"
#'     )
#'     xenium_gobject <- addFeatMetadata(
#'         gobject = xenium_gobject,
#'         new_metadata = feat_meta,
#'         by_column = TRUE,
#'         column_feat_ID = "feat_ID"
#'     )
#'
#'     return(xenium_gobject)
#' }
#'
#'
#'
#'
#' #' @title Read a structured xenium folder
#' #' @name .read_xenium_folder
#' #' @inheritParams createGiottoXeniumObject
#' #' @keywords internal
#' #' @returns path_list a list of xenium files discovered and their filepaths. NULL
#' #' values denote missing items
#' .read_xenium_folder <- function(
#'         xenium_dir,
#'         data_to_use = "subcellular",
#'         bounds_to_load = c("cell"),
#'         load_format = "csv",
#'         h5_expression = FALSE,
#'         verbose = TRUE
#' ) {
#'     # Check needed packages
#'     if (load_format == "parquet") {
#'         package_check(pkg_name = "arrow", repository = "CRAN")
#'         package_check(pkg_name = "dplyr", repository = "CRAN")
#'     }
#'     if (isTRUE(h5_expression)) {
#'         package_check(pkg_name = "hdf5r", repository = "CRAN")
#'     }
#'
#'     ch <- box_chars()
#'
#'
#'     # 0. test if folder structure exists and is as expected
#'
#'
#'     if (is.null(xenium_dir) | !dir.exists(xenium_dir))
#'         stop("The full path to a xenium directory must be given.")
#'     vmsg("A structured Xenium directory will be used\n", .v = verbose)
#'
#'     # find items (length = 1 if present, length = 0 if missing)
#'     dir_items <- list(
#'         `analysis info` = "*analysis*",
#'         `boundary info` = "*bound*",
#'         `cell feature matrix` = "*cell_feature_matrix*",
#'         `cell metadata` = "*cells*",
#'         `image info` = "*tif",
#'         `panel metadata` = "*panel*",
#'         `raw transcript info` = "*transcripts*",
#'         `experiment info (.xenium)` = "*.xenium"
#'     )
#'
#'     dir_items <- lapply(
#'         dir_items, function(x) Sys.glob(paths = file.path(xenium_dir, x)))
#'     dir_items_lengths <- lengths(dir_items)
#'
#'     if (isTRUE(verbose)) {
#'         message("Checking directory contents...")
#'         for (item in names(dir_items)) {
#'             # IF ITEM FOUND
#'
#'             if (dir_items_lengths[[item]] > 0) {
#'                 message(ch$s, "> ", item, " found")
#'                 for (item_i in seq_along(dir_items[[item]])) {
#'                     # print found item names
#'                     subItem <- gsub(pattern = ".*/", replacement = "",
#'                                     x = dir_items[[item]][[item_i]])
#'                     message(ch$s, ch$s, ch$l, ch$h, ch$h, subItem)
#'                 }
#'             } else {
#'                 # IF ITEM MISSING
#'                 # Based on workflow, determine if:
#'                 # necessary (error)
#'                 # optional (warning)
#'
#'                 if (data_to_use == "subcellular") {
#'                     # necessary items
#'                     if (item %in% c("boundary info", "raw transcript info"))
#'                         stop(item, " is missing")
#'                     # optional items
#'                     if (item %in% c(
#'                         "image info", "experiment info (.xenium)",
#'                         "panel metadata"))
#'                         warning(item, " is missing (optional)")
#'                     # items to ignore: analysis info, cell feature matrix,
#'                     # cell metadata
#'                 } else if (data_to_use == "aggregate") {
#'                     # necessary items
#'                     if (item %in% c("cell feature matrix", "cell metadata"))
#'                         stop(item, " is missing")
#'                     # optional items
#'                     if (item %in% c(
#'                         "image info", "experiment info (.xenium)",
#'                         "panel metadata", "analysis info"))
#'                         warning(item, " is missing (optional)")
#'                     # items to ignore: boundary info, raw transcript info
#'                 }
#'             }
#'         }
#'     }
#'
#'
#'     # 1. Select data to load
#'
#'
#'     # **** transcript info ****
#'     tx_path <- NULL
#'     tx_path <- dir_items$`raw transcript info`[grepl(
#'         pattern = load_format, dir_items$`raw transcript info`)]
#'     # **** cell metadata ****
#'     cell_meta_path <- NULL
#'     cell_meta_path <- dir_items$`cell metadata`[grepl(
#'         pattern = load_format, dir_items$`cell metadata`)]
#'
#'     # **** boundary info ****
#'     # Select bound load format
#'     if (load_format != "zarr") { # No zarr available for boundary info
#'         dir_items$`boundary info` <- dir_items$`boundary info`[grepl(
#'             pattern = load_format, dir_items$`boundary info`)]
#'     } else {
#'         dir_items$`boundary info` <- dir_items$`boundary info`[grepl(
#'             pattern = "csv", dir_items$`boundary info`)]
#'     }
#'
#'     # Organize bound paths by type of bound (bounds_to_load param)
#'     bound_paths <- NULL
#'     bound_names <- bounds_to_load
#'     bounds_to_load <- as.list(bounds_to_load)
#'     bound_paths <- lapply(bounds_to_load, function(x) dir_items$`boundary info`[
#'         grepl(pattern = x, dir_items$`boundary info`)])
#'     names(bound_paths) <- bound_names
#'
#'     # **** aggregated expression info ****
#'     agg_expr_path <- NULL
#'     if (isTRUE(h5_expression)) { # h5 expression matrix loading is default
#'         agg_expr_path <- dir_items$`cell feature matrix`[grepl(
#'             pattern = "h5", dir_items$`cell feature matrix`)]
#'     } else if (load_format == "zarr") {
#'         agg_expr_path <- dir_items$`cell feature matrix`[grepl(
#'             pattern = "zarr", dir_items$`cell feature matrix`)]
#'     } else { # No parquet for aggregated expression - default to normal 10x loading
#'         agg_expr_path <- dir_items$`cell feature matrix`[sapply(
#'             dir_items$`cell feature matrix`, function(x) file_test(op = "-d", x))]
#'         if (length(agg_expr_path) == 0) {
#'             stop(wrap_txt(
#'                 "Expression matrix cannot be loaded.\n
#'                 Has cell_feature_matrix(.tar.gz) been unpacked into a
#'                 directory?"
#'             ))
#'         }
#'     }
#'     if (data_to_use == "aggregate") {
#'         if (length(path_list$agg_expr_path) == 0) {
#'             stop(wrap_txt(
#'                 "Aggregated expression not found.\n
#'                 Please confirm h5_expression and load_format params are correct"
#'             ))
#'         }
#'     }
#'
#'     # **** panel info ****
#'     panel_meta_path <- NULL
#'     panel_meta_path <- dir_items$`panel metadata`
#'
#'
#'     vmsg("Directory check done", .v = verbose)
#'
#'     path_list <- list(
#'         "tx_path" = tx_path,
#'         "bound_paths" = bound_paths,
#'         "cell_meta_path" = cell_meta_path,
#'         "agg_expr_path" = agg_expr_path,
#'         "panel_meta_path" = panel_meta_path
#'     )
#'
#'     return(path_list)
#' }
#'
#' #' @title Load xenium data from folder
#' #' @name load_xenium_folder
#' #' @param path_list list of full filepaths from .read_xenium_folder
#' #' @inheritParams createGiottoXeniumObject
#' #' @returns list of loaded in xenium data
#' NULL
#'
#' #' @rdname load_xenium_folder
#' #' @keywords internal
#' .load_xenium_folder <- function(
#'         path_list,
#'         load_format = "csv",
#'         data_to_use = "subcellular",
#'         h5_expression = "FALSE",
#'         h5_gene_ids = "symbols",
#'         gene_column_index = 1,
#'         cores,
#'         verbose = TRUE
#' ) {
#'     data_list <- switch(load_format,
#'                         "csv" = .load_xenium_folder_csv(
#'                             path_list = path_list,
#'                             data_to_use = data_to_use,
#'                             h5_expression = h5_expression,
#'                             h5_gene_ids = h5_gene_ids,
#'                             gene_column_index = gene_column_index,
#'                             cores = cores,
#'                             verbose = verbose
#'                         ),
#'                         "parquet" = .load_xenium_folder_parquet(
#'                             path_list = path_list,
#'                             data_to_use = data_to_use,
#'                             h5_expression = h5_expression,
#'                             h5_gene_ids = h5_gene_ids,
#'                             gene_column_index = gene_column_index,
#'                             cores = cores,
#'                             verbose = verbose
#'                         ),
#'                         "zarr" = stop("load_format zarr:\n Not yet implemented", call. = FALSE)
#'     )
#'
#'     return(data_list)
#' }
#'
#'
#' #' @describeIn load_xenium_folder Load from csv files
#' #' @keywords internal
#' .load_xenium_folder_csv <- function(
#'         path_list,
#'         cores,
#'         data_to_use = "subcellular",
#'         h5_expression = FALSE,
#'         h5_gene_ids = "symbols",
#'         gene_column_index = 1,
#'         verbose = TRUE
#' ) {
#'     # initialize return vars
#'     feat_meta <- tx_dt <- bound_dt_list <- cell_meta <- agg_expr <- NULL
#'
#'     vmsg("Loading feature metadata...", .v = verbose)
#'     # updated for pipeline v1.6 json format
#'     fdata_path <- path_list$panel_meta_path[[1]]
#'     fdata_ext <- GiottoUtils::file_extension(fdata_path)
#'     if ("json" %in% fdata_ext) {
#'         feat_meta <- .load_xenium_panel_json(path = fdata_path,
#'                                              gene_ids = h5_gene_ids)
#'     } else {
#'         feat_meta <- data.table::fread(fdata_path, nThread = cores)
#'         colnames(feat_meta)[[1]] <- "feat_ID"
#'     }
#'
#'     # **** subcellular info ****
#'     if (data_to_use == "subcellular") {
#'         # append missing QC probe info to feat_meta
#'         if (isTRUE(h5_expression)) {
#'             h5 <- hdf5r::H5File$new(path_list$agg_expr_path)
#'             tryCatch({
#'                 root <- names(h5)
#'                 feature_id <- h5[[paste0(root, "/features/id")]][]
#'                 feature_info <- h5[[paste0(root, "/features/feature_type")]][]
#'                 feature_names <- h5[[paste0(root, "/features/name")]][]
#'                 features_dt <- data.table::data.table(
#'                     "id" = feature_id,
#'                     "name" = feature_names,
#'                     "feature_type" = feature_info
#'                 )
#'             }, finally = {
#'                 h5$close_all()
#'             })
#'         } else {
#'             features_dt <- data.table::fread(
#'                 paste0(path_list$agg_expr_path, "/features.tsv.gz"),
#'                 header = FALSE
#'             )
#'         }
#'         colnames(features_dt) <- c("id", "feat_ID", "feat_class")
#'         feat_meta <- merge(
#'             features_dt[, c(2, 3)], feat_meta, all.x = TRUE, by = "feat_ID")
#'
#'         GiottoUtils::vmsg("Loading transcript level info...", .v = verbose)
#'         tx_dt <- data.table::fread(path_list$tx_path[[1]], nThread = cores)
#'         data.table::setnames(
#'             x = tx_dt,
#'             old = c("feature_name", "x_location", "y_location"),
#'             new = c("feat_ID", "x", "y")
#'         )
#'
#'         GiottoUtils::vmsg("Loading boundary info...", .v = verbose)
#'         bound_dt_list <- lapply(
#'             path_list$bound_paths,
#'             function(x) data.table::fread(x[[1]], nThread = cores)
#'         )
#'     }
#'
#'     # **** aggregate info ****
#'     GiottoUtils::vmsg("loading cell metadata...", .v = verbose)
#'     cell_meta <- data.table::fread(
#'         path_list$cell_meta_path[[1]], nThread = cores)
#'
#'     if (data_to_use == "aggregate") {
#'         GiottoUtils::vmsg("Loading aggregated expression...", .v = verbose)
#'         if (isTRUE(h5_expression)) {
#'             agg_expr <- get10Xmatrix_h5(
#'                 path_to_data = path_list$agg_expr_path,
#'                 gene_ids = h5_gene_ids,
#'                 remove_zero_rows = TRUE,
#'                 split_by_type = TRUE
#'             )
#'         } else {
#'             agg_expr <- get10Xmatrix(
#'                 path_to_data = path_list$agg_expr_path,
#'                 gene_column_index = gene_column_index,
#'                 remove_zero_rows = TRUE,
#'                 split_by_type = TRUE
#'             )
#'         }
#'     }
#'
#'     data_list <- list(
#'         "feat_meta" = feat_meta,
#'         "tx_dt" = tx_dt,
#'         "bound_dt_list" = bound_dt_list,
#'         "cell_meta" = cell_meta,
#'         "agg_expr" = agg_expr
#'     )
#'
#'     return(data_list)
#' }
#'
#'
#'
#'
#' #' @describeIn load_xenium_folder Load from parquet files
#' #' @keywords internal
#' .load_xenium_folder_parquet <- function(
#'         path_list,
#'         cores,
#'         data_to_use = "subcellular",
#'         h5_expression = FALSE,
#'         h5_gene_ids = "symbols",
#'         gene_column_index = 1,
#'         verbose = TRUE
#' ) {
#'     # initialize return vars
#'     feat_meta <- tx_dt <- bound_dt_list <- cell_meta <- agg_expr <- NULL
#'     # dplyr variable
#'     cell_id <- NULL
#'
#'     vmsg("Loading feature metadata...", .v = verbose)
#'     # updated for pipeline v1.6 json format
#'     fdata_path <- path_list$panel_meta_path[[1]]
#'     fdata_ext <- GiottoUtils::file_extension(fdata_path)
#'     if ("json" %in% fdata_ext) {
#'         feat_meta <- .load_xenium_panel_json(
#'             path = fdata_path, gene_ids = h5_gene_ids)
#'     } else {
#'         feat_meta <- data.table::fread(fdata_path, nThread = cores)
#'         colnames(feat_meta)[[1]] <- "feat_ID"
#'     }
#'
#'     # **** subcellular info ****
#'     if (data_to_use == "subcellular") {
#'         # define for data.table
#'         transcript_id <- feature_name <- NULL
#'
#'         # append missing QC probe info to feat_meta
#'         if (isTRUE(h5_expression)) {
#'             h5 <- hdf5r::H5File$new(path_list$agg_expr_path)
#'             tryCatch({
#'                 root <- names(h5)
#'                 feature_id <- h5[[paste0(root, "/features/id")]][]
#'                 feature_info <- h5[[paste0(root, "/features/feature_type")]][]
#'                 feature_names <- h5[[paste0(root, "/features/name")]][]
#'                 features_dt <- data.table::data.table(
#'                     "id" = feature_id,
#'                     "name" = feature_names,
#'                     "feature_type" = feature_info
#'                 )
#'             }, finally = {
#'                 h5$close_all()
#'             })
#'         } else {
#'             features_dt <- arrow::read_tsv_arrow(paste0(
#'                 path_list$agg_expr_path, "/features.tsv.gz"),
#'                 col_names = FALSE
#'             ) %>%
#'                 data.table::setDT()
#'         }
#'         colnames(features_dt) <- c("id", "feat_ID", "feat_class")
#'         feat_meta <- merge(features_dt[
#'             , c(2, 3)], feat_meta, all.x = TRUE, by = "feat_ID")
#'
#'         vmsg("Loading transcript level info...", .v = verbose)
#'         tx_dt <- arrow::read_parquet(
#'             file = path_list$tx_path[[1]],
#'             as_data_frame = FALSE
#'         ) %>%
#'             dplyr::mutate(
#'                 transcript_id = cast(transcript_id, arrow::string())) %>%
#'             dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
#'             dplyr::mutate(
#'                 feature_name = cast(feature_name, arrow::string())) %>%
#'             as.data.frame() %>%
#'             data.table::setDT()
#'         data.table::setnames(
#'             x = tx_dt,
#'             old = c("feature_name", "x_location", "y_location"),
#'             new = c("feat_ID", "x", "y")
#'         )
#'         vmsg("Loading boundary info...", .v = verbose)
#'         bound_dt_list <- lapply(path_list$bound_paths, function(x) {
#'             arrow::read_parquet(file = x[[1]], as_data_frame = FALSE) %>%
#'                 dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
#'                 as.data.frame() %>%
#'                 data.table::setDT()
#'         })
#'     }
#'     # **** aggregate info ****
#'     if (data_to_use == "aggregate") {
#'         vmsg("Loading cell metadata...", .v = verbose)
#'         cell_meta <- arrow::read_parquet(
#'             file = path_list$cell_meta_path[[1]],
#'             as_data_frame = FALSE
#'         ) %>%
#'             dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
#'             as.data.frame() %>%
#'             data.table::setDT()
#'
#'         # NOTE: no parquet for agg_expr.
#'         vmsg("Loading aggregated expression...", .v = verbose)
#'         if (isTRUE(h5_expression)) {
#'             agg_expr <- get10Xmatrix_h5(
#'                 path_to_data = path_list$agg_expr_path,
#'                 gene_ids = h5_gene_ids,
#'                 remove_zero_rows = TRUE,
#'                 split_by_type = TRUE
#'             )
#'         } else {
#'             agg_expr <- get10Xmatrix(
#'                 path_to_data = path_list$agg_expr_path,
#'                 gene_column_index = gene_column_index,
#'                 remove_zero_rows = TRUE,
#'                 split_by_type = TRUE
#'             )
#'         }
#'     }
#'
#'     data_list <- list(
#'         "feat_meta" = feat_meta,
#'         "tx_dt" = tx_dt,
#'         "bound_dt_list" = bound_dt_list,
#'         "cell_meta" = cell_meta,
#'         "agg_expr" = agg_expr
#'     )
#'
#'     return(data_list)
#' }

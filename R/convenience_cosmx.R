# CLASS ####


setClass(
    "CosmxReader",
    slots = list(
        cosmx_dir = "character",
        protein_dir = "character",
        version = "character",
        slide = "numeric",
        fovs = "numeric",
        micron = "logical",
        px2um = "numeric",
        poly_pref = "character",
        offsets = "ANY",
        calls = "list",
        paths = "list"
    ),
    prototype = list(
        version = "default",
        slide = 1,
        micron = FALSE,
        px2um = 0.12028, # from cosmx output help files
        poly_pref = "mask",
        offsets = NULL,
        calls = list(),
        paths = list()
    )
)

# * show ####
setMethod("show", signature("CosmxReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "CosmxReader"))
    print_slots <- c("version", "dir", "slide", "fovs",
                    "micron", "poly_pref", "offsets", "funs")
    has_protein <- length(.cosmx_protein_dir(object)) > 0L
    if (has_protein) print_slots <- append(print_slots, "protein", after = 2L)
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
    if (has_protein) {
        cat(pre["protein"],
            GiottoUtils::str_abbreviate(object@protein_dir), "\n")
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
    
    # poly preference
    ppref <- object@poly_pref
    cat(pre["poly_pref"], ppref, "\n")

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
        if (length(x@fovs) == 0) {
            text(y ~ x, data = dat, labels = dat$fov, cex = cex, ...)
        } else {
            fov_sel <- dat$fov %in% x@fovs
            dat_grey <- dat[!fov_sel]
            dat_black <- dat[fov_sel]
            text(y ~ x, data = dat_grey, labels = dat_grey$fov, 
                cex = cex, col = "grey", ...)
            text(y ~ x, data = dat_black, labels = dat_black$fov, 
                cex = cex, col = "black", ...)
        }
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
#' based on `CosMx-ReadMe.html` info. May be different depending on dataset.
#' @param poly_pref character. Either "csv" (default) or "mask". Which format of
#' data to load as polygon info. "csv" will use vector polygons from the
#' `polygons.csv`. "mask" will load the mask images from `CellLabels` directory.
#' @param version character. Version of CosMx output. Current selections are
#' either "default", "v6", and "legacy" (for the NSCLC dataset).
#' @param backend (optional) a `gsource`-inheriting project backend (typically
#' produced by `GiottoDisk::sourceCreate()`). When provided, creates the
#' `giotto` object as a managed on-disk project.
#' @param protein_dir (optional) directory of the protein flat files of a
#' CosMx RNA + protein (multiomic) run: `exprMat_file`, `metadata_file` and
#' the `plex-*.txt` panel file. Kept apart from `cosmx_dir` because the two
#' exports use the same file names. Adds a `load_protein()` function and the
#' protein arguments of `create_gobject()`. With `backend`, it is passed to
#' `GiottoDisk::importCosMxDisk()`.
#' @details
#' Loading functions are generated after the `cosmx_dir` is added.
#' Transcripts, expression, and metadata loading are all expected to be done
#' from the top level of the directory. Loading of polys, and any image sets
#' are expected to be from specific subdirectories containing only those
#' images for the set of FOVs.
#'
#' \[**Protein**\] `load_protein()` reads the protein `exprMat_file`: target
#' antibodies become feature type `"protein"` and the isotype controls
#' (`CodeClass` "Negative" in the plex file; names containing `"IgG"` when
#' there is no plex file) become `"protein_neg"`. Values are the vendor's
#' per-cell mean intensities, stored unchanged as `"raw"`. The `Channel-*`
#' morphology columns are not loaded: they repeat the metadata `Mean.*`
#' columns.
#' @returns `CosmxReader` object, or `CosMxDiskReader` when `backend` is set
#' @examples
#' # Create a `CosmxReader` object
#' reader <- importCosMx()
#'
#' \dontrun{
#' # Set the cosmx_dir and fov parameters
#' path <- "path/to/cosmx/dir"
#' reader$cosmx_dir <- path
#' reader$fov <- c(1, 4)
#'
#' plot(reader) # displays FOVs (top left corner) in px scale.
#'
#' # Load polygons, transcripts, and images
#' polys <- reader$load_polys()
#' tx <- reader$load_transcripts()
#' imgs <- reader$load_images()
#' 
#' # polygons (mask) and images loading supports multiple filepaths
#' # This can be useful when loading AtoMx outputs
#' polys <- reader$load_polys(path = list.files(path,
#'     pattern = "CellLabels_F",
#'     recursive = TRUE,
#'     full.names = TRUE
#' ))
#' imgs <- reader$load_images(path = list.files(path,
#'     pattern = "Composite_F",
#'     recursive = TRUE,
#'     full.names = TRUE
#' ))
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
    micron = FALSE, px2um = 0.12028, poly_pref = "mask",
    version = "default", backend = NULL, protein_dir = NULL) {
    if (!is.null(backend)) {
        package_check(
            "GiottoDisk",
            repository = "github:giotto-suite/GiottoDisk"
        )
        disk_args <- list(
            cosmx_dir = cosmx_dir,
            backend = backend,
            slide = slide,
            fovs = fovs,
            version = version,
            micron = micron,
            px2um = px2um,
            poly_pref = poly_pref
        )
        if (!is.null(protein_dir)) {
            if (!"protein_dir" %in%
                names(formals(GiottoDisk::importCosMxDisk))) {
                stop(wrap_txt(
                    "[importCosMx] the installed GiottoDisk cannot read CosMx
                    protein data (`importCosMxDisk()` has no `protein_dir`).
                    Update GiottoDisk or leave out `protein_dir`."
                ), call. = FALSE)
            }
            disk_args$protein_dir <- protein_dir
        }
        return(do.call(GiottoDisk::importCosMxDisk, disk_args))
    }
    # get params
    a <- list(Class = "CosmxReader")
    if (!is.null(cosmx_dir)) {
        a$cosmx_dir <- cosmx_dir
    }
    if (!is.null(protein_dir)) {
        a$protein_dir <- protein_dir
    }
    if (!is.null(fovs)) {
        a$fovs <- fovs
    }
    a$slide <- slide
    a$micron <- micron
    a$px2um <- px2um
    a$poly_pref <- poly_pref

    do.call(new, args = a)
}

# * init ####
setMethod(
    "initialize", signature("CosmxReader"),
    function(.Object, cosmx_dir, version, slide, fovs, micron, px2um, poly_pref,
             protein_dir) {
        # provided params (if any)
        if (!missing(cosmx_dir)) {
            checkmate::assert_directory_exists(cosmx_dir)
            .Object@cosmx_dir <- cosmx_dir
        }
        if (!missing(protein_dir)) {
            checkmate::assert_directory_exists(protein_dir)
            .cosmx_check_protein_slot(.Object)
            .Object@protein_dir <- protein_dir
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
        if (!missing(poly_pref)) {
            .Object@poly_pref <- poly_pref
        }

        # NULL case
        if (length(.Object@cosmx_dir) == 0) {
            return(.Object) # return early if no path given
        }


        # detect paths and subdirs. Stored on @paths so subclass init
        # frames (e.g. CosMxDiskReader) can list2env() them in for their
        # own closures without re-running detection.
        .Object@paths <- c(
            .cosmx_detect_paths(.Object@cosmx_dir),
            .cosmx_protein_paths(.cosmx_protein_dir(.Object))
        )
        list2env(.Object@paths, envir = environment())

        v <- .Object@version
        if (v == "default") v <- "v6"


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
        col_classes_use <- switch(v,
            "v6" = c(
                "integer", # fov
                "integer", # cell_ID
                "character", # cell
                "integer", # x_local_px
                "integer", # y_local_px
                "double", # x_global_px
                "double", # y_global_px
                "integer", # z
                "character", # target
                "character" # CellComp
            ),
            "legacy" = c(
                "integer", # fov
                "integer", # cell_ID
                "double", # x_global_px
                "double", # y_global_px
                "double", # x_local_px
                "double", # y_local_px
                "integer", # z
                "character", # target
                "character" # CellComp
            )
        )
        tx_fun <- function(
        path = tx_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        col_classes = col_classes_use,
        dropcols = dcols,
        cores = determine_cores(),
        verbose = NULL) {
            .cosmx_transcript(
                path = path,
                fovs = .Object@fovs %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                col_classes = col_classes,
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
            "v6" = -1,
            "legacy" = FALSE
        )
        
        if (!is.null(mask_dir) && !is.null(poly_path)) {
            poly_use_path <- switch(.Object@poly_pref,
                "csv" = poly_path,
                "mask" = mask_dir
            )
        } else if (!is.null(poly_path)) {
            poly_use_path <- poly_path
        } else if (!is.null(mask_dir)) {
            poly_use_path <- mask_dir
        }
        
        poly_fun <- function(
        path = poly_use_path,
        file_fov_fmt = NULL,
        # VERTICAL FLIP + NO VERTICAL SHIFT
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        shift_vertical_step = vstep,
        shift_horizontal_step = FALSE,
        remove_background_polygon = TRUE,
        keep_ids = NULL,
        verbose = NULL) {
            .cosmx_poly(
                path = path,
                fovs = .Object@fovs %none% NULL,
                file_fov_fmt = file_fov_fmt,
                flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal,
                shift_vertical_step = shift_vertical_step,
                shift_horizontal_step = shift_horizontal_step,
                remove_background_polygon = remove_background_polygon,
                micron = .Object@micron,
                px2um = .Object@px2um,
                offsets = .Object@offsets,
                slide = .Object@slide,
                keep_ids = keep_ids,
                verbose = verbose
            )
        }
        .Object@calls$load_polys <- poly_fun


        # expression load call
        expr_fun <- function(
        path = expr_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb")) {
            .cosmx_expression(
                path = path,
                fovs = .Object@fovs %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                slide = .Object@slide
            )
        }
        .Object@calls$load_expression <- expr_fun


        # protein expression load call
        prot_fun <- function(
        path = protein_expr_path,
        feat_type = c("protein", "protein_neg"),
        split_keyword = NULL,
        drop_cols = "^Channel-",
        cores = determine_cores(),
        verbose = NULL) {
            if (is.null(path)) {
                stop(wrap_txt(
                    "[CosmxReader] no protein exprMat_file found.
                    Set `protein_dir` or pass `path`."
                ), call. = FALSE)
            }
            if (is.null(split_keyword)) {
                split_keyword <- if (length(feat_type) > 1L) {
                    .cosmx_protein_negatives(protein_plex_path, path)
                } else {
                    list()
                }
            }
            .cosmx_check_protein_header(.cosmx_csv_header(path),
                split_keyword, basename(path))
            if (length(feat_type) != length(split_keyword) + 1L) {
                stop(wrap_txt(sprintf(
                    "[CosmxReader] %d feat_type(s) given for %d feature
                    groups (the targets plus one per split_keyword).",
                    length(feat_type), length(split_keyword) + 1L
                )), call. = FALSE)
            }
            .cosmx_expression(
                path = path,
                fovs = .Object@fovs %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                slide = .Object@slide,
                drop_cols = drop_cols,
                cores = cores,
                verbose = verbose
            )
        }
        .Object@calls$load_protein <- prot_fun


        # images load call
        negy <- switch(v,
            "v6" = TRUE,
            "legacy" = FALSE,
        )
        
        img_fun <- function(
        path = composite_img_dir,
        file_fov_fmt = NULL,
        img_type = "composite",
        img_name_fmt = paste0(img_type, "_fov%03d"),
        negative_y = negy,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        verbose = NULL) {
            .cosmx_image(
                path = path,
                fovs = .Object@fovs %none% NULL,
                file_fov_fmt = file_fov_fmt,
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
        feat_type = "rna",
        select = NULL,
        rename = NULL,
        verbose = NULL) {
            .cosmx_cellmeta(
                path = path,
                slide = .Object@slide,
                fovs = .Object@fovs %none% NULL,
                dropcols = dropcols,
                feat_type = feat_type,
                select = select,
                rename = rename,
                cores = cores,
                verbose = verbose
            )
        }
        .Object@calls$load_cellmeta <- meta_fun


        # build gobject call
        gobject_fun <- function(
        transcript_path = tx_path,
        cell_labels_dir = poly_use_path,
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
        image_negative_y = NULL,
        load_expression = FALSE,
        load_cellmeta = TRUE,
        load_transcripts = TRUE,
        instructions = NULL,
        cores = determine_cores(),
        verbose = NULL,
        load_protein = !is.null(protein_expression_path),
        protein_expression_path = protein_expr_path,
        protein_metadata_path = protein_meta_path,
        protein_feat_type = c("protein", "protein_neg"),
        protein_split_keyword = NULL,
        protein_drop = "^Channel-",
        protein_meta_cols = c("nCount_RNA", "nFeature_RNA",
            "nCount_negprobes", "nFeature_negprobes", "Panel", "assay_type"),
        keep_cells = "intersect") {
            load_expression <- as.logical(load_expression)
            load_cellmeta <- as.logical(load_cellmeta)
            load_transcripts <- as.logical(load_transcripts)
            load_protein <- isTRUE(as.logical(
                load_protein %null% !is.null(protein_expression_path)
            ))
            if (load_protein) {
                keep_cells <- match.arg(keep_cells, choices = "intersect")
                checkmate::assert_character(protein_feat_type,
                    min.len = 1L, any.missing = FALSE)
            }
            if (!load_expression && !load_transcripts && !load_protein) {
                stop(wrap_txt("[CosmxReader] nothing to build the object
                    from. Set `load_expression = TRUE`, `load_transcripts =
                    TRUE` or load protein."), call. = FALSE)
            }

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
            if (isTRUE(load_transcripts)) {
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
            }

            # matrices before polygons: only cells with expression get polygons
            exlist <- plist <- keep <- cx <- NULL
            if (load_protein || load_expression) {
                rna_ids <- NULL
                if (load_expression) {
                    exlist <- funs$load_expression(
                        path = expression_path,
                        feat_type = feat_type,
                        split_keyword = split_keyword
                    )
                    rna_ids <- colnames(exlist[[1L]][])
                } else if (load_transcripts && length(metadata_path)) {
                    cx <- funs$load_cellmeta(
                        path = metadata_path,
                        cores = cores,
                        verbose = verbose
                    )
                    rna_ids <- cx[]$cell_ID
                }
                if (load_protein) {
                    plist <- funs$load_protein(
                        path = protein_expression_path,
                        feat_type = protein_feat_type,
                        split_keyword = protein_split_keyword,
                        drop_cols = protein_drop,
                        cores = cores,
                        verbose = verbose
                    )
                }
                keep <- .cosmx_match_cells(list(
                    rna = rna_ids,
                    protein = if (length(plist)) colnames(plist[[1L]][])
                ), keep_cells = keep_cells)
            }

            # polys
            poly_args <- list(
                path = cell_labels_dir,
                verbose = FALSE
            )
            if (!is.null(image_negative_y)) {
                # negative_y override
                if (isTRUE(image_negative_y)) {
                    poly_args$shift_vertical_step <- FALSE
                } else {
                    poly_args$shift_vertical_step <- 1
                }
            }
            if (!is.null(keep)) poly_args$keep_ids <- keep
            polys <- do.call(funs$load_polys, poly_args)
            if (!is.null(keep)) {
                # matrix order; the polygon csv is sorted by cell, then fov
                pid <- spatIDs(polys)
                .cosmx_report_cells(keep, pid)
                polys <- polys[match(keep[keep %in% pid], pid)]
            }
            g <- setGiotto(g, polys, verbose = verbose)

            # images
            if (!is.null(load_images)) {
                # replace convenient shortnames
                load_images[load_images == "composite"] <- composite_img_dir
                load_images[load_images == "overlay"] <- overlay_img_dir

                imglist <- list()
                dirnames <- names(load_images)
                for (imdir_i in seq_along(load_images)) {
                    img_args <- list(
                        path = load_images[[imdir_i]],
                        img_type = dirnames[[imdir_i]],
                        verbose = verbose
                    )
                    if (!is.null(image_negative_y)) {
                        # negative_y override
                        img_args$negative_y <- image_negative_y
                    }
                    dir_imgs <- do.call(funs$load_images, img_args)
                    imglist <- c(imglist, dir_imgs)
                }
                g <- addGiottoLargeImage(
                    g, 
                    largeImages = imglist, 
                    verbose = FALSE)
            }

            # expression & meta
            # Need to check that names agree for poly/expr/meta
            allowed_ids <- spatIDs(polys)

            # GiottoClass does not check that feature types share cell order
            for (ex in c(exlist, plist)) {
                idx <- match(allowed_ids, colnames(ex[]))
                if (anyNA(idx)) {
                    stop(wrap_txt("[CosmxReader] polygon cells are missing",
                        "from the", featType(ex), "expression."),
                        call. = FALSE)
                }
                ex[] <- ex[][, idx, drop = FALSE]
                g <- setGiotto(g, ex, verbose = verbose)
            }

            if (load_cellmeta) {
                if (is.null(cx)) {
                    cx <- funs$load_cellmeta(
                        path = metadata_path,
                        cores = cores,
                        verbose = verbose
                    )
                }

                if (!is.null(keep)) {
                    cx[] <- cx[][match(allowed_ids, cell_ID, nomatch = 0L), ]
                } else {
                    cx[] <- cx[][cell_ID %in% allowed_ids, ]
                }
                g <- setGiotto(g, cx, verbose = verbose)

                if (load_protein && !is.null(protein_metadata_path)) {
                    ren <- .cosmx_protein_meta_names(protein_meta_cols,
                        protein_feat_type)
                    px <- funs$load_cellmeta(
                        path = protein_metadata_path,
                        feat_type = protein_feat_type[[1L]],
                        select = c("Area.um2", names(ren)),
                        rename = ren,
                        cores = cores,
                        verbose = verbose
                    )
                    px[] <- px[][match(allowed_ids, cell_ID, nomatch = 0L), ]
                    g <- setGiotto(g, px, verbose = verbose)
                }
            }
            
            # create spatlocs
            g <- addSpatialCentroidLocations(g, verbose = FALSE)
            
            # add fovs
            g$fov <- as.integer(gsub("^c_\\d+_(\\d+)_\\d+$", "\\1",
                pDataDT(g)$cell_ID))
            if (load_protein) {
                g <- .cosmx_add_fov(g, feat_type = c(feat_type[[1L]],
                    protein_feat_type[[1L]]))
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
    basic_info <- c("cosmx_dir", "protein_dir", "version", "slide", "fovs",
                    "micron", "px2um", "poly_pref", "offsets")
    if (name %in% basic_info) {
        return(methods::slot(x, name))
    }

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("CosmxReader"), function(x, name, value) {
    basic_info <- c(
        "cosmx_dir", "protein_dir", "version", "slide", "fovs", "micron",
        "px2um", "poly_pref"
    )
    if (name %in% basic_info) {
        if (name == "protein_dir") {
            .cosmx_check_protein_slot(x)
            if (is.null(value)) value <- character(0L)
        }
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
    dn <- c("cosmx_dir", "protein_dir", "version", "slide", "fovs",
            "micron", "px2um", "poly_pref", "offsets")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}





# MODULAR ####

## paths ####

# Detect CosMx output file paths within a directory. Returns a named
# list of resolved paths. Callers (the reader's initialize and disk-
# subclass init frames) typically `list2env(paths, envir = environment())`
# into the frame whose closures reference these names via default-arg
# expressions.
.cosmx_detect_paths <- function(cosmx_dir) {
    if (length(cosmx_dir) == 0L) return(list())
    .cosmx_detect <- function(pattern) {
        .detect_in_dir(pattern = pattern, path = cosmx_dir,
                       platform = "CosMx")
    }
    list(
        shifts_path       = .cosmx_detect("fov_positions_file"),
        meta_path         = .cosmx_detect("metadata_file"),
        tx_path           = .cosmx_detect("tx_file"),
        mask_dir          = .cosmx_detect("CellLabels"),
        poly_path         = .cosmx_detect("polygons"),
        expr_path         = .cosmx_detect("exprMat_file"),
        composite_img_dir = .cosmx_detect("CellComposite"),
        overlay_img_dir   = .cosmx_detect("CellOverlay"),
        compart_img_dir   = .cosmx_detect("CompartmentLabels")
    )
}

# a subclass installed against an older CosmxReader has no protein_dir slot
.cosmx_protein_dir <- function(x) {
    if (methods::.hasSlot(x, "protein_dir")) x@protein_dir else character()
}

.cosmx_check_protein_slot <- function(x) {
    if (methods::.hasSlot(x, "protein_dir")) return(invisible(TRUE))
    stop(wrap_txt(sprintf(
        "[importCosMx] the %s class in use has no `protein_dir` slot: a
        package that extends CosmxReader (e.g. GiottoDisk) was installed
        against an older Giotto. Reinstall it to read CosMx protein data.",
        class(x)[[1L]]
    )), call. = FALSE)
}

# protein files share names with the RNA export, so they get their own folder
.cosmx_protein_paths <- function(protein_dir) {
    out <- list(
        protein_expr_path = NULL,
        protein_meta_path = NULL,
        protein_plex_path = NULL
    )
    if (length(protein_dir) == 0L) return(out)
    checkmate::assert_directory_exists(protein_dir)
    .detect <- function(pattern, warn = TRUE) {
        .detect_in_dir(pattern = pattern, path = protein_dir,
                       platform = "CosMx protein", warn = warn)
    }
    out["protein_expr_path"] <- list(.detect("exprMat_file"))
    out["protein_meta_path"] <- list(.detect("metadata_file"))
    out["protein_plex_path"] <- list(.detect("^plex-.*[.]txt$", warn = FALSE))
    out
}

# isotype names contain spaces and regex characters, so they are escaped
.cosmx_protein_negatives <- function(plex_path = NULL, expr_path = NULL) {
    neg <- character()
    if (length(plex_path) && file.exists(plex_path)) {
        plex <- data.table::fread(plex_path)
        if (all(c("DisplayName", "CodeClass") %in% colnames(plex))) {
            neg <- plex$DisplayName[plex$CodeClass == "Negative"]
        }
    }
    if (length(neg) && length(expr_path) && file.exists(expr_path)) {
        neg <- intersect(neg, .cosmx_csv_header(expr_path))
    }
    if (!length(neg)) return(list("IgG"))
    esc <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", neg)
    list(paste0("^(", paste(esc, collapse = "|"), ")$"))
}

.cosmx_check_protein_header <- function(header, split_keyword, file) {
    if (any(grepl("^(Negative|SystemControl)|NegPrb", header))) {
        warning(wrap_txt("[CosmxReader]", file, "has RNA negative-probe or
            system-control columns: is `protein_dir` the RNA export?"),
            call. = FALSE)
    }
    if (length(split_keyword) && !any(grepl(split_keyword[[1L]], header))) {
        warning(wrap_txt("[CosmxReader] no isotype control in", file,
            "matches", split_keyword[[1L]], "; every feature loads as a
            target"), call. = FALSE)
    }
    invisible(NULL)
}

# the vendor's protein metadata reuses RNA column names for protein totals
.cosmx_protein_meta_names <- function(cols,
    feat_type = c("protein", "protein_neg")) {
    if (!length(cols)) return(character())
    neg <- if (length(feat_type) > 1L) feat_type[[2L]] else
        paste0(feat_type[[1L]], "_neg")
    new <- sub("negprobes", neg, cols, fixed = TRUE)
    new <- sub("RNA", feat_type[[1L]], new, fixed = TRUE)
    same <- new == cols
    new[same] <- paste0(cols[same], "_", feat_type[[1L]])
    stats::setNames(new, cols)
}

.cosmx_match_cells <- function(ids, keep_cells = "intersect") {
    keep_cells <- match.arg(keep_cells, choices = "intersect")
    ids <- ids[lengths(ids) > 0L]
    if (!length(ids)) return(NULL)
    keep <- Reduce(function(a, b) a[a %in% b], ids)
    if (!length(keep)) {
        stop(wrap_txt(
            "[CosmxReader] the RNA and protein matrices share no cell; both
            should hold 'c_<slide>_<fov>_<cell>' ids of the same slide."
        ), call. = FALSE)
    }
    attr(keep, "dropped") <- lapply(ids, function(x) x[!x %in% keep])
    keep
}

.cosmx_report_cells <- function(keep, poly_ids) {
    d <- attr(keep, "dropped")
    dropped <- list(
        "with protein but no RNA expression" = d$protein,
        "with RNA but no protein expression" = d$rna,
        "with expression but no polygon" = keep[!keep %in% poly_ids]
    )
    dropped <- dropped[lengths(dropped) > 0L]
    if (!length(dropped)) return(invisible(NULL))
    n <- sum(lengths(dropped))
    why <- vapply(names(dropped), function(w) {
        sprintf("%d %s (e.g. %s)", length(dropped[[w]]), w,
            paste(utils::head(sort(dropped[[w]]), 2L), collapse = ", "))
    }, character(1L))
    warning(wrap_txt(sprintf(
        "[CosmxReader] %d %s not loaded: %s.",
        n, if (n == 1L) "cell is" else "cells are",
        paste(why, collapse = "; ")
    )), call. = FALSE)
    invisible(NULL)
}

.cosmx_add_fov <- function(g, feat_type) {
    have <- list_cell_metadata(g, spat_unit = "cell")
    for (ft in intersect(feat_type, have$feat_type)) {
        cx <- getCellMetadata(g, spat_unit = "cell", feat_type = ft,
            output = "data.table", copy_obj = FALSE)
        if ("fov" %in% colnames(cx)) next
        cx[, fov := as.integer(gsub("^c_\\d+_(\\d+)_\\d+$", "\\1", cell_ID))]
    }
    g
}

# first line only: fread() would decompress the whole file
.cosmx_csv_header <- function(path) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con))
    colnames(data.table::fread(text = readLines(con, n = 1L), sep = ",",
        header = TRUE))
}

# some exports write coordinates as `1e+05`
.cosmx_poly_col_classes <- function(path) {
    hdr <- .cosmx_csv_header(path)
    cls <- rep("character", length(hdr))
    cls[grepl("_px$", hdr)] <- "double"
    cls[hdr %in% c("fov", "cellID", "cell_ID")] <- "integer"
    cls
}


.cosmx_fov_shift <- function(path) {
    if (is.null(path)) return(NULL) # return early (empty case)
    fov_shifts <- data.table::fread(path)
    fs_colnames <- colnames(fov_shifts)
    # WTX datasets only have mm values for shifts.
    # mm->px conversion fragile if scalefactors not certain
    if ("X_mm" %in% fs_colnames) return(NULL) # return early (if WTX dataset)
    
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
        col_classes = NULL,
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
    
    if (is.null(fovs)) {
        tx <- data.table::fread(input = path, nThread = cores, drop = dropcols)
    } else {
        tx <- .cosmx_read_fovs(path, fovs = fovs, drop = dropcols,
            col_classes = col_classes, cores = cores)
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
    gpoints_params$verbose <- FALSE

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
    fov <- V1 <- NULL # NSE vars
    if (!missing(tx_dt) && !missing(meta_dt)) {
        stop("[.cosmx_infer_fov_shifts] Only one of tx_dt or meta_dt should be supplied\n",
             call. = FALSE)
    }
    
    if (!missing(tx_dt)) {
        tx_head <- tx_dt[, head(.SD, navg), by = fov]
        x <- tx_head[, mean(x_global_px - x_local_px), by = fov]
        if (flip_loc_y) {
            # test if flip is needed
            # Usual yshift variance / fov expected when correct is 0 to 1e-22
            # if var is too high for any fov, swap `flip_loc_y` value
            y <- tx_head[, var(y_global_px + y_local_px), by = fov]
            y <- subset(y, !is.na(V1)) # If FOV has only 1 entry, var returns NA and next step breaks
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
            y <- subset(y, !is.na(V1)) # If FOV has only 1 entry, var returns NA and next step breaks
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

# reads the FOV num of a cosmx image of style F(0)n where the number of prefixed
# 0s does not matter.
.cosmx_imgname_fovparser <- function(path) {
    if (length(path) == 1) im_names <- list.files(path)
    else im_names <- path
    fovs <- as.numeric(sub(".*F(\\d+)\\..*", "\\1", im_names))
    if (any(is.na(fovs))) {
        warning(wrap_txt(
            "Images to load should be sets of images/fov in subdirectories.
            No other files should be present."
        ), call. = FALSE)
        fovs <- fovs[!is.na(fovs)]
    }
    return(fovs)
}

# FOV width in image file names, read from the files: "F%03d" in older
# exports, "F%05d" in AtoMx ones; "F%03d" when mixed or unreadable.
.cosmx_file_fov_fmt <- function(path) {
    if (length(path) == 1) im_names <- list.files(path)
    else im_names <- basename(path)
    digits <- sub(".*F(\\d+)\\..*", "\\1", im_names)
    widths <- unique(nchar(digits[grepl("^\\d+$", digits)]))
    if (length(widths) == 1L) sprintf("F%%0%dd", widths) else "F%03d"
}

.cosmx_poly_maskimage <- function(path,
    slide = 1,
    fovs = NULL,
    file_fov_fmt = NULL,
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
    verbose = NULL,
    ...) {
    # NSE params
    f <- x <- y <- NULL
    
    path <- normalizePath(path)

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
    file_fov_fmt <- file_fov_fmt %null% .cosmx_file_fov_fmt(path)
    checkmate::assert_character(file_fov_fmt)
    with_pbar({
        p <- pbar(along = fovs)

        gpolys <- lapply(fovs, function(f) {
            if (length(path) == 1L) {
                search_fmt <- sprintf("%%s/*%s*", file_fov_fmt)
                segfile <- Sys.glob(paths = sprintf(search_fmt, path, f))
            } else {
                segfile <- path[grepl(pattern = sprintf(file_fov_fmt, f), path)]
            }
            
            if (length(segfile) == 0L) {
                fmt <- sprintf(
                    "No label images found matching fov [%s]", file_fov_fmt
                )
                msg <- sprintf(fmt, f)
                msg <- paste(msg, "\nCheck if `file_fov_fmt` =", file_fov_fmt,
                             "is appropriate and that image for fov", f,
                             "exists")
                stop(call. = FALSE, msg)
            }
            if (length(segfile) > 1L) {
                warning(sprintf(
                    "Multiple label images found for fov [%d]\n Using first one.", f
                ))
                segfile <- segfile[1L]
            }
            
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
            p(message = sprintf(file_fov_fmt, f))
            return(gpoly)
        })
    })

    if (length(gpolys) > 1L) {
        gpolys <- do.call(rbind, args = gpolys)
    } else {
        gpolys <- gpolys[[1]] # unlist
    }
    
    return(gpolys)
}

.cosmx_poly_csv <- function(
        path,
        slide = 1,
        fovs = NULL,
        name = "cell",
        micron = FALSE,
        px2um = 0.12028,
        make_valid = TRUE,
        # cellID is an integer identifier per FOV
        # x_local_px and y_local_px are local xy coords per FOV
        dropcols = c("x_local_px", "y_local_px"),
        keep_ids = NULL,
        verbose = NULL,
        ...) {
    
    GiottoUtils::vmsg(.v = verbose, "loading segmentation polygons...")
    vmsg(.v = verbose, .is_debug = TRUE, path)
    
    checkmate::assert_file_exists(path)
    
    if (is.null(fovs)) {
        dt <- data.table::fread(path, drop = dropcols)
    } else {
        dt <- .cosmx_read_fovs(path, fovs = as.integer(fovs), drop = dropcols,
            col_classes = .cosmx_poly_col_classes(path))
    }

    dt[, cell := sprintf("c_%d_%d_%d", as.integer(slide), fov, cellID)]
    dt[, cellID := NULL]
    if (!is.null(keep_ids)) dt <- dt[cell %in% keep_ids]

    data.table::setcolorder(dt, c("x_global_px", "y_global_px", "cell", "fov"))
    createGiottoPolygon(dt, make_valid = make_valid, verbose = FALSE)
}

.cosmx_poly <- function(
        path,
        slide = 1,
        fovs = NULL,
        file_fov_fmt = NULL,
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
        keep_ids = NULL,
        verbose = NULL,
        ...) {
    
    a <- GiottoUtils::get_args_list(...)

    if (missing(path)) {
        c("No path to mask image subdirectory or polygons csv", 
        "provided or auto-detected") %>%
            wrap_txt() %>%
            stop(call. = FALSE)
    }
    
    if (length(path) > 1L || dir.exists(path)) {
        a$keep_ids <- NULL
        gpolys <- do.call(.cosmx_poly_maskimage, args = a)
        if (!is.null(keep_ids)) {
            gpolys <- gpolys[gpolys$poly_ID %in% keep_ids]
        }
    } else if ("csv" %in% GiottoUtils::file_extension(path)) {
        a$file_fov_fmt <- NULL # not accepted by this function
        gpolys <- do.call(.cosmx_poly_csv, args = a)
    } else {
        "importCosMx - load_polys(): unrecognized path input"
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
        feat_type = "rna",
        select = NULL,
        rename = NULL,
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

    if (is.null(select)) {
        meta_dt <- data.table::fread(input = path, nThread = cores)
    } else {
        hdr <- .cosmx_csv_header(path)
        meta_dt <- data.table::fread(input = path, nThread = cores,
            select = intersect(hdr, c("fov", "cell_ID", "cell", select)))
    }

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

    if (length(rename)) {
        old <- intersect(names(rename), colnames(meta_dt))
        data.table::setnames(meta_dt, old = old, new = unname(rename[old]))
    }

    cx <- createCellMetaObj(
        metadata = meta_dt,
        spat_unit = "cell",
        feat_type = feat_type,
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
        drop_cols = NULL,
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to exprMat file provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, "loading expression matrix...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    hdr <- .cosmx_csv_header(path)
    drop <- NULL
    if (length(drop_cols)) {
        drop <- grep(paste(drop_cols, collapse = "|"),
            setdiff(hdr, c("fov", "cell_ID")), value = TRUE)
        if (!length(drop)) drop <- NULL
    }

    if (identical(hdr[[1L]], "fov")) {
        if (!is.null(fovs)) fovs <- as.integer(fovs)
        expr_mat <- .cosmx_read_exprmat(path, header = hdr, fovs = fovs,
            slide = slide, drop = drop, cores = cores)
    } else {
        if (length(drop)) {
            expr_dt <- data.table::fread(input = path, nThread = cores,
                drop = drop)
        } else {
            expr_dt <- data.table::fread(input = path, nThread = cores)
        }

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
        na <- !stats::complete.cases(expr_dt)
        .cosmx_warn_na_cells(expr_dt$cell_ID[na], path)
        expr_dt <- expr_dt[!na]

        # convert to Matrix
        expr_mat <- dt_to_matrix(expr_dt)
        expr_mat <- t_flex(expr_mat)
    }

    # split expression for rna / negprb if any split keywords provided.
    # Output of this chunk should always be a named list of 1 or more matrices
    if (length(split_keyword) > 0) {
        expr_list <- vector(mode = "list", length = length(feat_type))
        names(expr_list) <- feat_type
        # iterate through other expr types
        for (key_i in seq_along(split_keyword)) {
            feat_ids <- rownames(expr_mat)
            bool <- grepl(pattern = split_keyword[[key_i]], x = feat_ids)
            if (!any(bool)) next
            # subset and store split matrix
            expr_list[[key_i + 1L]] <- expr_mat[bool, , drop = FALSE]
            # remaining matrix
            expr_mat <- expr_mat[!bool, , drop = FALSE]
        }
        # assign the main expr
        expr_list[[1L]] <- expr_mat
        empty <- vapply(expr_list, is.null, logical(1L))
        if (any(empty)) {
            vmsg(.v = verbose, "no features for feat_type",
                paste(names(expr_list)[empty], collapse = ", "))
        }
        expr_list <- expr_list[!empty]
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

# read_colmatch() would parse the whole file
.cosmx_read_fovs <- function(path, fovs, drop = NULL, col_classes = NULL,
    cores = determine_cores(), batch_lines = 2e6) {
    j <- match("fov", .cosmx_csv_header(path))
    if (is.na(j)) {
        stop(wrap_txt("[CosmxReader] no `fov` column in", basename(path)),
            call. = FALSE)
    }
    con <- gzfile(path, open = "rt")
    on.exit(close(con))
    header_line <- readLines(con, n = 1L)
    kept <- list()
    repeat {
        lines <- readLines(con, n = batch_lines)
        if (!length(lines)) break
        lines <- lines[nzchar(lines)]
        if (!length(lines)) next
        f <- data.table::fread(text = lines, header = FALSE, select = j,
            sep = ",", nThread = cores)[[1L]]
        kept[[length(kept) + 1L]] <- lines[f %in% fovs]
    }
    data.table::fread(text = c(header_line, unlist(kept, use.names = FALSE)),
        sep = ",", header = TRUE, drop = drop, colClasses = col_classes,
        nThread = cores)
}

# batches of lines, one sparse matrix each: memory follows the selected cells
.cosmx_read_exprmat <- function(path, header, fovs = NULL, slide = 1,
    drop = NULL, cores = determine_cores(),
    batch_values = 2e7) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con))
    header_line <- readLines(con, n = 1L)
    batch_rows <- max(1L, as.integer(batch_values / length(header)))

    # one join at the end: cbind() of many parts peaks far above the result
    ids <- idx <- vals <- n_col <- list()
    feats <- NULL
    na_ids <- character(0L)
    repeat {
        lines <- readLines(con, n = batch_rows)
        if (!length(lines)) break
        if (!is.null(fovs)) {
            f1 <- substr(lines, 1L, regexpr(",", lines, fixed = TRUE) - 1L)
            f1 <- suppressWarnings(as.integer(gsub("\"", "", f1,
                fixed = TRUE)))
            lines <- lines[f1 %in% fovs]
        }
        lines <- lines[nzchar(lines)]
        if (!length(lines)) next
        dt <- data.table::fread(text = c(header_line, lines), sep = ",",
            header = TRUE, drop = drop, nThread = cores)
        dt <- dt[cell_ID != 0L, ]
        if (!nrow(dt)) next
        b <- length(ids) + 1L
        ids[[b]] <- sprintf("c_%d_%d_%d", slide, dt$fov, dt$cell_ID)
        dt[, c("fov", "cell_ID") := NULL]
        bad <- names(dt)[!vapply(dt, function(x) {
            is.numeric(x) || all(is.na(x))
        }, logical(1L))]
        if (length(bad)) {
            stop(wrap_txt("[CosmxReader]", basename(path),
                "has non-numeric values in",
                ngettext(length(bad), "column", "columns"),
                toString(utils::head(bad, 5L))), call. = FALSE)
        }
        na <- !stats::complete.cases(dt)
        if (any(na)) {
            na_ids <- c(na_ids, ids[[b]][na])
            ids[[b]] <- ids[[b]][!na]
            dt <- dt[!na]
            if (!nrow(dt)) {
                ids[[b]] <- NULL
                next
            }
        }
        m <- methods::as(t(as.matrix(dt)), "CsparseMatrix")
        feats <- rownames(m)
        idx[[b]] <- m@i
        vals[[b]] <- m@x
        n_col[[b]] <- diff(m@p)
    }
    .cosmx_warn_na_cells(na_ids, path)
    if (!length(ids)) {
        stop(wrap_txt(
            "[CosmxReader] no cells found in", basename(path),
            if (!is.null(fovs)) paste("for fovs", paste(fovs, collapse = ", "))
        ), call. = FALSE)
    }
    ids <- unlist(ids, use.names = FALSE)
    mat <- methods::new("dgCMatrix",
        i = unlist(idx, use.names = FALSE),
        p = c(0L, cumsum(unlist(n_col, use.names = FALSE))),
        x = unlist(vals, use.names = FALSE),
        Dim = c(length(feats), length(ids)),
        Dimnames = list(feats, ids)
    )
    if (length(mat) <= 2 * Matrix::nnzero(mat, na.counted = TRUE)) {
        mat <- Matrix::Matrix(as.matrix(mat))
    }
    mat
}

.cosmx_warn_na_cells <- function(ids, path) {
    if (!length(ids)) return(invisible(NULL))
    n <- length(ids)
    warning(wrap_txt(sprintf(
        "[CosmxReader] %d %s in %s %s NA values and %s left out, as cells
        missing from it (e.g. %s)", n, if (n == 1L) "cell" else "cells",
        basename(path),
        if (n == 1L) "has" else "have", if (n == 1L) "is" else "are",
        paste(utils::head(ids, 2L), collapse = ", ")
    )), call. = FALSE)
    invisible(NULL)
}

.cosmx_image <- function(path,
    fovs = NULL,
    img_type = "composite",
    img_name_fmt = paste(img_type, "_fov%03d"),
    file_fov_fmt = NULL,
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
    path <- normalizePath(path)

    GiottoUtils::vmsg(.v = verbose, sprintf("loading %s images...", img_type))
    vmsg(.v = verbose, .is_debug = TRUE, path)

    fovs <- fovs %null% .cosmx_imgname_fovparser(path) # ALL if NULL
    file_fov_fmt <- file_fov_fmt %null% .cosmx_file_fov_fmt(path)
    checkmate::assert_character(file_fov_fmt)
    verbose <- verbose %null% TRUE

    with_pbar({
        p <- pbar(along = fovs)

        gimg_list <- lapply(fovs, function(f) {
            if (length(path) == 1L) {
                search_fmt <- sprintf("%%s/*%s*", file_fov_fmt)
                imgfile <- Sys.glob(paths = sprintf(search_fmt, path, f))
            } else {
                imgfile <- path[grepl(pattern = sprintf(file_fov_fmt, f), path)]
            }
            
            if (length(imgfile) == 0L) {
                fmt <- sprintf(
                    "No images found matching fov [%s]", file_fov_fmt
                )
                msg <- sprintf(fmt, f)
                msg <- paste(msg, "\nCheck if `file_fov_fmt` =", file_fov_fmt,
                             "is appropriate and that image for fov", f,
                             "exists")
                stop(call. = FALSE, msg)
            }
            if (length(imgfile) > 1L) {
                warning(sprintf(
                    "Multiple images found for fov [%d]\n Using first one.", f
                ))
                imgfile <- imgfile[1L]
            }

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
            p(message = sprintf(file_fov_fmt, f))
            return(gimg)
        })
    })


    return(gimg_list)
}









# wrapper ####


#' @title Create Nanostring CosMx Giotto Object
#' @name createGiottoCosMxObject
#' @description Given the path to a CosMx experiment directory, creates a Giotto
#' object. For lower level control over loading, please see [importCosMx()]
#' @param cosmx_dir full path to the exported cosmx directory
#' @param version character. Version of CosMx output. Current selections are
#' either "default", "v6", and "legacy" (for the NSCLC dataset).
#' @param FOVs field of views to load.
#' @param slide numeric. Slide number. Defaults to 1. This must be correct so
#' that cell_IDs will match across polygons and (if loaded) expression matrix.
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
#' @param load_transcripts logical. (Default = TRUE) whether to load provided
#' transcript detections
#' @param poly_pref character. "mask" or "csv". Determines whether to load in
#' the polygons from the mask images (default) or the csv polygons file.
#' @param poly_read_fun Optional function applied to the polygons file as it is
#' read, for repairing a vendor export in flight. Requires a `backend`: the
#' in-memory reader has no polygon read hook.
#' @param meta_cols Optional. With a `backend`, the cell metadata columns to
#' load (default of the disk reader: `fov`, `Area.um2`, `nCount_RNA`,
#' `nFeature_RNA`). The in-memory reader loads every column, so it is ignored
#' there with a warning.
#' @param backend (optional) a `gsource`-inheriting project backend (typically
#' produced by `GiottoDisk::sourceCreate()`). When provided, the `giotto`
#' object is created as a managed on-disk project; see [importCosMx()].
#' Polygons are then read from the `polygons.csv` file (`poly_pref = "csv"`);
#' mask polygons (`CellLabels`) are not supported on this path.
#' @param image_negative_y Optional logical. Whether images are assumed to map
#' to positive or negative y values before fov shifts are applied. Affects
#' images (and polygons generated from masks). This overrides any settings from
#' selecting `version`.
#' @param protein_dir Optional. Directory of the protein flat files of a CosMx
#' RNA + protein run (`exprMat_file`, `metadata_file`, `plex-*.txt`). See
#' details and [importCosMx()].
#' @param load_protein logical. Whether to load the protein expression.
#' Default (`NULL`): `TRUE` when a protein `exprMat_file` is found or given.
#' @param protein_feat_type character. Feature types of the protein targets
#' and of the isotype controls. Default `c("protein", "protein_neg")`; one
#' name keeps the isotypes with the targets.
#' @param keep_cells character. Which cells to keep when protein is loaded.
#' `"intersect"` (only option): cells in the RNA export, protein expression
#' and a polygon. The RNA cells are those of the RNA matrix, or, when RNA is
#' loaded only as transcripts, those of the RNA cell metadata; with neither,
#' every protein cell with a polygon is kept.
#' @param fov_shifts_path Optional. Filepath to fov_positions_file
#' @param transcript_path Optional. Filepath to desired transcripts file to
#' load.
#' @param cell_labels_dir Optional. Path to directory containing CellLabels
#' images to load as polygons.
#' @param expression_path Optional. Filepath to cell feature matrix to load.
#' @param metadata_path Optional. Filepath to metadata file to load.
#' @param protein_expression_path,protein_metadata_path Optional. Filepaths
#' to the protein expression matrix and metadata, when they are not in
#' `protein_dir`.
#' @param cores nthreads to use when loading in cell metadata and transcripts
#' @param data_to_use deprecated. Not used
#' @param remove_background_polygon deprecated. Now always done
#' @param background_algo deprecated. Not used
#' @param remove_unvalid_polygons deprecated. Now always done
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @inheritDotParams importCosMx -cosmx_dir -fovs -slide
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
#'   \item{experimentname-\strong{polygons}.csv (file)}
#' }
#' 
#' The subdirectories should contain only the image files in order to be loaded
#' correctly.
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
#' 
#' \[**Polygons**\] Some outputs provide both the mask images and the 
#' `polygons.csv`. Giotto uses the mask images by default to convert to polygons
#' info. However, if only the `polygons.csv` file is present or 
#' `poly_pref = csv` is set, then the csv file will be used instead. One thing
#' to watch out for when loading from the `csv` is that in some datasets, the 
#' polygons preferably should not overlap each other. Giotto is not fully 
#' compatible with overlapping annotations. It also opens up the possibility 
#' of double counting transcripts. Polygon overlaps do not appear to be an
#' issue in the most recent outputs.
#' 
#' \[**Protein (RNA + protein runs)**\] With `protein_dir`, the protein
#' exprMat is loaded as two feature types: the targets (`"protein"`) and the
#' isotype controls (`"protein_neg"`, `CodeClass` "Negative" in the plex file).
#' Values are the vendor's per-cell mean intensities, stored unchanged as
#' `"raw"`; the `Channel-*` morphology columns are not loaded (they repeat the
#' metadata `Mean.*` columns). Cells are matched across modalities
#' (`keep_cells = "intersect"`): only cells in the RNA export (see
#' `keep_cells`), protein expression and a polygon are kept, every feature
#' type in the same cell order, and a warning gives the number dropped per
#' reason. Protein
#' metadata goes to feature type `"protein"`: the vendor protein file reuses
#' RNA column names for protein totals, so they are renamed
#' (`nCount_RNA` -> `nCount_protein`, which is the sum of the targets and the
#' `Channel-*` columns; `nFeature_RNA` -> `nFeature_protein`;
#' `nCount_negprobes` -> `nCount_protein_neg`, the isotype sum;
#' `nFeature_negprobes` -> `nFeature_protein_neg`; `Panel` ->
#' `Panel_protein`; `assay_type` -> `assay_type_protein`), plus `Area.um2`
#' and `fov`. Protein feature names can equal RNA gene names (`CD68`): pass
#' `feat_type = "protein"` to downstream functions.
#'
#' \[**Spatial Alignment Issues**\] Different versions of the CosMx output have
#' changed how images (and polygons generated from masks) should be aligned 
#' relative to the vector information (transcript detections). The `version`
#' param affects how the images are mapped by default. `image_negative_y` is
#' a toggle that overrides the image mapping style when an appropriate `version`
#' is difficult to determine.
#' \code{\link{showGiottoImageNames}} can be used to see the available images.
#' @md
#' @examples
#' \dontrun{
#' f <- "file path to cosmx flatfile output directory"
#' createGiottoCosMxObject(f)
#' 
#' # load older CosMx format
#' createGiottoCosMxObject(f, version = "legacy")
#' 
#' # force images and mask image polygons to shift up one image height
#' createGiottoCosMxObject(f, image_negative_y = FALSE)
#' 
#' # load only aggregated data
#' createGiottoCosMxObject(f,
#'     load_cellmeta = TRUE,
#'     load_expression = TRUE,
#'     load_transcripts = FALSE
#'     # data filepaths not needed unless they are not in expected locations
#' )
#' }
#' @export
createGiottoCosMxObject <- function(
        cosmx_dir,
        version = "default",
        FOVs = NULL,
        slide = 1,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        load_images = list(composite = "composite", overlay = "overlay"),
        load_expression = FALSE,
        load_cellmeta = TRUE,
        load_transcripts = TRUE,
        poly_pref = "mask",
        poly_read_fun = NULL,
        meta_cols = NULL,
        image_negative_y = NULL,
        protein_dir = NULL,
        load_protein = NULL,
        protein_feat_type = c("protein", "protein_neg"),
        keep_cells = "intersect",
        
        # optional filepaths
        fov_shifts_path = NULL,
        transcript_path = NULL,
        cell_labels_dir = NULL,
        expression_path = NULL,
        metadata_path = NULL,
        protein_expression_path = NULL,
        protein_metadata_path = NULL,
        cores = determine_cores(),
        verbose = NULL,
        instructions = NULL,
        backend = NULL,
        
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
    x <- importCosMx(cosmx_dir,
        fovs = FOVs,
        slide = slide,
        poly_pref = poly_pref,
        version = version,
        backend = backend,
        protein_dir = protein_dir,
        ...)
    if (is.null(backend)) x$version <- version
    if (!is.null(fov_shifts_path)) {
        checkmate::assert_file_exists(fov_shifts_path)
        x$offsets <- data.table::fread(fov_shifts_path)
    }

    if (!is.null(poly_read_fun) && is.null(backend)) {
        stop("[createGiottoCosMxObject] `poly_read_fun` needs a `backend`: ",
             "the in-memory reader has no polygon read hook.", call. = FALSE)
    }
    if (!is.null(meta_cols) && is.null(backend)) {
        warning("[createGiottoCosMxObject] `meta_cols` is ignored without a ",
                "`backend`: the in-memory reader loads every metadata column.",
                call. = FALSE)
    }

    # the disk reader's create_gobject takes the in-memory reader's arguments
    load_args <- list(
        feat_type = feat_type,
        split_keyword = split_keyword,
        load_images = load_images,
        image_negative_y = image_negative_y,
        load_expression = load_expression,
        load_cellmeta = load_cellmeta,
        load_transcripts = load_transcripts,
        instructions = instructions,
        cores = cores,
        verbose = verbose
    )

    if (!is.null(transcript_path)) load_args$transcript_path <- transcript_path
    if (!is.null(cell_labels_dir)) load_args$cell_labels_dir <- cell_labels_dir
    if (!is.null(expression_path)) load_args$expression_path <- expression_path
    if (!is.null(metadata_path)) load_args$metadata_path <- metadata_path
    if (!is.null(poly_read_fun)) load_args$poly_read_fun <- poly_read_fun
    if (!is.null(meta_cols) && !is.null(backend)) {
        load_args$meta_cols <- meta_cols
    }

    # RNA-only calls must still work with a GiottoDisk without protein support
    if (!is.null(protein_dir) || !is.null(protein_expression_path) ||
        isTRUE(load_protein)) {
        load_args$protein_feat_type <- protein_feat_type
        load_args$keep_cells <- keep_cells
        if (!is.null(load_protein)) load_args$load_protein <- load_protein
        if (!is.null(protein_expression_path)) {
            load_args$protein_expression_path <- protein_expression_path
        }
        if (!is.null(protein_metadata_path)) {
            load_args$protein_metadata_path <- protein_metadata_path
        }
    }

    do.call(x$create_gobject, load_args)
}




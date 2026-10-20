# CLASS ####



setClass(
    "MerscopeReader",
    slots = list(
        merscope_dir = "character",
        score_threshold = "ANY",
        polygon_format = "character",
        calls = "list",
        paths = "list"
    ),
    prototype = list(
        merscope_dir = NA_character_,
        score_threshold = NULL,
        polygon_format = "parquet",
        calls = list(),
        paths = list()
    )
)

# * show ####
setMethod("show", signature("MerscopeReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", class(object)[1L]))
    print_slots <- c("dir", "polygon_format", "score_cutoff", "funs")
    pre <- sprintf("%s :", format(print_slots))
    names(pre) <- print_slots

    d <- object@merscope_dir
    if (length(d) > 0L && !is.na(d)) {
        cat(pre["dir"], GiottoUtils::str_abbreviate(d), "\n")
    } else {
        cat(pre["dir"], "\n")
    }
    cat(pre["polygon_format"], object@polygon_format, "\n")
    sc <- object@score_threshold
    cat(pre["score_cutoff"], if (is.null(sc)) "none" else sc, "\n")
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("MerscopeReader"), function(x, ...) show(x))

# * init ####
setMethod(
    "initialize", signature("MerscopeReader"),
    function(.Object, ..., merscope_dir, score_threshold, polygon_format) {
        obj <- callNextMethod(.Object, ...)

        if (!missing(merscope_dir)) {
            checkmate::assert_directory_exists(merscope_dir)
            obj@merscope_dir <- normalizePath(path.expand(merscope_dir))
        }
        if (!missing(score_threshold)) {
            obj@score_threshold <- score_threshold
        }
        if (!missing(polygon_format)) {
            obj@polygon_format <- match.arg(
                polygon_format, c("parquet", "hdf5")
            )
        }

        obj@paths <- .merscope_detect_paths(obj@merscope_dir)
        list2env(obj@paths, envir = environment())
        score_th <- obj@score_threshold
        poly_fmt <- obj@polygon_format
        bound_path <- if (poly_fmt == "hdf5") cell_bound_dir else cell_bound_path

        # transcripts
        tx_fun <- function(
        path = tx_path,
        feat_type = c("rna", "neg_probe"),
        split_keyword = list("Blank"),
        FOVs = NULL,
        dropcols = c("x", "y", "barcode_id", "transcript_id"),
        score_threshold = score_th,
        flip_vertical = TRUE,
        output = c("giottoPoints", "data.table"),
        cores = determine_cores(),
        verbose = NULL) {
            .merscope_transcript(
                path = path,
                feat_type = feat_type,
                split_keyword = split_keyword,
                FOVs = FOVs,
                dropcols = dropcols,
                score_threshold = score_threshold,
                flip_vertical = flip_vertical,
                output = output,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_transcripts <- tx_fun

        # polygons
        poly_fun <- function(
        path = bound_path,
        name = "cell",
        poly_z_indices = NULL,
        polygon_format = poly_fmt,
        FOVs = NULL,
        metadata_path = cell_meta_path,
        flip_vertical = TRUE,
        calc_centroids = TRUE,
        verbose = NULL) {
            .merscope_poly(
                path = path,
                name = name,
                poly_z_indices = poly_z_indices,
                polygon_format = polygon_format,
                FOVs = FOVs,
                metadata_path = metadata_path,
                flip_vertical = flip_vertical,
                calc_centroids = calc_centroids,
                verbose = verbose
            )
        }
        obj@calls$load_polys <- poly_fun

        # expression
        expr_fun <- function(
        path = expr_path,
        spat_unit = "cell",
        feat_type = c("rna", "neg_probe"),
        split_keyword = "Blank",
        split_by_type = TRUE,
        remove_zero_rows = TRUE,
        cores = determine_cores(),
        verbose = NULL) {
            .merscope_expression(
                path = path,
                spat_unit = spat_unit,
                feat_type = feat_type,
                split_keyword = split_keyword,
                split_by_type = split_by_type,
                remove_zero_rows = remove_zero_rows,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_expression <- expr_fun

        # cell metadata
        cmeta_fun <- function(
        path = cell_meta_path,
        spat_unit = "cell",
        FOVs = NULL,
        cores = determine_cores(),
        verbose = NULL) {
            .merscope_cellmeta(
                path = path,
                spat_unit = spat_unit,
                FOVs = FOVs,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_cellmeta <- cmeta_fun

        # spatial locations (cell centers from cell metadata)
        sloc_fun <- function(
        path = cell_meta_path,
        spat_unit = "cell",
        FOVs = NULL,
        flip_vertical = TRUE,
        cores = determine_cores(),
        verbose = NULL) {
            .merscope_spatlocs(
                path = path,
                spat_unit = spat_unit,
                FOVs = FOVs,
                flip_vertical = flip_vertical,
                cores = cores,
                verbose = verbose
            )
        }
        obj@calls$load_spatlocs <- sloc_fun

        # images
        img_fun <- function(
        path,
        name = "image",
        transform_path = transform_path,
        negative_y = TRUE,
        verbose = NULL) {
            .merscope_image(
                path = path,
                name = name,
                transform_path = transform_path,
                negative_y = negative_y,
                verbose = verbose
            )
        }
        obj@calls$load_image <- img_fun

        # create giotto object
        gobject_fun <- function(
        transcript_path = tx_path,
        boundary_path = bound_path,
        expression_path = expr_path,
        metadata_path = cell_meta_path,
        data_to_use = c("subcellular", "aggregate"),
        feat_type = c("rna", "neg_probe"),
        split_keyword = list("Blank"),
        poly_z_indices = NULL,
        FOVs = NULL,
        dropcols = c("x", "y", "barcode_id", "transcript_id"),
        spat_unit = "cell",
        load_transcripts = TRUE,
        load_expression = TRUE,
        load_cellmeta = FALSE,
        load_images = NULL,
        image_transform_path = transform_path,
        instructions = NULL,
        verbose = NULL) {
            data_to_use <- match.arg(data_to_use, c("subcellular", "aggregate"))
            load_transcripts <- as.logical(load_transcripts)
            load_expression <- as.logical(load_expression)
            load_cellmeta <- as.logical(load_cellmeta)

            if (data_to_use == "aggregate") {
                load_transcripts <- FALSE
                load_expression <- TRUE
            }
            if (!load_transcripts && !load_expression) {
                warning(wrap_txt(
                    "One of either transcripts or expression info should be
                    loaded for a fully functioning object"
                ))
            }

            funs <- obj@calls

            g <- giotto(instructions = instructions)

            # transcripts
            if (load_transcripts) {
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    FOVs = FOVs,
                    dropcols = dropcols,
                    verbose = verbose
                )
                g <- setGiotto(g, tx_list, verbose = FALSE)
            }

            # polygons
            if (data_to_use == "subcellular" && is.na(boundary_path)) {
                warning(wrap_txt(sprintf(
                    "No cell boundaries found for polygon_format = '%s'.
                    The object is created without polygons.",
                    poly_fmt
                )), call. = FALSE)
            }
            if (data_to_use == "subcellular" && !is.na(boundary_path)) {
                b <- funs$load_polys(
                    path = boundary_path,
                    name = spat_unit,
                    poly_z_indices = poly_z_indices,
                    FOVs = FOVs,
                    metadata_path = metadata_path,
                    verbose = verbose
                )
                if (!is.null(b)) g <- setGiotto(g, b, verbose = FALSE)
            }

            # expression
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    split_keyword = unlist(split_keyword),
                    split_by_type = TRUE,
                    remove_zero_rows = TRUE,
                    verbose = verbose
                )
                g <- setGiotto(g, ex, verbose = FALSE)
            }

            # cell metadata
            if (load_cellmeta) {
                cx <- funs$load_cellmeta(
                    path = metadata_path,
                    spat_unit = spat_unit,
                    FOVs = FOVs,
                    verbose = verbose
                )
                g <- setGiotto(g, cx, verbose = FALSE)
            }

            # images
            if (!is.null(load_images)) {
                img_paths <- .merscope_image_paths(load_images, image_dir)
                imglist <- list()
                for (img_i in seq_along(img_paths)) {
                    im <- funs$load_image(
                        path = img_paths[[img_i]],
                        name = names(img_paths)[[img_i]],
                        transform_path = image_transform_path,
                        verbose = verbose
                    )
                    imglist <- c(imglist, im)
                }
                g <- setGiotto(g, imglist, verbose = FALSE)
            }

            # centroids from polygons, else cell centers from cell metadata
            spat_units_to_calc <- list_spatial_info_names(g)
            if (length(spat_units_to_calc) > 0L) {
                vmsg(.v = verbose, "calculating centroids")
                g <- addSpatialCentroidLocations(g,
                    poly_info = spat_units_to_calc,
                    provenance = as.list(spat_units_to_calc),
                    verbose = FALSE
                )
            } else if (!is.na(metadata_path)) {
                sl <- funs$load_spatlocs(
                    path = metadata_path,
                    spat_unit = spat_unit,
                    FOVs = FOVs,
                    verbose = verbose
                )
                g <- setGiotto(g, sl, verbose = FALSE)
            }

            vmsg(.v = verbose, "done")
            return(g)
        }
        obj@calls$create_gobject <- gobject_fun

        return(obj)
    }
)



# access ####

#' @export
setMethod("$", signature("MerscopeReader"), function(x, name) {
    basic_info <- c("merscope_dir", "polygon_format", "score_threshold")
    if (name %in% basic_info) {
        return(methods::slot(x, name))
    }
    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("MerscopeReader"), function(x, name, value) {
    basic_info <- c("merscope_dir", "polygon_format", "score_threshold")
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
`.DollarNames.MerscopeReader` <- function(x, pattern) {
    dn <- c("merscope_dir", "polygon_format", "score_threshold")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}



# CREATE READER ####

#' @title Import a Vizgen MERSCOPE Assay
#' @name importMerscope
#' @description
#' Giotto import functionalities for MERSCOPE datasets. This function creates
#' a `MerscopeReader` instance that has convenient reader functions for
#' converting individual pieces of MERSCOPE data into Giotto-compatible
#' representations.
#'
#' The reader loads the vendor outputs: transcripts
#' (`detected_transcripts.parquet` or `.csv`), cell boundaries
#' (`cell_boundaries.parquet`, or the per-FOV HDF5 files in `cell_boundaries/`),
#' the cell by gene matrix (`cell_by_gene.csv`) and cell metadata
#' (`cell_metadata.csv`). No overlap between transcripts and polygons is
#' calculated. Use [calculateOverlap()] and [overlapToMatrix()] on the created
#' object to aggregate transcripts by polygon.
#'
#' These functions should have all param values provided as defaults, but
#' can be flexibly modified to do things such as look in alternative
#' directories or paths.
#' @param merscope_dir MERSCOPE region directory (the one containing
#' `detected_transcripts.parquet`)
#' @param score_threshold numeric or `NULL` (default). Minimum
#' `transcript_score` for a transcript detection to be loaded. `NULL` loads
#' all transcripts.
#' @param polygon_format format of the cell boundaries, either `"parquet"`
#' (default, `cell_boundaries.parquet`) or `"hdf5"` (per-FOV files in
#' `cell_boundaries/`)
#' @param backend (optional) a `gsource`-inheriting project backend (typically
#' produced by `GiottoDisk::sourceCreate()`). When provided, creates the
#' `giotto` object as a managed on-disk project.
#' @returns `MerscopeReader` object, or `MerscopeDiskReader` when `backend` is
#' set
#' @examples
#' \dontrun{
#' x <- importMerscope("path/to/region_R1")
#' g <- x$create_gobject()
#' }
#' @export
importMerscope <- function(merscope_dir = NULL,
    score_threshold = NULL,
    polygon_format = c("parquet", "hdf5"),
    backend = NULL) {
    polygon_format <- match.arg(polygon_format, c("parquet", "hdf5"))
    if (!is.null(backend)) {
        package_check(
            "GiottoDisk",
            repository = "github:giotto-suite/GiottoDisk"
        )
        return(GiottoDisk::importMerscopeDisk(
            merscope_dir = merscope_dir,
            backend = backend,
            score_threshold = score_threshold,
            polygon_format = polygon_format
        ))
    }
    a <- list(
        Class = "MerscopeReader",
        score_threshold = score_threshold,
        polygon_format = polygon_format
    )
    if (!is.null(merscope_dir)) {
        a$merscope_dir <- merscope_dir
    }
    do.call(new, args = a)
}



# MODULAR ####

## paths ####

.merscope_detect_paths <- function(merscope_dir) {
    path_names <- c(
        "tx_path", "cell_bound_path", "cell_bound_dir", "expr_path",
        "cell_meta_path", "image_dir", "transform_path"
    )
    if (length(merscope_dir) == 0L || is.na(merscope_dir)) {
        out <- as.list(rep(NA_character_, length(path_names)))
        names(out) <- path_names
        return(out)
    }

    f <- function(...) {
        p <- file.path(merscope_dir, ...)
        if (file.exists(p)) normalizePath(p) else NULL
    }
    first <- function(...) {
        for (p in list(...)) if (!is.null(p)) return(p)
        NULL
    }

    paths <- list(
        tx_path = first(
            f("detected_transcripts.parquet"), f("detected_transcripts.csv")
        ),
        cell_bound_path = f("cell_boundaries.parquet"),
        cell_bound_dir = if (dir.exists(file.path(merscope_dir, "cell_boundaries"))) {
            f("cell_boundaries")
        },
        expr_path = f("cell_by_gene.csv"),
        cell_meta_path = f("cell_metadata.csv"),
        image_dir = f("images"),
        transform_path = f("images", "micron_to_mosaic_pixel_transform.csv")
    )

    if (is.null(paths$tx_path) && is.null(paths$expr_path)) {
        stop(wrap_txt(
            "[merscope] no `detected_transcripts.*` or `cell_by_gene.csv`
            found in:", merscope_dir
        ), call. = FALSE)
    }

    lapply(paths, function(p) if (is.null(p)) NA_character_ else p)
}



## transcripts ####

.merscope_transcript <- function(
        path,
        feat_type = c("rna", "neg_probe"),
        split_keyword = list("Blank"),
        FOVs = NULL,
        dropcols = c("x", "y", "barcode_id", "transcript_id"),
        score_threshold = NULL,
        flip_vertical = TRUE,
        output = c("giottoPoints", "data.table"),
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path) || is.na(path)) {
        stop(wrap_txt(
            "No path to transcripts file provided or auto-detected"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    output <- match.arg(output, choices = c("giottoPoints", "data.table"))

    vmsg("Loading transcript level info...", .v = verbose)
    tx <- .merscope_read_tx(
        path = path,
        FOVs = FOVs,
        dropcols = dropcols,
        cores = cores
    )

    transcript_score <- global_y <- cell_id <- NULL # NSE vars
    if (!is.null(score_threshold) && "transcript_score" %in% names(tx)) {
        n_before <- nrow(tx)
        tx <- tx[transcript_score >= score_threshold]
        vmsg(.v = verbose, sprintf(
            "score_threshold %s: %s of %s transcripts kept",
            score_threshold,
            format(nrow(tx), big.mark = ","),
            format(n_before, big.mark = ",")
        ))
    }
    if ("cell_id" %in% names(tx)) {
        tx[, cell_id := as.character(cell_id)]
    }

    if (isTRUE(flip_vertical)) {
        tx[, global_y := -global_y]
    }

    if (output == "data.table") return(tx)

    gpointslist <- createGiottoPoints(
        x = tx,
        x_colname = "global_x",
        y_colname = "global_y",
        feat_ID_colname = "gene",
        feat_type = feat_type,
        split_keyword = split_keyword,
        verbose = FALSE
    )
    if (!inherits(gpointslist, "list")) {
        gpointslist <- list(gpointslist)
    }
    return(gpointslist)
}

.merscope_read_tx <- function(path,
    FOVs = NULL,
    dropcols = character(0L),
    select = NULL,
    cores = determine_cores()) {
    fov <- cell_id <- NULL # NSE vars
    e <- tolower(file_extension(path)[[1L]])

    if (e == "parquet") {
        package_check("arrow")
        package_check("dplyr")
        sc <- arrow::schema(arrow::open_dataset(path))
        nms <- names(sc)
        unnamed <- !nzchar(nms)
        if (any(unnamed)) {
            nms[unnamed] <- paste0(".merscope_unnamed_", seq_len(sum(unnamed)))
        }
        sch <- do.call(arrow::schema, lapply(seq_along(nms), function(i) {
            arrow::field(nms[[i]], sc[[i]]$type)
        }))
        a <- arrow::open_dataset(sources = path, schema = sch)
        drop <- unique(c(nms[unnamed], as.character(dropcols)))
        if (!is.null(select)) drop <- setdiff(nms, select)
        if (length(drop) > 0L) {
            a <- dplyr::select(a, -dplyr::any_of(drop))
        }
        if (!is.null(FOVs)) {
            fovs_ <- as.integer(FOVs)
            a <- dplyr::filter(a, fov %in% fovs_)
        }
        if ("cell_id" %in% names(a)) {
            a <- dplyr::mutate(a, cell_id = arrow::cast(
                cell_id, arrow::string()
            ))
        }
        tx <- data.table::as.data.table(dplyr::collect(a))
    } else {
        hdr <- names(data.table::fread(path, nrows = 0L))
        first_line <- readLines(path, n = 1L)
        if (startsWith(first_line, ",")) hdr[[1L]] <- ".merscope_unnamed_1"
        keep <- if (is.null(select)) {
            setdiff(hdr, c(".merscope_unnamed_1", as.character(dropcols)))
        } else {
            intersect(hdr, select)
        }
        cc <- NULL
        if ("cell_id" %in% keep) cc <- c(cell_id = "character")
        if (is.null(FOVs)) {
            tx <- data.table::fread(
                path, select = keep, colClasses = cc, nThread = cores
            )
        } else {
            tx <- GiottoUtils::read_colmatch(
                file = path,
                col = "fov",
                values_to_match = FOVs,
                nThread = cores
            )
            tx <- tx[, intersect(keep, names(tx)), with = FALSE]
        }
    }
    tx
}

## polygons ####

.merscope_poly <- function(
        path,
        name = "cell",
        poly_z_indices = NULL,
        polygon_format = c("parquet", "hdf5"),
        FOVs = NULL,
        metadata_path = NULL,
        flip_vertical = TRUE,
        calc_centroids = TRUE,
        verbose = NULL) {
    if (missing(path) || is.na(path)) {
        stop(wrap_txt(
            "No path to cell boundaries provided or auto-detected"
        ), call. = FALSE)
    }
    checkmate::assert_character(name, len = 1L)
    polygon_format <- match.arg(polygon_format, c("parquet", "hdf5"))

    vmsg(.v = verbose, sprintf("Loading polygon info as '%s'...", name))

    gpoly <- switch(polygon_format,
        "parquet" = .merscope_poly_parquet(
            path = path,
            name = name,
            poly_z_indices = poly_z_indices,
            FOVs = FOVs,
            metadata_path = metadata_path,
            calc_centroids = calc_centroids,
            verbose = verbose
        ),
        "hdf5" = .merscope_poly_hdf5(
            path = path,
            name = name,
            poly_z_indices = poly_z_indices,
            FOVs = FOVs,
            flip_vertical = flip_vertical,
            calc_centroids = calc_centroids,
            verbose = verbose
        )
    )
    if (is.null(gpoly)) return(NULL)

    if (isTRUE(flip_vertical) && polygon_format == "parquet") {
        gpoly <- flip(gpoly, direction = "vertical", y0 = 0)
    }
    gpoly
}

.merscope_poly_parquet <- function(path,
    name = "cell",
    poly_z_indices = NULL,
    FOVs = NULL,
    metadata_path = NULL,
    calc_centroids = TRUE,
    verbose = NULL) {
    package_check("arrow")
    package_check("dplyr")
    checkmate::assert_file_exists(path)
    ZIndex <- EntityID <- Geometry <- NULL # NSE vars

    arch <- .merscope_zplane_architecture(
        path, poly_z_indices = poly_z_indices, verbose = verbose
    )

    q <- arrow::open_dataset(path)
    if (!is.na(arch$use_z) && "ZIndex" %in% names(q)) {
        use_z <- arch$use_z
        q <- dplyr::filter(q, ZIndex == use_z)
    }

    if (!is.null(FOVs)) {
        keep_ids <- .merscope_fov_cell_ids(metadata_path, FOVs)
        q <- dplyr::filter(q, EntityID %in% keep_ids)
    }

    x <- q |>
        dplyr::select(EntityID, Geometry) |>
        dplyr::collect()
    if (nrow(x) == 0L) {
        vmsg(.v = verbose, "[merscope] no polygons found")
        return(NULL)
    }

    sv <- terra::vect(as.list(x$Geometry))
    terra::values(sv) <- data.frame(
        poly_ID = as.character(x$EntityID),
        stringsAsFactors = FALSE
    )
    createGiottoPolygon(
        x = sv,
        name = name,
        calc_centroids = calc_centroids,
        verbose = FALSE
    )
}

.merscope_poly_hdf5 <- function(path,
    name = "cell",
    poly_z_indices = NULL,
    FOVs = NULL,
    flip_vertical = TRUE,
    calc_centroids = TRUE,
    verbose = NULL) {
    checkmate::assert_directory_exists(path)
    z <- if (is.null(poly_z_indices)) 1L:7L else as.integer(poly_z_indices)

    polys <- readPolygonFilesVizgenHDF5(
        boundaries_path = path,
        fovs = FOVs,
        z_indices = z,
        flip_y_axis = flip_vertical,
        calc_centroids = calc_centroids,
        verbose = verbose
    )
    if (length(polys) == 0L) return(NULL)

    if (length(polys) > 1L) {
        ids_a <- polys[[1L]]@unique_ID_cache
        same <- vapply(polys[-1L], function(p) {
            isTRUE(setequal(ids_a, p@unique_ID_cache)) &&
                isTRUE(all.equal(
                    utils::head(terra::crds(polys[[1L]]@spatVector), 100L),
                    utils::head(terra::crds(p@spatVector), 100L)
                ))
        }, logical(1L))
        if (!all(same)) {
            .merscope_stop_3d(z, reason = "polygons differ between z-planes")
        }
        vmsg(.v = verbose, sprintf(
            "[merscope] identical polygons across %d z-planes; using z-plane %s",
            length(z), z[[1L]]
        ))
    }

    gpoly <- polys[[1L]]
    objName(gpoly) <- name
    gpoly
}

.merscope_zplane_architecture <- function(path,
    poly_z_indices = NULL,
    n_sample = 200L,
    verbose = NULL) {
    ZIndex <- EntityID <- Geometry <- NULL # NSE vars
    ds <- arrow::open_dataset(path)
    if (!"ZIndex" %in% names(ds)) {
        return(list(kind = "2d", z_indices = NA_integer_, use_z = NA_integer_))
    }

    z_all <- ds |>
        dplyr::distinct(ZIndex) |>
        dplyr::collect() |>
        dplyr::pull(ZIndex) |>
        sort()
    if (!is.null(poly_z_indices)) {
        keep <- intersect(as.integer(poly_z_indices), as.integer(z_all))
        if (length(keep) == 0L) {
            stop(wrap_txt(
                "[merscope] none of `poly_z_indices` (",
                paste(poly_z_indices, collapse = ", "),
                ") are present. Available:", paste(z_all, collapse = ", ")
            ), call. = FALSE)
        }
        z_all <- keep
    }

    if (length(z_all) == 1L) {
        return(list(kind = "2d", z_indices = z_all, use_z = z_all[[1L]]))
    }

    za <- z_all[[1L]]
    zb <- z_all[[2L]]
    ids_a <- ds |>
        dplyr::filter(ZIndex == za) |>
        dplyr::select(EntityID) |>
        dplyr::collect() |>
        dplyr::pull(EntityID)
    ids_b <- ds |>
        dplyr::filter(ZIndex == zb) |>
        dplyr::select(EntityID) |>
        dplyr::collect() |>
        dplyr::pull(EntityID)
    if (!setequal(ids_a, ids_b)) {
        .merscope_stop_3d(z_all, reason = "cell ID sets differ between z-planes")
    }

    set.seed(1L)
    samp <- if (length(ids_a) > n_sample) sample(ids_a, n_sample) else ids_a
    g <- ds |>
        dplyr::filter(EntityID %in% samp) |>
        dplyr::select(EntityID, ZIndex, Geometry) |>
        dplyr::collect()
    sp <- split(g, g$EntityID)
    same <- vapply(sp, function(gi) {
        gi <- gi[order(gi$ZIndex), ]
        ref <- as.raw(gi$Geometry[[1L]])
        all(vapply(seq_len(nrow(gi)), function(i) {
            identical(as.raw(gi$Geometry[[i]]), ref)
        }, logical(1L)))
    }, logical(1L))
    if (!all(same)) {
        .merscope_stop_3d(z_all, reason = sprintf(
            "geometry differs across z-planes for %d of %d sampled cells",
            sum(!same), length(same)
        ))
    }

    vmsg(.v = verbose, sprintf(
        "[merscope] identical polygons across %d z-planes; using z-plane %s",
        length(z_all), za
    ))
    list(kind = "replicated_2d", z_indices = z_all, use_z = za)
}

.merscope_stop_3d <- function(z_indices, reason) {
    stop(wrap_txt(sprintf(
        "[merscope] true 3D segmentation detected (%s).
        Only 2D and replicated-2D boundaries are supported.
        Available z-indices: %s
        To load a single plane, pass e.g. poly_z_indices = %s",
        reason, paste(z_indices, collapse = ", "), z_indices[[1L]]
    )), call. = FALSE)
}

.merscope_fov_cell_ids <- function(metadata_path, FOVs) {
    if (is.null(metadata_path) || is.na(metadata_path)) {
        stop(wrap_txt(
            "[merscope] `FOVs` requires cell_metadata.csv to map cells to
            FOVs, but no metadata path was found."
        ), call. = FALSE)
    }
    cm <- .merscope_read_cellmeta(metadata_path)
    keep_ids <- cm[["cell_ID"]][cm[["fov"]] %in% as.integer(FOVs)]
    if (length(keep_ids) == 0L) {
        stop("[merscope] no cells found in FOV(s): ",
            paste(FOVs, collapse = ", "), call. = FALSE)
    }
    keep_ids
}



## expression ####

.merscope_feat_classes <- function(feat_ids,
    split_keyword = "Blank",
    feat_type = c("rna", "neg_probe")) {
    kw <- as.character(unlist(split_keyword))
    if (length(feat_type) != length(kw) + 1L) {
        stop("[merscope] `feat_type` needs one name more than `split_keyword`",
            call. = FALSE)
    }
    out <- rep(feat_type[[1L]], length(feat_ids))
    for (i in seq_along(kw)) {
        hit <- grepl(kw[[i]], feat_ids, fixed = TRUE)
        out[hit] <- feat_type[[i + 1L]]
    }
    out
}

.merscope_expression <- function(
        path,
        spat_unit = "cell",
        feat_type = c("rna", "neg_probe"),
        split_keyword = "Blank",
        split_by_type = TRUE,
        remove_zero_rows = TRUE,
        cores = determine_cores(),
        verbose = NULL) {
    if (missing(path) || is.na(path)) {
        stop(wrap_txt(
            "No path to cell_by_gene.csv provided or auto-detected"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    vmsg("Loading cell by gene matrix...", .v = verbose)

    dt <- data.table::fread(
        path, colClasses = c(cell = "character"), nThread = cores
    )
    cell_ids <- dt[["cell"]]
    feat_ids <- setdiff(names(dt), "cell")

    i_list <- vector("list", length(feat_ids))
    j_list <- vector("list", length(feat_ids))
    x_list <- vector("list", length(feat_ids))
    for (k in seq_along(feat_ids)) {
        v <- dt[[feat_ids[[k]]]]
        nz <- which(v != 0)
        i_list[[k]] <- rep.int(k, length(nz))
        j_list[[k]] <- nz
        x_list[[k]] <- as.numeric(v[nz])
    }
    mat <- Matrix::sparseMatrix(
        i = unlist(i_list, use.names = FALSE),
        j = unlist(j_list, use.names = FALSE),
        x = unlist(x_list, use.names = FALSE),
        dims = c(length(feat_ids), length(cell_ids)),
        dimnames = list(feat_ids, cell_ids)
    )
    rm(dt, i_list, j_list, x_list)

    feat_classes <- .merscope_feat_classes(
        feat_ids, split_keyword = split_keyword, feat_type = feat_type
    )

    if (isTRUE(remove_zero_rows)) {
        keep <- Matrix::rowSums(mat != 0) > 0
        mat <- mat[keep, , drop = FALSE]
        feat_classes <- feat_classes[keep]
    }

    uniq_classes <- unique(feat_classes)
    if (length(uniq_classes) > 1L && isTRUE(split_by_type)) {
        mat_list <- lapply(uniq_classes, function(fc) {
            mat[feat_classes == fc, , drop = FALSE]
        })
        names(mat_list) <- uniq_classes
    } else {
        nm <- if (length(uniq_classes) == 1L) uniq_classes else feat_type[[1L]]
        mat_list <- stats::setNames(list(mat), nm)
    }

    lapply(names(mat_list), function(ft) {
        createExprObj(
            expression_data = mat_list[[ft]],
            name = "raw",
            spat_unit = spat_unit,
            feat_type = ft,
            provenance = spat_unit
        )
    })
}



## cell metadata ####

.merscope_read_cellmeta <- function(path, cores = determine_cores()) {
    cell_ID <- NULL # NSE var
    first_line <- readLines(path, n = 1L)
    hdr <- names(data.table::fread(path, nrows = 0L))
    id_col <- if ("EntityID" %in% hdr) "EntityID" else hdr[[1L]]
    cc <- stats::setNames("character", id_col)
    dt <- data.table::fread(path, colClasses = cc, nThread = cores)
    if (startsWith(first_line, ",") && !"EntityID" %in% hdr) {
        data.table::setnames(dt, 1L, "cell_ID")
    } else {
        data.table::setnames(dt, id_col, "cell_ID")
    }
    dt[, cell_ID := as.character(cell_ID)]
    dt
}

.merscope_cellmeta <- function(path,
    spat_unit = "cell",
    FOVs = NULL,
    cores = determine_cores(),
    verbose = NULL) {
    if (missing(path) || is.na(path)) {
        stop(wrap_txt(
            "No path to cell_metadata.csv provided or auto-detected"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    vmsg("Loading cell metadata...", .v = verbose)

    dt <- .merscope_read_cellmeta(path, cores = cores)
    if (!is.null(FOVs) && "fov" %in% names(dt)) {
        dt <- dt[dt[["fov"]] %in% as.integer(FOVs)]
    }
    createCellMetaObj(
        metadata = dt,
        spat_unit = spat_unit,
        feat_type = "rna",
        provenance = spat_unit
    )
}

.merscope_spatlocs <- function(path,
    spat_unit = "cell",
    FOVs = NULL,
    flip_vertical = TRUE,
    cores = determine_cores(),
    verbose = NULL) {
    if (missing(path) || is.na(path)) {
        stop(wrap_txt(
            "No path to cell_metadata.csv provided or auto-detected"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    vmsg("Loading cell centers from cell metadata...", .v = verbose)
    sdimy <- NULL # NSE var

    dt <- .merscope_read_cellmeta(path, cores = cores)
    if (!is.null(FOVs) && "fov" %in% names(dt)) {
        dt <- dt[dt[["fov"]] %in% as.integer(FOVs)]
    }
    sl <- data.table::data.table(
        sdimx = dt[["center_x"]],
        sdimy = dt[["center_y"]],
        cell_ID = dt[["cell_ID"]]
    )
    if (isTRUE(flip_vertical)) {
        sl[, sdimy := -sdimy]
    }
    createSpatLocsObj(
        coordinates = sl,
        name = "raw",
        spat_unit = spat_unit,
        provenance = spat_unit,
        verbose = FALSE
    )
}



## images ####

.merscope_image_paths <- function(load_images, image_dir) {
    load_images <- as.list(load_images)
    nms <- names(load_images)
    if (is.null(nms)) nms <- rep("", length(load_images))
    paths <- character(length(load_images))
    for (i in seq_along(load_images)) {
        v <- as.character(load_images[[i]])
        if (file.exists(v)) {
            paths[[i]] <- normalizePath(v)
        } else {
            if (is.na(image_dir)) {
                stop(wrap_txt(
                    "[merscope] no `images` directory found for image", v
                ), call. = FALSE)
            }
            paths[[i]] <- file.path(image_dir, paste0("mosaic_", v, ".tif"))
            if (!file.exists(paths[[i]])) {
                stop(wrap_txt(
                    "[merscope] image not found:", paths[[i]]
                ), call. = FALSE)
            }
        }
        if (!nzchar(nms[[i]])) {
            nms[[i]] <- sub("^mosaic_", "", sub("(?i)\\.tiff?$", "",
                basename(paths[[i]]), perl = TRUE))
        }
    }
    stats::setNames(as.list(paths), nms)
}

.merscope_image <- function(path,
    name = "image",
    transform_path = NULL,
    negative_y = TRUE,
    verbose = NULL) {
    if (missing(path)) {
        stop(wrap_txt(
            "No path to image file provided"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    if (is.null(transform_path) || is.na(transform_path)) {
        stop(wrap_txt(
            "[merscope] no micron_to_mosaic_pixel_transform.csv found for
            image placement"
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(transform_path)
    vmsg(.v = verbose, sprintf("loading image as '%s'", name))

    tm <- as.matrix(data.table::fread(transform_path, header = FALSE))
    sx <- tm[1L, 1L]
    sy <- tm[2L, 2L]
    tx <- tm[1L, 3L]
    ty <- tm[2L, 3L]

    img <- createGiottoLargeImage(
        path,
        name = name,
        negative_y = negative_y,
        verbose = FALSE
    )
    img <- rescale(img, fx = 1 / sx, fy = 1 / sy, x0 = 0, y0 = 0)
    dy <- if (isTRUE(negative_y)) ty / sy else -ty / sy
    img <- spatShift(img, dx = -tx / sx, dy = dy)
    img
}

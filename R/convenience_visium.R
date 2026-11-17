# Standard 10x Visium (spot-level). Spots carry no vendor polygons or
# transcript coordinates; circular spot polygons are generated from the
# scalefactors.

# CLASS ####

setClass(
    "VisiumReader",
    slots = list(
        visium_dir = "character",
        expression_source = "character",
        feature_id_type = "character",
        calls = "list",
        paths = "list"
    ),
    prototype = list(
        visium_dir = NA_character_,
        expression_source = "raw",
        feature_id_type = "symbol",
        calls = list(),
        paths = list()
    )
)

# * show ####
setMethod("show", signature("VisiumReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", class(object)[1L]))
    print_slots <- c("dir", "expression_source", "feature_id_type", "funs")
    pre <- sprintf("%s :", format(print_slots))
    names(pre) <- print_slots

    d <- object@visium_dir
    if (length(d) > 0L && !is.na(d)) {
        cat(pre["dir"], GiottoUtils::str_abbreviate(d), "\n")
    } else {
        cat(pre["dir"], "\n")
    }
    cat(pre["expression_source"], object@expression_source, "\n")
    cat(pre["feature_id_type"], object@feature_id_type, "\n")
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("VisiumReader"), function(x, ...) show(x))

# * init ####
setMethod(
    "initialize", signature("VisiumReader"),
    function(.Object, ..., visium_dir, expression_source, feature_id_type) {
        obj <- callNextMethod(.Object, ...)

        if (!missing(visium_dir)) {
            checkmate::assert_directory_exists(visium_dir)
            obj@visium_dir <- normalizePath(path.expand(visium_dir))
        }
        if (!missing(expression_source)) {
            obj@expression_source <- expression_source
        }
        if (!missing(feature_id_type)) {
            obj@feature_id_type <- feature_id_type
        }
        checkmate::assert_choice(obj@expression_source, c("raw", "filtered"),
            .var.name = "expression_source"
        )
        checkmate::assert_choice(obj@feature_id_type, c("symbol", "ensembl"),
            .var.name = "feature_id_type"
        )

        obj@paths <- .visium_detect_paths(
            obj@visium_dir, expression_source = obj@expression_source
        )
        list2env(obj@paths, envir = environment())
        feat_id <- obj@feature_id_type

        # expression
        expr_fun <- function(
        path = expr_path,
        feature_id_type = feat_id,
        verbose = NULL) {
            .visium_expression(
                path = path,
                feature_id_type = feature_id_type,
                verbose = verbose
            )
        }
        obj@calls$load_expression <- expr_fun

        # spot locations
        tp_fun <- function(
        path = positions_path,
        barcodes = NULL,
        output = c("spatLocsObj", "data.table"),
        verbose = NULL) {
            output <- match.arg(output, c("spatLocsObj", "data.table"))
            x <- .visium_tissue_positions(path = path, barcodes = barcodes)
            sl <- x[, c("cell_ID", "sdimx", "sdimy")]
            if (output == "data.table") return(sl)
            createSpatLocsObj(sl,
                name = "raw", spat_unit = "cell", provenance = "cell",
                verbose = verbose
            )
        }
        obj@calls$load_tissue_position <- tp_fun

        # spot metadata
        cmeta_fun <- function(
        path = positions_path,
        barcodes = NULL,
        feat_type = "rna",
        verbose = NULL) {
            x <- .visium_tissue_positions(path = path, barcodes = barcodes)
            x <- x[, c("cell_ID", "in_tissue", "array_row", "array_col")]
            lapply(feat_type, function(ft) {
                createCellMetaObj(x,
                    spat_unit = "cell", feat_type = ft, provenance = "cell",
                    verbose = FALSE
                )
            })
        }
        obj@calls$load_cellmeta <- cmeta_fun

        # scalefactors
        sf_fun <- function(path = json_path) {
            .visium_read_scalefactors(json_path = path)
        }
        obj@calls$load_scalefactors <- sf_fun

        # image
        img_fun <- function(
        path = img_path,
        scalefactors = .visium_read_scalefactors(json_path),
        verbose = NULL) {
            .visium_image(
                image_path = path,
                json_info = scalefactors,
                verbose = verbose
            )
        }
        obj@calls$load_image <- img_fun

        # spot polygons
        poly_fun <- function(
        spatlocs = tp_fun(),
        scalefactors = .visium_read_scalefactors(json_path)) {
            if (is.null(scalefactors)) {
                stop("[Visium] spot polygons need scalefactors_json.json",
                    call. = FALSE
                )
            }
            .visium_spot_poly(
                spatlocs = spatlocs,
                json_scalefactors = scalefactors
            )
        }
        obj@calls$load_polygon <- poly_fun

        # create giotto object
        gobject_fun <- function(
        load_expression = TRUE,
        load_spatlocs = TRUE,
        load_metadata = TRUE,
        load_image = TRUE,
        load_polygons = TRUE,
        expression_path = expr_path,
        tissue_positions_path = positions_path,
        scalefactors_path = json_path,
        image_path = img_path,
        feature_id_type = feat_id,
        instructions = NULL,
        verbose = NULL) {
            funs <- obj@calls
            g <- giotto(instructions = instructions)

            barcodes <- NULL
            feat_types <- "rna"
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    feature_id_type = feature_id_type,
                    verbose = verbose
                )
                g <- setGiotto(g, ex, verbose = FALSE)
                barcodes <- colnames(ex[[1L]][])
                feat_types <- vapply(ex, featType, character(1L))
            }

            sl <- funs$load_tissue_position(
                path = tissue_positions_path, barcodes = barcodes
            )
            if (load_spatlocs) g <- setGiotto(g, sl, verbose = FALSE)

            if (load_metadata) {
                cx <- funs$load_cellmeta(
                    path = tissue_positions_path,
                    barcodes = barcodes,
                    feat_type = feat_types
                )
                g <- setGiotto(g, cx, verbose = FALSE)
            }

            json <- funs$load_scalefactors(path = scalefactors_path)

            if (load_image && !is.null(image_path)) {
                img <- funs$load_image(
                    path = image_path, scalefactors = json, verbose = verbose
                )
                g <- setGiotto(g, img, verbose = FALSE)
            }

            if (!is.null(json)) {
                if (load_polygons) {
                    poly <- funs$load_polygon(spatlocs = sl, scalefactors = json)
                    g <- setPolygonInfo(
                        gobject = g,
                        x = poly,
                        centroids_to_spatlocs = FALSE,
                        verbose = FALSE,
                        initialize = TRUE
                    )
                }
                instructions(g, "micron_scale") <-
                    65 / json$spot_diameter_fullres
            }

            return(g)
        }
        obj@calls$create_gobject <- gobject_fun

        return(obj)
    }
)



# access ####

.visium_editable_slots <- c("visium_dir", "expression_source", "feature_id_type")

#' @export
setMethod("$", signature("VisiumReader"), function(x, name) {
    if (name %in% .visium_editable_slots) {
        return(methods::slot(x, name))
    }
    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("VisiumReader"), function(x, name, value) {
    if (!name %in% .visium_editable_slots) {
        stop(sprintf(
            "Only items in '%s' can be set",
            toString(.visium_editable_slots)
        ))
    }
    methods::slot(x, name) <- value
    return(initialize(x))
})

#' @export
`.DollarNames.VisiumReader` <- function(x, pattern) {
    dn <- .visium_editable_slots
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}



# CREATE READER ####

#' @title Import a Visium assay
#' @name importVisium
#' @description Giotto import functionalities for standard 10x Visium
#' (spot-level) datasets. This function creates a `VisiumReader` instance that
#' has convenient reader functions for converting individual pieces of a Space
#' Ranger `outs` directory into Giotto-compatible representations.
#'
#' The reader loads the feature-barcode matrix (`<source>_feature_bc_matrix`
#' directory, or its `.h5`), the spot positions (`spatial/tissue_positions.csv`,
#' or `tissue_positions_list.csv` from Space Ranger 1), the scalefactors
#' (`spatial/scalefactors_json.json`) and a tissue image (the first `.png` in
#' `spatial/`). Visium has no vendor polygons; circular spot polygons are
#' created from `spot_diameter_fullres`.
#'
#' These functions should have all param values provided as defaults, but
#' can be flexibly modified to do things such as look in alternative
#' directories or paths.
#' @param visium_dir Space Ranger `outs` directory for a standard Visium run
#' @param expression_source character. `"raw"` (default) or `"filtered"`
#' feature-barcode matrix.
#' @param feature_id_type character. One of `"symbol"` or `"ensembl"`.
#'   Determines which to use as the feature identifiers.
#' @param backend (optional) a `gsource`-inheriting project backend (typically
#' produced by `GiottoDisk::sourceCreate()`). When provided, creates the
#' `giotto` object as a managed on-disk project.
#' @returns `VisiumReader` object, or `VisiumDiskReader` when `backend` is set
#' @examples
#' \dontrun{
#' x <- importVisium("path/to/outs", expression_source = "filtered")
#' sl <- x$load_tissue_position()
#' g <- x$create_gobject()
#' }
#' @export
importVisium <- function(visium_dir = NULL,
    expression_source = c("raw", "filtered"),
    feature_id_type = c("symbol", "ensembl"),
    backend = NULL) {
    expression_source <- match.arg(expression_source, c("raw", "filtered"))
    feature_id_type <- match.arg(feature_id_type, c("symbol", "ensembl"))
    if (!is.null(backend)) {
        package_check(
            "GiottoDisk",
            repository = "github:giotto-suite/GiottoDisk"
        )
        return(GiottoDisk::importVisiumDisk(
            visium_dir = visium_dir,
            backend = backend,
            expression_source = expression_source,
            feature_id_type = feature_id_type
        ))
    }
    a <- list(
        Class = "VisiumReader",
        expression_source = expression_source,
        feature_id_type = feature_id_type
    )
    if (!is.null(visium_dir)) {
        a$visium_dir <- visium_dir
    }
    do.call(new, args = a)
}



# MODULAR ####

## paths ####

.visium_detect_paths <- function(visium_dir, expression_source = "raw") {
    path_names <- c(
        "expr_path", "spatial_dir", "positions_path", "json_path", "img_path"
    )
    out <- vector("list", length(path_names))
    names(out) <- path_names
    if (length(visium_dir) == 0L || is.na(visium_dir)) {
        return(out)
    }

    # matrix directory first: folder reading has always used it
    stem <- paste0(expression_source, "_feature_bc_matrix")
    mtx_dir <- file.path(visium_dir, stem)
    out$expr_path <- if (dir.exists(mtx_dir)) {
        mtx_dir
    } else {
        .detect_in_dir(visium_dir,
            pattern = paste0(stem, "\\.h5$"),
            platform = "Visium", warn = FALSE
        )
    }
    if (is.null(out$expr_path)) {
        stop(wrap_txt(
            "[Visium] no", stem, "directory or .h5 found in:", visium_dir
        ), call. = FALSE)
    }

    spatial_dir <- file.path(visium_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
        stop(wrap_txt("[Visium] no spatial/ directory in:", visium_dir),
            call. = FALSE
        )
    }
    out$spatial_dir <- spatial_dir

    sp <- function(f) {
        p <- file.path(spatial_dir, f)
        if (file.exists(p)) p else NULL
    }
    # Space Ranger >= 2 writes a header; Space Ranger 1 the headerless list
    out$positions_path <- sp("tissue_positions.csv")
    if (is.null(out$positions_path)) {
        out$positions_path <- sp("tissue_positions_list.csv")
    }
    if (is.null(out$positions_path)) {
        stop(wrap_txt(
            "[Visium] no tissue_positions.csv or tissue_positions_list.csv in:",
            spatial_dir
        ), call. = FALSE)
    }

    out$json_path <- sp("scalefactors_json.json")

    pngs <- list.files(spatial_dir, pattern = "\\.png$", full.names = TRUE)
    if (length(pngs) > 0L) out$img_path <- pngs[[1L]]

    out
}



## expression ####

# exprObj list: "Gene Expression" as rna, "Antibody Capture" as protein
.visium_expression <- function(path,
    feature_id_type = c("symbol", "ensembl"),
    verbose = NULL) {
    feature_id_type <- match.arg(feature_id_type, c("symbol", "ensembl"))
    checkmate::assert_true(file.exists(path))
    vmsg(.v = verbose, "loading expression from", basename(path))

    if (dir.exists(path)) {
        expr_results <- get10Xmatrix(
            path_to_data = path,
            gene_column_index = switch(feature_id_type,
                symbol = 2L,
                ensembl = 1L
            )
        )
    } else {
        expr_results <- get10Xmatrix_h5(
            path_to_data = path,
            gene_ids = switch(feature_id_type,
                symbol = "symbols",
                ensembl = "ensembl"
            )
        )
    }
    if (!is.list(expr_results)) {
        expr_results <- list("Gene Expression" = expr_results)
    }

    mats <- list(rna = expr_results[["Gene Expression"]])
    if ("Antibody Capture" %in% names(expr_results)) {
        mats$protein <- expr_results[["Antibody Capture"]]
    }

    lapply(names(mats), function(ft) {
        createExprObj(mats[[ft]],
            name = "raw", spat_unit = "cell", feat_type = ft,
            provenance = "cell"
        )
    })
}



## spot positions ####

# spot table with `cell_ID`, `in_tissue`, `array_row`, `array_col` and the
# fullres pixel locations as `sdimx` / `sdimy` (y flipped). Columns are taken
# by position, so both the Space Ranger 1 headerless list and the
# Space Ranger >= 2 csv read the same. `barcodes` orders the rows to match.
.visium_tissue_positions <- function(path, barcodes = NULL) {
    # NSE vars
    cell_ID <- pxl_row <- pxl_col <- NULL

    checkmate::assert_file_exists(path)
    x <- data.table::fread(path)
    colnames(x) <- c(
        "cell_ID", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"
    )
    if (!is.null(barcodes)) x <- x[match(barcodes, cell_ID)]
    x[, `:=`(sdimx = pxl_col, sdimy = -pxl_row)]
    x[]
}



# CONVENIENCE ####

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
#' @param expression_matrix_class deprecated. See
#'   [GiottoClass::createExprObj] for details
#' @param h5_file optional path to create an on-disk h5 file
#' @param verbose be verbose
#' @param backend (optional) a `gsource`-inheriting project backend (typically
#' produced by `GiottoDisk::sourceCreate()`). When provided, creates the
#' `giotto` object as a managed on-disk project.
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
#' @seealso [importVisium()]
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
    expression_matrix_class = deprecated(),
    h5_file = NULL,
    cores = NA,
    verbose = NULL,
    backend = NULL) {
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

    if (!missing(expression_matrix_class)) {
        deprecate_warn(
            "4.2.5",
            "createGiottoVisiumObject(expression_matrix_class)",
            details = "See `?GiottoClass::createExprObj` for details."
        )
    }

    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    a <- list()
    if (!is.null(h5_visium_path)) {
        vmsg(
            .v = verbose,
            "A path to an .h5 10X file was provided and will be used"
        )
        if (!file.exists(h5_visium_path)) {
            .gstop("The provided path ", h5_visium_path, " does not exist")
        }
        if (is.null(h5_tissue_positions_path)) {
            .gstop("A path to the tissue positions (.csv) needs to be provided
                to h5_tissue_positions_path")
        }
        if (!file.exists(h5_tissue_positions_path)) {
            .gstop(
                "The provided path ", h5_tissue_positions_path,
                " does not exist"
            )
        }
        if (!is.null(h5_image_png_path) && !file.exists(h5_image_png_path)) {
            .gstop(
                "The provided h5 image path ", h5_image_png_path,
                "does not exist.
                Set to NULL to exclude or provide the correct path."
            )
        }
        h5_gene_ids <- match.arg(h5_gene_ids, c("symbols", "ensembl"))
        x <- importVisium(
            feature_id_type = c(symbols = "symbol", ensembl = "ensembl")[[
                h5_gene_ids
            ]],
            backend = backend
        )
        a$expression_path <- h5_visium_path
        a["tissue_positions_path"] <- list(h5_tissue_positions_path)
        a["scalefactors_path"] <- list(h5_json_scalefactors_path)
        a["image_path"] <- list(h5_image_png_path)
    } else {
        vmsg(.v = verbose, "A structured visium directory will be used")
        if (is.null(visium_dir)) {
            .gstop("visium_dir needs to be a path to a visium directory")
        }
        visium_dir <- path.expand(visium_dir)
        if (!dir.exists(visium_dir)) .gstop(visium_dir, " does not exist!")
        expr_data <- match.arg(expr_data, choices = c("raw", "filter"))
        if (!gene_column_index %in% c(1, 2)) {
            .gstop("gene_column_index must be 1 (ensembl) or 2 (symbol)")
        }
        x <- importVisium(
            visium_dir = visium_dir,
            expression_source = c(raw = "raw", filter = "filtered")[[
                expr_data
            ]],
            feature_id_type = c("ensembl", "symbol")[[gene_column_index]],
            backend = backend
        )
        if (!is.null(png_name)) {
            png_path <- file.path(visium_dir, "spatial", png_name)
            if (!file.exists(png_path)) .gstop(png_path, " does not exist!")
            a$image_path <- png_path
        }
    }

    a$instructions <- instructions
    a$verbose <- verbose
    do.call(x$create_gobject, args = a)
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

## CLASS ####
# ------- ###

setClass(
    "StereoSeqReader",
    slots = list(
        stereoseq_dir = "character",
        type = "character",      # "bin" or "cell"
        bin_size = "character",  # e.g. "bin100" (only used when type = "bin")
        gene_column = "character",
        negative_y = "logical",
        calls = "list"
    ),
    prototype = list(
        type = "bin",
        bin_size = "bin100",
        gene_column = "geneName",
        negative_y = TRUE,
        calls = list()
    )
)


# * show ####
setMethod("show", signature("StereoSeqReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", class(object)))
    print_slots <- c("dir", "type", "bin_size", "gene_column", "negative_y", "funs")
    pre <- sprintf("%s :", format(print_slots))
    names(pre) <- print_slots

    # dir
    d <- object@stereoseq_dir
    if (length(d) > 0L) {
        d <- GiottoUtils::str_abbreviate(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }

    # type
    cat(pre["type"], object@type, "\n")

    # bin_size (only relevant for bin type)
    if (object@type == "bin") {
        cat(pre["bin_size"], object@bin_size, "\n")
    }

    # gene_column
    cat(pre["gene_column"], object@gene_column, "\n")

    # negative_y
    cat(pre["negative_y"], object@negative_y, "\n")

    # funs
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("StereoSeqReader"), function(x, ...) show(x))



#' @title Import a Stereo-seq assay
#' @name importStereoSeq
#' @description Giotto import functionalities for Stereo-seq datasets. This
#' function generates a `StereoSeqReader` instance with convenient reader
#' functions for loading individual pieces of Stereo-seq data into
#' Giotto-compatible representations when the param `stereoseq_dir` is provided.
#' A function that creates the full `giotto` object is also available. These
#' functions should have all param values provided as defaults, but can be
#' flexibly modified to look in alternative directories or paths.
#' @param stereoseq_dir filepath to the Stereo-seq `outs` directory.
#' @param type character. One of `"bin"` (square bin, default) or `"cell"`
#'   (cell bin / cellbin).
#' @param bin_size character. Bin size level to load when `type = "bin"`.
#'   One of `"bin1"`, `"bin5"`, `"bin10"`, `"bin20"`, `"bin50"`,
#'   `"bin100"`, `"bin150"`, `"bin200"`. Default is `"bin100"`. List
#'   available bins in your sample with
#'   `rhdf5::h5ls("path/to/*.tissue.gef")`.
#' @param gene_column character. Which column to use as feature identifiers.
#'   One of `"geneName"` (default) or `"geneID"`.
#' @param negative_y logical. Map data to negative y spatial values (default
#'   `TRUE`). Origin is placed at the upper-left instead of lower-left.
#' @details Loading functions are generated after `stereoseq_dir` is set.
#' @returns `StereoSeqReader` object
#' @examples
#' \dontrun{
#' # Create a StereoSeqReader
#' reader <- importStereoSeq(
#'     stereoseq_dir = "path/to/outs",
#'     type = "bin",
#'     bin_size = "bin100"
#' )
#'
#' # Load expression or spatial locations individually
#' expr <- reader$load_expression()
#' sl <- reader$load_spatlocs()
#' img <- reader$load_image()
#'
#' # Create a full giotto object
#' g <- reader$create_gobject()
#'
#' # Change bin size and recreate
#' reader$bin_size <- "bin50"
#' g <- reader$create_gobject()
#' }
#' @export
importStereoSeq <- function(
    stereoseq_dir = NULL,
    type = c("bin", "cell"),
    bin_size = "bin100",
    gene_column = c("geneName", "geneID"),
    negative_y = TRUE) {

    a <- list(Class = "StereoSeqReader")

    if (!is.null(stereoseq_dir)) {
        a$stereoseq_dir <- stereoseq_dir
    }

    a$type <- match.arg(type, c("bin", "cell"))
    a$bin_size <- bin_size
    a$gene_column <- match.arg(gene_column, c("geneName", "geneID"))
    a$negative_y <- as.logical(negative_y)

    do.call(new, args = a)
}

#' @rdname importStereoSeq
#' @param x `StereoSeqReader`
#' @param ... additional params to pass
#' @export
setMethod("plot", signature("StereoSeqReader", "missing"),
    function(x, ...) {
    sl <- x$load_spatlocs()
    plot(hull(sl), ...)
})


# * init ####
setMethod("initialize", signature("StereoSeqReader"), function(.Object,
    stereoseq_dir,
    type,
    bin_size,
    gene_column,
    negative_y) {

    # provided params (if any) -------------------------------------------- #
    if (!missing(stereoseq_dir)) {
        checkmate::assert_directory_exists(stereoseq_dir)
        .Object@stereoseq_dir <- stereoseq_dir
    }
    if (!missing(type)) {
        .Object@type <- type
    }
    if (!missing(bin_size)) {
        .Object@bin_size <- bin_size
    }
    if (!missing(gene_column)) {
        .Object@gene_column <- gene_column
    }
    if (!missing(negative_y)) {
        .Object@negative_y <- negative_y
    }

    # NULL case: return early if no path given
    if (length(.Object@stereoseq_dir) == 0) {
        return(.Object)
    }

    # validate settings --------------------------------------------------- #
    type <- .Object@type
    bin_size <- .Object@bin_size
    gene_column <- .Object@gene_column

    if (!type %in% c("bin", "cell")) {
        stop("[StereoSeq] `$type` must be \"bin\" or \"cell\"\n", call. = FALSE)
    }

    .bin_size_options <- c("bin1", "bin5", "bin10", "bin20",
                           "bin50", "bin100", "bin150", "bin200")
    if (type == "bin" && !bin_size %in% .bin_size_options) {
        stop(sprintf("[StereoSeq] `$bin_size` must be one of: %s\n",
            toString(.bin_size_options)), call. = FALSE)
    }

    if (!gene_column %in% c("geneName", "geneID")) {
        stop("[StereoSeq] `$gene_column` must be \"geneName\" or \"geneID\"\n",
            call. = FALSE)
    }

    # detect paths -------------------------------------------------------- #
    p <- .Object@stereoseq_dir
    image_dir <- file.path(p, "image")

    gef_path <- if (type == "bin") {
        .stereoseq_find_tissue_gef(p)
    } else {
        .stereoseq_find_cellbin_gef(p)
    }

    spat_unit <- if (type == "bin") bin_size else "cell"

    # setup closures ------------------------------------------------------ #

    ## expression load call
    expression_fun <- function(
        path = gef_path,
        type = .Object@type,
        bin_size = .Object@bin_size,
        gene_column = .Object@gene_column,
        spat_unit = if (.Object@type == "bin") .Object@bin_size else "cell",
        verbose = NULL) {
        .stereoseq_expression(
            path = path,
            type = type,
            bin_size = bin_size,
            gene_column = gene_column,
            spat_unit = spat_unit,
            verbose = verbose
        )
    }
    .Object@calls$load_expression <- expression_fun

    ## spatlocs load call
    spatlocs_fun <- function(
        path = gef_path,
        type = .Object@type,
        bin_size = .Object@bin_size,
        gene_column = .Object@gene_column,
        negative_y = .Object@negative_y,
        spat_unit = if (.Object@type == "bin") .Object@bin_size else "cell",
        verbose = NULL) {
        .stereoseq_spatlocs(
            path = path,
            type = type,
            bin_size = bin_size,
            gene_column = gene_column,
            negative_y = negative_y,
            spat_unit = spat_unit,
            verbose = verbose
        )
    }
    .Object@calls$load_spatlocs <- spatlocs_fun

    ## image load call
    load_image_fun <- function(
        path = image_dir,
        negative_y = .Object@negative_y,
        verbose = NULL) {
        .stereoseq_image(path = path, negative_y = negative_y, verbose = verbose)
    }
    .Object@calls$load_image <- load_image_fun

    ## mask load call
    .default_mask_path <- .stereoseq_find_mask(p)
    load_mask_fun <- function(
        path = .default_mask_path,
        negative_y = .Object@negative_y,
        verbose = NULL) {
        .stereoseq_mask(path = path, negative_y = negative_y, verbose = verbose)
    }
    .Object@calls$load_mask <- load_mask_fun

    ## binpoints load call
    binpoints_fun <- function(
        path = gef_path,
        type = .Object@type,
        bin_size = .Object@bin_size,
        gene_column = .Object@gene_column,
        negative_y = .Object@negative_y,
        spat_unit = if (.Object@type == "bin") .Object@bin_size else "cell",
        verbose = NULL) {
        .stereoseq_binpoints(
            path = path,
            type = type,
            bin_size = bin_size,
            gene_column = gene_column,
            negative_y = negative_y,
            spat_unit = spat_unit,
            verbose = verbose
        )
    }
    .Object@calls$load_binpoints <- binpoints_fun

    ## create gobject call
    # Capture paths under different names to avoid recursive default arg error
    # (parameter name == captured variable name causes circular reference in R)
    .default_gef_path <- gef_path
    .default_image_dir <- image_dir
    .default_mask_path2 <- .stereoseq_find_mask(p)

    gobject_fun <- function(
        load_expression = TRUE,
        load_spatlocs = TRUE,
        load_binpoints = FALSE,
        load_image = TRUE,
        load_mask = TRUE,
        type = .Object@type,
        bin_size = .Object@bin_size,
        gene_column = .Object@gene_column,
        negative_y = .Object@negative_y,
        gef_path = .default_gef_path,
        image_path = .default_image_dir,
        mask_path = .default_mask_path2,
        instructions = NULL,
        verbose = NULL) {

        spat_unit <- if (type == "bin") bin_size else "cell"

        # spatLocsObj requires expression (or polygons) to be present first.
        # When only using binpoints, spatial info is embedded in the GBP itself.
        if (load_spatlocs && !load_expression) {
            vmsg(.v = verbose,
                "[StereoSeq] load_spatlocs = TRUE ignored: requires",
                "load_expression = TRUE. Spatial coordinates are embedded",
                "inside the giottoBinPoints object.")
            load_spatlocs <- FALSE
        }

        expr_obj <- sl <- gbp <- NULL

        if (load_expression || load_spatlocs || load_binpoints) {
            # read gef once, then build all requested representations.
            # Only request the sub-datasets each output actually needs.
            what_needed <- character(0)
            if (load_expression || load_binpoints) what_needed <- c(
                what_needed, "expression")
            if (load_spatlocs || load_binpoints) what_needed <- c(
                what_needed, "spatlocs")

            gef_data <- .stereoseq_read_gef(
                path = gef_path,
                type = type,
                bin_size = bin_size,
                gene_column = gene_column,
                what = what_needed,
                verbose = verbose
            )

            # For bin type, pre-compute unique (x,y) bin positions once so that
            # build_expression and build_spatlocs can reuse them without each
            # running their own unique() pass.
            if (type == "bin" && !is.null(gef_data$exprDT)) {
                b <- unique(gef_data$exprDT[, c("x", "y")], by = c("x", "y"))
                b[, bin_ID := .I]
                gef_data$bins <- b
            }

            if (load_expression) {
                expr_obj <- .stereoseq_build_expression(
                    gef_data = gef_data,
                    type = type,
                    gene_column = gene_column,
                    spat_unit = spat_unit,
                    verbose = verbose
                )
            }

            if (load_spatlocs) {
                sl <- .stereoseq_build_spatlocs(
                    gef_data = gef_data,
                    type = type,
                    negative_y = negative_y,
                    spat_unit = spat_unit,
                    verbose = verbose
                )
            }

            if (load_binpoints) {
                if (!is.null(expr_obj) && !is.null(sl)) {
                    # Both already built above — avoid redundant 
                    # matrix + unique() work.
                    # Keep in sync with .stereoseq_build_binpoints() if that 
                    # function changes.
                    vmsg(.v = verbose, 
                        "Building giottoBinPoints from pre-built objects...")
                    gbp <- createGiottoBinPoints(
                        expr_values = expr_obj[[1L]],
                        spatial_locs = sl,
                        feat_type = "rna"
                    )
                    spatUnit(gbp) <- spat_unit
                    vmsg(.v = verbose, "Finished giottoBinPoints")
                } else {
                    gbp <- .stereoseq_build_binpoints(
                        gef_data = gef_data,
                        type = type,
                        gene_column = gene_column,
                        negative_y = negative_y,
                        spat_unit = spat_unit,
                        verbose = verbose
                    )
                }
            }
        }

        gimg <- NULL
        if (load_image) {
            gimg <- .stereoseq_image(
                path = image_path,
                negative_y = negative_y,
                verbose = verbose
            )
        }

        gpoly <- NULL
        if (load_mask) {
            if (is.null(mask_path)) {
                warning(
                    "[StereoSeq] No mask file found. Skipping mask polygons.",
                    call. = FALSE)
            } else {
                gpoly <- .stereoseq_mask(
                    path = mask_path,
                    negative_y = negative_y,
                    verbose = verbose
                )
            }
        }

        # assemble giotto
        vmsg(.v = verbose, "assembling object ...")
        g <- giotto(instructions = instructions)
        if (!is.null(expr_obj)) {
            g <- setGiotto(g, expr_obj, verbose = verbose)
        }
        if (!is.null(sl)) {
            g <- setGiotto(g, sl, verbose = verbose)
        }
        if (!is.null(gbp)) {
            g <- setGiotto(g, gbp, verbose = verbose)
        }
        if (!is.null(gimg)) {
            g <- setGiotto(g, gimg, verbose = verbose)
        }
        if (!is.null(gpoly)) {
            g <- setGiotto(g, gpoly, verbose = verbose)
        }
        gc(verbose = FALSE)
        g
    }
    .Object@calls$create_gobject <- gobject_fun

    return(.Object)
})


# * access ####

.stereoseq_editable_slots <- c(
    "stereoseq_dir", "type", "bin_size", "gene_column", "negative_y"
)

#' @export
setMethod("$", signature("StereoSeqReader"), function(x, name) {
    if (name %in% .stereoseq_editable_slots) {
        return(methods::slot(x, name))
    }
    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("StereoSeqReader"), function(x, name, value) {
    if (!name %in% .stereoseq_editable_slots) {
        stop(sprintf(
            "Only items in '%s' can be set",
            toString(.stereoseq_editable_slots)
        ))
    }
    methods::slot(x, name) <- value
    return(initialize(x))
})

#' @export
.DollarNames.StereoSeqReader <- function(x, pattern) {
    dn <- .stereoseq_editable_slots
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}



# MODULAR ####
# ---------- #

# Find *.tissue.gef under path/feature_expression/
.stereoseq_find_tissue_gef <- function(path) {
    feat_dir <- file.path(path, "feature_expression")
    if (!dir.exists(feat_dir)) {
        stop(sprintf(
            "[StereoSeq] feature_expression directory not found at: %s",
            feat_dir), call. = FALSE)
    }
    files <- list.files(feat_dir)
    hit   <- files[grep("\\.tissue\\.gef$", files)]
    if (length(hit) == 0L) {
        stop("[StereoSeq] No *.tissue.gef file found in feature_expression 
            directory",
            call. = FALSE)
    }
    file.path(feat_dir, hit[[1L]])
}

# Find *.adjusted.cellbin.gef under path/feature_expression/
.stereoseq_find_cellbin_gef <- function(path) {
    feat_dir <- file.path(path, "feature_expression")
    if (!dir.exists(feat_dir)) {
        stop(sprintf(
            "[StereoSeq] feature_expression directory not found at: %s",
            feat_dir), call. = FALSE)
    }
    files <- list.files(feat_dir)
    hit   <- files[grep("\\.adjusted\\.cellbin\\.gef$", files)]
    if (length(hit) == 0L) {
        stop("[StereoSeq] No *.adjusted.cellbin.gef file found in 
            feature_expression directory",
            call. = FALSE)
    }
    file.path(feat_dir, hit[[1L]])
}

# Find *_HE_mask.tif (exact, not the _edm_dis_ variant) under path/image/
.stereoseq_find_mask <- function(path) {
    img_dir <- file.path(path, "image")
    if (!dir.exists(img_dir)) return(NULL)
    files <- list.files(img_dir)
    # match *_HE_mask.tif but NOT *_HE_mask_edm_dis_*.tif
    hit <- files[grep("_HE_mask\\.tif$", files)]
    if (length(hit) == 0L) return(NULL)
    file.path(img_dir, hit[[1L]])
}

# Read a .gef file and return raw data as a named list.
# `what` controls which sub-datasets are loaded:
#   "expression" — expression counts + gene names (needed for matrix/binpoints)
#   "spatlocs"   — spatial coordinates only:
#                    bin type : reads the expression compound dataset 
#                               (x,y live there) but skips the gene table
#                    cell type: reads only the small `cell` table (id, x, y);
#                               skips the large geneExp table entirely
# Pass both values (default) when multiple outputs will be built from one read.
.stereoseq_read_gef <- function(path, type, bin_size, gene_column = NULL,
                                what = c("expression", "spatlocs"),
                                verbose = NULL) {
    package_check(pkg_name = "rhdf5", repository = "Bioc")
    gene_idx <- NULL  # data.table var

    need_expr <- "expression" %in% what
    need_spat <- "spatlocs" %in% what

    if (type == "bin") {
        vmsg(.v = verbose, "Reading tissue.gef file...")
        # bin spatial coordinates live inside the expression compound dataset,
        # so it must be read for either purpose
        exprDT <- data.table::setDT(rhdf5::h5read(
            file = path,
            name = paste0("geneExp/", bin_size, "/expression")
        ))
        geneDT <- NULL
        if (need_expr) {
            geneDT <- data.table::setDT(rhdf5::h5read(
                file = path,
                name = paste0("geneExp/", bin_size, "/gene")
            ))
            exprDT[, gene_idx := rep(seq_len(nrow(geneDT)), geneDT$count)]
        }
        vmsg(.v = verbose, "Finished reading tissue.gef")
        list(type = "bin", exprDT = exprDT, geneDT = geneDT)
    } else {
        vmsg(.v = verbose, "Reading adjusted.cellbin.gef file...")
        exprDT <- geneDT <- cellDT <- NULL
        if (need_expr) {
            exprDT <- data.table::setDT(rhdf5::h5read(
                file = path, name = "cellBin/geneExp"
            ))
            geneDT <- data.table::setDT(rhdf5::h5read(
                file = path, name = "cellBin/gene"
            ))
            exprDT[, gene_idx := rep(seq_len(nrow(geneDT)), geneDT$cellCount)]
        }
        # cell table is small and needed both for spatlocs and for expression
        # matrix column ordering, so read it whenever either is requested
        if (need_expr || need_spat) {
            cellDT <- data.table::setDT(rhdf5::h5read(
                file = path, name = "cellBin/cell"
            ))
        }
        vmsg(.v = verbose, "Finished reading adjusted.cellbin.gef")
        list(type = "cell", exprDT = exprDT, geneDT = geneDT, cellDT = cellDT)
    }
}

# Build a giottoExprObj (wrapped in a list) from already-read gef data.
# Uses Matrix::sparseMatrix() directly from triplet data — avoids the
# memory-hungry dcast + Reduce(merge) approach entirely.
.stereoseq_build_expression <- function(
        gef_data, type, gene_column, spat_unit, verbose = NULL) {
    bin_ID <- gene_idx <- count <- x <- y <- cellID <- NULL  # data.table vars

    vmsg(.v = verbose, "Building expression matrix...")

    if (type == "bin") {
        exprDT <- gef_data$exprDT

        # assign integer bin IDs from unique (x, y) positions
        # use pre-computed bins from gef_data when available 
        # (set by gobject_fun)
        bins <- if (!is.null(gef_data$bins)) gef_data$bins else {
            b <- unique(exprDT[, c("x", "y")], by = c("x", "y"))
            b[, bin_ID := .I]
            b
        }
        dt <- merge(exprDT, bins, by = c("x", "y"))

        geneDT <- gef_data$geneDT
        all_names <- as.character(geneDT[[gene_column]]) # O(n_genes) string 
        # conversion derive gene_ids only from expressed genes (same behavior 
        # as the old sort(unique(exprDT$genes)) — zero-count genes are 
        # not included)
        expressed_names <- all_names[sort(unique(dt$gene_idx))]
        gene_ids    <- sort(unique(expressed_names))
        # map each geneDT index -> matrix row (two genes sharing a name 
        # share a row; non-expressed genes map to NA but their indices never 
        # appear in dt$gene_idx)
        name_to_row <- match(all_names, gene_ids)
        n_bins   <- nrow(bins)

        expMatrix <- Matrix::sparseMatrix(
            i = name_to_row[dt$gene_idx],
            j = dt$bin_ID,
            x = as.integer(dt$count),
            dims = c(length(gene_ids), n_bins),
            dimnames = list(gene_ids, paste0("bin_", seq_len(n_bins)))
        )
    } else {
        exprDT <- gef_data$exprDT
        cellDT <- gef_data$cellDT

        geneDT    <- gef_data$geneDT
        all_names <- as.character(geneDT[[gene_column]]) # O(n_genes) string 
        # conversion derive gene_ids only from expressed genes (same behavior 
        # as the old sort(unique(exprDT$genes)) — zero-count genes are not 
        # included)
        expressed_names <- all_names[sort(unique(exprDT$gene_idx))]
        gene_ids <- sort(unique(expressed_names))
        # map each geneDT index -> matrix row (non-expressed genes map to NA but
        # their indices never appear in exprDT$gene_idx)
        name_to_row <- match(all_names, gene_ids)

        cell_ids <- paste0("cell_", cellDT$id)
        # map raw cellID to 1-based column index
        cell_idx_map <- data.table::setattr(
            seq_len(nrow(cellDT)), "names", as.character(cellDT$id)
        )

        expMatrix <- Matrix::sparseMatrix(
            i = name_to_row[exprDT$gene_idx],
            j = cell_idx_map[as.character(exprDT$cellID)],
            x = as.integer(exprDT$count),
            dims = c(length(gene_ids), length(cell_ids)),
            dimnames = list(gene_ids, cell_ids)
        )
    }

    vmsg(.v = verbose, "Finished expression matrix")

    expr_obj <- createExprObj(
        expression_data = expMatrix,
        name = "raw",
        spat_unit = spat_unit,
        feat_type = "rna",
        provenance = spat_unit
    )
    list(expr_obj)
}

# Build a giottoBinPoints from already-read gef data.
# Builds a transient sparse matrix via Matrix::sparseMatrix() and a spatLocsObj,
# then delegates all internal slot construction to createGiottoBinPoints().
.stereoseq_build_binpoints <- function(
        gef_data, type, gene_column, negative_y, spat_unit,
        verbose = NULL) {
    vmsg(.v = verbose, "Building giottoBinPoints...")

    # Reuse the expression and spatlocs builders — the transient sparse matrix
    # is small (only non-zero entries) and createGiottoBinPoints() owns the
    # internal slot construction, keeping this function maintainable.
    expr_list <- .stereoseq_build_expression(
        gef_data = gef_data,
        type = type,
        gene_column = gene_column,
        spat_unit = spat_unit,
        verbose = NULL
    )
    sl <- .stereoseq_build_spatlocs(
        gef_data = gef_data,
        type = type,
        negative_y = negative_y,
        spat_unit = spat_unit,
        verbose = verbose
    )

    gbp <- createGiottoBinPoints(
        expr_values = expr_list[[1L]],
        spatial_locs = sl,
        feat_type = "rna"
    )
    spatUnit(gbp) <- spat_unit

    vmsg(.v = verbose, "Finished giottoBinPoints")
    gbp
}

# Build a spatLocsObj from already-read gef data.
.stereoseq_build_spatlocs <- function(
        gef_data, type, negative_y, spat_unit, verbose = NULL) {
    x <- y <- bin_ID <- cell_ID <- NULL  # data.table vars

    vmsg(.v = verbose, "Building spatial locations...")

    if (type == "bin") {
        exprDT <- gef_data$exprDT
        # use pre-computed bins from gef_data when available (set by 
        # gobject_fun);
        # copy() is required because `:=` below mutates in place (data.table 
        # semantics)
        spat_locs <- if (!is.null(gef_data$bins)) {
            data.table::copy(gef_data$bins)
        } else {
            b <- unique(exprDT[, c("x", "y")], by = c("x", "y"))
            b[, bin_ID := seq_len(nrow(b))]
            b
        }
        spat_locs[, cell_ID := paste0("bin_", bin_ID)]
        spat_locs <- spat_locs[, .(cell_ID, x, y)]
        vmsg(.v = verbose, nrow(spat_locs), " bins in total")
    } else {
        cellDT <- gef_data$cellDT
        cellDT[, cell_ID := paste0("cell_", id)]
        spat_locs <- cellDT[, .(cell_ID, x, y)]
        vmsg(.v = verbose, nrow(spat_locs), " cells in total")
    }

    spat_locs[, x := as.integer(x)]
    spat_locs[, y := as.integer(y)]

    if (isTRUE(negative_y)) {
        spat_locs[, y := 0L - y]
    }

    sl <- createSpatLocsObj(
        coordinates = spat_locs,
        name = "raw",
        spat_unit = spat_unit,
        provenance = spat_unit,
        verbose = FALSE
    )
    vmsg(.v = verbose, "Finished spatial locations")
    sl
}

# Public-facing expression loader: reads gef and returns list of giottoExprObj.
.stereoseq_expression <- function(
        path, type, bin_size, gene_column,
        spat_unit, verbose = NULL) {
    gef_data <- .stereoseq_read_gef(
        path = path,
        type = type,
        bin_size = bin_size,
        gene_column = gene_column,
        what = "expression",
        verbose = verbose
    )
    .stereoseq_build_expression(
        gef_data = gef_data,
        type = type,
        gene_column = gene_column,
        spat_unit = spat_unit,
        verbose = verbose
    )
}

# Public-facing spatlocs loader: reads gef and returns a spatLocsObj.
.stereoseq_spatlocs <- function(path, type, bin_size, gene_column,
                                negative_y, spat_unit, verbose = NULL) {
    gef_data <- .stereoseq_read_gef(
        path = path,
        type = type,
        bin_size = bin_size,
        gene_column = gene_column,
        what = "spatlocs",
        verbose = verbose
    )
    .stereoseq_build_spatlocs(
        gef_data = gef_data,
        type = type,
        negative_y = negative_y,
        spat_unit = spat_unit,
        verbose = verbose
    )
}

# Public-facing binpoints loader: reads gef and returns a giottoBinPoints.
# Most memory-efficient representation — no matrix is created.
.stereoseq_binpoints <- function(
        path, type, bin_size, gene_column,
        negative_y, spat_unit, verbose = NULL) {
    gef_data <- .stereoseq_read_gef(
        path = path,
        type = type,
        bin_size = bin_size,
        gene_column = gene_column,
        what = c("expression", "spatlocs"),
        verbose = verbose
    )
    .stereoseq_build_binpoints(
        gef_data = gef_data,
        type = type,
        gene_column = gene_column,
        negative_y = negative_y,
        spat_unit = spat_unit,
        verbose = verbose
    )
}

# Load the H&E registered image from the image directory.
# Returns a giottoLargeImage or NULL if no image is found.
.stereoseq_image <- function(path, negative_y = TRUE, verbose = NULL) {
    vmsg(.v = verbose, "Attaching HE image...")

    he_image_path <- list.files(
        path = path, pattern = "_regist", full.names = TRUE)

    if (length(he_image_path) == 0L) {
        warning("[StereoSeq] No *_regist image found in image directory. Skipping.",
            call. = FALSE)
        return(NULL)
    }

    gimg <- createGiottoLargeImage(
        he_image_path[[1L]],
        name = "image",
        negative_y = negative_y
    )
    vmsg(.v = verbose, "Finished attaching image")
    gimg
}


# Load the mask image and create cell polygons.
# `path` must be the full filepath to a *_HE_mask.tif file.
# Returns a giottoPolygon or NULL if path is NULL / file not found.
# When negative_y = TRUE the polygon y-coordinates are flipped (terra::flip)
# to match the negated-y convention used for spatial locations from the gef file
.stereoseq_mask <- function(path, negative_y = TRUE, verbose = NULL) {
    if (is.null(path) || !file.exists(path)) {
        warning(
            "[StereoSeq] No *_HE_mask.tif file found. Skipping mask polygons.",
            call. = FALSE)
        return(NULL)
    }
    vmsg(.v = verbose, "Creating polygons from mask...")
    poly <- createGiottoPolygonsFromMask(
        maskfile = path,
        calc_centroids = TRUE
    )
    if (isTRUE(negative_y)) {
        # The mask polygon comes out of terra::as.polygons in raster convention:
        # y = 0 at the BOTTOM of the image, y = nrows at the TOP.
        # The gef spatial locations use image convention (y = 0 at top) and
        # are negated (0 - y), placing them at y ∈ [-nrows, 0] with 0 at top.
        # Shifting the polygon down by nrows converts it to the same convention
        # without flipping orientation.
        # Concretely: a polygon vertex at terra-y T corresponds to image-row
        # (nrows - T), whose negated gef value is -(nrows - T) = T - nrows.
        #
        # IMPORTANT: must use the full image height (nrows of the source raster)
        # NOT the polygon bbox ymax. The tissue typically covers only a 
        # sub-region of the full slide image, so the polygon ymax < nrows. 
        # Using the polygon ymax would apply the wrong shift and misalign 
        # polygons with spatlocs.
        nrows <- nrow(terra::rast(path))
        poly <- spatShift(poly, dy = -nrows)
    }
    vmsg(.v = verbose, "Finished creating polygons from mask")
    poly
}


# -------------------------------------------------------------------------
# createGiottoStereoSeqObjectBin
# -------------------------------------------------------------------------

#' @title Create Stereo-seq Giotto Object from Square Bin Data
#' @name createGiottoStereoSeqObjectBin
#' @description Convenience function to create a Giotto object from a
#' Stereo-seq `outs` directory using square bin (squarebin) expression data.
#' Point `stereoseq_dir` directly at the `outs` directory produced by the SAW
#' pipeline. For lower-level loading of individual pieces of data, see
#' [importStereoSeq()].
#' @param stereoseq_dir filepath to the Stereo-seq `outs` directory.
#' @param bin_size character. Bin size level to load. One of `"bin1"`,
#'   `"bin5"`, `"bin10"`, `"bin20"`, `"bin50"`, `"bin100"` (default),
#'   `"bin150"`, `"bin200"`. List available bins with
#'   `rhdf5::h5ls("path/to/*.tissue.gef")`.
#' @param gene_column character. Feature identifier column. One of
#'   `"geneName"` (default) or `"geneID"`.
#' @param negative_y logical. Map data to negative y spatial values
#'   (default `TRUE`). Origin is placed at the upper-left.
#' @param load_expression logical. Whether to load the expression matrix.
#'   Uses `Matrix::sparseMatrix()` directly from the gef triplet data for
#'   memory efficiency. Set to `FALSE` to skip.
#' @param load_spatlocs logical. Whether to load spatial locations.
#' @param load_binpoints logical (default `FALSE`). Whether to load a
#'   [giottoBinPoints-class] object — the most memory-efficient representation.
#'   Data stays as integer triplets + `SpatVector`; no matrix is created.
#'   Useful for very large/fine bin datasets or as input to
#'   [calculateOverlap()] with custom polygons.
#' @param load_image logical. Whether to load the H&E registered image.
#' @param load_mask logical (default `TRUE`). Whether to create cell polygons
#'   from the `*_HE_mask.tif` file in `stereoseq_dir/image/`. Uses
#'   [createGiottoPolygonsFromMask()] with `calc_centroids = TRUE`.
#' @param gef_path (optional) direct filepath to the `*.tissue.gef` file.
#'   Auto-detected from `stereoseq_dir` when not provided.
#' @param image_path (optional) filepath or directory for the image.
#'   Auto-detected from `stereoseq_dir/image/` when not provided.
#' @param mask_path (optional) direct filepath to the `*_HE_mask.tif` file.
#'   Auto-detected from `stereoseq_dir/image/` when not provided.
#' @param instructions giotto instructions to apply.
#' @param verbose verbosity
#' @returns giotto object
#' @examples
#' if (FALSE) {
#' g <- createGiottoStereoSeqObjectBin(
#'     stereoseq_dir = "path/to/outs",
#'     bin_size = "bin100"
#' )
#' }
#' @seealso [importStereoSeq()] [createGiottoStereoSeqObjectCell()]
#' @export
createGiottoStereoSeqObjectBin <- function(
    stereoseq_dir,
    bin_size = "bin100",
    gene_column = c("geneName", "geneID"),
    negative_y = TRUE,
    load_expression = TRUE,
    load_spatlocs = TRUE,
    load_binpoints = FALSE,
    load_image = TRUE,
    load_mask = TRUE,
    gef_path = NULL,
    image_path = NULL,
    mask_path = NULL,
    instructions = NULL,
    verbose = NULL) {

    reader <- importStereoSeq(
        stereoseq_dir = stereoseq_dir,
        type = "bin",
        bin_size = bin_size,
        gene_column = match.arg(gene_column, c("geneName", "geneID")),
        negative_y = negative_y
    )

    read_args <- list(
        load_expression = load_expression,
        load_spatlocs = load_spatlocs,
        load_binpoints = load_binpoints,
        load_image = load_image,
        load_mask = load_mask,
        instructions = instructions,
        verbose = verbose
    )

    if (!is.null(gef_path)) read_args$gef_path <- gef_path
    if (!is.null(image_path)) read_args$image_path <- image_path
    if (!is.null(mask_path)) read_args$mask_path <- mask_path

    do.call(reader$create_gobject, read_args)
}


# -------------------------------------------------------------------------
# createGiottoStereoSeqObjectCell
# -------------------------------------------------------------------------

#' @title Create Stereo-seq Giotto Object from Cell Bin Data
#' @name createGiottoStereoSeqObjectCell
#' @description Convenience function to create a Giotto object from a
#' Stereo-seq `outs` directory using cell bin (cellbin) expression data.
#' Point `stereoseq_dir` directly at the `outs` directory produced by the SAW
#' pipeline. For lower-level loading of individual pieces of data, see
#' [importStereoSeq()].
#' @param stereoseq_dir filepath to the Stereo-seq `outs` directory.
#' @param gene_column character. Feature identifier column. One of
#'   `"geneName"` (default) or `"geneID"`.
#' @param negative_y logical. Map data to negative y spatial values
#'   (default `TRUE`). Origin is placed at the upper-left.
#' @param load_expression logical. Whether to load the expression matrix.
#'   Uses `Matrix::sparseMatrix()` directly from the gef triplet data for
#'   memory efficiency. Set to `FALSE` to skip.
#' @param load_spatlocs logical. Whether to load spatial locations.
#' @param load_binpoints logical (default `FALSE`). Whether to also load a
#'   [giottoBinPoints-class] object — the most memory-efficient representation.
#'   Data stays as integer triplets + `SpatVector`; no matrix is created.
#' @param load_image logical. Whether to load the H&E registered image.
#' @param load_mask logical (default `TRUE`). Whether to create cell polygons
#'   from the `*_HE_mask.tif` file in `stereoseq_dir/image/`. Uses
#'   [createGiottoPolygonsFromMask()] with `calc_centroids = TRUE`.
#' @param gef_path (optional) direct filepath to the `*.adjusted.cellbin.gef`
#'   file. Auto-detected from `stereoseq_dir` when not provided.
#' @param image_path (optional) filepath or directory for the image.
#'   Auto-detected from `stereoseq_dir/image/` when not provided.
#' @param mask_path (optional) direct filepath to the `*_HE_mask.tif` file.
#'   Auto-detected from `stereoseq_dir/image/` when not provided.
#' @param instructions giotto instructions to apply.
#' @param verbose verbosity
#' @returns giotto object
#' @examples
#' if (FALSE) {
#' g <- createGiottoStereoSeqObjectCell(
#'     stereoseq_dir = "path/to/outs"
#' )
#' }
#' @seealso [importStereoSeq()] [createGiottoStereoSeqObjectBin()]
#' @export
createGiottoStereoSeqObjectCell <- function(
    stereoseq_dir,
    gene_column = c("geneName", "geneID"),
    negative_y = TRUE,
    load_expression = TRUE,
    load_spatlocs = TRUE,
    load_binpoints = FALSE,
    load_image = TRUE,
    load_mask = TRUE,
    gef_path = NULL,
    image_path = NULL,
    mask_path = NULL,
    instructions = NULL,
    verbose = NULL) {

    reader <- importStereoSeq(
        stereoseq_dir = stereoseq_dir,
        type = "cell",
        gene_column = match.arg(gene_column, c("geneName", "geneID")),
        negative_y = negative_y
    )

    read_args <- list(
        load_expression = load_expression,
        load_spatlocs = load_spatlocs,
        load_binpoints = load_binpoints,
        load_image = load_image,
        load_mask = load_mask,
        instructions = instructions,
        verbose = verbose
    )

    if (!is.null(gef_path)) read_args$gef_path <- gef_path
    if (!is.null(image_path)) read_args$image_path <- image_path
    if (!is.null(mask_path)) read_args$mask_path <- mask_path

    do.call(reader$create_gobject, read_args)
}



# -------------------------------------------------------------------------
# createGiottoStereoSeqObject (deprecated)
# -------------------------------------------------------------------------

#' Create Stereo-seq Giotto Object
#'
#' @param stereoseq_dir filepath to the exported Stereo-seq directory.
#' @param type character. Use "squarebin" to read expression at different bin
#' levels (default), "cellbin" to read expression at cell resolution, or
#' "subcellular" to read individual transcripts and cell boundaries.
#' @param bin_size bin size. Choose a value from "bin1", "bin5", "bin10",
#' "bin20", "bin50", "bin100", "bin150", or "bin200". List the available bins
#' in your sample using rhdf5::h5ls("path/to/*.tissue.gef").
#' Only needed when using type = "squarebin".
#' @param gene_column (optional) character. Which column contains the gene names
#' within the geneExp information. Choose from "geneName" (default), or "geneID"
#' @param verbose logical. Be verbose.
#' @param h5_file (optional) name to create an on-disk HDF5 file.
#' @param instructions list of instructions or output result
#' @param negative_y logical. Map data to negative y spatial values during
#' automatic alignment (Default = TRUE). Meaning that origin is in upper left
#' instead of lower left.
#' from \code{\link{createGiottoInstructions}}
#' @returns Giotto Stereo-seq object
#' @details This function is deprecated. Use [createGiottoStereoSeqObjectBin()]
#' for square bin data or [createGiottoStereoSeqObjectCell()] for cell bin data.
#' The `"subcellular"` type is still only available through this function.
#' @seealso [importStereoSeq()] [createGiottoStereoSeqObjectBin()]
#' [createGiottoStereoSeqObjectCell()]
#' @export
createGiottoStereoSeqObject <- function(
        stereoseq_dir,
        type = "squarebin",
        bin_size = "bin100",
        gene_column = "geneName",
        negative_y = TRUE,
        shift_polygon_y = 0,
        verbose = TRUE,
        h5_file = NULL,
        instructions = NULL) {

    if (type != "subcellular") {
        .Deprecated(
            msg = paste0(
                "'createGiottoStereoSeqObject' is deprecated.\n",
                "Use 'createGiottoStereoSeqObjectBin()' for square bin outputs or\n",
                "'createGiottoStereoSeqObjectCell()' for cell bin outputs."
            )
        )
    }

    # data.table vars
    genes <- gene_idx <- x <- y <- sdimx <- sdimy <- cell_ID <- bin_ID <-
        count <- i.bin_ID <- NULL

    # package check
    package_check(pkg_name = "rhdf5", repository = "Bioc")

    # directory check
    if (!file.exists(stereoseq_dir)) stop(
        "Path to Stereo-seq directory does not exist")

    dir_files <- list.files(
        file.path(stereoseq_dir, "outs", "feature_expression"))

    # file reading type check
    if(!type %in% c("squarebin", "cellbin", "subcellular")) stop(
        "'type' should be either 'squarebin', 'cellbin', or 'subcellular'")

    # Read squarebin
    if(type == "squarebin") {
        expression_file <- file.path(
            stereoseq_dir, "outs", "feature_expression",
            dir_files[grep(".tissue.gef", dir_files)])

        if (!file.exists(expression_file)) stop(
            "Path to expression file ",
            file.path(
                stereoseq_dir, "outs", "feature_expression", ".tissue.gef"),
            " does not exist")

        # check if proper bin_size is selected. These are determined in SAW pipeline
        vmsg(.v = verbose, "Reading expression file... \n")

        bin_size_options <- c("bin1", "bin5", "bin10",
                              "bin20", "bin50", "bin100",
                              "bin150", "bin200")
        if (!(bin_size %in% bin_size_options)) {
            stop("Please select valid bin size, see
            ?createGiottoStereoSeqObject for details.")
        }

        # check valid gene_column value
        if(!gene_column %in% c("geneName", "geneID")) stop(
            "Provide a valid 'gene_column' value. It should be 'geneName' or
            'geneID'")

        # 1. read tissue.gef file at specific bin size
        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = paste0("geneExp/", bin_size)
        )

        exprDT <- data.table::setDT(geneExpData[["expression"]])
        geneDT <- data.table::setDT(geneExpData[["gene"]])
        exprDT[, `:=`(genes, rep(x = geneDT[[gene_column]], geneDT$count))]

        vmsg(.v = verbose, "Finished reading in tissue.gef")
    }

    # Read cellbin
    if(type == "cellbin") {
        expression_file <- file.path(
            stereoseq_dir, "outs", "feature_expression",
            dir_files[grep(".adjusted.cellbin.gef", dir_files)])

        if (!file.exists(expression_file)) stop(
            "Path to expression file ",
            file.path(
                stereoseq_dir, "outs", "feature_expression", ".adjusted.cellbin.gef"),
            " does not exist")

        # 1. read .adjusted.cellbin.gef file
        vmsg(.v = verbose, "Reading expression file... \n")

        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = "cellBin"
        )

        exprDT <- data.table::setDT(geneExpData[["geneExp"]])
        geneDT <- data.table::setDT(geneExpData[["gene"]])

        # check valid gene_column value
        if(!gene_column %in% c("geneName", "geneID")) stop(
            "Provide a valid 'gene_column' value. It should be 'geneName' or
            'geneID'")

        exprDT[, `:=`(genes, rep(x = geneDT[[gene_column]], geneDT$cellCount))]

        vmsg(.v = verbose, "Finished reading in adjusted.cellbin.gef")
    }

    # Read subcellular
    if(type == "subcellular") {

        vmsg(.v = verbose, "Reading transcripts file ...\n")

        expression_file <- file.path(
            stereoseq_dir, "outs", "feature_expression",
            dir_files[grep(".adjusted.cellbin.gef", dir_files)])

        if (!file.exists(expression_file)) stop(
            "Path to file ",
            file.path(
                stereoseq_dir, "outs", "feature_expression",
                ".adjusted.cellbin.gef"),
            " does not exist")

        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = "cellBin"
        )
        geneDT <- data.table::setDT(geneExpData[["gene"]])

        if (!file.exists(expression_file)) stop(
            "Path to transcripts file ",
            file.path(
                stereoseq_dir, "outs", "feature_expression",
                "_raw_barcode_gene_exp.txt"),
            " does not exist")

        transcripts_file <- file.path(
            stereoseq_dir, "outs", "feature_expression",
            dir_files[grep("_raw_barcode_gene_exp.txt", dir_files)])

        # read transcripts file with gene locations
        transcripts <- data.table::fread(transcripts_file)

        vmsg(.v = verbose, "Reading cell boundaries...\n")

        files_list <- list.files(file.path(stereoseq_dir, "outs", "image"))
        mask_file <- files_list[grep("_HE_mask.tif$", files_list)]

        # Check that mask file exists
        if(is.null(mask_file)) {
            stop("The expected mask file was not found in the image directory.
                 Make sure a *_HE_mask.tif files is located under outs/image/")
        }

        # Read polygons from mask
        mask_poly <- createGiottoPolygonsFromMask(
            maskfile = file.path(stereoseq_dir, "outs", "image", mask_file),
            calc_centroids = TRUE
        )

    }

    # 2. create spatial locations
    vmsg(.v = verbose, "Creating spatial_locations... \n")

    if(type == "squarebin") {
        spatial_locations <- unique(exprDT[, c("x", "y")], by = c("x", "y"))
        spatial_locations[, bin_ID := seq_len(nrow(spatial_locations))]

        vmsg(.v = verbose, nrow(spatial_locations), " bins in total \n")
    }

    if(type == "cellbin") {
        cellDT <- data.table::setDT(geneExpData[["cell"]])
        cellDT[, cell_ID := paste0("cell_", id)]
        spatial_locations <- cellDT[, .(cell_ID, x, y)]

        vmsg(.v = verbose, nrow(spatial_locations), " cells in total \n")
    }

    if(type == "subcellular") {

        if(isTRUE(negative_y)) {
            mask_poly <- terra::flip(mask_poly)
        }

        spatial_locations <- as.data.frame(
            terra::geom(terra::centroids(mask_poly)))
        spatial_locations$cell_ID <- paste0("cell_", spatial_locations$geom)
        spatial_locations <- spatial_locations[, c("x", "y", "cell_ID")]

        spat_locs <- createSpatLocsObj(coordinates = spatial_locations,
                                       spat_unit = "cell",
                                       name = "raw")

        vmsg(.v = verbose, nrow(spatial_locations),
            " polygons (cells) in total \n")

        # merge with geneID column
        transcript_locs <- merge(transcripts, geneDT[, .(geneID, geneName)])
        transcript_locs <- as.data.frame(transcript_locs)[
            , c("x", "y", gene_column, "readCount")]
        colnames(transcript_locs)[3] <- "feat_ID"
        colnames(transcript_locs)[4] <- "count"

        if(isTRUE(negative_y)) {
            transcript_locs$y <- 0 - transcript_locs$y
        } else {
            max_y <- max(transcript_locs$y)
            middle_point <- max_y/2
            transcript_locs$y <- 2*middle_point - transcript_locs$y
        }

        # Create giotto points
        g_points <- createGiottoPoints(x = transcript_locs)
        g_points <- crop(g_points, terra::ext(mask_poly))

        vmsg(.v = verbose, nrow(transcript_locs), " transcripts in total \n")
    }

    vmsg(.v = verbose, "Finished spatial_locations \n")

    # Create expression matrix
    if(type %in% c("squarebin", "cellbin")) {
        vmsg(.v = verbose, "Creating expression matrix... \n")

        if(type == "squarebin") {
            exprDT <- merge(exprDT, spatial_locations, by = c("x", "y"))
            exprDT <- exprDT[, .(bin_ID, genes, count)]
            exprDT[, count := as.integer(count)]
            exprDT[, bin_ID := as.integer(bin_ID)]

            genes_in_exprDT <- sort(unique(exprDT$genes))
            genes_in_exprDT_chunks <- split(genes_in_exprDT, ceiling(seq_along(genes_in_exprDT)/1000))

            w_list <- list()
            w_list <- lapply(genes_in_exprDT_chunks, function(x)
                exprDT[genes %in% x])
            w_list <- lapply(w_list , function(z)
                data.table::dcast(
                    z,
                    bin_ID ~ genes,
                    value.var = 'count',
                    fun.aggregate = sum))

            expMatrix <- Reduce(function( ... ) merge(
                ... ,
                by = "bin_ID",
                all = TRUE),
                w_list)

            rownames_matrix <- paste0("bin_", expMatrix$bin_ID)
            expMatrix <- as.matrix(expMatrix[, -1, with = FALSE])
            rownames(expMatrix) <- rownames_matrix
            expMatrix[is.na(expMatrix)] <- 0
            expMatrix <- Matrix::Matrix(expMatrix, sparse = TRUE)
            expMatrix <- t(expMatrix)

            spatial_locations[, cell_ID := paste0("bin_", bin_ID)]
            spatial_locations <- spatial_locations[, c("x", "y", "cell_ID")]
        }

        if(type == "cellbin") {
            exprDT <- exprDT[, .(cellID, genes, count)]
            exprDT[, cellID := as.integer(cellID)]
            exprDT[, count := as.integer(count)]

            genes_in_exprDT <- sort(unique(exprDT$genes))
            genes_in_exprDT_chunks <- split(genes_in_exprDT, ceiling(seq_along(genes_in_exprDT)/1000))

            w_list <- list()
            w_list <- lapply(genes_in_exprDT_chunks, function(x)
                exprDT[genes %in% x])
            w_list <- lapply(w_list , function(z)
                data.table::dcast(
                    z,
                    cellID ~ genes,
                    value.var = 'count',
                    fun.aggregate = sum))

            expMatrix <- Reduce(function( ... ) merge(
                ... ,
                by = "cellID",
                all = TRUE),
                w_list)

            rownames_matrix <- paste0("cell_", expMatrix$cellID)
            expMatrix <- as.matrix(expMatrix[, -1, with = FALSE])
            rownames(expMatrix) <- rownames_matrix
            expMatrix[is.na(expMatrix)] <- 0
            expMatrix <- Matrix::Matrix(expMatrix, sparse = TRUE)
            expMatrix <- t(expMatrix)
        }

        vmsg(.v = verbose, "finished expression matrix")
    }

    # Create Giotto object
    vmsg(.v = verbose, "Creating giotto object... \n")

    if(type %in% c("squarebin", "cellbin")) {
        # ensure first non-numerical col is cell_ID
        spatial_locations[, x := as.integer(x)]
        spatial_locations[, y := as.integer(y)]


        if(isTRUE(negative_y)) {
            spatial_locations[, y := 0 - y]
        }

        stereo <- createGiottoObject(
            expression = expMatrix,
            spatial_locs = spatial_locations,
            verbose = verbose,
            h5_file = h5_file,
            instructions = instructions
        )
    }

    if(type == "subcellular") {

        stereo <- giotto(instructions = instructions)

        # Add giotto points
        stereo <- setGiotto(stereo, g_points)

        # Add giotto polygons
        stereo <- setGiotto(stereo, mask_poly)

        # Add polygon centroids
        stereo <- setGiotto(stereo, spat_locs)

    }

    # 5. add image
    vmsg(.v = verbose, "Attaching HE image... \n")

    image_dir <- file.path(stereoseq_dir, "outs", "image")
    he_image_path <- list.files(
        path = image_dir, pattern = "_regist", full.names = TRUE)
    gimg <- createGiottoLargeImage(he_image_path,
                                   name = "image",
                                   negative_y = negative_y)

    stereo <- addGiottoLargeImage(
        gobject = stereo,
        largeImages = gimg,
        negative_y = negative_y,
        verbose = verbose
    )

    vmsg(.v = verbose, "Finished giotto object... \n")
    return(stereo)
}

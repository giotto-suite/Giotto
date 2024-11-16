#' @title normalizeGiotto
#' @name normalizeGiotto
#' @description fast normalize and/or scale expression values of Giotto object
#' @param gobject `giotto` object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param norm_methods normalization method to use
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param log_norm transform values to log-scale
#' @param log_offset offset value to add to expression matrix, default = 1
#' @param logbase log base to use to log normalize expression values
#' @param scale_feats z-score genes over all cells
#' @param scale_cells z-score cells over all genes
#' @param scale_order order to scale feats and cells
#' @param theta theta parameter for the pearson residual normalization step
#' @param name character. name to use for normalization results
#' @param verbose be verbose
#' @param scale_genes deprecated, use scale_feats
#' @param update_slot deprecated. Use `name` param instead
#' @md
#' @returns `giotto` object
#' @details Currently there are two 'methods' to normalize your raw counts data.
#'
#' A. The standard method follows the standard protocol which can be adjusted
#' using the provided parameters and follows the following order: \cr
#' \itemize{
#'   \item{1. Data normalization for total library size and scaling by a custom
#'   scale-factor.}
#'   \item{2. Log transformation of data.}
#'   \item{3. Z-scoring of data by genes and/or cells.}
#' }
#' B. The normalization method as provided by the osmFISH paper is also
#' implemented: \cr
#' \itemize{
#'   \item{1. First normalize genes, for each gene divide the counts by the
#'   total gene count and multiply by the total number of genes.}
#'   \item{2. Next normalize cells, for each cell divide the normalized gene
#'   counts by the total counts per cell and multiply by the total number of
#'   cells.}
#' }
#' C. The normalization method as provided by Lause/Kobak et al is also
#' implemented: \cr
#' \itemize{
#'   \item{1. First calculate expected values based on Pearson correlations.}
#'   \item{2. Next calculate z-scores based on observed and expected values.}
#' }
#' D. Quantile normalization across features
#' \itemize{
#'   \item{1. Rank feature expression}
#'   \item{2. Define a common distribution by sorting expression values per
#'   feature then finding the mean across all features per index}
#'   \item{3. Apply common distribution to expression information by using
#'   the ranks from step 1 as indices}
#' }
#' By default the latter two results will be saved in the Giotto slot for
#' scaled expression, this can be changed by changing the update_slot parameters
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' normalizeGiotto(g) # default is method A
#' @export
normalizeGiotto <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = "raw",
    norm_methods = c("standard", "pearson_resid", "osmFISH", "quantile"),
    library_size_norm = TRUE,
    scalefactor = 6e3,
    log_norm = TRUE,
    log_offset = 1,
    logbase = 2,
    scale_feats = TRUE,
    scale_genes = deprecated(),
    scale_cells = TRUE,
    scale_order = c("first_feats", "first_cells"),
    theta = 100,
    name = "scaled",
    update_slot = deprecated(),
    verbose = TRUE) {
    ## deprecated arguments
    scale_feats <- deprecate_param(
        scale_genes, scale_feats,
        fun = "normalizeGiotto",
        when = "3.0.0"
    )
    name <- deprecate_param(
        update_slot, name,
        fun = "normalizeGiotto",
        when = "4.1.3"
    )

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    ## default is to start from raw data
    values <- match.arg(expression_values, unique(c("raw", expression_values)))
    raw_expr <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )

    norm_methods <- match.arg(
        arg = norm_methods, choices = c(
            "standard", "pearson_resid", "osmFISH", "quantile"
        )
    )

    # normalization according to standard methods
    gobject <- switch(norm_methods,
        "standard" = .rna_standard_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            library_size_norm = library_size_norm,
            scalefactor = scalefactor,
            log_norm = log_norm,
            log_offset = log_offset,
            logbase = logbase,
            scale_feats = scale_feats,
            scale_cells = scale_cells,
            scale_order = scale_order,
            verbose = verbose
        ),
        "osmFISH" = .rna_osmfish_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            name = name,
            verbose = verbose
        ),
        "pearson_resid" = .rna_pears_resid_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            theta = theta,
            name = name,
            verbose = verbose
        ),
        "quantile" = .quantile_norm(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            name = name,
            verbose = verbose
        )
    )

    ## update parameters used ##

    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try(
        {
            upstream_func <- sys.call(-2)
            fname <- as.character(upstream_func[[1]])
            if (fname == "processGiotto") {
                return(gobject)
            }
        },
        silent = TRUE
    )


    # If this function call is not downstream of processGiotto, update normally
    gobject <- update_giotto_params(gobject, description = "_normalize")

    return(gobject)
}



# internals ####


#' @title Normalize expression matrix for library size
#' @param mymatrix matrix object
#' @param scalefactor scalefactor
#' @returns matrix
#' @keywords internal
#' @noRd
.lib_norm_giotto <- function(mymatrix, scalefactor) {
    libsizes <- colSums_flex(mymatrix)

    if (0 %in% libsizes) {
        warning(wrap_txt("Total library size or counts for individual spat
                    units are 0.
                    This will likely result in normalization problems.
                    filter (filterGiotto) or impute (imputeGiotto) spatial
                    units."))
    }

    norm_expr <- t_flex(t_flex(mymatrix) / libsizes) * scalefactor
    return(norm_expr)
}

#' @title Log normalize expression matrix
#' @returns matrix
#' @keywords internal
#' @noRd
.log_norm_giotto <- function(mymatrix, base, offset) {
    if (methods::is(mymatrix, "DelayedArray")) {
        mymatrix <- log(mymatrix + offset) / log(base)
        # } else if(methods::is(mymatrix, 'DelayedMatrix')) {
        #   mymatrix = log(mymatrix + offset)/log(base)
    } else if (methods::is(mymatrix, "dgCMatrix")) {
        mymatrix@x <- log(mymatrix@x + offset) / log(base)
        # replace with sparseMatrixStats
    } else if (methods::is(mymatrix, "Matrix")) {
        mymatrix@x <- log(mymatrix@x + offset) / log(base)
    } else if (methods::is(mymatrix, "dbMatrix")) {
        mymatrix[] <- dplyr::mutate(mymatrix[], x = x + offset)
        # workaround for lack of @x slot
        mymatrix <- log(mymatrix) / log(base)
    } else {
        mymatrix <- log(as.matrix(mymatrix) + offset) / log(base)
    }

    return(mymatrix)
}


#' @title compute_dbMatrix
#' @description saves dbMatrix to db if global option is set
#' @details
#' Set \code{options(giotto.dbmatrix_compute = FALSE)} if saving dbMatrix
#' after each step of normalization workflow is not desired.
#' @keywords internal
#' @noRd
.compute_dbMatrix <- function(dbMatrix, name, verbose = TRUE) {
    # input validation
    if (!inherits(dbMatrix, "dbMatrix")) {
        stop("dbMatrix must be of class dbMatrix")
    }

    if (!is.character(name)) {
        stop("name must be a character")
    }

    # TODO: update with dbData generic
    con <- dbMatrix:::get_con(dbMatrix)

    # overwrite table by default
    if (name %in% DBI::dbListTables(con)) {
        DBI::dbRemoveTable(con, name)
    }

    if (verbose) {
        msg <- glue::glue("Computing {name} expression matrix on disk...")
        cat(msg)
    }

    dbMatrix[] |>
        dplyr::compute(temporary = FALSE, name = name)

    # TODO: update below with proper setters from dbMatrix
    dbMatrix[] <- dplyr::tbl(con, name) # reassign to computed mat
    dbMatrix@name <- name

    if (verbose) cat("done \n")

    return(dbMatrix)
}

#' @title RNA standard normalization
#' @name .rna_standard_normalization
#' @description standard function for RNA normalization
#' @returns giotto object
#' @keywords internal
#' @noRd
.rna_standard_normalization <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    library_size_norm = TRUE,
    scalefactor = 6e3,
    log_norm = TRUE,
    log_offset = 1,
    logbase = 2,
    scale_feats = TRUE,
    scale_cells = TRUE,
    scale_order = c("first_feats", "first_cells"),
    verbose = TRUE) {
    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: Standard normalization was developed for RNA data \n")
    }

    # evaluate provenance before modifying raw_expr in case h5_file exists
    if (isS4(raw_expr)) {
        provenance <- raw_expr@provenance
    } else {
        provenance <- NULL
    }


    feat_names <- rownames(raw_expr[])
    col_names <- colnames(raw_expr[])

    ## 1. library size normalize
    if (isTRUE(library_size_norm)) {
        norm_expr <- .lib_norm_giotto(
            mymatrix = raw_expr[],
            scalefactor = scalefactor
        )
    } else {
        norm_expr <- raw_expr[]
    }

    ## 2. log normalize
    if (isTRUE(log_norm)) {
        norm_expr <- .log_norm_giotto(
            mymatrix = norm_expr,
            base = logbase,
            offset = log_offset
        )
    }

    ## 3. scale
    if (isTRUE(scale_feats) && isTRUE(scale_cells)) {
        scale_order <- match.arg(
            arg = scale_order, choices = c("first_feats", "first_cells")
        )

        if (scale_order == "first_feats") {
            if (isTRUE(verbose)) {
                vmsg(.v = verbose, "first scale feats and then cells")
            }

            norm_scaled_expr <- t_flex(standardise_flex(
                x = t_flex(norm_expr), center = TRUE, scale = TRUE
            ))
            norm_scaled_expr <- standardise_flex(
                x = norm_scaled_expr, center = TRUE, scale = TRUE
            )
        } else if (scale_order == "first_cells") {
            if (isTRUE(verbose)) {
                vmsg(.v = verbose, "first scale cells and then feats")
            }

            norm_scaled_expr <- standardise_flex(
                x = norm_expr, center = TRUE, scale = TRUE
            )
            norm_scaled_expr <- t_flex(standardise_flex(
                x = t_flex(norm_scaled_expr), center = TRUE, scale = TRUE
            ))
        } else {
            stop("\n scale order must be given \n")
        }
    } else if (isTRUE(scale_feats)) {
        norm_scaled_expr <- t_flex(standardise_flex(
            x = t_flex(norm_expr), center = TRUE, scale = TRUE
        ))
    } else if (isTRUE(scale_cells)) {
        norm_scaled_expr <- standardise_flex(
            x = norm_expr, center = TRUE, scale = TRUE
        )
    } else {
        norm_scaled_expr <- NULL
    }


    ## 4. add cell and gene names back
    if (!is.null(norm_expr)) {
        rownames(norm_expr) <- feat_names
        colnames(norm_expr) <- col_names
    }
    if (!is.null(norm_scaled_expr)) {
        rownames(norm_scaled_expr) <- feat_names
        colnames(norm_scaled_expr) <- col_names
    }

    ## 5. create and set exprObj
    # Save dbMatrix to db
    compute_mat <- getOption("giotto.dbmatrix_compute", default = FALSE)
    if (compute_mat && !is.null(norm_expr)) {
        norm_expr <- .compute_dbMatrix(
            dbMatrix = norm_expr,
            name = "normalized",
            verbose = verbose
        )
    }

    norm_expr <- create_expr_obj(
        name = "normalized",
        exprMat = norm_expr,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        misc = NULL
    )

    # Save dbMatrix to db
    if (compute_mat && !is.null(norm_scaled_expr)) {
        norm_scaled_expr <- .compute_dbMatrix(
            dbMatrix = norm_scaled_expr,
            name = "scaled",
            verbose = verbose
        )
    }

    norm_scaled_expr <- create_expr_obj(
        name = "scaled",
        exprMat = norm_scaled_expr,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        misc = NULL
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(gobject, norm_expr, initialize = FALSE)
    gobject <- setGiotto(gobject, norm_scaled_expr, initialize = FALSE)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    ## 6. return Giotto object
    return(initialize(gobject))
}



#' @title RNA osmfish normalization
#' @name .rna_osmfish_normalization
#' @description function for RNA normalization according to osmFISH paper
#' @returns giotto object
#' @keywords internal
#' @noRd
.rna_osmfish_normalization <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    name = "custom",
    verbose = TRUE) {
    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: osmFISH normalization was developed for RNA in situ
                data \n")
    }

    # 1. normalize per gene with scale-factor equal to number of genes
    norm_feats <- (raw_expr[] / rowSums_flex(raw_expr[])) * nrow(raw_expr[])
    # 2. normalize per cells with scale-factor equal to number of cells
    norm_feats_cells <- t_flex((t_flex(norm_feats) /
        colSums_flex(norm_feats)) * ncol(raw_expr[]))

    # return results to Giotto object
    if (verbose == TRUE) {
        message(
            "\n osmFISH-like normalized data will be returned to the",
            name, "Giotto slot \n"
        )
    }

    norm_feats_cells <- create_expr_obj(
        name = name,
        exprMat = norm_feats_cells,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = raw_expr@provenance
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(giotto, norm_feats_cells)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
}


#' @title RNA pearson residuals normalization
#' @name rna_pears_resid_normalization
#' @description function for RNA normalization according to Lause/Kobak et al
#' paper
#' Adapted from https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c#file-pearson_residuals-r
#' @returns giotto object
#' @keywords internal
#' @noRd
.rna_pears_resid_normalization <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    theta = 100,
    name = "scaled",
    verbose = TRUE) {
    # print message with information #
    if (verbose) {
        message("using 'Lause/Kobak' method to normalize count matrix If used in
      published research, please cite:
      Jan Lause, Philipp Berens, Dmitry Kobak (2020).
      'Analytic Pearson residuals for normalization of single-cell RNA-seq UMI
      data' ")
    }

    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: pearson residual normalization was developed for RNA
                count normalization \n")
    }

    if (methods::is(raw_expr[], "HDF5Matrix")) {
        .csums <- .csum_nodrop.HDF5Matrix
        .rsums <- .rsum_nodrop.HDF5Matrix
    } else {
        .csums <- .csum_nodrop.Matrix
        .rsums <- .rsum_nodrop.Matrix
    }

    z <- .prnorm(x = raw_expr[], theta, .csums = .csums, .rsums = .rsums)
    z <- create_expr_obj(
        name = name,
        exprMat = z,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = prov(raw_expr)
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(gobject, z, verbose = verbose)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
}

.quantile_norm <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    name = "quantile",
    verbose = TRUE) {
    z <- .qnorm(x = raw_expr[])
    z <- create_expr_obj(
        name = name,
        exprMat = z,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = prov(raw_expr)
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(gobject, z, verbose = verbose)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
}

# pearson residuals normalization
# x      : raw expression matrix
# .csums : function for colSums that does not drop to vector
# .rsums : function for rowSums that does not drop to vector
.prnorm <- function(x,
    theta = 100,
    .csums = .csum_nodrop.Matrix,
    .rsums = .rsum_nodrop.Matrix) {
    # find 1. colsums, 2. rowsums, 3. matrix sum
    counts_sum0 <- .csums(x)
    counts_sum1 <- .rsums(x)
    counts_sum <- sum(x)

    # get residuals
    mu <- (counts_sum1 %*% counts_sum0) / counts_sum
    z <- (x - mu) / sqrt(mu + mu^2 / theta)

    # clip to be within the range [-sqrt(n), sqrt(n)]
    # This is done to prevent extreme values from dominating the analysis.
    n <- ncol(x)
    z[z > sqrt(n)] <- sqrt(n)
    z[z < -sqrt(n)] <- -sqrt(n)
    return(z)
}



# quantile normalization
.qnorm <- function(x) {
    # apply on features by default
    x <- t_flex(x)
    # Rank the values within each column
    ranked_data <- t_flex(MatrixGenerics::colRanks(x, ties.method = "average"))

    # Calculate the mean of sorted values across all columns
    rank_means <- rowMeans(apply(x, 2, sort))

    # Replace the original values with the rank means
    # TODO revisit for large matrices
    normalized_data <- apply(ranked_data, 2, function(idx) {
        .qnorm_vector(idx, rank_means)
    }) |>
        methods::as("Matrix")

    # Retain the original column names
    colnames(normalized_data) <- colnames(x)
    normalized_data <- t_flex(normalized_data)
    return(normalized_data)
}

# create lookup value vector for quantile norm.
# .5 indices should pull the mean of the adjacent values
# indices: index values with some values being .5, designating ranking ties
# values: values to pull from with the indices
.qnorm_vector <- function(indices, values) {
    sorted_values <- sort(values)
    lower_indices <- floor(indices)
    upper_indices <- ceiling(indices)
    lower_values <- sorted_values[lower_indices]
    upper_values <- sorted_values[upper_indices]
    weights <- indices - lower_indices
    result <- (1 - weights) * lower_values + weights * upper_values
    return(result)
}

.csum_nodrop.Matrix <- function(x) {
    x |>
        Matrix::colSums() |>
        matrix(nrow = 1L) |>
        methods::as("Matrix")
}
.rsum_nodrop.Matrix <- function(x) {
    x |>
        Matrix::rowSums() |>
        matrix(ncol = 1L) |>
        methods::as("Matrix")
}
.csum_nodrop.HDF5Matrix <- function(x) {
    x |>
        MatrixGenerics::colSums2() |>
        matrix(nrow = 1L) |>
        methods::as("HDF5Matrix")
}
.rsum_nodrop.HDF5Matrix <- function(x) {
    x |>
        MatrixGenerics::rowSums2() |>
        matrix(ncol = 1L) |>
        methods::as("HDF5Matrix")
}

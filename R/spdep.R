#' Compute spatial auto correlation using spdep
#'
#' @param gobject Input a Giotto object.
#' @param method Specify a method name to compute auto correlation.
#' Available methods include
#' \code{"geary.test", "lee.test", "lm.morantest","moran.test"}.
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use, default = normalized
#' @param spatial_network_to_use spatial network to use,
#' default = spatial_network
#' @param verbose be verbose
#' @param return_gobject if FALSE, results are returned as data.table.
#' If TRUE, values will be appended to feature metadata
#' @returns A data table with computed values for each feature.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' spdepAutoCorr(g)
#' @export
spdepAutoCorr <- function(
        gobject,
        method = c("geary.test", "lee.test", "lm.morantest", "moran.test"),
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = "normalized",
        spatial_network_to_use = "spatial_network",
        return_gobject = FALSE,
        verbose = FALSE) {
    # Check and match the specified method argument
    method <- match.arg(method)

    # Check gobject and set spat_unit and feat_type
    if (!is.null(gobject)) {
        spat_unit <- set_default_spat_unit(
            gobject = gobject,
            spat_unit = spat_unit
        )
        feat_type <- set_default_feat_type(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        )
    } else {
        stop("gobject has not been provided")
    }

    # Evaluate spatial autocorrelation using Giotto
    resultSpdepCor <- .evaluate_autocor_input(
        gobject = gobject,
        use_ext_vals = FALSE,
        use_sn = TRUE,
        use_expr = TRUE,
        use_meta = FALSE,
        spat_unit = spat_unit,
        feat_type = feat_type,
        feats = NULL,
        data_to_use = "expression",
        expression_values = expression_values,
        meta_cols = NULL,
        spatial_network_to_use = spatial_network_to_use,
        wm_method = "distance",
        wm_name = "spat_weights",
        node_values = NULL,
        weight_matrix = NULL,
        verbose = verbose
    )


    # Extract feats and weight_matrix from the result
    feat <- resultSpdepCor$feats
    weight_matrix <- resultSpdepCor$weight_matrix
    use_values <- resultSpdepCor$use_values

    # progressr
    nfeats <- length(feat)
    step_size <- ceiling(nfeats / 10L)

    result_list <- list()
    with_pbar({
        if (step_size > 1) {
            pb <- pbar(
                steps = nfeats / step_size
            )
        }
        result_list <- lapply_flex(
            seq_along(feat),
            future.packages = c("data.table", "spdep"),
            function(feat_value) {
                callSpdepVar <- callSpdep(
                    method = method,
                    x = use_values[, feat_value],
                    listw = spdep::mat2listw(weight_matrix, style = "W")
                )
                # Extract the estimated value from the result
                result_value <- callSpdepVar$estimate[1]
                temp_dt <- data.table(
                    feat_ID = feat[feat_value], value = result_value
                )
                # increment progress
                if (exists("pb")) if (feat_value %% step_size == 0) pb()
                return(temp_dt)
            }
        )
    })
    # combine results
    result_dt <- data.table::rbindlist(result_list)

    # Return the resulting datatable

    if (isTRUE(return_gobject)) {
        if (isTRUE(verbose)) {
            wrap_msg(
                "Appending", method,
                "results to feature metadata: fDataDT()"
            )
        }
        gobject <- addFeatMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = result_dt,
            by_column = TRUE,
            column_feat_ID = "feat_ID"
        )

        return(gobject)
    } else {
        return(result_dt)
    }
}


#' Call the spdep function with required parameters
#'
#' @param method Specify method name to call from spdep with its required
#' parameters.
#' @param ... Additional parameters for the function. See spdep documentation
#' for relevant parameters.
#' @returns Computed statistics from the specified method.
#' @export
#' @seealso \url{https://cran.r-project.org/web/packages/spdep/index.html}
callSpdep <- function(method, ...) {
    # Load the 'spdep' package if not already installed
    package_check(pkg_name = "spdep", repository = "CRAN", optional = FALSE)

    # Check if 'method' argument is NULL, if so, stop with an error
    if (is.null(method)) {
        stop("The 'method' argument has not been provided. Please specify a
            valid method.")
    }

    # Check if 'method' exists in the 'spdep' package, if not, stop with an
    # error
    method <- try(eval(get(method, envir = loadNamespace("spdep"))),
        silent = TRUE
    )
    if (inherits(method, "try-error")) {
        stop(paste(
            "Invalid method name. Method", method,
            "is not available in the spdep package."
        ))
    }

    # Fetch the arguments of the 'method' from 'spdep'
    allArgs <- args(method) %>%
        as.list() %>%
        names()

    # Capture arguments provided by the user
    methodparam <- list(...)

    # Check if the user provided the listw argument
    if ("listw" %in% names(methodparam)) {
        listw_arg <- methodparam$listw

        # Check if listw_arg is a matrix
        if (is.matrix(listw_arg)) {
            # Convert the matrix to a listw object
            listw_arg <- spdep::mat2listw(listw_arg, style = "W")
        } else if (!inherits(listw_arg, "listw")) {
            stop("listw must be either a matrix or a listw object.")
        }

        # Update the listw argument in methodparam
        methodparam$listw <- listw_arg
    }


    # Check if all user-provided arguments are valid
    if (all(!(names(methodparam)) %in% allArgs)) {
        stop("Invalid or missing parameters.")
    }
    # A vector of specified arguments that trigger
    # 'spW <- spweights.constants()'
    requiredArgs <- c("n", "n1", "n2", "n3", "nn", "S0", "S1", "S2")

    # Check if any of the specified arguments are required by the method
    if (any(requiredArgs %in% allArgs)) {
        # Obtain arguments from 'spweights.constants'
        spW <- spdep::spweights.constants(listw = methodparam$listw)
        # Combine user-provided arguments and 'spW', checking only against
        # 'feats' value
        combinedParams <- append(methodparam, spW)
    } else {
        combinedParams <- methodparam
    }


    # Identify common parameters between user and 'spdep'
    commonParams <- intersect(names(combinedParams), allArgs)

    # Create a named list of common parameters
    combinedParams <- combinedParams[commonParams]

    # Call the function with its parameters
    do.call(method, combinedParams)
}

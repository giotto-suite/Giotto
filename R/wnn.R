#' Multi omics integration with WNN
#'
#' @param gobject A Giotto object with individual PCA feat_types pre-calculated
#' @param spat_unit spatial unit
#' @param feat_types feature types to integrate. Default = c("rna", "protein")
#' @param reduction_methods reduction methods for each feature type. Default = c("pca", "pca")
#' @param reduction_names names of the reduction methods to use. Default = c("rna.pca", "protein.pca")
#' @param k number of neighbors to calculate cell distances, default = 20
#' @param integrated_feat_type integrated feature type (e.g. 'rna_protein')
#' @param matrix_result_name Default = 'theta_weighted_matrix'
#' @param verbose be verbose
#' @param w_names optional. Names for the weighted matrices. If NULL, automatic names composed by w_feat_type will be created.
#'
#' @returns A Giotto object with a new multiomics slot containing the theta_weighted_matrix and individual weight matrices.
#' @export
runWNN <- function(
        gobject,
        spat_unit = "cell",
        feat_types = c("rna", "protein"),
        reduction_methods = c("pca", "pca"),
        reduction_names = c("rna.pca", "protein.pca"),
        k = 20,
        integrated_feat_type = NULL,
        matrix_result_name = NULL,
        w_names = c(NULL, NULL),
        verbose = FALSE) {
    # validate Giotto object
    if (!inherits(gobject, "giotto")) {
        stop("gobject needs to be a giotto object")
    }

    # validate feat_types
    for (feat_type in feat_types) {
        if (!feat_type %in% names(
            slot(gobject, "dimension_reduction")$cells[[spat_unit]]
        )) {
            stop(paste(feat_type, " and their dimension reductions must exist in the Giotto object"))
        }
    }

    # extract kNN and PCA
    kNN_list <- list()
    pca_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        kNN_list[[feat_type]] <- getNearestNetwork(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            nn_type = "kNN",
            output = "igraph"
        )

        pca_list[[feat_type]] <- getDimReduction(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = reduction_methods[i],
            name = reduction_names[i],
            output = "matrix"
        )
    }


    ## get cell names
    cell_names <- GiottoClass:::get_cell_id(gobject,
        spat_unit = spat_unit
    )

    ######################## distances calculation ############################

    cell_distances <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        cell_distances[[feat_type]] <- list()

        message("Calculating ", feat_type, "-", feat_type, " distance")

        cell_distances[[feat_type]][[feat_type]] <- list()

        for (cell_a in cell_names) {
            my_kNN <- kNN_list[[feat_type]][[cell_a]][[cell_a]]

            cell_distances[[feat_type]][[feat_type]][[cell_a]] <- rep(0, k)
            names(cell_distances[[feat_type]][[feat_type]][[cell_a]]) <- names(my_kNN)

            for (cell_i in names(my_kNN)) {
                dimensions_cell_a_i <- pca_list[[feat_type]][c(cell_a, cell_i), ]
                cell_distances[[feat_type]][[feat_type]][[cell_a]][cell_i] <- sqrt(sum((
                    dimensions_cell_a_i[1, ] - dimensions_cell_a_i[2, ])^2))
            }
        }
    }

    ########################### all cell-cell distances ########################

    all_cell_distances <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        all_cell_distances[[feat_type]] <- list()

        if (verbose) {
            message(paste(
                "Calculating low dimensional cell-cell distances for",
                feat_type
            ))
        }

        pca_distances <- dist(pca_list[[feat_type]])
        all_cell_distances[[feat_type]][[feat_type]] <- as.matrix(pca_distances)
    }

    ######################## within-modality prediction ########################

    if (verbose) message("Calculating within-modality prediction")

    ### predicted feat_type1 feat_type1
    predicted_values <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        predicted_values[[feat_type]] <- list()
        predicted_values[[feat_type]][[feat_type]] <- list()

        for (cell_a in cell_names) {
            dimensions_cell_a <- pca_list[[feat_type]][kNN_list[[feat_type]][[cell_a]][[cell_a]], ]

            predicted_values[[feat_type]][[feat_type]][[cell_a]] <- colSums(dimensions_cell_a) / k
        }
    }


    ######################## cross-modality prediction #########################

    if (verbose) message("Calculating cross-modality prediction")

    ## predicted feat_type1 modality2
    feat_type1 <- feat_types[1]
    feat_type2 <- feat_types[2]

    predicted_values[[feat_type1]][[feat_type2]] <- list()

    for (cell_a in cell_names) {
        dimensions_cell_a <- pca_list[[feat_type1]][kNN_list[[feat_type2]][[cell_a]][[cell_a]], ]

        predicted_values[[feat_type1]][[feat_type2]][[cell_a]] <- colSums(dimensions_cell_a) / k
    }

    ## predicted modality2 feat_type1
    predicted_values[[feat_type2]][[feat_type1]] <- list()

    for (cell_a in cell_names) {
        dimensions_cell_a <- pca_list[[feat_type2]][kNN_list[[feat_type1]][[cell_a]][[cell_a]], ]

        predicted_values[[feat_type2]][[feat_type1]][[cell_a]] <- colSums(dimensions_cell_a) / k
    }

    ###################### calculate jaccard similarities ######################

    if (verbose) message("Calculating Jaccard similarities")

    sNN_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        sNN_result <- createNearestNetwork(gobject,
            spat_unit = "cell",
            feat_type = feat_type,
            type = "sNN",
            dim_reduction_to_use = reduction_methods[i],
            dim_reduction_name = reduction_names[i],
            dimensions_to_use = 1:100,
            return_gobject = FALSE,
            minimum_shared = 1,
            k = 20
        )

        sNN_list[[feat_type]] <- igraph::as_data_frame(sNN_result)
    }

    if (verbose) message("Calculating kernel bandwidths")

    # cell-specific kernel bandwidth.
    sigma_i <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        sigma_i[[feat_type]] <- numeric()


        for (cell_a in cell_names) {
            ### 20 small jaccard values
            jaccard_values <- sNN_list[[feat_type]][sNN_list[[feat_type]]$from == cell_a, ]

            if (nrow(jaccard_values == 20)) {
                further_cell_cell_distances <- all_cell_distances[[feat_type]][[feat_type]][
                    cell_a, jaccard_values$to
                ]
            } else {
                further_cell_cell_distances <- tail(sort(all_cell_distances[[feat_type]][[feat_type]][
                    cell_a,
                ]), 20)
            }

            sigma_i[[feat_type]][cell_a] <- mean(further_cell_cell_distances)
            # cell-specific kernel bandwidth.
        }
    }

    ###################### cell-specific modality weights ######################

    if (verbose) message("Calculating modality weights")

    ## feat_type1 feat_type1
    theta_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        theta_list[[feat_type]] <- list()
        theta_list[[feat_type]][[feat_type]] <- list()

        for (cell_a in cell_names) {
            feat_type1_i <- pca_list[[feat_type]][cell_a, ] # profile of current cell

            d_feat_type1_i_modality2_predicted <- sqrt(sum((
                feat_type1_i - predicted_values[[feat_type]][[feat_type]][[cell_a]])^2))

            first_knn <- names(sort(cell_distances[[feat_type]][[feat_type]][[cell_a]]))[1]

            feat_type1_knn1 <- pca_list[[feat_type]][first_knn, ] # profile of the nearest neighbor

            d_feat_type1_i_feat_type1_knn1 <- sqrt(sum((
                feat_type1_i - feat_type1_knn1)^2))

            difference_distances <- d_feat_type1_i_modality2_predicted -
                d_feat_type1_i_feat_type1_knn1
            max_value <- max(c(difference_distances, 0))

            theta_list[[feat_type]][[feat_type]][[cell_a]] <- exp((
                -max_value) / (sigma_i[[feat_type]][cell_a] -
                d_feat_type1_i_feat_type1_knn1))
        }
    }


    ## feat_type1 modality2
    feat_type1 <- feat_types[1]
    feat_type2 <- feat_types[2]

    theta_list[[feat_type1]][[feat_type2]] <- list()

    for (cell_a in cell_names) {
        feat_type1_i <- pca_list[[feat_type1]][cell_a, ] # profile of current cell
        d_feat_type1_i_modality2_predicted <- sqrt(sum((
            feat_type1_i - predicted_values[[feat_type1]][[feat_type2]][[cell_a]])^2))
        first_knn <- names(sort(cell_distances[[feat_type1]][[feat_type1]][[cell_a]]))[1]

        feat_type1_knn1 <- pca_list[[feat_type1]][first_knn, ] # profile of the nearest neighbor
        d_feat_type1_i_feat_type1_knn1 <- sqrt(sum((
            feat_type1_i - feat_type1_knn1)^2))

        difference_distances <- d_feat_type1_i_modality2_predicted -
            d_feat_type1_i_feat_type1_knn1
        max_value <- max(c(difference_distances, 0))

        theta_list[[feat_type1]][[feat_type2]][[cell_a]] <- exp((
            -max_value) / (sigma_i[[feat_type1]][cell_a] -
            d_feat_type1_i_feat_type1_knn1))
    }


    ## modality2 feat_type1
    theta_list[[feat_type2]][[feat_type1]] <- list()

    for (cell_a in cell_names) {
        modality2_i <- pca_list[[feat_type2]][cell_a, ] # profile of current cell
        d_modality2_i_feat_type1_predicted <- sqrt(sum((
            modality2_i - predicted_values[[feat_type2]][[feat_type1]][[cell_a]])^2))
        first_knn <- names(sort(cell_distances[[feat_type2]][[feat_type2]][[cell_a]]))[1]

        modality2_knn1 <- pca_list[[feat_type2]][first_knn, ] # profile of the nearest neighbor
        d_modality2_i_modality2_knn1 <- sqrt(sum((
            modality2_i - modality2_knn1)^2))

        difference_distances <- d_modality2_i_feat_type1_predicted -
            d_modality2_i_modality2_knn1
        max_value <- max(c(difference_distances, 0))

        theta_list[[feat_type2]][[feat_type1]][[cell_a]] <- exp((
            -max_value) / (sigma_i[[feat_type2]][cell_a] -
            d_modality2_i_modality2_knn1))
    }


    ##################### ratio of affinities ##################################

    if (verbose) message("Calculating WNN")

    epsilon <- 10^-4

    ratio_list <- list()

    ## feat_type1
    ratio_list[[feat_type1]] <- list()

    for (cell_a in cell_names) {
        ratio_list[[feat_type1]][[cell_a]] <- theta_list[[feat_type1]][[feat_type1]][[cell_a]] /
            (theta_list[[feat_type1]][[feat_type2]][[cell_a]] + epsilon)
    }

    ## modality2
    ratio_list[[feat_type2]] <- list()

    for (cell_a in cell_names) {
        ratio_list[[feat_type2]][[cell_a]] <- theta_list[[feat_type2]][[feat_type2]][[cell_a]] /
            (theta_list[[feat_type2]][[feat_type1]][[cell_a]] + epsilon)
    }


    ########################### normalization ##################################

    if (verbose) message("Calculating WNN normalization")

    w_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        w_list[[feat_type]] <- rep(0, length(cell_names))
        names(w_list[[feat_type]]) <- cell_names

        for (cell_a in cell_names) {
            w_list[[feat_type]][cell_a] <- exp(ratio_list[[feat_type]][[cell_a]]) /
                (exp(ratio_list[[feat_type1]][[cell_a]]) + exp(ratio_list[[feat_type2]][[cell_a]]))
        }
    }

    ######################### Calculating a WNN graph ##########################

    if (verbose) message("Calculating WNN graph")

    theta_weighted <- matrix(rep(0, length(cell_names) * length(cell_names)),
        ncol = length(cell_names),
        nrow = length(cell_names)
    )

    colnames(theta_weighted) <- cell_names
    rownames(theta_weighted) <- cell_names

    kernelpower <- 1

    ## theta_feat_type1

    theta_cella_cellb <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        theta_cella_cellb[[feat_type]] <- exp(-1 * (all_cell_distances[[feat_type]][[feat_type]] /
            sigma_i[[feat_type]])**kernelpower)
    }

    ## theta_weighted
    theta_weighted <- w_list[[feat_type1]] * theta_cella_cellb[[feat_type1]] +
        w_list[[feat_type2]] * theta_cella_cellb[[feat_type2]]


    # save theta_weighted

    if (verbose) message("Saving WNN results")


    ## set integrated feat_type and result name if not provided
    if (is.null(integrated_feat_type)) {
        integrated_feat_type <- paste0(feat_type1, "_", feat_type2)
    }

    if (is.null(matrix_result_name)) {
        matrix_result_name <- "theta_weighted_matrix"
    }

    gobject <- setMultiomics(
        gobject = gobject,
        result = theta_weighted,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = "WNN",
        result_name = matrix_result_name,
        verbose = TRUE
    )

    # save feat_types weight
    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        if (is.null(w_names[i])) {
            w_names[i] <- paste0("w_", feat_type)
        }

        gobject <- setMultiomics(
            gobject = gobject,
            result = w_list[feat_type],
            spat_unit = spat_unit,
            feat_type = integrated_feat_type,
            integration_method = "WNN",
            result_name = w_names[i],
            verbose = TRUE
        )
    }

    return(gobject)
}


#' Run integrated UMAP
#'
#' @param gobject A giotto object
#' @param spat_unit spatial unit
#' @param feat_types feature types to integrate. Default = c("rna", "protein")
#' @param k k number
#' @param spread UMAP param: spread
#' @param min_dist UMAP param: min_dist
#' @param integrated_feat_type integrated feature type (e.g. 'rna_protein')
#' @param integration_method multiomics integration method used. Default = 'WNN'
#' @param matrix_result_name Default = 'theta_weighted_matrix'
#' @param force force calculation of integrated kNN. Default = FALSE
#' @param ... additional UMAP parameters
#'
#' @returns A Giotto object with integrated UMAP
#' @export
runIntegratedUMAP <- function(
        gobject,
        spat_unit = "cell",
        feat_types = c("rna", "protein"),
        integrated_feat_type = NULL,
        integration_method = "WNN",
        matrix_result_name = "theta_weighted_matrix",
        k = 20,
        spread = 5,
        min_dist = 0.01,
        force = FALSE,
        ...) {
    # validate feat_types
    for (feat_type in feat_types) {
        if (!feat_type %in% names(
            slot(gobject, "dimension_reduction")$cells[[spat_unit]]
        )) {
            stop(paste(feat_type, " and their dimension reductions must exist in the Giotto object"))
        }
    }

    feat_type1 <- feat_types[1]
    feat_type2 <- feat_types[2]

    if (is.null(integrated_feat_type)) {
        integrated_feat_type <- paste0(feat_type1, "_", feat_type2)
    }

    theta_weighted <- getMultiomics(gobject,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = integration_method,
        result_name = matrix_result_name
    )

    theta_weighted[is.na(theta_weighted)] <- 0

    if (is.null(gobject@nn_network[[spat_unit]][[
        feat_type1
    ]]$kNN$integrated_kNN) || force == TRUE) {
        ################# Calculate integrated Nearest Neighbors ###############

        message("Calculating integrated Nearest Neighbors")

        cell_names <- colnames(theta_weighted)

        nn_network <- dbscan::kNN(x = theta_weighted, k = k, sort = TRUE)
        from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <-
            shared <- NULL
        nn_network_dt <- data.table::data.table(
            from = rep(
                seq_len(nrow(nn_network$id)),
                k
            ),
            to = as.vector(nn_network$id),
            weight = 1 / (1 + as.vector(nn_network$dist)),
            distance = as.vector(nn_network$dist)
        )
        nn_network_dt[, `:=`(from_cell_ID, cell_names[from])]
        nn_network_dt[, `:=`(to_cell_ID, cell_names[to])]
        all_index <- unique(
            x = c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID)
        )

        ################################ Create igraph #########################

        message("Creating igraph")

        nn_network_igraph <- igraph::graph_from_data_frame(
            nn_network_dt[, .(from_cell_ID, to_cell_ID, weight, distance)],
            directed = TRUE,
            vertices = all_index
        )

        ## store igraph
        nnNetObj <- create_nn_net_obj(
            name = "integrated_kNN",
            nn_type = "kNN",
            igraph = nn_network_igraph,
            spat_unit = spat_unit,
            feat_type = feat_type1
        )

        gobject <- setGiotto(gobject, nnNetObj)

        ## store nn_network id
        gobject <- setMultiomics(
            gobject = gobject,
            result = nn_network$id,
            spat_unit = spat_unit,
            feat_type = integrated_feat_type,
            integration_method = "WNN",
            result_name = "integrated_kNN_id",
            verbose = TRUE
        )

        ## store nn_network dist
        gobject <- setMultiomics(
            gobject = gobject,
            result = nn_network$dist,
            spat_unit = spat_unit,
            feat_type = integrated_feat_type,
            integration_method = "WNN",
            result_name = "integrated_kNN_dist",
            verbose = TRUE
        )
    }

    ######################### Calculate integrated UMAP ########################

    message("Calculating integrated UMAP")

    nn_network_id <- getMultiomics(gobject,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = integration_method,
        result_name = "integrated_kNN_id"
    )

    nn_network_dist <- getMultiomics(gobject,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = integration_method,
        result_name = "integrated_kNN_dist"
    )


    #### using nn_network pre-calculation
    set.seed(4567)
    integrated_umap <- uwot::umap(
        X = theta_weighted,
        n_neighbors = k,
        nn_method = list(
            idx = nn_network_id,
            dist = nn_network_dist
        ),
        spread = spread,
        min_dist = min_dist,
        ...
    )

    colnames(integrated_umap) <- c("Dim.1", "Dim.2")

    ## add umap
    for (feat_type in feat_types) {
        umap_object <- createDimObj(
            coordinates = integrated_umap,
            spat_unit = spat_unit,
            feat_type = feat_type,
            method = "umap",
            name = "integrated.umap"
        )

        gobject <- setDimReduction(
            gobject = gobject,
            x = umap_object,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction_method = "umap",
            name = "integrated.umap"
        )
    }

    return(gobject)
}

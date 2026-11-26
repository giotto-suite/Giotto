## spatial gene detection ####

#' @title Fisher exact test
#' @name spat_fisher_exact
#' @description Perform fisher exact test
#' @returns A list with class "htest"
NULL

#' @rdname spat_fisher_exact
#' @keywords internal
.spat_fish_func <- function(
        feat,
        bin_matrix,
        spat_mat,
        calc_hub = FALSE,
        hub_min_int = 3) {
    feat_vector <- bin_matrix[rownames(bin_matrix) == feat, ]

    feat_vectorA <- feat_vector[names(feat_vector) %in% rownames(spat_mat)]
    feat_vectorA <- feat_vectorA[match(rownames(spat_mat), names(feat_vectorA))]

    feat_vectorB <- feat_vector[names(feat_vector) %in% colnames(spat_mat)]
    feat_vectorB <- feat_vectorB[match(colnames(spat_mat), names(feat_vectorB))]

    test1 <- spat_mat * feat_vectorA
    test2 <- t_flex(t_flex(spat_mat) * feat_vectorB)

    sourcevalues <- test1[spat_mat == 1]
    targetvalues <- test2[spat_mat == 1]

    # option 1
    test <- paste0(sourcevalues, "-", targetvalues)


    if (length(unique(test)) < 4) {
        possibs <- c("1-1", "0-1", "1-0", "0-0")
        missings_possibs <- possibs[!possibs %in% unique(test)]
        test <- c(test, missings_possibs)

        table_test <- table(test)
        table_test[names(table_test) %in% missings_possibs] <- 0
        table_matrix <- matrix(table_test, byrow = TRUE, nrow = 2)
    } else {
        table_matrix <- matrix(table(test), byrow = TRUE, nrow = 2)
    }

    if (calc_hub == TRUE) {
        high_cells <- names(feat_vector[feat_vector == 1])
        subset_spat_mat <- spat_mat[
            rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in%
                high_cells
        ]

        if (length(subset_spat_mat) == 1) {
            hub_nr <- 0
        } else {
            subset_spat_mat <- spat_mat[
                rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in%
                    high_cells
            ]
            rowhubs <- rowSums_flex(subset_spat_mat)
            colhubs <- colSums_flex(subset_spat_mat)
            hub_nr <- length(unique(c(
                names(colhubs[colhubs > hub_min_int]),
                names(rowhubs[colhubs > hub_min_int])
            )))
        }

        fish_res <- stats::fisher.test(table_matrix)[c("p.value", "estimate")]
        return(c(feats = list(feat), fish_res, hubs = list(hub_nr)))
    } else {
        fish_res <- stats::fisher.test(table_matrix)[c("p.value", "estimate")]
        return(c(feats = list(feat), fish_res))
    }
}

#' @describeIn spat_fisher_exact data.table implementation
#' @keywords internal
.spat_fish_func_dt <- function(
        bin_matrix_DTm,
        spat_netw_min,
        calc_hub = FALSE,
        hub_min_int = 3,
        cores = NA) {
    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    # data.table variables
    from_value <- to_value <- feat_ID <- N <- to <- from <- cell_ID <-
        V1 <- NULL

    # get binarized expression values for the neighbors
    spatial_network_min_ext <- data.table::merge.data.table(
        spat_netw_min,
        bin_matrix_DTm,
        by.x = "from",
        by.y = "variable",
        allow.cartesian = TRUE
    )
    data.table::setnames(spatial_network_min_ext, "value", "from_value")

    spatial_network_min_ext <- data.table::merge.data.table(
        spatial_network_min_ext,
        by.x = c("to", "feat_ID"),
        bin_matrix_DTm, by.y = c("variable", "feat_ID")
    )
    data.table::setnames(spatial_network_min_ext, "value", "to_value")


    # summarize the different combinations
    spatial_network_min_ext[, combn := paste0(from_value, "-", to_value)]
    freq_summary <- spatial_network_min_ext[, .N, by = .(feat_ID, combn)]
    data.table::setorder(freq_summary, feat_ID, combn)

    feats <- unique(freq_summary$feat_ID)
    all_combn <- c("0-0", "0-1", "1-0", "1-1")

    # create a zeroes DT to add missing observations
    freq_summary_zeroes <- data.table::data.table(
        feat_ID = rep(feats, each = 4),
        combn = rep(all_combn, length(feats)),
        N = 0
    )
    freq_summary2 <- rbind(freq_summary, freq_summary_zeroes)
    freq_summary2[, N := sum(N), by = .(feat_ID, combn)]
    freq_summary2 <- unique(freq_summary2)

    # sort the combinations and run fisher test
    data.table::setorder(freq_summary2, feat_ID, combn, -N)
    fish_results <- freq_summary2[, stats::fisher.test(
        matrix(N, nrow = 2)
    )[c(1, 3)], by = feat_ID]


    ## hubs ##
    if (calc_hub == TRUE) {
        double_pos <- spatial_network_min_ext[combn == "1-1"]

        double_pos_to <- double_pos[, .N, by = .(feat_ID, to)]
        data.table::setnames(double_pos_to, "to", "cell_ID")
        double_pos_from <- double_pos[, .N, by = .(feat_ID, from)]
        data.table::setnames(double_pos_from, "from", "cell_ID")

        double_pos_both <- rbind(double_pos_to, double_pos_from)
        double_pos_both <- double_pos_both[, sum(N), by = .(feat_ID, cell_ID)]
        data.table::setorder(double_pos_both, feat_ID, -V1)

        # get hubs and add 0's
        hub_DT <- double_pos_both[V1 > hub_min_int, .N, by = feat_ID]
        hub_DT_zeroes <- data.table::data.table(feat_ID = unique(
            spatial_network_min_ext$feat_ID
        ), N = 0)
        hub_DT2 <- rbind(hub_DT, hub_DT_zeroes)

        hub_DT2 <- hub_DT2[, sum(N), by = feat_ID]
        data.table::setnames(hub_DT2, "V1", "hub_nr")

        fish_results <- data.table::merge.data.table(
            fish_results, hub_DT2,
            by = "feat_ID"
        )
    }

    return(fish_results)
}




#' @title Spatial odds-ratio
#' @name spat_odds_ratio
#' @description calculate odds-ratio
#' @returns numeric
NULL

#' @rdname spat_odds_ratio
#' @keywords internal
.spat_or_func <- function(
        feat,
        bin_matrix,
        spat_mat,
        calc_hub = FALSE,
        hub_min_int = 3) {
    feat_vector <- bin_matrix[rownames(bin_matrix) == feat, ]

    feat_vectorA <- feat_vector[names(feat_vector) %in% rownames(spat_mat)]
    feat_vectorA <- feat_vectorA[match(rownames(spat_mat), names(feat_vectorA))]

    feat_vectorB <- feat_vector[names(feat_vector) %in% colnames(spat_mat)]
    feat_vectorB <- feat_vectorB[match(colnames(spat_mat), names(feat_vectorB))]

    test1 <- spat_mat * feat_vectorA
    test2 <- t_flex(t_flex(spat_mat) * feat_vectorB)

    sourcevalues <- test1[spat_mat == 1]
    targetvalues <- test2[spat_mat == 1]

    # option 1
    test <- paste0(sourcevalues, "-", targetvalues)


    if (length(unique(test)) < 4) {
        possibs <- c("1-1", "0-1", "1-0", "0-0")
        missings_possibs <- possibs[!possibs %in% unique(test)]
        test <- c(test, missings_possibs)

        table_test <- table(test)
        table_test[names(table_test) %in% missings_possibs] <- 0
        table_matrix <- matrix(table_test, byrow = TRUE, nrow = 2)
    } else {
        table_matrix <- matrix(table(test), byrow = TRUE, nrow = 2)
    }


    if (calc_hub == TRUE) {
        high_cells <- names(feat_vector[feat_vector == 1])
        subset_spat_mat <- spat_mat[
            rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in%
                high_cells
        ]

        if (length(subset_spat_mat) == 1) {
            hub_nr <- 0
        } else {
            rowhubs <- rowSums_flex(subset_spat_mat)
            colhubs <- colSums_flex(subset_spat_mat)
            hub_nr <- length(unique(c(names(
                colhubs[colhubs > hub_min_int]
            ), names(
                rowhubs[colhubs > hub_min_int]
            ))))
        }

        fish_matrix <- table_matrix
        fish_matrix <- fish_matrix / 1000
        OR <- ((fish_matrix[1] * fish_matrix[4]) /
            (fish_matrix[2] * fish_matrix[3]))

        return(c(feats = list(feat), OR, hubs = list(hub_nr)))
    }

    fish_matrix <- table_matrix
    fish_matrix <- fish_matrix / 1000
    OR <- ((fish_matrix[1] * fish_matrix[4]) / (
        fish_matrix[2] * fish_matrix[3]))
    return(c(feats = list(feat), OR))
}


#' @describeIn spat_odds_ratio data.table implementation
#' @keywords internal
.spat_or_func_dt <- function(
        bin_matrix_DTm,
        spat_netw_min,
        calc_hub = FALSE,
        hub_min_int = 3,
        cores = NA) {
    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    # data.table variables
    from_value <- to_value <- feat_ID <- N <- to <- from <- cell_ID <-
        V1 <- NULL

    # get binarized expression values for the neighbors
    spatial_network_min_ext <- data.table::merge.data.table(
        spat_netw_min, bin_matrix_DTm,
        by.x = "from", by.y = "variable",
        allow.cartesian = TRUE
    )
    data.table::setnames(spatial_network_min_ext, "value", "from_value")

    spatial_network_min_ext <- data.table::merge.data.table(
        spatial_network_min_ext,
        by.x = c("to", "feat_ID"),
        bin_matrix_DTm, by.y = c("variable", "feat_ID")
    )
    data.table::setnames(spatial_network_min_ext, "value", "to_value")


    # summarize the different combinations
    spatial_network_min_ext[, combn := paste0(from_value, "-", to_value)]
    freq_summary <- spatial_network_min_ext[, .N, by = .(feat_ID, combn)]
    data.table::setorder(freq_summary, feat_ID, combn)

    feats <- unique(freq_summary$feat_ID)
    all_combn <- c("0-0", "0-1", "1-0", "1-1")

    # create a zeroes DT to add missing observations
    freq_summary_zeroes <- data.table::data.table(
        feat_ID = rep(feats, each = 4),
        combn = rep(all_combn, length(feats)),
        N = 0
    )
    freq_summary2 <- rbind(freq_summary, freq_summary_zeroes)
    freq_summary2[, N := sum(N), by = .(feat_ID, combn)]
    freq_summary2 <- unique(freq_summary2)

    # sort the combinations and run fisher test
    setorder(freq_summary2, feat_ID, combn, -N)
    or_results <- freq_summary2[
        , .or_test_func(matrix(N, nrow = 2)),
        by = feat_ID
    ]


    ## hubs ##
    if (calc_hub == TRUE) {
        double_pos <- spatial_network_min_ext[combn == "1-1"]

        double_pos_to <- double_pos[, .N, by = .(feat_ID, to)]
        data.table::setnames(double_pos_to, "to", "cell_ID")
        double_pos_from <- double_pos[, .N, by = .(feat_ID, from)]
        data.table::setnames(double_pos_from, "from", "cell_ID")

        double_pos_both <- rbind(double_pos_to, double_pos_from)
        double_pos_both <- double_pos_both[, sum(N), by = .(feat_ID, cell_ID)]
        data.table::setorder(double_pos_both, feat_ID, -V1)

        # get hubs and add 0's
        hub_DT <- double_pos_both[V1 > hub_min_int, .N, by = feat_ID]
        hub_DT_zeroes <- data.table::data.table(
            feat_ID = unique(spatial_network_min_ext$feat_ID), N = 0
        )
        hub_DT2 <- rbind(hub_DT, hub_DT_zeroes)

        hub_DT2 <- hub_DT2[, sum(N), by = feat_ID]
        data.table::setnames(hub_DT2, "V1", "hub_nr")

        or_results <- data.table::merge.data.table(
            or_results, hub_DT2,
            by = "feat_ID"
        )
    }

    return(or_results)
}


#' @title Odds ratio test
#' @name .or_test_func
#' @description calculate odds-ratio from a 2x2 matrix
#' @returns list
#' @keywords internal
.or_test_func <- function(matrix) {
    OR <- ((matrix[1] * matrix[4]) / (matrix[2] * matrix[3]))
    list("estimate" = OR)
}



#' @title Calculate spatial enrichment
#' @name calculate_spatial_enrichment
#' @description Calculate spatial enrichment. Multiple methods are provided.
#' @returns spatial enrichment
NULL


#' @describeIn calculate_spatial_enrichment calculate using a 'simple' and
#' efficient for loop
#' @keywords internal
.calc_spatial_enrichment_minimum <- function(
        spatial_network,
        bin_matrix,
        adjust_method = "fdr",
        do_fisher_test = TRUE) {
    # data.table variables
    from <- to <- feats <- variable <- value <- p.value <- adj.p.value <-
        score <- estimate <- NULL

    spatial_network_min <- spatial_network[, .(from, to)]

    all_colindex <- seq_len(ncol(bin_matrix))
    names(all_colindex) <- colnames(bin_matrix)

    # code for possible combinations
    convert_code <- c(1, 2, 3, 4)
    names(convert_code) <- c("0-0", "0-1", "1-0", "1-1")

    # preallocate final matrix for results
    matrix_res <- matrix(
        data = NA, nrow = nrow(bin_matrix), ncol = nrow(spatial_network_min)
    )

    ## 1. summarize results for each edge in the network
    for (row_i in seq_len(nrow(spatial_network_min))) {
        from_id <- spatial_network_min[row_i][["from"]]
        to_id <- spatial_network_min[row_i][["to"]]

        sumres <- data.table::as.data.table(bin_matrix[
            , all_colindex[c(from_id, to_id)]
        ])
        sumres[, combn := paste0(get(from_id), "-", get(to_id))]

        code_res <- convert_code[sumres$combn]
        matrix_res[, row_i] <- code_res
    }

    rownames(matrix_res) <- rownames(bin_matrix)


    # preallocate matrix for table results
    table_res <- matrix(data = NA, nrow(matrix_res), ncol = 4)

    ## 2. calculate the frequencies of possible combinations ##
    # '0-0' = 1, '0-1' = 2, '1-0' = 3 and '1-1' = 4
    for (row_i in seq_len(nrow(matrix_res))) {
        x <- matrix_res[row_i, ]
        x <- factor(x, levels = c(1, 2, 3, 4))
        tabres <- as.vector(table(x))

        table_res[row_i, ] <- tabres
    }

    rownames(table_res) <- rownames(matrix_res)
    colnames(table_res) <- seq_len(4)

    rable_resDT <- data.table::as.data.table(table_res)
    rable_resDT[, feats := rownames(table_res)]

    rable_resDTm <- data.table::melt.data.table(rable_resDT, id.vars = "feats")
    data.table::setorder(rable_resDTm, feats, variable)

    ## run fisher test ##
    if (do_fisher_test == TRUE) {
        results <- rable_resDTm[, stats::fisher.test(matrix(
            value,
            nrow = 2
        ))[c(1, 3)], by = feats]

        # replace zero p-values with lowest p-value
        min_pvalue <- min(results$p.value[results$p.value > 0])
        results[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
        results[, adj.p.value := stats::p.adjust(
            p.value,
            method = adjust_method
        )]

        # sort feats based on p-value and estimate
        results[, score := -log(p.value) * estimate]
        data.table::setorder(results, -score)
    } else {
        results <- rable_resDTm[, .or_test_func(matrix(
            value,
            nrow = 2
        )), by = feats]
        data.table::setorder(results, -estimate)
    }

    return(results)
}

#' @describeIn calculate_spatial_enrichment calculate using 'matrix'
#' implementation
#' @keywords internal
.calc_spatial_enrichment_matrix <- function(
        spatial_network,
        bin_matrix,
        adjust_method = "fdr",
        do_fisher_test = TRUE,
        do_parallel = TRUE,
        cores = NA,
        calc_hub = FALSE,
        hub_min_int = 3,
        verbose = TRUE) {
    # data.table variables
    feats <- p.value <- estimate <- adj.p.value <- score <- NULL

    # convert spatial network data.table to spatial matrix. An edge is 1
    # whatever its distance, so edges between coincident cells (distance 0)
    # are kept.
    dc_spat_network <- data.table::dcast.data.table(
        data.table::copy(spatial_network)[, ".edge" := 1],
        formula = to ~ from, value.var = ".edge", fill = 0
    )
    spat_mat <- dt_to_matrix(dc_spat_network)
    spat_mat[spat_mat > 0] <- 1


    ## parallel
    if (do_parallel == TRUE) {
        if (do_fisher_test == TRUE) {
            save_list <- suppressMessages(lapply_flex(
                X = rownames(bin_matrix), cores = cores, fun = .spat_fish_func,
                bin_matrix = bin_matrix, spat_mat = spat_mat,
                calc_hub = calc_hub, hub_min_int = hub_min_int
            ))
        } else {
            save_list <- suppressMessages(lapply_flex(
                X = rownames(bin_matrix), cores = cores, fun = .spat_or_func,
                bin_matrix = bin_matrix, spat_mat = spat_mat,
                calc_hub = calc_hub, hub_min_int = hub_min_int
            ))
        }
    } else {
        ## serial
        save_list <- list()

        if (do_fisher_test == TRUE) {
            for (feat in rownames(bin_matrix)) {
                if (verbose == TRUE) print(feat)

                save_list[[feat]] <- suppressMessages(.spat_fish_func(
                    feat = feat, bin_matrix = bin_matrix, spat_mat = spat_mat,
                    calc_hub = calc_hub, hub_min_int = hub_min_int
                ))
            }
        } else {
            for (feat in rownames(bin_matrix)) {
                if (verbose == TRUE) print(feat)

                save_list[[feat]] <- suppressMessages(.spat_or_func(
                    feat = feat, bin_matrix = bin_matrix, spat_mat = spat_mat,
                    calc_hub = calc_hub, hub_min_int = hub_min_int
                ))
            }
        }
    }

    result <- data.table::as.data.table(do.call("rbind", save_list))
    result[, feats := unlist(feats)]


    if (do_fisher_test == TRUE) {
        result[, c("p.value", "estimate") := list(
            as.numeric(p.value), as.numeric(estimate)
        )]

        # convert p.value = 0 to lowest p-value
        min_pvalue <- min(result$p.value[result$p.value > 0])
        result[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
        result[, adj.p.value := stats::p.adjust(
            p.value,
            method = adjust_method
        )]

        result[, score := -log(p.value) * estimate]
        data.table::setorder(result, -score)
    } else {
        data.table::setnames(result, "V1", "estimate")
        data.table::setorder(result, -estimate)
    }

    return(result)
}


#' @describeIn calculate_spatial_enrichment calculate using 'data.table'
#' implementation
#' @keywords internal
.calc_spatial_enrichment_dt <- function(bin_matrix,
    spatial_network,
    calc_hub = FALSE,
    hub_min_int = 3,
    group_size = "automatic",
    do_fisher_test = TRUE,
    adjust_method = "fdr",
    cores = NA) {
    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    # data.table variables
    from <- to <- feat_ID <- p.value <- adj.p.value <- score <-
        estimate <- NULL

    # create minimum spatial network
    spat_netw_min <- spatial_network[, .(from, to)]

    # divide matrix in groups
    if (!is.na(group_size) & is.numeric(group_size)) {
        group_size <- group_size
        if (group_size > nrow(bin_matrix)) {
            stop("group_size is too big, it can not be greater than the
                number of feats")
        }
    } else if (group_size == "automatic") {
        test_number <- ceiling(nrow(bin_matrix) / 10)
        test_number <- max(2, test_number)
        group_size <- min(200, test_number)
    }

    groups <- ceiling(nrow(bin_matrix) / group_size)
    cut_groups <- cut(seq_len(nrow(bin_matrix)),
        breaks = groups,
        labels = seq_len(groups)
    )
    if (any(table(cut_groups) == 1)) {
        stop(
            "With group size = ", group_size,
            " you have a single gene in a group. Manually pick another group
            size"
        )
    }
    indexes <- seq_len(nrow(bin_matrix))
    names(indexes) <- cut_groups


    total_list <- list()
    for (group in unique(cut_groups)) {
        sel_indices <- indexes[names(indexes) == group]

        bin_matrix_DT <- data.table::as.data.table(bin_matrix[sel_indices, ])
        bin_matrix_DT[, feat_ID := rownames(bin_matrix[sel_indices, ])]
        bin_matrix_DTm <- data.table::melt.data.table(
            bin_matrix_DT,
            id.vars = "feat_ID"
        )

        if (do_fisher_test == TRUE) {
            test <- .spat_fish_func_dt(
                bin_matrix_DTm = bin_matrix_DTm,
                spat_netw_min = spat_netw_min,
                calc_hub = calc_hub,
                hub_min_int = hub_min_int,
                cores = cores
            )
        } else {
            test <- .spat_or_func_dt(
                bin_matrix_DTm = bin_matrix_DTm,
                spat_netw_min = spat_netw_min,
                calc_hub = calc_hub,
                hub_min_int = hub_min_int,
                cores = cores
            )
        }


        total_list[[group]] <- test
    }

    result <- do.call("rbind", total_list)

    if (do_fisher_test == TRUE) {
        min_pvalue <- min(result$p.value[result$p.value > 0])
        result[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
        result[, adj.p.value := stats::p.adjust(
            p.value,
            method = adjust_method
        )]

        result[, score := -log(p.value) * estimate]
        data.table::setorder(result, -score)
        data.table::setnames(result, "feat_ID", "feats")
    } else {
        data.table::setorder(result, -estimate)
        data.table::setnames(result, "feat_ID", "feats")
    }

    return(result)
}






#' @title binSpect
#' @name binSpect
#' @description Previously: `binGetSpatialGenes()`. \cr
#' BinSpect (Binary Spatial Extraction of genes) is a fast computational method
#' that identifies genes with a spatially coherent expression pattern. \cr
#' There are several functions documented together here, mainly differing in
#' how to provide expression and spatial connectivity/networks information.
#' When data is in a `giotto` object, use `binSpect()` which wraps
#' `binSpectSingle()` and `binSpectMulti()`.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param bin_method method to binarize gene expression
#' @param expression_values expression values to use
#' @param subset_feats only select a subset of features to test
#' @param spatial_network_name name of spatial network to use
#' (default = 'spatial_network')
#' @param spatial_network_k different k's for a spatial kNN to evaluate
#' @param reduce_network default uses the full network
#' @param kmeans_algo kmeans algorithm to use
#' (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use
#' (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high
#' expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation `"fast"` (default), or a legacy implementation:
#' `"data.table"`, `"simple"` or `"matrix"`. `"fast"` runs the join-count
#' test of [spatialAutoCor()] (`method = "binspect"`): an exact 2-means or
#' rank binarization and the closed-form randomisation null instead of a
#' Fisher test on edges, on in-memory and GiottoDisk-backed objects. Its
#' result has `statistic`, `z`, `p.value`, `adj.p.value`, `n_on` and
#' `odds_ratio`, ordered by `z`. The kmeans tuning arguments do not apply to
#' it. `spatial_network_k`, `bin_matrix`, `calc_hub = TRUE`,
#' `do_fisher_test = FALSE` and a non-default `kmeans_algo` need a legacy
#' implementation, so with any of them binSpect uses `"data.table"` and says
#' so. The legacy implementations return `p.value`, `estimate`, `score`,
#' `av_expr` and `high_expr`.
#' @param group_size number of genes to process together with data.table
#' implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if `do_parallel = TRUE`
#' @param verbose be verbose
#' @param knn_params list of parameters to create spatial kNN network
#' @param set.seed deprecated. Use \code{seed} param instead
#' @param seed seed for kmeans binarization. When \code{NULL}, no seed is set.
#' Otherwise, accepts a numeric input that will be used as seed.
#' @param bin_matrix a binarized matrix, when provided it will skip the
#' binarization process
#' @param summarize summarize the p-values or adjusted p-values
#' @param return_gobject whether to return values attached to the gobject or
#' separately (default)
#' @returns data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene
#' expression binarization.
#' Both methods are identicial except for how binarization is performed.
#'   1. **binarize:** Each gene is binarized (0 or 1) in each cell with
#'   **kmeans** (k = 2) or based on **rank** percentile
#'   2. **network:** All cells are connected through a spatial network based
#'   on the physical coordinates
#'   3. **contingency table:** A contingency table is calculated based on all
#'   edges of neighboring cells and the binarized expression
#'   (0-0, 0-1, 1-0 or 1-1)
#'   4. For each gene an odds-ratio (OR) and fisher.test (optional) is
#'   calculated
#'
#' Three different kmeans algorithms have been implemented:
#'   1. **kmeans:** default, see \code{\link[stats]{kmeans}}
#'   2. **kmeans_arma:** from ClusterR, see \code{\link[ClusterR]{KMeans_arma}}
#'   3. **kmeans_arma_subst:** from ClusterR, see
#'   \code{\link[ClusterR]{KMeans_arma}}, but randomly subsets the vector
#'   for each gene to increase speed. Change extreme_nr and sample_nr for
#'   control.
#'
#' Other statistics are provided (optional):
#'   * Number of cells with high expression (binary = 1)
#'   * Average expression of each gene within high expressing cells
#'   * Number of hub cells, these are high expressing cells that have a user
#'   defined number of high expressing neighbors
#'
#' By selecting a subset of likely spatial genes
#' (e.g. soft thresholding highly variable genes) can accelerate the speed.
#' The simple implementation is usually faster, but lacks the possibility to
#' run in parallel and to calculate hub cells.\cr
#' The data.table implementation might be more appropriate for large datasets
#' by setting the `group_size` (number of genes) parameter to divide the
#' workload.
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' binSpect(g)
#'
#' binSpectSingle(g)
#'
#' g_expression <- getExpression(g, output = "matrix")
#' g_spat_net <- getSpatialNetwork(g, output = "networkDT")
#'
#' binSpectSingleMatrix(
#'     expression_matrix = g_expression,
#'     spatial_network = g_spat_net
#' )
#'
NULL




#' @rdname binSpect
#' @export
binSpect <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    bin_method = c("kmeans", "rank"),
    expression_values = c("normalized", "scaled", "custom"),
    subset_feats = NULL,
    spatial_network_name = "Delaunay_network",
    spatial_network_k = NULL,
    reduce_network = FALSE,
    kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
    nstart = 3,
    iter_max = 10,
    extreme_nr = 50,
    sample_nr = 50,
    percentage_rank = 30,
    do_fisher_test = TRUE,
    adjust_method = "fdr",
    calc_hub = FALSE,
    hub_min_int = 3,
    get_av_expr = TRUE,
    get_high_expr = TRUE,
    implementation = c("fast", "data.table", "simple", "matrix"),
    group_size = "automatic",
    do_parallel = TRUE,
    cores = NA,
    verbose = TRUE,
    knn_params = NULL,
    set.seed = deprecated(),
    seed = 1234,
    bin_matrix = NULL,
    summarize = c("p.value", "adj.p.value"),
    return_gobject = FALSE) {
    # TODO align set.seed, set_seed, seed_number naming and usage across
    # packages
    # use only param seed. If NULL, set no seed. If !NULL set value as seed

    if (is_present(set.seed) && !is.function(set.seed)) {
        deprecate_warn(
            when = "4.0.3",
            what = "binSpect(set.seed)",
            with = "binSpect(seed)"
        )

        seed <- set.seed
        set.seed <- NULL
    }

    implementation <- match.arg(implementation)
    if (identical(implementation, "fast")) {
        # arguments only the legacy implementations honour: keep old calls
        # working by taking that path, and say so
        legacy_args <- c(
            spatial_network_k = !is.null(spatial_network_k),
            bin_matrix = !is.null(bin_matrix),
            calc_hub = isTRUE(calc_hub),
            do_fisher_test = isFALSE(do_fisher_test),
            kmeans_algo = !missing(kmeans_algo) && !identical(kmeans_algo[1L], "kmeans")
        )
        if (any(legacy_args)) {
            vmsg(.v = verbose, "binSpect: ",
                paste(names(legacy_args)[legacy_args], collapse = ", "),
                " needs the legacy implementation; using implementation = ",
                "\"data.table\"")
            implementation <- "data.table"
        }
    }
    if (identical(implementation, "fast")) {
        return(.binspect_fast(gobject,
            spat_unit = spat_unit, feat_type = feat_type,
            bin_method = bin_method, expression_values = expression_values,
            subset_feats = subset_feats,
            spatial_network_name = spatial_network_name,
            spatial_network_k = spatial_network_k,
            percentage_rank = percentage_rank, adjust_method = adjust_method,
            return_gobject = return_gobject
        ))
    }

    a <- get_args_list(keep = c(
        "spat_unit", "feat_type", "bin_method", "expression_values",
        "subset_feats", "reduce_network", "kmeans_algo",
        "nstart", "iter_max", "extreme_nr", "sample_nr",
        "percentage_rank", "do_fisher_test", "adjust_method",
        "calc_hub", "hub_min_int", "get_av_expr", "get_high_expr",
        "implementation", "group_size", "do_parallel", "cores", "seed",
        "verbose"
    ))

    if (!is.null(spatial_network_k)) {
        # list(), not c(): c() would flatten a vector or matrix argument
        # into one argument per element
        output <- do.call(binSpectMulti, args = c(a, list(
            gobject = gobject,
            spatial_network_k = spatial_network_k,
            knn_params = knn_params,
            summarize = summarize
        )))
    } else {
        output <- do.call(binSpectSingle, args = c(a, list(
            gobject = gobject,
            spatial_network_name = spatial_network_name,
            bin_matrix = bin_matrix
        )))
    }

    if (isTRUE(return_gobject)) {
        # Multi returns list(combined, simple); simple holds one p per feature
        result_dt <- if (!is.null(spatial_network_k)) {
            data.table::data.table(
                feats = output$simple$feats, pval = output$simple$p.val
            )
        } else {
            data.table::data.table(
                feats = output$feats, pval = output$adj.p.value
            )
        }
        data.table::setnames(result_dt, old = "pval", new = "binSpect.pval")
        gobject <- addFeatMetadata(
            gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = result_dt,
            by_column = TRUE,
            column_feat_ID = "feats"
        )
        return(gobject)
    } else {
        return(output)
    }
}





#' @describeIn binSpect binSpect for a single spatial network and a provided
#' expression matrix
#' @param expression_matrix expression matrix
#' @param spatial_network spatial network in data.table format
#' @export
binSpectSingleMatrix <- function(
        expression_matrix,
        spatial_network = NULL,
        bin_matrix = NULL,
        bin_method = c("kmeans", "rank"),
        subset_feats = NULL,
        kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
        nstart = 3,
        iter_max = 10,
        extreme_nr = 50,
        sample_nr = 50,
        percentage_rank = 30,
        do_fisher_test = TRUE,
        adjust_method = "fdr",
        calc_hub = FALSE,
        hub_min_int = 3,
        get_av_expr = TRUE,
        get_high_expr = TRUE,
        implementation = c("data.table", "simple", "matrix"),
        group_size = "automatic",
        do_parallel = TRUE,
        cores = NA,
        verbose = FALSE,
        set.seed = deprecated(),
        seed = 1234) {
    if (is_present(set.seed) && !is.function(set.seed)) {
        deprecate_warn(
            when = "4.0.3",
            what = "binSpectSingleMatrix(set.seed)",
            with = "binSpectSingleMatrix(seed)"
        )

        seed <- set.seed
        set.seed <- NULL
    }

    do_parallel <- as.logical(do_parallel)
    calc_hub <- as.logical(calc_hub)
    get_av_expr <- as.logical(get_av_expr)
    get_high_expr <- as.logical(get_high_expr)
    do_fisher_test <- as.logical(do_fisher_test)

    vmsg(.v = verbose, "\n This is the single parameter version of binSpect")


    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    # data.table: set global variable
    feats <- p.value <- estimate <- score <- NULL

    # set binarization method
    bin_method <- match.arg(bin_method, choices = c("kmeans", "rank"))

    # kmeans algorithm
    kmeans_algo <- match.arg(
        kmeans_algo,
        choices = c("kmeans", "kmeans_arma", "kmeans_arma_subset")
    )

    # implementation
    implementation <- match.arg(
        implementation,
        choices = c("data.table", "simple", "matrix")
    )


    # spatial network
    # TODO: verify binarization of spatial network
    if (is.null(spatial_network)) {
        stop("You need to provide a spatial network in data.table format to
            the 'spatial_network' parameter")
    }


    ## start binarization ##
    ## ------------------ ##

    if (!is.null(bin_matrix)) {
        # TODO: verify format of bin_matrix and compatibility with spatial
        # network
        bin_matrix <- bin_matrix
    } else {
        bin_matrix <- switch(bin_method,
            "kmeans" = kmeans_binarize_wrapper(
                expr_values = expression_matrix,
                subset_feats = subset_feats,
                kmeans_algo = kmeans_algo,
                nstart = nstart,
                iter_max = iter_max,
                extreme_nr = extreme_nr,
                sample_nr = sample_nr,
                # set.seed = set.seed,
                seed = seed
            ),
            "rank" = rank_binarize_wrapper(
                expr_values = expression_matrix,
                subset_feats = subset_feats,
                percentage_rank = percentage_rank
            )
        )
    }

    vmsg(.v = verbose, "\n 1. matrix binarization complete")

    ## start with enrichment ##
    ## --------------------- ##

    result <- switch(implementation,
        "simple" = {
            if (do_parallel) {
                warning("Parallel not yet implemented for simple.
                        Enrichment will default to serial.")
            }
            if (calc_hub) {
                warning("Hub calculation is not possible with the simple
                        implementation, change to matrix if required.")
            }

            .calc_spatial_enrichment_minimum(
                spatial_network = spatial_network,
                bin_matrix = bin_matrix,
                adjust_method = adjust_method,
                do_fisher_test = do_fisher_test
            )
        },
        "matrix" = .calc_spatial_enrichment_matrix(
            spatial_network = spatial_network,
            bin_matrix = bin_matrix,
            adjust_method = adjust_method,
            do_fisher_test = do_fisher_test,
            do_parallel = do_parallel,
            cores = cores,
            calc_hub = calc_hub,
            hub_min_int = hub_min_int,
            verbose = verbose
        ),
        "data.table" = .calc_spatial_enrichment_dt(
            bin_matrix = bin_matrix,
            spatial_network = spatial_network,
            calc_hub = calc_hub,
            hub_min_int = hub_min_int,
            group_size = group_size,
            do_fisher_test = do_fisher_test,
            adjust_method = adjust_method,
            cores = cores
        )
    )

    vmsg(.v = verbose, "\n 2. spatial enrichment test completed")



    ## start with average high expression ##
    ## ---------------------------------- ##

    if (get_av_expr) {
        # expression
        if (!is.null(subset_feats)) {
            expr_values <- expression_matrix[
                rownames(expression_matrix) %in% subset_feats,
            ]
        } else {
            expr_values <- expression_matrix
        }

        sel_expr_values <- expr_values * bin_matrix
        av_expr <- apply(sel_expr_values, MARGIN = 1, FUN = function(x) {
            mean(x[x > 0])
        })
        av_expr_DT <- data.table::data.table(
            feats = names(av_expr), av_expr = av_expr
        )
        result <- merge(result, av_expr_DT, by = "feats")

        vmsg(.v = verbose, "\n 3. (optional) average expression of high
            expressing cells calculated")
    }



    ## start with number of high expressing cells ##
    ## ------------------------------------------ ##

    if (get_high_expr) {
        high_expr <- rowSums(bin_matrix)
        high_expr_DT <- data.table::data.table(
            feats = names(high_expr), high_expr = high_expr
        )
        result <- merge(result, high_expr_DT, by = "feats")

        vmsg(.v = verbose, "\n 4. (optional) number of high expressing cells
            calculated")
    }


    # sort
    if (do_fisher_test) {
        data.table::setorder(result, -score)
    } else {
        data.table::setorder(result, -estimate)
    }

    return(result)
}



#' @describeIn binSpect binSpect for a single spatial network
#' @export
binSpectSingle <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        bin_method = c("kmeans", "rank"),
        expression_values = c("normalized", "scaled", "custom"),
        subset_feats = NULL,
        spatial_network_name = "Delaunay_network",
        reduce_network = FALSE,
        kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
        nstart = 3,
        iter_max = 10,
        extreme_nr = 50,
        sample_nr = 50,
        percentage_rank = 30,
        do_fisher_test = TRUE,
        adjust_method = "fdr",
        calc_hub = FALSE,
        hub_min_int = 3,
        get_av_expr = TRUE,
        get_high_expr = TRUE,
        implementation = c("data.table", "simple", "matrix"),
        group_size = "automatic",
        do_parallel = TRUE,
        cores = NA,
        verbose = TRUE,
        set.seed = deprecated(),
        seed = 1234,
        bin_matrix = NULL) {
    ## deprecated arguments

    if (is_present(set.seed) && !is.function(set.seed)) {
        deprecate_warn(
            when = "4.0.3",
            what = "binSpectSingle(set.seed)",
            with = "binSpectSingle(seed)"
        )

        seed <- set.seed
        set.seed <- NULL
    }

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

    ## 1. expression matrix
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )


    ## 2. spatial network
    spatial_network <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "networkDT"
    )
    if (is.null(spatial_network)) {
        stop(
            "spatial_network_name: ", spatial_network_name,
            " does not exist, create a spatial network first"
        )
    }

    # bidirectionalize canonical edges
    if (reduce_network == FALSE) {
        rev <- data.table::copy(spatial_network)
        data.table::setnames(rev, c("from", "to"), c("to", "from"))
        spatial_network <- unique(rbind(spatial_network, rev))
    }


    binSpectSingleMatrix(
        expression_matrix = expr_values,
        spatial_network = spatial_network,
        bin_matrix = bin_matrix,
        bin_method = bin_method,
        subset_feats = subset_feats,
        kmeans_algo = kmeans_algo,
        nstart = nstart,
        iter_max = iter_max,
        extreme_nr = extreme_nr,
        sample_nr = sample_nr,
        percentage_rank = percentage_rank,
        do_fisher_test = do_fisher_test,
        adjust_method = adjust_method,
        calc_hub = calc_hub,
        hub_min_int = hub_min_int,
        get_av_expr = get_av_expr,
        get_high_expr = get_high_expr,
        implementation = implementation,
        group_size = group_size,
        do_parallel = do_parallel,
        cores = cores,
        verbose = verbose,
        seed = seed
    )
}





#' @describeIn binSpect binSpect for multiple spatial kNN networks
#' @export
binSpectMulti <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        bin_method = c("kmeans", "rank"),
        expression_values = c("normalized", "scaled", "custom"),
        subset_feats = NULL,
        spatial_network_k = c(5, 10, 20),
        reduce_network = FALSE,
        kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
        nstart = 3,
        iter_max = 10,
        extreme_nr = 50,
        sample_nr = 50,
        percentage_rank = c(10, 30),
        do_fisher_test = TRUE,
        adjust_method = "fdr",
        calc_hub = FALSE,
        hub_min_int = 3,
        get_av_expr = TRUE,
        get_high_expr = TRUE,
        implementation = c("data.table", "simple", "matrix"),
        group_size = "automatic",
        do_parallel = TRUE,
        cores = NA,
        verbose = TRUE,
        knn_params = NULL,
        set.seed = deprecated(),
        seed = 1234,
        summarize = c("p.value", "adj.p.value")) {
    ## deprecated arguments
    if (is_present(set.seed) && !is.function(set.seed)) {
        deprecate_warn(
            when = "4.0.3",
            what = "binSpectMulti(set.seed)",
            with = "binSpectMulti(seed)"
        )

        seed <- set.seed
        set.seed <- NULL
    }

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

    if (verbose == TRUE) {
        message("This is the multi parameter version of binSpect")
    }

    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    # check bin_method
    bin_method <- match.arg(bin_method, choices = c("kmeans", "rank"))

    # summarization level
    summarize <- match.arg(summarize, choices = c("p.value", "adj.p.value"))

    ## bin method rank
    if (bin_method == "rank") {
        total_trials <- length(spatial_network_k) * length(percentage_rank)
        result_list <- vector(mode = "list", length = total_trials)
        i <- 1

        for (k in spatial_network_k) {
            if (is.null(knn_params)) {
                knn_params <- list(minimum_k = 1)
            }
            temp_gobject <- do.call("createSpatialKNNnetwork", c(
                gobject = gobject,
                spat_unit = spat_unit,
                name = "temp_knn_network",
                k = k,
                knn_params
            ))

            for (rank_i in percentage_rank) {
                if (verbose == TRUE) {
                    cat("Run for k = ", k, " and rank % = ", rank_i)
                }

                result <- binSpectSingle(
                    gobject = temp_gobject,
                    feat_type = feat_type,
                    spat_unit = spat_unit,
                    bin_method = bin_method,
                    expression_values = expression_values,
                    subset_feats = subset_feats,
                    spatial_network_name = "temp_knn_network",
                    reduce_network = reduce_network,
                    kmeans_algo = kmeans_algo,
                    percentage_rank = rank_i,
                    do_fisher_test = do_fisher_test,
                    adjust_method = adjust_method,
                    calc_hub = calc_hub,
                    hub_min_int = hub_min_int,
                    get_av_expr = get_av_expr,
                    get_high_expr = get_high_expr,
                    implementation = implementation,
                    group_size = group_size,
                    do_parallel = do_parallel,
                    cores = cores,
                    verbose = verbose,
                    # set.seed = set.seed,
                    seed = seed
                )

                result_list[[i]] <- result
                i <- i + 1
            }
        }
        combined_result <- data.table::rbindlist(result_list)
    } else if (bin_method == "kmeans") {
        ## bin method kmeans
        total_trials <- length(spatial_network_k)
        result_list <- vector(mode = "list", length = total_trials)
        i <- 1

        ## expression matrix
        values <- match.arg(
            expression_values,
            unique(c("normalized", "scaled", "custom", expression_values))
        )
        expr_values <- getExpression(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            values = values,
            output = "matrix"
        )


        # pre-calculate bin_matrix once
        bin_matrix <- kmeans_binarize_wrapper(
            expr_values = expr_values,
            subset_feats = subset_feats,
            kmeans_algo = kmeans_algo,
            nstart = nstart,
            iter_max = iter_max,
            extreme_nr = extreme_nr,
            sample_nr = sample_nr,
            # set.seed = set.seed,
            seed = seed
        )

        for (k in spatial_network_k) {
            if (is.null(knn_params)) {
                knn_params <- list(minimum_k = 1)
            }
            temp_gobject <- do.call("createSpatialKNNnetwork", c(
                gobject = gobject,
                spat_unit = spat_unit,
                name = "temp_knn_network",
                k = k, knn_params
            ))

            if (verbose == TRUE) cat("Run for k = ", k)

            result <- binSpectSingle(
                gobject = temp_gobject,
                feat_type = feat_type,
                spat_unit = spat_unit,
                bin_method = bin_method,
                expression_values = expression_values,
                subset_feats = subset_feats,
                spatial_network_name = "temp_knn_network",
                reduce_network = reduce_network,
                kmeans_algo = kmeans_algo,
                nstart = nstart,
                iter_max = iter_max,
                extreme_nr = extreme_nr,
                sample_nr = sample_nr,
                do_fisher_test = do_fisher_test,
                adjust_method = adjust_method,
                calc_hub = calc_hub,
                hub_min_int = hub_min_int,
                get_av_expr = get_av_expr,
                get_high_expr = get_high_expr,
                implementation = implementation,
                group_size = group_size,
                do_parallel = do_parallel,
                cores = cores,
                verbose = verbose,
                # set.seed = set.seed,
                seed = seed,
                bin_matrix = bin_matrix
            )

            result_list[[i]] <- result
            i <- i + 1
        }

        combined_result <- data.table::rbindlist(result_list)
    }


    # data.table variables
    feats <- V1 <- p.val <- NULL

    ## merge results into 1 p-value per feat ##
    simple_result <- combined_result[, sum(log(get(summarize))), by = feats]
    simple_result[, V1 := V1 * -2]
    simple_result[, p.val := stats::pchisq(
        q = V1, df = 2 * total_trials, log.p = FALSE, lower.tail = FALSE
    )]

    return(list(
        combined = combined_result, simple = simple_result[, .(feats, p.val)]
    ))
}



# not exported, so not linked to binSpect

#' @title binSpectMultiMatrix
#' @name binSpectMultiMatrix
#' @description binSpect for a single spatial network and a provided
#' expression matrix
#' @param expression_matrix expression matrix
#' @param spatial_networks list of spatial networks in data.table format
#' @param bin_method method to binarize gene expression
#' @param subset_feats only select a subset of features to test
#' @param kmeans_algo kmeans algorithm to use
#' (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank vector of percentages of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use
#' (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high
#' expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation enrichment implementation (data.table, simple, matrix)
#' @param group_size number of genes to process together with data.table
#' implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @param knn_params list of parameters to create spatial kNN network
#' @param set.seed deprecated
#' @param seed sets value as seed before kmeans binarization. If NULL, no seed
#' is set.
#' @param summarize summarize the p-values or adjusted p-values
#' @returns data.table with results
binSpectMultiMatrix <- function(
        expression_matrix,
        spatial_networks,
        bin_method = c("kmeans", "rank"),
        subset_feats = NULL,
        kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
        nstart = 3,
        iter_max = 10,
        extreme_nr = 50,
        sample_nr = 50,
        percentage_rank = c(10, 30),
        do_fisher_test = TRUE,
        adjust_method = "fdr",
        calc_hub = FALSE,
        hub_min_int = 3,
        get_av_expr = TRUE,
        get_high_expr = TRUE,
        implementation = c("data.table", "simple", "matrix"),
        group_size = "automatic",
        do_parallel = TRUE,
        cores = NA,
        verbose = TRUE,
        knn_params = NULL,
        set.seed = deprecated(),
        seed = 1234,
        summarize = c("p.value", "adj.p.value")) {
    if (is_present(set.seed) && !is.function(set.seed)) {
        deprecate_warn(
            when = "4.0.3",
            what = "binSpectMultiMatrix(set.seed)",
            with = "binSpectMultiMatrix(seed)"
        )

        seed <- set.seed
        set.seed <- NULL
    }


    if (verbose == TRUE) {
        message("This is the multi parameter version of binSpect")
    }

    # set number of cores automatically, but with limit of 10
    cores <- determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    # check bin_method
    bin_method <- match.arg(bin_method, choices = c("kmeans", "rank"))

    # summarization level
    summarize <- match.arg(summarize, choices = c("p.value", "adj.p.value"))


    ## bin method rank
    if (bin_method == "rank") {
        total_trials <- length(spatial_networks) * length(percentage_rank)
        result_list <- vector(mode = "list", length = total_trials)
        i <- 1

        for (k in seq_along(spatial_networks)) {
            for (rank_i in percentage_rank) {
                if (verbose == TRUE) {
                    cat("Run for spatial network ", k, " and rank % = ", rank_i)
                }

                result <- binSpectSingleMatrix(
                    expression_matrix = expression_matrix,
                    spatial_network = spatial_networks[[k]],
                    bin_method = bin_method,
                    subset_feats = subset_feats,
                    kmeans_algo = kmeans_algo,
                    percentage_rank = rank_i,
                    do_fisher_test = do_fisher_test,
                    adjust_method = adjust_method,
                    calc_hub = calc_hub,
                    hub_min_int = hub_min_int,
                    get_av_expr = get_av_expr,
                    get_high_expr = get_high_expr,
                    implementation = implementation,
                    group_size = group_size,
                    do_parallel = do_parallel,
                    cores = cores,
                    verbose = verbose,
                    # set.seed = set.seed,
                    seed = seed
                )

                result_list[[i]] <- result
                i <- i + 1
            }
        }
        combined_result <- data.table::rbindlist(result_list)
    } else if (bin_method == "kmeans") {
        ## bin method kmeans
        total_trials <- length(spatial_networks)
        result_list <- vector(mode = "list", length = total_trials)
        i <- 1


        # pre-calculate bin_matrix once
        bin_matrix <- kmeans_binarize_wrapper(
            expr_values = expression_matrix,
            subset_feats = subset_feats,
            kmeans_algo = kmeans_algo,
            nstart = nstart,
            iter_max = iter_max,
            extreme_nr = extreme_nr,
            sample_nr = sample_nr,
            # set.seed = set.seed,
            seed = seed
        )

        for (k in seq_along(spatial_networks)) {
            if (verbose == TRUE) cat("Run for spatial network = ", k)

            result <- binSpectSingleMatrix(
                expression_matrix = expression_matrix,
                bin_matrix = bin_matrix,
                spatial_network = spatial_networks[[k]],
                bin_method = bin_method,
                subset_feats = subset_feats,
                kmeans_algo = kmeans_algo,
                nstart = nstart,
                iter_max = iter_max,
                extreme_nr = extreme_nr,
                sample_nr = sample_nr,
                do_fisher_test = do_fisher_test,
                adjust_method = adjust_method,
                calc_hub = calc_hub,
                hub_min_int = hub_min_int,
                get_av_expr = get_av_expr,
                get_high_expr = get_high_expr,
                implementation = implementation,
                group_size = group_size,
                do_parallel = do_parallel,
                cores = cores,
                verbose = verbose,
                # set.seed = set.seed,
                seed = seed
            )

            result_list[[i]] <- result
            i <- i + 1
        }

        combined_result <- data.table::rbindlist(result_list)
    }


    # data.table variables
    feats <- V1 <- p.val <- NULL

    ## merge results into 1 p-value per feat ##
    simple_result <- combined_result[, sum(log(get(summarize))), by = feats]
    simple_result[, V1 := V1 * -2]
    simple_result[, p.val := stats::pchisq(
        q = V1, df = 2 * total_trials, log.p = FALSE, lower.tail = FALSE
    )]

    return(list(
        combined = combined_result, simple = simple_result[, .(feats, p.val)]
    ))
}





# * ####
## PCA spatial patterns ####

#' @title detectSpatialPatterns
#' @name detectSpatialPatterns
#' @description Identify spatial patterns through PCA on average expression
#' in a spatial grid.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param spatial_grid_name name of spatial grid to use
#' (default = 'spatial_grid')
#' @param min_cells_per_grid minimum number of cells in a grid to be considered
#' @param scale_unit scale features
#' @param ncp number of principal components to calculate
#' @param show_plot show plots
#' @param PC_zscore minimum z-score of variance explained by a PC
#' @returns spatial pattern object 'spatPatObj'
#' @details
#' Steps to identify spatial patterns:
#' \itemize{
#'   * 1. average gene expression for cells within a grid, see createSpatialGrid
#'   * 2. perform PCA on the average grid expression profiles
#'   * 3. convert variance of principal components (PCs) to z-scores and
#'   select PCs based on a z-score threshold
#' }
#' @export
detectSpatialPatterns <- function(
        gobject,
        expression_values = c("normalized", "scaled", "custom"),
        spatial_grid_name = "spatial_grid",
        min_cells_per_grid = 4,
        scale_unit = FALSE,
        ncp = 100,
        show_plot = TRUE,
        PC_zscore = 1.5) {
    ############################################################################
    stop(wrap_txt(
        "This function has not been updated for use with the current version
    of Giotto.
    See details:
    https://github.com/giotto-suite/Giotto/issues/666#issuecomment-1540447537",
        errWidth = TRUE
    ))
    ############################################################################
    # expression values to be used
    values <- match.arg(expression_values, c("normalized", "scaled", "custom"))
    expr_values <- getExpression(
        gobject = gobject,
        values = values,
        output = "matrix"
    )


    # spatial grid and spatial locations
    if (is.null(slot(gobject, "spatial_grid"))) {
        stop("you need to create a spatial grid, see createSpatialGrid(),
            for this function to work")
    }
    if (!spatial_grid_name %in% list_spatial_grids_names(gobject = gobject)) {
        stop("you need to provide an existing spatial grid name for this
            function to work")
    }

    spatial_grid <- getSpatialGrid(gobject, name = spatial_grid_name)

    # annotate spatial locations with spatial grid information
    spatial_locs <- getSpatialLocations(gobject,
        name = "raw",
        output = "data.table",
        copy_obj = TRUE
    )

    if (all(c("sdimx", "sdimy", "sdimz") %in% colnames(spatial_locs))) {
        spatial_locs <- annotate_spatlocs_with_spatgrid_3D(
            spatloc = spatial_locs, spatgrid = spatial_grid
        )
    } else if (all(c("sdimx", "sdimy") %in% colnames(spatial_locs))) {
        spatial_locs <- annotate_spatlocs_with_spatgrid_2D(
            spatloc = spatial_locs, spatgrid = spatial_grid
        )
    }


    # data.table variables
    gr_loc <- zscore <- variance.percent <- loc_ID <- gene_ID <- NULL

    # filter grid, minimum number of cells per grid
    cells_per_grid <- sort(table(spatial_locs$gr_loc))
    cells_per_grid <- cells_per_grid[cells_per_grid >= min_cells_per_grid]
    loc_names <- names(cells_per_grid)

    # average expression per grid
    loc_av_expr_list <- list()
    for (loc_name in loc_names) {
        loc_cell_IDs <- spatial_locs[gr_loc == loc_name]$cell_ID
        subset_expr <- expr_values[, colnames(expr_values) %in% loc_cell_IDs]
        if (is.vector(subset_expr) == TRUE) {
            loc_av_expr <- subset_expr
        } else {
            loc_av_expr <- rowMeans(subset_expr)
        }
        loc_av_expr_list[[loc_name]] <- loc_av_expr
    }
    loc_av_expr_matrix <- do.call("cbind", loc_av_expr_list)

    # START TEST
    loc_av_expr_matrix <- as.matrix(loc_av_expr_matrix)
    # STOP

    # perform pca on grid matrix
    mypca <- FactoMineR::PCA(
        X = t(loc_av_expr_matrix),
        scale.unit = scale_unit,
        ncp = ncp,
        graph = FALSE
    )

    # screeplot
    screeplot <- factoextra::fviz_eig(mypca, addlabels = TRUE, ylim = c(0, 50))
    if (show_plot == TRUE) {
        print(screeplot)
    }

    # select variable PCs
    eig.val <- factoextra::get_eigenvalue(mypca)
    eig.val_DT <- data.table::as.data.table(eig.val)
    eig.val_DT$names <- rownames(eig.val)
    eig.val_DT[, zscore := scale(variance.percent)]
    eig.val_DT[, rank := rank(variance.percent)]
    dims_to_keep <- eig.val_DT[zscore > PC_zscore]$names


    # if no dimensions are kept, return message
    if (is.null(dims_to_keep) | length(dims_to_keep) < 1) {
        message("no PC dimensions retained, lower the PC zscore")
    }

    # coordinates for cells
    pca_matrix <- mypca$ind$coord
    if (length(dims_to_keep) == 1) {
        pca_matrix_DT <- data.table::data.table(
            "dimkeep" = pca_matrix[, 1],
            loc_ID = colnames(loc_av_expr_matrix)
        )
        data.table::setnames(pca_matrix_DT, old = "dimkeep", dims_to_keep)
    } else {
        pca_matrix_DT <- data.table::as.data.table(pca_matrix[
            , seq_along(dims_to_keep)
        ])
        pca_matrix_DT[, loc_ID := colnames(loc_av_expr_matrix)]
    }


    # correlation of genes with PCs
    feat_matrix <- mypca$var$cor
    if (length(dims_to_keep) == 1) {
        feat_matrix_DT <- data.table::data.table(
            "featkeep" = feat_matrix[, 1],
            gene_ID = rownames(loc_av_expr_matrix)
        )
        data.table::setnames(feat_matrix_DT, old = "featkeep", dims_to_keep)
    } else {
        feat_matrix_DT <- data.table::as.data.table(feat_matrix[
            , seq_along(dims_to_keep)
        ])
        feat_matrix_DT[, gene_ID := rownames(loc_av_expr_matrix)]
    }


    spatPatObject <- list(
        pca_matrix_DT = pca_matrix_DT,
        feat_matrix_DT = feat_matrix_DT,
        spatial_grid = spatial_grid
    )

    class(spatPatObject) <- append("spatPatObj", class(spatPatObject))

    return(spatPatObject)
}



#' @title showPattern2D
#' @name showPattern2D
#' @description show patterns for 2D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of ggplot
#' @param point_size size of points
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see
#' \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns ggplot
#' @export
showPattern2D <- function(
        gobject,
        spatPatObj,
        dimension = 1,
        trim = c(0.02, 0.98),
        background_color = "white",
        grid_border_color = "grey",
        show_legend = TRUE,
        point_size = 1,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "showPattern2D") {
    if (!"spatPatObj" %in% class(spatPatObj)) {
        stop("spatPatObj needs to be the output from detectSpatialPatterns")
    }

    # select PC and subset data
    selected_PC <- paste0("Dim.", dimension)
    PC_DT <- spatPatObj$pca_matrix_DT
    if (!selected_PC %in% colnames(PC_DT)) {
        stop("This dimension was not found in the spatial pattern object")
    }
    PC_DT <- PC_DT[, c(selected_PC, "loc_ID"), with = FALSE]

    # annotate grid with PC values
    annotated_grid <- merge(
        spatPatObj$spatial_grid,
        by.x = "gr_name", PC_DT, by.y = "loc_ID"
    )

    # trim PC values
    if (!is.null(trim)) {
        boundaries <- stats::quantile(annotated_grid[[
            selected_PC
        ]], probs = trim)
        annotated_grid[[selected_PC]][annotated_grid[[
            selected_PC
        ]] < boundaries[1]] <- boundaries[1]
        annotated_grid[[selected_PC]][annotated_grid[[
            selected_PC
        ]] > boundaries[2]] <- boundaries[2]
    }

    # 2D-plot
    #


    dpl <- ggplot2::ggplot()
    dpl <- dpl + ggplot2::theme_bw()
    dpl <- dpl + ggplot2::geom_tile(
        data = annotated_grid,
        GiottoVisuals::aes_string2(
            x = "x_start", y = "y_start", fill = selected_PC),
        color = grid_border_color, 
        show.legend = show_legend
    )
    dpl <- dpl + ggplot2::scale_fill_gradient2(
        "low" = "darkblue", mid = "white", high = "darkred", midpoint = 0,
        guide = guide_legend(title = "")
    )
    dpl <- dpl + ggplot2::theme(
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        panel.background = element_rect(fill = background_color),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
    )
    dpl <- dpl + ggplot2::labs(x = "x coordinates", y = "y coordinates")


    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = dpl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}

#' @title showPattern
#' @name showPattern
#' @description show patterns for 2D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @inheritDotParams showPattern2D -gobject -spatPatObj
#' @returns ggplot
#' @seealso \code{\link{showPattern2D}}
#' @export
showPattern <- function(gobject, spatPatObj, ...) {
    showPattern2D(gobject = gobject, spatPatObj = spatPatObj, ...)
}

#' @title showPattern3D
#' @name showPattern3D
#' @description show patterns for 3D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of plot
#' @param point_size adjust the point size
#' @param axis_scale scale the axis
#' @param custom_ratio cutomize the scale of the axis
#' @param x_ticks the tick number of x_axis
#' @param y_ticks the tick number of y_axis
#' @param z_ticks the tick number of z_axis
#' @param show_plot show plot
#' @param return_plot return plot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see
#' \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns plotly
#' @export
showPattern3D <- function(
        gobject,
        spatPatObj,
        dimension = 1,
        trim = c(0.02, 0.98),
        background_color = "white",
        grid_border_color = "grey",
        show_legend = TRUE,
        point_size = 1,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "showPattern3D") {
    package_check("plotly", repository = "CRAN:plotly")

    # data.table variables
    center_x <- x_start <- x_end <- center_y <- y_start <- y_end <-
        center_z <- z_start <- z_end <- NULL

    if (!"spatPatObj" %in% class(spatPatObj)) {
        stop("spatPatObj needs to be the output from detectSpatialPatterns")
    }

    # select PC and subset data
    selected_PC <- paste0("Dim.", dimension)
    PC_DT <- spatPatObj$pca_matrix_DT
    if (!selected_PC %in% colnames(PC_DT)) {
        stop("This dimension was not found in the spatial pattern object")
    }
    PC_DT <- PC_DT[, c(selected_PC, "loc_ID"), with = FALSE]

    # annotate grid with PC values
    annotated_grid <- merge(
        spatPatObj$spatial_grid,
        by.x = "gr_name", PC_DT, by.y = "loc_ID"
    )

    # trim PC values
    if (!is.null(trim)) {
        boundaries <- stats::quantile(annotated_grid[[
            selected_PC
        ]], probs = trim)
        annotated_grid[[selected_PC]][annotated_grid[[
            selected_PC
        ]] < boundaries[1]] <- boundaries[1]
        annotated_grid[[selected_PC]][annotated_grid[[
            selected_PC
        ]] > boundaries[2]] <- boundaries[2]
    }


    annotated_grid <- data.table(annotated_grid)
    annotated_grid[, center_x := (x_start + x_end) / 2]
    annotated_grid[, center_y := (y_start + y_end) / 2]
    annotated_grid[, center_z := (z_start + z_end) / 2]


    axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))

    ratio <- plotly_axis_scale_3D(annotated_grid,
        sdimx = "center_x", sdimy = "center_y", sdimz = "center_z",
        mode = axis_scale, custom_ratio = custom_ratio
    )

    dpl <- plotly::plot_ly(
        type = "scatter3d",
        x = annotated_grid$center_x,
        y = annotated_grid$center_y, z = annotated_grid$center_z,
        color = annotated_grid[[selected_PC]], marker = list(size = point_size),
        mode = "markers", colors = c("darkblue", "white", "darkred")
    )
    dpl <- dpl %>% plotly::layout(scene = list(
        xaxis = list(title = "X", nticks = x_ticks),
        yaxis = list(title = "Y", nticks = y_ticks),
        zaxis = list(title = "Z", nticks = z_ticks),
        aspectmode = "manual",
        aspectratio = list(
            x = ratio[[1]],
            y = ratio[[2]],
            z = ratio[[3]]
        )
    ))
    dpl <- dpl %>% plotly::colorbar(
        title = paste(paste("dim.", dimension, sep = ""), "genes", sep = " ")
    )

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = dpl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}




#' @title showPatternGenes
#' @name showPatternGenes
#' @description show genes correlated with spatial patterns
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param point_size size of points
#' @param return_DT if TRUE, it will return the data.table used to generate
#' the plots
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see
#' \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns ggplot
#' @export
showPatternGenes <- function(
        gobject,
        spatPatObj,
        dimension = 1,
        top_pos_genes = 5,
        top_neg_genes = 5,
        point_size = 1,
        return_DT = FALSE,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "showPatternGenes") {
    # data.table variables
    gene_ID <- NULL

    if (!"spatPatObj" %in% class(spatPatObj)) {
        stop("spatPatObj needs to be the output from detectSpatialPatterns")
    }


    # select PC to use
    selected_PC <- paste0("Dim.", dimension)

    gene_cor_DT <- spatPatObj$feat_matrix_DT
    if (!selected_PC %in% colnames(gene_cor_DT)) {
        stop("This dimension was not found in the spatial pattern object")
    }
    gene_cor_DT <- gene_cor_DT[, c(selected_PC, "gene_ID"), with = FALSE]

    # order and subset
    gene_cor_DT <- gene_cor_DT[
        !is.na(get(selected_PC))
    ][order(get(selected_PC))]

    subset <- gene_cor_DT[
        c(seq_len(top_neg_genes), (nrow(
            gene_cor_DT
        ) - top_pos_genes):nrow(gene_cor_DT))
    ]
    subset[, gene_ID := factor(gene_ID, gene_ID)]

    ## return DT and make not plot ##
    if (return_DT == TRUE) {
        return(subset)
    }

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::geom_point(
        data = subset,
        GiottoVisuals::aes_string2(x = selected_PC, y = "gene_ID"), 
        size = point_size
    )
    pl <- pl + ggplot2::geom_vline(xintercept = 0, linetype = 2)
    pl <- pl + ggplot2::labs(x = "correlation", y = "", title = selected_PC)
    pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5))

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}


#' @title selectPatternGenes
#' @name selectPatternGenes
#' @description Select genes correlated with spatial patterns
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimensions dimensions to identify correlated genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param min_pos_cor Minimum positive correlation score to include a gene.
#' @param min_neg_cor Minimum negative correlation score to include a gene.
#' @param return_top_selection only return selection based on correlation
#' criteria (boolean)
#' @returns Data.table with genes associated with selected dimension (PC).
#' @details Description.
#' @export
selectPatternGenes <- function(
        spatPatObj,
        dimensions = 1:5,
        top_pos_genes = 10,
        top_neg_genes = 10,
        min_pos_cor = 0.5,
        min_neg_cor = -0.5,
        return_top_selection = FALSE) {
    if (!"spatPatObj" %in% class(spatPatObj)) {
        stop("spatPatObj needs to be the output from detectSpatialPatterns")
    }

    # data.table variables
    top_pos_rank <- value <- top_neg_rank <- topvalue <- gene_ID <-
        variable <- NULL


    # select PC to use
    selected_PCs <- paste0("Dim.", dimensions)
    gene_cor_DT <- spatPatObj$feat_matrix_DT
    if (any(selected_PCs %in% colnames(gene_cor_DT) == FALSE)) {
        stop("not all dimensions were found back")
    }
    gene_cor_DT <- gene_cor_DT[, c(selected_PCs, "gene_ID"), with = FALSE]

    # melt and select
    gene_cor_DT_m <- data.table::melt.data.table(
        gene_cor_DT,
        id.vars = "gene_ID"
    )
    gene_cor_DT_m[, top_pos_rank := rank(value), by = "variable"]
    gene_cor_DT_m[, top_neg_rank := rank(-value), by = "variable"]
    selection <- gene_cor_DT_m[
        top_pos_rank %in% seq_len(top_pos_genes) |
            top_neg_rank %in% seq_len(top_neg_genes)
    ]

    # filter on min correlation
    selection <- selection[value > min_pos_cor | value < min_neg_cor]

    # return all the top correlated genes + information
    if (return_top_selection == TRUE) {
        return(selection)
    }

    # remove duplicated genes by only retaining the most correlated dimension
    selection[, topvalue := max(abs(value)), by = "gene_ID"]
    uniq_selection <- selection[value == topvalue]

    # add other genes back
    output_selection <- uniq_selection[, .(gene_ID, variable)]
    other_genes <- gene_cor_DT[!gene_ID %in% output_selection$gene_ID][[
        "gene_ID"
    ]]
    other_genes_DT <- data.table::data.table(
        gene_ID = other_genes, variable = "noDim"
    )

    comb_output_genes <- rbind(output_selection, other_genes_DT)
    setnames(comb_output_genes, "variable", "patDim")

    return(comb_output_genes)
}

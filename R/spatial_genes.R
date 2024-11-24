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
    verbose <- feats <- p.value <- estimate <- adj.p.value <- score <- NULL

    # convert spatial network data.table to spatial matrix
    dc_spat_network <- data.table::dcast.data.table(
        spatial_network,
        formula = to ~ from, value.var = "distance", fill = 0
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
#' @param implementation enrichment implementation (data.table, simple, matrix)
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
    implementation = c("data.table", "simple", "matrix"),
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
        output <- do.call(binSpectMulti, args = c(a,
            gobject = gobject,
            spatial_network_k = spatial_network_k,
            knn_params = knn_params,
            summarize = summarize
        ))
    } else {
        output <- do.call(binSpectSingle, args = c(a,
            gobject = gobject,
            spatial_network_name = spatial_network_name,
            bin_matrix = bin_matrix
        ))
    }

    if (isTRUE(return_gobject)) {
        result_dt <- data.table::data.table(
            feats = output$feats, pval = output$adj.p.value
        )
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

    # convert to full network
    if (reduce_network == FALSE) {
        spatial_network <- convert_to_full_spatial_network(spatial_network)
        data.table::setnames(
            spatial_network, c("source", "target"), c("from", "to")
        )
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
        summarize = c("adj.p.value", "p.value")) {
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
    summarize <- match.arg(summarize, choices = c("adj.p.value", "p.value"))

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
        q = V1, df = total_trials, log.p = FALSE, lower.tail = FALSE
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
        summarize = c("adj.p.value", "p.value")) {
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
    summarize <- match.arg(summarize, choices = c("adj.p.value", "p.value"))


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
        q = V1, df = total_trials, log.p = FALSE, lower.tail = FALSE
    )]

    return(list(
        combined = combined_result, simple = simple_result[, .(feats, p.val)]
    ))
}





#' @title silhouetteRank
#' @name silhouetteRank
#' @description Previously: calculate_spatial_genes_python. This method
#' computes a silhouette score per gene based on the
#' spatial distribution of two partitions of cells
#' (expressed L1, and non-expressed L0).
#' Here, rather than L2 Euclidean norm, it uses a rank-transformed,
#' exponentially weighted
#' function to represent the local physical distance between two cells.
#' New multi aggregator implementation can be found at
#' \code{\link{silhouetteRankTest}}
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param metric distance metric to use
#' @param subset_genes only run on this subset of genes
#' @param rbp_p fractional binarization threshold
#' @param examine_top top fraction to evaluate with silhouette
#' @param python_path specify specific path to python if required
#' @returns data.table with spatial scores
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' silhouetteRank(g)
#' @export
silhouetteRank <- function(
        gobject,
        expression_values = c("normalized", "scaled", "custom"),
        metric = "euclidean",
        subset_genes = NULL,
        rbp_p = 0.95,
        examine_top = 0.3,
        python_path = NULL) {
    # expression values
    values <- match.arg(expression_values, c("normalized", "scaled", "custom"))
    expr_values <- getExpression(
        gobject = gobject,
        values = values,
        output = "matrix"
    )

    # subset genes
    if (!is.null(subset_genes)) {
        subset_genes <- subset_genes[subset_genes %in% gobject@feat_ID]
        expr_values <- expr_values[rownames(expr_values) %in% subset_genes, ]
    }


    # data.table variables
    sdimx <- sdimy <- NULL

    # spatial locations
    spatlocs <- getSpatialLocations(gobject,
        spat_unit = "cell",
        name = "raw",
        output = "data.table",
        copy_obj = TRUE
    )
    spatlocs <- as.matrix(spatlocs[, .(sdimx, sdimy)])

    # python path
    if (is.null(python_path)) {
        python_path <- readGiottoInstructions(gobject, param = "python_path")
    }

    ## prepare python path and louvain script
    reticulate::use_python(required = TRUE, python = python_path)
    python_silh_function <- system.file(
        "python", "python_spatial_genes.py",
        package = "Giotto"
    )
    reticulate::source_python(file = python_silh_function)

    output_python <- python_spatial_genes(
        spatial_locations = spatlocs,
        expression_matrix = as.data.frame(as.matrix(expr_values)),
        metric = metric,
        rbp_p = rbp_p,
        examine_top = examine_top
    )

    # unlist output
    genes <- unlist(lapply(output_python, FUN = function(x) {
        y <- x[1][[1]]
    }))
    scores <- unlist(lapply(output_python, FUN = function(x) {
        y <- x[2][[1]]
    }))

    spatial_python_DT <- data.table::data.table(genes = genes, scores = scores)

    return(spatial_python_DT)
}




#' @title silhouetteRankTest
#' @name silhouetteRankTest
#' @description Multi parameter aggregator version of
#' \code{\link{silhouetteRank}}
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param subset_genes only run on this subset of genes
#' @param overwrite_input_bin overwrite input bin
#' @param rbp_ps fractional binarization thresholds
#' @param examine_tops top fractions to evaluate with silhouette
#' @param matrix_type type of matrix
#' @param num_core number of cores to use
#' @param parallel_path path to GNU parallel function
#' @param output output directory
#' @param query_sizes size of query
#' @param verbose be verbose
#' @returns data.table with spatial scores
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' silhouetteRankTest(g)
#' @export
silhouetteRankTest <- function(
        gobject,
        expression_values = c("normalized", "scaled", "custom"),
        subset_genes = NULL,
        overwrite_input_bin = TRUE,
        rbp_ps = c(0.95, 0.99),
        examine_tops = c(0.005, 0.010, 0.050, 0.100, 0.300),
        matrix_type = "dissim",
        num_core = 4,
        parallel_path = "/usr/bin",
        output = NULL,
        query_sizes = 10L,
        verbose = FALSE) {
    # data.table variables
    cell_ID <- sdimx <- sdimy <- sdimz <- NULL

    ## test if R packages are installed
    # check envstats
    package_check(pkg_name = "EnvStats", repository = c("CRAN"))

    # check eva
    if ("eva" %in% rownames(installed.packages()) == FALSE) {
        stop(
            "\n package ", "eva", " is not yet installed \n",
            "To install: \n",
            "install.packages('eva')"
        )
    }

    ## test if python package is installed
    module_test <- reticulate::py_module_available("silhouetteRank")
    if (module_test == FALSE) {
        warning("silhouetteRank python module is not installed:
            install in the right environment or python path with:

            'pip install silhouetteRank'

            or from within R in the Giotto environment with:

            conda_path = reticulate::miniconda_path()
            conda_full_path = paste0(conda_path,'/','bin/conda')
            full_envname = paste0(conda_path,'/envs/giotto_env')
            reticulate::py_install(packages = 'silhouetteRank',
                       envname = full_envname,
                       method = 'conda',
                       conda = conda_full_path,
                       pip = TRUE,
                       python_version = '3.6')")
    }



    # expression values
    values <- match.arg(expression_values, c("normalized", "scaled", "custom"))
    expr_values <- getExpression(
        gobject = gobject,
        values = values,
        output = "matrix"
    )

    # subset genes
    if (!is.null(subset_genes)) {
        subset_genes <- subset_genes[subset_genes %in% gobject@gene_ID]
        expr_values <- expr_values[rownames(expr_values) %in% subset_genes, ]
    }

    # spatial locations
    spatlocs <- getSpatialLocations(gobject,
        name = "raw",
        output = "data.table",
        copy_obj = TRUE
    )

    ## save dir and log
    if (is.null(output)) {
        save_dir <- readGiottoInstructions(gobject, param = "save_dir")
        silh_output_dir <- paste0(save_dir, "/", "silhouetteRank_output/")
        if (!file.exists(silh_output_dir)) {
            dir.create(
                silh_output_dir,
                recursive = TRUE
            )
        }
    } else if (file.exists(output)) {
        silh_output_dir <- paste0(output, "/", "silhouetteRank_output/")
        if (!file.exists(silh_output_dir)) {
            dir.create(
                silh_output_dir,
                recursive = TRUE
            )
        }
    } else {
        silh_output_dir <- paste0(output, "/", "silhouetteRank_output/")
        if (!file.exists(silh_output_dir)) {
            dir.create(
                silh_output_dir,
                recursive = TRUE
            )
        }
    }

    # log directory
    log_dir <- paste0(silh_output_dir, "/", "logs/")
    if (!file.exists(log_dir)) dir.create(log_dir, recursive = TRUE)


    ## write spatial locations to .txt file
    if (ncol(spatlocs) == 3) {
        format_spatlocs <- spatlocs[, .(cell_ID, sdimx, sdimy)]
        colnames(format_spatlocs) <- c("ID", "x", "y")
    } else {
        format_spatlocs <- spatlocs[, .(cell_ID, sdimx, sdimy, sdimz)]
        colnames(format_spatlocs) <- c("ID", "x", "y", "z")
    }

    write.table(
        x = format_spatlocs, row.names = FALSE,
        file = paste0(silh_output_dir, "/", "format_spatlocs.txt"),
        quote = FALSE, sep = "\t"
    )

    spatlocs_path <- paste0(silh_output_dir, "/", "format_spatlocs.txt")


    silh_output_dir_norm <- normalizePath(silh_output_dir)
    expr_values_path_norm <- paste0(silh_output_dir_norm, "/", "expression.txt")

    data.table::fwrite(
        data.table::as.data.table(
            expr_values,
            keep.rownames = "gene"
        ),
        file = expr_values_path_norm,
        quote = FALSE,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
    )

    expr_values_path <- paste0(silh_output_dir, "/", "expression.txt")

    ## prepare python path and louvain script
    python_path <- readGiottoInstructions(gobject, param = "python_path")
    reticulate::use_python(required = TRUE, python = python_path)
    python_silh_function <- system.file(
        "python", "silhouette_rank_wrapper.py",
        package = "Giotto"
    )
    reticulate::source_python(file = python_silh_function)


    output_silh <- silhouette_rank(
        expr = expr_values_path,
        centroid = spatlocs_path,
        overwrite_input_bin = overwrite_input_bin,
        rbp_ps = rbp_ps,
        examine_tops = examine_tops,
        matrix_type = matrix_type,
        verbose = verbose,
        num_core = num_core,
        parallel_path = parallel_path,
        output = silh_output_dir,
        query_sizes = as.integer(query_sizes)
    )

    return(output_silh)
}






#' @title spatialDE
#' @name spatialDE
#' @description Compute spatial variable genes with spatialDE method
#' @param gobject Giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name for spatial locations
#' @param expression_values gene expression values to use
#' @param size size of plot
#' @param color low/medium/high color scheme for plot
#' @param sig_alpha alpha value for significance
#' @param unsig_alpha alpha value for unsignificance
#' @param python_path specify specific path to python if required
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see
#' \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns a list of data.frames with results and plot (optional)
#' @details This function is a wrapper for the SpatialDE method originally
#' implemented
#' in python. See publication \doi{10.1038/nmeth.4636}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' spatialDE(g)
#' @export
spatialDE <- function(
        gobject = NULL,
        feat_type = NULL,
        spat_unit = NULL,
        spat_loc_name = "raw",
        expression_values = c("raw", "normalized", "scaled", "custom"),
        size = c(4, 2, 1),
        color = c("blue", "green", "red"),
        sig_alpha = 0.5,
        unsig_alpha = 0.5,
        python_path = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "SpatialDE") {
    # test if SPARK is installed ##

    module_test <- reticulate::py_module_available("SpatialDE")
    if (module_test == FALSE) {
        warning("SpatialDE python module is not installed:
            install in the right environment or python path with:

            'pip install spatialde'

            or from within R in the Giotto environment with:

            conda_path = reticulate::miniconda_path()
            conda_full_path = paste0(conda_path,'/','bin/conda')
            full_envname = paste0(conda_path,'/envs/giotto_env')
            reticulate::py_install(packages = c('NaiveDE', 'patsy', 'SpatialDE'),
                                   envname = full_envname,
                                   method = 'conda',
                                   conda = conda_full_path,
                                   pip = TRUE,
                                   python_version = '3.6')")
    }


    # print message with information #
    message("using 'SpatialDE' for spatial gene/pattern detection. If used in
    published research, please cite:
    Svensson, Valentine, Sarah A. Teichmann, and Oliver Stegle.
    'SpatialDE: Identification of Spatially Variable Genes.'
    Nature Methods 15, no. 5 (May 2018): 343-46.
    https://doi.org/10.1038/nmeth.4636.")



    # data.table variables
    cell_ID <- NULL

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

    # expression
    values <- match.arg(
        expression_values, c("raw", "normalized", "scaled", "custom")
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )

    ## python path
    if (is.null(python_path)) {
        python_path <- readGiottoInstructions(gobject, param = "python_path")
    }

    ## source python file
    reticulate::use_python(required = TRUE, python = python_path)
    reader_path <- system.file(
        "python", "SpatialDE_wrapper.py",
        package = "Giotto"
    )
    reticulate::source_python(file = reader_path)

    ## get spatial locations
    spatial_locs <- getSpatialLocations(gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table"
    )
    spatial_locs <- as.data.frame(spatial_locs)
    rownames(spatial_locs) <- spatial_locs$cell_ID
    spatial_locs <- subset(spatial_locs, select = -cell_ID)

    ## run spatialDE
    Spatial_DE_results <- Spatial_DE(
        as.data.frame(t(as.matrix(expr_values))), spatial_locs
    )

    results <- as.data.frame(reticulate::py_to_r(Spatial_DE_results[[1]]))

    if (length(Spatial_DE_results) == 2) {
        ms_results <- as.data.frame(
            reticulate::py_to_r(Spatial_DE_results[[2]])
        )
        spatial_genes_results <- list(results, ms_results)
        names(spatial_genes_results) <- c("results", "ms_results")
    } else {
        spatial_genes_results <- results
        ms_results <- NULL
    }


    # print, return and save parameters
    show_plot <- ifelse(is.na(show_plot), readGiottoInstructions(
        gobject,
        param = "show_plot"
    ), show_plot)
    save_plot <- ifelse(is.na(save_plot), readGiottoInstructions(
        gobject,
        param = "save_plot"
    ), save_plot)
    return_plot <- ifelse(is.na(return_plot), readGiottoInstructions(
        gobject,
        param = "return_plot"
    ), return_plot)

    ## create plot
    if (isTRUE(show_plot) ||
        isTRUE(save_plot) ||
        isTRUE(return_plot)) {
        FSV_plot <- FSV_show(
            results = results,
            ms_results = ms_results,
            size = size,
            color = color,
            sig_alpha = sig_alpha,
            unsig_alpha = unsig_alpha
        )
    }

    ## print plot
    if (show_plot == TRUE) {
        print(FSV_plot)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call(
            "all_plots_save_function",
            c(list(
                gobject = gobject, plot_object = FSV_plot,
                default_save_name = default_save_name
            ), save_param)
        )
    }

    ## return results and plot (optional)
    if (return_plot == TRUE) {
        return(list(results = spatial_genes_results, plot = FSV_plot))
    } else {
        return(list(results = spatial_genes_results))
    }
}


#' @title spatialAEH
#' @name spatialAEH
#' @description Compute spatial variable genes with spatialDE method
#' @param gobject Giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name for spatial locations
#' @param SpatialDE_results results of \code{\link{spatialDE}} function
#' @param name_pattern name for the computed spatial patterns
#' @param expression_values gene expression values to use
#' @param pattern_num number of spatial patterns to look for
#' @param l lengthscale
#' @param python_path specify specific path to python if required
#' @param return_gobject show plot
#' @returns An updated giotto object
#' @details This function is a wrapper for the SpatialAEH method
#' implemented in the ...
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' spatialAEH(g)
#' @export
spatialAEH <- function(
        gobject = NULL,
        feat_type = NULL,
        spat_unit = NULL,
        spat_loc_name = "raw",
        SpatialDE_results = NULL,
        name_pattern = "AEH_patterns",
        expression_values = c("raw", "normalized", "scaled", "custom"),
        pattern_num = 6,
        l = 1.05,
        python_path = NULL,
        return_gobject = TRUE) {
    # data.table variables
    cell_ID <- NULL

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

    # expression
    values <- match.arg(
        expression_values, c("raw", "normalized", "scaled", "custom")
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values
    )

    ## python path
    if (is.null(python_path)) {
        python_path <- readGiottoInstructions(gobject, param = "python_path")
    }

    ## source python file
    reticulate::use_python(required = TRUE, python = python_path)
    reader_path <- system.file(
        "python", "SpatialDE_wrapper.py",
        package = "Giotto"
    )
    reticulate::source_python(file = reader_path)


    ## spatial locations
    spatial_locs <- getSpatialLocations(gobject,
        spat_unit = spat_unit,
        name = spat_loc_name
    )
    spatial_locs <- as.data.frame(spatial_locs)
    rownames(spatial_locs) <- spatial_locs$cell_ID
    spatial_locs <- subset(spatial_locs, select = -cell_ID)

    # extract results you need
    results <- SpatialDE_results[["results"]][["results"]]

    ## automatic expression histology
    AEH_results <- Spatial_DE_AEH(
        filterd_exprs = as.data.frame(t_flex(as.matrix(expr_values))),
        coordinates = spatial_locs,
        results = as.data.frame(results),
        pattern_num = pattern_num,
        l = l
    )
    histology_results <- as.data.frame(reticulate::py_to_r(AEH_results[[1]]))
    cell_pattern_score <- as.data.frame(reticulate::py_to_r(AEH_results[[2]]))

    spatial_pattern_results <- list(histology_results, cell_pattern_score)
    names(spatial_pattern_results) <- c(
        "histology_results", "cell_pattern_score"
    )


    if (return_gobject == TRUE) {
        dt_res <- data.table::as.data.table(
            spatial_pattern_results[["cell_pattern_score"]]
        )
        dt_res[["cell_ID"]] <- rownames(
            spatial_pattern_results[["cell_pattern_score"]]
        )
        gobject@spatial_enrichment[[name_pattern]] <- dt_res
        return(gobject)
    } else {
        return(list(results = spatial_pattern_results))
    }
}


#' @title FSV_show
#' @name FSV_show
#' @description Visualize spatial variable genes calculated by spatial_DE
#' @param results results calculated by spatial_DE
#' @param ms_results ms_results calculated by spatial_DE
#' @param size indicate different levels of qval
#' @param color indicate different SV features
#' @param sig_alpha transparency of significant genes
#' @param unsig_alpha transparency of unsignificant genes
#' @returns ggplot object
#' @keywords internal
FSV_show <- function(
        results,
        ms_results = NULL,
        size = c(4, 2, 1),
        color = c("blue", "green", "red"),
        sig_alpha = 0.5,
        unsig_alpha = 0.5) {
    results$FSV95conf <- 2 * sqrt(results$s2_FSV)
    results$intervals <- cut(
        results$FSV95conf, c(0, 1e-1, 1e0, Inf),
        label = FALSE
    )
    results$log_pval <- log10(results$pval)

    if (is.null(ms_results)) {
        results$model_bic <- results$model
    } else {
        results <- merge(results, ms_results[, c("g", "model")],
            by.x = "g", by.y = "g", all.x = TRUE,
            suffixes = (c(" ", "_bic"))
        )
    }

    results$model_bic <- factor(results$model_bic)
    results$intervals <- factor(results$intervals)


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()
    pl <- pl + ggplot2::geom_point(
        data = results[results$qval < 0.05, ],
        ggplot2::aes_string(
            x = "FSV", y = "log_pval", fill = "model_bic", size = "intervals"
        ),
        show.legend = TRUE, shape = 21, alpha = sig_alpha,
        stroke = 0.1, color = "black"
    ) +
        ggplot2::geom_point(
            data = results[results$qval > 0.05, ],
            ggplot2::aes_string(x = "FSV", y = "log_pval", size = "intervals"),
            show.legend = TRUE, shape = 21, alpha = unsig_alpha,
            fill = "black", # size = size[results_cp_ns$inftervals],
            stroke = 0.1, color = "black"
        ) +
        ggplot2::scale_size_manual(values = size, guide = FALSE) +
        ggplot2::scale_color_manual(values = color) +
        ggplot2::scale_fill_discrete(
            name = "Spatial Patterns",
            breaks = c("linear", "PER", "SE"),
            labels = c("linear", "periodical", "general")
        ) +
        ggplot2::geom_hline(yintercept = max(results[
            results$qval < 0.05,
        ]$log_pval), linetype = "dashed") +
        ggplot2::geom_text(ggplot2::aes(0.9, max(results[
            results$qval < 0.05,
        ]$log_pval),
        label = "FDR = 0.05", vjust = -1
        )) +
        ggplot2::scale_y_reverse()

    print(pl)
}




#' @title trendSceek
#' @name trendSceek
#' @description Compute spatial variable genes with trendsceek method
#' @param gobject Giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name for spatial locations
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to run trendsceek on
#' @param nrand An integer specifying the number of random resamplings of the
#' mark distribution as to create the null-distribution.
#' @param ncores An integer specifying the number of cores to be used by
#' BiocParallel
#' @param \dots Additional parameters to the
#' \code{\link[trendsceek]{trendsceek_test}} function
#' @returns data.frame with trendsceek spatial genes results
#' @details This function is a wrapper for the trendsceek_test method
#' implemented in the trendsceek package
#' Publication: \doi{10.1038/nmeth.4634}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' trendSceek(g)
#' @export
trendSceek <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        spat_loc_name = "raw",
        expression_values = c("normalized", "raw"),
        subset_genes = NULL,
        nrand = 100,
        ncores = 8,
        ...) {
    # verify if optional package is installed
    package_check(
        pkg_name = "trendsceek",
        repository = c("github"),
        github_repo = "edsgard/trendsceek"
    )

    # print message with information #
    message("using 'trendsceek' for spatial gene/pattern detection. If used
    in published research, please cite:
    Edsgard, Daniel, Per Johnsson, and Rickard Sandberg. 'Identification of
    Spatial Expression Trends in Single-Cell Gene Expression Data.'
    Nature Methods 15, no. 5 (May 2018): 339-42.
    https://doi.org/10.1038/nmeth.4634.")

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

    ## expression data
    values <- match.arg(expression_values, c("normalized", "raw"))
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values
    )

    ## normalization function
    if (values == "normalized") {
        log.fcn <- NA
    } else if (values == "raw") {
        log.fcn <- log10
    }

    ## subset genes
    if (!is.null(subset_genes)) {
        subset_genes <- subset_genes[subset_genes %in% gobject@gene_ID]
        expr_values <- expr_values[rownames(expr_values) %in% subset_genes, ]
    }


    ## initial locations

    # data.table variables
    cell_ID <- NULL

    spatial_locations <- getSpatialLocations(gobject,
        spat_unit = spat_unit,
        name = spat_loc_name
    )
    spatial_locations[, cell_ID := NULL]
    pp <- trendsceek::pos2pp(spatial_locations)

    ## initial gene counts
    pp <- trendsceek::set_marks(pp, as.matrix(expr_values), log.fcn = log.fcn)

    # eliminates running errors caused by too many zeros
    pp[["marks"]] <- pp[["marks"]] + 1e-7

    ## run trendsceek
    trendsceektest <- trendsceek::trendsceek_test(
        pp,
        nrand = nrand, ncores = ncores, ...
    )

    ## get final results
    trendsceektest <- trendsceektest$supstats_wide

    return(trendsceektest)
}




#' @title spark
#' @name spark
#' @description Compute spatially expressed genes with SPARK method
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name for spatial locations
#' @param percentage The percentage of cells that are expressed for analysis
#' @param min_count minimum number of counts for a gene to be included
#' @param expression_values type of values to use (raw by default)
#' @param num_core number of cores to use
#' @param covariates The covariates in experiments, i.e. confounding
#' factors/batch effect. Column name of giotto cell metadata.
#' @param return_object type of result to return (data.table or spark object)
#' @param \dots Additional parameters to the \code{\link[SPARK]{spark.vc}}
#' function
#' @returns data.table with SPARK spatial genes results or the SPARK object
#' @details This function is a wrapper for the method implemented in the
#' SPARK package:
#' \pkg{SPARK} package:
#'   1. **CreateSPARKObject** create a SPARK object from a giotto object
#'   2. **spark.vc** Fits the count-based spatial model to estimate the
#'   parameters, see \code{\link[SPARK]{spark.vc}} for additional parameters
#'   3. **spark.test** Testing multiple kernel matrices
#'
#' Publication: \doi{doi:10.1101/810903}
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' spark(g)
#' @export
spark <- function(
        gobject,
        spat_loc_name = "raw",
        feat_type = NULL,
        spat_unit = NULL,
        percentage = 0.1,
        min_count = 10,
        expression_values = "raw",
        num_core = 5,
        covariates = NULL,
        return_object = c("data.table", "spark"),
        ...) {
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

    # determine parameter
    return_object <- match.arg(return_object, c("data.table", "spark"))

    # data.table variables
    genes <- adjusted_pvalue <- combined_pvalue <- NULL

    ## test if SPARK is installed ##
    package_check(
        pkg_name = "SPARK",
        repository = c("github"),
        github_repo = "xzhoulab/SPARK"
    )


    # print message with information #
    message("using 'SPARK' for spatial gene/pattern detection. If used in
    published research, please cite:
    Sun, Shiquan, Jiaqiang Zhu, and Xiang Zhou. 'Statistical Analysis of
    Spatial Expression Pattern for Spatially Resolved Transcriptomic Studies.'
    BioRxiv, October 21, 2019, 810903. https://doi.org/10.1101/810903.")


    ## extract expression values from gobject
    expr <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values
    )

    ## extract coordinates from gobject
    locs <- getSpatialLocations(gobject,
        spat_unit = spat_unit,
        name = spat_loc_name
    )
    locs <- as.data.frame(locs)
    rownames(locs) <- colnames(expr)

    ## create SPARK object for analysis and filter out lowly expressed genes
    sobject <- SPARK::CreateSPARKObject(
        counts = expr,
        location = locs[, seq_len(2)],
        percentage = percentage,
        min_total_counts = min_count
    )

    ## total counts for each cell
    sobject@lib_size <- apply(sobject@counts, 2, sum)

    ## extract covariates ##
    if (!is.null(covariates)) {
        # first filter giotto object based on spark object
        filter_cell_ids <- colnames(sobject@counts)
        filter_gene_ids <- rownames(sobject@counts)
        tempgobject <- subsetGiotto(gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            cell_ids = filter_cell_ids,
            feat_ids = filter_gene_ids
        )

        metadata <- pDataDT(tempgobject)

        if (!covariates %in% colnames(metadata)) {
            warning(covariates, " was not found in the cell metadata of the
                    giotto object, will be set to NULL")
            covariates <- NULL
        } else {
            covariates <- metadata[[covariates]]
        }
    }

    ## Fit statistical model under null hypothesis
    sobject <- SPARK::spark.vc(sobject,
        covariates = covariates,
        lib_size = sobject@lib_size,
        num_core = num_core,
        verbose = FALSE,
        ...
    )

    ## test spatially expressed pattern genes
    ## calculating pval
    sobject <- SPARK::spark.test(sobject,
        check_positive = TRUE,
        verbose = FALSE
    )

    ## return results ##
    if (return_object == "spark") {
        return(sobject)
    } else if (return_object == "data.table") {
        DT_results <- data.table::as.data.table(sobject@res_mtest)
        gene_names <- rownames(sobject@counts)
        DT_results[, genes := gene_names]
        data.table::setorder(DT_results, adjusted_pvalue, combined_pvalue)
        return(DT_results)
    }
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
    https://github.com/drieslab/Giotto/issues/666#issuecomment-1540447537",
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
        aes_string(x = "x_start", y = "y_start", fill = selected_PC),
        color = grid_border_color, show.legend = show_legend
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
        aes_string(x = selected_PC, y = "gene_ID"), size = point_size
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







# ** ####
## Spatial co-expression ####
## ----------- ##

#' @title do_spatial_knn_smoothing
#' @name do_spatial_knn_smoothing
#' @description smooth gene expression over a kNN spatial network
#' @param expression_matrix gene expression values to use
#' @param subset_feats subset of features to use
#' @param spatial_network spatial network to use
#' @param b smoothing factor beteen 0 and 1 (default: automatic)
#' @returns matrix with smoothened gene expression values based on kNN
#' spatial network
#' @details This function will smoothen the gene expression values per cell
#' according to its neighbors in the selected spatial network. \cr
#' b is a smoothening factor that defaults to 1 - 1/k, where k is the median
#' number of k-neighbors in the selected spatial network. Setting b = 0 means
#' no smoothing and b = 1 means no contribution from its own expression.
#' @keywords internal
do_spatial_knn_smoothing <- function(
        expression_matrix,
        spatial_network,
        subset_feats = NULL,
        b = NULL) {
    # checks
    if (!is.null(b)) {
        if (b > 1 | b < 0) {
            stop("b needs to be between 0 (no spatial contribution) and 1
                (only spatial contribution)")
        }
    }

    # check if spatial network and expression matrix are compatible
    compatible_spatial_network(
        spatial_network = spatial_network,
        expression_matrix = expression_matrix
    )

    # matrix
    expr_values <- expression_matrix
    if (!is.null(subset_feats)) {
        expr_values <- expr_values[rownames(expr_values) %in% subset_feats, ]
    }

    # data.table variables
    feat_ID <- value <- NULL

    # merge spatial network with expression data
    expr_values_dt <- data.table::as.data.table(as.matrix(expr_values))
    expr_values_dt[, feat_ID := rownames(expr_values)]
    expr_values_dt_m <- data.table::melt.data.table(
        expr_values_dt,
        id.vars = "feat_ID", variable.name = "cell_ID"
    )


    # merge spatial network and matrix
    spatial_network_ext <- data.table::merge.data.table(
        spatial_network, expr_values_dt_m,
        by.x = "from", by.y = "cell_ID",
        allow.cartesian = TRUE
    )

    # calculate mean over all k-neighbours
    # exclude 0's?
    # trimmed mean?
    spatial_network_ext_smooth <- spatial_network_ext[
        , mean(value),
        by = c("to", "feat_ID")
    ]

    # convert back to matrix
    spatial_smooth_dc <- data.table::dcast.data.table(
        data = spatial_network_ext_smooth,
        formula = feat_ID ~ to,
        value.var = "V1"
    )
    spatial_smooth_matrix <- dt_to_matrix(spatial_smooth_dc)

    # if network was not fully connected, some cells might be missing and
    # are not smoothed
    # add the original values for those cells back
    all_cells <- colnames(expr_values)
    smoothed_cells <- colnames(spatial_smooth_matrix)
    missing_cells <- all_cells[!all_cells %in% smoothed_cells]

    if (length(missing_cells) > 0) {
        missing_matrix <- expr_values[, missing_cells]
        spatial_smooth_matrix <- cbind(spatial_smooth_matrix[
            rownames(expr_values),
        ], missing_matrix)
    }

    spatial_smooth_matrix <- spatial_smooth_matrix[
        rownames(expr_values), colnames(expr_values)
    ]

    # combine original and smoothed values according to smoothening b
    # create best guess for b if not given
    if (is.null(b)) {
        k <- stats::median(table(spatial_network$to))
        smooth_b <- 1 - 1 / k
    } else {
        smooth_b <- b
    }

    expr_b <- 1 - smooth_b
    spatsmooth_expr_values <- ((smooth_b * spatial_smooth_matrix) + (
        expr_b * expr_values))

    return(spatsmooth_expr_values)
}




#' @title Evaluate provided spatial locations
#' @name evaluate_provided_spatial_locations
#' @returns character
#' @keywords internal
evaluate_provided_spatial_locations <- function(spatial_locs) {
    if (!inherits(spatial_locs, "data.frame")) {
        stop("The spatial locations must be a data.frame(-like) object")
    }

    locs_names <- colnames(spatial_locs)
    required_cols <- c("sdimx", "sdimy", "cell_ID")
    missing_cols <- required_cols[!required_cols %in% locs_names]

    if (length(missing_cols) > 0) {
        stop("missing columns: ", list(missing_cols))
    } else {
        return(TRUE)
    }
}



#' @title do_spatial_grid_averaging
#' @name do_spatial_grid_averaging
#' @description smooth gene expression over a defined spatial grid
#' @returns matrix with smoothened gene expression values based on spatial grid
#' @keywords internal
do_spatial_grid_averaging <- function(
        expression_matrix,
        spatial_grid,
        spatial_locs,
        subset_feats = NULL,
        min_cells_per_grid = 4) {
    # matrix
    expr_values <- expression_matrix
    if (!is.null(subset_feats)) {
        expr_values <- expr_values[rownames(expr_values) %in% subset_feats, ]
    }

    # check spatial grid
    if (!inherits(x = spatial_grid, what = "spatialGridObj")) {
        stop("spatial_grid needs to be spatialGridObj")
    }

    # check spatial locations
    evaluate_provided_spatial_locations(spatial_locs = spatial_locs)

    # annoate spatial locations with spatial grid
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
    gr_loc <- NULL

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
    loc_av_expr_matrix <- as.matrix(loc_av_expr_matrix)

    return(loc_av_expr_matrix)
}



#' @title Detect spatially correlated features
#' @name detectSpatialCorFeats
#' @description Detect features that are spatially correlated. Functions for
#' starting from either a gobject (`detectSpatialCorFeats()`) or individual
#' pieces of data (`detectSpatialCorFeatsMatrix()`) are provided.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values gene expression values to use
#' @param expression_matrix provided expression matrix
#' @param spat_loc_name name for spatial locations
#' @param spatial_locs provided spatial locations
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network provided spatial network
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid provided spatial grid
#' @param method method to use for spatial averaging
#' @param subset_feats subset of features to use
#' @param network_smoothing smoothing factor between 0 and 1
#' (has automatic default, see details)
#' @param min_cells_per_grid minimum number of cells to consider a grid
#' @param cor_method correlation method
#' @returns returns a spatial correlation object: `spatCorObject`
#' @details
#' For `method = network`, it expects a fully connected spatial network.
#' You can make sure to create a
#' fully connected network by setting minimal_k > 0 in the
#' \code{\link{createSpatialNetwork}} function.
#'   1. **grid-averaging:** average gene expression values within a predefined
#'   spatial grid
#'   2. **network-averaging:** smoothens the gene expression matrix by
#'   averaging the expression within one cell by using the neighbours within
#'   the predefined spatial network. \eqn{b} is a smoothening factor passed by
#'   `network_smoothing` param that defaults to \eqn{1 - 1/k}, where \eqn{k}
#'   is the median number of k-neighbors in the selected spatial network.
#'   Setting \eqn{b = 0} means no smoothing and \eqn{b = 1} means no
#'   contribution from its own expression.
#'
#' The `spatCorObject` can be further explored with `showSpatialCorFeats()`
#' @seealso \code{\link{showSpatialCorFeats}}
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' # Perform with data in a gobject
#' detectSpatialCorFeats(g, method = "network")
#'
#' # This analysis can also be performed with data outside of the gobject
#' detectSpatialCorFeatsMatrix(
#'     expression_matrix = getExpression(
#'         g,
#'         output = "matrix"
#'     ),
#'     method = "network",
#'     spatial_network = getSpatialNetwork(g, output = "networkDT")
#' )
#'
NULL



#' @rdname detectSpatialCorFeats
#' @export
detectSpatialCorFeats <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    spat_loc_name = "raw",
    method = c("grid", "network"),
    expression_values = c("normalized", "scaled", "custom"),
    subset_feats = NULL,
    spatial_network_name = "Delaunay_network",
    network_smoothing = NULL,
    spatial_grid_name = "spatial_grid",
    min_cells_per_grid = 4,
    cor_method = c("pearson", "kendall", "spearman")) {
    # set default spat_unit and feat_type
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    ## correlation method to be used
    cor_method <- match.arg(
        cor_method,
        choices = c("pearson", "kendall", "spearman")
    )

    ## method to be used
    method <- match.arg(method, choices = c("grid", "network"))

    # get expression matrix
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

    if (!is.null(subset_feats)) {
        expr_values <- expr_values[rownames(expr_values) %in% subset_feats, ]
    }



    # get spatial locations
    spatial_locs <- getSpatialLocations(
        gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table",
        copy_obj = TRUE
    )

    ## spatial averaging or smoothing
    if (method == "grid") {
        # get spatial grid
        spatial_grid <- getSpatialGrid(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = spatial_grid_name,
            return_grid_Obj = FALSE
        )

        loc_av_expr_matrix <- do_spatial_grid_averaging(
            expression_matrix = as.matrix(expr_values),
            spatial_grid = spatial_grid,
            spatial_locs = spatial_locs,
            subset_feats = subset_feats,
            min_cells_per_grid = min_cells_per_grid
        )

        # data.table variables
        feat_ID <- variable <- NULL

        cor_spat_matrix <- cor_flex(t_flex(as.matrix(
            loc_av_expr_matrix
        )), method = cor_method)
        cor_spat_matrixDT <- data.table::as.data.table(cor_spat_matrix)
        cor_spat_matrixDT[, feat_ID := rownames(cor_spat_matrix)]
        cor_spat_DT <- data.table::melt.data.table(
            data = cor_spat_matrixDT,
            id.vars = "feat_ID", value.name = "spat_cor"
        )
    }

    if (method == "network") {
        # get spatial network
        spatial_network <- getSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            name = spatial_network_name,
            output = "networkDT"
        )

        knn_av_expr_matrix <- do_spatial_knn_smoothing(
            expression_matrix = as.matrix(expr_values),
            spatial_network = spatial_network,
            subset_feats = subset_feats,
            b = network_smoothing
        )




        cor_spat_matrix <- cor_flex(t_flex(as.matrix(
            knn_av_expr_matrix
        )), method = cor_method)
        cor_spat_matrixDT <- data.table::as.data.table(cor_spat_matrix)
        cor_spat_matrixDT[, feat_ID := rownames(cor_spat_matrix)]
        cor_spat_DT <- data.table::melt.data.table(
            data = cor_spat_matrixDT,
            id.vars = "feat_ID", value.name = "spat_cor"
        )
    }



    # data.table variables
    cordiff <- spat_cor <- expr_cor <- spatrank <- exprrank <- rankdiff <- NULL

    ## 2. perform expression correlation at single-cell level without
    ## spatial information
    cor_matrix <- cor_flex(t_flex(expr_values), method = cor_method)
    cor_matrixDT <- data.table::as.data.table(cor_matrix)
    cor_matrixDT[, feat_ID := rownames(cor_matrix)]
    cor_DT <- data.table::melt.data.table(
        data = cor_matrixDT,
        id.vars = "feat_ID", value.name = "expr_cor"
    )

    ## 3. merge spatial and expression correlation
    data.table::setorder(cor_spat_DT, feat_ID, variable)
    data.table::setorder(cor_DT, feat_ID, variable)
    doubleDT <- cbind(cor_spat_DT, expr_cor = cor_DT[["expr_cor"]])

    # difference in correlation scores
    doubleDT[, cordiff := spat_cor - expr_cor]

    # difference in rank scores
    doubleDT[, spatrank := frank(
        -spat_cor,
        ties.method = "first"
    ), by = feat_ID]
    doubleDT[, exprrank := frank(
        -expr_cor,
        ties.method = "first"
    ), by = feat_ID]
    doubleDT[, rankdiff := spatrank - exprrank]

    # sort data
    data.table::setorder(doubleDT, feat_ID, -spat_cor)

    spatCorObject <- list(
        cor_DT = doubleDT,
        feat_order = rownames(cor_spat_matrix),
        cor_hclust = list(),
        cor_clusters = list()
    )

    class(spatCorObject) <- append("spatCorObject", class(spatCorObject))

    return(spatCorObject)
}


#' @rdname detectSpatialCorFeats
#' @export
detectSpatialCorFeatsMatrix <- function(
        expression_matrix,
        method = c("grid", "network"),
        spatial_network,
        spatial_grid,
        spatial_locs,
        subset_feats = NULL,
        network_smoothing = NULL,
        min_cells_per_grid = 4,
        cor_method = c("pearson", "kendall", "spearman")) {
    ## correlation method to be used
    cor_method <- match.arg(
        cor_method,
        choices = c("pearson", "kendall", "spearman")
    )

    ## method to be used
    method <- match.arg(method, choices = c("grid", "network"))

    ## spatial averaging or smoothing
    if (method == "grid") {
        loc_av_expr_matrix <- do_spatial_grid_averaging(
            expression_matrix = as.matrix(expression_matrix),
            spatial_grid = spatial_grid,
            spatial_locs = spatial_locs,
            subset_feats = subset_feats,
            min_cells_per_grid = min_cells_per_grid
        )

        # data.table variables
        feat_ID <- variable <- NULL

        cor_spat_matrix <- cor_flex(t_flex(
            as.matrix(loc_av_expr_matrix)
        ), method = cor_method)
        cor_spat_matrixDT <- data.table::as.data.table(cor_spat_matrix)
        cor_spat_matrixDT[, feat_ID := rownames(cor_spat_matrix)]
        cor_spat_DT <- data.table::melt.data.table(
            data = cor_spat_matrixDT,
            id.vars = "feat_ID", value.name = "spat_cor"
        )
    }

    if (method == "network") {
        knn_av_expr_matrix <- do_spatial_knn_smoothing(
            expression_matrix = as.matrix(expression_matrix),
            spatial_network = spatial_network,
            subset_feats = subset_feats,
            b = network_smoothing
        )



        cor_spat_matrix <- cor_flex(t_flex(as.matrix(
            knn_av_expr_matrix
        )), method = cor_method)
        cor_spat_matrixDT <- data.table::as.data.table(cor_spat_matrix)
        cor_spat_matrixDT[, feat_ID := rownames(cor_spat_matrix)]
        cor_spat_DT <- data.table::melt.data.table(
            data = cor_spat_matrixDT,
            id.vars = "feat_ID", value.name = "spat_cor"
        )
    }



    # data.table variables
    cordiff <- spat_cor <- expr_cor <- spatrank <- exprrank <- rankdiff <- NULL

    ## 2. perform expression correlation at single-cell level without
    ## spatial information

    # matrix
    expr_values <- expression_matrix
    if (!is.null(subset_feats)) {
        expr_values <- expr_values[rownames(expr_values) %in% subset_feats, ]
    }

    cor_matrix <- cor_flex(t_flex(expr_values), method = cor_method)
    cor_matrixDT <- data.table::as.data.table(cor_matrix)
    cor_matrixDT[, feat_ID := rownames(cor_matrix)]
    cor_DT <- data.table::melt.data.table(
        data = cor_matrixDT,
        id.vars = "feat_ID", value.name = "expr_cor"
    )

    ## 3. merge spatial and expression correlation
    data.table::setorder(cor_spat_DT, feat_ID, variable)
    data.table::setorder(cor_DT, feat_ID, variable)
    doubleDT <- cbind(cor_spat_DT, expr_cor = cor_DT[["expr_cor"]])

    # difference in correlation scores
    doubleDT[, cordiff := spat_cor - expr_cor]

    # difference in rank scores
    doubleDT[, spatrank := data.table::frank(
        -spat_cor,
        ties.method = "first"
    ), by = feat_ID]
    doubleDT[, exprrank := data.table::frank(
        -expr_cor,
        ties.method = "first"
    ), by = feat_ID]
    doubleDT[, rankdiff := spatrank - exprrank]

    # sort data
    data.table::setorder(doubleDT, feat_ID, -spat_cor)

    spatCorObject <- list(
        cor_DT = doubleDT,
        feat_order = rownames(cor_spat_matrix),
        cor_hclust = list(),
        cor_clusters = list()
    )

    class(spatCorObject) <- append(class(spatCorObject), "spatCorObject")

    return(spatCorObject)
}






#' @title showSpatialCorFeats
#' @name showSpatialCorFeats
#' @description Shows and filters spatially correlated features
#' @param spatCorObject spatial correlation object
#' @param use_clus_name cluster information to show
#' @param selected_clusters subset of clusters to show
#' @param feats subset of features to show
#' @param min_spat_cor filter on minimum spatial correlation
#' @param min_expr_cor filter on minimum single-cell expression correlation
#' @param min_cor_diff filter on minimum correlation difference
#' (spatial vs expression)
#' @param min_rank_diff filter on minimum correlation rank difference
#' (spatial vs expression)
#' @param show_top_feats show top features per gene
#' @returns data.table with filtered information
#' @export
showSpatialCorFeats <- function(
        spatCorObject,
        use_clus_name = NULL,
        selected_clusters = NULL,
        feats = NULL,
        min_spat_cor = 0.5,
        min_expr_cor = NULL,
        min_cor_diff = NULL,
        min_rank_diff = NULL,
        show_top_feats = NULL) {
    # data.table variables
    clus <- feat_ID <- spat_cor <- cor_diff <- rankdiff <- NULL

    if (!"spatCorObject" %in% class(spatCorObject)) {
        stop("spatCorObject needs to be the output from
            detectSpatialCorfeats()")
    }

    filter_DT <- data.table::copy(spatCorObject[["cor_DT"]])

    if (!is.null(use_clus_name)) {
        clusters_part <- spatCorObject[["cor_clusters"]][[use_clus_name]]

        # combine spatial correlation info and clusters
        clusters <- clusters_part
        names_clusters <- names(clusters_part)
        clusters_DT <- data.table::data.table(
            "feat_ID" = names_clusters, "clus" = clusters
        )
        filter_DT <- data.table::merge.data.table(
            filter_DT, clusters_DT,
            by = "feat_ID"
        )
    }

    ## 0. subset clusters
    if (!is.null(selected_clusters)) {
        filter_DT <- filter_DT[clus %in% selected_clusters]
    }


    ## 1. subset feats
    if (!is.null(feats)) {
        filter_DT <- filter_DT[feat_ID %in% feats]
    }

    ## 2. select spatial correlation
    if (!is.null(min_spat_cor)) {
        filter_DT <- filter_DT[spat_cor >= min_spat_cor]
    }

    ## 3. minimum expression correlation
    if (!is.null(min_expr_cor)) {
        filter_DT <- filter_DT[spat_cor >= min_expr_cor]
    }

    ## 4. minimum correlation difference
    if (!is.null(min_cor_diff)) {
        filter_DT <- filter_DT[cor_diff >= min_cor_diff]
    }

    ## 5. minimum correlation difference
    if (!is.null(min_rank_diff)) {
        filter_DT <- filter_DT[rankdiff >= min_rank_diff]
    }

    ## 6. show only top feats
    if (!is.null(show_top_feats)) {
        filter_DT <- filter_DT[, head(.SD, show_top_feats), by = feat_ID]
    }

    return(filter_DT)
}



#' @title showSpatialCorGenes
#' @name showSpatialCorGenes
#' @description Shows and filters spatially correlated genes
#' @param spatCorObject spatial correlation object
#' @param use_clus_name cluster information to show
#' @param selected_clusters subset of clusters to show
#' @param genes subset of genes to show
#' @param min_spat_cor filter on minimum spatial correlation
#' @param min_expr_cor filter on minimum single-cell expression correlation
#' @param min_cor_diff filter on minimum correlation difference
#' (spatial vs expression)
#' @param min_rank_diff filter on minimum correlation rank difference
#' (spatial vs expression)
#' @param show_top_genes show top genes per gene
#' @returns data.table with filtered information
#' @export
showSpatialCorGenes <- function(
        spatCorObject,
        use_clus_name = NULL,
        selected_clusters = NULL,
        genes = NULL,
        min_spat_cor = 0.5,
        min_expr_cor = NULL,
        min_cor_diff = NULL,
        min_rank_diff = NULL,
        show_top_genes = NULL) {
    warning("Deprecated and replaced by showSpatialCorFeats")

    showSpatialCorFeats(
        spatCorObject = spatCorObject,
        use_clus_name = use_clus_name,
        selected_clusters = selected_clusters,
        feats = genes,
        min_spat_cor = min_spat_cor,
        min_expr_cor = min_expr_cor,
        min_cor_diff = min_cor_diff,
        min_rank_diff = min_rank_diff,
        show_top_feats = show_top_genes
    )
}







#' @title clusterSpatialCorFeats
#' @name clusterSpatialCorFeats
#' @description Cluster based on spatially correlated features
#' @param spatCorObject spatial correlation object
#' @param name name for spatial clustering results
#' @param hclust_method method for hierarchical clustering
#' @param k number of clusters to extract
#' @param return_obj return spatial correlation object (spatCorObject)
#' @returns spatCorObject or cluster results
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' clusterSpatialCorFeats(spatCorObject = detectSpatialCorFeats(
#'     g,
#'     method = "network"
#' ))
#' @export
clusterSpatialCorFeats <- function(
        spatCorObject,
        name = "spat_clus",
        hclust_method = "ward.D",
        k = 10,
        return_obj = TRUE) {
    # check input
    if (!"spatCorObject" %in% class(spatCorObject)) {
        stop("spatCorObject needs to be the output from
            detectSpatialCorfeats()")
    }

    # create correlation matrix
    cor_DT <- spatCorObject[["cor_DT"]]
    cor_DT_dc <- data.table::dcast.data.table(
        cor_DT,
        formula = feat_ID ~ variable, value.var = "spat_cor"
    )
    cor_matrix <- dt_to_matrix(cor_DT_dc)

    # re-ordering matrix
    my_feat_order <- spatCorObject[["feat_order"]]
    cor_matrix <- cor_matrix[my_feat_order, my_feat_order]

    # cluster
    cor_dist <- stats::as.dist(1 - cor_matrix)
    cor_h <- stats::hclust(d = cor_dist, method = hclust_method)
    cor_clus <- stats::cutree(cor_h, k = k)

    if (return_obj == TRUE) {
        spatCorObject[["cor_hclust"]][[name]] <- cor_h
        spatCorObject[["cor_clusters"]][[name]] <- cor_clus
        spatCorObject[["cor_coexpr_groups"]][[name]] <- NA

        return(spatCorObject)
    } else {
        return(list("hclust" = cor_h, "clusters" = cor_clus))
    }
}




#' @title clusterSpatialCorGenes
#' @name clusterSpatialCorGenes
#' @description Cluster based on spatially correlated genes
#' @param spatCorObject spatial correlation object
#' @param name name for spatial clustering results
#' @param hclust_method method for hierarchical clustering
#' @param k number of clusters to extract
#' @param return_obj return spatial correlation object (spatCorObject)
#' @returns spatCorObject or cluster results
#' @export
clusterSpatialCorGenes <- function(
        spatCorObject,
        name = "spat_clus",
        hclust_method = "ward.D",
        k = 10,
        return_obj = TRUE) {
    warning("Deprecated and replaced by clusterSpatialCorFeats")

    clusterSpatialCorFeats(
        spatCorObject = spatCorObject,
        name = name,
        hclust_method = hclust_method,
        k = k,
        return_obj = return_obj
    )
}





#' @title heatmSpatialCorFeats
#' @name heatmSpatialCorFeats
#' @description Create heatmap of spatially correlated features
#' @param gobject giotto object
#' @param spatCorObject spatial correlation object
#' @param use_clus_name name of clusters to visualize
#' (from clusterSpatialCorFeats())
#' @param show_cluster_annot show cluster annotation on top of heatmap
#' @param show_row_dend show row dendrogram
#' @param show_column_dend show column dendrogram
#' @param show_row_names show row names
#' @param show_column_names show column names
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see
#' \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @param \dots additional parameters to the
#' \code{\link[ComplexHeatmap]{Heatmap}} function from ComplexHeatmap
#' @returns Heatmap generated by ComplexHeatmap
#' @export
heatmSpatialCorFeats <- function(
        gobject,
        spatCorObject,
        use_clus_name = NULL,
        show_cluster_annot = TRUE,
        show_row_dend = TRUE,
        show_column_dend = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "heatmSpatialCorFeats",
        ...) {
    ## check input
    if (!"spatCorObject" %in% class(spatCorObject)) {
        stop("spatCorObject needs to be the output from
            detectSpatialCorFeats()")
    }

    ## package check for ComplexHeatmap
    package_check(pkg_name = "ComplexHeatmap", repository = "CRAN")

    ## create correlation matrix
    cor_DT <- spatCorObject[["cor_DT"]]
    cor_DT_dc <- data.table::dcast.data.table(
        cor_DT,
        formula = feat_ID ~ variable, value.var = "spat_cor"
    )
    cor_matrix <- dt_to_matrix(cor_DT_dc)

    # re-ordering matrix
    my_feat_order <- spatCorObject[["feat_order"]]
    cor_matrix <- cor_matrix[my_feat_order, my_feat_order]


    ## fix row and column names
    cor_matrix <- cor_matrix[rownames(cor_matrix), rownames(cor_matrix)]

    ## default top annotation
    ha <- NULL

    if (!is.null(use_clus_name)) {
        hclust_part <- spatCorObject[["cor_hclust"]][[use_clus_name]]

        if (is.null(hclust_part)) {
            cat(use_clus_name, " does not exist, make one with spatCorCluster")
            hclust_part <- TRUE
        } else {
            clusters_part <- spatCorObject[["cor_clusters"]][[use_clus_name]]

            if (show_cluster_annot) {
                uniq_clusters <- unique(clusters_part)

                # color vector
                mycolors <- getDistinctColors(length(uniq_clusters))
                names(mycolors) <- uniq_clusters
                ha <- ComplexHeatmap::HeatmapAnnotation(
                    bar = as.vector(clusters_part),
                    col = list(bar = mycolors),
                    annotation_legend_param = list(title = NULL)
                )
            }
        }
    } else {
        hclust_part <- TRUE
    }


    ## create heatmap
    heatm <- ComplexHeatmap::Heatmap(
        matrix = as.matrix(cor_matrix),
        cluster_rows = hclust_part,
        cluster_columns = hclust_part,
        show_row_dend = show_row_dend,
        show_column_dend = show_column_dend,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        top_annotation = ha, ...
    )

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = heatm,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}







#' @title rankSpatialCorGroups
#' @name rankSpatialCorGroups
#' @description Rank spatial correlated clusters according to correlation
#' structure
#' @param gobject giotto object
#' @param spatCorObject spatial correlation object
#' @param use_clus_name name of clusters to visualize
#' (from `clusterSpatialCorFeats()`)
#' @param show_plot logical. show plot
#' @param return_plot logical. return ggplot object
#' @param save_plot logical. directly save the plot
#' @param save_param list of saving parameters, see
#' \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns data.table with positive (within group) and negative
#' (outside group) scores
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' spatCorObject <- detectSpatialCorFeats(g, method = "network")
#' clusters <- clusterSpatialCorFeats(spatCorObject = spatCorObject)
#'
#' rankSpatialCorGroups(
#'     gobject = g, spatCorObject = clusters,
#'     use_clus_name = "spat_clus"
#' )
#' @md
#' @export
rankSpatialCorGroups <- function(
        gobject,
        spatCorObject,
        use_clus_name = NULL,
        show_plot = NULL,
        return_plot = FALSE,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "rankSpatialCorGroups") {
    ## check input
    if (!"spatCorObject" %in% class(spatCorObject)) {
        stop("spatCorObject needs to be the output from
            detectSpatialCorFeats()")
    }

    ## check if cluster exist
    if (is.null(use_clus_name)) {
        stop("use_clus_name does not exist")
    }
    clusters_part <- spatCorObject[["cor_clusters"]][[use_clus_name]]

    if (is.null(clusters_part)) {
        stop("use_clus_name does not exist")
    }

    ## create correlation matrix
    cor_DT <- spatCorObject[["cor_DT"]]
    cor_DT_dc <- data.table::dcast.data.table(
        cor_DT,
        formula = feat_ID ~ variable, value.var = "spat_cor"
    )
    cor_matrix <- dt_to_matrix(cor_DT_dc)

    # re-ordering matrix
    my_feat_order <- spatCorObject[["feat_order"]]
    cor_matrix <- cor_matrix[my_feat_order, my_feat_order]



    res_cor_list <- list()
    res_neg_cor_list <- list()
    nr_feats_list <- list()

    for (id in seq_along(unique(clusters_part))) {
        clus_id <- unique(clusters_part)[id]
        selected_feats <- names(clusters_part[clusters_part == clus_id])
        nr_feats_list[[id]] <- length(selected_feats)

        sub_cor_matrix <- cor_matrix[
            rownames(cor_matrix) %in% selected_feats,
            colnames(cor_matrix) %in% selected_feats
        ]
        mean_score <- mean_flex(sub_cor_matrix)
        res_cor_list[[id]] <- mean_score

        sub_neg_cor_matrix <- cor_matrix[
            rownames(cor_matrix) %in% selected_feats,
            !colnames(cor_matrix) %in% selected_feats
        ]
        mean_neg_score <- mean_flex(sub_neg_cor_matrix)
        res_neg_cor_list[[id]] <- mean_neg_score
    }


    # data.table variables
    cor_neg_adj <- cor_neg_score <- adj_cor_score <- cor_score <-
        clusters <- nr_feats <- NULL

    res_cor_DT <- data.table::data.table(
        "clusters" = unique(clusters_part),
        cor_score = unlist(res_cor_list),
        cor_neg_score = unlist(res_neg_cor_list),
        nr_feats = unlist(nr_feats_list)
    )

    res_cor_DT[, cor_neg_adj := 1 - (cor_neg_score - min(cor_neg_score))]
    res_cor_DT[, adj_cor_score := cor_neg_adj * cor_score]
    data.table::setorder(res_cor_DT, -adj_cor_score)
    res_cor_DT[, clusters := factor(x = clusters, levels = rev(clusters))]

    pl <- gg_simple_scatter(
        data = res_cor_DT,
        x = "clusters",
        y = "adj_cor_score",
        size = "nr_feats",
        xlab = "cluster",
        ylab = "pos r x (1 - (neg_r - min(neg_r)))"
    )

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = res_cor_DT
    ))
}




#' @title Get balanced spatial coexpression features
#' @name getBalancedSpatCoexpressionFeats
#' @description Extract features from spatial co-expression modules in a
#' balanced manner
#' @param spatCorObject spatial correlation object
#' @param maximum maximum number of genes to get from each spatial
#' co-expression module
#' @param rank ranking method (see details)
#' @param informed_ranking vector of ranked features
#' @param seed seed
#' @param verbose verbosity
#' @returns balanced vector with features for each co-expression module
#' @details There are 3 different ways of selecting features from the spatial
#' co-expression modules
#'   1. **weighted:** Features are ranked based on summarized pairwise
#'   co-expression scores
#'   2. **random:** A random selection of features, set seed for
#'   reproducibility
#'   3. **informed:** Features are selected based on prior information/ranking
#'
#' @md
#' @export
getBalancedSpatCoexpressionFeats <- function(
        spatCorObject,
        maximum = 50,
        rank = c("weighted", "random", "informed"),
        informed_ranking = NULL,
        seed = NA,
        verbose = TRUE) {
    # data.table vars
    feat_ID <- variable <- combo <- spat_cor <- rnk <- feat_id <- V1 <- NULL

    rank <- match.arg(rank, choices = c("weighted", "random", "informed"))

    clusters <- spatCorObject$cor_clusters$spat_netw_clus

    # rank = random
    if (rank == "random") {
        if (!is.na(seed) & is.numeric(seed)) {
            on.exit(random_seed(), add = TRUE)
            set.seed(seed)
            wrap_msg("Seed has been set for random")
        } else {
            wrap_msg("Random is selected, but no seed has been set \n
            Results might not be fully reproducible")
        }

        result_list <- list()
        for (clus in seq_along(unique(clusters))) {
            selected_cluster_features <- names(clusters[clusters == clus])

            feat_length <- length(selected_cluster_features)
            if (feat_length < maximum) {
                maximum_to_use <- feat_length
                wrap_msg(
                    "There are only ", feat_length, " features for cluster ",
                    clus, "\n",
                    "Maximum will be set to ", feat_length
                )
            } else {
                maximum_to_use <- maximum
            }

            selected_feats <- sample(
                x = selected_cluster_features,
                size = maximum_to_use,
                replace = FALSE
            )
            clus_id <- rep(clus, length(selected_feats))
            names(clus_id) <- selected_feats
            result_list[[clus]] <- clus_id
        }

        final_res <- do.call("c", result_list)
    }


    # rank = random
    if (rank == "weighted") {
        cor_data <- spatCorObject$cor_DT

        result_list <- list()
        for (clus in seq_along(unique(clusters))) {
            if (verbose) print(clus)

            # get all pairwise spatial feature correlations and rank them
            selected_cluster_features <- names(clusters[clusters == clus])
            subset_cor_data <- cor_data[
                feat_ID %in% selected_cluster_features &
                    variable %in% selected_cluster_features
            ]
            subset_cor_data <- subset_cor_data[feat_ID != variable]
            subset_cor_data <- dt_sort_combine_two_columns(
                DT = subset_cor_data,
                column1 = "feat_ID",
                column2 = "variable", myname = "combo"
            )
            subset_cor_data <- subset_cor_data[duplicated(combo)]
            data.table::setorder(subset_cor_data, -spat_cor)

            # create a ranked data.table
            rnk1DT <- data.table::data.table(
                feat_id = subset_cor_data$feat_ID,
                rnk = seq_along(subset_cor_data$feat_ID)
            )
            rnk2DT <- data.table::data.table(
                feat_id = subset_cor_data$variable,
                rnk = seq_along(subset_cor_data$variable)
            )
            rnkDT <- data.table::rbindlist(list(rnk1DT, rnk2DT))
            data.table::setorder(rnkDT, rnk)

            # summarize rank (weights)
            rnkcombined <- rnkDT[, sum(rnk), by = feat_id]
            data.table::setorder(rnkcombined, V1)

            feat_length <- nrow(rnkcombined)
            if (feat_length < maximum) {
                maximum_to_use <- feat_length
                wrap_msg(
                    "There are only ", feat_length, " features for cluster ",
                    clus, "\n",
                    "Maximum will be set to ", feat_length
                )
            } else {
                maximum_to_use <- maximum
            }

            selected_feats <- rnkcombined[seq_len(maximum_to_use)][["feat_id"]]

            clus_id <- rep(clus, length(selected_feats))
            names(clus_id) <- selected_feats
            result_list[[clus]] <- clus_id
        }


        final_res <- do.call("c", result_list)
    }


    # rank = random
    if (rank == "informed") {
        if (is.null(informed_ranking)) {
            stop("Informed has been selected, but no informed ranking
                vector has been provided")
        }

        # informed_ranking vector should be a ranked gene list
        informed_ranking_numerical <- seq_along(informed_ranking)
        names(informed_ranking_numerical) <- informed_ranking

        result_list <- list()
        for (clus in seq_along(unique(clusters))) {
            selected_cluster_features <- names(clusters[clusters == clus])

            feat_length <- length(selected_cluster_features)
            if (feat_length < maximum) {
                maximum_to_use <- feat_length
                wrap_msg(
                    "There are only ", feat_length, " features for cluster ",
                    clus, "\n",
                    "Maximum will be set to ", feat_length
                )
            } else {
                maximum_to_use <- maximum
            }


            informed_subset <- informed_ranking_numerical[
                names(informed_ranking_numerical) %in%
                    selected_cluster_features
            ]
            informed_subset <- sort(informed_subset)

            feat_length <- length(informed_subset)
            if (feat_length < maximum) {
                maximum_to_use <- feat_length
                wrap_msg(
                    "There are only ", feat_length, " features for cluster ",
                    clus, "\n",
                    "Maximum will be set to ", feat_length
                )
            } else {
                maximum_to_use <- maximum
            }

            selected_feats <- names(informed_subset[seq_len(maximum_to_use)])

            clus_id <- rep(clus, length(selected_feats))
            names(clus_id) <- selected_feats
            result_list[[clus]] <- clus_id
        }

        final_res <- do.call("c", result_list)
    }

    return(final_res)
}





# ** ####
## Simulate single-gene spatial patterns ####
## --------------------------------------- ##



#' @title simulateOneGenePatternGiottoObject
#' @name simulateOneGenePatternGiottoObject
#' @description Create a simulated spatial pattern for one selected gnee
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_name selected gene
#' @param spatial_prob probability for a high expressing gene value to be
#' part of the spatial pattern
#' @param gradient_direction direction of gradient
#' @param show_pattern show the discrete spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param normalization_params additional parameters for (re-)normalizing
#' @returns Reprocessed Giotto object for which one gene has a forced
#' spatial pattern
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' simulateOneGenePatternGiottoObject(
#'     gobject = g,
#'     pattern_cell_ids = c(
#'         "AAAGGGATGTAGCAAG-1", "TCAAACAACCGCGTCG-1",
#'         "ACGATCATACATAGAG-1", "TATGCTCCCTACTTAC-1"
#'     ),
#'     gene_name = "Gna12"
#' )
#' @export
simulateOneGenePatternGiottoObject <- function(
        gobject,
        pattern_name = "pattern",
        pattern_cell_ids = NULL,
        gene_name = NULL,
        spatial_prob = 0.95,
        gradient_direction = NULL,
        show_pattern = TRUE,
        pattern_colors = c("in" = "green", "out" = "red"),
        normalization_params = list()) {
    # data.table variables
    cell_ID <- sdimx_y <- sdimx <- sdimy <- NULL

    if (is.null(pattern_cell_ids)) {
        stop("pattern_cell_ids can not be NULL")
    }

    ## create and add annotation for pattern
    cell_meta <- pDataDT(gobject)
    cell_meta[, (pattern_name) := ifelse(
        cell_ID %in% pattern_cell_ids, "in", "out"
    )]

    newgobject <- addCellMetadata(
        gobject,
        new_metadata = cell_meta[, c("cell_ID", pattern_name), with = FALSE],
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )

    # show pattern
    if (show_pattern == TRUE) {
        spatPlot2D(
            gobject = newgobject,
            save_plot = FALSE,
            cell_color_code = pattern_colors,
            point_size = 2,
            cell_color = pattern_name
        )
    }


    ## merge cell metadata and cell coordinate data
    cell_meta <- pDataDT(newgobject)
    cell_coord <- getSpatialLocations(newgobject,
        name = "raw",
        output = "data.table",
        copy_obj = TRUE
    )
    cell_meta <- data.table::merge.data.table(
        cell_meta, cell_coord,
        by = "cell_ID"
    )

    ## get number of cells within pattern
    cell_number <- nrow(cell_meta[get(pattern_name) == "in"])


    ## normalized expression
    # expr_data <- newgobject@norm_expr
    expr_data <- getExpression(
        gobject = newgobject,
        values = "normalized",
        output = "matrix"
    )
    result_list <- list()

    ## raw expression
    # raw_expr_data <- newgobject@raw_exprs
    raw_expr_data <- getExpression(
        gobject = newgobject,
        values = "raw",
        output = "matrix"
    )
    raw_result_list <- list()


    ## create the spatial expression pattern for the specified gene
    # 1. rank all gene values from the cells from high to low
    # 2. move the highest expressing values to the spatial pattern using a
    # probability
    #     - 0.5 is the control = random
    #     - 1 is perfection: all the highest values go to the pattern
    #     - 0.5 to 1 is decreasing noise levels

    if (is.null(gene_name)) stop("a gene name needs to be provided")



    # rank genes
    gene_vector <- expr_data[rownames(expr_data) == gene_name, ]
    sort_expr_gene <- sort(gene_vector, decreasing = TRUE)

    # number of cells in and out the pattern
    total_cell_number <- length(sort_expr_gene)
    remaining_cell_number <- total_cell_number - cell_number

    # calculate outside probability
    outside_prob <- 1 - spatial_prob
    prob_vector <- c(
        rep(spatial_prob, cell_number),
        rep(outside_prob, remaining_cell_number)
    )

    # first get the 'in' pattern sample values randomly
    sample_values <- sample(
        sort_expr_gene,
        replace = FALSE, size = cell_number, prob = prob_vector
    )

    # then take the remaining 'out' pattern values randomly
    remain_values <- sort_expr_gene[
        !names(sort_expr_gene) %in% names(sample_values)
    ]
    remain_values <- sample(remain_values, size = length(remain_values))



    ## A. within pattern ##
    # ------------------- #
    in_cell_meta <- cell_meta[get(pattern_name) == "in"]

    # if gradient is wanted
    # does not work with 0.5!! is not random!!
    if (!is.null(gradient_direction)) {
        # sort in_ids according to x, y or  xy coordinates to create gradient
        in_cell_meta[, sdimx_y := abs(sdimx) + abs(sdimy)]
        # order according to gradient direction
        in_cell_meta <- in_cell_meta[order(get(gradient_direction))]
    }
    in_ids <- in_cell_meta$cell_ID

    # preparation for raw matrix
    sample_values_id_vector <- names(sample_values)
    names(sample_values_id_vector) <- in_ids


    ## B. outside pattern ##
    # -------------------- #
    out_ids <- cell_meta[get(pattern_name) == "out"]$cell_ID

    # preparation for raw matrix
    remain_values_id_vector <- names(remain_values)
    names(remain_values_id_vector) <- out_ids




    ## raw matrix
    # swap the cell ids #
    raw_gene_vector <- raw_expr_data[rownames(raw_expr_data) == gene_name, ]

    raw_new_sample_vector <- raw_gene_vector[sample_values_id_vector]
    names(raw_new_sample_vector) <- names(sample_values_id_vector)

    raw_new_remain_vector <- raw_gene_vector[remain_values_id_vector]
    names(raw_new_remain_vector) <- names(remain_values_id_vector)

    new_sim_raw_values <- c(raw_new_sample_vector, raw_new_remain_vector)
    new_sim_raw_values <- new_sim_raw_values[names(raw_gene_vector)]

    # change the original matrices
    raw_expr_data[rownames(raw_expr_data) == gene_name, ] <- new_sim_raw_values
    # newgobject@raw_exprs <- raw_expr_data
    newgobject <- setExpression(
        gobject = newgobject,
        x = createExprObj(
            expression_data = raw_expr_data,
            name = "raw"
        ),
        name = "raw",
        provenance = prov(getCellMetadata(newgobject))
    )

    # recalculate normalized values
    newgobject <- do.call(
        "normalizeGiotto",
        args = c(gobject = newgobject, normalization_params)
    )

    newgobject <- addStatistics(gobject = newgobject)

    return(newgobject)
}






#' @title run_spatial_sim_tests_one_rep
#' @name run_spatial_sim_tests_one_rep
#' @description runs all spatial tests for 1 probability and 1 rep
#' @returns data.table
#' @keywords internal
run_spatial_sim_tests_one_rep <- function(
        gobject,
        pattern_name = "pattern",
        pattern_cell_ids = NULL,
        gene_name = NULL,
        spatial_prob = 0.95,
        show_pattern = FALSE,
        spatial_network_name = "kNN_network",
        spat_methods = c(
            "binSpect_single", "binSpect_multi", "spatialDE",
            "spark", "silhouetteRank"
        ),
        spat_methods_params = list(NA, NA, NA, NA, NA),
        spat_methods_names = c(
            "binSpect_single", "binSpect_multi", "spatialDE",
            "spark", "silhouetteRank"
        ),
        save_plot = FALSE,
        save_raw = FALSE,
        save_norm = FALSE,
        save_dir = "~",
        save_name = "plot",
        run_simulations = TRUE,
        ...) {
    # data.table variables
    genes <- prob <- time <- adj.p.value <- method <- p.val <- sd <-
        qval <- pval <- g <- adjusted_pvalue <- feats <- NULL

    ## test if spat_methods, params and names have the same length
    if (length(spat_methods) != length(spat_methods_params)) {
        stop("number of spatial detection methods to test need to be equal
            to number of spatial methods parameters")
    }
    if (length(spat_methods) != length(spat_methods_names)) {
        stop("number of spatial detection methods to test need to be equal
            to number of spatial methods names")
    }


    ## simulate pattern ##
    simulate_patch <- simulateOneGenePatternGiottoObject(gobject,
        pattern_name = pattern_name,
        pattern_cell_ids = pattern_cell_ids,
        gene_name = gene_name,
        spatial_prob = spatial_prob,
        gradient_direction = NULL,
        show_pattern = show_pattern,
        ...
    )

    # save plot
    if (save_plot == TRUE) {
        spatFeatPlot2D(simulate_patch,
            expression_values = "normalized",
            feats = gene_name,
            point_shape = "border",
            point_border_stroke = 0.1,
            point_size = 2.5,
            cow_n_col = 1, show_plot = FALSE,
            save_plot = TRUE,
            save_param = list(
                save_dir = save_dir, save_folder = pattern_name,
                save_name = save_name,
                base_width = 9, base_height = 7, units = "cm"
            )
        )
    }

    # save raw data
    if (save_raw == TRUE) {
        folder_path <- paste0(save_dir, "/", pattern_name)
        if (!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

        write.table(
            x = as.matrix(getExpression(
                gobject = simulate_patch, values = "raw", output = "matrix"
            )),
            file = paste0(
                save_dir, "/", pattern_name, "/", save_name, "_raw_data.txt"
            ),
            sep = "\t"
        )
    }

    # save normalized data
    if (save_norm == TRUE) {
        folder_path <- paste0(save_dir, "/", pattern_name)
        if (!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

        write.table(
            x = as.matrix(getExpression(
                gobject = simulate_patch,
                values = "normalized",
                output = "matrix"
            )),
            file = paste0(
                save_dir, "/", pattern_name, "/", save_name, "_norm_data.txt"
            ),
            sep = "\t"
        )
    }



    ## do simulations ##
    if (run_simulations == TRUE) {
        result_list <- list()
        for (test in seq_along(spat_methods)) {
            # method
            selected_method <- spat_methods[test]
            if (!selected_method %in%
                c(
                    "binSpect_single", "binSpect_multi", "spatialDE", "spark",
                    "silhouetteRank"
                )) {
                stop(selected_method, " is not a know spatial method")
            }

            # params
            selected_params <- spat_methods_params[[test]]

            if (length(selected_params) == 1) {
                if (is.na(selected_params)) {
                    if (selected_method == "binSpect_single") {
                        selected_params <- list(
                            bin_method = "kmeans",
                            nstart = 3,
                            iter_max = 10,
                            expression_values = "normalized",
                            get_av_expr = FALSE,
                            get_high_expr = FALSE
                        )
                    } else if (selected_method == "binSpect_multi") {
                        selected_params <- list(
                            bin_method = "kmeans",
                            spatial_network_k = c(5, 10, 20),
                            nstart = 3,
                            iter_max = 10,
                            expression_values = "normalized",
                            get_av_expr = FALSE,
                            get_high_expr = FALSE,
                            summarize = "adj.p.value"
                        )
                    } else if (selected_method == "spatialDE") {
                        selected_params <- list(
                            expression_values = "raw",
                            sig_alpha = 0.5,
                            unsig_alpha = 0.5,
                            show_plot = FALSE,
                            return_plot = FALSE,
                            save_plot = FALSE
                        )
                    } else if (selected_method == "spark") {
                        selected_params <- list(
                            expression_values = "raw",
                            return_object = "data.table",
                            percentage = 0.1,
                            min_count = 10,
                            num_core = 5
                        )
                    } else if (selected_method == "silhouetteRank") {
                        selected_params <- list(
                            expression_values = "normalized",
                            overwrite_input_bin = FALSE,
                            rbp_ps = c(0.95, 0.99),
                            examine_tops = c(0.005, 0.010),
                            matrix_type = "dissim",
                            num_core = 4,
                            parallel_path = "/usr/bin",
                            output = NULL,
                            query_sizes = 10L
                        )
                    }
                }
            }

            # name
            selected_name <- spat_methods_names[test]


            ## RUN Spatial Analysis ##
            if (selected_method == "binSpect_single") {
                start <- proc.time()
                spatial_gene_results <- do.call("binSpectSingle", c(
                    gobject = simulate_patch,
                    selected_params
                ))

                spatial_gene_results <- spatial_gene_results[feats == gene_name]
                total_time <- proc.time() - start

                spatial_gene_results[, prob := spatial_prob]
                spatial_gene_results[, time := total_time[["elapsed"]]]

                spatial_gene_results <- spatial_gene_results[
                    , .(feats, adj.p.value, prob, time)
                ]
                colnames(spatial_gene_results) <- c(
                    "feats", "adj.p.value", "prob", "time"
                )

                spatial_gene_results[, method := selected_name]
            } else if (selected_method == "binSpect_multi") {
                start <- proc.time()
                spatial_gene_results <- do.call("binSpectMulti", c(
                    gobject = simulate_patch,
                    selected_params
                ))

                spatial_gene_results <- spatial_gene_results$simple
                spatial_gene_results <- spatial_gene_results[feats == gene_name]
                total_time <- proc.time() - start

                spatial_gene_results[, prob := spatial_prob]
                spatial_gene_results[, time := total_time[["elapsed"]]]

                spatial_gene_results <- spatial_gene_results[
                    , .(feats, p.val, prob, time)
                ]
                colnames(spatial_gene_results) <- c(
                    "feats", "adj.p.value", "prob", "time"
                )

                spatial_gene_results[, method := selected_name]
            } else if (selected_method == "spatialDE") {
                start <- proc.time()
                new_raw_sim_matrix <- getExpression(simulate_patch,
                    values = "raw",
                    output = "matrix"
                )
                sd_cells <- apply(new_raw_sim_matrix, 2, sd)
                sd_non_zero_cells <- names(sd_cells[sd_cells != 0])
                simulate_patch_fix <- subsetGiotto(
                    simulate_patch,
                    cell_ids = sd_non_zero_cells
                )

                spatial_gene_results <- do.call("spatialDE", c(
                    gobject = simulate_patch_fix,
                    selected_params
                ))

                spatialDE_spatialgenes_sim_res <- spatial_gene_results$results$results
                if (is.null(spatialDE_spatialgenes_sim_res)) {
                    spatialDE_spatialgenes_sim_res <- spatial_gene_results$results
                }

                spatialDE_spatialgenes_sim_res <- data.table::as.data.table(
                    spatialDE_spatialgenes_sim_res
                )
                data.table::setorder(spatialDE_spatialgenes_sim_res, qval, pval)
                spatialDE_result <- spatialDE_spatialgenes_sim_res[
                    g == gene_name
                ]

                spatialDE_time <- proc.time() - start

                spatialDE_result[, prob := spatial_prob]
                spatialDE_result[, time := spatialDE_time[["elapsed"]]]

                spatial_gene_results <- spatialDE_result[
                    , .(g, qval, prob, time)
                ]
                colnames(spatial_gene_results) <- c(
                    "feats", "adj.p.value", "prob", "time"
                )
                spatial_gene_results[, method := "spatialDE"]
            } else if (selected_method == "spark") {
                ## spark
                start <- proc.time()
                spark_spatialgenes_sim <- do.call("spark", c(
                    gobject = simulate_patch,
                    selected_params
                ))

                spark_result <- spark_spatialgenes_sim[genes == gene_name]
                spark_time <- proc.time() - start

                spark_result[, prob := spatial_prob]
                spark_result[, time := spark_time[["elapsed"]]]

                spatial_gene_results <- spark_result[
                    , .(genes, adjusted_pvalue, prob, time)
                ]
                colnames(spatial_gene_results) <- c(
                    "genes", "adj.p.value", "prob", "time"
                )
                spatial_gene_results[, method := "spark"]
            } else if (selected_method == "silhouetteRank") {
                ## silhouetterank
                start <- proc.time()

                spatial_gene_results <- do.call("silhouetteRankTest", c(
                    gobject = simulate_patch,
                    selected_params
                ))

                data.table::setnames(
                    spatial_gene_results,
                    old = "gene", new = "genes"
                )
                spatial_gene_results <- spatial_gene_results[genes == gene_name]
                silh_time <- proc.time() - start

                spatial_gene_results[, prob := spatial_prob]
                spatial_gene_results[, time := silh_time[["elapsed"]]]

                # silhrank uses qval by default
                spatial_gene_results <- spatial_gene_results[
                    , .(genes, qval, prob, time)
                ]
                colnames(spatial_gene_results) <- c(
                    "genes", "adj.p.value", "prob", "time"
                )
                spatial_gene_results[, method := "silhouette"]
            }

            result_list[[test]] <- spatial_gene_results
        }

        results <- data.table::rbindlist(l = result_list)
        return(results)
    } else {
        return(NULL)
    }
}





#' @title run_spatial_sim_tests_multi
#' @name run_spatial_sim_tests_multi
#' @description runs all spatial tests for multiple probabilities and
#' repetitions
#' @returns data.table
#' @keywords internal
run_spatial_sim_tests_multi <- function(
        gobject,
        pattern_name = "pattern",
        pattern_cell_ids = NULL,
        gene_name = NULL,
        spatial_probs = c(0.5, 1),
        reps = 2,
        spatial_network_name = "kNN_network",
        spat_methods = c(
            "binSpect_single", "binSpect_multi", "spatialDE",
            "spark", "silhouetteRank"
        ),
        spat_methods_params = list(NA, NA, NA, NA, NA),
        spat_methods_names = c(
            "binSpect_single", "binSpect_multi", "spatialDE",
            "spark", "silhouetteRank"
        ),
        save_plot = FALSE,
        save_raw = FALSE,
        save_norm = FALSE,
        save_dir = "~",
        verbose = TRUE,
        run_simulations = TRUE,
        ...) {
    prob_list <- list()
    for (prob_ind in seq_along(spatial_probs)) {
        prob_i <- spatial_probs[prob_ind]

        if (verbose) message("start with ", prob_i)

        rep_list <- list()
        for (rep_i in seq_len(reps)) {
            if (verbose) message("repetition = ", rep_i)


            plot_name <- paste0(
                "plot_", gene_name, "_prob",
                prob_i, "_rep", rep_i
            )


            rep_res <- run_spatial_sim_tests_one_rep(gobject,
                pattern_name = pattern_name,
                pattern_cell_ids = pattern_cell_ids,
                gene_name = gene_name,
                spatial_prob = prob_i,
                spatial_network_name = spatial_network_name,
                spat_methods = spat_methods,
                spat_methods_params = spat_methods_params,
                spat_methods_names = spat_methods_names,
                save_plot = save_plot,
                save_raw = save_raw,
                save_norm = save_norm,
                save_dir = save_dir,
                save_name = plot_name,
                run_simulations = run_simulations,
                ...
            )

            if (run_simulations == TRUE) {
                rep_res[, rep := rep_i]
                rep_list[[rep_i]] <- rep_res
            }
        }

        if (run_simulations == TRUE) {
            rep_list_res <- do.call("rbind", rep_list)
            prob_list[[prob_ind]] <- rep_list_res
        }
    }

    if (run_simulations == TRUE) {
        final_gene_results <- do.call("rbind", prob_list)
        return(final_gene_results)
    }
}




#' @title runPatternSimulation
#' @name runPatternSimulation
#' @description Creates a known spatial pattern for selected genes one-by-one
#' and runs the different spatial gene detection tests
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_names selected genes
#' @param spatial_probs probabilities to test for a high expressing gene
#' value to be part of the spatial pattern
#' @param reps number of random simulation repetitions
#' @param spatial_network_name which spatial network to use for binSpectSingle
#' @param spat_methods vector of spatial methods to test
#' @param spat_methods_params list of parameters list for each element in the
#' vector of spatial methods to test
#' @param spat_methods_names name for each element in the vector of spatial
#' elements to test
#' @param scalefactor library size scaling factor when re-normalizing dataset
#' @param save_plot save intermediate random simulation plots or not
#' @param save_raw save the raw expression matrix of the simulation
#' @param save_norm save the normalized expression matrix of the simulation
#' @param save_dir directory to save results to
#' @param max_col maximum number of columns for final plots
#' @param height height of final plots
#' @param width width of final plots
#' @param run_simulations run simulations (default = TRUE)
#' @param \dots additional parameters for renormalization
#' @returns data.table with results
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' runPatternSimulation(
#'     gobject = g, pattern_cell_ids = c(
#'         "AAAGGGATGTAGCAAG-1",
#'         "TCAAACAACCGCGTCG-1", "ACGATCATACATAGAG-1", "TATGCTCCCTACTTAC-1"
#'     ),
#'     spatial_network_name = "spatial_network",
#'     gene_names = c("Gna12", "Ccnd2")
#' )
#' @export
runPatternSimulation <- function(
        gobject,
        pattern_name = "pattern",
        pattern_colors = c("in" = "green", "out" = "red"),
        pattern_cell_ids = NULL,
        gene_names = NULL,
        spatial_probs = c(0.5, 1),
        reps = 2,
        spatial_network_name = "kNN_network",
        spat_methods = c(
            "binSpect_single", "binSpect_multi", "spatialDE",
            "spark", "silhouetteRank"
        ),
        spat_methods_params = list(NA, NA, NA, NA, NA),
        spat_methods_names = c(
            "binSpect_single", "binSpect_multi", "spatialDE",
            "spark", "silhouetteRank"
        ),
        scalefactor = 6000,
        save_plot = TRUE,
        save_raw = TRUE,
        save_norm = TRUE,
        save_dir = "~",
        max_col = 4,
        height = 7,
        width = 7,
        run_simulations = TRUE,
        ...) {
    # data.table variables
    prob <- method <- adj.p.value <- time <- NULL


    # plot pattern for first gene (the same for all)
    example_patch <- simulateOneGenePatternGiottoObject(gobject,
        pattern_name = pattern_name,
        pattern_cell_ids = pattern_cell_ids,
        gene_name = gene_names[1],
        spatial_prob = 1,
        normalization_params = list(
            scalefactor = scalefactor,
            verbose = TRUE
        )
    )

    spatPlot2D(example_patch,
        cell_color = pattern_name, cell_color_code = pattern_colors,
        save_plot = save_plot, save_param = list(
            save_dir = save_dir,
            save_folder = "original",
            save_name = paste0(pattern_name, "_pattern"),
            base_width = 9,
            base_height = 7,
            units = "cm"
        )
    )


    all_results <- list()
    for (gene_ind in seq_along(gene_names)) {
        gene <- gene_names[gene_ind]

        # plot original expression
        GiottoVisuals::spatFeatPlot2D(
            gobject = gobject,
            expression_values = "normalized",
            feats = gene,
            point_shape = "border",
            point_border_stroke = 0.1,
            show_network = FALSE,
            network_color = "lightgrey",
            point_size = 2.5,
            cow_n_col = 1,
            show_plot = FALSE,
            save_plot = save_plot,
            save_param = list(
                save_dir = save_dir, save_folder = "original",
                save_name = paste0(gene, "_original"),
                base_width = 9, base_height = 7,
                units = "cm"
            )
        )


        generesults <- run_spatial_sim_tests_multi(
            gobject,
            pattern_name = pattern_name,
            pattern_cell_ids = pattern_cell_ids,
            gene_name = gene,
            spatial_network_name = spatial_network_name,
            spat_methods = spat_methods,
            spat_methods_params = spat_methods_params,
            spat_methods_names = spat_methods_names,
            save_plot = save_plot,
            save_raw = save_raw,
            save_norm = save_norm,
            save_dir = save_dir,
            spatial_probs = spatial_probs,
            reps = reps,
            run_simulations = run_simulations,
            ...
        )

        if (run_simulations == TRUE) {
            generesults[, prob := as.factor(prob)]
            uniq_methods <- mixedsort(unique(generesults$method))
            generesults[, method := factor(method, levels = uniq_methods)]

            if (save_plot == TRUE) {
                subdir <- paste0(save_dir, "/", pattern_name, "/")
                if (!file.exists(subdir)) {
                    dir.create(
                        path = subdir, recursive = TRUE
                    )
                }
                # write results
                data.table::fwrite(
                    x = generesults,
                    file = paste0(subdir, "/", gene, "_results.txt"),
                    sep = "\t", quote = FALSE
                )
            }

            all_results[[gene_ind]] <- generesults
        }
    }


    ## create combined results and visuals
    if (run_simulations == TRUE) {
        results <- do.call("rbind", all_results)

        ## plot results ##

        if (save_plot == TRUE) {
            # 4 columns max
            nr_rows <- max(c(round(length(gene_names) / max_col), 1))

            # p-values
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::geom_boxplot(
                data = results,
                ggplot2::aes(x = method, y = adj.p.value, color = prob)
            )
            pl <- pl + ggplot2::geom_point(
                data = results,
                ggplot2::aes(x = method, y = adj.p.value, color = prob),
                size = 2, position = ggplot2::position_jitterdodge()
            )
            pl <- pl + ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(
                    angle = 90, vjust = 1, hjust = 1
                ))
            pl <- pl + ggplot2::facet_wrap(~genes, nrow = nr_rows)
            pl <- pl + ggplot2::geom_hline(
                yintercept = 0.05, color = "red", linetype = 2
            )

            grDevices::pdf(
                file = paste0(
                    save_dir, "/", pattern_name, "_boxplot_pvalues.pdf"
                ),
                width = width, height = height
            )
            print(pl)
            grDevices::dev.off()



            # -log10 p-values
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::geom_boxplot(
                data = results,
                ggplot2::aes(x = method, y = -log10(adj.p.value), color = prob)
            )
            pl <- pl + ggplot2::geom_point(
                data = results,
                ggplot2::aes(x = method, y = -log10(adj.p.value), color = prob),
                size = 2, position = ggplot2::position_jitterdodge()
            )
            pl <- pl + ggplot2::theme_bw() + ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    angle = 90, vjust = 1, hjust = 1
                )
            )
            pl <- pl + ggplot2::facet_wrap(~genes, nrow = nr_rows)

            grDevices::pdf(
                file = paste0(
                    save_dir, "/", pattern_name, "_boxplot_log10pvalues.pdf"
                ),
                width = width, height = height
            )
            print(pl)
            grDevices::dev.off()


            # time
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::geom_boxplot(
                data = results,
                ggplot2::aes(x = method, y = time, color = prob)
            )
            pl <- pl + ggplot2::geom_point(
                data = results,
                ggplot2::aes(x = method, y = time, color = prob), size = 2,
                position = ggplot2::position_jitterdodge()
            )
            pl <- pl + ggplot2::theme_bw() + ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    angle = 90, vjust = 1, hjust = 1
                )
            )

            grDevices::pdf(
                file = paste0(
                    save_dir, "/", pattern_name, "_boxplot_time.pdf"
                ),
                width = width, height = height
            )
            print(pl)
            grDevices::dev.off()
        }


        # write results
        data.table::fwrite(
            x = results,
            file = paste0(save_dir, "/", pattern_name, "_results.txt"),
            sep = "\t", quote = FALSE
        )
        return(results)
    } else {
        return(NULL)
    }
}

# cell type proximity enrichment ####

#' @title make_simulated_network
#' @name make_simulated_network
#' @description Simulate random network.
#' @returns data.table
#' @keywords internal
make_simulated_network <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    spatial_network_name = "Delaunay_network",
    cluster_column,
    number_of_simulations = 100,
    set_seed = TRUE,
    seed_number = 1234) {
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

    # data.table variables
    unified_cells <- NULL

    spatial_network_annot <- annotateSpatialNetwork(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )

    # remove double edges between same cells #
    spatial_network_annot <- dt_sort_combine_two_columns(spatial_network_annot,
        column1 = "from", column2 = "to",
        myname = "unified_cells"
    )
    spatial_network_annot <- spatial_network_annot[!duplicated(unified_cells)]

    # create a simulated network
    length_ints <- nrow(spatial_network_annot)
    s1_list <- list()
    s2_list <- list()

    all_cell_type <- c(
        spatial_network_annot$from_cell_type,
        spatial_network_annot$to_cell_type
    )
    middle_point <- length(all_cell_type) / 2

    for (sim in seq_len(number_of_simulations)) {
        if (set_seed == TRUE) {
            seed_number <- seed_number + sim
            set.seed(seed = seed_number)
        }

        reshuffled_all_cell_type <- sample(
            x = all_cell_type, size = length(all_cell_type), replace = FALSE
        )

        new_from_cell_type <- reshuffled_all_cell_type[seq_len(middle_point)]
        s1_list[[sim]] <- new_from_cell_type

        new_to_cell_type <- reshuffled_all_cell_type[
            (middle_point + 1):length(all_cell_type)
        ]
        s2_list[[sim]] <- new_to_cell_type
    }

    s1_vector <- do.call("c", s1_list)
    s2_vector <- do.call("c", s2_list)
    round_vector <- rep(x = seq_len(number_of_simulations), each = length_ints)
    round_vector <- paste0("sim", round_vector)

    # data.table variables
    s1 <- s2 <- unified_int <- type_int <- NULL

    sample_dt <- data.table::data.table(
        s1 = s1_vector, s2 = s2_vector, round = round_vector
    )
    uniq_sim_comb <- unique(sample_dt[, .(s1, s2)])
    uniq_sim_comb[, unified_int := paste(
        sort(c(s1, s2)),
        collapse = "--"
    ), by = seq_len(nrow(uniq_sim_comb))]
    sample_dt[uniq_sim_comb, unified_int := unified_int, on = c(
        s1 = "s1", s2 = "s2"
    )]
    sample_dt[, type_int := ifelse(s1 == s2, "homo", "hetero")]

    return(sample_dt)
}



#' @title cellProximityEnrichment
#' @name cellProximityEnrichment
#' @description Compute cell-cell interaction enrichment (observed vs expected)
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param number_of_simulations number of simulations to create expected
#' observations
#' @param adjust_method method to adjust p.values
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @returns List of cell Proximity scores (CPscores) in data.table format. The
#' first data.table (raw_sim_table) shows the raw observations of both the
#' original and simulated networks. The second data.table (enrichm_res) shows
#' the enrichment results.
#' @details Spatial proximity enrichment or depletion between pairs of cell
#' types is calculated by calculating the observed over the expected frequency
#' of cell-cell proximity interactions. The expected frequency is the average
#' frequency calculated from a number of spatial network simulations. Each
#' individual simulation is obtained by reshuffling the cell type labels of
#' each node (cell) in the spatial network.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' cellProximityEnrichment(g, cluster_column = "leiden_clus")
#' @export
cellProximityEnrichment <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    spatial_network_name = "Delaunay_network",
    cluster_column,
    number_of_simulations = 1000,
    adjust_method = c(
        "none", "fdr", "bonferroni", "BH",
        "holm", "hochberg", "hommel",
        "BY"
    ),
    set_seed = TRUE,
    seed_number = 1234) {
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

    # p.adj test
    sel_adjust_method <- match.arg(adjust_method, choices = c(
        "none", "fdr", "bonferroni", "BH",
        "holm", "hochberg", "hommel",
        "BY"
    ))

    spatial_network_annot <- annotateSpatialNetwork(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )

    # remove double edges between same cells #
    # a simplified network does not have double edges between cells #

    # data.table variables
    unified_cells <- type_int <- N <- NULL

    spatial_network_annot <- dt_sort_combine_two_columns(
        spatial_network_annot, "to", "from", "unified_cells"
    )
    spatial_network_annot <- spatial_network_annot[!duplicated(unified_cells)]

    sample_dt <- make_simulated_network(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column,
        number_of_simulations = number_of_simulations,
        set_seed = set_seed,
        seed_number = seed_number
    )

    # combine original and simulated network
    table_sim_results <- sample_dt[, .N, by = c(
        "unified_int", "type_int", "round"
    )]

    ## create complete simulations
    ## add 0 if no single interaction was found
    unique_ints <- unique(table_sim_results[, .(unified_int, type_int)])

    # data.table with 0's for all interactions
    minimum_simulations <- unique_ints[rep(
        seq_len(nrow(unique_ints)), number_of_simulations
    ), ]
    minimum_simulations[, round := rep(
        paste0("sim", seq_len(number_of_simulations)),
        each = nrow(unique_ints)
    )]
    minimum_simulations[, N := 0]

    table_sim_minimum_results <- rbind(table_sim_results, minimum_simulations)
    table_sim_minimum_results[, V1 := sum(N), by = c(
        "unified_int", "type_int", "round"
    )]
    table_sim_minimum_results <- unique(
        table_sim_minimum_results[, .(unified_int, type_int, round, V1)]
    )
    table_sim_results <- table_sim_minimum_results


    # data.table variables
    orig <- unified_int <- V1 <- original <- enrichm <- simulations <- NULL

    table_sim_results[, orig := "simulations"]
    spatial_network_annot[, round := "original"]

    table_orig_results <- spatial_network_annot[, .N, by = c(
        "unified_int", "type_int", "round"
    )]
    table_orig_results[, orig := "original"]
    data.table::setnames(table_orig_results, old = "N", new = "V1")

    table_results <- rbind(table_orig_results, table_sim_results)



    # add missing combinations from original or simulations
    # probably not needed anymore
    all_simulation_ints <- as.character(unique(table_results[
        orig == "simulations"
    ]$unified_int))
    all_original_ints <- as.character(unique(table_results[
        orig == "original"
    ]$unified_int))
    missing_in_original <- all_simulation_ints[
        !all_simulation_ints %in% all_original_ints
    ]
    missing_in_simulations <- all_original_ints[
        !all_original_ints %in% all_simulation_ints
    ]
    create_missing_for_original <- table_results[
        unified_int %in% missing_in_original
    ]
    create_missing_for_original <- unique(create_missing_for_original[
        , c("orig", "V1") := list("original", 0)
    ])
    create_missing_for_simulations <- table_results[
        unified_int %in% missing_in_simulations
    ]
    create_missing_for_simulations <- unique(
        create_missing_for_simulations[, c("orig", "V1") := list(
            "simulations", 0
        )]
    )

    table_results <- do.call(
        "rbind",
        list(
            table_results, create_missing_for_original,
            create_missing_for_simulations
        )
    )


    ## p-values
    combo_list <- rep(NA, length = length(unique(table_results$unified_int)))
    p_high <- rep(NA, length = length(unique(table_results$unified_int)))
    p_low <- rep(NA, length = length(unique(table_results$unified_int)))

    for (int_combo in seq_along(unique(table_results$unified_int))) {
        this_combo <- as.character(unique(table_results$unified_int)[int_combo])

        sub <- table_results[unified_int == this_combo]

        orig_value <- sub[orig == "original"]$V1
        sim_values <- sub[orig == "simulations"]$V1

        length_simulations <- length(sim_values)
        if (length_simulations != number_of_simulations) {
            additional_length_needed <- number_of_simulations -
                length_simulations
            sim_values <- c(sim_values, rep(0, additional_length_needed))
        }

        p_orig_higher <- 1 - (sum((orig_value + 1) > (sim_values + 1)) /
            number_of_simulations)
        p_orig_lower <- 1 - (sum((orig_value + 1) < (sim_values + 1)) /
            number_of_simulations)

        combo_list[[int_combo]] <- this_combo
        p_high[[int_combo]] <- p_orig_higher
        p_low[[int_combo]] <- p_orig_lower
    }
    res_pvalue_DT <- data.table::data.table(
        unified_int = as.vector(combo_list),
        p_higher_orig = p_high,
        p_lower_orig = p_low
    )


    # depletion or enrichment in barplot format
    table_mean_results <- table_results[, .(mean(V1)), by = c(
        "orig", "unified_int", "type_int"
    )]
    table_mean_results_dc <- data.table::dcast.data.table(
        data = table_mean_results, formula = type_int + unified_int ~ orig,
        value.var = "V1"
    )
    table_mean_results_dc[, original := ifelse(is.na(original), 0, original)]
    table_mean_results_dc[, enrichm := log2((original + 1) / (simulations + 1))]


    table_mean_results_dc <- merge(
        table_mean_results_dc, res_pvalue_DT,
        by = "unified_int"
    )
    data.table::setorder(table_mean_results_dc, enrichm)
    table_mean_results_dc[, unified_int := factor(unified_int, unified_int)]

    # adjust p-values for mht

    # data.table variables
    p.adj_higher <- p.adj_lower <- p_lower_orig <- p_higher_orig <-
        PI_value <- int_ranking <- NULL

    table_mean_results_dc[, p.adj_higher := stats::p.adjust(
        p_higher_orig,
        method = sel_adjust_method
    )]
    table_mean_results_dc[, p.adj_lower := stats::p.adjust(
        p_lower_orig,
        method = sel_adjust_method
    )]


    table_mean_results_dc[, PI_value := ifelse(p.adj_higher <= p.adj_lower,
        -log10(p.adj_higher + (1 / number_of_simulations)) * enrichm,
        -log10(p.adj_lower + (1 / number_of_simulations)) * enrichm
    )]
    data.table::setorder(table_mean_results_dc, PI_value)

    # order
    table_mean_results_dc <- table_mean_results_dc[order(-PI_value)]
    table_mean_results_dc[, int_ranking := seq_len(.N)]

    return(list(
        raw_sim_table = table_results,
        enrichm_res = table_mean_results_dc
    ))
}



#' @title addCellIntMetadata
#' @name addCellIntMetadata
#' @description Creates an additional metadata column with information about
#' interacting and non-interacting cell types of the
#' selected cell-cell interaction.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spatial_network name of spatial network to use
#' @param cluster_column column of cell types
#' @param cell_interaction cell-cell interaction to use
#' @param name name for the new metadata column
#' @param return_gobject return an updated giotto object
#' @returns Giotto object
#' @details This function will create an additional metadata column which
#' selects interacting cell types for a specific cell-cell
#' interaction. For example, if you want to color interacting astrocytes and
#' oligodendrocytes it will create a new metadata column with
#' the values "select_astrocytes", "select_oligodendrocytes",
#' "other_astrocytes", "other_oligodendroyctes" and "other". Where "other" is
#' all other cell types found within the selected cell type column.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addCellIntMetadata(g,
#'     cluster_column = "leiden_clus",
#'     cell_interaction = "custom_leiden"
#' )
#' @export
addCellIntMetadata <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    spatial_network = "spatial_network",
    cluster_column,
    cell_interaction,
    name = "select_int",
    return_gobject = TRUE) {
    # set spatial unit and feature type
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )


    if (is.null(spatial_network)) {
        stop("spatial_network must be provided, this must be an existing
            spatial network")
    }

    if (is.null(cluster_column)) {
        stop("cluster_column must be provided, this must be an existing cell
            metadata column, see pData(your_giotto_object)")
    }

    if (is.null(cell_interaction)) {
        stop("cell_interaction must be provided, this must be cell--cell
            interaction between cell types in cluster_column")
    }

    # create spatial network
    spatial_network_annot <- annotateSpatialNetwork(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        spatial_network_name = spatial_network,
        cluster_column = cluster_column
    )

    # selected vs other cells

    # data.table variables
    unified_int <- cell_ID <- NULL

    selected_cells <- unique(c(
        spatial_network_annot[unified_int == cell_interaction]$to,
        spatial_network_annot[unified_int == cell_interaction]$from
    ))

    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "cellMetaObj",
        copy_obj = TRUE
    )

    cell_type_1 <- strsplit(cell_interaction, split = "--")[[1]][1]
    cell_type_2 <- strsplit(cell_interaction, split = "--")[[1]][2]

    cell_metadata[][, c(name) := ifelse(!get(cluster_column) %in% c(
        cell_type_1, cell_type_2
    ), "other",
    ifelse(get(cluster_column) == cell_type_1 & cell_ID %in% selected_cells,
        paste0("select_", cell_type_1),
        ifelse(get(cluster_column) == cell_type_2 & cell_ID %in%
            selected_cells, paste0("select_", cell_type_2),
        ifelse(get(cluster_column) == cell_type_1,
            paste0("other_", cell_type_1),
            paste0("other_", cell_type_2)
        )
        )
    )
    )]

    if (return_gobject == TRUE) {
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, cell_metadata, verbose = FALSE)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        ## update parameters used ##
        gobject <- update_giotto_params(gobject,
            description = "_add_cell_int_info"
        )

        return(gobject)
    } else {
        return(cell_metadata)
    }
}





# * ####
# ICF cell proximity feature changes ####

#' @title Interaction changed features test methods
#' @name cell_proximity_tests
#' @description
#' Perform specified test on subsets of a matrix
#' @returns cell proximity values
NULL


#' @describeIn cell_proximity_tests t.test
#' @keywords internal
.do_ttest <- function(expr_values,
    select_ind,
    other_ind,
    adjust_method,
    mean_method,
    offset = 0.1) {
    vmsg(.is_debug = TRUE, ".do_ttest")

    # data.table variables
    p.value <- p.adj <- NULL

    mean_sel <- my_rowMeans(
        expr_values[, select_ind],
        method = mean_method, offset = offset
    )
    mean_all <- my_rowMeans(
        expr_values[, other_ind],
        method = mean_method, offset = offset
    )

    if (length(select_ind) == 1 | length(other_ind) == 1) {
        results <- NaN
    } else {
        results <- apply(expr_values, MARGIN = 1, function(x) {
            p.value <- stats::t.test(x[select_ind], x[other_ind])$p.value
        })
    }

    # other info
    log2fc <- log2((mean_sel + offset) / (mean_all + offset))
    diff <- mean_sel - mean_all

    resultsDT <- data.table(
        "feats" = rownames(expr_values),
        "sel" = mean_sel, "other" = mean_all,
        "log2fc" = log2fc, "diff" = diff,
        "p.value" = unlist(results)
    )
    resultsDT[, p.value := ifelse(is.nan(p.value), 1, p.value)]
    resultsDT[, p.adj := stats::p.adjust(p.value, method = adjust_method)]
    setorder(resultsDT, p.adj)

    return(resultsDT)
}




#' @describeIn cell_proximity_tests limma t.test
#' @keywords internal
.do_limmatest <- function(expr_values,
    select_ind,
    other_ind,
    mean_method,
    offset = 0.1) {
    vmsg(.is_debug = TRUE, ".do_limmatest")

    package_check("limma")

    # data.table variables
    sel <- other <- feats <- P.Value <- adj.P.Val <- p.adj <- NULL

    expr_values_subset <- cbind(
        expr_values[, select_ind], expr_values[, other_ind]
    )
    mygroups <- c(
        rep("sel", length(select_ind)),
        rep("other", length(other_ind))
    )
    mygroups <- factor(mygroups, levels = unique(mygroups))

    design <- stats::model.matrix(~ 0 + mygroups)
    colnames(design) <- levels(mygroups)
    fit <- limma::lmFit(expr_values_subset, design = design)


    cont.matrix <- limma::makeContrasts(
        sel_vs_other = sel - other,
        levels = design
    )

    fitcontrast <- limma::contrasts.fit(fit, cont.matrix)
    fitc_ebayes <- limma::eBayes(fitcontrast)

    # limma to DT
    limma_result <- limma::topTable(
        fitc_ebayes,
        coef = 1, number = 100000, confint = TRUE
    )
    limmaDT <- data.table::as.data.table(limma_result)
    limmaDT[, feats := rownames(limma_result)]

    # other info
    mean_sel <- my_rowMeans(
        expr_values[, select_ind],
        method = mean_method, offset = offset
    )
    mean_all <- my_rowMeans(
        expr_values[, other_ind],
        method = mean_method, offset = offset
    )

    log2fc <- log2((mean_sel + offset) / (mean_all + offset))
    diff <- mean_sel - mean_all

    tempDT <- data.table::data.table(
        "feats" = rownames(expr_values),
        "sel" = mean_sel,
        "other" = mean_all,
        "log2fc" = log2fc,
        "diff" = diff
    )
    limmaDT <- data.table::merge.data.table(limmaDT, tempDT, by = "feats")
    limmaDT <- limmaDT[
        , .(feats, sel, other, log2fc, diff, P.Value, adj.P.Val)
    ]
    colnames(limmaDT) <- c(
        "feats", "sel", "other", "log2fc", "diff",
        "p.value", "p.adj"
    )

    setorder(limmaDT, p.adj)

    return(limmaDT)
}




#' @describeIn cell_proximity_tests wilcoxon
#' @keywords internal
.do_wilctest <- function(expr_values,
    select_ind,
    other_ind,
    adjust_method,
    mean_method,
    offset = 0.1) {
    vmsg(.is_debug = TRUE, ".do_wilctest")

    # data.table variables
    p.value <- p.adj <- NULL

    mean_sel <- my_rowMeans(
        expr_values[, select_ind],
        method = mean_method, offset = offset
    )
    mean_all <- my_rowMeans(
        expr_values[, other_ind],
        method = mean_method, offset = offset
    )

    if (length(select_ind) == 1 | length(other_ind) == 1) {
        results <- NaN
    } else {
        results <- apply(expr_values, MARGIN = 1, function(x) {
            p.value <- stats::wilcox.test(x[select_ind], x[other_ind])$p.value
        })
    }

    # other info
    log2fc <- log2((mean_sel + offset) / (mean_all + offset))
    diff <- mean_sel - mean_all

    resultsDT <- data.table(
        "feats" = rownames(expr_values),
        "sel" = mean_sel,
        "other" = mean_all,
        "log2fc" = log2fc,
        "diff" = diff,
        "p.value" = unlist(results)
    )
    resultsDT[, p.value := ifelse(is.nan(p.value), 1, p.value)]
    resultsDT[, p.adj := stats::p.adjust(p.value, method = adjust_method)]
    setorder(resultsDT, p.adj)

    return(resultsDT)
}


# calculate original values
.do_permuttest_original <- function(expr_values,
    select_ind,
    other_ind,
    name = "orig",
    mean_method,
    offset = 0.1) {
    # data.table variables
    feats <- NULL

    mean_sel <- my_rowMeans(expr_values[
        , select_ind
    ], method = mean_method, offset = offset)
    mean_all <- my_rowMeans(expr_values[
        , other_ind
    ], method = mean_method, offset = offset)

    log2fc <- log2((mean_sel + offset) / (mean_all + offset))
    diff <- mean_sel - mean_all

    resultsDT <- data.table(
        "sel" = mean_sel, "other" = mean_all, "log2fc" = log2fc, "diff" = diff
    )
    resultsDT[, feats := rownames(expr_values)]
    resultsDT[, name := name]

    return(resultsDT)
}



# calculate random values
.do_permuttest_random <- function(expr_values,
    select_ind,
    other_ind,
    name = "perm_1",
    mean_method,
    offset = 0.1,
    set_seed = TRUE,
    seed_number = 1234) {
    # data.table variables
    feats <- NULL

    l_select_ind <- length(select_ind)
    l_other_ind <- length(other_ind)

    all_ind <- c(select_ind, other_ind)

    if (set_seed == TRUE) {
        set.seed(seed = seed_number)
    }
    random_select <- sample(all_ind, size = l_select_ind, replace = FALSE)
    random_other <- all_ind[!all_ind %in% random_select]

    # alternative
    mean_sel <- my_rowMeans(
        expr_values[, random_select],
        method = mean_method, offset = offset
    )
    mean_all <- my_rowMeans(
        expr_values[, random_other],
        method = mean_method, offset = offset
    )

    log2fc <- log2((mean_sel + offset) / (mean_all + offset))
    diff <- mean_sel - mean_all

    resultsDT <- data.table(
        "sel" = mean_sel, "other" = mean_all, "log2fc" = log2fc, "diff" = diff
    )
    resultsDT[, feats := rownames(expr_values)]
    resultsDT[, name := name]

    return(resultsDT)
}




# calculate multiple random values
.do_multi_permuttest_random <- function(expr_values,
    select_ind,
    other_ind,
    mean_method,
    offset = 0.1,
    n = 100,
    set_seed = TRUE,
    seed_number = 1234) {
    if (set_seed == TRUE) {
        seed_number_list <- seed_number:(seed_number + (n - 1))
    }

    result <- lapply(X = seq_len(n), FUN = function(x) {
        seed_number <- seed_number_list[x]

        perm_rand <- .do_permuttest_random(
            expr_values = expr_values,
            select_ind = select_ind,
            other_ind = other_ind,
            name = paste0("perm_", x),
            mean_method = mean_method,
            offset = offset,
            set_seed = set_seed,
            seed_number = seed_number
        )
    })

    final_result <- do.call("rbind", result)
}


#' @describeIn cell_proximity_tests random permutation
#' @keywords internal
.do_permuttest <- function(expr_values,
    select_ind, other_ind,
    n_perm = 1000,
    adjust_method = "fdr",
    mean_method,
    offset = 0.1,
    set_seed = TRUE,
    seed_number = 1234) {
    # data.table variables
    log2fc_diff <- log2fc <- sel <- other <- feats <- p_higher <- p_lower <-
        perm_sel <- NULL
    perm_other <- perm_log2fc <- perm_diff <- p.value <- p.adj <- NULL

    ## original data
    vmsg(.is_debug = TRUE, "ok1")
    original <- .do_permuttest_original(
        expr_values = expr_values,
        select_ind = select_ind, other_ind = other_ind,
        name = "orig",
        mean_method = mean_method,
        offset = offset
    )

    ## random permutations
    vmsg(.is_debug = TRUE, "ok2")
    random_perms <- .do_multi_permuttest_random(
        expr_values = expr_values,
        n = n_perm,
        select_ind = select_ind,
        other_ind = other_ind,
        mean_method = mean_method,
        offset = offset,
        set_seed = set_seed,
        seed_number = seed_number
    )

    ##
    random_perms[, log2fc_diff := rep(original$log2fc, n_perm) - log2fc]
    random_perms[,
        c("perm_sel", "perm_other", "perm_log2fc", "perm_diff") := list(
            mean(sel), mean(other), mean(log2fc), mean(diff)
        ),
        by = feats
    ]

    ## get p-values
    random_perms[, p_higher := sum(log2fc_diff > 0), by = feats]
    random_perms[, p_higher := 1 - (p_higher / n_perm)]
    random_perms[, p_lower := sum(log2fc_diff < 0), by = feats]
    random_perms[, p_lower := 1 - (p_lower / n_perm)]

    ## combine results permutation and original
    random_perms_res <- unique(random_perms[
        , .(
            feats, perm_sel, perm_other, perm_log2fc, perm_diff, p_higher,
            p_lower
        )
    ])
    results_m <- data.table::merge.data.table(
        random_perms_res, original[, .(feats, sel, other, log2fc, diff)],
        by = "feats"
    )

    # select lowest p-value and perform p.adj
    results_m[, p.value := ifelse(p_higher <= p_lower, p_higher, p_lower)]
    results_m[, p.adj := stats::p.adjust(p.value, method = adjust_method)]

    results_m <- results_m[
        , .(
            feats, sel, other, log2fc, diff, p.value, p.adj, perm_sel,
            perm_other, perm_log2fc, perm_diff
        )
    ]
    setorder(results_m, p.adj, -log2fc)

    return(results_m)
}



#' @title Do cell proximity test
#' @name .do_cell_proximity_test
#' @description Performs a selected differential test on subsets of a matrix
#' @param expr_values Matrix object
#' @returns differential test on subsets of a matrix
#' @keywords internal
#' @seealso [cell_proximity_tests]
.do_cell_proximity_test <- function(expr_values,
    select_ind, other_ind,
    diff_test = c("permutation", "limma", "t.test", "wilcox"),
    mean_method = c("arithmic", "geometric"),
    offset = 0.1,
    n_perm = 100,
    adjust_method = c(
        "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "fdr", "none"
    ),
    set_seed = TRUE,
    seed_number = 1234,
    verbose = FALSE) {
    # get parameters
    diff_test <- match.arg(
        diff_test,
        choices = c("permutation", "limma", "t.test", "wilcox")
    )
    adjust_method <- match.arg(adjust_method, choices = c(
        "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "fdr", "none"
    ))
    mean_method <- match.arg(mean_method, choices = c("arithmic", "geometric"))

    vmsg(.is_debug = TRUE, ".do_cell_proximity_test")

    if (diff_test == "permutation") {
        result <- .do_permuttest(
            expr_values = expr_values,
            select_ind = select_ind, other_ind = other_ind,
            n_perm = n_perm, adjust_method = adjust_method,
            mean_method = mean_method, offset = offset,
            set_seed = set_seed,
            seed_number = seed_number
        )
    } else if (diff_test == "limma") {
        result <- .do_limmatest(
            expr_values = expr_values,
            select_ind = select_ind, other_ind = other_ind,
            mean_method = mean_method, offset = offset
        )
    } else if (diff_test == "t.test") {
        result <- .do_ttest(
            expr_values = expr_values,
            select_ind = select_ind, other_ind = other_ind,
            mean_method = mean_method, offset = offset,
            adjust_method = adjust_method
        )
    } else if (diff_test == "wilcox") {
        result <- .do_wilctest(
            expr_values = expr_values,
            select_ind = select_ind, other_ind = other_ind,
            mean_method = mean_method, offset = offset,
            adjust_method = adjust_method
        )
    }

    return(result)
}




#' @title Find cell proximity features per interaction
#' @name .findCellProximityFeats_per_interaction
#' @description Identifies features that are differentially expressed due to
#' proximity to other cell types.
#' @param expr_values Matrix object
#' @returns data.table
#' @keywords internal
#' @seealso [.do_cell_proximity_test()] for specific tests
.findCellProximityFeats_per_interaction <- function(sel_int,
    expr_values,
    cell_metadata,
    annot_spatnetwork,
    cluster_column = NULL,
    minimum_unique_cells = 1,
    minimum_unique_int_cells = 1,
    exclude_selected_cells_from_test = TRUE,
    diff_test = c("permutation", "limma", "t.test", "wilcox"),
    mean_method = c("arithmic", "geometric"),
    offset = 0.1,
    adjust_method = "bonferroni",
    nr_permutations = 100,
    set_seed = TRUE,
    seed_number = 1234) {
    # data.table variables
    unified_int <- to_cell_type <- from_cell_type <- cell_type <-
        int_cell_type <- NULL
    nr_select <- int_nr_select <- nr_other <- int_nr_other <- unif_int <- NULL

    # select test to perform
    diff_test <- match.arg(
        arg = diff_test,
        choices = c("permutation", "limma", "t.test", "wilcox")
    )

    # select subnetwork
    sub_spatnetwork <- annot_spatnetwork[unified_int == sel_int]

    # unique cell types
    unique_cell_types <- unique(
        c(sub_spatnetwork$to_cell_type, sub_spatnetwork$from_cell_type)
    )

    if (length(unique_cell_types) == 2) {
        first_cell_type <- unique_cell_types[1]
        second_cell_type <- unique_cell_types[2]

        # first cell type ids
        to1 <- sub_spatnetwork[to_cell_type == first_cell_type][["to"]]
        from1 <- sub_spatnetwork[from_cell_type == first_cell_type][["from"]]
        cell1_ids <- unique(c(to1, from1))

        # second cell type ids
        to2 <- sub_spatnetwork[to_cell_type == second_cell_type][["to"]]
        from2 <- sub_spatnetwork[from_cell_type == second_cell_type][["from"]]
        cell2_ids <- unique(c(to2, from2))

        ## all cell ids
        all_cell1 <- cell_metadata[get(cluster_column) == first_cell_type][[
            "cell_ID"
        ]]
        all_cell2 <- cell_metadata[get(cluster_column) == second_cell_type][[
            "cell_ID"
        ]]

        ## exclude selected
        if (exclude_selected_cells_from_test == TRUE) {
            all_cell1 <- all_cell1[!all_cell1 %in% cell1_ids]
            all_cell2 <- all_cell2[!all_cell2 %in% cell2_ids]
        }

        ## FOR CELL TYPE 1
        sel_ind1 <- which(colnames(expr_values) %in% cell1_ids)
        all_ind1 <- which(colnames(expr_values) %in% all_cell1)

        ## FOR CELL TYPE 2
        sel_ind2 <- which(colnames(expr_values) %in% cell2_ids)
        all_ind2 <- which(colnames(expr_values) %in% all_cell2)


        ## do not continue if too few cells ##
        if (length(sel_ind1) < minimum_unique_cells |
            length(all_ind1) < minimum_unique_cells |
            length(sel_ind2) < minimum_unique_int_cells) {
            result_cell_1 <- NULL
        } else {
            result_cell_1 <- .do_cell_proximity_test(
                expr_values = expr_values,
                select_ind = sel_ind1,
                other_ind = all_ind1,
                diff_test = diff_test,
                n_perm = nr_permutations,
                mean_method = mean_method,
                offset = offset,
                adjust_method = adjust_method,
                set_seed = set_seed,
                seed_number = seed_number
            )
            result_cell_1[, cell_type := first_cell_type]
            result_cell_1[, int_cell_type := second_cell_type]
            result_cell_1[, nr_select := length(sel_ind1)]
            result_cell_1[, int_nr_select := length(sel_ind2)]
            result_cell_1[, nr_other := length(all_ind1)]
            result_cell_1[, int_nr_other := length(all_ind2)]
        }


        ## do not continue if too few cells ##
        if (length(sel_ind2) < minimum_unique_cells |
            length(all_ind2) < minimum_unique_cells |
            length(sel_ind1) < minimum_unique_int_cells) {
            result_cell_2 <- NULL
        } else {
            result_cell_2 <- .do_cell_proximity_test(
                expr_values = expr_values,
                select_ind = sel_ind2, other_ind = all_ind2,
                diff_test = diff_test,
                n_perm = nr_permutations,
                mean_method = mean_method,
                offset = offset,
                adjust_method = adjust_method,
                set_seed = set_seed,
                seed_number = seed_number
            )
            result_cell_2[, cell_type := second_cell_type]
            result_cell_2[, int_cell_type := first_cell_type]
            result_cell_2[, nr_select := length(sel_ind2)]
            result_cell_2[, int_nr_select := length(sel_ind1)]
            result_cell_2[, nr_other := length(all_ind2)]
            result_cell_2[, int_nr_other := length(all_ind1)]
        }


        ## COMBINE

        if (is.null(result_cell_1) & is.null(result_cell_2)) {
            return(NULL)
        } else {
            result_cells <- rbind(result_cell_1, result_cell_2)
        }
    } else if (length(unique_cell_types) == 1) {
        first_cell_type <- unique_cell_types[1]

        # first cell type ids
        to1 <- sub_spatnetwork[to_cell_type == first_cell_type][["to"]]
        from1 <- sub_spatnetwork[from_cell_type == first_cell_type][["from"]]
        cell1_ids <- unique(c(to1, from1))

        ## all cell ids
        all_cell1 <- cell_metadata[get(cluster_column) == first_cell_type][[
            "cell_ID"
        ]]

        ## exclude selected
        if (exclude_selected_cells_from_test == TRUE) {
            all_cell1 <- all_cell1[!all_cell1 %in% cell1_ids]
        }

        ## FOR CELL TYPE 1
        sel_ind1 <- which(colnames(expr_values) %in% cell1_ids)
        all_ind1 <- which(colnames(expr_values) %in% all_cell1)


        ## do not continue if too few cells ##
        if (length(sel_ind1) < minimum_unique_cells |
            length(all_ind1) < minimum_unique_cells) {
            return(NULL)
        }

        result_cells <- .do_cell_proximity_test(
            expr_values = expr_values,
            select_ind = sel_ind1, other_ind = all_ind1,
            diff_test = diff_test,
            n_perm = nr_permutations,
            mean_method = mean_method,
            offset = offset,
            adjust_method = adjust_method,
            set_seed = set_seed,
            seed_number = seed_number
        )

        result_cells[, cell_type := first_cell_type]
        result_cells[, int_cell_type := first_cell_type]
        result_cells[, nr_select := length(sel_ind1)]
        result_cells[, int_nr_select := length(sel_ind1)]
        result_cells[, nr_other := length(all_ind1)]
        result_cells[, int_nr_other := length(all_ind1)]
    }

    result_cells[, unif_int := sel_int]

    return(result_cells)
}




#' @title findInteractionChangedFeats
#' @name findInteractionChangedFeats
#' @aliases findICF
#' @description Identifies cell-to-cell Interaction Changed Features (ICF),
#' i.e. features that are differentially expressed due to proximity to other
#' cell types. This function is appropriate for single-cell level data. For
#' data from spot-based spatial assays or spatially binned data, see
#' [findICFSpot()], which runs on top of DWLS results or similar spot-level
#' cell-type enrichment outputs
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param selected_feats subset of selected features (optional)
#' @param cluster_column name of column to use for cell types
#' @param spatial_network_name name of spatial network to use
#' @param minimum_unique_cells minimum number of target cells required
#' @param minimum_unique_int_cells minimum number of interacting cells required
#' @param diff_test which differential expression test
#' @param mean_method method to use to calculate the mean
#' @param offset offset value to use when calculating log2 ratio
#' @param adjust_method which method to adjust p-values
#' @param nr_permutations number of permutations if diff_test = permutation
#' @param exclude_selected_cells_from_test exclude interacting cells other cells
#' @param do_parallel run calculations in parallel with mclapply
#' @param set_seed set a seed for reproducibility
#' @param seed_number seed number
#' @returns `icfObject` that contains the Interaction Changed differential
#' feature scores
#' @details Function to calculate if features are differentially expressed in
#' cell types when they interact (approximated by physical proximity) with
#' other cell types. The results data.table in the `icfObject` contains
#' - at least - the following columns:
#'   * **features:** All or selected list of tested features
#'   * **sel:** average feature expression in the interacting cells from the
#'   target cell type
#'   * **other:** average feature expression in the NOT-interacting cells from
#'   the target cell type
#'   * **log2fc:** log2 fold-change between sel and other
#'   * **diff:** spatial expression difference between sel and other
#'   * **p.value:** associated p-value
#'   * **p.adj:** adjusted p-value
#'   * **cell_type:** target cell type
#'   * **int_cell_type:** interacting cell type
#'   * **nr_select:** number of cells for selected target cell type
#'   * **int_nr_select:** number of cells for interacting cell type
#'   * **nr_other:** number of other cells of selected target cell type
#'   * **int_nr_other:** number of other cells for interacting cell type
#'   * **unif_int:** cell-cell interaction
#'
#' @seealso [filterInteractionChangedFeats()] [findICFSpot()]
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' icf1 <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"),
#'     nr_permutations = 10
#' )
#' force(icf1)
#' force(icf1$ICFscores)
#'
#' # this is just an alias with a shorter name
#' icf2 <- findICF(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"),
#'     nr_permutations = 10
#' )
#' @export
findInteractionChangedFeats <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = "normalized",
    selected_feats = NULL,
    cluster_column,
    spatial_network_name = "Delaunay_network",
    minimum_unique_cells = 1,
    minimum_unique_int_cells = 1,
    diff_test = c("permutation", "limma", "t.test", "wilcox"),
    mean_method = c("arithmic", "geometric"),
    offset = 0.1,
    adjust_method = c(
        "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "fdr", "none"
    ),
    nr_permutations = 1000,
    exclude_selected_cells_from_test = TRUE,
    do_parallel = TRUE,
    set_seed = TRUE,
    seed_number = 1234) {
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

    # expression values to be used
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
    ) # as matrix (never set back)

    ## test selected feats ##
    if (!is.null(selected_feats)) {
        expr_values <- expr_values[rownames(expr_values) %in% selected_feats, ]
    }

    ## stop test selected feats ##

    # difference test
    diff_test <- match.arg(
        diff_test,
        choices = c("permutation", "limma", "t.test", "wilcox")
    )

    # p.adj test
    adjust_method <- match.arg(adjust_method, choices = c(
        "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "fdr", "none"
    ))
    # how to calculate mean
    mean_method <- match.arg(mean_method, choices = c("arithmic", "geometric"))

    ## metadata
    cell_metadata <- pDataDT(
        gobject,
        spat_unit = spat_unit, feat_type = feat_type
    )


    ## annotated spatial network
    annot_spatnetwork <- annotateSpatialNetwork(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )

    all_interactions <- unique(annot_spatnetwork$unified_int)

    ## prepare function
    fcp_feats_per_i <- function(x) {
        .findCellProximityFeats_per_interaction(
            expr_values = expr_values,
            cell_metadata = cell_metadata,
            annot_spatnetwork = annot_spatnetwork,
            minimum_unique_cells = minimum_unique_cells,
            minimum_unique_int_cells = minimum_unique_int_cells,
            sel_int = x,
            cluster_column = cluster_column,
            exclude_selected_cells_from_test = exclude_selected_cells_from_test,
            diff_test = diff_test,
            mean_method = mean_method,
            offset = offset,
            adjust_method = adjust_method,
            nr_permutations = nr_permutations,
            set_seed = set_seed,
            seed_number = seed_number
        )
    }

    if (isTRUE(do_parallel)) { # parallel
        fin_result <- lapply_flex(
            X = all_interactions, future.seed = TRUE, FUN = fcp_feats_per_i
        )
    } else { # sequential
        fin_result <- list()

        for (i in seq_along(all_interactions)) {
            x <- all_interactions[i]
            tempres <- fcp_feats_per_i(x)
            fin_result[[i]] <- tempres
        }
    }

    final_result <- do.call("rbind", fin_result)

    # NSE variables
    spec_int <- cell_type <- int_cell_type <- type_int <- NULL

    final_result[, spec_int := paste0(cell_type, "--", int_cell_type)]
    final_result[, type_int := ifelse(
        cell_type == int_cell_type, "homo", "hetero"
    )]

    permutation_test <- ifelse(
        diff_test == "permutation", nr_permutations, "no permutations"
    )

    icfObject <- structure(
        .Data = list(
            ICFscores = final_result,
            Giotto_info = list(
                "values" = values,
                "cluster" = cluster_column,
                "spatial network" = spatial_network_name
            ),
            test_info = list(
                "test" = diff_test,
                "p.adj" = adjust_method,
                "min cells" = minimum_unique_cells,
                "min interacting cells" = minimum_unique_int_cells,
                "exclude selected cells" = exclude_selected_cells_from_test,
                "perm" = permutation_test
            )
        ),
        class = "icfObject"
    )
    return(icfObject)
}

#' @rdname findInteractionChangedFeats
#' @export
findICF <- findInteractionChangedFeats

#' @name print.icfObject
#' @title icfObject print method
#' @param x object to print
#' @param \dots additional params to pass (none implemented)
#' @keywords internal
#' @returns icfObject
#' @export
print.icfObject <- function(x, ...) {
    cat("An object of class", class(x), "\n")
    info <- list(
        dimensions = sprintf(
            "%d, %d (icfs, attributes)",
            nrow(x$ICFscores), ncol(x$ICFscores)
        )
    )
    print_list(info, pre = " -")
    cat("<giotto info>\n")
    print_list(x$Giotto_info, pre = " -")
    cat("<test info>\n")
    print_list(x$test_info, pre = " -")
}








#' @title filterInteractionChangedFeats
#' @name filterInteractionChangedFeats
#' @aliases filterICF
#' @description Filter Interaction Changed Feature scores.
#' @param icfObject ICF (interaction changed feature) score object
#' @param min_cells minimum number of source cell type
#' @param min_cells_expr minimum expression level for source cell type
#' @param min_int_cells minimum number of interacting neighbor cell type
#' @param min_int_cells_expr minimum expression level for interacting neighbor
#' cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum log2 fold-change
#' @param min_zscore minimum z-score change
#' @param zscores_column calculate z-scores over cell types or genes
#' @param direction differential expression directions to keep
#' @returns `icfObject` that contains the filtered differential feature scores
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' icf <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"),
#'     nr_permutations = 10
#' )
#' force(icf)
#' force(icf$ICFscores)
#'
#' icf_filter1 <- filterInteractionChangedFeats(icf, min_cells = 4)
#' force(icf_filter1)
#' force(icf_filter1$ICFscores)
#'
#' # filterICF is a simple alias with a shortened name
#' icf_filter2 <- filterICF(icf, min_cells = 4)
#' force(icf_filter2)
#'
#' @export
filterInteractionChangedFeats <- function(icfObject,
    min_cells = 4,
    min_cells_expr = 1,
    min_int_cells = 4,
    min_int_cells_expr = 1,
    min_fdr = 0.1,
    min_spat_diff = 0.2,
    min_log2_fc = 0.2,
    min_zscore = 2,
    zscores_column = c("cell_type", "feats"),
    direction = c("both", "up", "down")) {
    # NSE vars
    nr_select <- int_nr_select <- zscores <- log2fc <- sel <- other <-
        p.adj <- NULL

    if (!"icfObject" %in% class(icfObject)) {
        stop("icfObject needs to be the output from
            findInteractionChangedFeats() or findICF()")
    }

    zscores_column <- match.arg(
        zscores_column,
        choices = c("cell_type", "feats")
    )

    ICFscore <- copy(icfObject[["ICFscores"]])

    # other parameters
    direction <- match.arg(direction, choices = c("both", "up", "down"))


    ## sequential filter steps ##
    # 1. minimum number of source and target cells
    selection_scores <- ICFscore[nr_select >= min_cells &
        int_nr_select >= min_int_cells]

    # 2. create z-scores for log2fc per cell type
    selection_scores[, zscores := scale(log2fc), by = c(zscores_column)]

    # 3. filter based on z-scores and minimum levels
    comb_DT <- rbind(
        selection_scores[zscores >= min_zscore &
            abs(diff) >= min_spat_diff &
            log2fc >= min_log2_fc & sel >= min_cells_expr],
        selection_scores[zscores <= -min_zscore &
            abs(diff) >= min_spat_diff &
            log2fc <= -min_log2_fc &
            other >= min_int_cells_expr]
    )

    # 4. filter based on adjusted p-value (fdr)
    comb_DT <- comb_DT[p.adj < min_fdr]


    if (direction == "both") {
        selection_scores <- selection_scores
    } else if (direction == "up") {
        selection_scores <- selection_scores[log2fc >= min_log2_fc]
    } else if (direction == "down") {
        selection_scores <- selection_scores[log2fc <= -min_log2_fc]
    }


    newobj <- copy(icfObject)
    newobj[["ICFscores"]] <- comb_DT

    return(newobj)
}

#' @rdname filterInteractionChangedFeats
#' @export
filterICF <- filterInteractionChangedFeats








# * ####
# FTF feat-to-feat (pairs of ICF) ####

#' @title Combine ICF scores per interaction
#' @name .combineInteractionChangedFeatures_per_interaction
#' @description Combine ICF scores per interaction
#' @returns data.table
#' @keywords internal
.combineInteractionChangedFeatures_per_interaction <- function(icfObject,
    sel_int,
    selected_feats = NULL,
    specific_feats_1 = NULL,
    specific_feats_2 = NULL,
    min_cells = 5,
    min_int_cells = 3,
    min_fdr = 0.05,
    min_spat_diff = 0,
    min_log2_fc = 0.5) {
    # data.table variables
    unif_int <- feats <- cell_type <- p.adj <- nr_select <-
        int_nr_select <- log2fc <- sel <- NULL
    other <- p.value <- perm_sel <- perm_other <- perm_log2fc <-
        perm_diff <- NULL
    int_cell_type <- nr_other <- feats_combo <- feats_1 <- feats_2 <-
        type_int <- NULL

    if (!"icfObject" %in% class(icfObject)) {
        stop("icfObject needs to be the output from
            findInteractionChangedFeats() or findICF()")
    }

    ICFscore <- copy(icfObject[["ICFscores"]])

    test_used <- icfObject[["test_info"]][["test"]]

    ## subset on selected interaction
    subset <- ICFscore[unif_int == sel_int]

    ##  type of interactions
    type_interaction <- unique(subset[["type_int"]])

    ## first filtering ICFscores on feats
    if ((!is.null(specific_feats_1) & !is.null(specific_feats_2))) {
        if (length(specific_feats_1) != length(specific_feats_2)) {
            stop("specific_feats_1 must have the same length as
                specific_feats_2")
        }
        subset <- subset[feats %in% c(specific_feats_1, specific_feats_2)]
    } else if (!is.null(selected_feats)) {
        subset <- subset[feats %in% c(selected_feats)]
    }

    ## find number of unique cell types
    unique_cell_types <- unique(c(subset$cell_type, subset$int_cell_type))

    if (length(unique_cell_types) == 2) {
        ## CELL TYPE 1
        subset_cell_1 <- subset[cell_type == unique_cell_types[1]]

        if (nrow(subset_cell_1) == 0) {
            if (test_used == "permutation") {
                subset_cell_1 <- data.table::data.table(
                    "feats_1" = subset[["feats"]],
                    "sel_1" = NA,
                    "other_1" = NA,
                    "log2fc_1" = NA,
                    "diff_1" = NA,
                    "p.value_1" = NA,
                    "p.adj_1" = NA,
                    "perm_sel_1" = NA,
                    "perm_other_1" = NA,
                    "perm_log2fc_1" = NA,
                    "perm_diff_1" = NA,
                    "cell_type_1" = unique_cell_types[1],
                    "int_cell_type_1" = unique_cell_types[2],
                    "nr_select_1" = NA,
                    "nr_other_1" = NA,
                    "unif_int" = subset[["unif_int"]]
                )
            } else {
                subset_cell_1 <- data.table::data.table(
                    "feats_1" = subset[["feats"]],
                    "sel_1" = NA,
                    "other_1" = NA,
                    "log2fc_1" = NA,
                    "diff_1" = NA,
                    "p.value_1" = NA,
                    "p.adj_1" = NA,
                    "cell_type_1" = unique_cell_types[1],
                    "int_cell_type_1" = unique_cell_types[2],
                    "nr_select_1" = NA,
                    "nr_other_1" = NA,
                    "unif_int" = subset[["unif_int"]]
                )
            }
        } else {
            # filter on statistics
            subset_cell_1 <- subset_cell_1[p.adj <= min_fdr]
            subset_cell_1 <- subset_cell_1[nr_select >= min_cells]
            subset_cell_1 <- subset_cell_1[int_nr_select >= min_int_cells]
            subset_cell_1 <- subset_cell_1[abs(log2fc) >= min_log2_fc]
            subset_cell_1 <- subset_cell_1[abs(diff) >= min_spat_diff]

            if (test_used == "permutation") {
                # make it specific
                subset_cell_1 <- subset_cell_1[, .(
                    feats, sel, other, log2fc, diff, p.value, p.adj,
                    perm_sel, perm_other, perm_log2fc, perm_diff,
                    cell_type, int_cell_type, nr_select, nr_other, unif_int
                )]
                data.table::setnames(subset_cell_1,
                    old = c(
                        "feats", "sel", "other", "log2fc", "diff", "p.value",
                        "p.adj",
                        "perm_sel", "perm_other", "perm_log2fc", "perm_diff",
                        "cell_type", "int_cell_type", "nr_select", "nr_other"
                    ),
                    new = c(
                        "feats_1", "sel_1", "other_1", "log2fc_1", "diff_1",
                        "p.value_1", "p.adj_1",
                        "perm_sel_1", "perm_other_1", "perm_log2fc_1",
                        "perm_diff_1",
                        "cell_type_1", "int_cell_type_1", "nr_select_1",
                        "nr_other_1"
                    )
                )
            } else {
                # make it specific
                subset_cell_1 <- subset_cell_1[
                    , .(
                        feats, sel, other, log2fc, diff, p.value, p.adj,
                        cell_type, int_cell_type, nr_select, nr_other,
                        unif_int
                    )
                ]
                data.table::setnames(subset_cell_1,
                    old = c(
                        "feats", "sel", "other", "log2fc", "diff",
                        "p.value", "p.adj", "cell_type", "int_cell_type",
                        "nr_select", "nr_other"
                    ),
                    new = c(
                        "feats_1", "sel_1", "other_1", "log2fc_1",
                        "diff_1", "p.value_1", "p.adj_1", "cell_type_1",
                        "int_cell_type_1", "nr_select_1", "nr_other_1"
                    )
                )
            }
        }




        ## CELL TYPE 2
        subset_cell_2 <- subset[cell_type == unique_cell_types[2]]

        if (nrow(subset_cell_2) == 0) {
            if (test_used == "permutation") {
                subset_cell_2 <- data.table::data.table(
                    "feats_2" = subset[["feats"]],
                    "sel_2" = NA,
                    "other_2" = NA,
                    "log2fc_2" = NA,
                    "diff_2" = NA,
                    "p.value_2" = NA,
                    "p.adj_2" = NA,
                    "perm_sel_2" = NA,
                    "perm_other_2" = NA,
                    "perm_log2fc_2" = NA,
                    "perm_diff_2" = NA,
                    "cell_type_2" = unique_cell_types[2],
                    "int_cell_type_2" = unique_cell_types[1],
                    "nr_select_2" = NA,
                    "nr_other_2" = NA,
                    "unif_int" = subset[["unif_int"]]
                )
            } else {
                subset_cell_2 <- data.table::data.table(
                    "feats_2" = subset[["feats"]],
                    "sel_2" = NA,
                    "other_2" = NA,
                    "log2fc_2" = NA,
                    "diff_2" = NA,
                    "p.value_2" = NA,
                    "p.adj_2" = NA,
                    "cell_type_2" = unique_cell_types[2],
                    "int_cell_type_2" = unique_cell_types[1],
                    "nr_select_2" = NA,
                    "nr_other_2" = NA,
                    "unif_int" = subset[["unif_int"]]
                )
            }
        } else {
            # filter on statistics
            subset_cell_2 <- subset_cell_2[p.adj <= min_fdr]
            subset_cell_2 <- subset_cell_2[nr_select >= min_cells]
            subset_cell_2 <- subset_cell_2[int_nr_select >= min_int_cells]
            subset_cell_2 <- subset_cell_2[abs(log2fc) >= min_log2_fc]
            subset_cell_2 <- subset_cell_2[abs(diff) >= min_spat_diff]

            if (test_used == "permutation") {
                subset_cell_2 <- subset_cell_2[, .(
                    feats, sel, other, log2fc, diff, p.value, p.adj,
                    perm_sel, perm_other, perm_log2fc, perm_diff,
                    cell_type, int_cell_type, nr_select, nr_other, unif_int
                )]
                data.table::setnames(subset_cell_2,
                    old = c(
                        "feats", "sel", "other", "log2fc", "diff", "p.value",
                        "p.adj",
                        "perm_sel", "perm_other", "perm_log2fc", "perm_diff",
                        "cell_type", "int_cell_type", "nr_select", "nr_other"
                    ),
                    new = c(
                        "feats_2", "sel_2", "other_2", "log2fc_2", "diff_2",
                        "p.value_2", "p.adj_2",
                        "perm_sel_2", "perm_other_2", "perm_log2fc_2",
                        "perm_diff_2",
                        "cell_type_2", "int_cell_type_2", "nr_select_2",
                        "nr_other_2"
                    )
                )
            } else {
                subset_cell_2 <- subset_cell_2[
                    , .(
                        feats, sel, other, log2fc, diff, p.value, p.adj,
                        cell_type, int_cell_type, nr_select, nr_other,
                        unif_int
                    )
                ]
                data.table::setnames(subset_cell_2,
                    old = c(
                        "feats", "sel", "other", "log2fc", "diff",
                        "p.value", "p.adj", "cell_type", "int_cell_type",
                        "nr_select", "nr_other"
                    ),
                    new = c(
                        "feats_2", "sel_2", "other_2", "log2fc_2",
                        "diff_2", "p.value_2", "p.adj_2", "cell_type_2",
                        "int_cell_type_2", "nr_select_2", "nr_other_2"
                    )
                )
            }
        }

        merge_subsets <- data.table::merge.data.table(
            subset_cell_1, subset_cell_2,
            by = c("unif_int"),
            allow.cartesian = TRUE
        )
    } else if (length(unique_cell_types) == 1) {
        ## CELL TYPE 1
        subset_cell_1 <- subset[cell_type == unique_cell_types[1]]

        # filter on statistics
        subset_cell_1 <- subset_cell_1[p.adj <= min_fdr]
        subset_cell_1 <- subset_cell_1[nr_select >= min_cells]
        subset_cell_1 <- subset_cell_1[int_nr_select >= min_int_cells]
        subset_cell_1 <- subset_cell_1[abs(log2fc) >= min_log2_fc]
        subset_cell_1 <- subset_cell_1[abs(diff) >= min_spat_diff]

        # make it specific

        if (test_used == "permutation") {
            subset_cell_1A <- subset_cell_1[, .(
                feats, sel, other, log2fc, diff, p.value, p.adj,
                perm_sel, perm_other, perm_log2fc, perm_diff,
                cell_type, int_cell_type, nr_select, nr_other, unif_int
            )]
            data.table::setnames(subset_cell_1A,
                old = c(
                    "feats", "sel", "other", "log2fc", "diff", "p.value",
                    "p.adj",
                    "perm_sel", "perm_other", "perm_log2fc", "perm_diff",
                    "cell_type", "int_cell_type", "nr_select", "nr_other"
                ),
                new = c(
                    "feats_1", "sel_1", "other_1", "log2fc_1", "diff_1",
                    "p.value_1", "p.adj_1",
                    "perm_sel_1", "perm_other_1", "perm_log2fc_1",
                    "perm_diff_1",
                    "cell_type_1", "int_cell_type_1", "nr_select_1",
                    "nr_other_1"
                )
            )
        } else {
            subset_cell_1A <- subset_cell_1[
                , .(
                    feats, sel, other, log2fc, diff, p.value, p.adj,
                    cell_type, int_cell_type, nr_select, nr_other, unif_int
                )
            ]
            data.table::setnames(subset_cell_1A,
                old = c(
                    "feats", "sel", "other", "log2fc", "diff", "p.value",
                    "p.adj", "cell_type", "int_cell_type", "nr_select",
                    "nr_other"
                ),
                new = c(
                    "feats_1", "sel_1", "other_1", "log2fc_1", "diff_1",
                    "p.value_1", "p.adj_1", "cell_type_1",
                    "int_cell_type_1", "nr_select_1", "nr_other_1"
                )
            )
        }


        ## CELL TYPE 2

        if (test_used == "permutation") {
            subset_cell_1B <- subset_cell_1[, .(
                feats, sel, other, log2fc, diff, p.value, p.adj,
                perm_sel, perm_other, perm_log2fc, perm_diff,
                cell_type, int_cell_type, nr_select, nr_other, unif_int
            )]
            data.table::setnames(subset_cell_1B,
                old = c(
                    "feats", "sel", "other", "log2fc", "diff", "p.value",
                    "p.adj",
                    "perm_sel", "perm_other", "perm_log2fc", "perm_diff",
                    "cell_type", "int_cell_type", "nr_select", "nr_other"
                ),
                new = c(
                    "feats_2", "sel_2", "other_2", "log2fc_2", "diff_2",
                    "p.value_2", "p.adj_2",
                    "perm_sel_2", "perm_other_2", "perm_log2fc_2",
                    "perm_diff_2",
                    "cell_type_2", "int_cell_type_2", "nr_select_2",
                    "nr_other_2"
                )
            )
        } else {
            subset_cell_1B <- subset_cell_1[
                , .(
                    feats, sel, other, log2fc, diff, p.value, p.adj,
                    cell_type, int_cell_type, nr_select, nr_other, unif_int
                )
            ]
            data.table::setnames(subset_cell_1B,
                old = c(
                    "feats", "sel", "other", "log2fc", "diff", "p.value",
                    "p.adj", "cell_type", "int_cell_type", "nr_select",
                    "nr_other"
                ),
                new = c(
                    "feats_2", "sel_2", "other_2", "log2fc_2", "diff_2",
                    "p.value_2", "p.adj_2", "cell_type_2",
                    "int_cell_type_2", "nr_select_2", "nr_other_2"
                )
            )
        }

        merge_subsets <- data.table::merge.data.table(
            subset_cell_1A, subset_cell_1B,
            by = c("unif_int"),
            allow.cartesian = TRUE
        )
    }

    # restrict to feature combinations if needed
    if ((!is.null(specific_feats_1) & !is.null(specific_feats_2))) {
        merge_subsets[, feats_combo := paste0(feats_1, "--", feats_2)]
        all_combos <- c(
            paste0(specific_feats_1, "--", specific_feats_2),
            paste0(specific_feats_2, "--", specific_feats_1)
        )
        merge_subsets <- merge_subsets[feats_combo %in% all_combos]
        merge_subsets[, feats_combo := NULL]
    }

    merge_subsets[, type_int := type_interaction]
    return(merge_subsets)
}


#' @title combineInteractionChangedFeats
#' @name combineInteractionChangedFeats
#' @aliases combineICF
#' @description Combine ICF scores in a pairwise manner.
#' @param icfObject ICF (interaction changed feat) score object
#' @param selected_ints subset of selected cell-cell interactions (optional)
#' @param selected_feats subset of selected Features (optional)
#' @param specific_feats_1 specific Featureset combo
#' (need to position match specific_feats_2)
#' @param specific_feats_2 specific Featureset combo
#' (need to position match specific_feats_1)
#' @param min_cells minimum number of target cell type
#' @param min_int_cells minimum number of interacting cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum absolute log2 fold-change
#' @param do_parallel run calculations in parallel with mclapply
#' @param verbose verbose
#' @returns `combIcfObject` that contains the filtered differential feature
#' scores
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g_icf <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"),
#'     nr_permutations = 10
#' )
#'
#' cicf <- combineInteractionChangedFeats(g_icf)
#' force(cicf)
#' combineICF(g_icf) # this is a shortened alias
#' @export
combineInteractionChangedFeats <- function(icfObject,
    selected_ints = NULL,
    selected_feats = NULL,
    specific_feats_1 = NULL,
    specific_feats_2 = NULL,
    min_cells = 5,
    min_int_cells = 3,
    min_fdr = 0.05,
    min_spat_diff = 0,
    min_log2_fc = 0.5,
    do_parallel = TRUE,
    verbose = TRUE) {
    # NSE vars
    unif_int <- feat1_feat2 <- feats_1 <- feats_2 <- comb_logfc <-
        log2fc_1 <- log2fc_2 <- direction <- NULL

    ## check validity
    if (!"icfObject" %in% class(icfObject)) {
        stop("icfObject needs to be the output from
            findInteractionChangedFeats() or findICF()")
    }
    ICFscore <- copy(icfObject[["ICFscores"]])

    if (!is.null(selected_ints)) {
        ICFscore <- ICFscore[unif_int %in% selected_ints]
    }

    all_ints <- unique(ICFscore[["unif_int"]])

    # parallel
    if (do_parallel == TRUE) {
        FTFresults <- lapply_flex(X = all_ints, FUN = function(x) {
            tempres <- .combineInteractionChangedFeatures_per_interaction(
                icfObject = icfObject,
                sel_int = x,
                selected_feats = selected_feats,
                specific_feats_1 = specific_feats_1,
                specific_feats_2 = specific_feats_2,
                min_cells = min_cells,
                min_int_cells = min_int_cells,
                min_fdr = min_fdr,
                min_spat_diff = min_spat_diff,
                min_log2_fc = min_log2_fc
            )
        })
    } else {
        # for loop
        FTFresults <- list()

        for (i in seq_along(all_ints)) {
            x <- all_ints[[i]]

            if (verbose == TRUE) print(x)

            tempres <- .combineInteractionChangedFeatures_per_interaction(
                icfObject = icfObject,
                sel_int = x,
                selected_feats = selected_feats,
                specific_feats_1 = specific_feats_1,
                specific_feats_2 = specific_feats_2,
                min_cells = min_cells,
                min_int_cells = min_int_cells,
                min_fdr = min_fdr,
                min_spat_diff = min_spat_diff,
                min_log2_fc = min_log2_fc
            )
            FTFresults[[i]] <- tempres
        }
    }

    final_results <- do.call("rbind", FTFresults)

    final_results[, feat1_feat2 := paste0(feats_1, "--", feats_2)]

    final_results <- dt_sort_combine_two_columns(final_results,
        column1 = "feats_1", column2 = "feats_2",
        myname = "unif_feat_feat"
    )

    final_results[, comb_logfc := abs(log2fc_1) + abs(log2fc_2)]
    setorder(final_results, -comb_logfc)
    final_results[, direction := ifelse(log2fc_1 > 0 & log2fc_2 > 0, "both_up",
        ifelse(log2fc_1 < 0 & log2fc_2 < 0, "both_down", "mixed")
    )]

    combIcfObject <- structure(
        .Data = list(
            combICFscores = final_results,
            Giotto_info = list(
                "values" = icfObject[["Giotto_info"]][["values"]],
                "cluster" = icfObject[["Giotto_info"]][["cluster"]],
                "spatial network" =
                    icfObject[["Giotto_info"]][["spatial network"]]
            ),
            test_info = list(
                "test" = icfObject[["test_info"]][["test"]],
                "p.adj" = icfObject[["test_info"]][["p.adj"]],
                "min cells" = icfObject[["test_info"]][["min cells"]],
                "min interacting cells" = icfObject[["test_info"]][[
                    "min interacting cells"
                ]],
                "exclude selected cells" = icfObject[["test_info"]][[
                    "exclude selected cells"
                ]],
                "perm" = icfObject[["test_info"]][["perm"]]
            )
        ),
        class = "combIcfObject"
    )
    return(combIcfObject)
}

#' @rdname combineInteractionChangedFeats
#' @export
combineICF <- combineInteractionChangedFeats

#' @name print.combIcfObject
#' @title combIcfObject print method
#' @param x object to print
#' @param \dots additional params to pass (none implemented)
#' @keywords internal
#' @returns combIcfObject
#' @export
print.combIcfObject <- function(x, ...) {
    cat("An object of class", class(x), "\n")
    info <- list(
        dimensions = sprintf(
            "%d, %d (icf pairs, attributes)",
            nrow(x$combICFscores), ncol(x$combICFscores)
        )
    )
    print_list(info, pre = " -")
    cat("<giotto info>\n")
    print_list(x$Giotto_info, pre = " -")
    cat("<test info>\n")
    print_list(x$test_info, pre = " -")
}



# * ####
# cell communication ####


#' @title average_feat_feat_expression_in_groups
#' @name .average_feat_feat_expression_in_groups
#' @description calculate average expression per cluster
#' @param gobject giotto object to use
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param cluster_column cluster column with cell type information
#' @param feat_set_1 first specific feat set from feat pairs
#' @param feat_set_2 second specific feat set from feat pairs
#' @returns data.table with average expression scores for each cluster
#' @keywords internal
.average_feat_feat_expression_in_groups <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    cluster_column = "cell_types",
    feat_set_1,
    feat_set_2) {
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

    average_DT <- create_average_DT(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        meta_data_name = cluster_column
    )

    # change column names back to original
    new_colnames <- gsub(
        pattern = "cluster_", replacement = "", colnames(average_DT)
    )
    colnames(average_DT) <- new_colnames

    # keep order of colnames
    colnames_order <- sort(new_colnames)

    # feat_set_1 and feat_set_2 need to have same length and all feats need
    # to be present in data
    if (length(feat_set_1) != length(feat_set_2)) {
        stop("length of set1 needs to be the same as that of set2")
    }

    if (!all(c(feat_set_1, feat_set_2) %in% rownames(average_DT) == TRUE)) {
        stop("all selected feats from set 1 and 2 need to be present")
    }

    LR_pairs <- paste0(feat_set_1, "-", feat_set_2)

    # get ligand and receptor information
    ligand_match <- average_DT[
        match(feat_set_1, rownames(average_DT)), ,
        drop = FALSE
    ]
    receptor_match <- average_DT[
        match(feat_set_2, rownames(average_DT)), ,
        drop = FALSE
    ]

    # data.table variables
    ligand <- LR_comb <- receptor <- LR_expr <- lig_expr <- rec_expr <-
        lig_cell_type <- rec_cell_type <- NULL

    all_ligand_cols <- colnames(ligand_match)
    lig_test <- data.table::melt(
        data.table::as.data.table(ligand_match),
        measure.vars = all_ligand_cols
    )
    lig_test[, ligand := rep(rownames(ligand_match), ncol(ligand_match))]
    lig_test[, ligand := strsplit(ligand, "\\.")[[1]][1],
        by = seq_len(nrow(lig_test))
    ]
    lig_test[, LR_comb := rep(LR_pairs, ncol(ligand_match))]
    setnames(lig_test, "value", "lig_expr")
    setnames(lig_test, "variable", "lig_cell_type")

    all_receptor_cols <- colnames(receptor_match)
    rec_test <- data.table::melt(
        data.table::as.data.table(receptor_match),
        measure.vars = all_receptor_cols
    )
    rec_test[, receptor := rep(rownames(receptor_match), ncol(receptor_match))]
    rec_test[, receptor := strsplit(
        receptor, "\\."
    )[[1]][1], by = seq_len(nrow(rec_test))]
    rec_test[, LR_comb := rep(LR_pairs, ncol(receptor_match))]
    setnames(rec_test, "value", "rec_expr")
    setnames(rec_test, "variable", "rec_cell_type")

    lig_rec_test <- merge(
        lig_test, rec_test,
        by = "LR_comb", allow.cartesian = TRUE
    )
    lig_rec_test[, LR_expr := lig_expr + rec_expr]


    lig_rec_test[, lig_cell_type := factor(
        lig_cell_type,
        levels = colnames_order
    )]
    lig_rec_test[, rec_cell_type := factor(
        rec_cell_type,
        levels = colnames_order
    )]
    setorder(lig_rec_test, LR_comb, lig_cell_type, rec_cell_type)

    return(lig_rec_test)
}


#' @title exprCellCellcom
#' @name exprCellCellcom
#' @description Cell-Cell communication scores based on expression only
#' @param gobject giotto object to use
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param cluster_column cluster column with cell type information
#' @param random_iter number of iterations
#' @param feat_set_1 first specific feature set from feature pairs
#' @param feat_set_2 second specific feature set from feature pairs
#' @param log2FC_addendum addendum to add when calculating log2FC
#' @param detailed provide more detailed information
#' (random variance and z-score)
#' @param adjust_method which method to adjust p-values
#' @param adjust_target adjust multiple hypotheses at the cell or feature level
#' @param set_seed `logical`. set seed for random simulations (default = TRUE)
#' @param seed_number seed number
#' @param verbose verbose
#' @returns Cell-Cell communication scores for feature pairs based on
#' expression only
#' @details Statistical framework to identify if pairs of features
#' (such as ligand-receptor combinations)
#' are expressed at higher levels than expected based on a reshuffled null
#' distribution of feature expression values,
#' without considering the spatial position of cells.
#' More details will follow soon.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' res <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935",
#'     feat_set_2 = "9630013A20Rik"
#' )
#'
#' force(res)
#' @export
exprCellCellcom <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    cluster_column = "cell_types",
    random_iter = 1000,
    feat_set_1,
    feat_set_2,
    log2FC_addendum = 0.1,
    detailed = FALSE,
    adjust_method = c(
        "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "none"
    ),
    adjust_target = c("feats", "cells"),
    set_seed = TRUE,
    seed_number = 1234,
    verbose = TRUE) {
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

    # data.table variables
    lig_nr <- lig_cell_type <- rec_nr <- rec_cell_type <- rand_expr <-
        av_diff <- log2fc <- LR_expr <- pvalue <- NULL
    LR_cell_comb <- p.adj <- LR_comb <- PI <- sd_diff <- z_score <- NULL

    # get parameters
    adjust_method <- match.arg(adjust_method, choices = c(
        "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "none"
    ))
    adjust_target <- match.arg(adjust_target, choices = c("feats", "cells"))

    # get information about number of cells
    cell_metadata <- pDataDT(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit
    )
    nr_cell_types <- cell_metadata[, .N, by = c(cluster_column)]
    nr_cells <- nr_cell_types$N
    names(nr_cells) <- nr_cell_types$cluster_column


    comScore <- .average_feat_feat_expression_in_groups(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        cluster_column = cluster_column,
        feat_set_1 = feat_set_1,
        feat_set_2 = feat_set_2
    )

    comScore[, lig_nr := nr_cells[lig_cell_type]]
    comScore[, rec_nr := nr_cells[rec_cell_type]]

    # prepare for randomized scores
    total_av <- rep(0, nrow(comScore))

    if (detailed == FALSE) {
        total_sum <- rep(0, nrow(comScore))
    } else {
        total_sum <- matrix(nrow = nrow(comScore), ncol = random_iter)
    }

    total_bool <- rep(0, nrow(comScore))


    ## parallel option ##
    # not yet available


    with_pbar({
        pb <- pbar(steps = random_iter)

        for (sim in seq_len(random_iter)) {
            # create temporary giotto
            tempGiotto <- subsetGiotto(
                gobject = gobject,
                feat_type = feat_type,
                spat_unit = spat_unit
            )

            # randomize annoation
            cell_types <- cell_metadata[[cluster_column]]
            if (set_seed == TRUE) {
                seed_number <- seed_number + sim
                set.seed(seed = seed_number)
            }
            random_cell_types <- sample(
                x = cell_types, size = length(cell_types)
            )
            tempGiotto <- addCellMetadata(
                gobject = tempGiotto,
                feat_type = feat_type,
                spat_unit = spat_unit,
                new_metadata = random_cell_types,
                by_column = FALSE # on purpose since values are random
            )

            # get random communication scores
            randomScore <- .average_feat_feat_expression_in_groups(
                gobject = tempGiotto,
                feat_type = feat_type,
                spat_unit = spat_unit,
                cluster_column = "random_cell_types",
                feat_set_1 = feat_set_1,
                feat_set_2 = feat_set_2
            )

            # average random score
            total_av <- total_av + randomScore[["LR_expr"]]

            # difference between observed and random
            difference <- comScore[["LR_expr"]] - randomScore[["LR_expr"]]

            # calculate total difference
            if (detailed == FALSE) {
                total_sum <- total_sum + difference
            } else {
                total_sum[, sim] <- difference
            }

            # calculate p-values
            difference[difference > 0] <- 1
            difference[difference < 0] <- -1
            total_bool <- total_bool + difference

            pb(sprintf("simulation %d", sim))
        }
    })


    comScore[, rand_expr := total_av / random_iter]

    if (detailed == TRUE) {
        av_difference_scores <- rowMeans_flex(total_sum)
        sd_difference_scores <- apply(total_sum, MARGIN = 1, FUN = stats::sd)

        comScore[, av_diff := av_difference_scores]
        comScore[, sd_diff := sd_difference_scores]
        comScore[, z_score := (LR_expr - rand_expr) / sd_diff]
    } else {
        comScore[, av_diff := total_sum / random_iter]
    }

    comScore[, log2fc := log2((LR_expr + log2FC_addendum) / (
        rand_expr + log2FC_addendum))]
    comScore[, pvalue := total_bool / random_iter]
    comScore[, pvalue := ifelse(pvalue > 0, 1 - pvalue, 1 + pvalue)]
    comScore[, LR_cell_comb := paste0(lig_cell_type, "--", rec_cell_type)]

    if (adjust_target == "feats") {
        comScore[, p.adj := stats::p.adjust(pvalue, method = adjust_method),
            by = .(LR_cell_comb)
        ]
    } else if (adjust_target == "cells") {
        comScore[, p.adj := stats::p.adjust(pvalue, method = adjust_method),
            by = .(LR_comb)
        ]
    }


    # get minimum adjusted p.value that is not zero
    all_p.adj <- comScore[["p.adj"]]
    lowest_p.adj <- min(all_p.adj[all_p.adj != 0])
    comScore[, PI := ifelse(p.adj == 0, log2fc * (-log10(lowest_p.adj)),
        log2fc * (-log10(p.adj))
    )]

    data.table::setorder(comScore, LR_comb, -LR_expr)

    return(comScore)
}



#' @title Create randomized cell IDs within a selection of cell types
#' @name .create_cell_type_random_cell_IDs
#' @description creates randomized cell ids within a selection of cell types
#' @param gobject giotto object to use
#' @param feat_type feature type
#' @param cluster_column cluster column with cell type information
#' @param needed_cell_types vector of cell type names for which a random id
#' will be found
#' @param set_seed set a seed for reproducibility
#' @param seed_number seed number
#' @returns list of randomly sampled cell ids with same cell type composition
#' @keywords internal
.create_cell_type_random_cell_IDs <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    cluster_column = "cell_types",
    needed_cell_types,
    set_seed = FALSE,
    seed_number = 1234) {
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

    # subset metadata to choose from
    full_metadata <- pDataDT(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit
    )
    possible_metadata <- full_metadata[get(cluster_column) %in% unique(
        needed_cell_types
    )]

    sample_ids <- list()

    uniq_types <- unique(needed_cell_types)

    for (i in seq_along(uniq_types)) {
        uniq_type <- uniq_types[i]
        length_random <- length(needed_cell_types[
            needed_cell_types == uniq_type
        ])
        if (set_seed == TRUE) {
            set.seed(seed = seed_number)
        }
        sub_sample_ids <- possible_metadata[get(cluster_column) == uniq_type][
            sample(x = seq_len(.N), size = length_random)
        ][["cell_ID"]]
        sample_ids[[i]] <- sub_sample_ids
    }
    return(unlist(sample_ids))
}




#' @title Spatial cell cell communication scoring
#' @name spatCellCellcom
#' @description Spatial Cell-Cell communication scores based on spatial
#' expression of interacting cells
#' @param gobject giotto object to use
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spatial_network_name spatial network to use for identifying
#' interacting cells
#' @param cluster_column cluster column with cell type information
#' @param random_iter number of iterations
#' @param feat_set_1 first specific feature set from feature pairs
#' @param feat_set_2 second specific feature set from feature pairs
#' @param gene_set_1 deprecated, use feat_set_1
#' @param gene_set_2 deprecated, use feat_set_2
#' @param log2FC_addendum addendum to add when calculating log2FC
#' @param min_observations minimum number of interactions needed to be
#' considered
#' @param detailed provide more detailed information
#' (random variance and z-score)
#' @param adjust_method which method to adjust p-values
#' @param adjust_target adjust multiple hypotheses at the cell or feature level
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param set_seed set a seed for reproducibility
#' @param seed_number seed number
#' @param verbose verbose
#' @returns Cell-Cell communication scores for feature pairs based on spatial
#' interaction
#' @details Statistical framework to identify if pairs of genes
#' (such as ligand-receptor combinations)
#' are expressed at higher levels than expected based on a reshuffled null
#' distribution of feature expression values in cells that are spatially in
#' proximity to each other.
#' * **LR_comb:** Pair of ligand and receptor
#' * **lig_cell_type:** cell type to assess expression level of ligand
#' * **lig_expr:** average expression of ligand in lig_cell_type
#' * **ligand:** ligand name
#' * **rec_cell_type:** cell type to assess expression level of receptor
#' * **rec_expr:** average expression of receptor in rec_cell_type
#' * **receptor:** receptor name
#' * **LR_expr:** combined average ligand and receptor expression
#' * **lig_nr:** total number of cells from lig_cell_type that spatially
#' interact with cells from rec_cell_type
#' * **rec_nr:** total number of cells from rec_cell_type that spatially
#' interact with cells from lig_cell_type
#' * **rand_expr:** average combined ligand and receptor expression from
#' random spatial permutations
#' * **av_diff:** average difference between LR_expr and rand_expr over all
#' random spatial permutations
#' * **sd_diff:** (optional) standard deviation of the difference between
#' LR_expr and rand_expr over all random spatial permutations
#' * **z_score:** (optional) z-score
#' * **log2fc:** log2 fold-change (LR_expr/rand_expr)
#' * **pvalue:** p-value
#' * **LR_cell_comb:** cell type pair combination
#' * **p.adj:** adjusted p-value
#' * **PI:** significance score: \eqn{log2fc * -log10(p.adj)}
#'
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' res1 <- spatCellCellcom(
#'     gobject = g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935",
#'     feat_set_2 = "9630013A20Rik",
#'     verbose = "a lot",
#'     random_iter = 10
#' )
#' force(res1)
#'
#' res2 <- specificCellCellcommunicationScores(g,
#'     cluster_column = "leiden_clus",
#'     cell_type_1 = 1,
#'     cell_type_2 = 2,
#'     feat_set_1 = "Gm19935",
#'     feat_set_2 = "9630013A20Rik"
#' )
#'
#' force(res2)
#' @export
spatCellCellcom <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    spatial_network_name = "Delaunay_network",
    cluster_column = NULL,
    random_iter = 1000,
    feat_set_1,
    feat_set_2,
    gene_set_1 = NULL,
    gene_set_2 = NULL,
    log2FC_addendum = 0.1,
    min_observations = 2,
    detailed = FALSE,
    adjust_method = c(
        "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "none"
    ),
    adjust_target = c("feats", "cells"),
    do_parallel = TRUE,
    cores = NA,
    set_seed = TRUE,
    seed_number = 1234,
    verbose = c("a little", "a lot", "none")) {
    verbose <- match.arg(verbose, choices = c("a little", "a lot", "none"))

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

    ## check if spatial network exists ##
    spat_networks <- list_spatial_networks_names(gobject,
        spat_unit = spat_unit
    )

    if (!spatial_network_name %in% spat_networks) {
        stop(
            spatial_network_name, " is not an existing spatial network \n",
            "use showNetworks() to see the available networks \n",
            "or create a new spatial network with createSpatialNetwork()"
        )
    }

    ## deprecated arguments
    if (!is.null(gene_set_1)) {
        feat_set_1 <- gene_set_1
        warning("gene_set_1 is deprecated, use feat_set_1 in the future")
    }
    if (!is.null(gene_set_2)) {
        feat_set_2 <- gene_set_2
        warning("gene_set_2 is deprecated, use feat_set_2 in the future")
    }

    if (is.null(cluster_column)) {
        stop("Name of column in cell metadata with cell type info is needed")
    }


    cell_metadata <- pDataDT(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit
    )

    ## get all combinations between cell types
    all_uniq_values <- unique(cell_metadata[[cluster_column]])
    same_DT <- data.table::data.table(
        V1 = all_uniq_values, V2 = all_uniq_values
    )
    combn_DT <- data.table::as.data.table(t(combn(all_uniq_values, m = 2)))
    combn_DT <- rbind(same_DT, combn_DT)

    ## parallel option ##
    if (isTRUE(do_parallel)) {
        savelist <- lapply_flex(
            X = seq_len(nrow(combn_DT)), future.seed = TRUE,
            cores = cores, fun = function(row) {
                cell_type_1 <- combn_DT[row][["V1"]]
                cell_type_2 <- combn_DT[row][["V2"]]

                specific_scores <- specificCellCellcommunicationScores(
                    gobject = gobject,
                    feat_type = feat_type,
                    spat_unit = spat_unit,
                    cluster_column = cluster_column,
                    random_iter = random_iter,
                    cell_type_1 = cell_type_1,
                    cell_type_2 = cell_type_2,
                    feat_set_1 = feat_set_1,
                    feat_set_2 = feat_set_2,
                    spatial_network_name = spatial_network_name,
                    log2FC_addendum = log2FC_addendum,
                    min_observations = min_observations,
                    detailed = detailed,
                    adjust_method = adjust_method,
                    adjust_target = adjust_target,
                    set_seed = set_seed,
                    seed_number = seed_number,
                    verbose = verbose %in% c("a lot")
                )
            }
        )
    } else {
        ## for loop over all combinations ##
        savelist <- list()
        countdown <- nrow(combn_DT)

        for (row in seq_len(nrow(combn_DT))) {
            cell_type_1 <- combn_DT[row][["V1"]]
            cell_type_2 <- combn_DT[row][["V2"]]

            if (verbose == "a little" || verbose == "a lot") {
                cat(sprintf(
                    "[PROCESS nr  %d :  %d  and  %d]  ",
                    countdown, cell_type_1, cell_type_2
                ))
            }

            if (verbose %in% c("a little", "none")) {
                specific_verbose <- FALSE
            } else {
                specific_verbose <- TRUE
            }

            specific_scores <- specificCellCellcommunicationScores(
                gobject = gobject,
                feat_type = feat_type,
                spat_unit = spat_unit,
                cluster_column = cluster_column,
                random_iter = random_iter,
                cell_type_1 = cell_type_1,
                cell_type_2 = cell_type_2,
                feat_set_1 = feat_set_1,
                feat_set_2 = feat_set_2,
                spatial_network_name = spatial_network_name,
                log2FC_addendum = log2FC_addendum,
                min_observations = min_observations,
                detailed = detailed,
                adjust_method = adjust_method,
                adjust_target = adjust_target,
                set_seed = set_seed,
                seed_number = seed_number,
                verbose = specific_verbose
            )
            savelist[[row]] <- specific_scores
            countdown <- countdown - 1
        }
    }

    finalDT <- do.call("rbind", savelist)

    # data.table variables
    LR_comb <- LR_expr <- NULL

    data.table::setorder(finalDT, LR_comb, -LR_expr)

    return(finalDT)
}





#' @rdname spatCellCellcom
#' @param cell_type_1 character. First cell type
#' @param cell_type_2 character. Second cell type
#' @export
specificCellCellcommunicationScores <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        spatial_network_name = "Delaunay_network",
        cluster_column = NULL,
        random_iter = 100,
        cell_type_1 = NULL,
        cell_type_2 = NULL,
        feat_set_1,
        feat_set_2,
        gene_set_1 = NULL,
        gene_set_2 = NULL,
        log2FC_addendum = 0.1,
        min_observations = 2,
        detailed = FALSE,
        adjust_method = c(
            "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel",
            "BY", "none"
        ),
        adjust_target = c("feats", "cells"),
        set_seed = FALSE,
        seed_number = 1234,
        verbose = TRUE) {
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

    ## deprecated arguments
    if (!is.null(gene_set_1)) {
        feat_set_1 <- gene_set_1
        warning("gene_set_1 is deprecated, use feat_set_1 in the future")
    }
    if (!is.null(gene_set_2)) {
        feat_set_2 <- gene_set_2
        warning("gene_set_2 is deprecated, use feat_set_2 in the future")
    }

    if (is.null(cluster_column)) {
        stop("Name of column in cell metadata with cell type info is needed")
    }

    if (is.null(cell_type_1) || is.null(cell_type_2)) {
        stop(sprintf(
            "`%s` and `%s` in `%s` must be given",
            "cell_type_1", "cell_type_2", "cluster_column"
        ))
    }


    # data.table variables
    from_to <- cell_ID <- lig_cell_type <- rec_cell_type <- lig_nr <-
        rec_nr <- rand_expr <- NULL
    av_diff <- log2fc <- LR_expr <- pvalue <- LR_cell_comb <- p.adj <-
        LR_comb <- PI <- NULL
    sd_diff <- z_score <- NULL

    # get parameters
    adjust_method <- match.arg(adjust_method, choices = c(
        "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel",
        "BY", "none"
    ))
    adjust_target <- match.arg(adjust_target, choices = c("feats", "cells"))

    # metadata
    cell_metadata <- pDataDT(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit
    )

    # get annotated spatial network
    annot_network <- annotateSpatialNetwork(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )

    cell_direction_1 <- paste0(cell_type_1, "-", cell_type_2)
    cell_direction_2 <- paste0(cell_type_2, "-", cell_type_1)

    subset_annot_network <- annot_network[from_to %in% c(
        cell_direction_1, cell_direction_2
    )]

    # make sure that there are sufficient observations
    if (nrow(subset_annot_network) <= min_observations) {
        return(NULL)
    } else {
        # subset giotto object to only interacting cells
        subset_ids <- unique(c(
            subset_annot_network$to, subset_annot_network$from
        ))
        subsetGiotto <- subsetGiotto(
            gobject = gobject,
            cell_ids = subset_ids,
            feat_type = feat_type,
            spat_unit = spat_unit
        )

        # get information about number of cells
        temp_meta <- pDataDT(subsetGiotto,
            feat_type = feat_type,
            spat_unit = spat_unit
        )
        nr_cell_types <- temp_meta[cell_ID %in% subset_ids][
            , .N,
            by = c(cluster_column)
        ]
        nr_cells <- nr_cell_types$N
        names(nr_cells) <- nr_cell_types$cell_types

        # get average communication scores
        comScore <- .average_feat_feat_expression_in_groups(
            gobject = subsetGiotto,
            feat_type = feat_type,
            spat_unit = spat_unit,
            cluster_column = cluster_column,
            feat_set_1 = feat_set_1,
            feat_set_2 = feat_set_2
        )
        comScore <- comScore[(lig_cell_type == cell_type_1 &
            rec_cell_type == cell_type_2) |
            (lig_cell_type == cell_type_2 & rec_cell_type == cell_type_1)]

        comScore[, lig_nr := nr_cells[lig_cell_type]]
        comScore[, rec_nr := nr_cells[rec_cell_type]]

        # prepare for randomized scores
        total_av <- rep(0, nrow(comScore))

        if (detailed == FALSE) {
            total_sum <- rep(0, nrow(comScore))
        } else {
            total_sum <- matrix(nrow = nrow(comScore), ncol = random_iter)
        }

        total_bool <- rep(0, nrow(comScore))

        # identify which cell types you need
        subset_metadata <- cell_metadata[cell_ID %in% subset_ids]
        needed_cell_types <- subset_metadata[[cluster_column]]

        if (verbose) cat("simulations:")

        ## simulations ##
        for (sim in seq_len(random_iter)) {
            if (verbose) {
                cat(sprintf(" %s ", sim))
            }

            # get random ids and subset
            if (set_seed == TRUE) {
                seed_number <- seed_number + sim
                set.seed(seed = seed_number)
            }
            random_ids <- .create_cell_type_random_cell_IDs(
                gobject = gobject,
                feat_type = feat_type,
                spat_unit = spat_unit,
                cluster_column = cluster_column,
                needed_cell_types = needed_cell_types,
                set_seed = set_seed,
                seed_number = seed_number
            )
            tempGiotto <- subsetGiotto(
                gobject = gobject,
                cell_ids = random_ids,
                feat_type = feat_type,
                spat_unit = spat_unit
            )

            # get random communication scores
            randomScore <- .average_feat_feat_expression_in_groups(
                gobject = tempGiotto,
                feat_type = feat_type,
                spat_unit = spat_unit,
                cluster_column = cluster_column,
                feat_set_1 = feat_set_1,
                feat_set_2 = feat_set_2
            )
            randomScore <- randomScore[(lig_cell_type == cell_type_1 &
                rec_cell_type == cell_type_2) |
                (lig_cell_type == cell_type_2 & rec_cell_type == cell_type_1)]




            # average random score
            total_av <- total_av + randomScore[["LR_expr"]]

            # difference between observed and random
            difference <- comScore[["LR_expr"]] - randomScore[["LR_expr"]]

            # calculate total difference
            if (detailed == FALSE) {
                total_sum <- total_sum + difference
            } else {
                total_sum[, sim] <- difference
            }

            # calculate p-values
            difference[difference > 0] <- 1
            difference[difference < 0] <- -1
            total_bool <- total_bool + difference
        }

        comScore[, rand_expr := total_av / random_iter]

        if (detailed == TRUE) {
            av_difference_scores <- rowMeans_flex(total_sum)
            sd_difference_scores <- apply(
                total_sum,
                MARGIN = 1, FUN = stats::sd
            )

            comScore[, av_diff := av_difference_scores]
            comScore[, sd_diff := sd_difference_scores]
            comScore[, z_score := (LR_expr - rand_expr) / sd_diff]
        } else {
            comScore[, av_diff := total_sum / random_iter]
        }


        comScore[, log2fc := log2((LR_expr + log2FC_addendum) / (
            rand_expr + log2FC_addendum))]
        comScore[, pvalue := total_bool / random_iter]
        comScore[, pvalue := ifelse(pvalue > 0, 1 - pvalue, 1 + pvalue)]
        comScore[, LR_cell_comb := paste0(lig_cell_type, "--", rec_cell_type)]

        if (adjust_target == "feats") {
            comScore[, p.adj := stats::p.adjust(
                pvalue,
                method = adjust_method
            ), by = .(LR_cell_comb)]
        } else if (adjust_target == "cells") {
            comScore[, p.adj := stats::p.adjust(
                pvalue,
                method = adjust_method
            ), by = .(LR_comb)]
        }

        # get minimum adjusted p.value that is not zero
        all_p.adj <- comScore[["p.adj"]]
        nonzero_p.adj <- all_p.adj[all_p.adj != 0]
        if (length(nonzero_p.adj) == 0L) {
            warning(
                call. = FALSE,
                "no adjusted p.values that are not zero; returning Inf"
            )
            if (verbose) cat("<- Inf returned")
            lowest_p.adj <- Inf
        } else {
            lowest_p.adj <- min(nonzero_p.adj)
        }

        comScore[, PI := ifelse(p.adj == 0, log2fc * (
            -log10(lowest_p.adj)), log2fc * (-log10(p.adj)))]

        if (verbose) cat("\n")

        return(comScore)
    }
}




#' @title Combine cell cell communication tables
#' @name combCCcom
#' @description Combine spatial and expression based cell-cell communication
#' data.tables
#' @param spatialCC spatial cell-cell communication scores
#' @param exprCC expression cell-cell communication scores
#' @param min_lig_nr minimum number of ligand cells
#' @param min_rec_nr minimum number of receptor cells
#' @param min_padj_value minimum adjusted p-value
#' @param min_log2fc minimum log2 fold-change
#' @param min_av_diff minimum average expression difference
#' @param detailed detailed option used with \code{\link{spatCellCellcom}}
#' (default = FALSE)
#' @returns combined data.table with spatial and expression communication data
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' exprCC <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935",
#'     feat_set_2 = "9630013A20Rik"
#' )
#'
#' spatialCC <- spatCellCellcom(
#'     gobject = g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935",
#'     feat_set_2 = "9630013A20Rik",
#'     verbose = "a lot",
#'     random_iter = 10
#' )
#'
#' combCC <- combCCcom(spatialCC = spatialCC, exprCC = exprCC)
#' force(combCC)
#' @export
combCCcom <- function(spatialCC,
    exprCC,
    min_lig_nr = 3,
    min_rec_nr = 3,
    min_padj_value = 1,
    min_log2fc = 0,
    min_av_diff = 0,
    detailed = FALSE) {
    # data.table variables
    lig_nr <- rec_nr <- p.adj <- log2fc <- av_diff <- NULL

    spatialCC <- spatialCC[lig_nr >= min_lig_nr & rec_nr >= min_rec_nr &
        p.adj <= min_padj_value & abs(log2fc) >= min_log2fc &
        abs(av_diff) >= min_av_diff]


    if (detailed == TRUE) {
        old_detailed <- c("sd_diff", "z_score")
        new_detailed <- c("sd_diff_spat", "z_score_spat")
    } else {
        old_detailed <- NULL
        new_detailed <- NULL
    }


    data.table::setnames(
        x = spatialCC,
        old = c(
            "lig_expr", "rec_expr", "LR_expr", "lig_nr", "rec_nr",
            "rand_expr", "av_diff", old_detailed,
            "log2fc", "pvalue", "p.adj", "PI"
        ),
        new = c(
            "lig_expr_spat", "rec_expr_spat", "LR_expr_spat", "lig_nr_spat",
            "rec_nr_spat",
            "rand_expr_spat", "av_diff_spat", new_detailed, "log2fc_spat",
            "pvalue_spat", "p.adj_spat", "PI_spat"
        )
    )

    merge_DT <- data.table::merge.data.table(spatialCC, exprCC, by = c(
        "LR_comb", "LR_cell_comb",
        "lig_cell_type", "rec_cell_type",
        "ligand", "receptor"
    ))

    # data.table variables
    LR_expr_rnk <- LR_expr <- LR_comb <- LR_spat_rnk <- LR_expr_spat <-
        exprPI_rnk <- PI <- spatPI_rnk <- PI_spat <- NULL

    # rank for expression levels
    merge_DT[, LR_expr_rnk := rank(-LR_expr), by = LR_comb]
    merge_DT[, LR_spat_rnk := rank(-LR_expr_spat), by = LR_comb]

    # rank for differential activity levels
    merge_DT[, exprPI_rnk := rank(-PI), by = LR_comb]
    merge_DT[, spatPI_rnk := rank(-PI_spat), by = LR_comb]

    return(merge_DT)
}







# DEPRECATED ####

#' @title deprecated
#' @name findInteractionChangedGenes
#' @description Identifies cell-to-cell Interaction Changed Features (ICF),
#' i.e. genes that are differentially expressed due to interactions with other
#' cell types.
#' @param ... params to pass to \code{findInteractionChangedFeats}
#' @seealso \code{\link{findInteractionChangedFeats}}
#' @returns interaction changed genes
#' @export
findInteractionChangedGenes <- function(...) {
    .Deprecated(new = "findInteractionChangedFeats")

    findInteractionChangedFeats(...)
}



#' @title deprecated
#' @name findCellProximityGenes
#' @description Identifies cell-to-cell Interaction Changed Features (ICF),
#' i.e. genes that are differentially expressed due to proximity to other cell
#' types.
#' @inheritDotParams findInteractionChangedFeats
#' @seealso \code{\link{findInteractionChangedFeats}}
#' @returns cell-cell interaction changed genes
#' @export
findCellProximityGenes <- function(...) {
    .Deprecated(new = "findInteractionChangedFeats")

    findInteractionChangedFeats(...)
}



#' @title deprecated
#' @name findICG
#' @description Identifies cell-to-cell Interaction Changed Features (ICF),
#' i.e. genes that are differentially expressed due to interaction with other
#' cell types.
#' @inheritDotParams findICF
#' @seealso \code{\link{findICF}}
#' @returns cell-cell interaction changed features
#' @export
findICG <- function(...) {
    .Deprecated(new = "findICF")

    findICF(...)
}



#' @title deprecated
#' @name findCPG
#' @description Identifies cell-to-cell Interaction Changed Features (ICF),
#' i.e. genes that are differentially expressed due to proximity to other cell
#' types.
#' @inheritDotParams findICF
#' @returns cell-to-cell Interaction Changed Genes
#' @seealso \code{\link{findICF}}
#' @export
findCPG <- function(...) {
    .Deprecated(new = "findICF")

    findICF(...)
}

#' @title deprecated
#' @name filterInteractionChangedGenes
#' @description Filter Interaction Changed Feature scores.
#' @inheritDotParams filterInteractionChangedFeats
#' @seealso \code{\link{filterInteractionChangedFeats}}
#' @returns filtered interaction changed feature scores
#' @export
filterInteractionChangedGenes <- function(...) {
    .Deprecated(new = "filterInteractionChangedFeats")

    filterInteractionChangedFeats(...)
}


#' @title deprecated
#' @name filterCellProximityGenes
#' @description Filter Interaction Changed Feature scores.
#' @inheritDotParams filterInteractionChangedFeats
#' @seealso \code{\link{filterInteractionChangedFeats}}
#' @returns proximity genes
#' @export
filterCellProximityGenes <- function(...) {
    .Deprecated(new = "filterInteractionChangedFeats")

    filterInteractionChangedFeats(...)
}





#' @title deprecated
#' @name filterICG
#' @description Filter Interaction Changed Gene scores.
#' @inheritDotParams filterICF
#' @seealso \code{\link{filterICF}}
#' @returns filtered interaction changed gene scores
#' @export
filterICG <- function(...) {
    .Deprecated(new = "filterICF")

    filterICF(...)
}



#' @title deprecated
#' @name filterCPG
#' @description Filter Interaction Changed Gene scores.
#' @inheritDotParams filterICF
#' @seealso \code{\link{filterICF}}
#' @returns filtered interaction changed gene scores
#' @export
filterCPG <- function(...) {
    .Deprecated(new = "filterICF")

    filterICF(...)
}


#' @title deprecated
#' @name combineInteractionChangedGenes
#' @description Combine ICF scores in a pairwise manner.
#' @inheritDotParams combineInteractionChangedFeats
#' @returns ICF scores
#' @seealso \code{\link{combineInteractionChangedFeats}}
#' @export
combineInteractionChangedGenes <- function(...) {
    .Deprecated(new = "combineInteractionChangedFeats")

    combineInteractionChangedFeats(...)
}


#' @title deprecated
#' @name combineCellProximityGenes
#' @description Combine ICF scores in a pairwise manner.
#' @inheritDotParams combineInteractionChangedFeats
#' @returns ICF scores
#' @seealso \code{\link{combineInteractionChangedFeats}}
#' @export
combineCellProximityGenes <- function(...) {
    .Deprecated(new = "combineInteractionChangedFeats")

    combineInteractionChangedFeats(...)
}

#' @title deprecated
#' @name combineICG
#' @description Combine ICF scores in a pairwise manner.
#' @inheritDotParams combineICF
#' @returns ICF scores
#' @seealso \code{\link{combineICF}}
#' @export
combineICG <- function(...) {
    .Deprecated(new = "combineICF")

    combineICF(...)
}

#' @title deprecated
#' @name combineCPG
#' @description Combine ICF scores in a pairwise manner.
#' @inheritDotParams combineICF
#' @returns ICF scores
#' @seealso \code{\link{combineICF}}
#' @export
combineCPG <- function(...) {
    .Deprecated(new = "combineICF")

    combineICF(...)
}

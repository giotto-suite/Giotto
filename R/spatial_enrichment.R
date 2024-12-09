## create spatial enrichment matrix ####

#' @title PAGE feature enrichment
#' @name enrichment_PAGE
#' @description
#' Expression feature-based enrichment scoring of labels. \cr
#' A binary matrix of signature features (e.g. for cell types or processes) can
#' either be directly provided or converted from a list using
#' `makeSignMatrixPAGE()`. This matrix is then used with `runPAGEEnrich()` in
#' order to calculate feature signature enrichment scores per spatial position
#' using PAGE.
#' @param sign_names `character` vector with names (labels) for each provided
#' feat signature
#' @param sign_list list of feats in signature
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param sign_matrix binary matrix of signature feats for each cell type /
#' process. Alternatively a list of signature feats can be provided to
#' `makeSignMatrixPAGE()`, which will create the matrix for you.
#' @param expression_values expression values to use
#' @param min_overlap_genes minimum number of overlapping feats in
#' `sign_matrix` required to calculate enrichment
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @param p_value logical. Default = `FALSE`. calculate p-values
#' @param include_depletion calculate both enrichment and depletion
#' @param n_times number of permutations to calculate for p_value
#' @param max_block number of lines to process together (default = 20e6)
#' @param name to give to spatial enrichment results, default = PAGE
#' @param verbose be verbose
#' @param return_gobject return giotto object
#' @returns `matrix` (`makeSignMatrixPAGE()`) and
#' `giotto` (`runPAGEEnrich(return_gobject = TRUE)`) or
#' `data.table` (`runPAGEEnrich(return_gobject = FALSE)`)
#' @details
#' The enrichment Z score is calculated by using method (PAGE) from
#' Kim SY et al., BMC bioinformatics, 2005 as \cr
#' \eqn{Z = ((Sm – mu)*m^(1/2)) / delta}. \cr
#' For each gene in each spot, mu is the fold change values versus the mean
#' expression and delta is the standard deviation. Sm is the mean fold change
#' value of a specific marker gene set and  m is the size of a given marker
#' gene set.
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' sign_list <- list(
#'     cell_type1 = c(
#'         "Bcl11b", "Lmo1", "F3", "Cnih3", "Ppp1r3c",
#'         "Rims2", "Gfap", "Gjc3", "Chrna4", "Prkcd"
#'     ),
#'     cell_type2 = c(
#'         "Prr18", "Grb14", "Tprn", "Clic1", "Olig2", "Hrh3",
#'         "Tmbim1", "Carhsp1", "Tmem88b", "Ugt8a"
#'     ),
#'     cell_type2 = c(
#'         "Arpp19", "Lamp5", "Galnt6", "Hlf", "Hs3st2",
#'         "Tbr1", "Myl4", "Cygb", "Ttc9b", "Ipcef1"
#'     )
#' )
#'
#' sm <- makeSignMatrixPAGE(
#'     sign_names = c("cell_type1", "cell_type2", "cell_type3"),
#'     sign_list = sign_list
#' )
#'
#' g <- runPAGEEnrich(
#'     gobject = g,
#'     sign_matrix = sm,
#'     min_overlap_genes = 2
#' )
#'
#' spatPlot2D(g,
#'     cell_color = "cell_type2",
#'     spat_enr_names = "PAGE",
#'     color_as_factor = FALSE
#' )
NULL

#' @rdname enrichment_PAGE
#' @export
makeSignMatrixPAGE <- function(
        sign_names,
        sign_list) {
    ## check input
    if (!inherits(sign_list, "list")) {
        stop("sign_list needs to be a list of signatures for each cell type /
            process")
    }
    if (length(sign_names) != length(sign_list)) {
        stop("the number of names needs to match the number of signatures
            provided")
    }

    ## create genes and signatures
    genes <- do.call("c", sign_list)
    types <- lapply(seq_along(sign_names), FUN = function(x) {
        subset <- sign_list[[x]]
        name_subset <- sign_names[[x]]

        res <- rep(x = name_subset, length(subset))
    })
    mydt <- data.table::data.table(
        genes = genes, types = unlist(types), value = 1
    )

    # convert data.table to signature matrix
    dtmatrix <- data.table::dcast.data.table(
        mydt,
        formula = genes ~ types, value.var = "value", fill = 0
    )
    final_sig_matrix <- Matrix::as.matrix(dtmatrix[, -1])
    rownames(final_sig_matrix) <- dtmatrix$genes

    return(final_sig_matrix)
}



## create spatialDWLS matrix ####

#' @title makeSignMatrixDWLSfromMatrix
#' @name makeSignMatrixDWLSfromMatrix
#' @description Function to convert a single-cell RNAseq matrix into a format
#'  that can be used with \code{\link{runDWLSDeconv}}.
#' @param matrix scRNA-seq matrix
#' @param sign_gene genes to use (e.g. marker genes)
#' @param cell_type_vector vector with cell types (length = ncol(matrix))
#' @returns matrix
#' @seealso \code{\link{runDWLSDeconv}}
#' @examples
#' sign_gene <- c(
#'     "Bcl11b", "Lmo1", "F3", "Cnih3", "Ppp1r3c", "Rims2", "Gfap",
#'     "Gjc3", "Chrna4", "Prkcd", "Prr18", "Grb14", "Tprn", "Clic1", "Olig2",
#'     "Hrh3", "Tmbim1", "Carhsp1", "Tmem88b", "Ugt8a", "Arpp19", "Lamp5",
#'     "Galnt6", "Hlf", "Hs3st2", "Tbr1", "Myl4", "Cygb", "Ttc9b", "Ipcef1"
#' )
#'
#' sign_matrix <- matrix(rnorm(length(sign_gene) * 3), nrow = length(sign_gene))
#' rownames(sign_matrix) <- sign_gene
#' colnames(sign_matrix) <- c("cell_type1", "cell_type2", "cell_type3")
#'
#' makeSignMatrixDWLSfromMatrix(
#'     matrix = sign_matrix, sign_gene = sign_gene,
#'     cell_type_vector = c("cell_type1", "cell_type2", "cell_type3")
#' )
#' @export
makeSignMatrixDWLSfromMatrix <- function(
        matrix,
        sign_gene,
        cell_type_vector) {
    # 1. check if cell_type_vector and matrix are compatible
    if (ncol(matrix) != length(cell_type_vector)) {
        stop("ncol(matrix) needs to be the same as length(cell_type_vector)")
    }

    # check input for sign_gene
    if (!is.character(sign_gene)) {
        stop("sign_gene needs to be a character vector of cell type specific
            genes")
    }


    # 2. get the common genes from the matrix and vector of signature genes
    intersect_sign_gene <- intersect(rownames(matrix), sign_gene)
    matrix_subset <- matrix[intersect_sign_gene, ]


    # 3. for each cell type
    # calculate average expression for all signature genes
    signMatrix <- matrix(
        data = NA,
        nrow = nrow(matrix_subset),
        ncol = length(unique(cell_type_vector))
    )

    for (cell_type_i in seq_along(unique(cell_type_vector))) {
        cell_type <- unique(cell_type_vector)[cell_type_i]
        selected_cells <- colnames(matrix_subset)[cell_type_vector == cell_type]
        mean_expr_in_selected_cells <- rowMeans_flex(matrix_subset[
            , selected_cells
        ])

        signMatrix[, cell_type_i] <- mean_expr_in_selected_cells
    }

    rownames(signMatrix) <- rownames(matrix_subset)
    colnames(signMatrix) <- unique(cell_type_vector)

    return(signMatrix)
}



#' @title makeSignMatrixDWLS
#' @description Function to convert a matrix within a Giotto object into a
#' format that can be used with \code{\link{runDWLSDeconv}} for deconvolution.
#' A vector of cell types for parameter \code{cell_type_vector} can be created
#' from the cell metadata (\code{\link{pDataDT}}).
#' @param gobject Giotto object of single cell
#' @param spat_unit spatial unit
#' @param feat_type feature type to use
#' @param expression_values expression values to use
#' @param reverse_log reverse a log-normalized expression matrix
#' @param log_base the logarithm base (default = 2)
#' @param sign_gene all of DE genes (signature)
#' @param cell_type_vector vector with cell types (length = ncol(matrix))
#' @param cell_type deprecated, use \code{cell_type_vector}
#' @returns matrix
#' @seealso \code{\link{runDWLSDeconv}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' sign_gene <- c(
#'     "Bcl11b", "Lmo1", "F3", "Cnih3", "Ppp1r3c", "Rims2", "Gfap",
#'     "Gjc3", "Chrna4", "Prkcd", "Prr18", "Grb14", "Tprn", "Clic1", "Olig2",
#'     "Hrh3", "Tmbim1", "Carhsp1", "Tmem88b", "Ugt8a", "Arpp19", "Lamp5",
#'     "Galnt6", "Hlf", "Hs3st2", "Tbr1", "Myl4", "Cygb", "Ttc9b", "Ipcef1"
#' )
#'
#' makeSignMatrixDWLS(
#'     gobject = g, sign_gene = sign_gene,
#'     cell_type_vector = pDataDT(g)[["leiden_clus"]]
#' )
#' @export
makeSignMatrixDWLS <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        reverse_log = TRUE,
        log_base = 2,
        sign_gene,
        cell_type_vector,
        cell_type = NULL) {
    ## deprecated arguments
    if (!is.null(cell_type)) {
        warning("cell_type is deprecated, use cell_type_vector in the future")
        cell_type_vector <- cell_type
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
        output = "exprObj"
    )

    ## 2. reverse log-normalization
    if (reverse_log == TRUE) {
        expr_values[] <- log_base^(expr_values[]) - 1
    }

    ## 3. run signature matrix function
    res <- makeSignMatrixDWLSfromMatrix(
        matrix = expr_values[],
        sign_gene = sign_gene,
        cell_type_vector = cell_type_vector
    )

    return(res)
}





#' @title makeSignMatrixRank
#' @description Function to convert a single-cell count matrix
#' and a corresponding single-cell cluster vector into
#' a rank matrix that can be used with the Rank enrichment option.
#' @param sc_matrix matrix of single-cell RNAseq expression data
#' @param sc_cluster_ids vector of cluster ids
#' @param ties_method how to handle rank ties
#' @param gobject if giotto object is given then only genes present in both
#' datasets will be considered
#' @returns matrix
#' @seealso \code{\link{rankEnrich}}
#' @examples
#' sign_gene <- c(
#'     "Bcl11b", "Lmo1", "F3", "Cnih3", "Ppp1r3c", "Rims2", "Gfap",
#'     "Gjc3", "Chrna4", "Prkcd", "Prr18", "Grb14", "Tprn", "Clic1", "Olig2",
#'     "Hrh3", "Tmbim1", "Carhsp1", "Tmem88b", "Ugt8a", "Arpp19", "Lamp5",
#'     "Galnt6", "Hlf", "Hs3st2", "Tbr1", "Myl4", "Cygb", "Ttc9b", "Ipcef1"
#' )
#'
#' sign_matrix <- matrix(rnorm(length(sign_gene) * 3), nrow = length(sign_gene))
#' rownames(sign_matrix) <- sign_gene
#' colnames(sign_matrix) <- c("cell_type1", "cell_type2", "cell_type3")
#'
#' makeSignMatrixRank(
#'     sc_matrix = sign_matrix,
#'     sc_cluster_ids = c("cell_type1", "cell_type2", "cell_type3")
#' )
#' @export
makeSignMatrixRank <- function(
        sc_matrix,
        sc_cluster_ids,
        ties_method = c("random", "max"),
        gobject = NULL) {
    if (inherits(sc_matrix, "exprObj")) {
        sc_matrix <- sc_matrix[]
    }
    if (!inherits(sc_matrix, c("matrix", "Matrix", "DelayedArray"))) {
        stop("'sc_matrix' must be a matrix input")
    }

    # select ties_method
    ties_method <- match.arg(ties_method, choices = c("random", "max"))

    # check input
    if (length(sc_cluster_ids) != ncol(sc_matrix)) {
        stop("Number of columns of the scRNAseq matrix needs to have the same
            length as the cluster ids")
    }

    mean_list <- list()
    group_list <- list()
    total_nr_genes <- nrow(sc_matrix)

    # calculate means for each cluster group
    for (group in seq_along(unique(sc_cluster_ids))) {
        group_id <- unique(sc_cluster_ids)[group]
        cell_ind <- which(sc_cluster_ids == group_id)
        cluster_rowmeans <- rowMeans_flex(sc_matrix[, cell_ind])
        mean_list[[group_id]] <- cluster_rowmeans
        group_list[[group]] <- rep(group_id, total_nr_genes)
    }

    mean_list_res <- data.table::as.data.table(do.call("c", mean_list))
    group_list_res <- data.table::as.data.table(do.call("c", group_list))

    # average expression for all cells
    av_expression <- rowMeans_flex(sc_matrix)
    av_expression_res <- rep(av_expression, length(unique(sc_cluster_ids)))

    gene_names <- rownames(sc_matrix)
    gene_names_res <- rep(gene_names, length(unique(sc_cluster_ids)))

    # create data.table with genes, mean expression per cluster, mean
    # expression overall and cluster ids
    comb_dt <- data.table::data.table(
        genes = gene_names_res,
        mean_expr = mean_list_res[[1]],
        av_expr = av_expression_res,
        clusters = group_list_res[[1]]
    )

    # data.table variables
    fold <- mean_expr <- av_expr <- rankFold <- clusters <- NULL

    # calculate fold change and rank of fold-change
    comb_dt[, fold := log2(mean_expr + 1) - log2(av_expr + 1)]
    comb_dt[, rankFold := data.table::frank(
        -fold,
        ties.method = ties_method
    ), by = clusters]

    # create matrix
    comb_rank_mat <- data.table::dcast.data.table(
        data = comb_dt, genes ~ clusters, value.var = "rankFold"
    )
    comb_rank_matrix <- dt_to_matrix(comb_rank_mat)
    comb_rank_matrix <- comb_rank_matrix[
        rownames(sc_matrix), unique(sc_cluster_ids)
    ]

    return(comb_rank_matrix)
}



# * ####
## spatial enrichment functions ####


#' @title PAGE permutation
#' @description creates permutation for the PAGEEnrich test
#' @returns PAGEEnrich test
#' @keywords internal
.do_page_permutation <- function(
        gobject,
        sig_gene,
        ntimes) {
    # check available gene
    available_ct <- c()
    for (i in colnames(sig_gene)) {
        gene_i <- rownames(sig_gene)[which(sig_gene[, i] == 1)]
        overlap_i <- intersect(
            gene_i, rownames(gobject@expression$rna$normalized)
        )
        if (length(overlap_i) <= 5) {
            output <- paste0(
                "Warning, ", i, " only has ", length(overlap_i),
                " overlapped genes. Will remove it."
            )
        } else {
            available_ct <- c(available_ct, i)
        }
    }
    if (length(available_ct) == 1) {
        stop("Only one cell type available.")
    }
    # only continue with genes present in both datasets
    interGene <- intersect(
        rownames(sig_gene), rownames(gobject@expression$rna$normalized)
    )
    sign_matrix <- sig_gene[interGene, available_ct]

    ct_gene_counts <- NULL
    for (i in seq_len(dim(sign_matrix)[2])) {
        a <- length(which(sign_matrix[, i] == 1))
        ct_gene_counts <- c(ct_gene_counts, a)
    }
    uniq_ct_gene_counts <- unique(ct_gene_counts)
    background_mean_sd <- matrix(
        data = NA, nrow = length(uniq_ct_gene_counts) + 1, ncol = 3
    )
    for (i in seq_along(uniq_ct_gene_counts)) {
        gene_num <- uniq_ct_gene_counts[i]
        all_sample_names <- NULL
        all_sample_list <- NULL
        for (j in seq_len(ntimes)) {
            set.seed(j)
            random_gene <- sample(rownames(
                gobject@expression$rna$normalized
            ), gene_num, replace = FALSE)
            ct_name <- paste("ct", j, sep = "")
            all_sample_names <- c(all_sample_names, ct_name)
            all_sample_list <- c(all_sample_list, list(random_gene))
        }
        random_sig <- makeSignMatrixPAGE(all_sample_names, all_sample_list)
        random_DT <- runPAGEEnrich(
            gobject,
            sign_matrix = random_sig,
            p_value = FALSE
        )
        background <- unlist(random_DT[, 2:dim(random_DT)[2]])
        df_row_name <- paste("gene_num_", uniq_ct_gene_counts[i], sep = "")
        list_back_i <- c(df_row_name, mean(background), stats::sd(background))
        background_mean_sd[i, ] <- list_back_i
    }
    background_mean_sd[length(uniq_ct_gene_counts) + 1, ] <- c("temp", "0", "1")
    df_back <- data.frame(background_mean_sd, row.names = 1)
    colnames(df_back) <- c("mean", "sd")
    return(df_back)
}





#' @title PAGE data.table method
#' @param expr_values matrix of expression values
#' @returns data.table
#' @keywords internal
.page_dt_method <- function(
        sign_matrix,
        expr_values,
        min_overlap_genes = 5,
        logbase = 2,
        reverse_log_scale = TRUE,
        output_enrichment = c("original", "zscore"),
        p_value = FALSE,
        include_depletion = FALSE,
        n_times = 1000,
        max_block = 20e6,
        verbose = TRUE) {
    # data.table variables
    Var1 <- value <- Var2 <- V1 <- marker <- nr_markers <- fc <- cell_ID <-
        zscore <- colmean <- colSd <- pval <- NULL
    mean_zscore <- sd_zscore <- pval_score <- NULL

    # output enrichment
    output_enrichment <- match.arg(
        output_enrichment,
        choices = c("original", "zscore")
    )

    ## identify available cell types
    all_genes <- rownames(expr_values)
    sign_matrix <- as.matrix(sign_matrix)
    sign_matrix_DT <- melt_matrix(sign_matrix)
    sign_matrix_DT <- sign_matrix_DT[Var1 %in% all_genes]
    detected_DT <- sign_matrix_DT[, sum(value), by = Var2]

    lost_cell_types_DT <- detected_DT[V1 <= min_overlap_genes]
    if (nrow(lost_cell_types_DT) > 0) {
        for (row in seq_len(nrow(lost_cell_types_DT))) {
            output <- paste0(
                "Warning, ", lost_cell_types_DT[row][["Var2"]], " only has ",
                lost_cell_types_DT[row][["V1"]],
                " overlapping genes. Will be removed."
            )
            if (verbose) print(output)
        }
    }
    available_ct <- as.character(detected_DT[V1 > min_overlap_genes][["Var2"]])

    if (length(available_ct) == 1) {
        stop("Only one cell type available.")
    }

    # create subset of signature matrix
    interGene <- intersect(rownames(sign_matrix), rownames(expr_values))
    filterSig <- sign_matrix[interGene, available_ct]

    # create fold expression for each gene in each spot
    # calculate mean gene expression
    if (reverse_log_scale == TRUE) {
        mean_gene_expr <- log(rowMeans(logbase^expr_values - 1, dims = 1) + 1)
    } else {
        mean_gene_expr <- rowMeans(expr_values)
    }
    geneFold <- expr_values - mean_gene_expr

    # calculate sample/spot mean and sd
    cellColMean <- colMeans(geneFold)
    cellColSd <- apply(geneFold, 2, stats::sd)
    cellColMeanSd <- data.table::data.table(
        cell_ID = names(cellColMean),
        colmean = cellColMean,
        colSd = cellColSd
    )

    filterSig_DT <- melt_matrix(filterSig)
    colnames(filterSig_DT) <- c("gene", "cell_type", "marker")
    sub_ct_DT <- filterSig_DT[marker == 1]
    sub_ct_DT[, nr_markers := .N, by = cell_type]

    ## reshape gene fold-expression
    geneFold_DT <- melt_matrix(geneFold)
    colnames(geneFold_DT) <- c("gene", "cell_ID", "fc")

    mergetest <- data.table::merge.data.table(
        sub_ct_DT, geneFold_DT,
        by = "gene"
    )
    mergetest <- mergetest[, mean(fc), by = .(cell_type, cell_ID, nr_markers)]
    if (is.integer(mergetest$cell_ID) && is.character(cellColMeanSd$cell_ID)) {
        mergetest$cell_ID <- as.character(mergetest$cell_ID)
    }
    mergetest <- data.table::merge.data.table(
        mergetest, cellColMeanSd,
        by = "cell_ID"
    )
    mergetest[, zscore := ((V1 - colmean) * nr_markers^(1 / 2)) / colSd]

    if (output_enrichment == "zscore") {
        mergetest[, zscore := scale(zscore), by = "cell_type"]
    }




    ## return p-values based on permutations ##
    if (p_value == TRUE) {
        ## 1. get number of markers instructions ##
        sample_intrs <- unique(sub_ct_DT[, .(cell_type, nr_markers)])


        ## 2. first create the random samples all together ##
        cell_type_list <- list()
        perm_type_list <- list()
        for (row in seq_len(nrow(sample_intrs))) {
            cell_type <- sample_intrs[row][["cell_type"]]
            nr_genes <- as.numeric(sample_intrs[row][["nr_markers"]])

            gene_list <- list()
            perm_list <- list()
            for (i in seq_len(n_times)) {
                sampled_genes <- sample(rownames(expr_values), size = nr_genes)
                gene_list[[i]] <- sampled_genes
                perm_list[[i]] <- rep(paste0("p_", i), nr_genes)
            }

            gene_res <- unlist(gene_list)
            names(gene_res) <- rep(cell_type, length(gene_res))
            cell_type_list[[row]] <- gene_res

            perm_res <- unlist(perm_list)
            perm_type_list[[row]] <- perm_res
        }

        cell_type_perm <- unlist(cell_type_list)
        perm_round <- unlist(perm_type_list)

        cell_type_perm_DT <- data.table::data.table(
            cell_type = names(cell_type_perm),
            gene = cell_type_perm,
            round = perm_round
        )

        sample_intrs_vec <- sample_intrs$nr_markers
        names(sample_intrs_vec) <- sample_intrs$cell_type
        cell_type_perm_DT[, nr_markers := sample_intrs_vec[cell_type]]


        ## 3. decide on number of blocks to process ##
        nr_perm_lines <- as.numeric(nrow(cell_type_perm_DT))
        nr_spots <- as.numeric(ncol(expr_values))
        total_lines <- nr_spots * nr_perm_lines
        nr_groups <- round(total_lines / max_block)

        ## 4. create groups
        all_perms <- unique(perm_round)
        all_perms_num <- seq_along(all_perms)
        names(all_perms_num) <- all_perms
        group_labels <- paste0("group_", seq_len(nr_groups))
        groups_vec <- cut(
            all_perms_num,
            breaks = nr_groups, labels = group_labels
        )
        names(all_perms) <- groups_vec


        ## 5. do random enrichment per block
        res_list <- list()
        for (group_i in seq_along(group_labels)) {
            group <- group_labels[group_i]
            sub_perms <- all_perms[names(all_perms) == group]
            cell_type_perm_DT_sub <- cell_type_perm_DT[round %in% sub_perms]

            mergetest_perm_sub <- data.table::merge.data.table(
                cell_type_perm_DT_sub, geneFold_DT,
                allow.cartesian = TRUE
            )
            mergetest_perm_sub <- mergetest_perm_sub[
                , mean(fc),
                by = .(cell_type, cell_ID, nr_markers, round)
            ]
            if (is.integer(mergetest_perm_sub$cell_ID) && is.character(
                cellColMeanSd$cell_ID
            )) {
                mergetest_perm_sub$cell_ID <- as.character(
                    mergetest_perm_sub$cell_ID
                )
            }
            mergetest_perm_sub <- data.table::merge.data.table(
                mergetest_perm_sub, cellColMeanSd,
                by = "cell_ID"
            )
            mergetest_perm_sub[, zscore := ((
                V1 - colmean) * nr_markers^(1 / 2)) / colSd]

            res_list[[group_i]] <- mergetest_perm_sub
        }

        res_list_comb <- do.call("rbind", res_list)
        res_list_comb_average <- res_list_comb[
            , .(mean_zscore = mean(zscore), sd_zscore = stats::sd(zscore)),
            by = c("cell_ID", "cell_type")
        ]
        mergetest_final <- data.table::merge.data.table(
            mergetest, res_list_comb_average,
            by = c("cell_ID", "cell_type")
        )

        ## calculate p.values based on normal distribution
        if (include_depletion == TRUE) {
            mergetest_final[, pval := stats::pnorm(
                abs(zscore),
                mean = mean_zscore, sd = sd_zscore,
                lower.tail = FALSE, log.p = FALSE
            )]
        } else {
            mergetest_final[, pval := stats::pnorm(
                zscore,
                mean = mean_zscore, sd = sd_zscore,
                lower.tail = FALSE, log.p = FALSE
            )]
        }

        data.table::setorder(mergetest_final, pval)

        ## calculate pval_score
        if (include_depletion == TRUE) {
            mergetest_final[, pval_score := sign(zscore) * -log10(pval)]
        } else {
            mergetest_final[, pval_score := -log10(pval)]
        }


        resultmatrix <- data.table::dcast(
            mergetest_final,
            formula = cell_ID ~ cell_type,
            value.var = "pval_score"
        )
        return(list(DT = mergetest_final, matrix = resultmatrix))
    } else {
        resultmatrix <- data.table::dcast(
            mergetest,
            formula = cell_ID ~ cell_type, value.var = "zscore"
        )
        return(list(DT = mergetest, matrix = resultmatrix))
    }
}



#' @rdname enrichment_PAGE
#' @export
runPAGEEnrich <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        sign_matrix,
        expression_values = c("normalized", "scaled", "custom"),
        min_overlap_genes = 5,
        reverse_log_scale = TRUE,
        logbase = 2,
        output_enrichment = c("original", "zscore"),
        p_value = FALSE,
        include_depletion = FALSE,
        n_times = 1000,
        max_block = 20e6,
        name = NULL,
        verbose = TRUE,
        return_gobject = TRUE) {
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
        unique(c("normalized", "scaled", "custom"), expression_values)
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )

    # check parameters
    if (is.null(name)) name <- "PAGE"

    PAGE_results <- .page_dt_method(
        sign_matrix = sign_matrix,
        expr_values = as.matrix(expr_values[]),
        min_overlap_genes = min_overlap_genes,
        logbase = logbase,
        reverse_log_scale = reverse_log_scale,
        output_enrichment = c("original", "zscore"),
        p_value = p_value,
        include_depletion = include_depletion,
        n_times = n_times,
        max_block = max_block,
        verbose = verbose
    )

    # create spatial enrichment object
    enrObj <- create_spat_enr_obj(
        name = name,
        method = "PAGE",
        enrichDT = PAGE_results[["matrix"]],
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = expr_values@provenance,
        misc = list(
            expr_values_used = expression_values,
            reverse_log_scale = reverse_log_scale,
            logbase = logbase,
            p_values_calculated = p_value,
            output_enrichment_scores = output_enrichment,
            include_depletion = include_depletion,
            nr_permutations = n_times
        )
    )

    ## return object or results ##
    if (return_gobject == TRUE) {
        spenr_names <- list_spatial_enrichments_names(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        )

        if (name %in% spenr_names) {
            cat(name, " has already been used, will be overwritten")
        }

        ## update parameters used ##
        parameters_list <- gobject@parameters
        number_of_rounds <- length(parameters_list)
        update_name <- paste0(number_of_rounds, "_spatial_enrichment")

        # parameters to include
        parameters_list[[update_name]] <- c(
            "method used" = "PAGE",
            "enrichment name" = name,
            "expression values" = expression_values,
            "reverse log scale" = reverse_log_scale,
            "logbase" = logbase,
            "output enrichment scores" = output_enrichment,
            "p values calculated" = p_value,
            "include depletion" = include_depletion,
            "nr permutations" = n_times
        )

        gobject@parameters <- parameters_list

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, enrObj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


        return(gobject)
    } else {
        PAGE_results[["matrix"]] <- enrObj
        return(PAGE_results)
    }
}









#' @title Rank permutation
#' @description creates permutation for the rankEnrich test
#' @returns data.frame
#' @keywords internal
.do_rank_permutation <- function(sc_gene, n) {
    random_df <- data.frame(matrix(ncol = n, nrow = length(sc_gene)))
    for (i in seq_len(n)) {
        set.seed(i)
        random_rank <- sample(
            seq_along(sc_gene), length(sc_gene),
            replace = FALSE
        )
        random_df[, i] <- random_rank
    }
    rownames(random_df) <- sc_gene
    return(random_df)
}


#' @title runRankEnrich
#' @description Function to calculate gene signature enrichment scores per
#' spatial position using a rank based approach.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @param ties_method how to handle rank ties
#' @param p_value calculate p-values (boolean, default = FALSE)
#' @param n_times number of permutations to calculate for p_value
#' @param rbp_p fractional binarization threshold (default = 0.99)
#' @param num_agg number of top genes to aggregate (default = 100)
#' @param name to give to spatial enrichment results, default = rank
#' @param return_gobject return giotto object
#' @returns data.table with enrichment results
#' @details
#' sign_matrix: a rank-fold matrix with genes as row names and cell-types as
#' column names.
#' Alternatively a scRNA-seq matrix and vector with clusters can be provided
#' to makeSignMatrixRank, which will create the matrix for you. \cr
#'
#' First a new rank is calculated as R = (R1*R2)^(1/2), where R1 is the rank of
#' fold-change for each gene in each spot and R2 is the rank of each marker in
#' each cell type.
#' The Rank-Biased Precision is then calculated as:
#' RBP = (1 - 0.99) * (0.99)^(R - 1)
#' and the final enrichment score is then calculated as the sum of top 100 RBPs.
#' @seealso \code{\link{makeSignMatrixRank}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- findMarkers_one_vs_all(g,
#'     cluster_column = "leiden_clus", min_feats = 20
#' )
#' sign_gene <- x$feats
#'
#' sign_matrix <- matrix(rnorm(length(sign_gene) * 8, mean = 10),
#'     nrow = length(sign_gene)
#' )
#' rownames(sign_matrix) <- sign_gene
#' colnames(sign_matrix) <- paste0("celltype_", unique(x$cluster))
#'
#' runRankEnrich(
#'     gobject = g, sign_matrix = sign_matrix,
#'     expression_values = "normalized"
#' )
#' @export
runRankEnrich <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        sign_matrix,
        expression_values = c("normalized", "raw", "scaled", "custom"),
        reverse_log_scale = TRUE,
        logbase = 2,
        output_enrichment = c("original", "zscore"),
        ties_method = c("average", "max"),
        p_value = FALSE,
        n_times = 1000,
        rbp_p = 0.99,
        num_agg = 100,
        name = NULL,
        return_gobject = TRUE) {
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

    # determine ties.method
    ties_method <- match.arg(ties_method, choices = c("average", "max"))

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
        output = "exprObj"
    )

    if (values == "raw") {
        expr_values[] <- Matrix::as.matrix(expr_values[])
    }

    # check parameters
    if (is.null(name)) name <- "rank"

    # check gene list
    interGene <- intersect(rownames(sign_matrix), rownames(expr_values[]))
    if (length(interGene) < 100) {
        stop("Please check the gene numbers or names of scRNA-seq. The names
            of scRNA-seq should be consistent with spatial data.")
    }

    # output enrichment
    output_enrichment <- match.arg(
        output_enrichment,
        choices = c("original", "zscore")
    )

    enrichment <- matrix(
        data = NA,
        nrow = dim(sign_matrix)[2],
        ncol = dim(expr_values[])[2]
    )

    # calculate mean gene expression
    if (reverse_log_scale == TRUE) {
        mean_gene_expr <- log(Matrix::rowMeans(
            logbase^expr_values[] - 1,
            dims = 1
        ) + 1)
    } else {
        mean_gene_expr <- Matrix::rowMeans(expr_values[])
    }

    # fold change and ranking

    ties_1 <- ties_method
    ties_2 <- ties_method
    if (ties_method == "max") {
        ties_1 <- "min"
        ties_2 <- "max"
    }
    # else ties_1=ties_2 is equal to random
    geneFold <- expr_values[]
    geneFold <- sparseMatrixStats::rowRanks(geneFold, ties.method = ties_1)
    rankFold <- t(sparseMatrixStats::colRanks(-geneFold, ties.method = ties_2))

    rownames(rankFold) <- rownames(expr_values[])
    colnames(rankFold) <- colnames(expr_values[])

    for (i in seq_len(dim(sign_matrix)[2])) {
        signames <- rownames(sign_matrix)[which(sign_matrix[, i] > 0)]
        interGene <- intersect(signames, rownames(rankFold))
        filterSig <- sign_matrix[interGene, ]
        filterRankFold <- rankFold[interGene, ]

        multiplyRank <- (filterRankFold * filterSig[, i])^(1 / 2)
        rpb <- (1.0 - rbp_p) * (rbp_p^(multiplyRank - 1))

        vectorX <- rep(NA, dim(filterRankFold)[2])

        for (j in seq_len(dim(filterRankFold)[2])) {
            toprpb <- sort(rpb[, j], decreasing = TRUE)
            zscore <- sum(toprpb[seq_len(num_agg)])
            vectorX[j] <- zscore
        }
        enrichment[i, ] <- vectorX
    }

    rownames(enrichment) <- colnames(sign_matrix)
    colnames(enrichment) <- colnames(rankFold)

    enrichment <- t(enrichment)

    if (output_enrichment == "zscore") {
        enrichment <- scale(enrichment)
    }

    enrichmentDT <- data.table::data.table(cell_ID = rownames(enrichment))
    enrichmentDT <- cbind(enrichmentDT, data.table::as.data.table(enrichment))


    # default name for page enrichment

    if (isTRUE(p_value)) {
        random_rank <- .do_rank_permutation(
            sc_gene = rownames(sign_matrix),
            n = n_times
        )

        random_DT <- runRankEnrich(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sign_matrix = random_rank,
            expression_values = expression_values,
            reverse_log_scale = reverse_log_scale,
            logbase = logbase,
            output_enrichment = output_enrichment,
            p_value = FALSE
        )

        background <- unlist(random_DT[, 2:dim(random_DT)[2]])
        fit.gamma <- fitdistrplus::fitdist(
            background,
            distr = "gamma", method = "mle"
        )
        pvalue_DT <- enrichmentDT
        enrichmentDT[, 2:dim(enrichmentDT)[2]] <- lapply(
            enrichmentDT[, 2:dim(enrichmentDT)[2]], function(x) {
                stats::pgamma(
                    x, fit.gamma$estimate[1],
                    rate = fit.gamma$estimate[2],
                    lower.tail = FALSE, log.p = FALSE
                )
            }
        )
    }

    # create spatial enrichment object
    enrObj <- create_spat_enr_obj(
        name = name,
        method = "rank",
        enrichDT = enrichmentDT,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = expr_values@provenance,
        misc = list(
            expr_values_used = expression_values,
            reverse_log_scale = reverse_log_scale,
            logbase = logbase,
            p_values_calculated = p_value,
            output_enrichment_scores = output_enrichment,
            nr_permutations = n_times
        )
    )

    ## return object or results ##
    if (return_gobject == TRUE) {
        spenr_names <- list_spatial_enrichments_names(
            gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
        )


        if (name %in% spenr_names) {
            cat(name, " has already been used, will be overwritten")
        }

        ## update parameters used ##
        parameters_list <- gobject@parameters
        number_of_rounds <- length(parameters_list)
        update_name <- paste0(number_of_rounds, "_spatial_enrichment")

        # parameters to include
        parameters_list[[update_name]] <- c(
            "method used" = "rank",
            "enrichment name" = name,
            "expression values" = expression_values,
            "reverse log scale" = reverse_log_scale,
            "logbase" = logbase,
            "p-values calculated" = p_value,
            "output enrichment scores" = output_enrichment,
            "p values calculated" = p_value,
            "nr permutations" = n_times
        )
        gobject@parameters <- parameters_list

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, enrObj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        return(gobject)
    } else {
        return(enrObj)
    }
}




#' @title runHyperGeometricEnrich
#' @description Function to calculate gene signature enrichment scores per
#' spatial position using a hypergeometric test.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param top_percentage percentage of cells that will be considered to have
#' gene expression with matrix binarization
#' @param output_enrichment how to return enrichment output
#' @param p_value calculate p-values (boolean, default = FALSE)
#' @param name to give to spatial enrichment results, default = hypergeometric
#' @param return_gobject return giotto object
#' @returns data.table with enrichment results
#' @details The enrichment score is calculated based on the p-value from the
#' hypergeometric test, -log10(p-value).
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- findMarkers_one_vs_all(g,
#'     cluster_column = "leiden_clus", min_feats = 20
#' )
#' sign_gene <- x$feats
#'
#' sign_matrix <- matrix(rnorm(length(sign_gene) * 8, mean = 10),
#'     nrow = length(sign_gene)
#' )
#' rownames(sign_matrix) <- sign_gene
#' colnames(sign_matrix) <- paste0("celltype_", unique(x$cluster))
#'
#' runHyperGeometricEnrich(gobject = g, sign_matrix = sign_matrix)
#' @export
runHyperGeometricEnrich <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        sign_matrix,
        expression_values = c("normalized", "scaled", "custom"),
        reverse_log_scale = TRUE,
        logbase = 2,
        top_percentage = 5,
        output_enrichment = c("original", "zscore"),
        p_value = FALSE,
        name = NULL,
        return_gobject = TRUE) {
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

    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )


    # check parameters
    if (is.null(name)) name <- "hypergeometric"

    # output enrichment
    output_enrichment <- match.arg(
        output_enrichment,
        choices = c("original", "zscore")
    )

    # calculate mean gene expression
    if (reverse_log_scale == TRUE) {
        expr_values[] <- logbase^expr_values[] - 1
    }

    interGene <- intersect(rownames(expr_values[]), rownames(sign_matrix))

    inter_sign_matrix <- sign_matrix[interGene, ]

    aveExp <- log2(2 * (Matrix::rowMeans(2^(expr_values[] - 1), dims = 1)) + 1)

    foldChange <- expr_values[] - aveExp

    top_q <- 1 - top_percentage / 100
    quantilecut <- apply(
        foldChange, 2, stats::quantile,
        probs = top_q, na.rm = TRUE
    )
    expbinary <- t_flex(1 * t_flex(foldChange > quantilecut))

    markerGenes <- rownames(inter_sign_matrix)
    expbinaryOverlap <- expbinary[markerGenes, ]
    total <- length(markerGenes)
    enrichment <- matrix(
        data = NA,
        nrow = dim(inter_sign_matrix)[2],
        ncol = dim(expbinaryOverlap)[2]
    )

    for (i in seq_len(dim(inter_sign_matrix)[2])) {
        signames <- rownames(inter_sign_matrix)[
            which(inter_sign_matrix[, i] == 1)
        ]
        vectorX <- NULL

        for (j in seq_len(dim(expbinaryOverlap)[2])) {
            cellsiggene <- names(expbinaryOverlap[
                which(expbinaryOverlap[, j] == 1), j
            ])
            x <- length(intersect(cellsiggene, signames))
            m <- length(rownames(inter_sign_matrix)[which(
                inter_sign_matrix[, i] == 1
            )])
            n <- total - m
            k <- length(intersect(cellsiggene, markerGenes))
            enrich <- (0 - log10(stats::phyper(
                x, m, n, k,
                log.p = FALSE, lower.tail = FALSE
            )))
            vectorX <- append(vectorX, enrich)
        }
        enrichment[i, ] <- vectorX
    }

    rownames(enrichment) <- colnames(inter_sign_matrix)
    colnames(enrichment) <- colnames(expbinaryOverlap)

    enrichment <- t(enrichment)

    if (output_enrichment == "zscore") {
        enrichment <- scale(enrichment)
    }

    enrichmentDT <- data.table::data.table(cell_ID = rownames(enrichment))
    enrichmentDT <- cbind(enrichmentDT, data.table::as.data.table(enrichment))


    ## calculate p-values ##
    if (p_value == TRUE) {
        enrichmentDT[, 2:dim(enrichmentDT)[2]] <- lapply(
            enrichmentDT[, 2:dim(enrichmentDT)[2]], function(x) {
                10^(-x)
            }
        )
    }

    # create spatial enrichment object
    enrObj <- create_spat_enr_obj(
        name = name,
        method = "hypergeometric",
        enrichDT = enrichmentDT,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = expr_values@provenance,
        misc = list(
            expr_values_used = expression_values,
            reverse_log_scale = reverse_log_scale,
            logbase = logbase,
            top_percentage = top_percentage,
            p_values_calculated = p_value,
            output_enrichment_scores = output_enrichment
        )
    )

    ## return object or results ##
    if (return_gobject == TRUE) {
        spenr_names <- list_spatial_enrichments_names(
            gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
        )

        if (name %in% spenr_names) {
            cat(name, " has already been used, will be overwritten")
        }

        ## update parameters used ##
        parameters_list <- gobject@parameters
        number_of_rounds <- length(parameters_list)
        update_name <- paste0(number_of_rounds, "_spatial_enrichment")

        # parameters to include
        parameters_list[[update_name]] <- c(
            "method used" = "hypergeometric",
            "enrichment name" = name,
            "expression values" = expression_values,
            "reverse log scale" = reverse_log_scale,
            "logbase" = logbase,
            "top percentage" = top_percentage,
            "p-values calculated" = p_value,
            "output enrichment scores" = output_enrichment,
            "p values calculated" = p_value
        )
        gobject@parameters <- parameters_list

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, enrObj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        return(gobject)
    } else {
        return(enrObj)
    }
}





#' @title runSpatialEnrich
#' @description Function to calculate gene signature enrichment scores per
#' spatial position using an enrichment test.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param enrich_method method for gene signature enrichment calculation
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param min_overlap_genes minimum number of overlapping genes in sign_matrix
#' required to calculate enrichment (PAGE)
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param p_value calculate p-value (default = FALSE)
#' @param n_times (page/rank) number of permutation iterations to calculate
#' p-value
#' @param rbp_p (rank) fractional binarization threshold (default = 0.99)
#' @param num_agg (rank) number of top genes to aggregate (default = 100)
#' @param max_block number of lines to process together (default = 20e6)
#' @param top_percentage (hyper) percentage of cells that will be considered
#' to have gene expression with matrix binarization
#' @param output_enrichment how to return enrichment output
#' @param name to give to spatial enrichment results, default = PAGE
#' @param verbose be verbose
#' @param return_gobject return giotto object
#' @returns Giotto object or enrichment results if return_gobject = FALSE
#' @details For details see the individual functions:
#'   * **PAGE:** \code{\link{runPAGEEnrich}}
#'   * **Rank:** \code{\link{runRankEnrich}}
#'   * **Hypergeometric:** \code{\link{runHyperGeometricEnrich}}
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- findMarkers_one_vs_all(g,
#'     cluster_column = "leiden_clus", min_feats = 20
#' )
#' sign_gene <- x$feats
#'
#' sign_matrix <- matrix(rnorm(length(sign_gene) * 8, mean = 10),
#'     nrow = length(sign_gene)
#' )
#' rownames(sign_matrix) <- sign_gene
#' colnames(sign_matrix) <- paste0("celltype_", unique(x$cluster))
#'
#' runSpatialEnrich(gobject = g, sign_matrix = sign_matrix)
#' @export
runSpatialEnrich <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        enrich_method = c("PAGE", "rank", "hypergeometric"),
        sign_matrix,
        expression_values = c("normalized", "scaled", "custom"),
        min_overlap_genes = 5,
        reverse_log_scale = TRUE,
        logbase = 2,
        p_value = FALSE,
        n_times = 1000,
        rbp_p = 0.99,
        num_agg = 100,
        max_block = 20e6,
        top_percentage = 5,
        output_enrichment = c("original", "zscore"),
        name = NULL,
        verbose = TRUE,
        return_gobject = TRUE) {
    enrich_method <- match.arg(
        enrich_method,
        choices = c("PAGE", "rank", "hypergeometric")
    )

    output_enrichment <- match.arg(
        output_enrichment,
        choices = c("original", "zscore")
    )


    if (enrich_method == "PAGE") {
        results <- runPAGEEnrich(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sign_matrix = sign_matrix,
            expression_values = expression_values,
            reverse_log_scale = reverse_log_scale,
            logbase = logbase,
            output_enrichment = output_enrichment,
            p_value = p_value,
            n_times = n_times,
            name = name,
            return_gobject = return_gobject
        )
    } else if (enrich_method == "rank") {
        results <- runRankEnrich(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sign_matrix = sign_matrix,
            expression_values = expression_values,
            reverse_log_scale = reverse_log_scale,
            logbase = logbase,
            output_enrichment = output_enrichment,
            p_value = p_value,
            n_times = n_times,
            rbp_p = rbp_p,
            num_agg = num_agg,
            name = name,
            return_gobject = return_gobject
        )
    } else if (enrich_method == "hypergeometric") {
        results <- runHyperGeometricEnrich(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sign_matrix = sign_matrix,
            expression_values = expression_values,
            reverse_log_scale = reverse_log_scale,
            logbase = logbase,
            top_percentage = top_percentage,
            output_enrichment = output_enrichment,
            p_value = p_value,
            name = name,
            return_gobject = return_gobject
        )
    }

    return(results)
}





# * ####
## spatial autocorrelation functions ####


#' @title Spatial autocorrelation
#' @name spatialAutoCor
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param feats features (expression) on which to run autocorrelation.
#' (leaving as NULL means that all features will be tested)
#' @param method method of autocorrelation. See details (default = 'moran')
#' @param data_to_use if using data from gobject, whether to test using
#' expression ('expression') or cell metadata ('cell_meta')
#' @param expression_values name of expression information to use
#' @param meta_cols columns in cell metadata to test
#' @param spatial_network_to_use spatial network to use
#' @param wm_method type of weight matrix to generate from spatial network
#' if no weight matrix is found attached to the spatial network
#' @param wm_name name of attached weight matrix to use
#' @param node_values alternative method of directly supplying a set of node
#' values
#' @param weight_matrix alternative method of directly supplying a spatial
#' weight matrix
#' @param test_method method to test values for significance (default is no
#' testing)
#' @param verbose be verbose
#' @returns spatial autocorrelation
#' @description Find spatial autocorrelation. Note that
#' \code{spatialAutoCorGlobal} will return values as a data.table instead of
#' appending information to the gobject.
#' \code{spatialAutoCorLocal} will append the results as a spatial enrichment
#' object by default. \cr
#' If providing external data using either the \code{node_values} and/or
#' \code{weight_matrix} params, the order of values provided should be the
#' same as the ordering of the columns and rows of the weight matrix.
NULL


# internals for spatial autocorrelation using terra
#' @keywords internal
.spat_autocor_terra_numeric <- function(x, w, method) {
    return(terra::autocor(x = x, w = w, method = method))
}

#' @keywords internal
.spat_autocor_terra_raster <- function(x, w, global = TRUE, method) {
    return(terra::autocor(x = x, w = w, global = global, method = method))
}





#' @describeIn spatialAutoCor Global autocorrelation (single value returned)
#'
#' @param mc_nsim when \code{test_method = 'monte_carlo'} this is number of
#' simulations to perform
#' @param cor_name name to assign the results in global autocorrelation output
#' @param return_gobject (default = FALSE) whether to return results appended to
#' metadata in the giotto object or as a data.table
#' @details
#' \strong{Global Methods:}
#' \itemize{
#'   \item{\emph{Moran's I} 'moran'}
#'   \item{\emph{Geary's C} 'geary'}
#' }
#' @export
spatialAutoCorGlobal <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        feats = NULL,
        method = c("moran", "geary"),
        data_to_use = c("expression", "cell_meta"),
        expression_values = c("normalized", "scaled", "custom"),
        meta_cols = NULL,
        spatial_network_to_use = "kNN_network",
        wm_method = c("distance", "adjacency"),
        wm_name = "spat_weights",
        node_values = NULL,
        weight_matrix = NULL,
        test_method = c("none", "monte_carlo"),
        mc_nsim = 99,
        cor_name = NULL,
        return_gobject = FALSE,
        verbose = TRUE) {
    # 0. determine inputs
    method <- match.arg(method, choices = c("moran", "geary"))
    test_method <- match.arg(test_method, choices = c("none", "monte_carlo"))
    data_to_use <- match.arg(
        data_to_use,
        choices = c("expression", "cell_meta")
    )
    if (is.null(cor_name)) cor_name <- method
    if (!is.null(node_values)) {
        if (is.numeric(node_values)) {
            stop(wrap_txt('External "node_values" must be type numeric.',
                errWidth = TRUE
            ))
        }
    }

    use_ext_vals <- data.table::fifelse(
        !is.null(node_values),
        yes = TRUE, no = FALSE
    )
    use_sn <- data.table::fifelse(
        !is.null(weight_matrix),
        yes = FALSE, no = TRUE
    )

    use_expr <- data.table::fcase(
        isTRUE(use_ext_vals), FALSE,
        data_to_use != "expression", FALSE,
        default = TRUE
    )

    use_meta <- data.table::fcase(
        isTRUE(use_ext_vals), FALSE,
        data_to_use != "cell_meta", FALSE,
        default = TRUE
    )

    if (data_to_use == "cell_meta") {
        if (is.null(meta_cols)) {
            stop(wrap_txt(
                'If "data_to_use" is "cell_meta" then a character vector of
                cell metadata',
                'columns to use must be provided in "meta_cols"',
                errWidth = TRUE
            ))
        }
    }
    if (isTRUE(return_gobject)) {
        if (data_to_use == "cell_meta" | isTRUE(use_ext_vals)) {
            stop(wrap_txt(
                "Global spatial autocorrelations on cell_meta or external data
                should not",
                "be returned to the gobject.
         > Please set return_gobject = FALSE",
                errWidth = TRUE
            ))
        }
    }

    # 1. setup
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
    } else { # if null
        if (any(!use_ext_vals, use_sn, return_gobject)) {
            stop("gobject has not been provided")
        }
    }

    # select and format input
    data_list <- .evaluate_autocor_input(
        gobject = gobject,
        use_ext_vals = use_ext_vals,
        use_sn = use_sn,
        use_expr = use_expr,
        use_meta = use_meta,
        spat_unit = spat_unit,
        feat_type = feat_type,
        feats = feats,
        data_to_use = data_to_use,
        expression_values = expression_values,
        meta_cols = meta_cols,
        spatial_network_to_use = spatial_network_to_use,
        wm_method = wm_method,
        wm_name = wm_name,
        node_values = node_values,
        weight_matrix = weight_matrix,
        verbose = verbose
    )
    # unpack formatted data
    use_values <- data_list$use_values
    feats <- data_list$feats
    weight_matrix <- data_list$weight_matrix


    # 2. perform autocor
    res_dt <- .run_spat_autocor_global(
        use_values = use_values,
        feats = feats,
        weight_matrix = weight_matrix,
        method = method,
        test_method = test_method,
        mc_nsim = mc_nsim,
        cor_name = cor_name
    )

    # return info
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
            new_metadata = res_dt,
            by_column = TRUE,
            column_feat_ID = "feat_ID"
        )

        return(gobject)
    } else {
        return(res_dt)
    }
}


#' @describeIn spatialAutoCor Local autocorrelation
#' (values generated for each spatial ID)
#'
#' @param enrich_name name to assign local autocorrelation spatial enrichment
#' results
#' @param return_gobject (default = FALSE) whether to return results appended to
#' @param output 'spatEnrObj' or 'data.table'
#' metadata in the giotto object or as a data.table
#' @details
#' \strong{Local Methods:}
#' \itemize{
#'   \item{\emph{Local Moran's I} 'moran'}
#'   \item{\emph{Getis-Ord Gi} 'Gi'}
#'   \item{\emph{Getis-Ord Gi*} 'Gi*'}
#'   \item{\emph{Local mean} 'mean'}
#' }
#' @export
spatialAutoCorLocal <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        feats = NULL,
        method = c("moran", "gi", "gi*", "mean"),
        data_to_use = c("expression", "cell_meta"),
        expression_values = c("normalized", "scaled", "custom"),
        meta_cols = NULL,
        spatial_network_to_use = "kNN_network",
        wm_method = c("distance", "adjacency"),
        wm_name = "spat_weights",
        node_values = NULL,
        weight_matrix = NULL,
        test_method = c("none"),
        # cor_name = NULL,
        enrich_name = NULL,
        return_gobject = TRUE,
        output = c("spatEnrObj", "data.table"),
        verbose = TRUE) {
    # 0. determine inputs
    method_select <- match.arg(
        method,
        choices = c("moran", "gi", "gi*", "mean")
    )
    data_to_use <- match.arg(
        data_to_use,
        choices = c("expression", "cell_meta")
    )
    output <- match.arg(output, choices = c("spatEnrObj", "data.table"))
    # if(is.null(cor_name)) cor_name = method

    if (method_select == "moran") {
        method <- "locmor"
    } else {
        method <- method_select
    }

    if (!is.null(node_values)) {
        if (is.numeric(node_values)) {
            stop(wrap_txt('External "node_values" must be type numeric',
                errWidth = TRUE
            ))
        }
    }


    use_ext_vals <- data.table::fifelse(
        !is.null(node_values),
        yes = TRUE, no = FALSE
    )
    use_sn <- data.table::fifelse(
        !is.null(weight_matrix),
        yes = FALSE, no = TRUE
    )

    use_expr <- data.table::fcase(
        isTRUE(use_ext_vals), FALSE,
        data_to_use != "expression", FALSE,
        default = TRUE
    )

    use_meta <- data.table::fcase(
        isTRUE(use_ext_vals), FALSE,
        data_to_use != "cell_meta", FALSE,
        default = TRUE
    )

    if (is.null(enrich_name)) { # name of spatEnrObj
        enrich_name <- data.table::fcase(
            isTRUE(use_ext_vals), method_select,
            data_to_use == "expression", paste0("expr_", method_select),
            data_to_use == "cell_meta", paste0("meta_", method_select),
            default = method_select
        )
    }



    # 1. setup
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
    } else { # if null
        if (any(!use_ext_vals, use_sn, return_gobject)) {
            stop("gobject has not been provided")
        }
    }


    # select and format input
    data_list <- .evaluate_autocor_input(
        gobject = gobject,
        use_ext_vals = use_ext_vals,
        use_sn = use_sn,
        use_expr = use_expr,
        use_meta = use_meta,
        spat_unit = spat_unit,
        feat_type = feat_type,
        feats = feats,
        data_to_use = data_to_use,
        expression_values = expression_values,
        meta_cols = meta_cols,
        spatial_network_to_use = spatial_network_to_use,
        wm_method = wm_method,
        wm_name = wm_name,
        node_values = node_values,
        weight_matrix = weight_matrix,
        verbose = verbose
    )
    # unpack formatted input
    use_values <- data_list$use_values
    feats <- data_list$feats
    weight_matrix <- data_list$weight_matrix
    provenance <- data_list$provenance
    values <- data_list$expr_values
    IDs <- data_list$IDs

    # spatIDs to use when returning autocor results
    # Provide default spatIDs if missing
    if (is.null(IDs)) {
        IDs <- seq(nrow(use_values))
    }

    # 2. perform autocor
    res_dt <- .run_spat_autocor_local(
        use_values = use_values,
        feats = feats,
        weight_matrix = weight_matrix,
        method = method,
        test_method = test_method,
        IDs = IDs
    )

    # create spatial enrichment object
    enr <- create_spat_enr_obj(
        name = enrich_name,
        method = method_select,
        enrichDT = res_dt,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        misc = if (use_expr) list(expr_values_used = values)
    )


    # return info
    if (isTRUE(return_gobject)) {
        if (isTRUE(verbose)) {
            wrap_msg("Attaching ", method_select,
                ' results as spatial enrichment: "',
                enrich_name, '"',
                sep = ""
            )
        }
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, enr)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        return(gobject)
    } else {
        if (output == "spatEnrObj") {
            return(enr)
        }
        if (output == "data.table") {
            return(res_dt)
        }
    }
}







#' .run_spat_autocor_global
#' @returns data.table
#' @keywords internal
.run_spat_autocor_global <- function(
        use_values,
        feats,
        weight_matrix,
        method,
        test_method,
        mc_nsim,
        cor_name) {
    # data.table vars
    cell_ID <- nsim <- NULL

    nfeats <- length(feats)
    if (test_method != "none") {
        step_size <- ceiling(nfeats / 100L)
    } else {
        step_size <- step_size <- ceiling(nfeats / 10L)
    }

    with_pbar({
        if (step_size > 1) {
            pb <- pbar(
                steps = nfeats / step_size
            )
        }
        res_list <- lapply_flex(
            seq_along(feats),
            # future.packages = c('terra', 'data.table'),
            function(feat_i) {
                feat <- feats[feat_i]
                if (inherits(use_values, "data.table")) {
                    feat_vals <- eval(call("[", use_values, j = as.name(feat)))
                } else {
                    feat_vals <- use_values[, feat]
                }


                spat_ac <- .spat_autocor_terra_numeric(
                    x = feat_vals,
                    w = weight_matrix,
                    method = method
                )


                # test
                if (test_method != "none") {
                    if (test_method == "monte_carlo") {
                        mc <- sapply(seq(mc_nsim), function(i) {
                            .spat_autocor_terra_numeric(
                                x = sample(feat_vals),
                                w = weight_matrix,
                                method = method
                            )
                        })
                        P <- 1 - sum((spat_ac > mc) / (nsim + 1))
                    }
                    if (test_method == "spdep") {
                        wrap_msg("spdep not yet implemented")
                    }
                }
                # increment progress
                if (exists("pb")) if (feat_i %% step_size == 0) pb()


                if (test_method == "none") {
                    return(data.table::data.table(feat, spat_ac))
                } else {
                    return(data.table::data.table(feat, spat_ac, P))
                }
            }
        )
    })
    res_dt <- do.call("rbind", res_list)
    if (test_method == "none") {
        colnames(res_dt) <- c("feat_ID", cor_name)
    } else {
        colnames(res_dt) <- c("feat_ID", cor_name, paste0(
            cor_name, "_", test_method
        ))
    }
    return(res_dt)
}

#' .run_spat_autocor_local
#' @returns data.table
#' @keywords internal
.run_spat_autocor_local <- function(
        use_values,
        feats,
        weight_matrix,
        method,
        test_method,
        IDs) {
    cell_ID <- NULL

    nfeats <- length(feats)
    if (test_method != "none") {
        step_size <- ceiling(nfeats / 100L)
    } else {
        step_size <- step_size <- ceiling(nfeats / 10L)
    }

    with_pbar({
        if (step_size > 1) {
            pb <- pbar(
                steps = nfeats / step_size
            )
        }
        res_list <- lapply_flex(
            seq_along(feats),
            # future.packages = c('terra', 'data.table'),
            function(feat_i) {
                feat <- feats[feat_i]
                if (inherits(use_values, "data.table")) {
                    feat_vals <- eval(call("[", use_values, j = as.name(feat)))
                } else {
                    feat_vals <- use_values[, feat]
                }

                spat_ac <- .spat_autocor_terra_numeric(
                    x = feat_vals,
                    w = weight_matrix,
                    method = method
                )

                # increment progress
                if (exists("pb")) if (feat_i %% step_size == 0) pb()


                out_dt <- data.table::data.table(spat_ac)
                colnames(out_dt) <- feat
                return(out_dt)
            }
        )
    })
    res_dt <- do.call("cbind", res_list)
    # append cell_ID column
    res_dt[, cell_ID := IDs]
    return(res_dt)
}



# Determine which information to retrieve and how to format the information
# Vars from upstream:
# use_sn - if true, extracts spatial network from gobject. Otherwise use
# externally provided info
# use_expr - if true, extracts expression information from gobject to use as
# node values
# use_meta - if true, extracts cell metadata information from gobject to use as
# node values
# use_ext_vals - directly use externally provided node value information

# Expected input:
# 1. source of data per spatial ID, whether that be expression information,
# cell metadata annotations, or external data
# 2. a spatial weight matrix for defining how important spatial interactions
# should
# be considered. This information can either be extracted spatial networks in
# the gobject with a pre-generated spatial weight matrix or generated during
# this call.

# Expected output:
# list of the following...
# 1. use_values - data per spatial ID. Formatted to be spatial ID (rows) by
# feats (cols)
# 2. feats - character vector of features in use_values to iterate through for
# autocor
# 3. weight_matrix - weight matrix (ordering checked to match with use_values
# if possible)
# 4, IDs - cell_IDs if available
# Some additional information about information used in specific workflows are
# also returned
.evaluate_autocor_input <- function(
        gobject,
        use_ext_vals,
        use_sn,
        use_expr,
        use_meta,
        spat_unit,
        feat_type,
        feats,
        data_to_use,
        expression_values,
        meta_cols,
        spatial_network_to_use,
        wm_method,
        wm_name,
        node_values,
        weight_matrix,
        verbose = TRUE) {
    cell_ID <- NULL

    # 1. Get spatial network to either get or generate a spatial weight matrix
    # End output is weight_matrix
    if (isTRUE(use_sn)) {
        # SPATNET==============================================================#
        sn <- getSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            name = spatial_network_to_use,
            output = "spatialNetworkObj"
        )
        weight_matrix <- slot(sn, "misc")$weight_matrix[[wm_name]]

        # if no weight_matrix already generated...
        if (is.null(weight_matrix)) {
            wm_method <- match.arg(
                wm_method,
                choices = c("distance", "adjacency")
            )
            if (isTRUE(verbose)) {
                wrap_msg(
                    "No spatial weight matrix found in selected spatial network
                    Generating", wm_method, "matrix from",
                    spatial_network_to_use
                )
            }
            weight_matrix <- createSpatialWeightMatrix(
                gobject = gobject,
                spat_unit = spat_unit,
                spatial_network_to_use = spatial_network_to_use,
                wm_name = wm_name,
                method = wm_method,
                return_gobject = FALSE,
                verbose = FALSE
            )
        }
        wm_colnames <- colnames(weight_matrix)
        # SPATNET==============================================================#
    }
    if (!isTRUE(use_sn)) {
        # EXTSPATNET===========================================================#
        if (!is.null(colnames(weight_matrix))) {
            wm_colnames <- colnames(weight_matrix)
            if (isTRUE(verbose)) {
                wrap_msg(
                    "colnames of externally provided weight matrix will be
                    matched to"
                )
            }
        }
        # EXTSPATNET===========================================================#
    }


    # 2. Get and format node values for use with autocorrelation function.
    # End outputs are:
    #  - use_values for a spatID (rows) by features (cols) table or matrix
    #  - feats the names of selected features to use that will be iterated
    # through downstream
    if (isTRUE(use_expr)) {
        # EXPR=================================================================#
        values <- match.arg(
            expression_values,
            unique(c("normalized", "scaled", "custom", expression_values))
        )
        use_values <- getExpression(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            values = values,
            output = "matrix"
        )
        use_values <- t_flex(use_values)

        # ensure identical ordering with giotto weight matrix
        if (exists("wm_colnames")) use_values <- use_values[wm_colnames, ]

        if (is.null(feats)) feats <- colnames(use_values)
        IDs <- rownames(use_values)
        # EXPR=================================================================#
    }
    if (isTRUE(use_meta)) {
        # META=================================================================#
        if (is.null(meta_cols)) {
            stop(wrap_txt("Metadata columns to autocorrelate must be given",
                errWidth = TRUE
            ))
        }
        use_values <- getCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "data.table",
            copy_obj = TRUE
        )

        # ensure identical ordering with giotto weight matrix
        if (exists("wm_colnames")) {
            new_order <- data.table::chmatch(wm_colnames, use_values$cell_ID)
            GiottoUtils::dt_set_row_order(use_values, new_order)
        }

        feats <- meta_cols
        IDs <- use_values[, cell_ID]
        # META=================================================================#
    }
    if (isTRUE(use_ext_vals)) {
        # EXTDATA==============================================================#
        use_values <- data.table::as.data.table(values = node_values)

        feats <- "values"
        # EXTDATA==============================================================#
    }


    # 3. general formatting and checking
    ## weight matrix type
    if (!inherits(weight_matrix, c("Matrix", "matrix", "listw", "nb"))) {
        stop(wrap_txt("weight_matrix must be a matrix, Matrix, or listw",
            errWidth = TRUE
        ))
    }

    ## check if weight matrix dimensions match use_values
    if ((nrow(use_values) != ncol(weight_matrix)) |
        (nrow(use_values) != nrow(weight_matrix))) {
        stop(wrap_txt("Number of values to correlate do not match number of
                      weight matrix entries",
            errWidth = TRUE
        ))
    }


    # return formatted values
    # provenance included if available
    return(list(
        use_values = use_values,
        feats = feats,
        weight_matrix = weight_matrix,
        # method specific items:
        expr_values = if (use_expr) values else NULL,
        provenance = if (use_sn) prov(sn) else NULL,
        IDs = if (use_expr | use_meta) IDs else NULL
    ))
}









# * ####
## spatial deconvolution functions ####


#' @title enrich_deconvolution
#' @description Rui to fill in
#' @returns matrix
#' @keywords internal
enrich_deconvolution <- function(
        expr,
        log_expr,
        cluster_info,
        ct_exp,
        cutoff) {
    ##### generate enrich 0/1 matrix based on expression matrix
    ct_exp <- ct_exp[rowSums(ct_exp) > 0, , drop = FALSE]
    enrich_matrix <- matrix(0, nrow = dim(ct_exp)[1], ncol = dim(ct_exp)[2])
    rowmax_col <- Rfast::rowMaxs(ct_exp)
    for (i in seq_along(rowmax_col)) {
        enrich_matrix[i, rowmax_col[i]] <- 1
    }
    colsum_ct_binary <- colSums(enrich_matrix)
    for (i in seq_along(colsum_ct_binary)) {
        if (colsum_ct_binary[i] <= 2) {
            rank <- rank(-ct_exp[, i])
            enrich_matrix[rank <= 2, i] <- 1
        }
    }
    rownames(enrich_matrix) <- rownames(ct_exp)
    colnames(enrich_matrix) <- colnames(ct_exp)

    ##### page enrich
    enrich_result <- enrich_analysis(log_expr, enrich_matrix)
    ##### initialize dwls matrix
    dwls_results <- matrix(0, nrow = dim(enrich_matrix)[2], ncol = dim(expr)[2])
    rownames(dwls_results) <- colnames(enrich_matrix)
    colnames(dwls_results) <- colnames(expr)
    cluster_sort <- mixedsort(unique(cluster_info))
    cluster_info <- cluster_info
    for (i in seq_along(cluster_sort)) {
        cluster_i_enrich <- enrich_result[
            , which(cluster_info == cluster_sort[i])
        ]
        row_i_max <- Rfast::rowMaxs(cluster_i_enrich, value = TRUE)
        ct <- rownames(enrich_result)[which(row_i_max > cutoff)]
        if (length(ct) < 2) {
            sort_rank <- sort(row_i_max, decreasing = TRUE)
            ct <- rownames(enrich_result)[which(row_i_max >= sort_rank[2])]
        }
        ct_gene <- c()
        for (j in seq_along(ct)) {
            sig_gene_j <- rownames(enrich_matrix)[
                which(enrich_matrix[, ct[j]] == 1)
            ]
            ct_gene <- c(ct_gene, sig_gene_j)
        }
        uniq_ct_gene <- intersect(rownames(expr), unique(ct_gene))
        select_sig_exp <- ct_exp[uniq_ct_gene, ct, drop = FALSE]
        cluster_i_cell <- which(cluster_info == cluster_sort[i])
        cluster_cell_exp <- expr[uniq_ct_gene, cluster_i_cell, drop = FALSE]

        cluster_i_dwls <- optimize_deconvolute_dwls(
            cluster_cell_exp, select_sig_exp
        )
        dwls_results[ct, cluster_i_cell] <- cluster_i_dwls
    }
    ##### remove negative values
    for (i in dim(dwls_results)[1]) {
        negtive_index <- which(dwls_results[i, ] < 0)
        dwls_results[i, negtive_index] <- 0
        ####* fixed a typo, with "==" the negative values weren't removed
    }
    return(dwls_results)
}


#' @title spot_deconvolution
#' @description Rui to fill in
#' @returns matrix
#' @keywords internal
spot_deconvolution <- function(
        expr,
        cluster_info,
        ct_exp,
        binary_matrix) {
    ##### generate enrich 0/1 matrix based on expression matrix
    enrich_matrix <- matrix(0, nrow = dim(ct_exp)[1], ncol = dim(ct_exp)[2])
    rowmax_col <- Rfast::rowMaxs(ct_exp)
    for (i in seq_along(rowmax_col)) {
        enrich_matrix[i, rowmax_col[i]] <- 1
    }
    rownames(enrich_matrix) <- rownames(ct_exp)
    colnames(enrich_matrix) <- colnames(ct_exp)

    cluster_sort <- mixedsort(unique(cluster_info))
    #### initialize dwls matrix
    dwls_results <- matrix(0, nrow = dim(ct_exp)[2], ncol = dim(expr)[2])
    rownames(dwls_results) <- colnames(ct_exp)
    colnames(dwls_results) <- colnames(expr)

    for (i in seq_along(cluster_sort)) {
        cluster_i_matrix <- binary_matrix[
            , which(cluster_info == cluster_sort[i])
        ]
        row_i_max <- Rfast::rowMaxs(cluster_i_matrix, value = TRUE)
        ct_i <- rownames(cluster_i_matrix)[which(row_i_max == 1)]
        ######## calculate proportion based on binarized deconvolution
        ######## results at first step
        if (length(ct_i) == 1) {
            dwls_results[ct_i[1], which(cluster_info == cluster_sort[i])] == 1
        } else {
            ct_gene <- c()
            for (j in seq_along(ct_i)) {
                sig_gene_j <- rownames(enrich_matrix)[
                    which(enrich_matrix[, ct_i[j]] == 1)
                ]
                ct_gene <- c(ct_gene, sig_gene_j)
            }
            uniq_ct_gene <- intersect(rownames(expr), unique(ct_gene))
            select_sig_exp <- ct_exp[uniq_ct_gene, ct_i, drop = FALSE]
            cluster_i_cell <- which(cluster_info == cluster_sort[i])
            cluster_cell_exp <- expr[uniq_ct_gene, cluster_i_cell, drop = FALSE]
            ###### calculate
            ###### overlap signature with spatial genes
            all_exp <- Matrix::rowMeans(cluster_cell_exp)
            solution_all_exp <- solve_OLS_internal(select_sig_exp, all_exp)
            constant_J <- find_dampening_constant(
                select_sig_exp, all_exp, solution_all_exp
            )
            ###### deconvolution for each spot
            for (k in seq_len(dim(cluster_cell_exp)[2])) {
                B <- Matrix::as.matrix(cluster_cell_exp[, k])
                ct_spot_k <- rownames(cluster_i_matrix)[
                    which(cluster_i_matrix[, k] == 1)
                ]
                if (sum(B) == 0 || length(ct_spot_k) == 0) {
                    ####* must include the case where all genes are 0
                    dwls_results[, colnames(cluster_cell_exp)[k]] <- NA
                    ####* will produce NAs for some spots in the output
                    next ####* no need to look into this spot any more
                }
                if (length(ct_spot_k) == 1) {
                    dwls_results[
                        ct_spot_k[1], colnames(cluster_cell_exp)[k]
                    ] <- 1
                } else {
                    ct_k_gene <- c()
                    for (m in seq_along(ct_spot_k)) {
                        sig_gene_k <- rownames(enrich_matrix)[which(
                            enrich_matrix[, ct_spot_k[m]] == 1
                        )]
                        ct_k_gene <- c(ct_k_gene, sig_gene_k)
                    }
                    uniq_ct_k_gene <- intersect(
                        rownames(ct_exp), unique(ct_k_gene)
                    )
                    S_k <- Matrix::as.matrix(ct_exp[uniq_ct_k_gene, ct_spot_k])
                    if (sum(B[uniq_ct_k_gene, ]) == 0) {
                        ####* must include the case all genes are 0
                        dwls_results[, colnames(cluster_cell_exp)[k]] <- NA
                        ####* will produce NAs for some spots in the output
                    } else {
                        solDWLS <- optimize_solveDampenedWLS(S_k, B[
                            uniq_ct_k_gene,
                        ], constant_J)
                        dwls_results[
                            names(solDWLS),
                            colnames(cluster_cell_exp)[k]
                        ] <- solDWLS
                    }
                }
            }
        }
    }
    ##### remove negative values
    for (i in dim(dwls_results)[1]) {
        negtive_index <- which(dwls_results[i, ] < 0)
        dwls_results[i, negtive_index] <- 0
    }
    return(dwls_results)
}



#' @title cluster_enrich_analysis
#' @description Rui to fill in
#' @returns enrichment values
#' @keywords internal
cluster_enrich_analysis <- function(
        exp_matrix,
        cluster_info,
        enrich_sig_matrix) {
    uniq_cluster <- mixedsort(unique(cluster_info))
    if (length(uniq_cluster) == 1) {
        stop("Only one cluster identified, need at least two.")
    }
    cluster_exp <- NULL
    for (i in uniq_cluster) {
        cluster_exp <- cbind(
            cluster_exp,
            (apply(
                exp_matrix, 1,
                function(y) mean(y[which(cluster_info == i)])
            ))
        )
    }
    log_cluster_exp <- log2(cluster_exp + 1)
    colnames(log_cluster_exp) <- uniq_cluster
    cluster_enrich <- enrich_analysis(log_cluster_exp, enrich_sig_matrix)
    return(cluster_enrich)
}


#' @title enrich_analysis
#' @description Rui to fill in
#' @returns enrichment matrix
#' @keywords internal
enrich_analysis <- function(
        expr_values,
        sign_matrix) {
    # output enrichment
    # only continue with genes present in both datasets
    interGene <- intersect(rownames(sign_matrix), rownames(expr_values))
    filterSig <- sign_matrix[interGene, ]
    signames <- rownames(filterSig)[which(filterSig[, 1] == 1)]
    # calculate mean gene expression
    # mean_gene_expr = rowMeans(expr_values)
    mean_gene_expr <- log2(Matrix::rowMeans(2^expr_values - 1, dims = 1) + 1)
    geneFold <- expr_values - mean_gene_expr
    # calculate sample/spot mean and sd
    cellColMean <- apply(geneFold, 2, mean)
    cellColSd <- apply(geneFold, 2, stats::sd)

    # get enrichment scores
    enrichment <- matrix(
        data = NA, nrow = dim(filterSig)[2], ncol = length(cellColMean)
    )
    for (i in seq_len(dim(filterSig)[2])) {
        signames <- rownames(filterSig)[which(filterSig[, i] == 1)]
        sigColMean <- apply(geneFold[signames, ], 2, mean)
        m <- length(signames)
        vectorX <- NULL
        for (j in (seq_along(cellColMean))) {
            Sm <- sigColMean[j]
            u <- cellColMean[j]
            sigma <- cellColSd[j]
            zscore <- (Sm - u) * m^(1 / 2) / sigma
            vectorX <- append(vectorX, zscore)
        }
        enrichment[i, ] <- vectorX
    }
    rownames(enrichment) <- colnames(filterSig)
    colnames(enrichment) <- names(cellColMean)
    return(enrichment)
}



#' @title optimize_deconvolute_dwls
#' @description Rui to fill in
#' @returns matrix
#' @keywords internal
optimize_deconvolute_dwls <- function(
        exp,
        Signature) {
    ###### overlap signature with spatial genes
    Genes <- intersect(rownames(Signature), rownames(exp))
    S <- Signature[Genes, ]
    S <- Matrix::as.matrix(S)
    Bulk <- Matrix::as.matrix(exp)
    subBulk <- Bulk[Genes, ]
    allCounts_DWLS <- NULL
    all_exp <- Matrix::rowMeans(exp)

    solution_all_exp <- solve_OLS_internal(S, all_exp[Genes])

    constant_J <- find_dampening_constant(S, all_exp[Genes], solution_all_exp)
    for (j in seq_len(dim(subBulk)[2])) {
        B <- subBulk[, j]
        if (sum(B) > 0) {
            solDWLS <- optimize_solveDampenedWLS(S, B, constant_J)
        } else {
            # solDWLS <- rep(0, length(B))  ####*  wrong dimension,
            # causes warnings
            solDWLS <- rep(0, ncol(S)) ####* corrected dim
            # names(solDWLS) <- names(B)  ####*  wrong dim name
            names(solDWLS) <- colnames(S) ####* corrected dim name
        }
        allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
    }
    colnames(allCounts_DWLS) <- colnames(exp)
    return(allCounts_DWLS)
}


#' @title optimize_solveDampenedWLS
#' @returns numeric
#' @keywords internal
optimize_solveDampenedWLS <- function(
        S,
        B,
        constant_J) {
    # first solve OLS, use this solution to find a starting point for the
    # weights
    solution <- solve_OLS_internal(S, B)
    # now use dampened WLS, iterate weights until convergence
    iterations <- 0
    changes <- c()
    # find dampening constant for weights using cross-validation
    j <- constant_J
    change <- 1

    while (change > .01 & iterations < 1000) {
        newsolution <- solve_dampened_WLSj(S, B, solution, j)
        # decrease step size for convergence
        solutionAverage <- rowMeans(cbind(
            newsolution,
            matrix(solution,
                nrow = length(solution),
                ncol = 4
            )
        ))
        change <- norm(Matrix::as.matrix(solutionAverage - solution))
        solution <- solutionAverage
        iterations <- iterations + 1
        changes <- c(changes, change)
    }


    return(solution / sum(solution))
}


#' @title find_dampening_constant
#' @description find a dampening constant for the weights using cross-validation
#' @returns numeric
#' @keywords internal
find_dampening_constant <- function(
        S,
        B,
        goldStandard) {
    solutionsSd <- NULL

    # goldStandard is used to define the weights
    sol <- goldStandard
    ws <- as.vector((1 / (S %*% sol))^2)
    wsScaled <- ws / min(ws)
    wsScaledMinusInf <- wsScaled

    # ignore infinite weights
    if (max(wsScaled) == "Inf") {
        wsScaledMinusInf <- wsScaled[-which(wsScaled == "Inf")]
    }

    # try multiple values of the dampening constant (multiplier)
    # for each, calculate the variance of the dampened weighted solution for
    # a subset of genes
    for (j in seq_len(ceiling(log2(max(wsScaledMinusInf))))) {
        multiplier <- 1 * 2^(j - 1)
        wsDampened <- wsScaled
        wsDampened[which(wsScaled > multiplier)] <- multiplier
        solutions <- NULL
        seeds <- seq_len(100)
        for (i in seq_len(100)) {
            set.seed(seeds[i]) # make nondeterministic
            subset <- sample(length(ws), size = length(ws) * 0.5)
            # randomly select half of gene set

            # solve dampened weighted least squares for subset
            fit <- stats::lm(
                B[subset] ~ -1 + S[subset, , drop = FALSE],
                weights = wsDampened[subset]
            )
            sol <- fit$coef * sum(goldStandard) / sum(fit$coef)
            solutions <- cbind(solutions, sol)
        }
        solutionsSd <- cbind(solutionsSd, apply(solutions, 1, stats::sd))
    }

    # choose dampening constant that results in least cross-validation variance
    j <- which.min(colMeans(solutionsSd^2))
    return(j)
}


#' @title solve_OLS_internal
#' @description basic functions for dwls
#' @returns numeric
#' @keywords internal
solve_OLS_internal <- function(
        S,
        B) {
    D <- t(S) %*% S
    d <- t(S) %*% B
    A <- cbind(diag(dim(S)[2]))
    bzero <- c(rep(0, dim(S)[2]))


    out <- tryCatch(
        expr = {
            quadprog::solve.QP(
                Dmat = D,
                dvec = d,
                Amat = A,
                bvec = bzero
            )$solution
        },
        error = function(cond) {
            message("Original error message:")
            message(cond)

            message("Will try to fix error with Matrix::nearPD()")
            return(NULL)
        }
    )


    if (is.null(out)) {
        D <- Matrix::nearPD(D)
        D <- as.matrix(D$mat)

        out <- tryCatch(
            expr = {
                quadprog::solve.QP(
                    Dmat = D,
                    dvec = d,
                    Amat = A,
                    bvec = bzero
                )$solution
            },
            error = function(cond) {
                message("Original error message:")
                message(cond)

                message("nearPD() did not fix the error")
                return(NULL)
            }
        )

        if (is.null(out)) {
            stop("Errors could not be fixed")
        } else {
            names(out) <- colnames(S)
        }
    } else {
        names(out) <- colnames(S)
    }

    return(out)
}

#
#' @title solve_dampened_WLSj
#' @description solve WLS given a dampening constant
#' @returns matrix
#' @keywords internal
solve_dampened_WLSj <- function(
        S,
        B,
        goldStandard,
        j) {
    multiplier <- 1 * 2^(j - 1)
    sol <- goldStandard
    ws <- as.vector((1 / (S %*% sol))^2)
    wsScaled <- ws / min(ws)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    W <- diag(wsDampened)
    D <- t(S) %*% W %*% S
    d <- t(S) %*% W %*% B
    A <- cbind(diag(dim(S)[2]))
    bzero <- c(rep(0, dim(S)[2]))
    sc <- norm(D, "2")

    D_positive_definite <- Matrix::nearPD(x = D / sc)

    solution <- quadprog::solve.QP(
        Dmat = as.matrix(D_positive_definite$mat),
        dvec = d / sc,
        Amat = A,
        bvec = bzero
    )$solution

    names(solution) <- colnames(S)
    return(solution)
}


#' @title runDWLSDeconv
#' @description Function to perform DWLS deconvolution based on single cell
#' expression data
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param logbase base used for log normalization
#' @param cluster_column name of cluster column
#' @param sign_matrix sig matrix for deconvolution
#' @param n_cell number of cells per spot
#' @param cutoff cut off (default = 2)
#' @param name name to give to spatial deconvolution results, default = DWLS
#' @param return_gobject return giotto object
#' @returns giotto object or deconvolution results
#' @seealso \url{https://github.com/dtsoucas/DWLS} for the \emph{DWLS} bulk
#' deconvolution method, and \doi{10.1186/s13059-021-02362-7} for
#' \emph{spatialDWLS}, the spatial implementation used here.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- findMarkers_one_vs_all(g,
#'     cluster_column = "leiden_clus", min_feats = 20
#' )
#' sign_gene <- x$feats
#'
#' sign_matrix <- matrix(rnorm(length(sign_gene) * 8, mean = 10),
#'     nrow = length(sign_gene)
#' )
#' rownames(sign_matrix) <- sign_gene
#' colnames(sign_matrix) <- paste0("celltype_", unique(x$cluster))
#'
#' runDWLSDeconv(gobject = g, sign_matrix = sign_matrix)
#' @export
runDWLSDeconv <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized"),
        logbase = 2,
        cluster_column = "leiden_clus",
        sign_matrix,
        n_cell = 50,
        cutoff = 2,
        name = NULL,
        return_gobject = TRUE) {
    # verify if optional package is installed
    package_check(pkg_name = "quadprog", repository = "CRAN")
    package_check(pkg_name = "Rfast", repository = "CRAN")

    ## check parameters ##
    if (is.null(name)) name <- "DWLS"


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

    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )

    # #transform expression data to no log data
    nolog_expr <- logbase^(expr_values[]) - 1

    # cluster column
    cell_metadata <- pDataDT(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    if (!cluster_column %in% colnames(cell_metadata)) {
        stop("cluster column not found")
    }
    cluster <- cell_metadata[[cluster_column]]


    ##### getting overlapped gene lists
    sign_matrix <- as.matrix(sign_matrix)
    intersect_gene <- intersect(rownames(sign_matrix), rownames(nolog_expr))
    filter_Sig <- sign_matrix[intersect_gene, ]
    filter_expr <- nolog_expr[intersect_gene, ]
    filter_log_expr <- expr_values[][intersect_gene, ]

    ##### first round spatial deconvolution ##spot or cluster
    enrich_spot_proportion <- enrich_deconvolution(
        expr = filter_expr,
        log_expr = filter_log_expr,
        cluster_info = cluster,
        ct_exp = filter_Sig,
        cutoff = cutoff
    )

    #### re-deconvolution based on spatial resolution
    resolution <- (1 / n_cell)
    binarize_proportion <- ifelse(enrich_spot_proportion >= resolution, 1, 0)
    spot_proportion <- spot_deconvolution(
        expr = filter_expr,
        cluster_info = cluster,
        ct_exp = filter_Sig,
        binary_matrix = binarize_proportion
    )
    deconvolutionDT <- data.table::data.table(
        cell_ID = colnames(spot_proportion)
    )
    deconvolutionDT <- cbind(
        deconvolutionDT, data.table::as.data.table(t(spot_proportion))
    )

    # create spatial enrichment object
    enrObj <- create_spat_enr_obj(
        name = name,
        method = "DWLS",
        enrichDT = deconvolutionDT,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = expr_values@provenance,
        misc = list(
            expr_values_used = expression_values,
            logbase = logbase,
            cluster_column_used = cluster_column,
            number_of_cells_per_spot = n_cell,
            used_cut_off = cutoff
        )
    )

    ## return object or results ##
    if (return_gobject == TRUE) {
        spenr_names <- list_spatial_enrichments_names(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        )

        if (name %in% spenr_names) {
            cat(name, " has already been used, will be overwritten")
        }

        ## update parameters used ##
        parameters_list <- gobject@parameters
        number_of_rounds <- length(parameters_list)
        update_name <- paste0(number_of_rounds, "_spatial_deconvolution")

        # parameters to include
        parameters_list[[update_name]] <- c(
            "method used" = "DWLS",
            "deconvolution name" = name,
            "expression values" = expression_values,
            "logbase" = logbase,
            "cluster column used" = cluster_column,
            "number of cells per spot" = n_cell,
            "used cut off" = cutoff
        )

        gobject@parameters <- parameters_list

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, enrObj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        return(gobject)
    } else {
        return(enrObj)
    }
}



#' @title runSpatialDeconv
#' @name runSpatialDeconv
#' @description Function to perform deconvolution based on single cell
#' expression data
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param deconv_method method to use for deconvolution
#' @param expression_values expression values to use
#' @param logbase base used for log normalization
#' @param cluster_column name of cluster column
#' @param sign_matrix signature matrix for deconvolution
#' @param n_cell number of cells per spot
#' @param cutoff cut off (default = 2)
#' @param name name to give to spatial deconvolution results
#' @param return_gobject return giotto object
#' @returns giotto object or deconvolution results
#' @seealso \code{\link{runDWLSDeconv}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- findMarkers_one_vs_all(g,
#'     cluster_column = "leiden_clus", min_feats = 20
#' )
#' sign_gene <- x$feats
#'
#' sign_matrix <- matrix(rnorm(length(sign_gene) * 8, mean = 10),
#'     nrow = length(sign_gene)
#' )
#' rownames(sign_matrix) <- sign_gene
#' colnames(sign_matrix) <- paste0("celltype_", unique(x$cluster))
#'
#' runSpatialDeconv(gobject = g, sign_matrix = sign_matrix)
#' @export
runSpatialDeconv <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        deconv_method = c("DWLS"),
        expression_values = c("normalized"),
        logbase = 2,
        cluster_column = "leiden_clus",
        sign_matrix,
        n_cell = 50,
        cutoff = 2,
        name = NULL,
        return_gobject = TRUE) {
    deconv_method <- match.arg(deconv_method, choices = c("DWLS"))


    if (deconv_method == "DWLS") {
        results <- runDWLSDeconv(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            expression_values = expression_values,
            logbase = logbase,
            cluster_column = cluster_column,
            sign_matrix = sign_matrix,
            n_cell = n_cell,
            cutoff = cutoff,
            name = name,
            return_gobject = return_gobject
        )
    }

    return(results)
}

# General Giotto Helper Functions

# statistical helpers belong here as opposed to GiottoUtils


#' @title mygini_fun
#' @description calculate gini coefficient
#' @keywords internal
#' @returns gini coefficient
mygini_fun <- function(
        x,
        weights = rep(1, length(x))) {
    # adapted from R package GiniWegNeg
    dataset <- cbind(x, weights)
    ord_x <- order(x)
    dataset_ord <- dataset[ord_x, ]
    x <- dataset_ord[, 1]
    weights <- dataset_ord[, 2]
    N <- sum(weights)
    xw <- x * weights
    C_i <- cumsum(weights)
    num_1 <- sum(xw * C_i)
    num_2 <- sum(xw)
    num_3 <- sum(xw * weights)
    G_num <- (2 / N^2) * num_1 - (1 / N) * num_2 - (1 / N^2) * num_3
    t_neg <- subset(xw, xw <= 0)
    T_neg <- sum(t_neg)
    T_pos <- sum(xw) + abs(T_neg)
    n_RSV <- (2 * (T_pos + (abs(T_neg))) / N)
    mean_RSV <- (n_RSV / 2)
    G_RSV <- (1 / mean_RSV) * G_num
    return(G_RSV)
}


#' @title extended_gini_fun
#' @description calculate gini coefficient on a minimum length vector
#' @keywords internal
#' @returns gini coefficient
extended_gini_fun <- function(
        x,
        weights = rep(1, length = length(x)),
        minimum_length = 16) {
    if (length(x) < minimum_length) {
        difference <- minimum_length - length(x)
        min_value <- min(x)
        x <- c(x, rep(min_value, difference))
    }

    result <- mygini_fun(x = x, weights = weights)
    return(result)
}


# matrix binarization methods ####

#' @title .kmeans_binarize
#' @name .kmeans_binarize
#' @description create binarized scores from a vector using kmeans
#' @returns numeric
#' @keywords internal
.kmeans_binarize <- function(
        x,
        nstart = 3,
        iter.max = 10,
        seed = NULL) {
    if (!is.null(seed)) {
        on.exit(random_seed(), add = TRUE)
        set.seed(seed)
    }
    sel_gene_km <- stats::kmeans(
        x,
        centers = 2, nstart = nstart, iter.max = iter.max
    )$cluster
    mean_1 <- mean(x[sel_gene_km == 1])
    mean_2 <- mean(x[sel_gene_km == 2])

    if (mean_1 > mean_2) {
        mean_1_value <- 1
        mean_2_value <- 0
    } else {
        mean_1_value <- 0
        mean_2_value <- 1
    }

    sel_gene_bin <- x
    sel_gene_bin[sel_gene_km == 1] <- mean_1_value
    sel_gene_bin[sel_gene_km == 2] <- mean_2_value

    return(sel_gene_bin)
}

#' @title .kmeans_arma_binarize
#' @name .kmeans_arma_binarize
#' @description create binarized scores from a vector using kmeans_arma
#' @returns numeric
#' @keywords internal
.kmeans_arma_binarize <- function(x, n_iter = 5, seed = NULL) {
    if (!is.null(seed)) {
        on.exit(random_seed(), add = TRUE)
        set.seed(seed)
    }
    sel_gene_km_res <- ClusterR::KMeans_arma(
        data = as.matrix(x),
        clusters = 2,
        n_iter = n_iter
    )
    sel_gene_km <- ClusterR::predict_KMeans(
        data = as.matrix(x),
        CENTROIDS = sel_gene_km_res
    )

    mean_1 <- mean(x[sel_gene_km == 1])
    mean_2 <- mean(x[sel_gene_km == 2])

    if (mean_1 > mean_2) {
        mean_1_value <- 1
        mean_2_value <- 0
    } else {
        mean_1_value <- 0
        mean_2_value <- 1
    }

    sel_gene_bin <- x
    sel_gene_bin[sel_gene_km == 1] <- mean_1_value
    sel_gene_bin[sel_gene_km == 2] <- mean_2_value

    return(sel_gene_bin)
}

#' @title .kmeans_arma_subset_binarize
#' @name .kmeans_arma_subset_binarize
#' @description create binarized scores from a subsetted vector using
#' kmeans_arma
#' @returns numeric
#' @keywords internal
.kmeans_arma_subset_binarize <- function(
        x,
        n_iter = 5,
        extreme_nr = 20,
        sample_nr = 200,
        seed = NULL) {
    length_x <- length(x)

    vector_x <- sort(x)
    first_set <- vector_x[seq_len(extreme_nr)]
    last_set <- vector_x[(length_x - (extreme_nr - 1)):length_x]
    random_set <- sample(
        vector_x[(extreme_nr + 1):(length_x - extreme_nr)],
        size = sample_nr
    )
    testset <- c(first_set, last_set, random_set)

    if (!is.null(seed)) {
        on.exit(random_seed(), add = TRUE)
        set.seed(seed)
    }
    sel_gene_km_res <- ClusterR::KMeans_arma(
        data = as.matrix(testset),
        clusters = 2,
        n_iter = n_iter
    )
    sel_gene_km <- ClusterR::predict_KMeans(
        data = as.matrix(x),
        CENTROIDS = sel_gene_km_res
    )

    mean_1 <- mean(x[sel_gene_km == 1])
    mean_2 <- mean(x[sel_gene_km == 2])

    if (mean_1 > mean_2) {
        mean_1_value <- 1
        mean_2_value <- 0
    } else {
        mean_1_value <- 0
        mean_2_value <- 1
    }

    sel_gene_bin <- x
    sel_gene_bin[sel_gene_km == 1] <- mean_1_value
    sel_gene_bin[sel_gene_km == 2] <- mean_2_value

    return(sel_gene_bin)
}



#' @title kmeans_binarize_wrapper
#' @name kmeans_binarize_wrapper
#' @description wrapper for different binarization functions
#' @returns matrix
#' @keywords internal
kmeans_binarize_wrapper <- function(expr_values,
    subset_feats = NULL,
    kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
    nstart = 3,
    iter_max = 10,
    extreme_nr = 50,
    sample_nr = 50,
    seed = NULL) {
    # expression values
    if (!is.null(subset_feats)) {
        expr_values <- expr_values[rownames(expr_values) %in% subset_feats, ]
    }

    # check parameter
    kmeans_algo <- match.arg(
        arg = kmeans_algo,
        choices = c("kmeans", "kmeans_arma", "kmeans_arma_subset")
    )

    bin_matrix <- switch(kmeans_algo,
        "kmeans" = t_flex(apply(
            X = expr_values, MARGIN = 1, FUN = .kmeans_binarize,
            nstart = nstart, iter.max = iter_max, seed = seed
        )),
        "kmeans_arma" = t_flex(apply(
            X = expr_values, MARGIN = 1, FUN = .kmeans_arma_binarize,
            n_iter = iter_max, seed = seed
        )),
        "kmeans_arma_subset" = t_flex(apply(
            X = expr_values, MARGIN = 1, FUN = .kmeans_arma_subset_binarize,
            n_iter = iter_max,
            extreme_nr = extreme_nr,
            sample_nr = sample_nr,
            seed = seed
        ))
    )

    return(bin_matrix)
}



#' @title Rank binarize
#' @name .rank_binarize
#' @description create binarized scores from a vector using arbitrary rank
#' @returns data.table
#' @keywords internal
.rank_binarize <- function(x, max_rank = 200) {
    sel_gene_rank <- rank(-x, ties.method = "average")

    sel_gene_bin <- x
    sel_gene_bin[sel_gene_rank <= max_rank] <- 1
    sel_gene_bin[sel_gene_rank > max_rank] <- 0

    return(sel_gene_bin)
}



#' @title rank_binarize_wrapper
#' @name rank_binarize_wrapper
#' @description wrapper for rank binarization function
#' @returns matrix
#' @keywords internal
rank_binarize_wrapper <- function(
        expr_values,
        subset_feats = NULL,
        percentage_rank = 30) {
    # expression values
    if (!is.null(subset_feats)) {
        expr_values <- expr_values[rownames(expr_values) %in% subset_feats, ]
    }

    max_rank <- (ncol(expr_values) / 100) * percentage_rank
    bin_matrix <- t_flex(apply(
        X = expr_values, MARGIN = 1, FUN = .rank_binarize, max_rank = max_rank
    ))

    return(bin_matrix)
}


## chatgpt queries ####

#' @title writeChatGPTqueryDEG
#' @name writeChatGPTqueryDEG
#' @description This function writes a query as a .txt file that can be used
#' with ChatGPT or a similar LLM service to find the most likely cell types
#' based on the top differential expressed genes (DEGs) between identified
#' clusters.
#' @param DEG_output the output format from the differential expression
#' functions
#' @param top_n_genes number of genes for each cluster
#' @param tissue_type tissue type
#' @param folder_name path to the folder where you want to save the .txt file
#' @param file_name name of .txt file
#' @returns writes a .txt file to the desired location
#' @details This function does not run any LLM service. It simply creates the
#' .txt file that can then be used any LLM service (e.g. OpenAI, Gemini, ...)
#' @export
writeChatGPTqueryDEG <- function(
        DEG_output,
        top_n_genes = 10,
        tissue_type = "human breast cancer",
        folder_name = getwd(),
        file_name = "chatgpt_query.txt") {
    chatgpt_query <- paste0(
        "Identify cell types of ", tissue_type,
        " tissue using the following markers. Identify one cell type for each
        row. Only provide the cell type name and the marker genes used for cell
        type identification."
    )

    selected_DEG_output <- DEG_output[, head(.SD, top_n_genes), by = "cluster"]

    finallist <- list()
    finallist[[1]] <- chatgpt_query

    for (clus in unique(selected_DEG_output$cluster)) {
        x <- selected_DEG_output[cluster == clus][["feats"]]
        x <- c(clus, x)
        finallist[[as.numeric(clus) + 1]] <- x
    }

    outputdt <- data.table::data.table(finallist)

    cat("\n start writing \n")
    data.table::fwrite(
        x = outputdt,
        file = paste0(folder_name, "/", file_name),
        sep2 = c("", " ", ""), col.names = FALSE
    )
}



# IDs ####


#' @title convertEnsemblToGeneSymbol
#' @name convertEnsemblToGeneSymbol
#' @description This function convert ensembl gene IDs from a matrix to
#' official gene symbols
#' @param matrix an expression matrix with ensembl gene IDs as rownames
#' @param species species to use for gene symbol conversion
#' @returns expression matrix with gene symbols as rownames
#' @details This function requires that the biomaRt library is installed
#' @export
convertEnsemblToGeneSymbol <- function(
        matrix,
        species = c("mouse", "human")) {
    # data.table: set global variable
    dupes <- mgi_symbol <- gene_symbol <- ensembl_gene_id <- hgnc_symbol <- NULL

    package_check("biomaRt", repository = "Bioc")

    species <- match.arg(species, choices = c("mouse", "human"))

    if (species == "mouse") {
        # ensembl IDs to change
        ensemblsIDS <- rownames(matrix)

        # prepare ensembl database
        ensembl <- biomaRt::useMart("ensembl",
            dataset = "mmusculus_gene_ensembl"
        )
        gene_names <- biomaRt::getBM(
            attributes = c("mgi_symbol", "ensembl_gene_id"),
            filters = "ensembl_gene_id",
            values = ensemblsIDS,
            mart = ensembl
        )
        gene_names_DT <- data.table::as.data.table(gene_names)
        gene_names_DT[, dupes := duplicated(mgi_symbol)]
        gene_names_DT[, gene_symbol := ifelse(any(dupes) == FALSE, mgi_symbol,
            ifelse(mgi_symbol == "", ensembl_gene_id, "temporary")
        ), by = mgi_symbol]
        gene_names_DT[, gene_symbol := ifelse(
            mgi_symbol == "", ensembl_gene_id, gene_symbol
        )]
        gene_names_DT[, gene_symbol := ifelse(
            gene_symbol == "temporary",
            paste0(mgi_symbol, "--", seq_len(.N)), gene_symbol
        ),
        by = mgi_symbol
        ]

        # filter
        matrix <- matrix[rownames(matrix) %in% gene_names_DT$ensembl_gene_id, ]

        # create swapping vector
        new_symbols <- gene_names_DT[["gene_symbol"]]
        names(new_symbols) <- gene_names_DT[["ensembl_gene_id"]]

        # replace
        new_rownames <- new_symbols[rownames(matrix)]
        rownames(matrix) <- new_rownames

        return(matrix)
    }

    if (species == "human") {
        # ensembl IDs to change
        ensemblsIDS <- rownames(matrix)

        # prepare ensembl database
        ensembl <- biomaRt::useMart("ensembl",
            dataset = "hsapiens_gene_ensembl"
        )
        gene_names <- biomaRt::getBM(
            attributes = c("hgnc_symbol", "ensembl_gene_id"),
            filters = "ensembl_gene_id",
            values = ensemblsIDS,
            mart = ensembl
        )
        gene_names_DT <- data.table::as.data.table(gene_names)
        gene_names_DT[, dupes := duplicated(hgnc_symbol)]
        gene_names_DT[, gene_symbol := ifelse(any(dupes) == FALSE, hgnc_symbol,
            ifelse(hgnc_symbol == "", ensembl_gene_id, "temporary")
        ), by = hgnc_symbol]
        gene_names_DT[, gene_symbol := ifelse(
            hgnc_symbol == "", ensembl_gene_id, gene_symbol
        )]
        gene_names_DT[, gene_symbol := ifelse(
            gene_symbol == "temporary",
            paste0(hgnc_symbol, "--", seq_len(.N)), gene_symbol
        ),
        by = hgnc_symbol
        ]

        # filter
        matrix <- matrix[rownames(matrix) %in% gene_names_DT$ensembl_gene_id, ]

        # create swapping vector
        new_symbols <- gene_names_DT[["gene_symbol"]]
        names(new_symbols) <- gene_names_DT[["ensembl_gene_id"]]

        # replace
        new_rownames <- new_symbols[rownames(matrix)]
        rownames(matrix) <- new_rownames

        return(matrix)
    }
}






# I/O Helpers ####



## Parallelized giottoPolygon creation workflows ####

# Internal function to create a giottoPolygon object, smooth it, then wrap it so
# that results are portable/possible to use with parallelization.
# dotparams are passed to smoothGiottoPolygons
#' @title Polygon creation and smoothing for parallel
#' @name gpoly_from_dfr_smoothed_wrapped
#' @returns giottoPolygon
#' @keywords internal
gpoly_from_dfr_smoothed_wrapped <- function(segmdfr,
    name = "cell",
    calc_centroids = FALSE,
    smooth_polygons = FALSE,
    vertices = 20L,
    k = 3L,
    set_neg_to_zero = TRUE,
    skip_eval_dfr = FALSE,
    copy_dt = TRUE,
    verbose = TRUE) {
    gpoly <- createGiottoPolygonsFromDfr(
        segmdfr = segmdfr,
        name = name,
        calc_centroids = FALSE,
        skip_eval_dfr = skip_eval_dfr,
        copy_dt = copy_dt,
        verbose = verbose
    )
    if (isTRUE(smooth_polygons)) {
        gpoly <- smoothGiottoPolygons(
            gpolygon = gpoly,
            vertices = vertices,
            k = k,
            set_neg_to_zero = set_neg_to_zero
        )
    }
    if (isTRUE(calc_centroids)) {
        gpoly <- centroids(
            gpoly,
            append_gpolygon = TRUE
        )
    }

    slot(gpoly, "spatVector") <- terra::wrap(slot(gpoly, "spatVector"))
    if (isTRUE(calc_centroids)) {
        slot(gpoly, "spatVectorCentroids") <- terra::wrap(
            slot(gpoly, "spatVectorCentroids")
        )
    }
    return(gpoly)
}



## 10X ####

#' @title get10Xmatrix
#' @name get10Xmatrix
#' @description This function creates an expression matrix from a 10X
#' structured folder
#' @param path_to_data path to the 10X folder
#' @param gene_column_index which column from the features or genes .tsv file
#' to use for row ids
#' @param remove_zero_rows removes rows with sum equal to zero
#' @param split_by_type split into multiple matrices based on 3rd column of
#' features.tsv(.gz)
#' @returns sparse expression matrix from 10X
#' @details A typical 10X folder is named raw_feature_bc_matrix or
#' filtered_feature_bc_matrix and it has 3 files:
#' \itemize{
#'   \item{barcodes.tsv(.gz)}
#'   \item{features.tsv(.gz) or genes.tsv(.gz)}
#'   \item{matrix.mtx(.gz)}
#' }
#' By default the first column of the features or genes .tsv file will be used,
#' however if multiple
#' annotations are provided (e.g. ensembl gene ids and gene symbols) the user
#' can select another column.
#' @export
get10Xmatrix <- function(
        path_to_data,
        gene_column_index = 1,
        remove_zero_rows = TRUE,
        split_by_type = TRUE) {
    # data.table variables
    total <- gene_symbol <- gene_id <- gene_id_num <- cell_id <-
        cell_id_num <- sort_gene_id_num <- NULL

    # data directory
    files_10X <- list.files(path_to_data)

    # get barcodes and create vector
    barcodes_file <- grep(files_10X, pattern = "barcodes", value = TRUE)
    barcodesDT <- data.table::fread(
        input = paste0(path_to_data, "/", barcodes_file), header = FALSE
    )
    barcodes_vec <- barcodesDT$V1
    names(barcodes_vec) <- seq_len(nrow(barcodesDT))

    # get features and create vector
    features_file <- grep(files_10X, pattern = "features|genes", value = TRUE)
    featuresDT <- data.table::fread(
        input = paste0(path_to_data, "/", features_file), header = FALSE
    )

    g_name <- colnames(featuresDT)[gene_column_index]
    ## convert ensembl gene id to gene symbol ##
    ## TODO

    featuresDT[, total := .N, by = get(g_name)]
    featuresDT[, gene_symbol := ifelse(
        total > 1, paste0(get(g_name), "--", seq_len(.N)),
        get(g_name)
    ), by = get(g_name)]
    features_vec <- featuresDT$gene_symbol
    names(features_vec) <- seq_len(nrow(featuresDT))

    # get matrix
    matrix_file <- grep(files_10X, pattern = "matrix", value = TRUE)
    MMmatrix <- Matrix::readMM(paste0(path_to_data, "/", matrix_file))
    rownames(MMmatrix) <- features_vec
    colnames(MMmatrix) <- barcodes_vec


    # Split by type of feature (features.tsv 3rd col)
    feat_classes <- unique(featuresDT$V3)
    if (length(feat_classes) > 1 && isTRUE(split_by_type)) {
        result_list <- list()

        for (fclass in feat_classes) {
            result_list[[fclass]] <- MMmatrix[featuresDT$V3 == fclass, ]
        }

        if (isTRUE(remove_zero_rows)) {
            result_list <- lapply(result_list, function(MMmatrix) {
                rowsums_result <- rowSums_flex(MMmatrix)
                rowsums_bool <- rowsums_result != 0
                MMmatrix <- MMmatrix[rowsums_bool, ]
            })
        }

        return(result_list)
    } else {
        if (remove_zero_rows == TRUE) {
            rowsums_result <- rowSums_flex(MMmatrix)
            rowsums_bool <- rowsums_result != 0
            MMmatrix <- MMmatrix[rowsums_bool, ]
        }

        return(MMmatrix)
    }
}





#' @title get10Xmatrix_h5
#' @name get10Xmatrix_h5
#' @description This function creates an expression matrix from a 10X h5 file
#' path
#' @param path_to_data path to the 10X .h5 file
#' @param gene_ids use gene symbols (default) or ensembl ids for the gene
#' expression matrix
#' @inheritParams get10Xmatrix
#' @returns (list of) sparse expression matrix from 10X
#' @details If the .h5 10x file has multiple classes of features
#' (e.g. expression vs QC probes) or modalities (e.g. RNA and protein), and
#' \code{split_by_type} param is \code{TRUE}, multiple matrices will be returned
#' @export
get10Xmatrix_h5 <- function(
        path_to_data,
        gene_ids = c("symbols", "ensembl"),
        remove_zero_rows = TRUE,
        split_by_type = TRUE) {
    ## function inspired by and modified from the VISION package
    ## see read_10x_h5_v3 in
    ## https://github.com/YosefLab/VISION/blob/master/R/Utilities.R

    # verify if optional package is installed
    package_check(pkg_name = "hdf5r", repository = "CRAN")

    # select parameter
    gene_ids <- match.arg(gene_ids, choices = c("symbols", "ensembl"))

    h5 <- hdf5r::H5File$new(path_to_data)

    tryCatch(
        {
            # list objects part of the h5 file
            # hdf5r::list.objects(h5)

            # get root folder name e.g. 'matrix'
            root <- names(h5)
            root <- root[1]

            # extraction information
            data <- h5[[paste0(root, "/data")]][]
            data <- as.numeric(data)

            barcodes <- h5[[paste0(root, "/barcodes")]][]
            feature_tag_keys <- h5[[paste0(root, "/features/_all_tag_keys")]][]
            feature_types <- h5[[paste0(root, "/features/feature_type")]][]
            genome <- unique(h5[[paste0(root, "/features/genome")]][])
            feature_id <- h5[[paste0(root, "/features/id")]][]
            feature_names <- h5[[paste0(root, "/features/name")]][]
            indices <- h5[[paste0(root, "/indices")]][]
            indptr <- h5[[paste0(root, "/indptr")]][]
            data_shape <- h5[[paste0(root, "/shape")]][]

            # create a feature data.table
            features_dt <- data.table::data.table(
                "id" = feature_id,
                "name" = feature_names,
                "feature_type" = feature_types,
                "genome" = genome
            )
        },
        finally = {
            h5$close_all()
        }
    )

    # create uniq name symbols
    # duplicate gene symbols will be given a suffix '_1', '_2', ...

    # data.table variables
    nr_name <- name <- uniq_name <- NULL

    features_dt[, nr_name := seq_len(.N), by = name]
    features_dt[, uniq_name := ifelse(
        nr_name == 1, name, paste0(name, "_", (nr_name - 1))
    )]


    # dimension names
    dimnames <- list(feature_id, barcodes)

    sparsemat <- Matrix::sparseMatrix(
        i = indices + 1,
        p = indptr,
        x = data,
        dims = data_shape,
        dimnames = dimnames
    )

    # multiple feature classes (e.g. gene vs diagnostic or even modalities?)
    if (isTRUE(split_by_type)) {
        result_list <- list()

        for (fclass in unique(feature_types)) {
            result_list[[fclass]] <- sparsemat[
                features_dt$feature_type == fclass,
            ]

            # change names to gene symbols if it's expression
            if (fclass == "Gene Expression" && gene_ids == "symbols") {
                conv_vector <- features_dt$uniq_name
                names(conv_vector) <- features_dt$id

                current_names <- rownames(result_list[[fclass]])
                new_names <- conv_vector[current_names]
                rownames(result_list[[fclass]]) <- new_names
            }
        }

        if (isTRUE(remove_zero_rows)) {
            result_list <- lapply(result_list, function(sparsemat) {
                rowsums_result <- rowSums_flex(sparsemat)
                rowsums_bool <- rowsums_result != 0
                sparsemat <- sparsemat[rowsums_bool, ]
            })
        }

        return(result_list)
    } else {
        if (isTRUE(remove_zero_rows)) {
            rowsums_result <- rowSums_flex(sparsemat)
            rowsums_bool <- rowsums_result != 0
            sparsemat <- sparsemat[rowsums_bool, ]
        }

        return(list("Gene Expression" = sparsemat))
    }
}

#' read10xAffineImage
#' @name read10xAffineImage
#' @description Read a 10x image that is provided with an affine matrix
#' transform. Loads the image in with an orientation that matches the dataset
#' points and polygons vector information
#' @param file filepath to image
#' @param imagealignment_path filepath to alignment file which contains
#' an affine transformation matrix. Usually a `.csv` file
#' @param name character. Name to assign the image. Default is 'image'.
#' @param micron micron scaling. Directly used if a numeric is supplied.
#' Also prefers a filepath to the `experiment.xenium` file which contains this
#' info. A default of 0.2125 is provided.
#' @param \dots additional params to pass to
#' `[GiottoClass::createGiottoLargeImage]`
#' @md
#' @returns 10xAffineImage
#' @export
read10xAffineImage <- function(file, imagealignment_path,
    name = "aligned_image", micron = 0.2125, ...) {
    checkmate::assert_file_exists(file)
    checkmate::assert_file_exists(imagealignment_path)
    if (!is.numeric(micron)) {
        checkmate::assert_file_exists(micron)
        micron <- read_json(micron)$pixel_size
    }

    aff <- data.table::fread(imagealignment_path) %>%
        as.matrix()

    img <- createGiottoLargeImage(file, name = name, ...)

    aff_img <- .tenx_img_affine(x = img, affine = aff, micron = micron)

    return(aff_img)
}

.tenx_img_affine <- function(x, affine, micron) {
    x %>%
        affine(affine[seq(2), seq(2)]) %>%
        rescale(micron, x0 = 0, y0 = 0) %>%
        spatShift(dx = affine[1, 3] * micron, dy = -affine[2, 3] * micron)
}





## Vizgen ####



#' @title readPolygonFilesVizgenHDF5
#' @name readPolygonFilesVizgenHDF5_old
#' @description Read and create polygons for all cells, or for only selected
#' FOVs.
#' @param boundaries_path path to the cell_boundaries folder
#' @param fovs subset of fovs to use
#' @param custom_polygon_names a character vector to provide custom polygon
#' names (optional)
#' @param polygon_feat_types a vector containing the polygon feature types
#' @param flip_x_axis flip x axis of polygon coordinates (multiply by -1)
#' @param flip_y_axis flip y axis of polygon coordinates (multiply by -1)
#' @param smooth_polygons smooth polygons (default = TRUE)
#' @param smooth_vertices number of vertices for smoothing
#' @param set_neg_to_zero set negative values to zero when smoothing
#' @param H5Fopen_flags see \code{\link[rhdf5]{H5Fopen}} for more details
#' @param cores cores to use
#' @param verbose be verbose
#' @seealso \code{\link{smoothGiottoPolygons}}
#' @returns data.table
#' @details Set H5Fopen_flags to "H5F_ACC_RDONLY" if you encounter permission
#' issues.
#' @export
readPolygonFilesVizgenHDF5_old <- function(
        boundaries_path,
        fovs = NULL,
        polygon_feat_types = 0:6,
        custom_polygon_names = NULL,
        flip_x_axis = FALSE,
        flip_y_axis = FALSE,
        smooth_polygons = TRUE,
        smooth_vertices = 60,
        set_neg_to_zero = FALSE,
        H5Fopen_flags = "H5F_ACC_RDWR",
        cores = NA,
        verbose = TRUE) {
    # necessary pkgs
    package_check(pkg_name = "rhdf5", repository = "Bioc")

    cores <- determine_cores(cores)

    # data.table vars
    x <- y <- cell_id <- file_id <- my_id <- NULL

    # prepare poly feat names
    poly_feat_names <- paste0("z", polygon_feat_types)
    poly_feat_indexes <- paste0("zIndex_", polygon_feat_types)

    # provide your own custom names
    if (!is.null(custom_polygon_names)) {
        if (!is.character(custom_polygon_names)) {
            stop(wrap_txt("If custom_polygon_names are provided, it needs to
                        be a character vector"))
        }

        if (length(custom_polygon_names) != length(poly_feat_names)) {
            stop(wrap_txt("length of custom names need to be same as
                        polygon_feat_types"))
        } else {
            poly_feat_names <- custom_polygon_names
        }
    }

    if (isTRUE(verbose)) wrap_msg("Reading from:", boundaries_path)
    # list all files in the folder
    hdf5_boundary_list <- list.files(full.names = TRUE, boundaries_path)
    # only load subset of files if fov is given
    if (!is.null(fovs)) {
        selected_hdf5s <- paste0("feature_data_", fovs, ".hdf5")
        selected_hdf5s_concatenated <- paste0(selected_hdf5s, collapse = "|")
        hdf5_boundary_selected_list <- grep(
            selected_hdf5s_concatenated,
            x = hdf5_boundary_list, value = TRUE
        )
    } else {
        hdf5_boundary_selected_list <- hdf5_boundary_list
    }

    if (isTRUE(verbose)) {
        wrap_msg("finished listing .hdf5 files start extracting .hdf5
                information")
    }

    # open selected polygon files
    hdf5_list_length <- length(hdf5_boundary_selected_list)

    # append data from all FOVs to single list
    init <- proc.time()
    with_pbar({
        pb <- pbar(along = hdf5_boundary_selected_list)
        read_list <- lapply_flex(seq_along(hdf5_boundary_selected_list),
            cores = cores,
            future.packages = c("rhdf5", "Rhdf5lib"),
            function(bound_i) {
                # get feature data
                read_file <- rhdf5::H5Fopen(
                    hdf5_boundary_selected_list[[bound_i]][[1]],
                    flags = H5Fopen_flags
                )
                fov_info <- read_file$featuredata

                # update progress
                if (verbose) {
                    print(basename(hdf5_boundary_selected_list[[bound_i]]))
                }
                elapsed <- (proc.time() - init)[[3L]]
                step_time <- elapsed / bound_i
                est <- (hdf5_list_length * step_time) - elapsed
                pb(message = c(
                    "// E:", time_format(elapsed), "| R:", time_format(est)
                ))
                rhdf5::H5Fclose(read_file)
                return(fov_info)
            }
        )
    })

    # # combine to FOV data single list
    read_list <- Reduce("append", read_list)
    cell_names <- names(read_list)


    # extract values for each z index and cell from read_list
    result_list <- lapply_flex(
        seq_along(poly_feat_indexes),
        cores = cores, function(z_i) {
            lapply_flex(seq_along(read_list), cores = cores, function(cell_i) {
                singlearray <- read_list[[cell_i]][[
                    poly_feat_indexes[z_i]
                ]]$p_0$coordinates
                cell_name <- cell_names[[cell_i]]
                if (!is.null(singlearray)) {
                    singlearraydt <- data.table::as.data.table(t_flex(
                        as.matrix(singlearray[, , 1])
                    ))
                    data.table::setnames(
                        singlearraydt,
                        old = c("V1", "V2"), new = c("x", "y")
                    )
                    if (flip_x_axis) singlearraydt[, x := -1 * x]
                    if (flip_y_axis) singlearraydt[, y := -1 * y]

                    singlearraydt[, cell_id := cell_name]
                }
            })
        }
    )
    result_list_rbind <- lapply_flex(
        seq_along(result_list),
        cores = cores, function(z_i) {
            data.table::rbindlist(result_list[[z_i]])
        }
    )



    if (isTRUE(verbose)) {
        wrap_msg("finished extracting .hdf5 files start creating polygons")
    }


    # create Giotto polygons and add them to gobject
    with_pbar({
        pb <- pbar(along = result_list_rbind)
        smooth_cell_polygons_list <- lapply_flex(seq_along(result_list_rbind),
            cores = cores, function(i) {
                dfr_subset <- result_list_rbind[[i]][, .(x, y, cell_id)]
                cell_polygons <- createGiottoPolygonsFromDfr(
                    segmdfr = dfr_subset,
                    name = poly_feat_names[i],
                    verbose = verbose
                )

                pb(message = poly_feat_names[i])

                if (smooth_polygons == TRUE) {
                    return(smoothGiottoPolygons(cell_polygons,
                        vertices = smooth_vertices,
                        set_neg_to_zero = set_neg_to_zero
                    ))
                } else {
                    return(cell_polygons)
                }
            }
        )
    })


    # TODO: add spatial centroids
    # needs to happen after smoothing to be correct

    return(smooth_cell_polygons_list)
}




#' @title readPolygonFilesVizgenHDF5
#' @name readPolygonFilesVizgenHDF5
#' @description Read polygon info for all cells or for only selected FOVs from
#' Vizgen HDF5 files. Data is returned as a list of giottoPolygons or
#' data.tables of the requested z indices.
#' @param boundaries_path path to the cell_boundaries folder
#' @param fovs subset of fovs to use
#' @param z_indices z indices of polygons to use
#' @param segm_to_use segmentation results to use (usually = 1. Depends on if
#' alternative segmentations were generated)
#' @param custom_polygon_names a character vector to provide custom polygon
#' names (optional)
#' @param polygon_feat_types deprecated. Use \code{z_indices}
#' @param flip_x_axis flip x axis of polygon coordinates (multiply by -1)
#' @param flip_y_axis flip y axis of polygon coordinates (multiply by -1)
#' @param smooth_polygons smooth polygons (default = TRUE)
#' @param smooth_vertices number of vertices for smoothing
#' @param set_neg_to_zero set negative values to zero when smoothing
#' @param calc_centroids calculate centroids (default = FALSE)
#' @param H5Fopen_flags see \code{\link[rhdf5]{H5Fopen}} for more details
#' @param cores cores to use
#' @param create_gpoly_parallel (default = TRUE) Whether to run gpoly creation
#' in parallel
#' @param create_gpoly_bin (Optional, default = FALSE) Parallelization option.
#' Accepts integer values as an binning size when generating giottoPolygon
#' objects
#' @param verbose be verbose
#' @param output whether to return as list of giottoPolygon or data.table
#' @seealso \code{\link{smoothGiottoPolygons}}
#' @returns list of giottoPolygon or data.table
#' @details Set H5Fopen_flags to "H5F_ACC_RDONLY" if you encounter permission
#' issues.
#' @export
readPolygonFilesVizgenHDF5 <- function(
        boundaries_path,
        fovs = NULL,
        z_indices = 1L:7L,
        segm_to_use = 1L,
        custom_polygon_names = NULL,
        flip_x_axis = FALSE,
        flip_y_axis = TRUE,
        calc_centroids = FALSE,
        smooth_polygons = TRUE,
        smooth_vertices = 60L,
        set_neg_to_zero = FALSE,
        H5Fopen_flags = "H5F_ACC_RDWR",
        cores = determine_cores(),
        create_gpoly_parallel = TRUE,
        create_gpoly_bin = FALSE,
        verbose = TRUE,
        output = c("giottoPolygon", "data.table"),
        polygon_feat_types = NULL) {
    # necessary pkgs
    package_check(pkg_name = "rhdf5", repository = "Bioc")

    output <- match.arg(output, choices = c("giottoPolygon", "data.table"))

    # deprecation
    if (!is.null(polygon_feat_types)) {
        warning("polygon_feat_types is deprecated.\n Use z_indices instead")
        z_indices <- polygon_feat_types + 1L
    }

    segm_to_use <- paste0("p_", (segm_to_use - 1L))

    # data.table vars
    x <- y <- z <- NULL

    # provide your own custom names
    if (!is.null(custom_polygon_names)) {
        if (!is.character(custom_polygon_names)) {
            stop(wrap_txt("If custom_polygon_names are provided, it needs to
                        be a character vector"))
        }

        if (length(custom_polygon_names) != length(z_indices)) {
            stop(wrap_txt(
                "length of custom names need to be same as z_indices"
            ))
        }
    }

    if (isTRUE(verbose)) wrap_msg("Reading from:", boundaries_path)
    # list all files in the folder
    hdf5_boundary_list <- list.files(full.names = TRUE, boundaries_path)
    # only load subset of files if fov is given
    if (!is.null(fovs)) {
        selected_hdf5s <- paste0("feature_data_", fovs, ".hdf5")
        selected_hdf5s_concatenated <- paste0(selected_hdf5s, collapse = "|")
        hdf5_boundary_selected_list <- grep(
            selected_hdf5s_concatenated,
            x = hdf5_boundary_list, value = TRUE
        )
    } else {
        hdf5_boundary_selected_list <- hdf5_boundary_list
    }

    if (isTRUE(verbose)) {
        message("finished listing .hdf5 files start extracting .hdf5
                information")
    }

    # open selected polygon files

    # append data from all FOVs to single list
    init <- Sys.time()
    with_pbar({
        pb <- pbar(length(hdf5_boundary_selected_list) / 5)
        read_list <- lapply_flex(seq_along(hdf5_boundary_selected_list),
            future.packages = c("rhdf5", "Rhdf5lib"),
            function(init, z_indices, segm_to_use, bound_i) {
                read_file <- .h5_read_vizgen(
                    h5File = hdf5_boundary_selected_list[[bound_i]][[1]],
                    z_indices = z_indices,
                    segm_to_use = segm_to_use,
                    H5Fopen_flags = H5Fopen_flags
                )

                # update progress
                if (verbose) {
                    print(basename(hdf5_boundary_selected_list[[bound_i]]))
                }
                if (bound_i %% 5 == 0) {
                    pb()
                }

                return(read_file)
            },
            cores = cores,
            init = init,
            z_indices = z_indices,
            segm_to_use = segm_to_use
        )
    })

    # combine to FOV data single list
    read_DT <- data.table::rbindlist(read_list)

    # perform any necessary flips
    if (flip_x_axis) read_DT[, x := -1 * x]
    if (flip_y_axis) read_DT[, y := -1 * y]

    # separate polygons by z index
    zvals <- read_DT[, unique(z)]
    z_names <- paste0("z", zvals)
    z_read_DT <- lapply(seq_along(zvals), function(z_idx) {
        read_DT[z == zvals[z_idx], ]
    })
    names(z_read_DT) <- z_names
    if (!is.null(custom_polygon_names)) {
        poly_names <- custom_polygon_names
    } else {
        poly_names <- z_names
    }

    if (isTRUE(verbose)) message("finished extracting .hdf5 files")

    # outputs
    switch(output,
        "giottoPolygon" = .create_giotto_polygons_vizgen(
            z_read_DT = z_read_DT,
            poly_names = poly_names,
            set_neg_to_zero = set_neg_to_zero,
            calc_centroids = calc_centroids,
            smooth_polygons = smooth_polygons,
            smooth_vertices = smooth_vertices,
            create_gpoly_parallel = create_gpoly_parallel,
            create_gpoly_bin = create_gpoly_bin,
            verbose = verbose
        ),
        "data.table" = z_read_DT
    )
}




#' @keywords internal
#' @noRd
.create_giotto_polygons_vizgen <- function(
        z_read_DT,
        poly_names = names(z_read_DT),
        set_neg_to_zero = FALSE,
        calc_centroids = FALSE,
        smooth_polygons = TRUE,
        smooth_vertices = 60L,
        create_gpoly_parallel = TRUE,
        create_gpoly_bin = FALSE,
        verbose = TRUE) {
    checkmate::assert_list(z_read_DT)
    checkmate::assert_numeric(smooth_vertices)

    # data.table vars
    x <- y <- cell_id <- poly_ID <- NULL

    if (isTRUE(verbose)) message("start creating polygons")

    # **** sequential method ****
    if (!isTRUE(create_gpoly_parallel)) {
        with_pbar({
            pb <- pbar(along = z_read_DT)
            smooth_cell_polygons_list <- lapply(
                seq_along(z_read_DT), function(i) {
                    dfr_subset <- z_read_DT[[i]][, .(x, y, cell_id)]
                    data.table::setnames(
                        dfr_subset,
                        old = "cell_id", new = "poly_ID"
                    )
                    cell_polygons <- createGiottoPolygonsFromDfr(
                        segmdfr = dfr_subset,
                        name = poly_names[i],
                        calc_centroids = FALSE,
                        skip_eval_dfr = TRUE,
                        copy_dt = FALSE,
                        verbose = verbose
                    )
                    if (isTRUE(smooth_polygons)) {
                        cell_polygons <- smoothGiottoPolygons(
                            gpolygon = cell_polygons,
                            vertices = smooth_vertices,
                            k = 3L,
                            set_neg_to_zero = set_neg_to_zero
                        )
                    }
                    if (isTRUE(calc_centroids)) {
                        # NOTE: won't recalculate if centroids are already
                        # attached
                        cell_polygons <- centroids(
                            cell_polygons,
                            append_gpolygon = TRUE
                        )
                    }
                    pb(message = c(
                        poly_names[i], " (", i, "/", length(z_read_DT), ")"
                    ))
                    return(cell_polygons)
                }
            )
        })
        return(smooth_cell_polygons_list)
    }


    # **** parallel methods ****
    # no binning
    if (!is.numeric(create_gpoly_bin)) {
        with_pbar({
            pb <- pbar(along = z_read_DT)
            smooth_cell_polygons_list <- lapply_flex(
                seq_along(z_read_DT),
                future.packages = c("terra", "stats", "data.table"),
                function(i) {
                    dfr_subset <- z_read_DT[[i]][, .(x, y, cell_id)]
                    data.table::setnames(
                        dfr_subset,
                        old = "cell_id", new = "poly_ID"
                    )
                    cell_polygons <- gpoly_from_dfr_smoothed_wrapped(
                        segmdfr = dfr_subset,
                        name = poly_names[i],
                        skip_eval_dfr = TRUE,
                        copy_dt = FALSE,
                        smooth_polygons = smooth_polygons,
                        vertices = smooth_vertices,
                        set_neg_to_zero = set_neg_to_zero,
                        calc_centroids = calc_centroids,
                        verbose = verbose
                    )

                    pb(message = c(
                        poly_names[i], " (", i, "/", length(z_read_DT), ")"
                    ))
                    return(cell_polygons)
                }
            )
        })

        # unwrap results
        smooth_cell_polygons_list <- lapply(
            smooth_cell_polygons_list, function(x) {
                slot(x, "spatVector") <- terra::vect(slot(x, "spatVector"))
                if (isTRUE(calc_centroids)) {
                    slot(x, "spatVectorCentroids") <- terra::vect(
                        slot(x, "spatVectorCentroids")
                    )
                }
                return(x)
            }
        )
    } else {
        # with binning

        dfr_subset <- lapply(z_read_DT, function(bin, DT) {
            DT <- DT[, .(x, y, cell_id)]
            data.table::setnames(DT, old = "cell_id", new = "poly_ID")
            pid <- DT[, unique(poly_ID)]

            bin_pid <- data.table::data.table(
                "poly_ID" = pid,
                "bin_ID" = as.numeric(
                    cut(
                        x = seq_along(pid),
                        breaks = ceiling(length(pid) / bin)
                    )
                )
            )
            DT <- data.table::merge.data.table(
                DT, bin_pid,
                by = "poly_ID", all.x = TRUE
            )
            DT <- split(DT, DT$bin_ID)
        }, bin = create_gpoly_bin)

        bin_steps <- sum(unlist(lapply(dfr_subset, length)))

        with_pbar({
            pb <- pbar(steps = bin_steps)
            smooth_cell_polygons_list <- lapply( # sequential across z index
                seq_along(dfr_subset),
                function(i) {
                    lapply_flex( # parallelize across bins
                        dfr_subset[[i]],
                        future.packages = c("terra", "stats", "data.table"),
                        function(bin_DT) {
                            cell_polygons <- gpoly_from_dfr_smoothed_wrapped(
                                segmdfr = bin_DT,
                                name = poly_names[i],
                                skip_eval_dfr = TRUE,
                                copy_dt = FALSE,
                                smooth_polygons = smooth_polygons,
                                vertices = smooth_vertices,
                                set_neg_to_zero = set_neg_to_zero,
                                calc_centroids = calc_centroids,
                                verbose = verbose
                            )

                            pb(message = c(
                                poly_names[i], " (", i, "/",
                                length(dfr_subset), ")"
                            ))
                            return(cell_polygons)
                        }
                    )
                }
            )
        })

        # unwrap results
        smooth_cell_polygons_list <- lapply(
            seq_along(smooth_cell_polygons_list), function(i) {
                p_list <- lapply(smooth_cell_polygons_list[[i]], function(x) {
                    slot(x, "spatVector") <- terra::vect(slot(x, "spatVector"))
                    if (isTRUE(calc_centroids)) {
                        slot(x, "spatVectorCentroids") <- terra::vect(
                            slot(x, "spatVectorCentroids")
                        )
                    }
                    return(x)
                })
                # rbind results
                names(p_list) <- NULL
                return(do.call("rbind", p_list))
            }
        )
    }


    return(smooth_cell_polygons_list)
}







#' @title Read MERSCOPE polygons from parquet
#' @name readPolygonVizgenParquet
#' @description
#' Read Vizgen exported cell boundary parquet files as giottoPolyons. The z
#' level can be selected.
#' @param file parquet file to load
#' @param z_index either 'all' or a numeric vector of z_indices to get polygons
#' for
#' @param calc_centroids calculate centroids for the polygons (default = TRUE)
#' @param verbose be verbose
#' @returns giottoPolygons
#' @export
readPolygonVizgenParquet <- function(file,
    z_index = "all",
    calc_centroids = TRUE,
    verbose = TRUE) {
    # package checks
    package_check("arrow")
    package_check("sf")
    package_check("dplyr")


    checkmate::assert_file_exists(file)
    if (!setequal(z_index, "all")) {
        checkmate::assert_numeric(z_index)
    } else {
        checkmate::assert_true(identical(z_index, "all"))
    }

    # NSE vars
    ZIndex <- Geometry <- NULL

    # 1. determine z indices to get
    avail_z_idx <- arrow::open_dataset(file) %>%
        dplyr::distinct(ZIndex) %>%
        dplyr::pull() %>%
        # dplyr::pull(as_vector = TRUE) %>% # switch to this in future and add
        # arrow version requirement
        sort()

    get_z_idx <- if (setequal(z_index, "all")) {
        avail_z_idx
    } else if (is.numeric(z_index)) {
        z_index <- as.integer(z_index)
        if (!all(z_index %in% avail_z_idx)) {
            stop(paste("Not all z indices found in cell boundaries.\n
                    Existing indices are:", paste(avail_z_idx, collapse = " ")))
        }
        z_index
    }
    if (isTRUE(verbose)) {
        message("loading poly z_indices: ", paste(get_z_idx, collapse = " "))
    }


    # 2. collect by z index filter and convert WKB to multipolygon
    multipolygons <- lapply_flex(
        get_z_idx,
        function(z_idx) {
            # set schema
            schema <- arrow::open_dataset(file)$schema
            schema$EntityID <- arrow::string()

            # read and convert
            arrow::open_dataset(file, schema = schema) %>%
                dplyr::filter(ZIndex == z_idx) %>%
                dplyr::collect() %>%
                dplyr::mutate(Geometry = sf::st_as_sfc(Geometry))
        },
        future.seed = TRUE
    )
    names(multipolygons) <- lapply(
        multipolygons, function(x) paste0("z", unique(x$ZIndex))
    )


    # 3. convert to giottoPolygons and append meta
    out <- lapply(seq_along(multipolygons), function(i) {
        # append poly IDs and meta
        poly_table <- multipolygons[[i]]
        sv <- terra::vect(poly_table$Geometry)
        sv$poly_ID <- poly_table$EntityID
        sv$z_level <- poly_table$ZLevel

        gpoly <- giottoPolygon(
            name = names(multipolygons)[[i]],
            spatVector = sv,
            unique_ID_cache = poly_table$EntityID
        )

        if (isTRUE(calc_centroids)) {
            # NOTE: will not recalculate if centroids are already attached
            gpoly <- GiottoClass::centroids(x = gpoly, append_gpolygon = TRUE)
        }
    })

    return(out)
}









#' @title readPolygonFilesVizgen
#' @name readPolygonFilesVizgen
#' @description Read selected polygon files for the FOVs present in the Giotto
#' object and add the smoothed polygons to the object
#' @param gobject giotto object
#' @param boundaries_path path to the cell_boundaries folder
#' @param fovs selected fovs, if NULL select all fovs within Giotto object
#' @param polygon_feat_types a vector containing the polygon feature types
#' @param flip_x_axis flip x axis of polygon coordinates (multiply by -1)
#' @param flip_y_axis flip y axis of polygon coordinates (multiply by -1)
#' @param smooth_polygons smooth polygons (default = TRUE)
#' @param smooth_vertices number of vertices for smoothing
#' @param set_neg_to_zero set negative values to zero when smoothing
#' @param return_gobject return giotto object
#' @param verbose be verbose
#' @returns giotto object or cell polygons list
#' @seealso \code{\link{smoothGiottoPolygons}}
#' @export
readPolygonFilesVizgen <- function(
        gobject,
        boundaries_path,
        fovs = NULL,
        polygon_feat_types = 0:6,
        flip_x_axis = FALSE,
        flip_y_axis = FALSE,
        smooth_polygons = TRUE,
        smooth_vertices = 60,
        set_neg_to_zero = FALSE,
        return_gobject = TRUE,
        verbose = TRUE) {
    # define names
    poly_feat_names <- paste0("z", polygon_feat_types)
    poly_feat_indexes <- paste0("zIndex_", polygon_feat_types)

    # select FOVs present in the subset
    if (is.null(fovs)) {
        subset_metadata <- pDataDT(gobject)
        fovs <- unique(subset_metadata$fov)
    }



    smooth_cell_polygons_list <- readPolygonFilesVizgenHDF5(
        boundaries_path = boundaries_path,
        fovs = fovs,
        polygon_feat_types = polygon_feat_types,
        flip_x_axis = flip_x_axis,
        flip_y_axis = flip_y_axis,
        smooth_polygons = smooth_polygons,
        smooth_vertices = smooth_vertices,
        set_neg_to_zero = set_neg_to_zero,
        verbose = verbose
    )


    if (return_gobject) {
        # add cell polygons to Giotto object
        names(smooth_cell_polygons_list) <- poly_feat_names
        gobject <- addGiottoPolygons(
            gobject = gobject,
            gpolygons = smooth_cell_polygons_list
        )
        return(gobject)
    } else {
        return(smooth_cell_polygons_list)
    }
}




#' @describeIn readPolygonFilesVizgen (internal) Optimized .hdf5 reading for
#' vizgen merscope output. Returns a data.table of xyz coords and cell_id
#' @keywords internal
.h5_read_vizgen <- function(
        h5File,
        z_indices = 1L:7L,
        segm_to_use = "p_0",
        H5Fopen_flags = "H5F_ACC_RDWR") {
    # data.table vars
    group <- name <- cell <- z_name <- otype <- d_name <- cell_id <- NULL

    h5_ls <- data.table::setDT(
        rhdf5::h5ls(h5File, recursive = 5, datasetinfo = FALSE)
    )
    cell_names <- as.character(h5_ls[group == "/featuredata", name])
    z_names <- h5_ls[grep("zIndex", name), unique(name)]

    dset_names <- h5_ls[otype == "H5I_DATASET" & name == "coordinates", ]
    # subset by segm_to_use
    dset_names <- dset_names[grep(segm_to_use, group), ]
    # tag cellnames
    dset_names[, cell := gsub(
        pattern = "/featuredata/|/zIndex.*$", replacement = "", x = group
    )]
    # tag z_names
    dset_names[, z_name := gsub(
        pattern = "^.*/(zIndex_\\d*).*$", replacement = "\\1", x = group
    )]
    # subset by z_indices
    dset_names <- dset_names[z_name %in% z_names[z_indices], ]
    # create full file location
    dset_names[, d_name := paste0(group, "/", name)]

    fid <- rhdf5::H5Fopen(h5File, flags = H5Fopen_flags)
    dapl <- rhdf5::H5Pcreate("H5P_DATASET_ACCESS")

    contents <- lapply(cell_names, function(fid, dapl, cell_name) {
        zvals <- .h5_read_bare(
            file = fid,
            name = paste0(
                c("/featuredata", cell_name, "z_coordinates"),
                collapse = "/"
            ),
            dapl = dapl
        )
        names(zvals) <- z_names

        # subset to datasets related to cell
        cell_dsets <- dset_names[cell == cell_name, ]

        cell_data <- lapply(
            seq(nrow(cell_dsets)), function(fid, dapl, zvals, d_i) {
                res <- .h5_read_bare(
                    file = fid, name = cell_dsets[d_i, d_name], dapl = dapl
                )
                res <- t_flex(res[, , 1L])
                res <- cbind(res, zvals[cell_dsets[d_i, z_name]])
                colnames(res) <- c("x", "y", "z")
                res
            },
            fid = fid, dapl = dapl, zvals = zvals
        )
        cell_data <- data.table::as.data.table(do.call("rbind", cell_data))
        cell_data[, cell_id := cell_name]
        cell_data
    }, fid = fid, dapl = dapl)

    rhdf5::H5Pclose(dapl)
    rhdf5::H5Fclose(fid)
    contents <- data.table::rbindlist(contents)
    contents
}



#' @name .h5_read_bare
#' @title Read dataset from opened HDF5 with C functions
#' @param file opened HDF5 file
#' @param name dataset name within
#' @param dapl HDF5 property list (H5Pcreate('H5P_DATASET_ACCESS'))
#' @returns HDF5 contents
#' @keywords internal
.h5_read_bare <- function(file, name = "", dapl) {
    did <- .Call("_H5Dopen", file@ID, name, dapl@ID, PACKAGE = "rhdf5")
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, FALSE,
        PACKAGE = "rhdf5"
    )
    invisible(.Call("_H5Dclose", did, PACKAGE = "rhdf5"))

    res
}







## stereo-seq ####

#' @title getGEFtxCoords
#' @name getGEFtxCoords
#' @description Converts .gef file (output stereo-seq pipeline) into
#' transcript with coordinates
#' @param gef_file path to .gef file
#' @param bin_size bin size to select from .gef file
#' @returns transcript with coordinates
#' @export
getGEFtxCoords <- function(
        gef_file,
        bin_size = "bin100") {
    # data.table vars
    genes <- NULL

    # package check
    package_check(pkg_name = "rhdf5", repository = "Bioc")
    if (!file.exists(gef_file)) stop("File path to .gef file does not exist")

    # Step 1: Parse tx coords
    exprDT <- rhdf5::h5read(
        file = gef_file,
        name = paste0("geneExp/", bin_size, "/expression")
    )
    setDT(exprDT)

    # Step 2: Parse gene expression info using index
    geneDT <- rhdf5::h5read(
        file = gef_file,
        name = paste0("geneExp/", bin_size, "/gene")
    )
    setDT(geneDT)

    # Step 3: Combine read expression and gene data by repeating count
    # (match offset index)
    # See STOMICS file format manual for more information about exprDT and
    # geneDT
    exprDT[, genes := rep(x = geneDT$gene, geneDT$count)]

    return(exprDT)
}

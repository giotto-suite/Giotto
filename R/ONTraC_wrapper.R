#' @title installGiottoONTraCEnvironment
#' @description Installs a conda environment contains ONTraC. This
#' includes a miniconda installation and also a set of python packages that
#' Giotto may often use. See details for further information
#' @param packages_to_install python modules (packages) to install for Giotto.
#' @param python_version python version to use within the giotto conda
#' environment. Default is v3.11.9
#' @param ontrac_version ONTraC version to install. Default is "latest"
#' @param mini_install_path (optional) desired miniconda installation location.
#' Default is chosen by `reticulate::install_miniconda()`
#' @param confirm whether to pause for confirmation of conda environment
#' install location (default = TRUE)
#' @param envname name to assign environment. Default = "giotto_ontrac_env"
#' @param conda either "auto" (default) to allow reticulate to handle it, or
#' the full filepath to the conda executable. You can also set the option
#' "reticulate.conda_binary" or `Sys.setenv()` "RETICULATE_CONDA" to tell
#' reticulate where to look.
#' @param force_miniconda force reinstallation of miniconda
#' @param force_environment force reinstallation of the giotto environment
#' @param verbose be verbose
#' @returns installs a giotto environment using the reticulate miniconda system
#' @details This function will install a local conda environment using
#' the miniconda system as implemented by \pkg{reticulate}. Which contains
#' ONTraC and a set of python packages that Giotto may often use.
#'
#' @examples
#' installGiottoONTraCEnvironment()
#'
#' @export
installGiottoONTraCEnvironment <- function(python_version = "3.11.9",
    ontrac_version = "latest",
    mini_install_path = NULL,
    confirm = TRUE,
    envname = "giotto_ontrac_env",
    conda = "auto",
    force_miniconda = FALSE,
    force_environment = FALSE,
    verbose = NULL) {
    # handle ontrac version
    if (ontrac_version == "latest") {
        ontrac <- "ONTraC"
    } else {
        ontrac <- paste0("ONTraC==", ontrac_version)
    }

    # install conda env
    installGiottoEnvironment(
        packages_to_install = c(
            "pandas==2.2.1",
            "networkx==2.8.8",
            "python-igraph==0.10.2",
            "leidenalg==0.9.0",
            "python-louvain==0.16",
            "python.app==1.4",
            "scikit-learn==1.1.3",
            "smfishhmrf",
            "session-info",
            ontrac
        ),
        pip_packages = c(
            "python-louvain",
            "smfishhmrf",
            "session-info",
            "ONTraC"
        ),
        python_version = python_version,
        mini_install_path = mini_install_path,
        confirm = confirm,
        envname = envname,
        conda = "auto",
        force_miniconda = force_miniconda,
        force_environment = force_environment,
        verbose = verbose
    )
}


#' @title getONTraCv1Input
#' @name getONTraCv1Input
#' @description generate the input data for ONTraC v1
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param output_path the path to save the output file
#' @param cell_type the cell type column name in the metadata
#' @returns data.table with columns: Cell_ID, Sample, x, y, Cell_Type
#' @details This function generate the input data for ONTraC v1
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' getONTraCv1Input(
#'     gobject = g,
#'     cell_type = "custom_leiden"
#' )
#' @export
getONTraCv1Input <- function(
        gobject,
        cell_type,
        output_path = getwd(),
        spat_unit = NULL,
        feat_type = NULL,
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

    pos_df <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        output = "data.table"
    )
    meta_df <- pDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    output_df <- merge(x = pos_df, y = meta_df, by = "cell_ID")

    # check if the cell_type column exits
    if (!cell_type %in% colnames(output_df)) {
        vmsg(.v = verbose, paste(
            "Given",
            cell_type,
            "do not exist in giotto object's metadata!"
        ))
        return(NULL)
    }

    # add default sample name for one sample obj
    if (!"list_ID" %in% colnames(output_df)) {
        output_df$list_ID <- "ONTraC"
    }

    output_df <- output_df[, .SD, .SDcols = c(
        "cell_ID",
        "list_ID",
        "sdimx",
        "sdimy",
        cell_type
    )]
    colnames(output_df) <- c("Cell_ID", "Sample", "x", "y", "Cell_Type")
    file_path <- file.path(output_path, "ONTraC_dataset_input.csv")
    write.csv(output_df, file = file_path, quote = FALSE, row.names = FALSE)
    vmsg(.v = verbose, paste("ONTraC input file was saved as", file_path))

    return(output_df)
}


#' @title runONTraCV1
#' @name runONTraCV1
#' @description run ONTraC
#' @param dataset the path to the input data file
#' @param preprocessing_dir the directory to save the preprocessing results
#' @param GNN_dir the directory to save the GNN results
#' @param NTScore_dir the directory to save the NTScore results
#' @param n_cpu the number of CPUs used for niche network constructing. Default
#' is 4L
#' @param n_neighbors the number of neighbors used for ONTraC in niche network
#' construction. Default is 50L
#' @param n_local the index of local neighbor used for ONTraC in niche network
#' construction for normalization. Default is 20L
#' @param device the device used for ONTraC running GNN model. Default is "cpu"
#' @param epochs the maximum number of epochs for model training. Default is
#' 1000L
#' @param patience the number of epochs wait for better result. Default is 100L
#' @param min_delta the minimum change of loss to be considered as improvement.
#' Default is 0.001
#' @param min_epochs the minimum number of epochs to train. Default is 50L
#' @param batch_size the batch size for training. Default is 0L for whole
#' dataset
#' @param seed the random seed for reproducibility. Default is 42L
#' @param lr the learning rate for training. Default is 0.03
#' @param hidden_feats the number of hidden features for GNN model. Default is
#' 4L
#' @param k the number of neighbors for GNN model. Default is 6L
#' @param modularity_loss_weight the weight of modularity loss. Default is 0.3
#' @param purity_loss_weight the weight of purity loss. Default is 300.0
#' @param regularization_loss_weight the weight of regularization loss. Default
#' is 0.1
#' @param beta the weight of entropy loss. Default is 0.03
#' @param python_path, path to python executable within a conda/miniconda
#' environment. Default is "giotto_ontrac_env"
#' @returns none
#' @details This function runs ONTraC
#' @examples
#' runONTraCV1(
#'     dataset = "ONTraC_dataset_input.csv",
#'     preprocessing_dir = "preprocessing_dir",
#'     GNN_dir = "GNN_dir",
#'     NTScore_dir = "NTScore_dir",
#'     envname = "giotto_ontrac_env"
#' )
#' @export
runONTraCV1 <- function(ONTraC_input,
    dataset,
    preprocessing_dir,
    GNN_dir,
    NTScore_dir,
    n_cpu = 4L,
    n_neighbors = 50L,
    n_local = 20L,
    device = c("cpu", "cuda"),
    epochs = 1000L,
    patience = 100L,
    min_delta = 0.001,
    min_epochs = 50L,
    batch_size = 0L,
    seed = 42L,
    lr = 0.03,
    hidden_feats = 4L,
    k = 6L,
    modularity_loss_weight = 0.3,
    purity_loss_weight = 300.0,
    regularization_loss_weight = 0.1,
    beta = 0.03,
    python_path = "giotto_ontrac_env") {
    # parameters check
    device <- match.arg(device)

    # 1. identify operating system
    my_os <- get_os()

    # handle conda env
    set_giotto_python_path(python_path)
    python_path <- reticulate::conda_python(envname = python_path)
    if (my_os == "windows") {
        ONTraC_path <- file.path(dirname(python_path), "Scripts", "ONTraC")
    } else {
        ONTraC_path <- file.path(dirname(python_path), "ONTraC")
    }

    # run ONTraC
    command <- paste(
        ONTraC_path,
        "-d", dataset,
        "--preprocessing-dir", preprocessing_dir,
        "--GNN-dir", GNN_dir,
        "--NTScore-dir", NTScore_dir,
        "--n-cpu", n_cpu,
        "--n-neighbors", n_neighbors,
        "--n-local", n_local,
        "--device", device,
        "--epochs", epochs,
        "--patience", patience,
        "--min-delta", min_delta,
        "--min-epochs", min_epochs,
        "--batch-size", batch_size,
        "--seed", seed,
        "--lr", lr,
        "--hidden-feats", hidden_feats,
        "--k", k,
        "--modularity-loss-weight", modularity_loss_weight,
        "--purity-loss-weight", purity_loss_weight,
        "--regularization-loss-weight", regularization_loss_weight,
        "--beta", beta
    )
    wrap_msg(paste0("+", command))
    system(command)
}


#' @title load_cell_NT_score
#' @name load_cell_NT_score
#' @description load cell-level NT score
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved.
#' Default is getwd()
#' @param NTScore_dir the directory to save the NTScore results. Default is
#' file.path(ontrac_results_dir, "NTScore_dir")
#' @param NTScore_reverse whether to reverse the NTScore. Default is FALSE
#' @returns gobject with cell-level NT score
#' @details This function loads the ONTraC outputed cell-level NT score
load_cell_NT_score <- function(
        gobject,
        ontrac_results_dir = getwd(),
        NTScore_dir = file.path(
            ontrac_results_dir,
            "NTScore_dir"
        ),
        NTScore_reverse = FALSE) {
    NT_score_df <- read.csv(file = file.path(
        NTScore_dir, "NTScore.csv.gz"
    ))[c("Cell_ID", "Cell_NTScore")]
    colnames(NT_score_df) <- c("cell_ID", "NTScore")
    if (NTScore_reverse) {
        NT_score_df$NTScore <- 1 - NT_score_df$NTScore
    }
    gobject <- addCellMetadata(gobject,
        new_metadata = NT_score_df,
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )

    return(gobject)
}


#' @title load_cell_niche_cluster_prob
#' @name load_cell_niche_cluster_prob
#' @description load cell-niche cluster probability
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved.
#' Default is getwd()
#' @param GNN_dir the directory to save the GNN results. Default is
#' file.path(ontrac_results_dir, "GNN_dir")
#' @param name name for the probability matrix
#' @returns gobject with cell-niche cluster probability matrix
#' @details This function loads the ONTraC outputed cell-niche cluster
#' probability as an exprObj into the giotto object.
load_cell_niche_cluster_prob <- function(
        gobject,
        ontrac_results_dir = getwd(),
        GNN_dir = file.path(
            ontrac_results_dir,
            "GNN_dir"
        ),
        spat_unit = "cell",
        feat_type = "niche cluster",
        name = "prob") {
    niche_cluster_prob_df <- read.csv(file = file.path(
        GNN_dir, "cell_level_niche_cluster.csv.gz"
    ))
    rownames(niche_cluster_prob_df) <- niche_cluster_prob_df$Cell_ID
    niche_cluster_prob_df$Cell_ID <- NULL
    expobj <- createExprObj(t(niche_cluster_prob_df),
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    gobject <- GiottoClass::setExpression(
        gobject = gobject,
        x = expobj,
        spat_unit = spat_unit,
        feat_type = feat_type,
        name = name
    )

    return(gobject)
}


#' @title load_nc_connectivity
#' @name load_nc_connectivity
#' @description load niche cluster connectivity
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved.
#' Default is getwd()
#' @param GNN_dir the directory to save the GNN results. Default is
#' file.path(ontrac_results_dir, "GNN_dir")
#' @param name name for the connectivity matrix
#' @returns gobject with niche cluster connectivity matrix
#' @details This function loads the ONTraC outputed niche cluster connectivity
#' matrix as an exprObj into the giotto object.
load_nc_connectivity <- function(
        gobject,
        ontrac_results_dir = getwd(),
        GNN_dir = file.path(
            ontrac_results_dir,
            "GNN_dir"
        ),
        spat_unit = "niche cluster",
        feat_type = "connectivity",
        name = "normalized") {
    connectivity_df <- read.csv(file = file.path(
        GNN_dir, "consolidate_out_adj.csv.gz"
    ), header = FALSE)
    rownames(connectivity_df) <- paste0(
        "NicheCluster_",
        seq_len(dim(connectivity_df)[1]) - 1
    )
    colnames(connectivity_df) <- paste0(
        "NicheCluster_",
        seq_len(dim(connectivity_df)[2]) - 1
    )
    expobj <- createExprObj(t(connectivity_df),
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    gobject <- GiottoClass::setExpression(
        gobject = gobject,
        x = expobj,
        spat_unit = spat_unit,
        feat_type = feat_type,
        name = name
    )

    return(gobject)
}


#' @title load_niche_cluster_nt_score
#' @name load_niche_cluster_nt_score
#' @description load niche cluster NT score
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved.
#' Default is getwd()
#' @param NTScore_dir the directory to save the NTScore results. Default is
#' file.path(ontrac_results_dir, "NTScore_dir")
#' @param NTScore_reverse whether to reverse the NTScore. Default is FALSE
#' @returns gobject with niche cluster NT score
#' @details This function loads the ONTraC outputed niche cluster NT score
#' into the giotto object.
load_niche_cluster_nt_score <- function(
        gobject,
        ontrac_results_dir = getwd(),
        NTScore_dir = file.path(
            ontrac_results_dir,
            "NTScore_dir"
        ),
        NTScore_reverse = FALSE) {
    niche_cluster_df <- read.csv(file = file.path(
        NTScore_dir, "niche_cluster_score.csv.gz"
    ), header = FALSE)
    colnames(niche_cluster_df) <- c("NTScore")
    niche_cluster_df$feat_ID <- paste0(
        "NicheCluster_",
        seq_len(dim(niche_cluster_df)[1]) - 1
    )
    if (NTScore_reverse) {
        niche_cluster_df$NTScore <- 1 - niche_cluster_df$NTScore
    }
    gobject <- GiottoClass::addCellMetadata(
        gobject = gobject,
        spat_unit = "niche cluster",
        feat_type = "connectivity",
        new_metadata = niche_cluster_df,
        by_column = TRUE,
        column_cell_ID = "feat_ID"
    )
    niche_cluster_meta_obj <- GiottoClass::createFeatMetaObj(niche_cluster_df)
    gobject <- GiottoClass::setFeatureMetadata(
        gobject = gobject,
        x = niche_cluster_meta_obj,
        spat_unit = "cell",
        feat_type = "niche cluster"
    )

    return(gobject)
}


#' @title cal_cell_niche_cluster_bin
#' @name cal_cell_niche_cluster_bin
#' @description calculate binarized cell-level niche cluster assignment
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @returns gobject with binarized cell-level niche cluster assignment
cal_cell_niche_cluster_bin <- function(gobject,
    spat_unit = "cell",
    feat_type = "niche cluster") {
    # calculate the binarized cell-level niche cluster assignment
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "matrix"
    )
    niche_cluster_bin <- rownames(expr_values)[apply(
        expr_values,
        2,
        function(x) which.max(x)
    )]
    ori_meta_df <- pDataDT(gobject)
    new_meta_df <- data.frame(
        cell_ID = ori_meta_df$cell_ID,
        new_feat = niche_cluster_bin
    )
    colnames(new_meta_df) <- c("cell_ID", feat_type)
    # get NTScore for each niche cluster
    nc_meta_df <- fDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    # order the niche cluster by NT score
    sorted_nc <- nc_meta_df$feat_ID[order(nc_meta_df$NTScore)]
    new_meta_df[[feat_type]] <- factor(
        new_meta_df[[feat_type]],
        levels = sorted_nc
    )
    # add the new metadata to the giotto object
    gobject <- addCellMetadata(gobject,
        new_metadata = new_meta_df,
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )
    return(gobject)
}


#' @title loadOntraCResults
#' @name loadOntraCResults
#' @description load ONTraC results
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved
#' @param preprocessing_dir the directory to save the preprocessing results.
#' Default is file.path(ontrac_results_dir, "preprocessing_dir")
#' @param GNN_dir the directory to save the GNN results. Default is
#' file.path(ontrac_results_dir, "GNN_dir")
#' @param NTScore_dir the directory to save the NTScore results. Default is
#' file.path(ontrac_results_dir, "NTScore_dir")
#' @param NTScore_reverse whether to reverse the NTScore. Default is FALSE
#' @returns gobject with ONTraC results
#' @details This function loads the ONTraC results into the giotto object.
#' @export
loadOntraCResults <- function(
        gobject,
        ontrac_results_dir = getwd(),
        preprocessing_dir = file.path(
            ontrac_results_dir,
            "preprocessing_dir"
        ),
        GNN_dir = file.path(
            ontrac_results_dir,
            "GNN_dir"
        ),
        NTScore_dir = file.path(
            ontrac_results_dir,
            "NTScore_dir"
        ),
        NTScore_reverse = FALSE) {
    gobject <- load_cell_NT_score(
        gobject = gobject,
        ontrac_results_dir = ontrac_results_dir,
        NTScore_dir = NTScore_dir,
        NTScore_reverse = NTScore_reverse
    )
    gobject <- load_cell_niche_cluster_prob(
        gobject = gobject,
        ontrac_results_dir = ontrac_results_dir,
        GNN_dir = GNN_dir
    )
    gobject <- GiottoClass::addCellMetadata(
        gobject = gobject,
        spat_unit = "cell",
        feat_type = "niche cluster",
        new_metadata = pDataDT(gobject),
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )
    gobject <- load_nc_connectivity(
        gobject = gobject,
        ontrac_results_dir = ontrac_results_dir,
        GNN_dir = GNN_dir
    )
    gobject <- load_niche_cluster_nt_score(
        gobject = gobject,
        ontrac_results_dir = ontrac_results_dir,
        NTScore_dir = NTScore_dir,
        NTScore_reverse = NTScore_reverse
    )
    gobject <- cal_cell_niche_cluster_bin(
        gobject = gobject
    )

    return(gobject)
}



#' @title plotSpatNicheClusterProb
#' @name plotSpatNicheClusterProb
#' @description plot spatial niche cluster probability
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_unit name of spatial unit niche stored cluster features
#' @param feat_type name of the feature type stored probability matrix
#' @param expression_values name of the expression matrix stored probability
#' values
#' @param ... additional arguments to be passed to the spatFeatPlot2D function
#' @details This function plots the spatial niche cluster probability
#' @returns ggplot
#' @export
plotSpatNicheClusterProb <- function(gobject,
    spat_unit = "cell",
    feat_type = "niche cluster",
    expression_values = "prob",
    ...,
    default_save_name = "spatNicheClusterProb") {
    nc_meta_df <- fDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    sorted_nc <- nc_meta_df$feat_ID[order(nc_meta_df$NTScore)]

    spatFeatPlot2D(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = expression_values,
        feats = sorted_nc,
        ...
    )
}


#' @title plotSpatNicheClusterBin
#' @name plotSpatNicheClusterBin
#' @description plot spatial niche cluster binarized
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_unit name of spatial unit niche stored cluster features
#' @param feat_type name of the feature type stored binarized niche cluster
#' @param niche_cluster_label name of the niche cluster label
#' @param ... additional arguments to be passed to the spatFeatPlot2D function
#' @details This function plots the spatial niche cluster binarized
#' @returns ggplot
#' @export
plotSpatNicheClusterBin <- function(gobject,
    spat_unit = "cell",
    feat_type = "niche cluster",
    ...,
    default_save_name = "spatNicheClusterBin") {
    # determine the color code
    nc_meta_df <- fDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    # order the niche cluster by NT score
    sorted_nc <- nc_meta_df$feat_ID[order(nc_meta_df$NTScore)]
    nc_num <- length(sorted_nc)
    # cell color code is a named vector with the niche cluster label
    cell_color_code <- setNames(
        viridis::turbo(n = nc_num + 2)[1:nc_num + 1],
        sorted_nc
    )

    spatPlot2D(
        gobject = gobject,
        spat_unit = spat_unit,
        cell_color = feat_type,
        cell_color_code = cell_color_code,
        ...
    )
}


#' @title plotNicheClusterConnectivity
#' @name plotNicheClusterConnectivity
#' @description plot niche cluster connectivity
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_unit name of spatial unit niche stored cluster features
#' @param feat_type name of the feature type stored niche cluster connectivities
#' @param values name of the expression matrix stored connectivity values
#' @details This function plots the niche cluster connectivity matrix
#' @returns ggplot
#' @export
plotNicheClusterConnectivity <- function(gobject,
    spat_unit = "niche cluster",
    feat_type = "connectivity",
    values = "normalized",
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "NicheClusterConnectivity") {
    # load `guide_edge_colourbar` function in ggraph,
    # otherwise it will raise an error when using `scale_edge_colour_gradientn`
    library(ggraph)

    # get the niche cluster connectivity matrix
    niche_cluster_connectivites <- getExpression(
        gobject = gobject,
        values = values,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "matrix"
    )
    nc_num <- dim(niche_cluster_connectivites)[1]

    # transform the matrix to data.frame for constructing igraph object
    niche_cluster_connectivites <- cbind(
        expand.grid(dimnames(niche_cluster_connectivites)),
        value = as.vector(as.matrix(
            niche_cluster_connectivites
        ))
    )
    colnames(niche_cluster_connectivites) <- c("from", "to", "connectivites")

    # construct igraph object
    igd <- igraph::graph_from_data_frame(
        d = niche_cluster_connectivites[, c("from", "to", "connectivites")],
        directed = FALSE
    )
    igd <- igraph::simplify(
        graph = igd,
        remove.loops = TRUE,
        remove.multiple = FALSE
    )
    # set edge attributes
    edges_sizes <- igraph::edge_attr(igd, "connectivites")
    edges_colors <- edges_sizes
    igd <- igraph::set_edge_attr(
        graph = igd,
        index = igraph::E(igd),
        name = "color",
        value = edges_colors
    )
    igd <- igraph::set_edge_attr(
        graph = igd,
        index = igraph::E(igd),
        name = "size",
        value = edges_sizes
    )
    # set node attributes
    nc_meta_df <- pDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    igraph::V(igd)$NTScore <- nc_meta_df$NTScore[
        match(
            igraph::V(igd)$name,
            nc_meta_df$cell_ID
        )
    ]

    # plot
    ## layout
    coords <- igraph::layout_with_drl(
        graph = igd,
        weights = edges_sizes,
        use.seed = TRUE
    )
    gpl <- ggraph::ggraph(graph = igd, layout = coords)

    ## edges
    gpl <- gpl + ggraph::geom_edge_link(
        ggplot2::aes(
            colour = edges_sizes,
            edge_width = 5,
            edge_alpha = size # nolint: object_usage_linter.
        ),
        show.legend = FALSE
    )
    gpl <- gpl + ggraph::scale_edge_alpha(range = c(0.1, 1))
    gpl <- gpl + ggraph::scale_edge_colour_gradientn(
        colours = getColors("Reds", 9, src = "RColorBrewer"),
        name = "Value"
    )

    ## node
    gpl <- gpl + ggraph::geom_node_point(
        ggplot2::aes(colour = NTScore), # nolint: object_usage_linter.
        size = 10
    )
    gpl <- gpl + ggplot2::scale_colour_gradientn(
        colours = viridis::turbo(n = nc_num + 2)[1:nc_num + 1]
    )
    gpl <- gpl + ggraph::geom_node_text(
        ggplot2::aes(label = name), # nolint: object_usage_linter.
        repel = TRUE
    )

    ## theme
    gpl <- gpl + ggplot2::theme_bw() + ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
    )
    gpl

    # return or save
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = gpl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}

#' @title plotCTCompositionInNicheCluster
#' @name plotCTCompositionInNicheCluster
#' @description plot cell type composition within each niche cluster
#' @param cell_type the cell type column name in the metadata
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_unit name of spatial unit niche stored cluster features
#' @param feat_type name of the feature type stored probability matrix
#' @param normalization normalization method for the cell type composition
#' @param values name of the expression matrix stored probability of each cell
#' assigned to each niche cluster
#' @details This function plots the cell type composition within each niche
#' cluster
#' @returns ggplot
#' @export
plotCTCompositionInNicheCluster <- function(gobject,
    cell_type,
    values = "prob",
    spat_unit = "cell",
    feat_type = "niche cluster",
    normalization = c("by_niche_cluster", "by_cell_type", NULL),
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "CellTypeCompositionInNicheCluster") {
    normalization <- match.arg(normalization)

    # Get the cell type composition within each niche cluster
    ## extract the cell-level niche cluster probability matrix
    exp <- getExpression(
        gobject = gobject,
        values = values,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "matrix"
    )
    prob_df <- as.data.frame(t(as.matrix(exp)))
    prob_df$cell_ID <- rownames(prob_df)
    ## combine the cell type and niche cluster probability matrix
    combined_df <- merge(
        as.data.frame(pDataDT(gobject))[, c(
            "cell_ID",
            cell_type
        )],
        prob_df,
        by = "cell_ID"
    )

    # Calculate the normalized cell type composition within each niche cluster
    cell_type_counts_df <- combined_df %>%
        tidyr::pivot_longer(
            cols = dplyr::starts_with("NicheCluster_"),
            names_to = "Cluster",
            values_to = "Probability"
        ) %>%
        dplyr::group_by(
            !!rlang::sym(cell_type),
            Cluster # nolint: object_usage_linter.
        ) %>%
        dplyr::summarise(Sum = sum(Probability, # nolint: object_usage_linter.
            na.rm = TRUE
        )) %>%
        tidyr::spread(key = "Cluster", value = "Sum", fill = 0)
    cell_type_counts_df <- as.data.frame(cell_type_counts_df)
    rownames(cell_type_counts_df) <- cell_type_counts_df[[cell_type]]
    cell_type_counts_df[[cell_type]] <- NULL

    if (normalization == "by_cell_type") {
        normalized_df <- as.data.frame(
            cell_type_counts_df / rowSums(cell_type_counts_df)
        )
    } else if (normalization == "by_niche_cluster") {
        normalized_df <- as.data.frame(t(
            t(cell_type_counts_df) / colSums(cell_type_counts_df)
        ))
    } else if (normalization == NULL) {
        normalized_df <- cell_type_counts_df
    }

    # Reshape the data frame into long format
    normalized_df[[cell_type]] <- rownames(normalized_df)
    df_long <- normalized_df %>%
        tidyr::pivot_longer(
            cols = -!!rlang::sym(cell_type), # nolint: object_usage_linter.
            names_to = "Cluster",
            values_to = "Composition"
        )

    # Order the niche clusters by NTscore
    nc_meta_df <- fDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    df_long$Cluster <- factor(df_long$Cluster,
        levels = nc_meta_df$feat_ID[order(nc_meta_df$NTScore)]
    )

    # Order the cell types by the average NTScore
    data_df <- pDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    avg_scores <- data_df %>%
        dplyr::group_by(!!rlang::sym(cell_type)) %>%
        dplyr::summarise(Avg_NTScore = mean(NTScore))
    # nolint: object_usage_linter.
    df_long[[cell_type]] <- factor(df_long[[cell_type]],
        levels = avg_scores[[cell_type]][order(avg_scores$Avg_NTScore)]
    )

    # Create the heatmap using ggplot2
    pl <- ggplot(df_long, aes(
        x = !!rlang::sym(cell_type), # nolint: object_usage_linter.
        y = Cluster, # nolint: object_usage_linter.
        fill = Composition # nolint: object_usage_linter.
    )) +
        geom_tile() +
        viridis::scale_fill_viridis(option = "inferno", limits = c(0, 1)) +
        theme_minimal() +
        labs(
            title = "Normalized cell type compositions within each niche
            cluster",
            x = "Cell_Type",
            y = "Cluster"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # return or save
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


#' @title plotCellTypeNTScore
#' @name plotCellTypeNTScore
#' @description plot NTScore by cell type
#' @param cell_type the cell type column name in the metadata
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @returns ggplot
#' @export
plotCellTypeNTScore <- function(
        gobject,
        cell_type,
        values = "NTScore",
        spat_unit = "cell",
        feat_type = "niche cluster",
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "CellTypeNTScore") {
    # Get the cell type composition within each niche cluster
    data_df <- pDataDT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    avg_scores <- data_df %>%
        dplyr::group_by(!!rlang::sym(cell_type)) %>%
        dplyr::summarise(Avg_NTScore = mean(!!rlang::sym(values)))
    data_df[[cell_type]] <- factor(data_df[[cell_type]],
        levels = avg_scores[[cell_type]][order(avg_scores$Avg_NTScore)]
    )

    pl <- ggplot(data_df, aes(
        x = !!rlang::sym(values),
        y = !!rlang::sym(cell_type),
        fill = !!rlang::sym(cell_type)
    )) +
        geom_violin() +
        theme_minimal() +
        labs(
            title = "Violin Plot of NTScore by Cell Type",
            x = values,
            y = "Cell Type"
        ) +
        ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # return or save
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

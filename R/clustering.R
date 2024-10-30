#' @title doLeidenCluster
#' @name doLeidenCluster
#' @description cluster cells using a NN-network and the Leiden community
#' detection algorithm
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name for cluster, default to "leiden_clus"
#' @param nn_network_to_use type of NN network to use (kNN vs sNN), default to
#' "sNN"
#' @param network_name name of NN network to use, default to "sNN.pca"
#' @param python_path specify specific path to python if required
#' @param resolution resolution, default = 1
#' @param weight_col weight column to use for edges, default to "weight"
#' @param partition_type The type of partition to use for optimization.
#' (e.g. "RBConfigurationVertexPartition", "ModularityVertexPartition")
#' @param init_membership initial membership of cells for the partition
#' @param n_iterations number of interactions to run the Leiden algorithm.
#' If the number of iterations is negative, the Leiden algorithm is run until
#' an iteration in which there was no improvement.
#' @param return_gobject logical. return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @returns giotto object with new clusters appended to cell metadata
#' @details
#' This function is a wrapper for the Leiden algorithm implemented in python,
#' which can detect communities in graphs of millions of nodes (cells),
#' as long as they can fit in memory. See the
#' [leidenalg](https://github.com/vtraag/leidenalg)
#' github page or the
#' [readthedocs](https://leidenalg.readthedocs.io/en/stable/index.html)
#' page for more information.
#'
#' Partition types available and information:
#' * **RBConfigurationVertexPartition:** Implements Reichardt and Bornholdt’s
#' Potts model with a configuration null model. This quality function is
#' well-defined only for positive edge weights. This quality function uses a
#' linear resolution parameter.
#' * **ModularityVertexPartition:** Implements modularity. This quality
#' function is well-defined only for positive edge weights. It does \emph{not}
#' use the resolution parameter
#'
#' Set \emph{weight_col = NULL} to give equal weight (=1) to each edge.
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' doLeidenCluster(g)
#' @export
doLeidenCluster <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    name = "leiden_clus",
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    python_path = NULL,
    resolution = 1,
    weight_col = "weight",
    partition_type = c(
        "RBConfigurationVertexPartition",
        "ModularityVertexPartition"
    ),
    init_membership = NULL,
    n_iterations = 1000,
    return_gobject = TRUE,
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

    ## get cell IDs ##
    cell_ID_vec <- gobject@cell_ID[[spat_unit]]

    ## select network to use
    igraph_object <- getNearestNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        nn_type = nn_network_to_use,
        name = network_name,
        output = "igraph"
    )



    ## select partition type
    partition_type <- match.arg(partition_type,
        choices = c(
            "RBConfigurationVertexPartition", "ModularityVertexPartition"
        )
    )

    ## check or make paths
    # python path
    if (is.null(python_path)) {
        python_path <- readGiottoInstructions(gobject, param = "python_path")
    }

    ## prepare python path and louvain script
    reticulate::use_python(required = TRUE, python = python_path)
    python_leiden_function <- system.file("python", "python_leiden.py",
        package = "Giotto"
    )
    reticulate::source_python(file = python_leiden_function)

    ## set seed
    if (isTRUE(set_seed)) {
        seed_number <- as.integer(seed_number)
    } else {
        seed_number <- as.integer(sample(x = seq(10000), size = 1))
    }

    ## extract NN network
    network_edge_dt <- data.table::as.data.table(
        igraph::as_data_frame(x = igraph_object, what = "edges")
    )

    # data.table variables
    weight <- NULL

    # add weight for edges or set to 1 for all
    if (!is.null(weight_col)) {
        if (!weight_col %in% colnames(network_edge_dt)) {
            stop("weight column is not an igraph attribute")
        } else {
            # weight is defined by attribute of igraph object
            network_edge_dt <- network_edge_dt[
                , c("from", "to", weight_col),
                with = FALSE
            ]
            data.table::setnames(network_edge_dt, weight_col, "weight")
        }
    } else {
        # weight is the same
        network_edge_dt <- network_edge_dt[, c("from", "to"), with = FALSE]
        network_edge_dt[, weight := 1]
    }




    ## do python leiden clustering
    reticulate::py_set_seed(
        seed = seed_number,
        disable_hash_randomization = TRUE
    )
    pyth_leid_result <- python_leiden(
        df = network_edge_dt,
        partition_type = partition_type,
        initial_membership = init_membership,
        weights = "weight",
        n_iterations = n_iterations,
        seed = seed_number,
        resolution_parameter = resolution
    )

    ident_clusters_DT <- data.table::data.table(
        cell_ID = pyth_leid_result[[1]], "name" = pyth_leid_result[[2]]
    )
    data.table::setnames(ident_clusters_DT, "name", name)



    ## add clusters to metadata ##
    if (return_gobject == TRUE) {
        cluster_names <- names(pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        ))

        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- getCellMetadata(
                gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "cellMetaObj",
                copy_obj = TRUE
            )

            cell_metadata[][, eval(name) := NULL]

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setCellMetadata(
                gobject,
                x = cell_metadata,
                verbose = FALSE,
                initialize = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }

        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = ident_clusters_DT[
                , c("cell_ID", name), with = FALSE],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )

        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_cluster")
        return(gobject)
    } else {
        # else return clustering result
        return(ident_clusters_DT)
    }
}







#' @title doLeidenClusterIgraph
#' @name doLeidenClusterIgraph
#' @description cluster cells using a NN-network and the Leiden community
#' detection algorithm as implemented in igraph
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name for cluster, default to "leiden_clus"
#' @param nn_network_to_use type of NN network to use (kNN vs sNN), default to
#' "sNN"
#' @param network_name name of NN network to use, default to "sNN.pca"
#' @param objective_function objective function for the leiden algo
#' @param weights weights of edges
#' @param resolution resolution, default = 1
#' @param resolution_parameter deprecated. Use `resolution` instead
#' @param beta leiden randomness
#' @param initial_membership initial membership of cells for the partition
#' @param n_iterations number of interations to run the Leiden algorithm.
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @inheritDotParams igraph::cluster_leiden -graph -objective_function
#' -resolution_parameter -beta -weights -initial_membership -n_iterations
#' -resolution
#' @returns giotto object with new clusters appended to cell metadata
#' @details
#' This function is a wrapper for the Leiden algorithm implemented in igraph,
#' which can detect communities in graphs of millions of nodes (cells),
#' as long as they can fit in memory. See \code{\link[igraph]{cluster_leiden}}
#' for more information.
#'
#' Set \emph{weights = NULL} to use the vertices weights associated with the
#' igraph network.
#' Set \emph{weights = NA} if you don't want to use vertices weights
#'
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' doLeidenClusterIgraph(g)
#' @export
doLeidenClusterIgraph <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    name = "leiden_clus",
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    objective_function = c("modularity", "CPM"),
    weights = NULL,
    resolution = 1,
    resolution_parameter = deprecated(),
    beta = 0.01,
    initial_membership = NULL,
    n_iterations = 1000,
    return_gobject = TRUE,
    set_seed = TRUE,
    seed_number = 1234,
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
    
    resolution <- deprecate_param(
        x = resolution_parameter, 
        y = resolution,
        fun = "doLeidenClusterIgraph",
        when = "4.1.4"
    )

    ## get cell IDs ##
    cell_ID_vec <- gobject@cell_ID[[spat_unit]]

    ## select network to use
    igraph_object <- getNearestNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        nn_type = nn_network_to_use,
        name = network_name,
        output = "igraph"
    )

    ## select partition type
    objective_function <- match.arg(objective_function,
        choices = c("modularity", "CPM")
    )

    ## set seed
    if (isTRUE(set_seed)) {
        seed_number <- as.integer(seed_number)
        set.seed(seed_number)
        on.exit(expr = {
            GiottoUtils::random_seed(set.seed = TRUE)
        }, add = TRUE)
    }

    # make igraph network undirected
    graph_object_undirected <- igraph::as.undirected(igraph_object)
    leiden_clusters <- igraph::cluster_leiden(
        graph = graph_object_undirected,
        objective_function = objective_function,
        resolution = resolution,
        beta = beta,
        weights = weights,
        initial_membership = initial_membership,
        n_iterations = n_iterations,
        ...
    )

    # summarize results
    ident_clusters_DT <- data.table::data.table(
        "cell_ID" = leiden_clusters$names, "name" = leiden_clusters$membership
    )
    data.table::setnames(ident_clusters_DT, "name", name)



    ## add clusters to metadata ##
    if (isTRUE(return_gobject)) {
        cluster_names <- names(pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        ))

        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- getCellMetadata(gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "cellMetaObj",
                copy_obj = TRUE
            )

            cell_metadata[][, eval(name) := NULL]

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setCellMetadata(gobject,
                x = cell_metadata,
                verbose = FALSE,
                initialize = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }

        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = ident_clusters_DT[
                , c("cell_ID", name), with = FALSE],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )

        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_cluster")
        return(gobject)
    } else {
        # else return clustering result
        return(ident_clusters_DT)
    }
}







#' @title doGiottoClustree
#' @name doGiottoClustree
#' @description cluster cells using leiden methodology to visualize different
#' resolutions
#' @param gobject giotto object
#' @param res_vector vector of different resolutions to test
#' @param res_seq list of float numbers indicating start, end, and step size
#' for resolution testing, i.e. (0.1, 0.6, 0.1)
#' @param return_gobject default FALSE. See details for more info.
#' @param show_plot by default, pulls from provided gobject instructions
#' @param save_plot by default, pulls from provided gobject instructions
#' @param return_plot by default, pulls from provided gobject instructions
#' @param save_param list of saving parameters from
#' [GiottoVisuals::all_plots_save_function()]
#' @param default_save_name name of saved plot, default "clustree"
#' @param verbose be verbose
#' @inheritDotParams clustree::clustree -x
#' @returns a plot object (default), OR a giotto object (if specified)
#' @details This function tests different resolutions for Leiden clustering and
#' provides a visualization of cluster sizing as resolution varies.
#'
#' By default, the tested leiden clusters are NOT saved to the Giotto object,
#' and a plot is returned.
#'
#' If return_gobject is set to TRUE, and a giotto object with *all* tested
#' leiden cluster information
#' will be returned.
#' @seealso \code{\link{doLeidenCluster}}
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' doGiottoClustree(
#'     gobject = g, res_vector = c(0.5, 0.8), return_plot = FALSE,
#'     show_plot = FALSE, save_plot = FALSE
#' )
#' @export
doGiottoClustree <- function(gobject,
    res_vector = NULL,
    res_seq = NULL,
    return_gobject = FALSE,
    show_plot = NULL,
    save_plot = NULL,
    return_plot = NULL,
    save_param = list(),
    default_save_name = "clustree",
    verbose = TRUE,
    ...) {
    package_check(pkg_name = "clustree", repository = "CRAN")
    ## setting resolutions to use
    if (is.null(res_vector)) {
        if (!is.null(res_seq)) {
            res_vector <- seq(res_seq[1], res_seq[2], res_seq[3])
        } else {
            stop("Please input res_vector or res_seq parameters")
        }
    }

    ## performing multiple leiden clusters at resolutions specified
    for (i in res_vector) {
        if (isTRUE(verbose)) wrap_msg("Calculating leiden res:", i)
        gobject <- doLeidenCluster(
            gobject = gobject,
            resolution = i,
            name = paste0("leiden_clustree_", i)
        )
    }

    ## plotting clustree graph
    pl <- clustree::clustree(
        x = pDataDT(gobject),
        prefix = "leiden_clustree_",
        ...
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
        else_return = NULL
    ))
}


#' @title doLouvainCluster community
#' @name .doLouvainCluster_community
#' @description cluster cells using a NN-network and the Louvain algorithm
#' from the community module in Python
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name for cluster, default to "louvain_clus"
#' @param nn_network_to_use type of NN network to use (kNN vs sNN), default to
#' "sNN"
#' @param network_name name of NN network to use, default to "sNN.pca"
#' @param python_path specify specific path to python if required
#' @param resolution resolution, default = 1
#' @param weight_col weight column to use for edges
#' @param louv_random Will randomize the node evaluation order and the
#' community evaluation order to get different partitions at each call
#' (default = FALSE)
#' @param return_gobject logical. return giotto object (default = TRUE)
#' @param set_seed set seed (default = FALSE)
#' @param seed_number number for seed
#' @param \dots additional params to pass
#' @returns giotto object with new clusters appended to cell metadata
#' @details This function is a wrapper for the Louvain algorithm implemented in
#' Python, which can detect communities in graphs of nodes (cells).
#' See the
#' [readthedocs](https://python-louvain.readthedocs.io/en/latest/index.html)
#' page for more information.
#'
#' Set \emph{weight_col = NULL} to give equal weight (=1) to each edge.
#' @md
#' @keywords internal
.doLouvainCluster_community <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    name = "louvain_clus",
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    python_path = NULL,
    resolution = 1,
    weight_col = NULL,
    louv_random = FALSE,
    return_gobject = TRUE,
    set_seed = FALSE,
    seed_number = 1234,
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


    ## get cell IDs ##
    cell_ID_vec <- gobject@cell_ID[[spat_unit]]

    ## select network to use
    igraph_object <- getNearestNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        nn_type = nn_network_to_use,
        name = network_name,
        output = "igraph"
    )

    ## check or make paths
    # python path
    if (is.null(python_path)) {
        python_path <- readGiottoInstructions(gobject, param = "python_path")
    }

    # prepare python path and louvain script
    reticulate::use_python(required = TRUE, python = python_path)
    python_louvain_function <- system.file(
        "python", "python_louvain.py",
        package = "Giotto"
    )
    reticulate::source_python(file = python_louvain_function)

    # set seed
    if (isTRUE(set_seed)) {
        seed_number <- as.integer(seed_number)
    } else {
        seed_number <- as.integer(sample(x = seq(10000), size = 1))
    }

    network_edge_dt <- data.table::as.data.table(igraph::as_data_frame(
        x = igraph_object, what = "edges"
    ))

    # data.table variables
    weight <- NULL

    if (!is.null(weight_col)) {
        if (!weight_col %in% colnames(network_edge_dt)) {
            stop("weight column is not an igraph attribute")
        } else {
            # weight is defined by attribute of igraph object
            network_edge_dt <- network_edge_dt[
                , c("from", "to", weight_col),
                with = FALSE
            ]
            setnames(network_edge_dt, weight_col, "weight")
        }
    } else {
        # weight is the same
        network_edge_dt <- network_edge_dt[, c("from", "to"), with = FALSE]
        network_edge_dt[, weight := 1]
    }

    # do python louvain clustering
    if (louv_random == FALSE) {
        reticulate::py_set_seed(
            seed = seed_number, disable_hash_randomization = TRUE
        )
        pyth_louv_result <- python_louvain(
            df = network_edge_dt, resolution = resolution, randomize = FALSE
        )
    } else {
        reticulate::py_set_seed(
            seed = seed_number, disable_hash_randomization = TRUE
        )
        pyth_louv_result <- python_louvain(
            df = network_edge_dt,
            resolution = resolution,
            random_state = seed_number
        )
    }
    ident_clusters_DT <- data.table::data.table(
        cell_ID = rownames(pyth_louv_result), "name" = pyth_louv_result[[1]]
    )
    data.table::setnames(ident_clusters_DT, "name", name)


    ## return
    if (isTRUE(return_gobject)) {
        # get cell metadata names
        cluster_names <- names(pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        ))


        # set name/cluster results to NULL if already exist
        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- getCellMetadata(gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "cellMetaObj",
                copy_obj = TRUE
            )

            cell_metadata[][, eval(name) := NULL]

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setCellMetadata(gobject,
                x = cell_metadata,
                verbose = FALSE,
                initialize = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }


        # add new metadata information
        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = ident_clusters_DT[
                , c("cell_ID", name), with = FALSE],
            by_column = TRUE, column_cell_ID = "cell_ID"
        )

        ## update parameters used ##

        # 1. get parent function name
        cl <- sys.call(-1)

        # 2. check if this function call is within doLouvainCluster
        if (is.null(cl)) {
            gobject <- update_giotto_params(gobject, description = "_cluster")
        } else {
            fname <- as.character(cl[[1]])
            if (fname == "doLouvainCluster") {
                gobject <- update_giotto_params(gobject,
                    description = "_cluster",
                    toplevel = 3
                )
            } else {
                gobject <- update_giotto_params(gobject,
                    description = "_cluster"
                )
            }
        }

        return(gobject)
    } else {
        # else return clustering result
        return(ident_clusters_DT)
    }
}


#' @title doLouvainCluster multinet
#' @name .doLouvainCluster_multinet
#' @description cluster cells using a NN-network and the Louvain algorithm from
#' the multinet package in R.
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name for cluster, default to "louvain_clus"
#' @param nn_network_to_use type of NN network to use (kNN vs sNN), default to
#' "sNN"
#' @param network_name name of NN network to use, default to "sNN.pca"
#' @param gamma Resolution parameter for modularity in the generalized louvain
#' method. default  = 1
#' @param omega Inter-layer weight parameter in the generalized louvain method.
#' default = 1
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed (default = FALSE)
#' @param seed_number number for seed
#' @returns giotto object with new clusters appended to cell metadata
#' @details See \code{\link[multinet]{glouvain_ml}} from the multinet package
#' in R for more information.
#'
#' @keywords internal
.doLouvainCluster_multinet <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    name = "louvain_clus",
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    gamma = 1,
    omega = 1,
    return_gobject = TRUE,
    set_seed = FALSE,
    seed_number = 1234) {
    if ("multinet" %in% rownames(installed.packages()) == FALSE) {
        stop(
            "package 'multinet' is not yet installed \n",
            "To install: \n",
            "install.packages('multinet')"
        )
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

    ## get cell IDs ##
    cell_ID_vec <- gobject@cell_ID[[spat_unit]]

    ## select network to use
    igraph_object <- getNearestNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        nn_type = nn_network_to_use,
        name = network_name,
        output = "igraph"
    )

    # create mlnetworkobject
    mln_object <- multinet::ml_empty()
    # multinet::add_vertices_ml(
    #     n = mln_object, vertices = igraph::V(igraph_object))
    multinet::add_igraph_layer_ml(
        n = mln_object, g = igraph_object, name = name
    )

    # start seed
    if (isTRUE(set_seed)) {
        set.seed(seed = as.integer(seed_number))
    }

    # data.table variables
    cell_ID <- actor <- weight_col <- NULL

    louvain_clusters <- multinet::glouvain_ml(
        n = mln_object, gamma = gamma, omega = omega
    )
    ident_clusters_DT <- data.table::as.data.table(louvain_clusters)
    ident_clusters_DT[, cell_ID := actor]
    data.table::setnames(ident_clusters_DT, "cid", name)

    # exit seed
    if (isTRUE(set_seed)) {
        set.seed(Sys.time())
    }


    ## return
    if (isTRUE(return_gobject)) {
        # get cell metadata names
        cluster_names <- names(pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        ))


        # set name/cluster results to NULL if already exist
        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- getCellMetadata(gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "cellMetaObj",
                copy_obj = TRUE
            )

            cell_metadata[][, eval(name) := NULL]

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setCellMetadata(gobject,
                x = cell_metadata,
                verbose = FALSE,
                initialize = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }


        # add new metadata information
        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = ident_clusters_DT[
                , c("cell_ID", name), with = FALSE],
            by_column = TRUE, column_cell_ID = "cell_ID"
        )

        ## update parameters used ##
        # 1. get parent function name
        cl <- sys.call(-1)

        # 2. check if this function call is within doLouvainCluster
        if (is.null(cl)) {
            gobject <- update_giotto_params(gobject, description = "_cluster")
        } else {
            fname <- as.character(cl[[1]])
            if (fname == "doLouvainCluster") {
                gobject <- update_giotto_params(gobject,
                    description = "_cluster",
                    toplevel = 3
                )
            } else {
                gobject <- update_giotto_params(gobject,
                    description = "_cluster"
                )
            }
        }
        return(gobject)
    } else {
        # else return clustering result
        return(ident_clusters_DT)
    }
}


#' @title doLouvainCluster
#' @name doLouvainCluster
#' @description cluster cells using a NN-network and the Louvain algorithm.
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param version implemented version of Louvain clustering to use
#' @param name name for cluster, default to "louvain_clus"
#' @param nn_network_to_use type of NN network to use (kNN vs sNN), default to
#' "sNN"
#' @param network_name name of NN network to use, default to "sNN.pca"
#' @param python_path [community] specify specific path to python if required
#' @param resolution [community] resolution, default = 1
#' @param louv_random [community] Will randomize the node evaluation order and
#' the community evaluation order to get different partitions at each call
#' (default = FALSE)
#' @param weight_col weight column name
#' @param gamma [multinet] Resolution parameter for modularity in the
#' generalized louvain method, default = 1
#' @param omega [multinet] Inter-layer weight parameter in the generalized
#' louvain method, default = 1
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed (default = FALSE)
#' @param ... arguments passed to \code{\link{.doLouvainCluster_community}} or
#' \code{\link{.doLouvainCluster_multinet}}
#' @param seed_number number for seed
#'
#' @returns giotto object with new clusters appended to cell metadata
#' @details Louvain clustering using the community or multinet implementation
#' of the louvain clustering algorithm.
#' @seealso \code{\link{.doLouvainCluster_community}} and
#' \code{\link{.doLouvainCluster_multinet}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' doLouvainCluster(g)
#' @export
doLouvainCluster <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    version = c("community", "multinet"),
    name = "louvain_clus",
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    python_path = NULL,
    resolution = 1,
    weight_col = NULL,
    gamma = 1,
    omega = 1,
    louv_random = FALSE,
    return_gobject = TRUE,
    set_seed = FALSE,
    seed_number = 1234,
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



    ## louvain clustering version to use
    version <- match.arg(version, c("community", "multinet"))

    # python community implementation
    if (version == "community") {
        result <- .doLouvainCluster_community(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = name,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            python_path = python_path,
            resolution = resolution,
            weight_col = weight_col,
            louv_random = louv_random,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number,
            ...
        )

        return(result)

        ## r multinet implementation
    } else if (version == "multinet") {
        result <- .doLouvainCluster_multinet(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = name,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            gamma = gamma,
            omega = omega,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )

        return(result)
    }
}



#' @title doRandomWalkCluster
#' @name doRandomWalkCluster
#' @description Cluster cells using a random walk approach.
#' @param gobject giotto object
#' @param name name for cluster, default to "random_walk_clus"
#' @param nn_network_to_use type of NN network to use (kNN vs sNN), default to
#' "sNN"
#' @param network_name name of NN network to use, default to "sNN.pca"
#' @param walk_steps number of walking steps, default = 4
#' @param walk_clusters number of final clusters, default =  10
#' @param walk_weights cluster column defining the walk weights
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed (default = FALSE)
#' @param seed_number number for seed
#' @returns giotto object with new clusters appended to cell metadata
#' @details See \code{\link[igraph]{cluster_walktrap}} function from the igraph
#' package in R for more information.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' g <- doRandomWalkCluster(g)
#' pDataDT(g)
#' @export
doRandomWalkCluster <- function(gobject,
    name = "random_walk_clus",
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    walk_steps = 4,
    walk_clusters = 10,
    walk_weights = NA,
    return_gobject = TRUE,
    set_seed = FALSE,
    seed_number = 1234) {
    ## get cell IDs ##
    cell_ID_vec <- gobject@cell_ID

    ## select network to use
    igraph_object <- getNearestNetwork(
        gobject,
        nn_type = nn_network_to_use,
        name = network_name,
        output = "igraph"
    )


    # start seed
    if (isTRUE(set_seed)) {
        set.seed(seed = seed_number)
    }

    randomwalk_clusters <- igraph::cluster_walktrap(
        graph = igraph_object, steps = walk_steps, weights = walk_weights
    )
    randomwalk_clusters <- as.factor(igraph::cut_at(
        communities = randomwalk_clusters, no = walk_clusters
    ))

    ident_clusters_DT <- data.table::data.table(
        "cell_ID" = igraph::V(igraph_object)$name,
        "name" = randomwalk_clusters
    )
    data.table::setnames(ident_clusters_DT, "name", name)

    # exit seed
    if (isTRUE(set_seed)) {
        set.seed(Sys.time())
    }

    ## return
    if (return_gobject == TRUE) {
        gobject <- addCellMetadata(
            gobject = gobject,
            new_metadata = ident_clusters_DT[, c("cell_ID", name),
                with = FALSE
            ],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )


        ## update parameters used ##
        gobject <- update_giotto_params(gobject,
            description = "_randomwalk_cluster"
        )
        return(gobject)
    } else {
        # else return clustering result
        return(ident_clusters_DT)
    }
}


#' @title doSNNCluster
#' @name doSNNCluster
#' @description Cluster cells using a SNN cluster approach.
#' @param gobject giotto object
#' @param name name for cluster, default to "sNN_clus"
#' @param nn_network_to_use type of NN network to use (only works on kNN),
#' default to "kNN"
#' @param network_name name of kNN network to use, default to "kNN.pca"
#' @param k Neighborhood size for nearest neighbor sparsification to create
#' the shared NN graph, default = 20
#' @param eps Two objects are only reachable from each other if they share at
#' least eps nearest neighbors, default = 4
#' @param minPts minimum number of points that share at least eps nearest
#' neighbors for a point to be considered a core points, default = 16
#' @param borderPoints should borderPoints be assigned to clusters like in
#' DBSCAN? (default = TRUE)
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed (default = FALSE)
#' @param seed_number number for seed
#' @returns giotto object with new clusters appended to cell metadata
#' @details See \code{\link[dbscan]{sNNclust}} from dbscan package
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- createNearestNetwork(g, type = "kNN")
#'
#' doSNNCluster(g)
#' @export
doSNNCluster <- function(gobject,
    name = "sNN_clus",
    nn_network_to_use = "kNN",
    network_name = "kNN.pca",
    k = 20,
    eps = 4,
    minPts = 16,
    borderPoints = TRUE,
    return_gobject = TRUE,
    set_seed = FALSE,
    seed_number = 1234) {
    ## get cell IDs ##
    cell_ID_vec <- gobject@cell_ID

    ## select network to use
    igraph_object <- getNearestNetwork(
        gobject,
        nn_type = nn_network_to_use,
        name = network_name,
        output = "igraph"
    )


    if (nn_network_to_use == "sNN") {
        stop("sNNclust can only be used with kNN-network")
    }


    # start seed
    if (isTRUE(set_seed)) {
        set.seed(seed = seed_number)
    }

    # data.table variables
    from <- from_T <- to_T <- to <- weight <- distance <- NULL

    ## SNN clust
    igraph_DT <- data.table::as.data.table(igraph::as_data_frame(
        igraph_object,
        what = "edges"
    ))
    igraph_DT <- igraph_DT[order(from)]

    cell_id_numeric <- unique(x = c(igraph_DT$from, igraph_DT$to))
    names(cell_id_numeric) <- seq_along(cell_id_numeric)
    igraph_DT[, from_T := as.numeric(names(cell_id_numeric[
        cell_id_numeric == from
    ])), by = 1:nrow(igraph_DT)]
    igraph_DT[, to_T := as.numeric(names(cell_id_numeric[
        cell_id_numeric == to
    ])), by = 1:nrow(igraph_DT)]
    temp_igraph_DT <- igraph_DT[, .(from_T, to_T, weight, distance)]
    data.table::setnames(
        temp_igraph_DT,
        old = c("from_T", "to_T"), new = c("from", "to")
    )

    kNN_object <- nnDT_to_kNN(nnDT = temp_igraph_DT)
    sNN_clusters <- dbscan::sNNclust(
        x = kNN_object, k = k, eps = eps,
        minPts = minPts, borderPoints = borderPoints
    )

    ident_clusters_DT <- data.table::data.table(
        "cell_ID" = cell_id_numeric[seq_len(nrow(kNN_object$dist))],
        "name" = sNN_clusters$cluster
    )
    data.table::setnames(ident_clusters_DT, "name", name)

    # exit seed
    if (isTRUE(set_seed)) {
        set.seed(Sys.time())
    }

    ## add clusters to metadata ##
    if (return_gobject == TRUE) {
        cluster_names <- names(gobject@cell_metadata)
        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- gobject@cell_metadata
            cell_metadata[, eval(name) := NULL]
            gobject@cell_metadata <- cell_metadata
        }

        gobject <- addCellMetadata(
            gobject = gobject,
            new_metadata = ident_clusters_DT[, c("cell_ID", name),
                with = FALSE
            ],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )


        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_SNN_cluster")
        return(gobject)
    } else {
        # else return clustering result
        return(ident_clusters_DT)
    }
}





#' @title doKmeans
#' @name doKmeans
#' @description cluster cells using kmeans algorithm
#' @param gobject giotto object
#' @param feat_type feature type (e.g. "cell")
#' @param spat_unit spatial unit (e.g. "rna", "dna", "protein")
#' @param expression_values expression values from list_expression()
#' (e.g. "normalized", "scaled", "custom")
#' @param feats_to_use (optional) subset of features to use
#' @param dim_reduction_to_use dimension reduction from list_dim_reductions()
#' (e.g. "pca", "umap", "tsne")
#' @param dim_reduction_name dimensions reduction name, default to "pca"
#' @param dimensions_to_use dimensions to use, default = 1:10
#' @param distance_method distance method (e.g. "original", "pearson",
#' "spearman", "euclidean", "maximum", "manhattan", "canberra", "binary",
#' "minkowski")
#' @param centers number of final clusters, default = 10
#' @param iter_max kmeans maximum iterations, default = 100
#' @param nstart kmeans nstart, default = 1000
#' @param algorithm kmeans algorithm, default to "Hartigan-Wong"
#' @param name name for kmeans clustering, default to "kmeans"
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed (default = TRUE)
#' @param seed_number number for seed
#' @returns if return_gobject = TRUE: giotto object with new clusters appended
#' to cell metadata
#' @details The default settings will use dimension reduction results as input.
#' Set dim_reduction_to_use = NULL if you want to directly use expression
#' values as input.
#' By providing a feature vector to feats_to_use you can subset the expression
#' matrix.
#' @seealso  \code{\link[stats]{kmeans}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' doKmeans(g)
#' @export
doKmeans <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    feats_to_use = NULL,
    dim_reduction_to_use = c("cells", "pca", "umap", "tsne"),
    dim_reduction_name = "pca",
    dimensions_to_use = 1:10,
    distance_method = c(
        "original", "pearson", "spearman",
        "euclidean", "maximum", "manhattan",
        "canberra", "binary", "minkowski"
    ),
    centers = 10,
    iter_max = 100,
    nstart = 1000,
    algorithm = "Hartigan-Wong",
    name = "kmeans",
    return_gobject = TRUE,
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

    distance_method <- match.arg(distance_method, choices = c(
        "original", "pearson", "spearman",
        "euclidean", "maximum", "manhattan",
        "canberra", "binary", "minkowski"
    ))
    dim_reduction_to_use <- match.arg(
        dim_reduction_to_use, c("cells", "pca", "umap", "tsne")
    )
    expression_values <- match.arg(
        expression_values, c("normalized", "scaled", "custom")
    )

    ## using dimension reduction ##
    if (dim_reduction_to_use != "cells" && !is.null(dim_reduction_to_use)) {
        # use only available dimensions if dimensions < dimensions_to_use
        dim_coord <- getDimReduction(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "dimObj"
        )

        dimensions_to_use <- dimensions_to_use[
            dimensions_to_use %in% seq_len(ncol(dim_coord[]))
        ]
        matrix_to_use <- dim_coord[][, dimensions_to_use]
    } else {
        vmsg(.is_debug = TRUE, "clustering from expression values")
        ## using original matrix ##
        expr_values <- getExpression(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            values = expression_values,
            output = "exprObj"
        )

        # subset expression matrix
        if (!is.null(feats_to_use)) {
            expr_values[] <- expr_values[][
                rownames(expr_values[]) %in% feats_to_use,
            ]
        }

        # features as columns
        # cells as rows
        matrix_to_use <- t_flex(expr_values[])
    }


    ## distance
    if (distance_method == "original") {
        celldist <- matrix_to_use
    } else if (distance_method %in% c("spearman", "pearson")) {
        celldist <- stats::as.dist(1 - cor_flex(
            x = t_flex(matrix_to_use), method = distance_method
        ))
    } else if (distance_method %in% c(
        "euclidean", "maximum", "manhattan",
        "canberra", "binary", "minkowski"
    )) {
        celldist <- stats::dist(x = matrix_to_use, method = distance_method)
    }

    ## kmeans clustering
    # start seed
    if (isTRUE(set_seed)) {
        GiottoUtils::local_seed(seed_number)
    }

    # start clustering
    kclusters <- stats::kmeans(
        x = celldist,
        centers = centers,
        iter.max = iter_max,
        nstart = nstart,
        algorithm = algorithm
    )

    ident_clusters_DT <- data.table::data.table(
        "cell_ID" = names(kclusters[["cluster"]]),
        "name" = kclusters[["cluster"]]
    )
    data.table::setnames(ident_clusters_DT, "name", name)


    ## add clusters to metadata ##
    if (isTRUE(return_gobject)) {
        cluster_names <- names(pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        ))

        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- getCellMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "cellMetaObj",
                copy_obj = TRUE
            )

            cell_metadata[][, eval(name) := NULL]

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setGiotto(
                gobject, cell_metadata,
                verbose = FALSE, initialize = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }

        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = ident_clusters_DT[
                , c("cell_ID", name), with = FALSE],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )

        ## update parameters used ##
        gobject <- update_giotto_params(gobject,
            description = "_kmeans_cluster"
        )
        return(gobject)
    } else {
        return(ident_clusters_DT)
    }
}





#' @title doHclust
#' @name doHclust
#' @description cluster cells using hierarchical clustering algorithm
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param feats_to_use subset of features to use
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimensions reduction name
#' @param dimensions_to_use dimensions to use
#' @param distance_method distance method
#' @param agglomeration_method agglomeration method for hclust
#' @param k number of final clusters
#' @param h cut hierarchical tree at height = h
#' @param name name for hierarchical clustering
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @returns giotto object with new clusters appended to cell metadata
#' @details Description on how to use Kmeans clustering method.
#' @seealso  \code{\link[stats]{hclust}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' doHclust(g)
#' @export
doHclust <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    feats_to_use = NULL,
    dim_reduction_to_use = c("cells", "pca", "umap", "tsne"),
    dim_reduction_name = "pca",
    dimensions_to_use = 1:10,
    distance_method = c(
        "pearson", "spearman", "original",
        "euclidean", "maximum", "manhattan",
        "canberra", "binary", "minkowski"
    ),
    agglomeration_method = c(
        "ward.D2", "ward.D", "single",
        "complete", "average", "mcquitty",
        "median", "centroid"
    ),
    k = 10,
    h = NULL,
    name = "hclust",
    return_gobject = TRUE,
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


    distance_method <- match.arg(
        distance_method,
        choices = c(
            "pearson", "spearman", "original",
            "euclidean", "maximum", "manhattan",
            "canberra", "binary", "minkowski"
        )
    )
    agglomeration_method <- match.arg(
        agglomeration_method,
        choices = c(
            "ward.D2", "ward.D", "single",
            "complete", "average", "mcquitty",
            "median", "centroid"
        )
    )
    values <- match.arg(expression_values, c("normalized", "scaled", "custom"))
    dim_reduction_to_use <- match.arg(
        dim_reduction_to_use, c("cells", "pca", "umap", "tsne")
    )


    ## using dimension reduction ##
    if (dim_reduction_to_use != "cells" && !is.null(dim_reduction_to_use)) {

        # use only available dimensions if dimensions < dimensions_to_use
        dim_coord <- getDimReduction(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "data.table"
        )

        dimensions_to_use <- dimensions_to_use[
            dimensions_to_use %in% seq_len(ncol(dim_coord))
        ]
        matrix_to_use <- dim_coord[, dimensions_to_use]
    } else {
        vmsg(.is_debug = TRUE, "clustering from expression values")
        ## using original matrix ##
        expr_values <- getExpression(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            values = values,
            output = "matrix"
        )

        # subset expression matrix
        if (!is.null(feats_to_use)) {
            expr_values <- expr_values[
                rownames(expr_values) %in% feats_to_use,
            ]
        }

        # features as columns
        # cells as rows
        matrix_to_use <- t_flex(expr_values)
    }

    ## distance
    if (distance_method == "original") {
        celldist <- matrix_to_use
    } else if (distance_method %in% c("spearman", "pearson")) {
        celldist <- stats::as.dist(1 - cor_flex(x = t_flex(
            matrix_to_use
        ), method = distance_method))
    } else if (distance_method %in% c(
        "euclidean", "maximum", "manhattan",
        "canberra", "binary", "minkowski"
    )) {
        celldist <- stats::dist(x = matrix_to_use, method = distance_method)
    }

    ## hierarchical clustering
    # start seed
    if (isTRUE(set_seed)) {
        set.seed(seed = as.integer(seed_number))
    }

    # start clustering
    hclusters <- stats::hclust(d = celldist, method = agglomeration_method)
    hclusters_cut <- stats::cutree(tree = hclusters, k = k, h = h)

    # exit seed
    if (isTRUE(set_seed)) {
        set.seed(seed = Sys.time())
    }

    ident_clusters_DT <- data.table::data.table(
        cell_ID = names(hclusters_cut),
        "name" = hclusters_cut
    )
    data.table::setnames(ident_clusters_DT, "name", name)


    ## add clusters to metadata ##
    if (return_gobject == TRUE) {
        cluster_names <- names(pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        ))

        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- getCellMetadata(gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "cellMetaObj",
                copy_obj = TRUE
            )

            cell_metadata[][, eval(name) := NULL]

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setCellMetadata(gobject,
                x = cell_metadata,
                verbose = FALSE,
                initialize = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }

        gobject <- addCellMetadata(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            new_metadata = ident_clusters_DT[
                , c("cell_ID", name), with = FALSE],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )



        ## update parameters used ##
        gobject <- update_giotto_params(gobject,
            description = "_hierarchical_cluster"
        )
        return(gobject)
    } else {
        return(list("hclust" = hclusters, "DT" = ident_clusters_DT))
    }
}





#' @title clusterCells
#' @name clusterCells
#' @description cluster cells using a variety of different methods
#' @param gobject giotto object
#' @param name name for new clustering result
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param cluster_method community cluster method to use
#' @param pyth_leid_resolution resolution for leiden
#' @param pyth_leid_weight_col column to use for weights
#' @param pyth_leid_part_type partition type to use
#' @param pyth_leid_init_memb initial membership
#' @param pyth_leid_iterations number of iterations
#' @param pyth_louv_resolution resolution for louvain
#' @param pyth_louv_weight_col python louvain param: weight column
#' @param python_louv_random python louvain param: random
#' @param python_path specify specific path to python if required
#' @param louvain_gamma louvain param: gamma or resolution
#' @param louvain_omega louvain param: omega
#' @param walk_steps randomwalk: number of steps
#' @param walk_clusters randomwalk: number of clusters
#' @param walk_weights randomwalk: weight column
#' @param sNNclust_k SNNclust: k neighbors to use
#' @param sNNclust_eps SNNclust: epsilon
#' @param sNNclust_minPts SNNclust: min points
#' @param borderPoints SNNclust: border points
#' @param expression_values expression values to use
#' @param feats_to_use features to use in clustering,
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name name of reduction 'pca',
#' @param dimensions_to_use dimensions to use
#' @param distance_method distance method
#' @param km_centers kmeans centers
#' @param km_iter_max kmeans iterations
#' @param km_nstart kmeans random starting points
#' @param km_algorithm kmeans algorithm
#' @param hc_agglomeration_method hierarchical clustering method
#' @param hc_k hierachical number of clusters
#' @param hc_h hierarchical tree cutoff
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @returns giotto object with new clusters appended to cell metadata
#' @details Wrapper for the different clustering methods.
#' @seealso \code{\link{doLeidenCluster}},
#' \code{\link{.doLouvainCluster_community}},
#' \code{\link{.doLouvainCluster_multinet}},
#' \code{\link{doLouvainCluster}}, \code{\link{doRandomWalkCluster}},
#' \code{\link{doSNNCluster}},
#' \code{\link{doKmeans}}, \code{\link{doHclust}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' clusterCells(g)
#' @export
clusterCells <- function(gobject,
    cluster_method = c(
        "leiden",
        "louvain_community", "louvain_multinet",
        "randomwalk", "sNNclust",
        "kmeans", "hierarchical"
    ),
    name = "cluster_name",
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    pyth_leid_resolution = 1,
    pyth_leid_weight_col = "weight",
    pyth_leid_part_type = c(
        "RBConfigurationVertexPartition",
        "ModularityVertexPartition"
    ),
    pyth_leid_init_memb = NULL,
    pyth_leid_iterations = 1000,
    pyth_louv_resolution = 1,
    pyth_louv_weight_col = NULL,
    python_louv_random = FALSE,
    python_path = NULL,
    louvain_gamma = 1,
    louvain_omega = 1,
    walk_steps = 4,
    walk_clusters = 10,
    walk_weights = NA,
    sNNclust_k = 20,
    sNNclust_eps = 4,
    sNNclust_minPts = 16,
    borderPoints = TRUE,
    expression_values = c("normalized", "scaled", "custom"),
    feats_to_use = NULL,
    dim_reduction_to_use = c("cells", "pca", "umap", "tsne"),
    dim_reduction_name = "pca",
    dimensions_to_use = 1:10,
    distance_method = c(
        "original", "pearson", "spearman",
        "euclidean", "maximum", "manhattan",
        "canberra", "binary", "minkowski"
    ),
    km_centers = 10,
    km_iter_max = 100,
    km_nstart = 1000,
    km_algorithm = "Hartigan-Wong",
    hc_agglomeration_method = c(
        "ward.D2", "ward.D", "single",
        "complete", "average", "mcquitty",
        "median", "centroid"
    ),
    hc_k = 10,
    hc_h = NULL,
    return_gobject = TRUE,
    set_seed = TRUE,
    seed_number = 1234) {
    ## select cluster method
    cluster_method <- match.arg(
        arg = cluster_method,
        choices = c(
            "leiden",
            "louvain_community",
            "louvain_multinet",
            "randomwalk",
            "sNNclust",
            "kmeans",
            "hierarchical"
        )
    )


    if (cluster_method == "leiden") {
        result <- doLeidenCluster(
            gobject = gobject,
            name = name,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            python_path = python_path,
            resolution = pyth_leid_resolution,
            weight_col = pyth_leid_weight_col,
            partition_type = pyth_leid_part_type,
            init_membership = pyth_leid_init_memb,
            n_iterations = pyth_leid_iterations,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )
    } else if (cluster_method == "louvain_community") {
        result <- .doLouvainCluster_community(
            gobject = gobject,
            name = name,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            python_path = python_path,
            resolution = pyth_louv_resolution,
            weight_col = pyth_louv_weight_col,
            louv_random = python_louv_random,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )
    } else if (cluster_method == "louvain_multinet") {
        result <- .doLouvainCluster_multinet(
            gobject = gobject,
            name = name,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            gamma = louvain_gamma,
            omega = louvain_omega,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )
    } else if (cluster_method == "randomwalk") {
        result <- doRandomWalkCluster(
            gobject = gobject,
            name = name,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            walk_steps = walk_steps,
            walk_clusters = walk_clusters,
            walk_weights = walk_weights,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )
    } else if (cluster_method == "sNNclust") {
        result <- doSNNCluster(
            gobject = gobject,
            name = name,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            k = sNNclust_k,
            eps = sNNclust_eps,
            minPts = sNNclust_minPts,
            borderPoints = borderPoints,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )
    } else if (cluster_method == "kmeans") {
        result <- doKmeans(
            gobject = gobject,
            name = name,
            expression_values = expression_values,
            feats_to_use = feats_to_use,
            dim_reduction_to_use = dim_reduction_to_use,
            dim_reduction_name = dim_reduction_name,
            dimensions_to_use = dimensions_to_use,
            distance_method = distance_method,
            centers = km_centers,
            iter_max = km_iter_max,
            nstart = km_nstart,
            algorithm = km_algorithm,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )
    } else if (cluster_method == "hierarchical") {
        result <- doHclust(
            gobject = gobject,
            name = name,
            expression_values = expression_values,
            feats_to_use = feats_to_use,
            dim_reduction_to_use = dim_reduction_to_use,
            dim_reduction_name = dim_reduction_name,
            dimensions_to_use = dimensions_to_use,
            distance_method = distance_method,
            agglomeration_method = hc_agglomeration_method,
            k = hc_k,
            h = hc_h,
            return_gobject = return_gobject,
            set_seed = set_seed,
            seed_number = seed_number
        )
    }


    return(result)
}




# subclustering ####


#' @title Cell subclustering
#' @name subClusterCells
#' @description Perform cell subclustering by taking an annotated group of
#' cells and performing another round of clustering on just that subset.
#' Several methods are implemented. `subClusterCells()` is the main wrapper
#' function. `doLeidenSubCluster()` and `doLouvainSubCluster()` are more
#' specific implementations.
#' @param gobject `giotto` object
#' @param name name for new clustering result
#' @param cluster_method clustering method to use. Currently one of "leiden"
#' (default), "louvain_community", "louvain_multinet"
#' @param cluster_column cluster column to subcluster
#' @param selected_clusters only do subclustering on these clusters
#' @param hvf_param list of parameters for [calculateHVF()]
#' @param hvg_param deprecated
#' @param hvf_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_min_perc_cells deprecated
#' @param hvf_mean_expr_det threshold for mean expression level in cells with
#' detection
#' @param hvg_mean_expr_det deprecated
#' @param use_all_feats_as_hvf forces all features to be HVF and to be used as
#' input for PCA
#' @param use_all_genes_as_hvg deprecated
#' @param min_nr_of_hvf minimum number of HVF, or all features will be used as
#' input for PCA
#' @param min_nr_of_hvg deprecated
#' @param pca_param list of parameters for [runPCA()]
#' @param nn_param list of parameters for [createNearestNetwork()]
#' @param k_neighbors number of k for [createNearestNetwork()]
#' @param resolution resolution for community algorithm
#' @param n_iterations number of iterations to run the Leiden algorithm.
#' @param gamma gamma
#' @param omega omega
#' @param python_path specify specific path to python if required
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param return_gobject logical. return `giotto` object (default = TRUE)
#' @param verbose verbose
#' @returns `giotto` object with new subclusters appended to cell metadata
#' @details This function performs subclustering on selected clusters.
#' The systematic steps are:
#'   1. subset Giotto object
#'   2. identify highly variable genes
#'   3. run PCA
#'   4. create nearest neighbouring network
#'   5. do clustering
#'
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' # Run some subclusterings based on "leiden_clus" annotations that already
#' # exist in the visium mini object
#'
#' # default method is leiden subclustering
#' subClusterCells(g, cluster_column = "leiden_clus")
#'
#' # use louvain instead
#' subClusterCells(g,
#'     cluster_column = "leiden_clus",
#'     cluster_method = "louvain_community"
#' )
#'
#' # directly call the more specific functions
#' doLeidenSubCluster(g, cluster_column = "leiden_clus")
#'
#' doLouvainSubCluster(g, cluster_column = "leiden_clus")
#' @md
NULL




#' @rdname subClusterCells
#' @export
subClusterCells <- function(
        gobject,
        name = "sub_clus",
        cluster_method = c(
            "leiden",
            "louvain_community",
            "louvain_multinet"
        ),
        cluster_column = NULL,
        selected_clusters = NULL,
        hvg_param = deprecated(),
        hvf_param = list(
            reverse_log_scale = TRUE, difference_in_cov = 1,
            expression_values = "normalized"
        ),
        hvg_min_perc_cells = deprecated(),
        hvf_min_perc_cells = 5,
        hvg_mean_expr_det = deprecated(),
        hvf_mean_expr_det = 1,
        use_all_genes_as_hvg = deprecated(),
        use_all_feats_as_hvf = FALSE,
        min_nr_of_hvg = deprecated(),
        min_nr_of_hvf = 5,
        pca_param = list(expression_values = "normalized", scale_unit = TRUE),
        nn_param = list(dimensions_to_use = 1:20),
        k_neighbors = 10,
        resolution = 1,
        n_iterations = 1000,
        gamma = 1,
        omega = 1,
        python_path = NULL,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        return_gobject = TRUE,
        verbose = TRUE) {
    ## select cluster method
    cluster_method <- match.arg(arg = cluster_method, choices = c(
        "leiden",
        "louvain_community",
        "louvain_multinet"
    ))

    # deprecations
    .dep_param <- function(...) {
        GiottoUtils::deprecate_param(
            ...,
            fun = "subClusterCells", when = "4.0.9"
        )
    }

    hvf_param <- .dep_param(hvg_param, hvf_param)
    hvf_min_perc_cells <- .dep_param(hvg_min_perc_cells, hvf_min_perc_cells)
    hvf_mean_expr_det <- .dep_param(hvg_mean_expr_det, hvf_mean_expr_det)
    use_all_feats_as_hvf <- .dep_param(use_all_genes_as_hvg,
                                       use_all_feats_as_hvf)
    min_nr_of_hvf <- .dep_param(min_nr_of_hvg, min_nr_of_hvf)

    # gather common args
    common_args <- get_args_list(keep = c(
        "gobject",
        "cluster_column",
        "selected_clusters",
        "hvf_param",
        "hvf_min_perc_cells",
        "hvf_mean_expr_det",
        "use_all_feats_as_hvf",
        "min_nr_of_hvf",
        "pca_param",
        "nn_param",
        "k_neighbors",
        "nn_network_to_use",
        "network_name",
        "name",
        "return_gobject",
        "verbose"
    ))

    result <- switch(cluster_method,
        "leiden" = {
            do.call(doLeidenSubCluster, args = c(
                common_args,
                list(
                    resolution = resolution,
                    n_iterations = n_iterations,
                    python_path = python_path,
                    toplevel = 4
                )
            ))
        },
        "louvain_community" = {
            do.call(.doLouvainSubCluster_community, args = c(
                common_args,
                list(
                    resolution = resolution,
                    python_path = python_path
                )
            ))
        },
        "louvain_multinet" = {
            do.call(.doLouvainSubCluster_multinet, args = c(
                common_args,
                list(
                    gamma = gamma,
                    omega = omega
                )
            ))
        }
    )

    return(result)
}





#' @describeIn subClusterCells Further subcluster cells using a NN-network and
#' the Leiden algorithm
#' @param toplevel do not use
#' @param feat_type feature type
#' @export
doLeidenSubCluster <- function(gobject,
    feat_type = NULL,
    name = "sub_leiden_clus",
    cluster_column = NULL,
    selected_clusters = NULL,
    hvf_param = list(
        reverse_log_scale = TRUE, difference_in_cov = 1,
        expression_values = "normalized"
    ),
    hvg_param = deprecated(),
    hvf_min_perc_cells = 5,
    hvg_min_perc_cells = deprecated(),
    hvf_mean_expr_det = 1,
    hvg_mean_expr_det = deprecated(),
    use_all_feats_as_hvf = FALSE,
    use_all_genes_as_hvg = deprecated(),
    min_nr_of_hvf = 5,
    min_nr_of_hvg = deprecated(),
    pca_param = list(expression_values = "normalized", scale_unit = TRUE),
    nn_param = list(dimensions_to_use = 1:20),
    k_neighbors = 10,
    resolution = 0.5,
    n_iterations = 500,
    python_path = NULL,
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    return_gobject = TRUE,
    toplevel = 2,
    verbose = TRUE) {
    # specify feat_type
    if (is.null(feat_type)) {
        feat_type <- gobject@expression_feat[[1]]
    }

    # deprecated arguments
    .dep_param <- function(x, y) {
        GiottoUtils::deprecate_param(
            x, y,
            fun = "doLeidenSubCluster", when = "4.0.9"
        )
    }

    hvf_param <- .dep_param(hvg_param, hvf_param)
    hvf_min_perc_cells <- .dep_param(hvg_min_perc_cells, hvf_min_perc_cells)
    hvf_mean_expr_det <- .dep_param(hvg_mean_expr_det, hvf_mean_expr_det)
    use_all_feats_as_hvf <- .dep_param(use_all_genes_as_hvg,
                                       use_all_feats_as_hvf)
    min_nr_of_hvf <- .dep_param(min_nr_of_hvg, min_nr_of_hvf)


    iter_list <- list()

    cell_metadata <- pDataDT(gobject,
        feat_type = feat_type
    )

    if (is.null(cluster_column)) {
        stop("You need to provide a cluster column to subcluster on")
    }
    unique_clusters <- mixedsort(unique(cell_metadata[[cluster_column]]))


    # data.table variables
    hvf <- perc_cells <- mean_expr_det <- parent_cluster <- comb <-
        tempclus <- NULL


    for (cluster in unique_clusters) {
        vmsg(.v = verbose, "start with cluster: ", cluster, "\n")

        ## get subset
        subset_cell_IDs <- cell_metadata[
            get(cluster_column) == cluster
        ][["cell_ID"]]
        temp_giotto <- subsetGiotto(
            gobject = gobject,
            feat_type = feat_type,
            cell_ids = subset_cell_IDs
        )

        ## if cluster is not selected
        if (!is.null(selected_clusters) & !cluster %in% selected_clusters) {
            temp_cluster <- data.table(
                "cell_ID" = subset_cell_IDs,
                "tempclus" = 1,
                "parent_cluster" = cluster
            )
            iter_list[[cluster]] <- temp_cluster
        } else {
            # continue for selected clusters or all clusters if there is no
            # selection

            ## calculate stats
            temp_giotto <- addStatistics(
                gobject = temp_giotto,
                feat_type = feat_type,
                verbose = FALSE
            )

            ## calculate variable feats
            hvf_param$verbose <- FALSE
            temp_giotto <- do.call(
                "calculateHVF", c(gobject = temp_giotto, hvf_param)
            )

            ## get hvg
            feat_metadata <- fDataDT(temp_giotto,
                feat_type = feat_type
            )
            usefeats <- feat_metadata[
                hvf == "yes" & perc_cells >= hvf_min_perc_cells &
                    mean_expr_det >= hvf_mean_expr_det
            ]$feat_ID

            ## catch too low number of hvg
            if (use_all_feats_as_hvf == TRUE) {
                usefeats == feat_metadata$feat_ID
            } else {
                if (verbose == TRUE) {
                    cat(
                        length(usefeats),
                        "highly variable feats have been selected\n"
                    )
                }
                if (length(usefeats) <= min_nr_of_hvf) {
                    message("too few feats, will continue with all feats
                            instead")
                    usefeats <- feat_metadata$feat_ID
                }
            }

            ## run PCA
            pca_param$verbose <- FALSE
            temp_giotto <- do.call(
                "runPCA",
                c(
                    gobject = temp_giotto, feats_to_use = list(usefeats),
                    pca_param
                )
            )

            ## nearest neighbor and clustering
            nn_param$verbose <- FALSE
            temp_giotto <- do.call(
                "createNearestNetwork",
                c(gobject = temp_giotto, k = k_neighbors, nn_param)
            )

            ## Leiden Cluster
            ## TO DO: expand to all clustering options
            temp_cluster <- doLeidenCluster(
                gobject = temp_giotto,
                feat_type = feat_type,
                resolution = resolution,
                n_iterations = n_iterations,
                python_path = python_path,
                name = "tempclus",
                return_gobject = FALSE
            )

            temp_cluster[, parent_cluster := cluster]

            iter_list[[cluster]] <- temp_cluster
        }
    }

    together <- do.call("rbind", iter_list)
    together[, comb := paste0(parent_cluster, ".", tempclus)]

    # rename with subcluster of original name
    setnames(together, "comb", name)


    if (return_gobject == TRUE) {
        cluster_names <- names(gobject@cell_metadata[[feat_type]])
        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- gobject@cell_metadata[[feat_type]]
            cell_metadata[, eval(name) := NULL]
            gobject@cell_metadata[[feat_type]] <- cell_metadata
        }

        gobject <- addCellMetadata(gobject,
            feat_type = feat_type,
            new_metadata = together[, c("cell_ID", name), with = FALSE],
            by_column = TRUE, column_cell_ID = "cell_ID"
        )

        ## update parameters used ##
        gobject <- update_giotto_params(
            gobject,
            description = "_sub_cluster", toplevel = toplevel
        )
        return(gobject)
    } else {
        return(together)
    }
}


# subcluster cells using a NN-network and the Louvain community
# detection algorithm
.doLouvainSubCluster_community <- function(gobject,
    name = "sub_louvain_comm_clus",
    cluster_column = NULL,
    selected_clusters = NULL,
    hvf_param = list(
        reverse_log_scale = TRUE,
        difference_in_cov = 1,
        expression_values = "normalized"
    ),
    hvf_min_perc_cells = 5,
    hvf_mean_expr_det = 1,
    use_all_feats_as_hvf = FALSE,
    min_nr_of_hvf = 5,
    pca_param = list(expression_values = "normalized", scale_unit = TRUE),
    nn_param = list(dimensions_to_use = 1:20),
    k_neighbors = 10,
    resolution = 0.5,
    python_path = NULL,
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    return_gobject = TRUE,
    verbose = TRUE) {
    iter_list <- list()

    cell_metadata <- pDataDT(gobject)

    if (is.null(cluster_column)) {
        stop("You need to provide a cluster column to subcluster on")
    }
    unique_clusters <- mixedsort(unique(cell_metadata[[cluster_column]]))

    ## if clusters start at 0, then add +1 for the index ##
    index_offset <- ifelse(0 %in% unique_clusters, 1, 0)

    for (cluster in unique_clusters) {
        if (verbose == TRUE) cat("start with cluster: ", cluster, "\n")

        ## get subset
        subset_cell_IDs <- cell_metadata[
            get(cluster_column) == cluster
        ][["cell_ID"]]
        temp_giotto <- subsetGiotto(
            gobject = gobject,
            cell_ids = subset_cell_IDs
        )

        ## if cluster is not selected
        if (!is.null(selected_clusters) & !cluster %in% selected_clusters) {
            temp_cluster <- data.table(
                "cell_ID" = subset_cell_IDs,
                "tempclus" = 1,
                "parent_cluster" = cluster
            )
            iter_list[[cluster + index_offset]] <- temp_cluster
        } else {
            # continue for selected clusters or all clusters if there is no
            # selection

            ## calculate stats
            temp_giotto <- addStatistics(
                gobject = temp_giotto, verbose = FALSE
            )

            ## calculate variable genes
            hvf_param$verbose <- FALSE
            temp_giotto <- do.call(
                "calculateHVF", c(gobject = temp_giotto, hvf_param)
            )

            ## get hvf
            feat_metadata <- fDataDT(temp_giotto)

            # NSE variables
            hvf <- perc_cells <- mean_expr_det <- NULL

            usefeats <- feat_metadata[
                hvf == "yes" &
                    perc_cells >= hvf_min_perc_cells &
                    mean_expr_det >= hvf_mean_expr_det
            ]$feat_ID

            ## catch too low number of hvf
            if (isTRUE(use_all_feats_as_hvf)) {
                usefeats == feat_metadata$feat_ID
            } else {
                if (isTRUE(verbose)) {
                    cat(
                        length(usefeats),
                        "highly variable features have been selected\n"
                    )
                }
                if (length(usefeats) <= min_nr_of_hvf) {
                    wrap_msg("too few features
                            will continue with all features instead")
                    usefeats <- feat_metadata$feat_ID
                }
            }

            ## run PCA
            pca_param$verbose <- FALSE
            temp_giotto <- do.call(
                "runPCA",
                c(
                    gobject = temp_giotto, feats_to_use = list(usefeats),
                    pca_param
                )
            )

            ## nearest neighbor and clustering
            nn_param$verbose <- FALSE
            temp_giotto <- do.call(
                "createNearestNetwork",
                c(gobject = temp_giotto, k = k_neighbors, nn_param)
            )

            ## TO DO: expand to all clustering options
            temp_cluster <- .doLouvainCluster_community(
                gobject = temp_giotto,
                resolution = resolution,
                python_path = python_path,
                name = "tempclus",
                return_gobject = FALSE
            )

            # data.table variables
            parent_cluster <- NULL

            temp_cluster[, parent_cluster := cluster]

            iter_list[[cluster + index_offset]] <- temp_cluster
        }
    }

    together <- do.call("rbind", iter_list)

    # data.table variables
    comb <- tempclus <- NULL

    together[, comb := paste0(parent_cluster, ".", tempclus)]

    # rename with subcluster of original name
    # new_cluster_column = paste0(cluster_column,'_sub')
    setnames(together, "comb", name)


    if (return_gobject == TRUE) {
        cluster_names <- names(gobject@cell_metadata)
        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- gobject@cell_metadata
            cell_metadata[, eval(name) := NULL]
            gobject@cell_metadata <- cell_metadata
        }

        gobject <- addCellMetadata(gobject,
            new_metadata = together[, c("cell_ID", name), with = FALSE],
            by_column = TRUE, column_cell_ID = "cell_ID"
        )

        ## update parameters used ##
        parameters_list <- gobject@parameters
        number_of_rounds <- length(parameters_list)
        update_name <- paste0(number_of_rounds, "_sub_cluster")

        # parameters to include
        parameters_list[[update_name]] <- c(
            "subclus name" = name,
            "resolution " = resolution,
            "k neighbors " = k_neighbors
        )

        gobject@parameters <- parameters_list

        return(gobject)
    } else {
        return(together)
    }
}




# subcluster cells using a NN-network and the Louvain multinet
# detection algorithm
.doLouvainSubCluster_multinet <- function(gobject,
    name = "sub_louvain_mult_clus",
    cluster_column = NULL,
    selected_clusters = NULL,
    hvf_param = list(
        reverse_log_scale = TRUE, difference_in_cov = 1,
        expression_values = "normalized"
    ),
    hvf_min_perc_cells = 5,
    hvf_mean_expr_det = 1,
    use_all_feats_as_hvf = FALSE,
    min_nr_of_hvf = 5,
    pca_param = list(expression_values = "normalized", scale_unit = TRUE),
    nn_param = list(dimensions_to_use = 1:20),
    k_neighbors = 10,
    gamma = 1,
    omega = 1,
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    return_gobject = TRUE,
    verbose = TRUE) {
    if ("multinet" %in% rownames(installed.packages()) == FALSE) {
        stop(
            "package 'multinet' is not yet installed \n",
            "To install: \n",
            "install.packages('multinet')"
        )
    }

    iter_list <- list()

    cell_metadata <- pDataDT(gobject)

    if (is.null(cluster_column)) {
        stop("You need to provide a cluster column to subcluster on")
    }
    unique_clusters <- mixedsort(unique(cell_metadata[[cluster_column]]))

    ## if clusters start at 0, then add +1 for the index ##
    index_offset <- ifelse(0 %in% unique_clusters, 1, 0)


    # data.table variables
    hvf <- perc_cells <- mean_expr_det <- parent_cluster <- cell_ID <-
        comb <- tempclus <- NULL


    for (cluster in unique_clusters) {
        if (verbose == TRUE) cat("start with cluster: ", cluster, "\n")

        ## get subset
        subset_cell_IDs <- cell_metadata[
            get(cluster_column) == cluster
        ][["cell_ID"]]
        temp_giotto <- subsetGiotto(
            gobject = gobject,
            cell_ids = subset_cell_IDs
        )

        ## if cluster is not selected
        if (!is.null(selected_clusters) & !cluster %in% selected_clusters) {
            temp_cluster <- data.table(
                "cell_ID" = subset_cell_IDs,
                "tempclus" = 1,
                "parent_cluster" = cluster
            )
            iter_list[[cluster + index_offset]] <- temp_cluster
        } else {
            # continue for selected clusters or all clusters if there is no
            # selection

            ## calculate stats
            temp_giotto <- addStatistics(
                gobject = temp_giotto, verbose = FALSE
            )

            ## calculate variable genes
            hvf_param$verbose <- FALSE
            temp_giotto <- do.call(
                "calculateHVF", c(gobject = temp_giotto, hvf_param)
            )

            ## get hvf
            feat_metadata <- fDataDT(temp_giotto)
            usefeats <- feat_metadata[
                hvf == "yes" & perc_cells >= hvf_min_perc_cells &
                    mean_expr_det >= hvf_mean_expr_det
            ]$feat_ID

            ## catch too low number of hvf
            if (use_all_feats_as_hvf == TRUE) {
                usefeats == feat_metadata$feat_ID
            } else {
                if (verbose == TRUE) {
                    cat(
                        length(usefeats),
                        "highly variable features have been selected\n"
                    )
                }
                if (length(usefeats) <= min_nr_of_hvf) {
                    message("too few features, will continue with all features
                            instead")
                    usefeats <- feat_metadata$feat_ID
                }
            }

            ## run PCA
            pca_param$verbose <- FALSE
            temp_giotto <- do.call(
                "runPCA",
                c(
                    gobject = temp_giotto, feats_to_use = list(usefeats),
                    pca_param
                )
            )

            ## nearest neighbor and clustering
            nn_param$verbose <- FALSE
            temp_giotto <- do.call(
                "createNearestNetwork",
                c(gobject = temp_giotto, k = k_neighbors, nn_param)
            )

            ## TO DO: expand to all clustering options
            temp_cluster <- .doLouvainCluster_multinet(
                gobject = temp_giotto,
                gamma = gamma,
                omega = omega,
                name = "tempclus",
                return_gobject = FALSE
            )

            temp_cluster[, parent_cluster := cluster]
            temp_cluster <- temp_cluster[, .(cell_ID, tempclus, parent_cluster)]

            iter_list[[cluster + index_offset]] <- temp_cluster
        }
    }

    together <- do.call("rbind", iter_list)
    together[, comb := paste0(parent_cluster, ".", tempclus)]

    # rename with subcluster of original name
    setnames(together, "comb", name)


    if (return_gobject == TRUE) {
        cluster_names <- names(gobject@cell_metadata)
        if (name %in% cluster_names) {
            cat(name, " has already been used, will be overwritten")
            cell_metadata <- gobject@cell_metadata
            cell_metadata[, eval(name) := NULL]
            gobject@cell_metadata <- cell_metadata
        }

        gobject <- addCellMetadata(
            gobject,
            new_metadata = together[, c("cell_ID", name), with = FALSE],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )

        ## update parameters used ##
        parameters_list <- gobject@parameters
        number_of_rounds <- length(parameters_list)
        update_name <- paste0(number_of_rounds, "_sub_cluster")

        # parameters to include
        parameters_list[[update_name]] <- c(
            "subclus name" = name,
            "k neighbors " = k_neighbors,
            "gamma" = gamma,
            "omega" = omega
        )

        gobject@parameters <- parameters_list

        return(gobject)
    } else {
        return(together)
    }
}




#' @describeIn subClusterCells subcluster cells using a NN-network and the
#' Louvain algorithm
#' @param version version of Louvain algorithm to use. One of "community" or
#' "multinet", with the default being "community"
#' @export
doLouvainSubCluster <- function(gobject,
    name = "sub_louvain_clus",
    version = c("community", "multinet"),
    cluster_column = NULL,
    selected_clusters = NULL,
    hvg_param = deprecated(),
    hvf_param = list(
        reverse_log_scale = TRUE, difference_in_cov = 1,
        expression_values = "normalized"
    ),
    hvg_min_perc_cells = deprecated(),
    hvf_min_perc_cells = 5,
    hvg_mean_expr_det = deprecated(),
    hvf_mean_expr_det = 1,
    use_all_genes_as_hvg = deprecated(),
    use_all_feats_as_hvf = FALSE,
    min_nr_of_hvg = deprecated(),
    min_nr_of_hvf = 5,
    pca_param = list(expression_values = "normalized", scale_unit = TRUE),
    nn_param = list(dimensions_to_use = 1:20),
    k_neighbors = 10,
    resolution = 0.5,
    gamma = 1,
    omega = 1,
    python_path = NULL,
    nn_network_to_use = "sNN",
    network_name = "sNN.pca",
    return_gobject = TRUE,
    verbose = TRUE) {
    ## louvain clustering version to use
    version <- match.arg(version, c("community", "multinet"))

    # deprecations
    .dep_param <- function(x, y) {
        GiottoUtils::deprecate_param(
            x, y,
            fun = "doLouvainSubCluster", when = "4.0.9"
        )
    }

    hvf_param <- .dep_param(hvg_param, hvf_param)
    hvf_min_perc_cells <- .dep_param(hvg_min_perc_cells, hvf_min_perc_cells)
    hvf_mean_expr_det <- .dep_param(hvg_mean_expr_det, hvf_mean_expr_det)
    use_all_feats_as_hvf <- .dep_param(use_all_genes_as_hvg,
                                       use_all_feats_as_hvf)
    min_nr_of_hvf <- .dep_param(min_nr_of_hvg, min_nr_of_hvf)

    # get common args
    common_args <- get_args_list(keep = c(
        "gobject",
        "cluster_column",
        "selected_clusters",
        "hvf_param",
        "hvf_min_perc_cells",
        "hvf_mean_expr_det",
        "use_all_feats_as_hvf",
        "min_nr_of_hvf",
        "pca_param",
        "nn_param",
        "k_neighbors",
        "nn_network_to_use",
        "network_name",
        "name",
        "return_gobject",
        "verbose"
    ))

    result <- switch(version,
        "community" = {
            do.call(.doLouvainSubCluster_community, args = c(
                common_args,
                list(
                    resolution = resolution,
                    python_path = python_path
                )
            ))
        },
        "multinet" = {
            do.call(.doLouvainSubCluster_multinet, args = c(
                common_args,
                list(
                    gamma = gamma,
                    omega = omega
                )
            ))
        }
    )

    return(result)
}








# cluster manipulation ####


#' @title getClusterSimilarity
#' @name getClusterSimilarity
#' @description Creates data.table with pairwise correlation scores between
#' each cluster.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @returns data.table
#' @details Creates data.table with pairwise correlation scores between each
#' cluster and the group size (# of cells) for each cluster. This information
#' can be used together with mergeClusters to combine very similar or small
#' clusters into bigger clusters.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' getClusterSimilarity(g, cluster_column = "leiden_clus")
#' @export
getClusterSimilarity <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    cluster_column,
    cor = c("pearson", "spearman")) {
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
    group1 <- group2 <- unified_group <- value <- NULL

    cor <- match.arg(cor, c("pearson", "spearman"))
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )

    metadata <- pDataDT(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit
    )

    # get clustersize
    clustersize <- metadata[, .N, by = cluster_column]
    colnames(clustersize) <- c("clusters", "size")

    # data.table variables
    clusters <- NULL

    clustersize[, clusters := as.character(clusters)]

    # scores per cluster
    metatable <- calculateMetaTable(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        metadata_cols = cluster_column
    )
    dcast_metatable <- data.table::dcast.data.table(
        metatable,
        formula = variable ~ uniq_ID, value.var = "value"
    )
    testmatrix <- dt_to_matrix(x = dcast_metatable)

    # correlation matrix
    cormatrix <- cor_flex(x = testmatrix, method = cor)
    cor_table <- melt_matrix(cormatrix)
    data.table::setnames(
        cor_table,
        old = c("Var1", "Var2"), c("group1", "group2")
    )
    cor_table[, c("group1", "group2") := list(
        as.character(group1), as.character(group2)
    )]
    cor_table[, unified_group := paste(
        sort(c(group1, group2)),
        collapse = "--"
    ),
    by = 1:nrow(cor_table)
    ]
    cor_table <- cor_table[!duplicated(cor_table[, .(value, unified_group)])]

    cor_table <- merge(
        cor_table,
        by.x = "group1", clustersize, by.y = "clusters"
    )
    setnames(cor_table, "size", "group1_size")
    cor_table <- merge(
        cor_table,
        by.x = "group2", clustersize, by.y = "clusters"
    )
    setnames(cor_table, "size", "group2_size")

    return(cor_table)
}




#' @title mergeClusters
#' @name mergeClusters
#' @description Merge selected clusters based on pairwise correlation scores
#' and size of cluster.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @param new_cluster_name new name for merged clusters
#' @param min_cor_score min correlation score to merge pairwise clusters
#' @param max_group_size max cluster size that can be merged
#' @param force_min_group_size size of clusters that will be merged with their
#' most similar neighbor(s)
#' @param max_sim_clusters maximum number of clusters to potentially merge to
#' reach force_min_group_size
#' @param return_gobject return giotto object
#' @param verbose be verbose
#' @returns Giotto object
#' @details Merge selected clusters based on pairwise correlation scores and
#' size of cluster.
#' To avoid large clusters to merge the max_group_size can be lowered. Small
#' clusters can be forcibly merged with their most similar pairwise cluster by
#' adjusting the force_min_group_size parameter. Clusters smaller than this
#' value will be merged independent on the provided min_cor_score value. The
#' force_min_group_size might not always be reached if clusters have already
#' been merged before \cr
#' A giotto object is returned by default, if FALSE then the merging vector
#' will be returned.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' mergeClusters(g, cluster_column = "leiden_clus")
#' @export
mergeClusters <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    cluster_column,
    cor = c("pearson", "spearman"),
    new_cluster_name = "merged_cluster",
    min_cor_score = 0.8,
    max_group_size = 20,
    force_min_group_size = 10,
    max_sim_clusters = 10,
    return_gobject = TRUE,
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

    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )

    # correlation score to be used
    cor <- match.arg(cor, c("pearson", "spearman"))

    # calculate similarity data.table
    similarityDT <- getClusterSimilarity(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        cluster_column = cluster_column,
        cor = cor
    )

    ## get clusters that can be merged
    # 1. clusters with high correlation

    # data.table variables
    group1 <- group2 <- group1_size <- value <- cumsum_val <- group2_size <-
        min_reached <- cumsum_reached <- NULL

    filter_set_first <- similarityDT[group1 != group2][
        group1_size < max_group_size
    ][value >= min_cor_score]

    # 2. small clusters
    minimum_set <- similarityDT[group1 != group2][
        group1_size < force_min_group_size
    ][order(-value)][
        , head(.SD, max_sim_clusters),
        by = group1
    ]

    # 2.1 take all clusters necessary to reach force_min_group_size
    minimum_set[, cumsum_val := cumsum(group2_size) + group1_size, by = group1]
    minimum_set[, min_reached := ifelse(cumsum_val > max_group_size, 1, 0)]
    minimum_set[, cumsum_reached := cumsum(min_reached), by = group1]
    minimum_set <- minimum_set[cumsum_reached <= 1]
    minimum_set[, c("cumsum_val", "min_reached", "cumsum_reached") := NULL]

    filter_set <- unique(do.call("rbind", list(filter_set_first, minimum_set)))

    ## get list of correlated groups
    finallist <- list()
    start_i <- 1
    for (row in seq_len(nrow(filter_set))) {
        first_clus <- filter_set[row][["group1"]]
        second_clus <- filter_set[row][["group2"]]

        res <- lapply(finallist, function(x) {
            any(x %in% c(first_clus, second_clus))
        })

        if (all(res == FALSE)) {
            finallist[[start_i]] <- c(first_clus, second_clus)
            start_i <- start_i + 1
        } else if (length(res[res == TRUE]) == 2) {
            NULL
        } else {
            who <- which(res == TRUE)[[1]]
            finallist[[who]] <- unique(
                c(finallist[[who]], first_clus, second_clus)
            )
        }
    }


    ## update metadata
    metadata <- data.table::copy(pDataDT(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    ))

    finalvec <- NULL
    for (ll in seq_along(finallist)) {
        tempvec <- finallist[[ll]]
        names(tempvec) <- rep(paste0("m_", ll), length(tempvec))
        finalvec <- c(finalvec, tempvec)
    }

    metadata[, eval(new_cluster_name) := ifelse(
        as.character(get(cluster_column)) %in% finalvec,
        names(finalvec[finalvec == as.character(get(cluster_column))]),
        as.character(get(cluster_column))
    ), by = 1:nrow(metadata)]




    if (return_gobject == TRUE) {
        cluster_names <- names(pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        ))


        if (new_cluster_name %in% cluster_names) {
            cat(new_cluster_name, " has already been used, will be overwritten")
            cell_metadata <- gobject@cell_metadata[[feat_type]][[spat_unit]]
            cell_metadata[, eval(new_cluster_name) := NULL]
            gobject@cell_metadata[[feat_type]][[spat_unit]] <- cell_metadata
        }

        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = metadata[
                , c("cell_ID", new_cluster_name),
                with = FALSE
            ],
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )


        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_merge_cluster")
        return(gobject)
    } else {
        return(list(mergevector = finalvec, metadata = metadata))
    }
}





#' @title Split dendrogram in two
#' @name .split_dendrogram_in_two
#' @description Merge selected clusters based on pairwise correlation scores
#' and size of cluster.
#' @param dend dendrogram object
#' @returns list of two dendrograms and height of node
#' @keywords internal
.split_dendrogram_in_two <- function(dend) {
    top_height <- attributes(dend)$height
    divided_leaves_labels <- dendextend::cut_lower_fun(dend, h = top_height)

    # this works for both numericala nd character leave names
    all_leaves <- dendextend::get_leaves_attr(dend = dend, attribute = "label")
    selected_labels_ind_1 <- all_leaves %in% divided_leaves_labels[[1]]
    selected_labels_ind_2 <- all_leaves %in% divided_leaves_labels[[2]]
    numerical_leaves <- unlist(dend)
    names(numerical_leaves) <- all_leaves

    dend_1 <- dendextend::find_dendrogram(
        dend = dend,
        selected_labels = names(numerical_leaves[selected_labels_ind_1])
    )
    dend_2 <- dendextend::find_dendrogram(
        dend = dend,
        selected_labels = names(numerical_leaves[selected_labels_ind_2])
    )

    return(list(theight = top_height, dend1 = dend_1, dend2 = dend_2))
}



#' @title Node clusters
#' @name .node_clusters
#' @description Merge selected clusters based on pairwise correlation scores
#' and size of cluster.
#' @param hclus_obj hclus object
#' @param verbose be verbose
#' @returns list of splitted dendrogram nodes from high to low node height
#' @keywords internal
.node_clusters <- function(hclus_obj, verbose = TRUE) {
    heights <- sort(hclus_obj[["height"]], decreasing = TRUE)
    mydend <- stats::as.dendrogram(hclus_obj)


    result_list <- list()
    j <- 1

    dend_list <- list()
    i <- 1
    dend_list[[i]] <- mydend

    ## create split at each height ##
    for (n_height in heights) {
        if (verbose == TRUE) cat("height ", n_height, "\n")

        # only use dendrogram objects
        ind <- lapply(dend_list, FUN = function(x) class(x) == "dendrogram")
        dend_list <- dend_list[unlist(ind)]

        # check which heights are available
        available_h <- as.numeric(unlist(lapply(
            dend_list,
            FUN = function(x) attributes(x)$height
        )))

        # get dendrogram associated with height and split in two
        select_dend_ind <- which.min(abs(available_h - n_height))
        select_dend <- dend_list[[select_dend_ind]]
        tempres <- .split_dendrogram_in_two(dend = select_dend)

        # find leave labels
        toph <- tempres[[1]]
        first_group <- dendextend::get_leaves_attr(
            tempres[[2]],
            attribute = "label"
        )
        second_group <- dendextend::get_leaves_attr(
            tempres[[3]],
            attribute = "label"
        )

        result_list[[j]] <- list(
            "height" = toph,
            "first" = first_group,
            "sec" = second_group
        )
        j <- j + 1



        ## add dendrograms to list
        ind <- lapply(tempres, FUN = function(x) class(x) == "dendrogram")
        tempres_dend <- tempres[unlist(ind)]

        dend_list[[i + 1]] <- tempres_dend[[1]]
        dend_list[[i + 2]] <- tempres_dend[[2]]

        i <- i + 2
    }

    return(list(dend_list, result_list))
}



#' @title getDendrogramSplits
#' @name getDendrogramSplits
#' @description Split dendrogram at each node and keep the leave (label)
#' information.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @param distance distance method to use for hierarchical clustering
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param show_dend show dendrogram
#' @param verbose be verbose
#' @returns data.table object
#' @details Creates a data.table with three columns and each row represents a
#' node in the dendrogram. For each node the height of the node is given
#' together with the two subdendrograms. This information can be used to
#' determine in a hierarchical manner differentially expressed marker genes at
#' each node.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' getDendrogramSplits(g, cluster_column = "leiden_clus")
#' @export
getDendrogramSplits <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    cluster_column,
    cor = c("pearson", "spearman"),
    distance = "ward.D",
    h = NULL,
    h_color = "red",
    show_dend = TRUE,
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

    # package check for dendextend
    package_check(pkg_name = "dendextend", repository = "CRAN")

    # data.table variables
    nodeID <- NULL

    cor <- match.arg(cor, c("pearson", "spearman"))
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )

    # create average expression matrix per cluster
    metatable <- calculateMetaTable(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        metadata_cols = cluster_column
    )
    dcast_metatable <- data.table::dcast.data.table(
        metatable,
        formula = variable ~ uniq_ID, value.var = "value"
    )
    testmatrix <- dt_to_matrix(x = dcast_metatable)

    # correlation
    cormatrix <- cor_flex(x = testmatrix, method = cor)
    cordist <- stats::as.dist(1 - cormatrix, diag = TRUE, upper = TRUE)
    corclus <- stats::hclust(d = cordist, method = distance)

    cordend <- stats::as.dendrogram(object = corclus)


    if (show_dend == TRUE) {
        # plot dendrogram
        plot(cordend)

        # add horizontal lines
        if (!is.null(h)) {
            abline(h = h, col = h_color)
        }
    }


    splitList <- .node_clusters(hclus_obj = corclus, verbose = verbose)

    splitDT <- data.table::as.data.table(t_flex(
        data.table::as.data.table(splitList[[2]])
    ))
    colnames(splitDT) <- c("node_h", "tree_1", "tree_2")
    splitDT[, nodeID := paste0("node_", seq_len(.N))]

    return(splitDT)
}






# projection ####


# * labelTransfer ####

#' @name labelTransfer
#' @title Transfer labels/annotations between sets of data via similarity
#' voting
#' @description
#' When two sets of data share an embedding space, transfer the labels from
#' one of the sets to the other based on KNN similarity voting in that space.
#' @param x target object
#' @param y source object
#' @param source_cell_ids cell/spatial IDs with the source labels to transfer
#' @param target_cell_ids cell/spatial IDs to transfer the labels to.
#' IDs from `source_cell_ids` are always included as well.
#' @param labels metadata column in source with labels to transfer
#' @param k number of k-neighbors to train a KNN classifier
#' @param name metadata column in target to apply the full set of labels to
#' @param prob output knn probabilities together with label predictions
#' @param reduction reduction on cells or features (default = "cells")
#' @param reduction_method shared reduction method (default = "pca" space)
#' @param reduction_name name of shared reduction space (default name = "pca")
#' @param dimensions_to_use dimensions to use in shared reduction space
#' (default = 1:10)
#' @returns object `x` with new transferred labels added to metadata
#' @inheritDotParams FNN::knn -train -test -cl -k -prob
#' @details
#' This function trains a KNN classifier with [FNN::knn()].
#' The training data is from object `y` or `source_cell_ids` subset in `x` and
#' uses existing annotations within the cell metadata.
#' Cells without annotation/labels from `x` or `target_cell_ids` subset in `x`
#' will receive predicted labels (and optional probabilities when
#' `prob = TRUE`).
#'
#' **IMPORTANT** This projection assumes that you're using the same dimension
#' reduction space (e.g. PCA) and number of dimensions (e.g. first 10 PCs) to
#' train the KNN classifier as you used to create the initial
#' annotations/labels in the source Giotto object.
#'
#' This function can allow you to work with very big data as you can predict
#' cell labels on a smaller & subsetted Giotto object and then project the cell
#' labels to the remaining cells in the target Giotto object. It can also be
#' used to transfer labels from one set of annotated data to another dataset
#' based on expression similarity after joining and integrating.
#'
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' id_subset <- sample(spatIDs(g), 300)
#' n_pred <- nrow(pDataDT(g)) - 300
#'
#' # transfer labels from one object to another ###################
#' g_small <- g[, id_subset]
#' # additional steps to get labels to transfer on smaller object...
#' g <- labelTransfer(g, g_small, labels = "leiden_clus")
#'
#' # transfer labels between subsets of a single object ###########
#' g <- labelTransfer(g,
#'     label = "leiden_clus", source_cell_ids = id_subset, name = "knn_leiden2"
#' )
#' @md
NULL

setGeneric("labelTransfer",
           function(x, y, ...) standardGeneric("labelTransfer"))

#' @rdname labelTransfer
#' @export
setMethod("labelTransfer", signature(x = "giotto", y = "giotto"), function(x, y,
    spat_unit = NULL,
    feat_type = NULL,
    labels,
    k = 10,
    name = paste0("trnsfr_", labels),
    prob = TRUE,
    reduction = "cells",
    reduction_method = "pca",
    reduction_name = "pca",
    dimensions_to_use = 1:10,
    return_gobject = TRUE,
    ...) {
    # NSE vars
    temp_name <- cell_ID <- temp_name_prob <- NULL

    package_check(pkg_name = "FNN", repository = "CRAN")
    spat_unit <- set_default_spat_unit(x, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(x,
        spat_unit = spat_unit, feat_type = feat_type
    )

    # get data
    cx_src <- getCellMetadata(y,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )
    cx_tgt <- getCellMetadata(x,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )
    dim_coord <- getDimReduction(x,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = reduction,
        reduction_method = reduction_method,
        name = reduction_name,
        output = "matrix"
    )

    # source annotation vector #
    # names : cell_ID
    # values: label
    source_annot_vec <- cx_src[[labels]]
    names(source_annot_vec) <- cx_src[["cell_ID"]]

    # create the matrix from the target object that you want to use for the
    # kNN classifier
    # the matrix should be the same for the source and target objects
    # (e.g. same PCA space)
    dimensions_to_use <- dimensions_to_use[
        # ensure dims to use exist
        dimensions_to_use %in% seq_len(ncol(dim_coord))
    ]
    matrix_to_use <- dim_coord[, dimensions_to_use]

    ## create the training and testset from the matrix

    # the training set is those spatial IDs that are in the source
    # (w/ labels) AND target giotto object
    in_common <- rownames(matrix_to_use) %in% names(source_annot_vec)
    train <- matrix_to_use[in_common, ]
    train <- train[match(names(source_annot_vec), rownames(train)), ]

    # the test set are the remaining cell_IDs that need a label
    test <- matrix_to_use[!in_common, ]

    # make prediction
    knnprediction <- FNN::knn(
        train = train,
        test = test,
        cl = source_annot_vec,
        k = k,
        prob = prob,
        ...
    )

    # get prediction results
    knnprediction_vec <- as.vector(knnprediction)
    names(knnprediction_vec) <- rownames(test)

    # add probability information
    if (isTRUE(prob)) {
        probs <- attr(knnprediction, "prob")
        names(probs) <- rownames(test)
    }

    # create annotation vector for all cell IDs (from source and predicted)
    all_vec <- c(source_annot_vec, knnprediction_vec)
    cx_tgt[, temp_name := all_vec[cell_ID]]

    if (isTRUE(prob)) {
        cx_tgt[, temp_name_prob := probs[cell_ID]]
        cx_tgt <- cx_tgt[, .(cell_ID, temp_name, temp_name_prob)]
        cx_tgt[, temp_name_prob := ifelse(
            is.na(temp_name_prob), 1, temp_name_prob
        )]

        data.table::setnames(cx_tgt,
            old = c("temp_name", "temp_name_prob"),
            new = c(name, paste0(name, "_prob"))
        )
    } else {
        cx_tgt <- cx_tgt[, .(cell_ID, temp_name)]
        data.table::setnames(cx_tgt, old = "temp_name", new = name)
    }


    if (return_gobject) {
        x <- addCellMetadata(x,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = cx_tgt,
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )
        return(x)
    } else {
        return(cx_tgt)
    }
})

#' @rdname labelTransfer
#' @export
setMethod("labelTransfer", signature(x = "giotto", y = "missing"), function(x,
    spat_unit = NULL,
    feat_type = NULL,
    source_cell_ids,
    target_cell_ids,
    labels,
    k = 10,
    name = paste0("trnsfr_", labels),
    prob = TRUE,
    reduction = "cells",
    reduction_method = "pca",
    reduction_name = "pca",
    dimensions_to_use = 1:10,
    return_gobject = TRUE,
    ...) {
    # NSE vars
    temp_name <- cell_ID <- temp_name_prob <- NULL

    package_check(pkg_name = "FNN", repository = "CRAN")
    spat_unit <- set_default_spat_unit(x, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(x,
        spat_unit = spat_unit, feat_type = feat_type
    )

    # get data
    cx <- getCellMetadata(x,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )
    dim_coord <- getDimReduction(x,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = reduction,
        reduction_method = reduction_method,
        name = reduction_name,
        output = "matrix"
    )

    # source annotation vector #
    # names : cell_ID
    # values: label
    source_annot_vec <- cx[[labels]]
    names(source_annot_vec) <- cx[["cell_ID"]]
    source_annot_vec <- source_annot_vec[source_cell_ids]

    # target cell IDs (if not provided) are everything not in the source cell
    # IDs
    if (missing(target_cell_ids)) {
        sids <- cx[["cell_ID"]]
        target_cell_ids <- sids[!sids %in% source_cell_ids]
    }

    # create the matrix from the target object that you want to use for the
    # kNN classifier
    # the matrix should be the same for the source and target objects
    # (e.g. same PCA space)
    dimensions_to_use <- dimensions_to_use[
        # ensure dims to use exist
        dimensions_to_use %in% seq_len(ncol(dim_coord))
    ]
    matrix_to_use <- dim_coord[, dimensions_to_use]

    ## create the training and testset from the matrix

    # the training set is those spatial IDs that are in the source
    # (w/ labels) AND target giotto object
    train <- matrix_to_use[source_cell_ids, ]
    train <- train[match(names(source_annot_vec), rownames(train)), ]

    # the test set are the remaining cell_IDs that need a label
    test <- matrix_to_use[target_cell_ids, ]

    # make prediction
    knnprediction <- FNN::knn(
        train = train,
        test = test,
        cl = source_annot_vec,
        k = k,
        prob = prob,
        ...
    )

    # get prediction results
    knnprediction_vec <- as.vector(knnprediction)
    names(knnprediction_vec) <- rownames(test)

    # add probability information
    if (isTRUE(prob)) {
        probs <- attr(knnprediction, "prob")
        names(probs) <- rownames(test)
    }

    # create annotation vector for all cell IDs (from source and predicted)
    all_vec <- c(source_annot_vec, knnprediction_vec)
    cx[, temp_name := all_vec[cell_ID]]

    if (isTRUE(prob)) {
        cx[, temp_name_prob := probs[cell_ID]]
        cx <- cx[, .(cell_ID, temp_name, temp_name_prob)]
        cx[, temp_name_prob := ifelse(
            is.na(temp_name_prob), 1, temp_name_prob
        )]

        data.table::setnames(cx,
            old = c("temp_name", "temp_name_prob"),
            new = c(name, paste0(name, "_prob"))
        )
    } else {
        cx <- cx[, .(cell_ID, temp_name)]
        data.table::setnames(cx, old = "temp_name", new = name)
    }


    if (return_gobject) {
        x <- addCellMetadata(x,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = cx,
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )
        return(x)
    } else {
        return(cx)
    }
})





#' @title Projection of cluster labels
#' @name doClusterProjection
#' @description Use a fast KNN classifier to predict labels from a smaller
#' giotto object
#' @param target_gobject target giotto object
#' @param target_cluster_label_name name for predicted clusters
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param source_gobject source giotto object with annotation data
#' @param source_cluster_labels annotation/labels to use to train KNN classifier
#' @param reduction reduction on cells or features (default = cells)
#' @param reduction_method shared reduction method (default = pca space)
#' @param reduction_name name of shared reduction space (default name = 'pca')
#' @param dimensions_to_use dimensions to use in shared reduction space
#' (default = 1:10)
#' @param knn_k number of k-neighbors to train a KNN classifier
#' @param prob output probabilities together with label predictions
#' @param algorithm nearest neighbor search algorithm
#' @param return_gobject return giotto object
#' @returns giotto object (default) or data.table with cell metadata
#'
#' @details Function to train a KNN with \code{\link[FNN]{knn}}. The training
#' data is obtained from the source giotto object (source_gobject) using
#' existing annotations within the cell metadata. Cells without
#' annotation/labels from the target giotto object (target_gobject) will
#' receive predicted labels (and optional probabilities with prob = TRUE).
#'
#' **IMPORTANT** This projection assumes that you're using the same dimension
#' reduction space (e.g. PCA) and number of dimensions (e.g. first 10 PCs) to
#' train the KNN classifier as you used to create the initial
#' annotations/labels in the source Giotto object.
#'
#' Altogether this is a convenience function that allow you to work with very
#' big data as you can predict cell labels on a smaller & subsetted Giotto
#' object and then project the cell labels to the remaining cells in the target
#' Giotto object.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- pDataDT(g)
#' g_small <- subsetGiotto(g, cell_ids = sample(x$cell_ID, 300))
#' doClusterProjection(
#'     target_gobject = g, source_gobject = g_small,
#'     source_cluster_labels = "leiden_clus"
#' )
#' @export
doClusterProjection <- function(target_gobject,
    target_cluster_label_name = "knn_labels",
    spat_unit = NULL,
    feat_type = NULL,
    source_gobject,
    source_cluster_labels = NULL,
    reduction = "cells",
    reduction_method = "pca",
    reduction_name = "pca",
    dimensions_to_use = 1:10,
    knn_k = 10,
    prob = FALSE,
    algorithm = c(
        "kd_tree",
        "cover_tree", "brute"
    ),
    return_gobject = TRUE) {
    deprecate_warn(
        when = "4.1.2",
        what = "doClusterProjection()",
        with = "labelTransfer()"
    )

    # NSE vars
    cell_ID <- temp_name_prob <- NULL

    package_check(pkg_name = "FNN", repository = "CRAN")

    spat_unit <- set_default_spat_unit(
        gobject = target_gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = target_gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # identify clusters from source object and create annotation vector
    cell_meta_source <- getCellMetadata(
        gobject = source_gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )
    source_annot_vec <- cell_meta_source[[source_cluster_labels]]
    names(source_annot_vec) <- cell_meta_source[["cell_ID"]]

    # create the matrix from the target object that you want to use for the
    # kNN classifier
    # the matrix should be the same for the source and target objects
    # (e.g. same PCA space)
    dim_obj <- getDimReduction(
        gobject = target_gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = reduction,
        reduction_method = reduction_method,
        name = reduction_name,
        output = "dimObj"
    )

    dim_coord <- dim_obj[]
    dimensions_to_use <- dimensions_to_use[
        dimensions_to_use %in% seq_len(ncol(dim_coord))
    ]
    matrix_to_use <- dim_coord[, dimensions_to_use]

    ## create the training and testset from the matrix

    # the training set is the set of cell IDs that are in both the source
    # (w/ labels)
    # and target giotto object
    train <- matrix_to_use[
        rownames(matrix_to_use) %in% names(source_annot_vec),
    ]
    train <- train[match(names(source_annot_vec), rownames(train)), ]

    # the test set are the remaining cell_IDs that need a label
    test <- matrix_to_use[
        !rownames(matrix_to_use) %in% names(source_annot_vec),
    ]
    cl <- source_annot_vec

    # make prediction
    knnprediction <- FNN::knn(
        train = train, test = test,
        cl = cl, k = knn_k, prob = prob,
        algorithm = algorithm
    )

    # get prediction results
    knnprediction_vec <- as.vector(knnprediction)
    names(knnprediction_vec) <- rownames(test)

    # add probability information
    if (isTRUE(prob)) {
        probs <- attr(knnprediction, "prob")
        names(probs) <- rownames(test)
    }



    # create annotation vector for all cell IDs (from source and predicted)
    all_vec <- c(source_annot_vec, knnprediction_vec)

    cell_meta_target <- getCellMetadata(
        gobject = target_gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )
    # data.table variables
    temp_name <- NULL
    cell_meta_target[, temp_name := all_vec[cell_ID]]

    if (isTRUE(prob)) {
        cell_meta_target[, temp_name_prob := probs[cell_ID]]
        cell_meta_target <- cell_meta_target[
            , .(cell_ID, temp_name, temp_name_prob)
        ]
        cell_meta_target[, temp_name_prob := ifelse(
            is.na(temp_name_prob), 1, temp_name_prob
        )]

        data.table::setnames(cell_meta_target,
            old = c("temp_name", "temp_name_prob"),
            new = c(
                target_cluster_label_name,
                paste0(target_cluster_label_name, "_prob")
            )
        )
    } else {
        cell_meta_target <- cell_meta_target[, .(cell_ID, temp_name)]
        data.table::setnames(cell_meta_target,
            old = "temp_name",
            new = target_cluster_label_name
        )
    }


    if (return_gobject) {
        if (isTRUE(prob)) {
            prob_label <- paste0(target_cluster_label_name, "_prob")

            target_gobject <- addCellMetadata(
                gobject = target_gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                new_metadata = cell_meta_target[
                    , c("cell_ID", target_cluster_label_name, prob_label),
                    with = FALSE
                ],
                by_column = TRUE,
                column_cell_ID = "cell_ID"
            )
        } else {
            target_gobject <- addCellMetadata(
                gobject = target_gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                new_metadata = cell_meta_target[
                    , c("cell_ID", target_cluster_label_name),
                    with = FALSE
                ],
                by_column = TRUE,
                column_cell_ID = "cell_ID"
            )
        }
    } else {
        return(cell_meta_target)
    }
}

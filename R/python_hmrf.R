#' @title doHMRF
#' @name doHMRF
#' @description Run HMRF
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param spatial_network_name name of spatial network to use for HMRF
#' @param spat_loc_name name of spatial locations
#' @param spatial_genes spatial genes to use for HMRF
#' @param spatial_dimensions select spatial dimensions to use, default is all
#' possible dimensions
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name name of HMRF run
#' @param k  number of HMRF domains
#' @param seed seed to fix random number generator
#' (for creating initialization of HMRF) (-1 if no fixing)
#' @param betas betas to test for. three numbers: start_beta, beta_increment,
#' num_betas e.g. c(0, 2.0, 50)
#' @param tolerance tolerance
#' @param zscore zscore
#' @param numinit number of initializations
#' @param python_path python path to use
#' @param output_folder output folder to save results
#' @param overwrite_output overwrite output folder
#' @returns Creates a directory with results that can be viewed with
#' viewHMRFresults
#' @details Description of HMRF parameters ...
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' spat_genes <- binSpect(g)
#'
#' doHMRF(g,
#'     spatial_genes = spat_genes[seq_len(10)]$feats,
#'     output_folder = tempdir()
#' )
#' @export
doHMRF <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    spatial_network_name = "Delaunay_network",
    spat_loc_name = "raw",
    spatial_genes = NULL,
    spatial_dimensions = c("sdimx", "sdimy", "sdimz"),
    dim_reduction_to_use = NULL,
    dim_reduction_name = "pca",
    dimensions_to_use = 1:10,
    seed = 100,
    name = "test",
    k = 10,
    betas = c(0, 2, 50),
    tolerance = 1e-10,
    zscore = c("none", "rowcol", "colrow"),
    numinit = 100,
    python_path = NULL,
    output_folder = NULL,
    overwrite_output = TRUE) {
    package_check("smfishHmrf", repository = "pip")

    # data.table set global variable
    to <- from <- NULL

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

    ## check or make paths
    # python path
    if (is.null(python_path)) {
        python_path <- readGiottoInstructions(gobject, param = "python_path")
    }

    ## reader.py and get_result.py paths
    reader_path <- system.file("python", "reader2.py", package = "Giotto")

    ## output folder
    # no folder path specified
    if (is.null(output_folder)) {
        output_folder <- paste0(getwd(), "/", "HMRF_output")
        if (!file.exists(output_folder)) {
            dir.create(
                path = paste0(getwd(), "/", "HMRF_output"), recursive = TRUE
            )
        }
    }
    # folder path specified
    else if (!is.null(output_folder)) {
        output_folder <- path.expand(output_folder)
        if (!file.exists(output_folder)) {
            dir.create(path = output_folder, recursive = TRUE)
        }
    }


    ## first write necessary txt files to output folder ##
    # cell location / spatial network / expression data and selected spatial
    # genes

    ## 1. expression values
    if (!is.null(dim_reduction_to_use)) {
        expr_values <- getDimReduction(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "data.table"
        )
        expr_values <- expr_values[, dimensions_to_use]
        expr_values <- t_flex(expr_values)
    } else {
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
    }

    if (!"matrix" %in% class(expr_values)) {
        warning("\n this matrix will be converted to a dense and memory intensive
                base matrix ...\n")
        expr_values <- as.matrix(expr_values)
    }


    expression_file <- paste0(output_folder, "/", "expression_matrix.txt")

    # overwrite if exists
    if (file.exists(expression_file) & overwrite_output == TRUE) {
        message("\n expression_matrix.txt already exists at this location,
        will be overwritten")
        data.table::fwrite(
            data.table::as.data.table(expr_values, keep.rownames = "gene"),
            file = expression_file, quote = FALSE, col.names = TRUE,
            row.names = FALSE, sep = " "
        )
    } else if (file.exists(expression_file) & overwrite_output == FALSE) {
        message("\n expression_matrix.txt already exists at this location,
        will be used again")
    } else {
        data.table::fwrite(
            data.table::as.data.table(expr_values, keep.rownames = "gene"),
            file = expression_file, quote = FALSE, col.names = TRUE,
            row.names = FALSE, sep = " "
        )
    }






    ## 2. spatial genes
    if (!is.null(dim_reduction_to_use)) {
        dimred_rownames <- rownames(expr_values)
        spatial_genes_detected <- dimred_rownames[dimensions_to_use]
        spatial_genes_detected <- spatial_genes_detected[
            !is.na(spatial_genes_detected)
        ]
    } else {
        if (is.null(spatial_genes)) {
            stop("\n you need to provide a vector of spatial genes (~500)")
        }
        spatial_genes_detected <- spatial_genes[
            spatial_genes %in% rownames(expr_values)
        ]
    }
    spatial_genes_file <- paste0(output_folder, "/", "spatial_genes.txt")

    # overwrite if exists
    if (file.exists(spatial_genes_file) & overwrite_output == TRUE) {
        message("\n spatial_genes.txt already exists at this location, will be
                overwritten")
        write.table(spatial_genes_detected,
            file = spatial_genes_file,
            quote = FALSE, col.names = FALSE, row.names = FALSE
        )
    } else if (file.exists(spatial_genes_file) & overwrite_output == FALSE) {
        message("\n spatial_genes.txt already exists at this location, will be
                used again")
    } else {
        write.table(spatial_genes_detected,
            file = spatial_genes_file,
            quote = FALSE, col.names = FALSE, row.names = FALSE
        )
    }




    ## 3. spatial network
    spatial_network <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "networkDT"
    )
    spatial_network <- spatial_network[, .(to, from)]
    spatial_network_file <- paste0(output_folder, "/", "spatial_network.txt")

    if (file.exists(spatial_network_file) & overwrite_output == TRUE) {
        message("\n spatial_network.txt already exists at this location, will be
                overwritten")
        write.table(spatial_network,
            file = spatial_network_file,
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
        )
    } else if (file.exists(spatial_network_file) & overwrite_output == FALSE) {
        message("\n spatial_network.txt already exists at this location, will be
                used again")
    } else {
        write.table(spatial_network,
            file = spatial_network_file,
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
        )
    }




    ## 4. cell location
    spatial_location <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table",
        copy_obj = TRUE
    )

    # select spatial dimensions that are available #
    spatial_dimensions <- spatial_dimensions[
        spatial_dimensions %in% colnames(spatial_location)
    ]
    spatial_location <- spatial_location[
        , c(spatial_dimensions, "cell_ID"),
        with = FALSE
    ]
    spatial_location_file <- paste0(
        output_folder, "/", "spatial_cell_locations.txt"
    )

    if (file.exists(spatial_location_file) & overwrite_output == TRUE) {
        message("\n spatial_cell_locations.txt already exists at this location,
                will be overwritten")
        write.table(spatial_location,
            file = spatial_location_file,
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
        )
    } else if (file.exists(spatial_location_file)) {
        message("\n spatial_cell_locations.txt already exists at this location,
                will be used again")
    } else {
        write.table(spatial_location,
            file = spatial_location_file,
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
        )
    }




    # prepare input paths
    cell_location <- paste0(output_folder, "/", "spatial_cell_locations.txt")
    spatial_genes <- paste0(output_folder, "/", "spatial_genes.txt")
    spatial_network <- paste0(output_folder, "/", "spatial_network.txt")
    expression_data <- paste0(output_folder, "/", "expression_matrix.txt")

    # create output subfolder for HMRF
    output_data <- paste0(output_folder, "/", "result.spatial.zscore")
    if (!file.exists(output_data)) dir.create(output_data)

    # encapsulate to avoid path problems
    # python code also needs to be updated internally
    cell_location <- paste0('"', cell_location, '"')
    spatial_genes <- paste0('"', spatial_genes, '"')
    spatial_network <- paste0('"', spatial_network, '"')
    expression_data <- paste0('"', expression_data, '"')
    output_data <- paste0('"', output_data, '"')

    # process other params
    zscore <- match.arg(zscore, c("none", "rowcol", "colrow"))
    betas_param <- c("-b", betas)
    betas_final <- paste(betas_param, collapse = " ")

    ## reader part ##
    reader_command <- paste0(
        python_path, " ", reader_path,
        " -l ", cell_location,
        " -g ", spatial_genes,
        " -n ", spatial_network,
        " -e ", expression_data,
        " -o ", output_data,
        " -a ", name,
        " -k ", k,
        " ", betas_final,
        " -t ", tolerance,
        " -z ", zscore,
        " -s ", seed,
        " -i ", numinit
    )

    print(reader_command)
    system(command = reader_command)


    # store parameter results in HMRF S3 object
    HMRFObj <- list(
        name = name,
        feat_type = feat_type,
        output_data = output_data,
        k = k,
        betas = betas,
        python_path = python_path
    )

    class(HMRFObj) <- append(class(HMRFObj), "HMRFoutput")


    return(HMRFObj)
}



#' @title loadHMRF
#' @name loadHMRF
#' @description load previous HMRF
#' @param name_used name of HMRF that was run
#' @param k_used  number of HMRF domains that was tested
#' @param betas_used betas that were tested
#' @param python_path_used python path that was used
#' @param output_folder_used output folder that was used
#' @returns reloads a previous ran HMRF from doHRMF
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- tempdir()
#' doHMRF(g,
#'     spatial_genes = c("Gna12", "Ccnd2"), output_folder = x,
#'     betas = c(0, 2, 50)
#' )
#'
#' loadHMRF(
#'     output_folder_used = x, betas_used = c(0, 2, 50),
#'     python_path_used = NULL
#' )
#'
#' @export
loadHMRF <- function(
        name_used = "test",
        output_folder_used,
        k_used = 10,
        betas_used,
        python_path_used) {
    output_data <- paste0(output_folder_used, "/", "result.spatial.zscore")
    if (!file.exists(output_data)) {
        stop("\n doHMRF was not run in this output directory")
    }

    # check if it indeed exists

    HMRFObj <- list(
        name = name_used,
        output_data = output_data,
        k = k_used,
        betas = betas_used,
        python_path = python_path_used
    )

    class(HMRFObj) <- append(class(HMRFObj), "HMRFoutput")


    return(HMRFObj)
}



#' @title viewHMRFresults
#' @name viewHMRFresults
#' @description View results from doHMRF.
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF
#' @param k number of HMRF domains
#' @param betas_to_view results from different betas that you want to view
#' @param third_dim 3D data (boolean)
#' @param \dots additional paramters (see details)
#' @returns spatial plots with HMRF domains
#' @seealso \code{\link{spatPlot2D}} and \code{\link{spatPlot3D}}
#' @export
viewHMRFresults <- function(
        gobject,
        HMRFoutput,
        k = NULL,
        betas_to_view = NULL,
        third_dim = FALSE,
        ...) {
    if (!"HMRFoutput" %in% class(HMRFoutput)) {
        stop("\n HMRFoutput needs to be output from doHMRFextend")
    }

    ## reader.py and get_result.py paths
    # TODO: part of the package
    get_result_path <- system.file(
        "python", "get_result2.py",
        package = "Giotto"
    )

    # paths and name
    name <- HMRFoutput$name
    output_data <- HMRFoutput$output_data
    python_path <- HMRFoutput$python_path

    # k-values
    if (is.null(k)) {
        stop("you need to select a k that was used with doHMRFextend")
    }
    k <- HMRFoutput$k

    # betas
    betas <- HMRFoutput$betas
    possible_betas <- seq(betas[1],
        to = betas[1] + (betas[2] * (betas[3] - 1)),
        by = betas[2]
    )

    betas_to_view_detected <- betas_to_view[betas_to_view %in% possible_betas]

    # plot betas
    for (b in betas_to_view_detected) {
        ## get results part ##
        result_command <- paste0(
            python_path, " ", get_result_path,
            " -r ", output_data,
            " -a ", name,
            " -k ", k,
            " -b ", b
        )

        print(result_command)

        output <- system(command = result_command, intern = TRUE)


        title_name <- paste0("k = ", k, " b = ", b)

        spatPlot2D(
            gobject = gobject,
            cell_color = output,
            show_plot = TRUE,
            title = title_name,
            ...
        )

        if (third_dim == TRUE) {
            spatPlot3D(
                gobject = gobject,
                cell_color = output,
                show_plot = TRUE,
                ...
            )
        }
    }
}



#' @title writeHMRFresults
#' @name writeHMRFresults
#' @description write results from doHMRF to a data.table.
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF
#' @param k k to write results for
#' @param betas_to_view results from different betas that you want to view
#' @param print_command see the python command
#' @returns data.table with HMRF results for each b and the selected k
#' @export
writeHMRFresults <- function(
        gobject,
        HMRFoutput,
        k = NULL,
        betas_to_view = NULL,
        print_command = FALSE) {
    if (!"HMRFoutput" %in% class(HMRFoutput)) {
        stop("\n HMRFoutput needs to be output from doHMRFextend")
    }

    ## reader.py and get_result.py paths
    # TODO: part of the package
    get_result_path <- system.file(
        "python", "get_result2.py",
        package = "Giotto"
    )

    # paths and name
    name <- HMRFoutput$name
    output_data <- HMRFoutput$output_data
    python_path <- HMRFoutput$python_path

    # k-values
    if (is.null(k)) {
        stop("\n you need to select a k that was used with doHMRFextend")
    }
    k <- HMRFoutput$k

    # betas
    betas <- HMRFoutput$betas
    possible_betas <- seq(betas[1],
        to = betas[1] + (betas[2] * (betas[3] - 1)),
        by = betas[2]
    )

    betas_to_view_detected <- betas_to_view[betas_to_view %in% possible_betas]

    result_list <- list()

    # plot betas
    for (i in seq_along(betas_to_view_detected)) {
        b <- betas_to_view_detected[i]

        ## get results part ##
        result_command <- paste0(
            python_path, " ", get_result_path,
            " -r ", output_data,
            " -a ", name,
            " -k ", k,
            " -b ", b
        )

        if (print_command == TRUE) {
            print(result_command)
        }

        output <- system(command = result_command, intern = TRUE)
        title_name <- paste0("k.", k, ".b.", b)
        result_list[[title_name]] <- output
    }

    result_DT <- data.table::as.data.table(do.call("cbind", result_list))
    result_DT <- cbind(data.table::data.table(
        "cell_ID" = gobject@cell_ID
    ), result_DT)
    return(result_DT)
}




#' @title addHMRF
#' @name addHMRF
#' @description Add selected results from doHMRF to the giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param HMRFoutput HMRF output from doHMRF()
#' @param k number of domains
#' @param betas_to_add results from different betas that you want to add
#' @param hmrf_name specify a custom name
#' @returns giotto object
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' spat_genes <- binSpect(g)
#'
#' output_folder <- file.path(tempdir(), "HMRF")
#' if (!file.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
#'
#' out <- doHMRF(
#'     g,
#'     spatial_genes = spat_genes[seq_len(20)]$feats,
#'     expression_values = "scaled",
#'     spatial_network_name = "Delaunay_network",
#'     k = 6, betas = c(0, 10, 5),
#'     output_folder = output_folder
#' )
#'
#' g <- addHMRF(
#'     gobject = g,
#'     HMRFoutput = out,
#'     k = 6,
#'     betas_to_add = 20,
#'     hmrf_name = "HMRF"
#' )
#'
#' spatPlot(
#'     gobject = g, cell_color = "HMRF_k6_b.20",
#' )
#' @export
addHMRF <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        HMRFoutput,
        k = NULL,
        betas_to_add = NULL,
        hmrf_name = NULL) {
    if (!"HMRFoutput" %in% class(HMRFoutput)) {
        stop("\n HMRFoutput needs to be output from doHMRFextend")
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

    # get feat_type
    feat_type <- HMRFoutput$feat_type

    ## reader.py and get_result.py paths
    # TODO: part of the package
    get_result_path <- system.file(
        "python", "get_result2.py",
        package = "Giotto"
    )

    # paths and name
    name <- HMRFoutput$name
    output_data <- HMRFoutput$output_data
    python_path <- HMRFoutput$python_path

    # k-values
    if (is.null(k)) {
        stop("\n you need to select a k that was used with doHMRFextend")
    }
    k <- HMRFoutput$k

    # betas
    betas <- HMRFoutput$betas
    possible_betas <- seq(betas[1],
        to = betas[1] + (betas[2] * (betas[3] - 1)),
        by = betas[2]
    )

    betas_to_add_detected <- betas_to_add[betas_to_add %in% possible_betas]


    # get cell metadata for object
    cell_metadata <- pDataDT(gobject, feat_type = feat_type)


    # plot betas
    for (b in betas_to_add_detected) {
        ## get results part ##
        result_command <- paste0(
            python_path, " ", get_result_path,
            " -r ", output_data,
            " -a ", name,
            " -k ", k,
            " -b ", b
        )
        print(result_command)
        output <- system(command = result_command, intern = TRUE)

        # get cell_IDs
        cids <- gsub(basename(output_data), "", output_data) %>%
            gsub(pattern = "\"", replacement = "") %>%
            list.files(
                pattern = "expression_matrix",
                full.names = TRUE
            ) %>%
            data.table::fread(nrows = 0L, header = TRUE) %>%
            colnames()
        cids <- cids[-1] # gene colname is also included

        # create unique name
        annot_DT <- data.table::data.table(
            cell_ID = cids,
            temp_name = output
        )

        if (!is.null(hmrf_name)) {
            annot_name <- paste0(hmrf_name, "_k", k, "_b.", b)
            setnames(annot_DT, old = "temp_name", new = annot_name)
        } else {
            annot_name <- paste0("hmrf_k.", k, "_b.", b)
            data.table::setnames(annot_DT, old = "temp_name", new = annot_name)
        }


        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            column_cell_ID = "cell_ID",
            new_metadata = annot_DT,
            by_column = TRUE
        )
    }

    return(gobject)
}





#' @title viewHMRFresults2D
#' @name viewHMRFresults2D
#' @description View results from doHMRF.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param HMRFoutput HMRF output from doHMRF
#' @param k number of HMRF domains
#' @param betas_to_view results from different betas that you want to view
#' @param \dots additional parameters to spatPlot2D()
#' @returns spatial plots with HMRF domains
#' @seealso \code{\link{spatPlot2D}}
#' @export
viewHMRFresults2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        HMRFoutput,
        k = NULL,
        betas_to_view = NULL,
        ...) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    if (!"HMRFoutput" %in% class(HMRFoutput)) {
        stop("\n HMRFoutput needs to be output from doHMRFextend")
    }

    ## reader.py and get_result.py paths
    # TODO: part of the package
    get_result_path <- system.file(
        "python", "get_result2.py",
        package = "Giotto"
    )

    # paths and name
    name <- HMRFoutput$name
    output_data <- HMRFoutput$output_data
    python_path <- HMRFoutput$python_path

    # k-values
    if (is.null(k)) {
        stop("\n please select a valid k used with doHMRFextend")
    }
    k <- HMRFoutput$k

    # betas
    betas <- HMRFoutput$betas
    possible_betas <- seq(betas[1],
        to = betas[1] + (betas[2] * (betas[3] - 1)),
        by = betas[2]
    )

    betas_to_view_detected <- betas_to_view[betas_to_view %in% possible_betas]

    # plot betas
    for (b in betas_to_view_detected) {
        ## get results part ##
        result_command <- paste0(
            python_path, " ", get_result_path,
            " -r ", output_data,
            " -a ", name,
            " -k ", k,
            " -b ", b
        )

        print(result_command)
        output <- system(command = result_command, intern = TRUE)

        # get cell_IDs
        cids <- gsub(basename(output_data), "", output_data) %>%
            gsub(pattern = "\"", replacement = "") %>%
            list.files(
                pattern = "expression_matrix",
                full.names = TRUE
            ) %>%
            data.table::fread(nrows = 0L, header = TRUE) %>%
            colnames()
        cids <- cids[-1] # gene colname is also included

        # create unique name
        annot_DT <- data.table::data.table(
            cell_ID = cids,
            temp_name = output
        )

        annot_name <- paste0("hmrf_k.", k, "_b.", b)
        data.table::setnames(annot_DT, old = "temp_name", new = annot_name)


        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            column_cell_ID = "cell_ID",
            new_metadata = annot_DT,
            by_column = TRUE
        )

        title_name <- paste0("k = ", k, " b = ", b)

        spatPlot2D(
            gobject = gobject,
            cell_color = annot_name,
            show_plot = TRUE,
            save_plot = FALSE,
            title = title_name,
            ...
        )
    }
}


#' @title viewHMRFresults3D
#' @name viewHMRFresults3D
#' @description View results from doHMRF.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param HMRFoutput HMRF output from doHMRF
#' @param k number of HMRF domains
#' @param betas_to_view results from different betas that you want to view
#' @param \dots additional parameters to spatPlot3D()
#' @returns spatial plots with HMRF domains
#' @seealso \code{\link{spatPlot3D}}
#' @export
viewHMRFresults3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        HMRFoutput,
        k = NULL,
        betas_to_view = NULL,
        ...) {
    if (!"HMRFoutput" %in% class(HMRFoutput)) {
        stop("\n HMRFoutput needs to be output from doHMRFextend")
    }

    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    ## reader.py and get_result.py paths
    # TODO: part of the package
    get_result_path <- system.file(
        "python", "get_result2.py",
        package = "Giotto"
    )

    # paths and name
    name <- HMRFoutput$name
    output_data <- HMRFoutput$output_data
    python_path <- HMRFoutput$python_path

    # k-values
    if (is.null(k)) {
        stop("\n please select a valid k used with doHMRFextend")
    }
    k <- HMRFoutput$k

    # betas
    betas <- HMRFoutput$betas
    possible_betas <- seq(betas[1],
        to = betas[1] + (betas[2] * (betas[3] - 1)),
        by = betas[2]
    )

    betas_to_view_detected <- betas_to_view[betas_to_view %in% possible_betas]

    # plot betas
    for (b in betas_to_view_detected) {
        ## get results part ##
        result_command <- paste0(
            python_path, " ", get_result_path,
            " -r ", output_data,
            " -a ", name,
            " -k ", k,
            " -b ", b
        )

        print(result_command)
        output <- system(command = result_command, intern = TRUE)

        # get cell_IDs
        cids <- gsub(basename(output_data), "", output_data) %>%
            gsub(pattern = "\"", replacement = "") %>%
            list.files(
                pattern = "expression_matrix",
                full.names = TRUE
            ) %>%
            data.table::fread(nrows = 0L, header = TRUE) %>%
            colnames()
        cids <- cids[-1] # gene colname is also included

        # create unique name
        annot_DT <- data.table::data.table(
            cell_ID = cids,
            temp_name = output
        )

        annot_name <- paste0("hmrf_k.", k, "_b.", b)
        data.table::setnames(annot_DT, old = "temp_name", new = annot_name)


        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            column_cell_ID = "cell_ID",
            new_metadata = annot_DT,
            by_column = TRUE
        )

        title_name <- paste0("k = ", k, " b = ", b)

        spatPlot3D(
            gobject = gobject,
            cell_color = annot_name,
            show_plot = TRUE,
            save_plot = FALSE,
            title = title_name,
            ...
        )
    }
}



################################################################################
### new HMRF functions V2

#' @title sampling_sp_genes
#' @name sampling_sp_genes
#' @description function to select a set of spatial genes
#' @param clust spatial gene clusters
#' @param sample_rate sampling rate, takes values equal or greater than 1
#' @param target target length of gene list
#' @param seed random seed
#' @returns list
#' @details
#' This function samples a subset of spatial genes among different clusters,
#' with size n = target.
#' Number of samples from each cluster denpends on the relative proportion of
#' each cluster.
#' Changing from equal size by setting sample_rate = 1 to with exact proportion
#' of each cluster by setting sample_rate = +Inf
#' @keywords internal
sampling_sp_genes <- function(
        clust,
        sample_rate = 2,
        target = 500,
        seed = 10) {
    tot <- 0
    num_cluster <- length(unique(clust))
    gene_list <- list()

    for (i in seq(1, num_cluster)) {
        gene_list[[i]] <- colnames(t(clust[which(clust == i)]))
    }
    for (i in seq(1, num_cluster)) {
        num_g <- length(gene_list[[i]])
        tot <- tot + num_g / (num_g^(1 / sample_rate))
    }
    factor <- target / tot
    num_sample <- c()
    genes <- c()
    for (i in seq(1, num_cluster)) {
        num_g <- length(gene_list[[i]])
        genes[i] <- num_g
        num_sample[i] <- round(num_g / (num_g^(1 / sample_rate)) * factor)
    }
    set.seed(seed)
    samples <- list()
    union_genes <- c()
    for (i in seq(1, num_cluster)) {
        if (length(gene_list[[i]]) < num_sample[i]) {
            samples[[i]] <- gene_list[[i]]
        } else {
            samples[[i]] <- sample(gene_list[[i]], num_sample[i])
        }
        union_genes <- union(union_genes, samples[[i]])
    }
    union_genes <- unique(union_genes)

    return(list(
        union_genes = union_genes, num_sample = num_sample,
        num_gene = genes, gene_list = gene_list
    ))
}


#' @title numPts_below_line
#' @name numPts_below_line
#' @description function to calculate the number of data points below a given
#' line
#' @param myVector input sequence of sorted positive values from smallest to
#' greatest
#' @param slope slope to compare
#' @param x location point of the line to compare, integer from 1 to length of
#' myVector
#' @returns numeric
#' @details
#' This function calculates the number of data points in a sorted sequence
#' below a line with given slope through a certain point on this sequence.
#' @keywords internal
numPts_below_line <- function(
        myVector,
        slope,
        x) {
    yPt <- myVector[x]
    b <- yPt - (slope * x)
    xPts <- seq_along(myVector)
    return(sum(myVector <= (xPts * slope + b)))
}


#' @title filterSpatialGenes
#' @name filterSpatialGenes
#' @description function to filter gene list with existing spatial gene sets
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param spatial_genes input gene list
#' @param max max number of genes selected from spatial test
#' @param name name of spatial gene test for the filtering
#' @param method method of spatial gene selection
#' @returns list
#' @details
#' This function filters given gene list with the gene sets of selected
#' spatial gene test in Giotto,
#' also controls the total size of the gene set with given max number.
#' @keywords external
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' filterSpatialGenes(g, spatial_genes = "Gm19935")
#' @export
filterSpatialGenes <- function(gobject, spat_unit = NULL, feat_type = NULL,
    spatial_genes, max = 2500,
    name = c("binSpect", "silhouetteRank", "silhouetteRankTest"),
    method = c("none", "elbow")) {
    name <- match.arg(
        name,
        unique(c("binSpect", "silhouetteRank", "silhouetteRankTest", name))
    )
    method <- match.arg(method, unique(c("none", "elbow", method)))


    # NSE vars
    binSpect.pval <- silhouetteRank.score <- silhouetteRankTest.pval <-
        feat_ID <- NULL

    # first determine how many spatial genes in this dataset
    gx <- fDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)

    if (name == "binSpect") {
        gx <- gx[!is.na(binSpect.pval) & binSpect.pval < 1]
        gx_sorted <- gx[order(gx$binSpect.pval, decreasing = FALSE), ]
    } else if (name == "silhouetteRank") {
        gx <- gx[!is.na(silhouetteRank.score) & silhouetteRank.score > 0]
        gx_sorted <- gx[order(gx$silhouetteRank.score, decreasing = TRUE), ]
    } else if (name == "silhouetteRankTest") {
        gx <- gx[!is.na(silhouetteRankTest.pval) & silhouetteRankTest.pval < 1]
        gx_sorted <- gx[order(gx$silhouetteRankTest.pval, decreasing = FALSE), ]
    }

    # print(gx_sorted)
    if (method == "none") {
        if (name == "binSpect") {
            gx_sorted <- gx_sorted[binSpect.pval < 0.01]
        } else if (name == "silhouetteRankTest") {
            gx_sorted <- gx_sorted[silhouetteRankTest.pval < 0.01]
        }
        gx_sorted <- head(gx_sorted, n = max)
    } else if (method == "elbow") {
        y0 <- c()
        if (name == "binSpect") {
            y0 <- -log10(gx_sorted$binSpect.pval)
        } else if (name == "silhouetteRankTest") {
            y0 <- -log10(gx_sorted$silhouetteRankTest.pval)
        } else if (name == "silhouetteRank") {
            y0 <- gx_sorted$silhouetteRank.score
        }
        x0 <- seq(1, nrow(gx_sorted))

        y0s <- sort(y0)
        y0s[y0s < 0] <- 0 # strictly positive
        # plot(x0, y0)
        slope <- (max(y0s) - min(y0s)) / length(y0s) # This is the slope of the
        # line we want to slide. This is the diagonal.
        xPt <- floor(optimize(
            numPts_below_line,
            lower = 1, upper = length(y0s),
            myVector = y0s, slope = slope
        )$minimum)
        xPt <- length(y0s) - xPt
        y_cutoff <- y0[xPt] # The y-value at this x point. This is our y_cutoff.
        gx_sorted <- head(gx_sorted, n = xPt)
        message("\n Elbow method chosen to determine number of spatial genes.")
        cat(paste0(
            "\n Elbow point determined to be at x=", xPt, " genes",
            " y=", y_cutoff
        ))
    }

    # filter user's gene list (spatial_genes)
    gx_sorted <- gx_sorted[feat_ID %in% spatial_genes]

    num_genes_removed <- length(spatial_genes) - nrow(gx_sorted)

    return(list(
        genes = gx_sorted$feat_ID, num_genes_removed = num_genes_removed
    ))
}



#' @title chooseAvailableSpatialGenes
#' @name chooseAvailableSpatialGenes
#' @description function to find the test name for existing spatial gene sets
#' in Giotto
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @returns character
#' @details
#' This function outputs the available test name for existing spatial gene sets
#' in Giotto,
#' which could be used in parameter ‘name’ in `filterSpatialGenes`.
#' Priorities for showing the spatial gene test names are ‘binSpect’ >
#' ‘silhouetteRankTest’ > ‘silhouetteRank’.
#' @keywords internal
chooseAvailableSpatialGenes <- function(gobject,
    spat_unit = NULL, feat_type = NULL) {
    gx <- fDataDT(gobject, spat_unit = NULL, feat_type = NULL)
    eval1 <- "binSpect.pval" %in% names(gx)
    eval2 <- "silhouetteRankTest.pval" %in% names(gx)
    eval3 <- "silhouetteRank.score" %in% names(gx)
    if (eval1 == TRUE) {
        return("binSpect")
    } else if (eval2 == TRUE) {
        return("silhouetteRankTest")
    } else if (eval3 == TRUE) {
        return("silhouetteRank")
    } else {
        stop(paste0("\n No available spatial genes. Please run binSpect or
                    silhouetteRank\n"), call. = FALSE)
    }
}


#' @title checkAndFixSpatialGenes
#' @name checkAndFixSpatialGenes
#' @description function to check the selected test name for spatial gene set
#' in Giotto object
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param use_spatial_genes test name of spatial gene set to check
#' @param use_score logical variable to select silhouetteRank score
#' @returns character
#' @details
#' This function checks the user specified test name of spatial gene set in
#' Giotto object.
#' SilhouetteRank works only with score, and SilhouetteRankTest works only
#' with pval. Use parameter use_score to specify.
#' @keywords internal
checkAndFixSpatialGenes <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        use_spatial_genes,
        use_score = FALSE) {
    gx <- fDataDT(gobject, spat_unit = NULL, feat_type = NULL)

    if (use_spatial_genes == "silhouetteRank") {
        if (use_score == TRUE) {
            use_spatial_genes <- "silhouetteRank"
        } else {
            eval1 <- "silhouetteRank.score" %in% names(gx)
            eval2 <- "silhouetteRankTest.pval" %in% names(gx)
            if (eval1 == TRUE && eval2 == TRUE) {
                # if both evaluate to true, then decide by use_score.
                # silhouetteRank works only with score, silhouetteRankTest
                # works only with pval
                if (use_score == TRUE) {
                    use_spatial_genes <- "silhouetteRank"
                } else {
                    use_spatial_genes <- "silhouetteRankTest"
                }
            } else if (eval1 == TRUE) {
                use_spatial_genes <- "silhouetteRank"
            } else if (eval2 == TRUE) {
                use_spatial_genes <- "silhouetteRankTest"
            } else {
                stop(paste0("\n use_spatial_genes is set to silhouetteRank,
                            but it has not been run yet. Run silhouetteRank
                            first.\n"), call. = FALSE)
            }
        }
        return(use_spatial_genes)
    } else if (use_spatial_genes == "binSpect") {
        eval1 <- "binSpect.pval" %in% names(gx)
        if (eval1 == FALSE) {
            stop(paste0("\n use_spatial_genes is set to binSpect, but it has
                        not been run yet. Run binSpect first."),
                call. = FALSE
            )
        }
        return(use_spatial_genes)
    } else {
        stop("use_spatial_genes is set to one that is not supported.",
            call. = FALSE
        )
    }
}




#' @title initHMRF_V2
#' @name initHMRF_V2
#' @description Run initialization for HMRF model
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param spatial_network_name name of spatial network to use for HMRF
#' @param use_spatial_genes which of Giotto's spatial genes to use
#' @param use_score use score as gene selection criterion
#' (applies when use_spatial_genes=silhouetteRank)
#' @param gene_list_from_top total spatial genes before sampling
#' @param filter_method filter genes by top or by elbow method, prior to
#' sampling
#' @param user_gene_list user-specified genes (optional)
#' @param use_pca if PCA is used on the spatial gene expression value for
#' clustering
#' @param use_pca_dim dimensions of the PCs of the selected expression
#' @param gene_samples number of spatial gene subset to use for HMRF
#' @param gene_sampling_rate parameter (1-50) controlling proportion of gene
#' samples from different module when sampling, 1 corresponding to equal gene
#' samples between different modules; 50 corresponding to gene samples
#' proportional to module size.
#' @param gene_sampling_seed random number seed to sample spatial genes
#' @param use_metagene if metagene expression is used for clustering
#' @param cluster_metagene number of metagenes to use
#' @param top_metagene = number of genes in each cluster for the metagene
#' calculation
#' @param existing_spatial_enrichm_to_use name of existing spatial enrichment
#' result to use
#' @param use_neighborhood_composition if neighborhood composition is used for
#' hmrf
#' @param spatial_network_name_for_neighborhood spatial network used to
#' calculate neighborhood composition
#' @param metadata_to_use metadata used to calculate neighborhood composition
#' @param hmrf_seed random number seed to generate initial mean vector of HMRF
#' model
#' @param cl.method clustering method to calculate the initial mean vector,
#' selecting from 'km', 'leiden', or 'louvain'
#' @param resolution.cl resolution of Leiden or Louvain clustering
#' @param k number of HMRF domains
#' @param tolerance error tolerance threshold
#' @param zscore type of zscore to use
#' @param nstart number of Kmeans initializations from which to select the
#' best initialization
#' @param factor_step dampened factor step
#' @param python_path python_path
#' @returns initialized HMRF
#' @details
#' This function is the initialization step of HMRF domain clustering. First,
#' user specify which of Giotto's spatial genes to run,
#' through use_spatial_genes. Spatial genes have been stored in the gene
#' metadata table. A first pass of genes will filter genes that
#' are not significantly spatial, as determined by filter_method. If
#' filter_method is none, then top 2500 (gene_list_from_top) genes
#' ranked by pvalue are considered spatial. If filter_method is elbow, then the
#' exact cutoff is determined by the elbow in
#' the -log10 P-value vs. gene rank plot. Second, users have a few options to
#' decrease the dimension of the spatial genes for
#' clustering, listed with selection priority:
#'    1. use PCA of the spatial gene expressions (selected by use_pca)
#'    2. use metagene expressions (selected by use_metagene)
#'    3. sampling to select 500 spatial genes (controlled by gene_samples).
#' Third, once spatial genes are finalized, we are using clustering method to
#' initialize HMRF.
#' Instead of select spatial genes for domain clustering, HMRF method could
#' also applied on unit neighborhood composition of any group
#' membership(such as cell types), specified by parameter:
#' use_neighborhood_composition,  spatial_network_name_for_neighborhood and
#' metadata_to_use. Also HMRF provides the opportunity for user to do
#' clustering by any customized spatial enrichment matrix
#' (existing_spatial_enrichm_to_use).
#' There are 3 clustering algorithm: K-means, Leiden, and Louvain to determine
#' initial centroids of HMRF. The initialization is
#' then finished. This function returns a list containing y (expression),
#' nei (neighborhood structure), numnei (number of neighbors),
#' blocks (graph colors), damp (dampened factor), mu (mean),
#' sigma (covariance), k, genes, edgelist, init.cl (initial clusters),
#' spat_unit, feat_type. This information is needed for the second step, doHMRF.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- binSpect(g, return_gobject = TRUE)
#'
#' initHMRF_V2(gobject = g, cl.method = "km")
#' @export
initHMRF_V2 <-
    function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("scaled", "normalized", "custom"),
        spatial_network_name = "Delaunay_network",
        use_spatial_genes = c("binSpect", "silhouetteRank"),
        use_score = FALSE,
        gene_list_from_top = 2500,
        filter_method = c("none", "elbow"),
        user_gene_list = NULL,
        use_pca = FALSE,
        use_pca_dim = 1:20,
        gene_samples = 500,
        gene_sampling_rate = 2,
        gene_sampling_seed = 10,
        use_metagene = FALSE,
        cluster_metagene = 50,
        top_metagene = 20,
        existing_spatial_enrichm_to_use = NULL,
        use_neighborhood_composition = FALSE,
        spatial_network_name_for_neighborhood = NULL,
        metadata_to_use = NULL,
        hmrf_seed = 100,
        cl.method = c("km", "leiden", "louvain"),
        resolution.cl = 1,
        k = 10,
        tolerance = 1e-05,
        zscore = c("none", "rowcol", "colrow"),
        nstart = 1000,
        factor_step = 1.05,
        python_path = NULL) {
        wrap_msg(
            "\n If used in published research, please cite:
        Q Zhu, S Shah, R Dries, L Cai, GC Yuan.
        'Identification of spatially associated subpopulations by combining
        scRNAseq and sequential fluorescence in situ hybridization data'
        Nature biotechnology 36 (12), 1183-1190. 2018\n"
        )

        package_check(
            pkg_name = c(
                "tidygraph", "ggraph", "dplyr", "smfishHmrf", "graphcoloring"
            ),
            repository = c(
                "CRAN:tidygraph",
                "CRAN:ggraph",
                "CRAN:dplyr",
                "bitbucket:qzhudfci/smfishHmrf-r",
                "bitbucket:qzhudfci/graphcoloring"
            )
        )


        # DT vars
        to <- from <- clus <- NULL

        spat_unit <- set_default_spat_unit(
            gobject = gobject,
            spat_unit = spat_unit
        )
        feat_type <- set_default_feat_type(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        )

        gx <- fDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
        cx <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)

        spatial_network <- getSpatialNetwork(
            gobject,
            spat_unit = spat_unit,
            name = spatial_network_name,
            output = "networkDT",
            copy_obj = FALSE
        )
        spatial_network <- spatial_network[, .(to, from)]

        if (use_neighborhood_composition) {
            if (is.null(spatial_network_name_for_neighborhood)) {
                stop("\n spatial network is required to define neighborhood,
                    set with \'spatial_network_name_for_neighborhood\' \n",
                    call. = FALSE
                )
            } else if (is.null(metadata_to_use)) {
                stop("\n please specify the cluster in meta data, set with
                    \'metadata_to_use\' \n",
                    call. = FALSE
                )
            } else if (is.null(cx[[metadata_to_use]])) {
                stop("\n please provide a valid index in meta data, set with
                    \'metadata_to_use\'",
                    call. = FALSE
                )
            }

            cat(paste0(
                "\n use spatial network composition of \'",
                metadata_to_use, "\' for domain clustering"
            ))

            name.cl <- as.character(sort(unique(cx[[metadata_to_use]])))

            spatial_network_for_neighborhood <- getSpatialNetwork(
                gobject,
                spat_unit = spat_unit,
                name = spatial_network_name_for_neighborhood,
                output = "networkDT",
                copy_obj = FALSE
            )

            from.all <- c(
                spatial_network_for_neighborhood$from,
                spatial_network_for_neighborhood$to
            )
            to.all <- c(
                spatial_network_for_neighborhood$to,
                spatial_network_for_neighborhood$from
            )

            ct.tab <- aggregate(
                cx[[metadata_to_use]][match(
                    to.all, cx[["cell_ID"]]
                )],
                by = list(cell_ID = from.all), function(y) {
                    table(y)[name.cl]
                }
            )

            y0 <- ct.tab[, -1]
            y0[is.na(y0)] <- 0
            rownames(y0) <- ct.tab$cell_ID
            y0 <- y0 / rowSums(y0)
        } else if (!is.null(existing_spatial_enrichm_to_use)) {
            y0 <- getSpatialEnrichment(
                gobject,
                spat_unit = spat_unit, feat_type = feat_type,
                name = existing_spatial_enrichm_to_use,
                output = "data.table"
            )
            cell_ID_enrich <- y0$cell_ID
            y0 <- as.data.frame(y0[, -"cell_ID"])
            rownames(y0) <- cell_ID_enrich

            cat(paste0(
                "\n Spatial enrichment result: \'",
                existing_spatial_enrichm_to_use, "\' is used."
            ))

            if (sum(!rownames(y0) %in% cx$cell_ID) > 0) {
                stop("\n Rownames of selected spatial enrichment result do not
                    match to (a subset of) Cell IDs, please fix them.",
                    call. = FALSE
                )
            }
        } else {
            zscore <- match.arg(zscore, unique(c(
                "none", "rowcol", "colrow",
                zscore
            )))
            use_spatial_genes <- match.arg(use_spatial_genes, unique(c(
                "binSpect",
                "silhouetteRank", use_spatial_genes
            )))
            filter_method <- match.arg(filter_method, unique(c(
                "none",
                "elbow", filter_method
            )))
            values <- match.arg(expression_values, unique(c(
                "scaled",
                "normalized", "custom", expression_values
            )))
            expr_values <- getExpression(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                values = values, output = "matrix"
            )
            if (zscore != "none") {
                zscore <- match.arg(zscore, c("none", "colrow", "rowcol"))
                expr_values <- getExpression(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    values = "normalized",
                    output = "matrix"
                )
                if (zscore == "colrow") {
                    expr_values <- t(scale(t(scale(expr_values))))
                }
                if (zscore == "rowcol") {
                    expr_values <- scale(t(scale(t(expr_values))))
                }
            }

            spatial_genes <- c()

            if (!"binSpect.pval" %in% names(gx) &&
                !"silhouetteRank.score" %in% names(gx) &&
                !"silhouetteRankTest.pval" %in% names(gx)) {
                stop(paste0("\n Giotto spatial gene detection has not been run.
                    Please run spatial gene detection first: binSpect,
                    silhouetteRank."),
                    call. = FALSE
                )
            }

            if (!is.null(user_gene_list)) {
                message("\n User supplied gene list detected.")
                message("\n Checking user gene list is spatial...")

                use_spatial_genes <- chooseAvailableSpatialGenes(gobject)
                filtered <- filterSpatialGenes(
                    gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    spatial_genes = user_gene_list,
                    max = gene_list_from_top,
                    name = use_spatial_genes,
                    method = filter_method
                )
                if (filtered$num_genes_removed > 0) {
                    cat(paste0(
                        "\n Removed ", filtered$num_genes_removed,
                        " from user's input gene list due to being absent or
                        non-spatial genes."
                    ))
                    cat(paste0(
                        "\n Kept ", length(filtered$genes),
                        " spatial genes for next step"
                    ))
                }
                spatial_genes <- filtered$genes

                if (length(spatial_genes) == 0) {
                    stop("\n No genes are remaining to do HMRF. Please give a
                        larger gene list.",
                        call. = FALSE
                    )
                }
            } else {
                cat(paste0(
                    "\n Choosing spatial genes from the results of ",
                    use_spatial_genes
                ))
                use_spatial_genes <- checkAndFixSpatialGenes(
                    gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    use_spatial_genes = use_spatial_genes,
                    use_score = use_score
                )
                all_genes <- gx$feat_ID
                filtered <- filterSpatialGenes(
                    gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    spatial_genes = all_genes,
                    max = gene_list_from_top,
                    name = use_spatial_genes,
                    method = filter_method
                )
                cat(paste0(
                    "\n Kept ", length(filtered$genes),
                    " top spatial genes for next step"
                ))
                spatial_genes <- filtered$genes
            }

            if (use_pca == TRUE) {
                expr_values <- expr_values[spatial_genes, ]
                pc.expr <- prcomp(expr_values)[[2]]
                use_pca_dim <- use_pca_dim[
                    use_pca_dim %in% seq_len(ncol(pc.expr))
                ]
                y0 <- (pc.expr[, use_pca_dim])
            } else {
                message("\n Computing spatial coexpression modules...")
                spat_cor_netw_DT <- detectSpatialCorFeats(
                    gobject = gobject,
                    feat_type = feat_type,
                    spat_unit = spat_unit,
                    expression_values = values,
                    method = "network",
                    spatial_network_name = spatial_network_name,
                    subset_feats = spatial_genes,
                    network_smoothing = 0
                )

                if (use_metagene == FALSE) {
                    n <- min(gene_samples, 500, length(spatial_genes))
                    if (n < length(spatial_genes)) {
                        spat_cor_netw_DT <- clusterSpatialCorFeats(
                            spat_cor_netw_DT,
                            name = "spat_netw_clus",
                            k = 20
                        )
                        message("\n Sampling spatial genes from coexpression
                                modules...")
                        sample_genes <- sampling_sp_genes(
                            spat_cor_netw_DT$cor_clusters$spat_netw_clus,
                            sample_rate = gene_sampling_rate,
                            target = n,
                            seed = gene_sampling_seed
                        )
                        spatial_genes_selected <- sample_genes$union_genes
                        cat(paste0(
                            "\n Sampled ", length(spatial_genes_selected),
                            " genes."
                        ))
                    } else {
                        spatial_genes_selected <- spatial_genes
                    }
                    cat(paste0(
                        "\n Will use ", length(spatial_genes_selected),
                        " spatial genes for initialization of HMRF."
                    ))
                    expr_values <- expr_values[spatial_genes_selected, ]
                } else {
                    k.sp <- min(
                        ceiling(length(spatial_genes) / 20), cluster_metagene
                    )
                    if (k.sp < cluster_metagene) {
                        cat(paste0(
                            "\n construct ", k.sp,
                            " coexpression modules due to limited gene size..."
                        ))
                    }
                    spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT,
                        name = "spat_netw_clus", k = k.sp
                    )

                    cluster_genes_DT <- showSpatialCorFeats(spat_cor_netw_DT,
                        use_clus_name = "spat_netw_clus",
                        show_top_feats = 1
                    )

                    cat(paste0("\n Collecting top spatial genes and calculating
                        metagenes from ", k.sp, " coexpression modules..."))

                    top_per_module <- cluster_genes_DT[
                        , head(.SD, top_metagene),
                        by = clus
                    ]
                    cluster_genes <- top_per_module$clus
                    names(cluster_genes) <- top_per_module$feat_ID

                    meta.genes <- createMetafeats(
                        gobject,
                        spat_unit = spat_unit,
                        feat_type = feat_type,
                        expression_values = values,
                        feat_clusters = cluster_genes,
                        return_gobject = FALSE
                    )

                    expr_values <- t(meta.genes@enrichDT[, seq_len(k.sp)])
                    colnames(expr_values) <- unlist(
                        meta.genes@enrichDT[, "cell_ID"]
                    )
                    rownames(expr_values) <- paste0(
                        "metagene_", rownames(expr_values)
                    )
                }

                y0 <- t(as.matrix(expr_values))
            }
        }

        cell.rm <- setdiff(rownames(y0), unique(c(
            spatial_network$to,
            spatial_network$from
        )))

        if (length(cell.rm) > 0) {
            y0 <- y0[-match(cell.rm, rownames(y0)), ]
        }

        ##############################
        ## scale y matrix on each sample
        y <- t(scale(t(y0)))

        ## do not scale y if using neighborhood composition
        if (use_neighborhood_composition) {
            y <- y0
        }


        numcell <- dim(y)[1]
        m <- dim(y)[2]
        ncol.nei <- max(table(c(spatial_network$to, spatial_network$from)))
        nei <- matrix(-1, ncol = ncol.nei, nrow = numcell)
        rownames(nei) <- rownames(y)
        for (i in seq_len(numcell)) {
            nei.i <- c(spatial_network$from[spatial_network$to ==
                rownames(nei)[i]], spatial_network$to[spatial_network$from ==
                rownames(nei)[i]])
            if (length(nei.i) > 0) {
                nei[i, seq_along(nei.i)] <- sort(match(nei.i, rownames(y)))
            }
        }
        numnei <- as.integer(rowSums(nei != (-1)))
        nn <- nei
        numedge <- 0
        for (i in seq_len(numcell)) {
            numedge <- numedge + length(nn[i, nn[i, ] != -1])
        }
        edgelist <- matrix(0, nrow = numedge, ncol = 2)
        edge_ind <- 1
        for (i in seq_len(numcell)) {
            neighbors <- nn[i, nn[i, ] != -1]
            for (j in seq_along(neighbors)) {
                edgelist[edge_ind, ] <- c(i, neighbors[j])
                edge_ind <- edge_ind + 1
            }
        }
        message("\n Parsing neighborhood graph...")
        pp <- tidygraph::tbl_graph(
            edges = as.data.frame(edgelist), directed = FALSE
        )
        yy <- pp %>% dplyr::mutate(
            color = as.factor(graphcoloring::color_dsatur())
        )
        colors <- as.list(yy)$nodes$color
        cl_color <- sort(unique(colors))
        blocks <- lapply(cl_color, function(cl) {
            which(colors == cl)
        })

        cl.method <- tolower(cl.method)
        if (!cl.method %in% c("km", "leiden", "louvain")) {
            cl.method <- "km"
            message("clustering method not specified, use kmeans as default...")
        }

        if (cl.method == "km") {
            message("\n Kmeans initialization...")
            kk <- smfishHmrf::smfishHmrf.generate.centroid(
                y = y, par_k = k, par_seed = hmrf_seed,
                nstart = nstart
            )
            mu <- t(kk$centers)
            lclust <- lapply(seq_len(k), function(x) which(kk$cluster == x))
        } else {
            ##### need to double check leiden and louvain cluster functions
            gobject@dimension_reduction$cells$spatial <- NULL
            gobject@dimension_reduction$cells$spatial$spatial_feat <- NULL
            gobject@dimension_reduction$cells$spatial$spatial_feat$name <- "spatial_feat"
            gobject@dimension_reduction$cells$spatial$spatial_feat$reduction_method <- "spatial"
            gobject@dimension_reduction$cells$spatial$spatial_feat$coordinates <- y

            gobject <- createNearestNetwork(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                dim_reduction_to_use = "spatial",
                dim_reduction_name = "spatial_feat",
                dimensions_to_use = seq_len(ncol(y)),
                name = "sNN.initHMRF"
            )

            if (cl.method == "leiden") {
                message("\n Leiden clustering initialization...")
                leiden.cl <- doLeidenCluster(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    nn_network_to_use = "sNN",
                    network_name = "sNN.initHMRF",
                    set_seed = hmrf_seed,
                    return_gobject = FALSE,
                    python_path = python_path,
                    resolution = resolution.cl
                )
                cl.match <- leiden.cl$leiden_clus[
                    match(rownames(y), leiden.cl$cell_ID)
                ]
                mu <- aggregate(y, by = list(cl.match), FUN = mean)
            } else if (cl.method == "louvain") {
                message("\n Louvain clustering initialization...")
                louvain.cl <- doLouvainCluster(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    nn_network_to_use = "sNN",
                    network_name = "sNN.initHMRF",
                    set_seed = hmrf_seed,
                    return_gobject = FALSE,
                    python_path = python_path,
                    resolution = resolution.cl
                )
                cl.match <- louvain.cl$louvain_clus[
                    match(rownames(y), louvain.cl$cell_ID)
                ]
                mu <- aggregate(y, by = list(cl.match), FUN = mean)
            }

            rownames(mu) <- mu[, 1]
            mu <- t(mu[, -1])
            k <- dim(mu)[2]
            lclust <- lapply(colnames(mu), function(x) which(cl.match == x))
            cat(paste0("k is automatically identified as ", k))
        }


        damp <- array(0, c(k))
        sigma <- array(0, c(m, m, k))
        for (i in seq_len(k)) {
            sigma[, , i] <- cov(y[lclust[[i]], ])
            di <- smfishHmrf::findDampFactor(
                sigma[, , i],
                factor = factor_step,
                d_cutoff = tolerance, startValue = 1e-04
            )
            damp[i] <- ifelse(is.null(di), 0, di)
        }
        message("\n Done")
        list(
            y = y, nei = nei, numnei = numnei, blocks = blocks,
            damp = damp, mu = mu, sigma = sigma, k = k, genes = colnames(y),
            edgelist = edgelist, init.cl = lclust, spat_unit = spat_unit,
            feat_type = feat_type
        )
    }





#' @title doHMRF_V2
#' @name doHMRF_V2
#' @description function to run HMRF model
#' @param HMRF_init_obj initialization object list returned from initHMRF()
#' function
#' @param betas beta value of the HMRF model, controlling the smoothness of
#' clustering. NULL value of beta will provide default values based on feature
#' numbers, otherwise, a vector of three values: initial beta, beta increment,
#' and number of betas
#' @returns HMRF model
#' @details
#' This function will run a HMRF model after initialization of HMRF. Of note
#' is the beta parameter, the smoothing parameter.
#' If the users are interested in selecting results from different smoothness,
#' we recommend running a range of betas,
#' hence betas specify what this range is. For example, betas=c(0,10,5) will
#' run for the following betas: 0, 10, 20, 30, 40.
#' betas=c(0,5,2) will run for betas: 0, 5, 10. Setting the beta can use the
#' following guideline:
#' If number of features N is 10<N<50, set betas=c(0, 1, 20)
#' For 50<N<100, set betas=c(0, 2, 25)
#' For 100<N<500, set betas=c(0, 5, 20)
#' If user didn't specify the number or sequence of beta, the function will
#' provide a default value based on the number of features.
#' Returns a list of results for betas, spat_unit and feat_type. Result for
#' each beta is a list with probability(normalized or non-normalized), class,
#' and model log-likelihood value.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- binSpect(g, return_gobject = TRUE)
#' HMRF_init_obj <- initHMRF_V2(gobject = g, cl.method = "km")
#'
#' doHMRF_V2(HMRF_init_obj = HMRF_init_obj, betas = c(0, 5, 2))
#' @export
doHMRF_V2 <- function(HMRF_init_obj, betas = NULL) {
    message(
        "\n If used in published research, please cite:
    Q Zhu, S Shah, R Dries, L Cai, GC Yuan.
    'Identification of spatially associated subpopulations by combining
    scRNAseq and sequential fluorescence in situ hybridization data'
    Nature biotechnology 36 (12), 1183-1190. 2018"
    )

    message("\n Please find more explanation and instruction of the HMRF function
    on \n https://bitbucket.org/qzhudfci/smfishhmrf-r/src/master/TRANSITION.md")
    if (!"y" %in% names(HMRF_init_obj)) {
        stop("\n expression matrix 'y' not in the intialization object")
    }
    if (!"nei" %in% names(HMRF_init_obj)) {
        stop("\n neighbor matrix 'nei' not in the intialization object")
    }
    if (!"numnei" %in% names(HMRF_init_obj)) {
        stop("\n number of neighbors 'numnei' not in the intialization object")
    }
    if (!"blocks" %in% names(HMRF_init_obj)) {
        stop("\n iteration groups 'blocks' not in the intialization object")
    }
    if (!"damp" %in% names(HMRF_init_obj)) {
        stop("\n dampen factors 'damp' not in the intialization object")
    }
    if (!"mu" %in% names(HMRF_init_obj)) {
        stop("\n initial mean vector 'mu' not in the intialization object")
    }
    if (!"sigma" %in% names(HMRF_init_obj)) {
        stop("\n initial covariance matrix 'sigma' not in the intialization
            object")
    }
    if (!"k" %in% names(HMRF_init_obj)) {
        stop("\n cluster number 'k' not in the intialization object")
    }
    if (!"spat_unit" %in% names(HMRF_init_obj)) {
        HMRF_init_obj[["spat_unit"]] <- NULL
    }
    if (!"feat_type" %in% names(HMRF_init_obj)) {
        HMRF_init_obj[["feat_type"]] <- NULL
    }

    y <- HMRF_init_obj$y
    nei <- HMRF_init_obj$nei
    numnei <- HMRF_init_obj$numnei
    blocks <- HMRF_init_obj$blocks
    damp <- HMRF_init_obj$damp
    mu <- HMRF_init_obj$mu
    sigma <- HMRF_init_obj$sigma
    k <- HMRF_init_obj$k
    spat_unit <- HMRF_init_obj$spat_unit
    feat_type <- HMRF_init_obj$feat_type

    if (is.null(betas)) {
        beta_seq <- max(ceiling(ncol(y) / 10), 2)
        cat(paste0("\n Default value beta = ", beta_seq, " is used..."))
    } else if (length(betas) != 3 || (sum(betas[seq_len(3)] < 0) > 0)) {
        stop(wrap_txt(
            "\n please provide betas as a vector of 3 non-negative numbers
            (initial value, increment, total iteration number)",
            errWidth = TRUE
        ))
    } else {
        beta_init <- betas[1]
        beta_increment <- betas[2]
        beta_num_iter <- betas[3]
        beta_seq <- (seq_len(beta_num_iter) - 1) * beta_increment + beta_init
        beta_seq <- sort(unique(c(0, beta_seq)))
    }

    res <- c()
    for (beta_current in beta_seq) {
        print(sprintf("Doing beta=%.3f", beta_current))
        tc.hmrfem <- smfishHmrf::smfishHmrf.hmrfem.multi(
            y = y, neighbors = nei,
            beta = beta_current,
            numnei = numnei,
            blocks = blocks,
            mu = mu,
            sigma = sigma,
            verbose = TRUE,
            err = 1e-07,
            maxit = 50,
            dampFactor = damp
        )

        t_key <- sprintf("k=%d b=%.2f", k, beta_current)
        tc.hmrfem$sigma <- NULL
        tc.hmrfem$mu <- NULL
        rownames(tc.hmrfem$prob) <- rownames(y)
        rownames(tc.hmrfem$unnormprob) <- rownames(y)
        # names(tc.hmrfem$class) <- rownames(y)
        res[[t_key]] <- tc.hmrfem
    }
    result.hmrf <- res
    result.hmrf[["spat_unit"]] <- spat_unit
    result.hmrf[["feat_type"]] <- feat_type
    class(result.hmrf) <- append(class(result.hmrf), "HMRFoutput")
    return(result.hmrf)
}


#' @title addHMRF_V2
#' @name addHMRF_V2
#' @description function to add HMRF Domain Type to cell meta data
#' @param gobject giotto object
#' @param HMRFoutput result object from HMRF model
#' @param name name of HMRF models
#' @returns giotto object
#' @details
#' This function appends HMRF domain clusters to corresponding cell meta data
#' for all the beta values, with the given HMRF model names. For example, if
#' name = ‘hmrf1’ and name of result in HMRFoutput is ‘k=8 b=0.00’, the
#' appended cell meta data column will be named with ‘hmrf1 k=8 b=0.00’
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- binSpect(g, return_gobject = TRUE)
#' HMRF_init_obj <- initHMRF_V2(gobject = g, cl.method = "km")
#' HMRFoutput <- doHMRF_V2(HMRF_init_obj = HMRF_init_obj, betas = c(0, 5, 2))
#'
#' addHMRF_V2(gobject = g, HMRFoutput = HMRFoutput)
#' @export
addHMRF_V2 <- function(gobject, HMRFoutput, name = "hmrf") {
    if (!"HMRFoutput" %in% class(HMRFoutput)) {
        stop("\n HMRFoutput needs to be output from doHMRF_V2()")
    }
    if (!"spat_unit" %in% names(HMRFoutput)) {
        HMRFoutput[["spat_unit"]] <- NULL
    }
    if (!"feat_type" %in% names(HMRFoutput)) {
        HMRFoutput[["feat_type"]] <- NULL
    }
    spat_unit <- HMRFoutput$spat_unit
    feat_type <- HMRFoutput$feat_type

    cx <- getCellMetadata(
        gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    ordered_cell_IDs <- cx$cell_ID

    for (i in seq_along(grep("k", names(HMRFoutput)))) {
        gobject <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = HMRFoutput[[i]]$class[match(
                ordered_cell_IDs,
                rownames(HMRFoutput[[i]]$prob)
            )],
            vector_name = paste(name, names(HMRFoutput)[i])
            # ,column_cell_ID = 'cell_ID',
            # by_column = TRUE
        )
    }
    return(gobject)
}



#' @title viewHMRFresults_V2
#' @name viewHMRFresults_V2
#' @description function to view HMRF results with multiple betas
#' @param gobject giotto object
#' @param k number of clusters in hmrf results
#' @param betas values of beta to plot
#' @param hmrf_name name of hmrf models
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param third_dim index of 3D plot
#' @param cow_n_col number of columns to show in figure
#' @param cow_rel_h relative height
#' @param cow_rel_w relative width
#' @param cow_align cowplot alignment parameter
#' @param show_plot option to show plot
#' @param save_plot option to save plot
#' @param return_plot if function return a plot
#' @param default_save_name names of figure file to save
#' @param save_param other saving parameters
#' @param \dots additional params to pass to plotting
#' @returns spatial plots with HMRF domains
#' @details
#' This function plots spatial map of HMRF domain clusters for multiple beta
#' with the name (hmrf_name),
#' matching the first part of the cell meta column names with HMRF clusters
#' (for example name of ‘hmrf1 k=8 b=0.00’ is ‘hmrf1’)
#' @export
viewHMRFresults_V2 <-
    function(
        gobject, k, betas,
        hmrf_name,
        spat_unit = NULL,
        feat_type = NULL,
        third_dim = FALSE,
        cow_n_col = 2,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        show_plot = TRUE,
        save_plot = TRUE,
        return_plot = TRUE,
        default_save_name = "HMRF_result",
        save_param = list(),
        ...) {
        # beta_seq = round(betas,digits = 2)
        # t_key = paste0(hmrf_name,'_k', k, '_b.',beta_seq)
        t_key <- paste(hmrf_name, sprintf("k=%d b=%.2f", k, betas))

        meta_names <- colnames(combineMetadata(
            gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
        ))

        if (length(setdiff(t_key, meta_names)) > 0) {
            beta_null <- paste(betas[which(!t_key %in% meta_names)],
                collapse = ","
            )
            stop(paste0(
                '\n HMRF result "', hmrf_name, '" of k = ', k,
                ", beta = ", beta_null,
                " was not found in the Giotto object."
            ))
        }

        savelist <- list()
        for (kk in seq_along(t_key))
        {
            if (third_dim == TRUE) {
                pl <- spatPlot3D(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    cell_color = t_key[kk],
                    show_plot = FALSE,
                    save_plot = FALSE,
                    title = t_key[kk],
                    default_save_name = "HMRF_result",
                    return_plot = TRUE,
                    ...
                )
            } else {
                pl <- spatPlot2D(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    cell_color = t_key[kk],
                    show_plot = FALSE,
                    save_plot = FALSE,
                    title = t_key[kk],
                    cow_n_col = 1,
                    cow_rel_h = 1,
                    cow_rel_w = 1,
                    cow_align = "h",
                    default_save_name = "HMRF_result",
                    return_plot = TRUE,
                    ...
                )
            }
            savelist[[kk]] <- pl
        }

        if (cow_n_col > length(savelist)) {
            cow_n_col <- length(savelist)
        }

        # combine plots with cowplot
        combo_plot <- plot_grid(
            plotlist = savelist,
            ncol = cow_n_col,
            rel_heights = cow_rel_h,
            rel_widths = cow_rel_w,
            align = cow_align
        )

        # assign default ncol and nrow if not in save_param
        if (!"ncol" %in% names(save_param)) {
            save_param[["ncol"]] <- cow_n_col
        }
        if (!"nrow" %in% names(save_param)) {
            save_param[["nrow"]] <- ceiling(length(savelist) / cow_n_col)
        }

        # output plot
        return(GiottoVisuals::plot_output_handler(
            gobject = gobject,
            plot_object = combo_plot,
            save_plot = save_plot,
            return_plot = return_plot,
            show_plot = show_plot,
            default_save_name = default_save_name,
            save_param = save_param,
            else_return = NULL
        ))
    }

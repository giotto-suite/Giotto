#' Create Stereo-seq Giotto Object
#'
#' @param stereoseq_dir filepath to the exported Stereo-seq directory.
#' @param type character. Use "squarebin" to read expression at different bin 
#' levels (default), "cellbin" to read expression at cell resolution, or 
#' "subcellular" to read individual transcripts and cell boundaries.
#' @param bin_size bin size to select from *.tissue.gef file. Choose a value 
#' from "bin1", "bin5", "bin10", "bin20", "bin50", "bin100", "bin150", or 
#' "bin200". Only needed when using type = "squarebin".
#' @param gene_column (optional) character. Which column contains the gene names
#' within the geneExp information. Choose from "geneName" (default), or "geneID"
#' @param verbose logical. Be verbose.
#' @param h5_file (optional) name to create an on-disk HDF5 file.
#' @param instructions list of instructions or output result
#' @param negative_y logical. Map data to negative y spatial values during 
#' automatic alignment (Default = TRUE). Meaning that origin is in upper left 
#' instead of lower left.
#' from \code{\link{createGiottoInstructions}}
#' @returns Giotto Stereo-seq object
#' @export
createGiottoStereoSeqObject <- function(
        stereoseq_dir,
        type = "squarebin", 
        bin_size = "bin100",
        gene_column = "geneName",
        negative_y = TRUE,
        shift_polygon_y = 0,
        verbose = TRUE,
        h5_file = NULL,
        instructions = NULL) {
    # data.table vars
    genes <- gene_idx <- x <- y <- sdimx <- sdimy <- cell_ID <- bin_ID <-
        count <- i.bin_ID <- NULL
    
    # package check
    package_check(pkg_name = "rhdf5", repository = "Bioc")
    
    # directory check
    if (!file.exists(stereoseq_dir)) stop(
        "Path to Stereo-seq directory does not exist")
    
    dir_files <- list.files(
        file.path(stereoseq_dir, "outs", "feature_expression"))
    
    # file reading type check
    if(!type %in% c("squarebin", "cellbin", "subcellular")) stop(
        "'type' should be either 'squarebin', 'cellbin', or 'subcellular'")
    
    # Read squarebin
    if(type == "squarebin") {
        expression_file <- file.path(
            stereoseq_dir, "outs", "feature_expression", 
            dir_files[grep(".tissue.gef", dir_files)])
        
        if (!file.exists(expression_file)) stop(
            "Path to expression file ", 
            file.path(
                stereoseq_dir, "outs", "feature_expression", ".tissue.gef"),
            " does not exist")
        
        # check if proper bin_size is selected. These are determined in SAW pipeline
        vmsg(.v = verbose, "Reading expression file... \n")
        
        bin_size_options <- c("bin1", "bin10", "bin20", "bin50", "bin100", "bin200")
        if (!(bin_size %in% bin_size_options)) {
            stop("Please select valid bin size, see ?createGiottoSTOmicsObject for 
            details.")
        }
        
        # check valid gene_column value
        if(!gene_column %in% c("geneName", "geneID")) stop(
            "Provide a valid 'gene_column' value. It should be 'geneName' or
            'geneID'")
        
        # 1. read tissue.gef file at specific bin size
        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = paste0("geneExp/", bin_size)
        )
        
        exprDT <- data.table::setDT(geneExpData[["expression"]])
        geneDT <- data.table::setDT(geneExpData[["gene"]])
        exprDT[, `:=`(genes, rep(x = geneDT[[gene_column]], geneDT$count))]
        
        vmsg(.v = verbose, "Finished reading in tissue.gef")
    }
    
    # Read cellbin
    if(type == "cellbin") {
        expression_file <- file.path(
            stereoseq_dir, "outs", "feature_expression", 
            dir_files[grep(".adjusted.cellbin.gef", dir_files)])
        
        if (!file.exists(expression_file)) stop(
            "Path to expression file ", 
            file.path(
                stereoseq_dir, "outs", "feature_expression", ".adjusted.cellbin.gef"),
            " does not exist")
        
        # 1. read .adjusted.cellbin.gef file 
        vmsg(.v = verbose, "Reading expression file... \n")
        
        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = "cellBin"
        )
        
        exprDT <- data.table::setDT(geneExpData[["geneExp"]])
        geneDT <- data.table::setDT(geneExpData[["gene"]])
        
        # check valid gene_column value
        if(!gene_column %in% c("geneName", "geneID")) stop(
            "Provide a valid 'gene_column' value. It should be 'geneName' or
            'geneID'")
        
        exprDT[, `:=`(genes, rep(x = geneDT[[gene_column]], geneDT$cellCount))]
        
        vmsg(.v = verbose, "Finished reading in adjusted.cellbin.gef")
    }
    
    # Read subcellular
    if(type == "subcellular") {
        
        vmsg(.v = verbose, "Reading transcripts file ...\n")
        
        expression_file <- file.path(
            stereoseq_dir, "outs", "feature_expression", 
            dir_files[grep(".adjusted.cellbin.gef", dir_files)])
        
        if (!file.exists(expression_file)) stop(
            "Path to file ", 
            file.path(
                stereoseq_dir, "outs", "feature_expression", 
                ".adjusted.cellbin.gef"),
            " does not exist")
        
        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = "cellBin"
        )
        geneDT <- data.table::setDT(geneExpData[["gene"]])
        
        if (!file.exists(expression_file)) stop(
            "Path to transcripts file ", 
            file.path(
                stereoseq_dir, "outs", "feature_expression", 
                "_raw_barcode_gene_exp.txt"),
            " does not exist")
        
        transcripts_file <- file.path(
            stereoseq_dir, "outs", "feature_expression", 
            dir_files[grep("_raw_barcode_gene_exp.txt", dir_files)])
        
        # read transcripts file with gene locations
        transcripts <- data.table::fread(transcripts_file)
        
        vmsg(.v = verbose, "Reading cell boundaries...\n")
        
        files_list <- list.files(file.path(stereoseq_dir, "outs", "image"))
        mask_file <- files_list[grep("_HE_mask.tif$", files_list)]
        
        # Check that mask file exists
        if(is.null(mask_file)) {
            stop("The expected mask file was not found in the image directory.
                 Make sure a *_HE_mask.tif files is located under outs/image/")
        }
        
        # Read polygons from mask
        mask_poly <- createGiottoPolygonsFromMask(
            maskfile = file.path(stereoseq_dir, "outs", "image", mask_file),
            calc_centroids = TRUE
        )
        
    }
    
    # 2. create spatial locations
    vmsg(.v = verbose, "Creating spatial_locations... \n")

    if(type == "squarebin") {
        spatial_locations <- unique(exprDT[, c("x", "y")], by = c("x", "y"))
        spatial_locations[, bin_ID := seq_len(nrow(spatial_locations))]
        
        vmsg(.v = verbose, nrow(spatial_locations), " bins in total \n")
    }
    
    if(type == "cellbin") {
        cellDT <- data.table::setDT(geneExpData[["cell"]])
        cellDT[, cell_ID := paste0("cell_", id)]
        spatial_locations <- cellDT[, .(cell_ID, x, y)]
        
        vmsg(.v = verbose, nrow(spatial_locations), " cells in total \n")
    }
    
    if(type == "subcellular") {
        
        if(isTRUE(negative_y)) {
            mask_poly <- spatShift(
                mask_poly,
                dy = -(terra::ext(mask_poly)[3] + terra::ext(mask_poly)[4]))
        }
        
        mask_poly <- spatShift(mask_poly, dy = shift_polygon_y)
        
        spatial_locations <- as.data.frame(
            terra::geom(terra::centroids(mask_poly)))
        spatial_locations$cell_ID <- paste0("cell_", spatial_locations$geom)
        spatial_locations <- spatial_locations[, c("x", "y", "cell_ID")]
        
        spat_locs <- createSpatLocsObj(coordinates = spatial_locations,
                                       spat_unit = "cell",
                                       name = "raw")
        
        vmsg(.v = verbose, nrow(spatial_locations), 
            " polygons (cells) in total \n")
        
        # merge with geneID column
        transcript_locs <- merge(transcripts, geneDT[, .(geneID, geneName)])
        transcript_locs <- as.data.frame(transcript_locs)[
            , c("x", "y", gene_column, "readCount")]
        colnames(transcript_locs)[3] <- "feat_ID"
        
        # Create giotto points
        g_points <- createGiottoPoints(x = transcript_locs)
        
        vmsg(.v = verbose, nrow(transcript_locs), " transcripts in total \n")
        
        if(isTRUE(negative_y)) {
            g_points <- flip(g_points, direction = "vertical")
        }
        
        g_points <- crop(g_points, terra::ext(mask_poly))
        
    }
    
    vmsg(.v = verbose, "Finished spatial_locations \n")

    # Create expression matrix
    if(type %in% c("squarebin", "cellbin")) {
        vmsg(.v = verbose, "Creating expression matrix... \n")
        
        if(type == "squarebin") {
            exprDT <- merge(exprDT, spatial_locations, by = c("x", "y"))
            exprDT <- exprDT[, .(bin_ID, genes, count)]
            exprDT[, bin_ID := as.integer(bin_ID)]
            
            expMatrix <- data.table::dcast(
                exprDT, 
                formula = bin_ID ~ genes, 
                value.var = 'count', 
                fun.aggregate = sum)
            
            rownames_matrix <- paste0("cell_", expMatrix$bin_ID)
            expMatrix <- as.matrix(expMatrix[, -1, with = FALSE])
            rownames(expMatrix) <- rownames_matrix
            expMatrix <- t(expMatrix)
            
            spatial_locations[, cell_ID := paste0("cell_", bin_ID)]
            spatial_locations <- spatial_locations[, c("x", "y", "cell_ID")]
        }
        
        if(type == "cellbin") {
            exprDTmat <- exprDT[, .(cellID, genes, count)]
            exprDTmat[, cellID := as.integer(cellID)]
            expMatrix <- data.table::dcast(
                exprDTmat, 
                formula = cellID ~ genes, 
                value.var = 'count', 
                fun.aggregate = sum)
            
            rownames_matrix <- paste0("cell_", expMatrix$cellID)
            expMatrix <- as.matrix(expMatrix[, -1, with = FALSE])
            rownames(expMatrix) <- rownames_matrix
            expMatrix <- t(expMatrix)
        }
        
        vmsg(.v = verbose, "finished expression matrix")
    }
    
    # Create Giotto object
    vmsg(.v = verbose, "Creating giotto object... \n")
    
    if(type %in% c("squarebin", "cellbin")) {
        # ensure first non-numerical col is cell_ID
        spatial_locations[, x := as.integer(x)]
        spatial_locations[, y := as.integer(y)]
        
        if(isTRUE(negative_y)) {
            spatial_locations[, y := 0 - y]
        }
        
        stereo <- createGiottoObject(
            expression = expMatrix,
            spatial_locs = spatial_locations,
            verbose = verbose,
            h5_file = h5_file,
            instructions = instructions
        )
    }
    
    if(type == "subcellular") {
        
        stereo <- giotto()
        
        # Add giotto points 
        stereo <- setGiotto(stereo, g_points)
        
        # Add giotto polygons 
        stereo <- setGiotto(stereo, mask_poly)
        
        # Add polygon centroids
        stereo <- setGiotto(stereo, spat_locs)
        
    }

    # 5. add image
    vmsg(.v = verbose, "Attaching HE image... \n")

    image_dir <- file.path(stereoseq_dir, "outs", "image")
    he_image_path <- list.files(
        path = image_dir, pattern = "HE_regist", full.names = TRUE)
    gimg <- createGiottoLargeImage(he_image_path, name = "image")
    
    # if(type == "subcellular") {
    #     terra::ext(gimg) <- terra::ext(stereo)
    # }

    stereo <- addGiottoLargeImage(
        gobject = stereo,
        largeImages = gimg,
        negative_y = negative_y,
        verbose = verbose
    )
    
    vmsg(.v = verbose, "Finished giotto object... \n")
    return(stereo)
}

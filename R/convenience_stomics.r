#' Create STOmics Giotto Object
#'
#' @param stomics_dir filepath to the exported STOmics directory.
#' @param type character. Use "squarebin" to read expression from *.tissue.gef 
#' file (default) or "cellbin" to read expression from *.adjusted.cellbin.gef 
#' file.
#' @param bin_size bin size to select from *.tissue.gef file. Choose a value 
#' from "bin1", "bin5", "bin10", "bin20", "bin50", "bin100", "bin150", or 
#' "bin200". Only needed when using type = "squarebin".
#' @param gene_column (optional) character. Which column contains the gene names
#' within the geneExp information. Choose from "geneName" (default), or "geneID"
#' @param verbose logical. Be verbose.
#' @param h5_file (optional) name to create an on-disk HDF5 file.
#' @param instructions list of instructions or output result
#' @param flip_spatial_locs Flip spatial locations in the y axis
#' from \code{\link{createGiottoInstructions}}
#' @returns Giotto STOmics object
#' @export
createGiottoSTOmicsObject <- function(
        stomics_dir,
        type = c("squarebin", "cellbin"),
        bin_size = "bin100",
        gene_column = c("geneName", "geneID"),
        flip_spatial_locs = TRUE,
        verbose = TRUE,
        h5_file = NULL,
        instructions = NULL) {
    # data.table vars
    genes <- gene_idx <- x <- y <- sdimx <- sdimy <- cell_ID <- bin_ID <-
        count <- i.bin_ID <- NULL
    
    # package check
    package_check(pkg_name = "rhdf5", repository = "Bioc")
    
    # directory check
    if (!file.exists(stomics_dir)) stop(
        "Path to stomics directory does not exist")
    
    dir_files <- list.files(file.path(stomics_dir, "outs", "feature_expression"))
    
    # file reading type check
    if(!type %in% c("squarebin", "cellbin")) stop(
        "'type' value should be either 'squarebin' or 'cellbin'")
    
    # Read squarebin
    if(type == "squarebin") {
        expression_file <- file.path(
            stomics_dir, "outs", "feature_expression", 
            dir_files[grep(".tissue.gef", dir_files)])
        
        if (!file.exists(expression_file)) stop(
            "Path to expression file ", expression_file, " does not exist")
        
        # check if proper bin_size is selected. These are determined in SAW pipeline
        wrap_msg("1. Reading expression file... \n")
        bin_size_options <- c("bin1", "bin10", "bin20", "bin50", "bin100", "bin200")
        if (!(bin_size %in% bin_size_options)) {
            stop("Please select valid bin size, see ?createGiottoSTOmicsObject for 
            details.")
        }
        
        # 1. read tissue.gef file at specific bin size
        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = paste0("geneExp/", bin_size)
        )
        
        exprDT <- data.table::as.data.table(geneExpData[["expression"]])
        exprDT[, count := lapply(.SD, as.integer), .SDcols = "count"]
        data.table::setorder(exprDT, x, y) # sort by x, y coords (ascending)
        geneDT <- data.table::as.data.table(geneExpData[["gene"]])
        
        # check valid gene_column value
        if(!gene_column %in% c("geneName", "geneID")) stop(
            "Provide a valid 'gene_column' value. It should be 'geneName' or
            'geneID'")
        
        # process duplicated gene symbol
        if (any(duplicated(geneDT[[gene_column]]))) {
            duplicated_genes <- unique(
                geneDT[[gene_column]][duplicated(geneDT[[gene_column]])])
            cat("Ops!!! Duplicated_genes, processing: sum(count), mean(offset)\n")
            # merge
            for (gene in duplicated_genes) {
                # indices
                idx <- which(geneDT[[gene_column]] == gene)
                
                cat("Original count values for", gene, ":", geneDT$count[idx], "\n")
                
                # update
                geneDT$count[idx[1]] <- sum(geneDT$count[idx]) # 对重复的 count 求和
                geneDT$offset[idx[1]] <- mean(geneDT$offset[idx]) # 对重复的 offset 求平均
                
                #
                cat("Updated count for", gene, ":", geneDT$count[idx[1]], "\n")
                cat("Updated offset for", gene, ":", geneDT$offset[idx[1]], "\n")
                
                #
                geneDT <- geneDT[-idx[-1], ]
            }
        }
        
        if (isTRUE(verbose)) wrap_msg(
            "Finished reading in tissue.gef", bin_size)
    }
    
    if(type == "cellbin") {
        expression_file <- file.path(
            stomics_dir, "outs", "feature_expression", 
            dir_files[grep(".adjusted.cellbin.gef", dir_files)])
        
        if (!file.exists(expression_file)) stop(
            "Path to expression file ", expression_file, " does not exist")
        
        # check if proper bin_size is selected. These are determined in SAW pipeline
        wrap_msg("1. Reading expression file... \n")
        
        # 1. read tissue.gef file at specific bin size
        geneExpData <- rhdf5::h5read(
            file = expression_file,
            name = "cellBin"
        )
        
        exprDT <- data.table::as.data.table(geneExpData[["cellExp"]])
        exprDT_cell <- data.table::as.data.table(geneExpData[["geneExp"]])
        exprDT[, cell_ID := exprDT_cell[["cellID"]]]
        data.table::setorder(exprDT, geneID) # sort by geneID (ascending)
        geneDT <- data.table::as.data.table(geneExpData[["gene"]])
        
        # check valid gene_column value
        if(!gene_column %in% c("geneName", "geneID")) stop(
            "Provide a valid 'gene_column' value. It should be 'geneName' or
            'geneID'")
        
        # process duplicated gene symbol
        if (any(duplicated(geneDT[[gene_column]]))) {
            duplicated_genes <- unique(
                geneDT[[gene_column]][duplicated(geneDT[[gene_column]])])
            cat("Ops!!! Duplicated_genes, processing: sum(count), mean(offset)\n")
            # merge
            for (gene in duplicated_genes) {
                # indices
                idx <- which(geneDT[[gene_column]] == gene)
                
                cat("Original count values for", gene, ":", geneDT$cellCount[idx], "\n")
                
                # update
                geneDT$cellCount[idx[1]] <- sum(geneDT$cellCount[idx]) # 对重复的 count 求和
                geneDT$offset[idx[1]] <- mean(geneDT$offset[idx]) # 对重复的 offset 求平均
                
                #
                cat("Updated count for", gene, ":", geneDT$cellCount[idx[1]], "\n")
                cat("Updated offset for", gene, ":", geneDT$offset[idx[1]], "\n")
                
                #
                geneDT <- geneDT[-idx[-1], ]
            }
        }
        
        if (isTRUE(verbose)) wrap_msg(
            "Finished reading in adjusted.cellbin.gef")
    }
    
    # 2. create spatial locations
    if (isTRUE(verbose)) wrap_msg("2. create spatial_locations... \n")
    
    if(type == "squarebin") {
        cell_locations <- unique(exprDT[, c("x", "y")], by = c("x", "y"))
        cell_locations[, bin_ID := as.factor(seq_len(nrow(cell_locations)))]
        cell_locations[, cell_ID := paste0("cell_", bin_ID)]
        data.table::setcolorder(cell_locations, c("x", "y", "cell_ID", "bin_ID"))
        # ensure first non-numerical col is cell_ID
        if (isTRUE(verbose)) wrap_msg(nrow(cell_locations), " bins in total \n")
    }
    
    if(type == "cellbin") {
        cell_locations <- unique(
            geneExpData[["cell"]][, c("x", "y", "id")], by = c("x", "y"))
        cell_locations <- data.table::as.data.table(cell_locations)
        cell_locations[, cell_ID := paste0("cell_", id)]
        cell_locations <- cell_locations[, c("x", "y", "cell_ID")]
        
        # ensure first non-numerical col is cell_ID
        data.table::setcolorder(cell_locations, c("x", "y", "cell_ID"))
        cell_locations[, x := lapply(.SD, as.integer), .SDcols = "x"]
        cell_locations[, y := lapply(.SD, as.integer), .SDcols = "y"]
        
        if (isTRUE(verbose)) wrap_msg(nrow(cell_locations), " cells in total \n")
    }
    
    if(isTRUE(flip_spatial_locs)) {
        cell_locations[, y := 0 - y]
    }
    
    if (isTRUE(verbose)) wrap_msg("finished spatial_locations \n")
    
    # 3. create expression matrix
    if (isTRUE(verbose)) wrap_msg("3. create expression matrix... \n")
    
    if(type == "squarebin") {
        exprDT[, genes := as.character(
            rep(x = geneDT[[gene_column]], geneDT$count))]
        
        exprDT[, gene_idx := as.integer(factor(exprDT$genes,
                                               levels = unique(exprDT$genes)
        ))]
        
        # merge on x,y and populate based on bin_ID values in cell_locations
        exprDT[cell_locations, cell_ID := i.bin_ID, on = .(x, y)]
        # exprDT$cell_ID <- as.integer(exprDT$cell_ID)
        exprDT[, cell_ID := paste0("cell_", cell_ID)]
        exprDT <- exprDT[, c("genes", "cell_ID", "count")]

        exprDT_wide <- data.table::dcast(
            exprDT, genes ~ cell_ID, value.var = "count",
            fun.aggregate = sum)
        
        expMatrix <-  Matrix::Matrix(Matrix::as.matrix(exprDT_wide[,-1]),
                                     sparse = TRUE)
        
        colnames(expMatrix) <- colnames(exprDT_wide)[-1]
        rownames(expMatrix) <- exprDT_wide$genes
        

        # expMatrix <- Matrix::sparseMatrix(
        #     i = exprDT$gene_idx,
        #     j = exprDT$cell_ID,
        #     x = exprDT$count
        # )
        # 
        # colnames(expMatrix) <- cell_locations$cell_ID
        # rownames(expMatrix) <- geneDT[[gene_column]]
    }
    
    if(type == "cellbin") {
        exprDT[, genes := as.character(
            rep(x = geneDT[[gene_column]], geneDT$cellCount))]
        
        exprDT[, cell_ID := paste0("cell_", cell_ID)]
        
        exprDT_wide <- data.table::dcast(
            exprDT, genes ~ cell_ID, value.var = "count",
            fun.aggregate = sum)
        
        expMatrix <-  Matrix::Matrix(Matrix::as.matrix(exprDT_wide[,-1]),
                                     sparse = TRUE)
        
        colnames(expMatrix) <- colnames(exprDT_wide)[-1]
        rownames(expMatrix) <- exprDT_wide$genes
    }
    
    rm(exprDT)
    if (isTRUE(verbose)) wrap_msg("finished expression matrix")
    
    # 4. create minimal giotto object
    if (isTRUE(verbose)) wrap_msg("4. create giotto object... \n")
    stereo <- createGiottoObject(
        expression = expMatrix,
        spatial_locs = cell_locations,
        verbose = verbose,
        h5_file = h5_file,
        instructions = instructions
    )
    
    # 5. add image
    if (isTRUE(verbose)) wrap_msg("5. attaching HE image... \n")
    image_dir <- file.path(stomics_dir, "outs", "image")
    he_image_path <- list.files(path = image_dir, pattern = "HE_regist", full.names = TRUE)
    gimg <- createGiottoLargeImage(he_image_path, name = "HE_regist")
    stereo <- addGiottoLargeImage(
        gobject = stereo,
        largeImages = gimg,
        negative_y = TRUE,
        verbose = verbose
    )

    if (isTRUE(verbose)) wrap_msg("finished giotto object... \n")
    return(stereo)
}

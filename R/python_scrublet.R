#' @title doScrubletDetect
#' @name doScrubletDetect
#' @description Run *scrublet* doublet detection for raw expression. Intended
#' for single cell data
#' @param gobject giotto object containing expression data
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param expected_doublet_rate expected transcriptomes that are doublets. 0.06
#' is from 10x Chromium guide.
#' @param min_counts scrublet internal data filtering, min counts found to be
#' considered a cell
#' @param min_cells scrublet internal data filtering. min cells expressed to be
#' considered a feat
#' @param min_gene_variability_pctl scrublet internal PCA generation. highly
#' variable gene percentile cutoff
#' @param n_prin_comps number of PCs to use in PCA for detection
#' @param return_gobject return as gobject if TRUE, data.frame with cell_ID if
#' FALSE
#' @param seed If a numeric is provided, then it will be used as a seed. If
#' NULL, no seed will be set.
#' @seealso This function wraps the python package scrublet
#' \doi{10.1016/j.cels.2018.11.005}
#' @returns if `return_gobject = FALSE`, a `data.table` cell_ID, doublet scores,
#' and classifications are returned. If `TRUE`, that information is appended
#' into the input `giotto` object's metadata and the `giotto` object is
#' returned.
#' @md
#' @examples
#' # Should only be done with single cell data, but this is just a
#' # convenient example.
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' g <- doScrubletDetect(g)
#'
#' pDataDT(g) # doublet_scores and doublet cols are added
#' dimPlot2D(g, cell_color = "doublet_scores", color_as_factor = FALSE)
#' @export
doScrubletDetect <- function(
        gobject,
        feat_type = NULL,
        spat_unit = "cell",
        expression_values = "raw",
        expected_doublet_rate = 0.06,
        min_counts = 1,
        min_cells = 1,
        min_gene_variability_pctl = 85,
        n_prin_comps = 30,
        return_gobject = TRUE,
        seed = 1234) {
    # verify if optional package is installed
    package_check(
        pkg_name = "scrublet",
        repository = "pip"
    )

    # print message with information #
    message("using 'scrublet' to detect doublets. If used in published
    research, please cite: \n
    Wolock, S. L., Lopez, R. & Klein, A. M.
    Scrublet: Computational Identification of Cell Doublets in Single-Cell
    Transcriptomic Data. Cell Syst. 8, 281-291.e9 (2019).
    https://doi.org/10.1016/j.cels.2018.11.005")

    # prepare python path and scrublet_script
    python_path <- readGiottoInstructions(gobject, param = "python_path")
    reticulate::use_python(required = TRUE, python = python_path)
    python_scrublet_function <- system.file(
        "python", "python_scrublet.py",
        package = "Giotto"
    )
    reticulate::source_python(file = python_scrublet_function, convert = TRUE)

    # set seed
    if (!is.null(seed)) {
        seed_number <- as.numeric(seed)
        reticulate::py_set_seed(
            seed = seed_number,
            disable_hash_randomization = TRUE
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


    # 1. convert input to char for python inputs that must be type int
    min_counts <- as.character(min_counts)
    min_cells <- as.character(min_cells)
    min_gene_variability_pctl <- as.character(min_gene_variability_pctl)
    n_prin_comps <- as.character(n_prin_comps)

    # 2. get expression data
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "matrix"
    )

    # input is a sparse matrix with cells as rows and genes as columns
    scr_input <- to_scipy_sparse(expr_values, format = "C", transpose = TRUE)

    scrublet_out <- python_scrublet(
        counts_matrix = scr_input,
        expected_doublet_rate = expected_doublet_rate,
        min_counts = min_counts,
        min_cells = min_cells,
        min_gene_variability_pctl = min_gene_variability_pctl,
        n_prin_comps = n_prin_comps
    )

    scrublet_out <- data.table::data.table(
        cell_ID = colnames(expr_values),
        doublet_scores = scrublet_out[[1]],
        doublet = scrublet_out[[2]]
    )


    if (isTRUE(return_gobject)) {
        # Add to metadata
        gobject <- addCellMetadata(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            new_metadata = scrublet_out,
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )
        return(gobject)
    }

    return(scrublet_out)
}

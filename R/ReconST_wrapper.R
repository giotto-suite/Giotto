#' ReconST integration for Giotto
#'
#' Run gene-panel selection with \code{runReconST(gobject_sc = reference,
#' gobject_spatial = spatial)}. Export and result annotation are automatic.
#' The separate export and loading functions remain available for file-based
#' workflows, including running Python on another machine.
#'
#' @name ReconST
NULL

.reconst_use_environment <- function(envname) {
    package_check("reticulate")
    if (!is.null(envname)) {
        reticulate::use_condaenv(envname, required = TRUE)
        # Giotto object creation also consults this option before Python starts.
        options(giotto.py_path = envname)
    }
    invisible(NULL)
}

.reconst_integer <- function(x, name, lower = 0L, nullable = FALSE) {
    if (nullable && is.null(x)) {
        return(NULL)
    }
    checkmate::assert_int(x,
        lower = lower, upper = .Machine$integer.max,
        .var.name = name
    )
    as.integer(x)
}

.reconst_python <- function() {
    if (!reticulate::py_module_available("reconst")) {
        stop(
            "ReconST is not installed in the active Python environment. ",
            "Run installGiottoReconSTEnvironment() before creating objects."
        )
    }
    adapter <- new.env(parent = baseenv())
    reticulate::source_python(
        system.file("python", "python_reconst.py",
            package = "Giotto", mustWork = TRUE
        ),
        envir = adapter
    )
    adapter
}


# ---------------------------------------------------------------------------
# 1. Environment setup
# ---------------------------------------------------------------------------

#' Install a conda environment for ReconST
#'
#' Creates a reticulate-managed conda environment with ReconST and its
#' Python dependencies. Re-run with \code{force = TRUE} to rebuild.
#'
#' @param envname Name of the conda environment.
#' @param python_version Python version for the env.
#' @param reconst_version Optional PyPI version or pip installation URL.
#'   The default \code{NULL} installs a pinned revision of the public
#'   \href{https://github.com/haoranlustat/ReconST}{ReconST repository}.
#'   Use \code{package_path} for a local package.
#' @param use_gpu If \code{FALSE}, use the CPU torch index on Linux/Windows.
#'   If \code{TRUE}, use the default PyPI build; GPU support depends on the
#'   platform and installed drivers. macOS uses the PyPI CPU/MPS build.
#' @param force Re-create the env if it already exists.
#' @param package_path Local ReconST source directory or wheel. Use this to
#'   install an unpublished package; overrides \code{reconst_version}.
#'
#' @return Invisible \code{NULL}. Side effect: a conda env named
#'   \code{envname} is created and selected as the active reticulate env.
#'
#' @export
installGiottoReconSTEnvironment <- function(envname = "giotto_reconst_env",
                                            python_version = "3.10",
                                            reconst_version = NULL,
                                            use_gpu = FALSE,
                                            force = FALSE,
                                            package_path = NULL) {
    package_check("reticulate")

    package_spec <- if (!is.null(package_path)) {
        normalizePath(package_path, mustWork = TRUE)
    } else if (is.null(reconst_version)) {
        paste0(
            "git+https://github.com/haoranlustat/ReconST.git@",
            "7479571b9e579a00ff31d7b703ac0cb74ad1494c"
        )
    } else if (grepl("^(git\\+|https?://)", reconst_version)) {
        reconst_version
    } else {
        paste0("reconst==", reconst_version)
    }

    if (force && reticulate::condaenv_exists(envname)) {
        reticulate::conda_remove(envname)
    }

    if (!reticulate::condaenv_exists(envname)) {
        reticulate::conda_create(envname, python_version = python_version)
    }

    # macOS wheels (CPU/MPS) come from PyPI; the CPU index is for Linux/Windows.
    torch_options <- if (!use_gpu && Sys.info()[["sysname"]] != "Darwin") {
        c("--index-url", "https://download.pytorch.org/whl/cpu")
    } else {
        character()
    }
    reticulate::conda_install(
        envname,
        packages = "torch", pip = TRUE, pip_options = torch_options
    )
    reticulate::conda_install(
        envname,
        packages = package_spec,
        pip = TRUE
    )

    .reconst_use_environment(envname)
    invisible(NULL)
}


# ---------------------------------------------------------------------------
# 2. Giotto -> h5ad export
# ---------------------------------------------------------------------------

#' Export paired Giotto objects to h5ad for ReconST
#'
#' ReconST needs paired scRNA-seq + spatial counts. This helper extracts
#' the raw count matrix and feature IDs from each Giotto object and writes
#' them as AnnData (\code{.h5ad}) files via reticulate/anndata.
#'
#' @param gobject_sc A Giotto object holding the reference scRNA-seq data.
#' @param gobject_spatial A Giotto object holding the spatial data
#'   (e.g. MERFISH). Must share the feature ID convention with
#'   \code{gobject_sc}.
#' @param output_dir Directory to write the two h5ad files into.
#' @param spat_unit,feat_type Giotto slot identifiers; \code{NULL} uses each
#'   object's active spatial unit and feature type.
#' @param expression_values Which expression matrix to export. ReconST
#'   expects raw counts.
#' @param spat_unit_sc,feat_type_sc Reference-object slots. Default to the
#'   spatial object's slot arguments; use these when the objects differ.
#'
#' @return A named list with \code{sc_h5ad} and \code{spatial_h5ad} paths.
#'
#' @export
getReconSTInput <- function(gobject_sc,
                            gobject_spatial,
                            output_dir,
                            spat_unit = NULL,
                            feat_type = NULL,
                            expression_values = "raw",
                            spat_unit_sc = spat_unit,
                            feat_type_sc = feat_type) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ad <- reticulate::import("anndata", convert = FALSE)

    .write_one <- function(gobject, path, label, unit, type) {
        expr <- GiottoClass::getExpression(
            gobject = gobject,
            spat_unit = unit,
            feat_type = type,
            values = expression_values,
            output = "matrix"
        )
        # Giotto stores genes-by-cells; AnnData expects cells-by-genes.
        mat <- t_flex(expr)
        if (inherits(mat, "sparseMatrix")) {
            mat <- methods::as(mat, "dgCMatrix")
        } else {
            mat <- as.matrix(mat)
        }
        cells <- rownames(mat)
        genes <- colnames(mat)
        if (nrow(mat) == 0L || ncol(mat) == 0L ||
            is.null(cells) || is.null(genes) ||
            anyNA(cells) || anyNA(genes) ||
            any(!nzchar(cells)) || any(!nzchar(genes)) ||
            anyDuplicated(cells) || anyDuplicated(genes)) {
            stop(label, " expression needs nonempty, unique IDs.")
        }
        obs <- data.frame(row.names = cells, cell_id = cells)
        var <- data.frame(row.names = genes, gene_symbol = genes)

        adata <- ad$AnnData(
            X = reticulate::r_to_py(mat)$astype("float32"),
            obs = reticulate::r_to_py(obs),
            var = reticulate::r_to_py(var)
        )
        adata$write_h5ad(path)
        message(sprintf(
            "[ReconST] wrote %s: %d cells x %d genes (%s)",
            label, nrow(mat), ncol(mat), path
        ))
        path
    }

    list(
        sc_h5ad = .write_one(
            gobject_sc, file.path(output_dir, "sc.h5ad"),
            "sc", spat_unit_sc, feat_type_sc
        ),
        spatial_h5ad = .write_one(
            gobject_spatial, file.path(output_dir, "spatial.h5ad"),
            "spatial", spat_unit, feat_type
        )
    )
}


# ---------------------------------------------------------------------------
# 3. Run ReconST
# ---------------------------------------------------------------------------

#' Run the ReconST gene-panel selection pipeline
#'
#' Supply two Giotto objects to export, run ReconST, and annotate the spatial
#' object in one call. Alternatively, supply \code{input} to run existing h5ad
#' files and return a summary, as in the original file-based API.
#'
#' @param input Output of \code{\link{getReconSTInput}} (a list with
#'   \code{sc_h5ad} and \code{spatial_h5ad}), or any list/character vector
#'   with those names.
#' @param output_dir Directory ReconST writes results into. Created if
#'   missing.
#' @param embedding_size,num_epochs,batch_size Model dimensions and training
#'   schedule, forwarded to the ReconST Python adapter.
#' @param lr,weight_decay,l_lambda Optimizer and L1 penalty settings.
#' @param dropout,leaky_slope Autoencoder activation settings.
#' @param threshold Select genes whose learned weight meets this threshold.
#' @param seed Random seed. Numeric values such as \code{42} are converted
#'   to Python integers automatically.
#' @param device Optional torch device string ("cuda", "mps", "cpu").
#'   Auto-detected when \code{NULL}.
#' @param envname Conda env to activate before any Python import. Use
#'   \code{NULL} to use an already configured Python environment.
#' @param gobject_sc,gobject_spatial Reference scRNA-seq and spatial Giotto
#'   objects. Both must use the same feature IDs; cell IDs need not match.
#' @param spat_unit,feat_type Spatial-object slots.
#' @param spat_unit_sc,feat_type_sc Reference-object slots.
#' @param subset_to_panel Keep only selected features in the returned object.
#' @param train_split Fraction of reference cells used for training.
#' @param min_genes_per_cell,min_cells_per_gene Reference-data filters.
#'   Set either to \code{NULL} to skip that filter.
#' @param normalize_target_sum Reference-data normalization target; \code{NULL}
#'   skips normalization. Spatial expression is evaluated as supplied.
#' @param verbose Print training progress.
#' @param num_threads Torch CPU threads during training. Defaults to one for
#'   compatibility with R native libraries; the previous setting is restored.
#'
#' @return With Giotto objects, the spatial object with
#'   \code{reconst_importance}
#'   and \code{reconst_selected} feature metadata. With \code{input}, a list
#'   of artifact paths and summary statistics. Result files are saved in
#'   \code{output_dir} in both modes.
#'
#' @export
runReconST <- function(input = NULL,
                       output_dir = "reconst_results",
                       embedding_size = 200L,
                       num_epochs = 20L,
                       batch_size = 256L,
                       lr = 1e-3,
                       weight_decay = 1e-5,
                       l_lambda = 1e-4,
                       dropout = 0.05,
                       leaky_slope = 0.3,
                       threshold = 0.001,
                       seed = 42L,
                       device = NULL,
                       envname = "giotto_reconst_env",
                       gobject_sc = NULL,
                       gobject_spatial = NULL,
                       spat_unit = NULL,
                       feat_type = NULL,
                       spat_unit_sc = spat_unit,
                       feat_type_sc = feat_type,
                       subset_to_panel = FALSE,
                       train_split = 0.8,
                       min_genes_per_cell = 200L,
                       min_cells_per_gene = 100L,
                       normalize_target_sum = 1e4,
                       verbose = TRUE,
                       num_threads = 1L) {
    object_mode <- !is.null(gobject_sc) || !is.null(gobject_spatial)
    if (object_mode && (!is.null(input) ||
        is.null(gobject_sc) || is.null(gobject_spatial))) {
        stop("Supply both Giotto objects, or input with two h5ad paths.")
    }
    if (!object_mode && (is.null(input) ||
        !all(c("sc_h5ad", "spatial_h5ad") %in% names(input)))) {
        stop(
            "Supply both Giotto objects, or input with ",
            "sc_h5ad and spatial_h5ad."
        )
    }

    .reconst_use_environment(envname)
    reconst <- .reconst_python()
    if (object_mode) {
        io_dir <- tempfile("reconst_input_")
        dir.create(io_dir)
        on.exit(unlink(io_dir, recursive = TRUE), add = TRUE)
        input <- getReconSTInput(
            gobject_sc, gobject_spatial,
            output_dir = io_dir,
            spat_unit = spat_unit, feat_type = feat_type,
            spat_unit_sc = spat_unit_sc, feat_type_sc = feat_type_sc
        )
    }
    input <- as.list(input)
    for (name in c("sc_h5ad", "spatial_h5ad")) {
        if (length(input[[name]]) != 1L || is.na(input[[name]]) ||
            !file.exists(input[[name]])) {
            stop("Input file not found: ", name)
        }
    }

    summary <- reconst$run_reconst_pipeline(
        sc_h5ad = input$sc_h5ad,
        spatial_h5ad = input$spatial_h5ad,
        output_dir = output_dir,
        embedding_size = .reconst_integer(
            embedding_size, "embedding_size",
            lower = 1L
        ),
        num_epochs = .reconst_integer(num_epochs, "num_epochs", lower = 1L),
        batch_size = .reconst_integer(batch_size, "batch_size", lower = 1L),
        lr = lr,
        weight_decay = weight_decay,
        l_lambda = l_lambda,
        dropout = dropout,
        leaky_slope = leaky_slope,
        threshold = threshold,
        seed = .reconst_integer(seed, "seed", lower = 0L),
        device = device,
        train_split = train_split,
        min_genes_per_cell = .reconst_integer(
            min_genes_per_cell, "min_genes_per_cell",
            nullable = TRUE
        ),
        min_cells_per_gene = .reconst_integer(
            min_cells_per_gene, "min_cells_per_gene",
            nullable = TRUE
        ),
        normalize_target_sum = normalize_target_sum,
        verbose = verbose,
        num_threads = .reconst_integer(num_threads, "num_threads", lower = 1L)
    )

    if (!object_mode) {
        return(summary)
    }
    loadReconSTResults(
        gobject_spatial, summary,
        spat_unit = spat_unit, feat_type = feat_type,
        subset_to_panel = subset_to_panel
    )
}


# ---------------------------------------------------------------------------
# 4. Load results back into Giotto
# ---------------------------------------------------------------------------

#' Load ReconST results into a Giotto object
#'
#' Reads \code{reconst_results.csv} and writes \code{reconst_importance} /
#' \code{reconst_selected} columns into the spatial Giotto object's
#' feature metadata. Genes that exist in the Giotto object but are not in
#' the ReconST common-gene set are filled with \code{NA} (not 0) so that
#' "not scored" stays distinguishable from "scored and rejected".
#'
#' @param gobject_spatial The spatial Giotto object to annotate.
#' @param reconst_output Either the path to \code{reconst_results.csv} or
#'   the list returned by \code{\link{runReconST}}.
#' @param spat_unit,feat_type Giotto slot identifiers.
#' @param subset_to_panel If \code{TRUE}, subset the object to selected
#'   genes only (analogous to running with the targeted panel).
#'
#' @return The Giotto object with importance/selection annotations
#'   written into feature metadata.
#'
#' @export
loadReconSTResults <- function(gobject_spatial,
                               reconst_output,
                               spat_unit = NULL,
                               feat_type = NULL,
                               subset_to_panel = FALSE) {
    results_csv <- if (is.list(reconst_output)) {
        reconst_output$results_csv
    } else {
        reconst_output
    }
    if (is.null(results_csv) || !file.exists(results_csv)) {
        stop("Could not locate reconst_results.csv.")
    }
    # Read identifiers as text so numeric-looking IDs keep leading zeroes.
    res <- utils::read.csv(results_csv,
        colClasses = "character", check.names = FALSE
    )
    key <- if ("gene_id" %in% names(res)) "gene_id" else "gene_symbol"
    required <- c(key, "importance", "selected")
    if (!all(required %in% names(res)) || nrow(res) == 0L) {
        stop("Results need gene_id (or gene_symbol), importance, and selected.")
    }
    if (anyNA(res[[key]]) || any(!nzchar(res[[key]])) ||
        anyDuplicated(res[[key]])) {
        stop("Result feature IDs must be nonempty and unique.")
    }
    res$importance <- suppressWarnings(as.numeric(res$importance))
    res$selected <- suppressWarnings(as.logical(res$selected))
    if (any(!is.finite(res$importance)) || anyNA(res$selected)) {
        stop("Results contain invalid importance scores or selection flags.")
    }

    fmeta <- GiottoClass::getFeatureMetadata(
        gobject = gobject_spatial,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )
    if (!"feat_ID" %in% colnames(fmeta)) {
        stop("Feature metadata must contain feat_ID.")
    }
    feat_col <- "feat_ID"

    importance <- rep(NA_real_, nrow(fmeta))
    selected <- rep(NA, nrow(fmeta))
    idx <- match(fmeta[[feat_col]], res[[key]])
    has <- !is.na(idx)
    if (!any(has)) {
        stop("No result feature IDs match Giotto feat_ID; check identifiers.")
    }
    importance[has] <- res$importance[idx[has]]
    selected[has] <- as.logical(res$selected[idx[has]])
    if (isTRUE(subset_to_panel) && !any(selected, na.rm = TRUE)) {
        stop(
            "No genes were selected; lower threshold or ",
            "use subset_to_panel = FALSE."
        )
    }

    gobject_spatial <- GiottoClass::addFeatMetadata(
        gobject = gobject_spatial,
        spat_unit = spat_unit,
        feat_type = feat_type,
        new_metadata = data.frame(
            feat_ID = fmeta[[feat_col]],
            reconst_importance = importance,
            reconst_selected = selected
        ),
        by_column = TRUE,
        column_feat_ID = "feat_ID"
    )

    if (isTRUE(subset_to_panel)) {
        keep <- fmeta[[feat_col]][which(selected)]
        gobject_spatial <- GiottoClass::subsetGiotto(
            gobject = gobject_spatial,
            feat_ids = keep,
            spat_unit = spat_unit,
            feat_type = feat_type
        )
    }

    gobject_spatial
}

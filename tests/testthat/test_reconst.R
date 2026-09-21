.reconst_test_object <- function(ids = c("002", "001", "unscored"), n = 8L) {
    values <- matrix(
        rep(seq_along(ids), n),
        nrow = length(ids),
        dimnames = list(ids, paste0("cell", seq_len(n)))
    )
    createGiottoObject(
        expression = createExprObj(
            Matrix::Matrix(values, sparse = TRUE),
            name = "raw"
        ),
        spatial_locs = data.frame(
            cell_ID = colnames(values), sdimx = seq_len(n), sdimy = 0
        ),
        cores = 1L, verbose = FALSE
    )
}

.reconst_test_csv <- function(rows = NULL) {
    if (is.null(rows)) {
        rows <- data.frame(
            gene_id = c("001", "002"),
            gene_symbol = c("symbolA", "symbolB"),
            importance = c(0.2, 0.0001), selected = c(TRUE, FALSE)
        )
    }
    path <- tempfile(fileext = ".csv")
    write.csv(rows, path, row.names = FALSE)
    path
}

test_that("ReconST annotates by feature ID and preserves unscored genes", {
    g <- .reconst_test_object()
    original <- fDataDT(g)
    path <- .reconst_test_csv()
    on.exit(unlink(path))
    annotated <- loadReconSTResults(g, path)
    meta <- fDataDT(annotated)
    expect_equal(meta$feat_ID, original$feat_ID)
    expect_equal(meta$reconst_importance, c(0.0001, 0.2, NA_real_))
    expect_identical(meta$reconst_selected, c(FALSE, TRUE, NA))
    expect_false("reconst_importance" %in% names(fDataDT(g)))
    panel <- loadReconSTResults(g, path, subset_to_panel = TRUE)
    expect_equal(featIDs(panel), "001")
})

test_that("ReconST validates results and supports the older CSV schema", {
    g <- .reconst_test_object()
    path <- .reconst_test_csv(data.frame(
        gene_symbol = "001", importance = 0.3, selected = TRUE
    ))
    on.exit(unlink(path))
    meta <- fDataDT(loadReconSTResults(g, path))
    expect_true(meta$reconst_selected[meta$feat_ID == "001"])
    write.csv(
        data.frame(
            gene_id = c("001", "001"), importance = 1, selected = TRUE
        ),
        path,
        row.names = FALSE
    )
    expect_error(loadReconSTResults(g, path), "unique")
    write.csv(data.frame(gene_id = "missing", importance = 1, selected = TRUE),
        path,
        row.names = FALSE
    )
    expect_error(loadReconSTResults(g, path), "No result feature IDs")
    write.csv(data.frame(gene_id = "001", importance = "bad", selected = TRUE),
        path,
        row.names = FALSE
    )
    expect_error(loadReconSTResults(g, path), "invalid importance")
    write.csv(data.frame(gene_id = "001", importance = 0, selected = FALSE),
        path,
        row.names = FALSE
    )
    expect_error(
        loadReconSTResults(g, path, subset_to_panel = TRUE), "No genes"
    )
})

test_that("one-call ReconST annotates objects and cleans inputs on failure", {
    g <- .reconst_test_object()
    state <- new.env(parent = emptyenv())
    state$fail <- FALSE
    path <- .reconst_test_csv()
    on.exit(unlink(path))
    local_mocked_bindings(
        .reconst_use_environment = function(envname) NULL,
        .reconst_python = function() {
            list(run_reconst_pipeline = function(...) {
                args <- list(...)
                expect_type(args$seed, "integer")
                expect_type(args$num_epochs, "integer")
                expect_null(args$min_genes_per_cell)
                if (state$fail) stop("simulated training failure")
                list(results_csv = path)
            })
        },
        getReconSTInput = function(gobject_sc, gobject_spatial,
                                   output_dir, ...) {
            state$paths <- file.path(output_dir, c("sc.h5ad", "spatial.h5ad"))
            file.create(state$paths)
            setNames(as.list(state$paths), c("sc_h5ad", "spatial_h5ad"))
        }
    )
    result <- runReconST(
        gobject_sc = g, gobject_spatial = g,
        num_epochs = 1, seed = 42, min_genes_per_cell = NULL
    )
    expect_s4_class(result, "giotto")
    expect_false(any(file.exists(state$paths)))
    expect_true("reconst_selected" %in% names(fDataDT(result)))
    state$fail <- TRUE
    expect_error(runReconST(
        gobject_sc = g, gobject_spatial = g, min_genes_per_cell = NULL
    ), "simulated training failure")
    expect_false(any(file.exists(state$paths)))
    expect_error(runReconST(gobject_sc = g), "Supply both")
})

test_that("file-based ReconST accepts a named vector and retains input files", {
    paths <- c(sc_h5ad = tempfile(), spatial_h5ad = tempfile())
    file.create(paths)
    on.exit(unlink(paths))
    local_mocked_bindings(
        .reconst_use_environment = function(envname) NULL,
        .reconst_python = function() {
            list(run_reconst_pipeline = function(...) {
                args <- list(...)
                expect_equal(args$sc_h5ad, unname(paths[["sc_h5ad"]]))
                list(results_csv = "scores.csv")
            })
        }
    )
    expect_equal(runReconST(input = paths)$results_csv, "scores.csv")
    expect_true(all(file.exists(paths)))
})

test_that("ReconST installer accepts public and local sources", {
    installs <- list()
    old_options <- options(giotto.py_path = NULL)
    on.exit(options(old_options), add = TRUE)
    local_mocked_bindings(
        condaenv_exists = function(...) TRUE,
        conda_install = function(...) {
            installs[[length(installs) + 1L]] <<- list(...)
        },
        use_condaenv = function(...) NULL,
        .package = "reticulate"
    )
    installGiottoReconSTEnvironment()
    expect_identical(getOption("giotto.py_path"), "giotto_reconst_env")
    expect_identical(installs[[1]]$packages, "torch")
    expect_match(
        installs[[2]]$packages,
        "^git\\+https://github.com/haoranlustat/ReconST.git@[a-f0-9]{40}$"
    )
    installs <- list()
    installGiottoReconSTEnvironment(package_path = tempdir())
    expect_identical(installs[[2]]$packages, normalizePath(tempdir()))
})

test_that("ReconST exports sparse h5ad and trains with real Python", {
    skip_if(
        Sys.getenv("GIOTTO_TEST_RECONST") != "true",
        "Set GIOTTO_TEST_RECONST=true in an environment with ReconST installed."
    )
    skip_if_not(reticulate::py_module_available("reconst"))
    sc <- .reconst_test_object(c("001", "002", "003", "004"), n = 12L)
    spatial <- .reconst_test_object(c("004", "001", "002", "unscored"), n = 6L)
    folder <- tempfile("reconst_integration_")
    dir.create(folder)
    on.exit(unlink(folder, recursive = TRUE))
    inputs <- getReconSTInput(sc, spatial, file.path(folder, "inputs"))
    ad <- reticulate::import("anndata", convert = FALSE)
    scipy <- reticulate::import("scipy.sparse", convert = TRUE)
    exported <- ad$read_h5ad(inputs$sc_h5ad)
    expect_true(scipy$issparse(exported$X))
    expect_equal(as.integer(reticulate::py_to_r(exported$shape)), c(12L, 4L))
    expect_equal(
        reticulate::py_to_r(exported$var_names$to_list()),
        c("001", "002", "003", "004")
    )
    annotated <- runReconST(
        gobject_sc = sc, gobject_spatial = spatial,
        output_dir = folder, envname = NULL, device = "cpu",
        num_epochs = 1, batch_size = 4, embedding_size = 2, seed = 42,
        min_genes_per_cell = NULL, min_cells_per_gene = NULL,
        normalize_target_sum = NULL, verbose = FALSE
    )
    meta <- fDataDT(annotated)
    scores <- read.csv(file.path(folder, "reconst_results.csv"),
        colClasses = "character"
    )
    expect_setequal(scores$gene_id, c("001", "002", "004"))
    expect_equal(
        meta$reconst_importance,
        as.numeric(scores$importance[match(meta$feat_ID, scores$gene_id)])
    )
    expect_true(is.na(meta$reconst_selected[meta$feat_ID == "unscored"]))
    expect_true(all(file.exists(file.path(folder, c(
        "autoencoder_model.pth", "training_history.csv",
        "evaluation_metrics.json"
    )))))
})


test_that("ReconST integer parameters reject fractional or missing values", {
    expect_error(.reconst_integer(1.5, "num_epochs", lower = 1L), "num_epochs")
    expect_error(.reconst_integer(NA_real_, "seed"), "seed")
    expect_identical(.reconst_integer(42, "seed"), 42L)
    expect_null(.reconst_integer(NULL, "min_cells", nullable = TRUE))
})

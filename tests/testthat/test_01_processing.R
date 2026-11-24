
# DATA TO USE
# visium mini expression
g <- GiottoData::loadGiottoMini("visium")


g <- test_data$vis

# normalize ####

test_that("expression is normalized", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    # remove normalized and scaled expr matrices
    g <- setExpression(g, NULL, spat_unit = "cell", feat_type = "rna", name = "normalized", verbose = FALSE)
    g <- setExpression(g, NULL, spat_unit = "cell", feat_type = "rna", name = "scaled", verbose = FALSE)

    expect_false(any(c("normalized", "scaled") %in% list_expression(g)$name))

    g <- normalizeGiotto(g, verbose = FALSE)

    expect_true(any(c("normalized", "scaled") %in% list_expression(g)$name))
})



# hvf ####

test_that("highly variable gene detection", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    # reset feature metadata
    g <- setFeatureMetadata(g, NULL, spat_unit = "cell", feat_type = "rna", verbose = FALSE)
    expect_false("hvf" %in% names(fDataDT(g)))

    g <- calculateHVF(g)
    expect_true("hvf" %in% names(fDataDT(g)))

    # character "yes" and "no" expected
    checkmate::expect_character(fDataDT(g)$hvf)
})

test_that("highly variable gene detections - pearson resid", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    # reset feature metadata
    g <- setFeatureMetadata(g, NULL, spat_unit = "cell", feat_type = "rna", verbose = FALSE)
    expect_false(any(c("var", "hvf") %in% names(fDataDT(g))))

    g <- normalizeGiotto(
        gobject = g,
        feat_type = "rna",
        scalefactor = 5000,
        verbose = FALSE,
        norm_methods = "pearson_resid",
        update_slot = "pearson"
    )

    # var_p_resid builds the residuals itself from raw counts, so the
    # pearson_resid normalization above is no longer a prerequisite and
    # expression_values is not consulted.
    g <- calculateHVF(
        g,
        method = "var_p_resid",
        var_threshold = 7
    )

    expect_true(any(c("var", "hvf") %in% names(fDataDT(g))))

    # character "yes" and "no" expected
    checkmate::expect_numeric(fDataDT(g)$var)
    checkmate::expect_character(fDataDT(g)$hvf)
})

test_that("var_p_resid warns that expression_values is ignored", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    g <- setFeatureMetadata(g, NULL, spat_unit = "cell", feat_type = "rna",
                            verbose = FALSE)
    expect_warning(
        calculateHVF(g, method = "var_p_resid", expression_values = "normalized"),
        "ignored for method"
    )
})

test_that("highly variable features by plain variance", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    g <- setFeatureMetadata(g, NULL, spat_unit = "cell", feat_type = "rna",
        verbose = FALSE)
    g <- calculateHVF(g, method = "var", n_top_feats = 50,
        show_plot = FALSE, return_plot = FALSE, save_plot = FALSE)
    fm <- fDataDT(g)
    x <- getExpression(g, values = "normalized", output = "matrix")
    v <- apply(as.matrix(x), 1, stats::var)
    top <- names(sort(v[v > 0], decreasing = TRUE))[1:50]
    expect_setequal(fm[hvf == "yes"]$feat_ID, top)
    expect_equal(fm$var[match(top, fm$feat_ID)], unname(v[top]),
        tolerance = 1e-8)
})

test_that("var_p_resid row variance matches base variance on Giotto matrices", {
    rlang::local_options(lifecycle_verbosity = "quiet")

    g_var <- normalizeGiotto(
        gobject = g,
        feat_type = "rna",
        scalefactor = 5000,
        verbose = FALSE,
        norm_methods = "pearson_resid",
        update_slot = "pearson"
    )

    pearson_mat <- getExpression(
        g_var,
        spat_unit = "cell",
        feat_type = "rna",
        values = "pearson",
        output = "matrix"
    )

    expect_equal(
        .rowVars_flex(pearson_mat),
        apply(pearson_mat, 1, stats::var),
        tolerance = 1e-12
    )
})


# statistics ####

test_that("statistics are added", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    # reset cell metadata
    g <- setCellMetadata(g, NULL, spat_unit = "cell", feat_type = "rna", verbose = FALSE)
    expect_false(any(c("nr_feats", "perc_feats", "total_expr") %in% names(pDataDT(g))))

    g <- addStatistics(g)
    expect_true(any(c("nr_feats", "perc_feats", "total_expr") %in% names(pDataDT(g))))
})

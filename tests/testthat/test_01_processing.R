
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

    g <- calculateHVF(
        g,
        method = "var_p_resid",
        var_threshold = 7,
        expression_values = "pearson"
    )

    expect_true(any(c("var", "hvf") %in% names(fDataDT(g))))

    # character "yes" and "no" expected
    checkmate::expect_numeric(fDataDT(g)$var)
    checkmate::expect_character(fDataDT(g)$hvf)
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

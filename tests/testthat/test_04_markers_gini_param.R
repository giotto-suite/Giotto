# analyzeData(x, markersParam(method = "gini")) -- the verb the gini machinery
# now lives in. findGiniMarkers() and findGiniMarkers_one_vs_all() are thin
# gobject wrappers over it, so these cover the statistic and those cover the
# plumbing.

g <- test_data$vis
EX <- getExpression(g, values = "normalized", output = "matrix")
CM <- getCellMetadata(g, output = "data.table")
GRP <- stats::setNames(as.character(CM$leiden_clus), CM$cell_ID)

.gp <- function(...) markersParam(method = "gini", ...)


test_that("markersParam(method = 'gini') builds the documented defaults", {
    p <- .gp()
    expect_s4_class(p, "giniMarkersParam")
    expect_true(is(p, "markersParam"))
    expect_identical(p$comparison, "pairwise")
    expect_identical(p$min_expression, 0.2)
    expect_identical(p$min_detection, 0.2)
    expect_identical(p$min_expression_gini, -Inf)
    expect_identical(p$min_detection_gini, -Inf)
    expect_identical(p$detection_threshold, 0)
    expect_identical(p$min_length, 0)
    expect_identical(p$rank_score, Inf)
    expect_identical(p$min_feats, 5)

    # explicit values win
    expect_identical(.gp(min_feats = 2, rank_score = 3)$min_feats, 2)
    expect_identical(.gp(min_feats = 2, rank_score = 3)$rank_score, 3)
})

test_that("the factory still routes scran, and rejects unknown methods", {
    expect_s4_class(markersParam(method = "scran"), "scranMarkersParam")
    expect_error(markersParam(method = "nope"))
})


# the statistic ####

test_that("the verb reproduces findGiniMarkers exactly", {
    # findGiniMarkers()'s defaults are the param's defaults, so this is a
    # like-for-like comparison with nothing passed on either side.
    expect_identical(
        analyzeData(EX, .gp(), groups = GRP),
        findGiniMarkers(g,
            cluster_column = "leiden_clus", expression_values = "normalized"
        )
    )
})

test_that("the verb reproduces findGiniMarkers_one_vs_all exactly", {
    # that wrapper's defaults differ from the param's -- 0.5/0.5/4 rather than
    # 0.2/0.2/5 -- and `comparison` deliberately does not switch them, so they
    # are passed explicitly here.
    got <- analyzeData(EX,
        .gp(comparison = "one_vs_rest", min_expression = 0.5,
            min_detection = 0.5, min_feats = 4),
        groups = GRP, verbose = FALSE
    )
    ref <- findGiniMarkers_one_vs_all(g,
        cluster_column = "leiden_clus", expression_values = "normalized",
        verbose = FALSE
    )
    expect_identical(got, ref)
})

test_that("gates and knobs reach the statistic through the param", {
    loose <- analyzeData(EX, .gp(min_expression = 0, min_detection = 0),
        groups = GRP)
    tight <- analyzeData(EX, .gp(min_expression = 5, min_detection = 0.9),
        groups = GRP)
    expect_gt(nrow(loose), nrow(tight))

    # gates are OR'd with min_feats, so tightening them cannot go below
    # min_feats per group
    n_groups <- length(unique(loose$cluster))
    expect_gte(nrow(tight), .gp()$min_feats * n_groups)

    # min_length changes the coefficients, being a different statistic
    padded <- analyzeData(EX, .gp(min_length = 20), groups = GRP)
    expect_false(isTRUE(all.equal(
        padded$expression_gini,
        analyzeData(EX, .gp(), groups = GRP)$expression_gini
    )))
})

test_that("the column contract is stated, not incidental", {
    expect_named(analyzeData(EX, .gp(), groups = GRP), c(
        "feats", "cluster", "expression", "expression_gini",
        "detection", "detection_gini", "expression_rank",
        "detection_rank", "comb_score", "comb_rank"
    ))
})


# contract ####

test_that("groups is required", {
    expect_error(analyzeData(EX, .gp()), "`groups` is required")
})

test_that("NA groups narrow without subsetting the object", {
    keep <- unique(GRP)[1:2]
    narrowed <- GRP
    narrowed[!narrowed %in% keep] <- NA
    got <- analyzeData(EX, .gp(), groups = narrowed)
    expect_setequal(unique(got$cluster), keep)
})

test_that("a named grouping is matched on cell ID, not position", {
    expect_false(identical(names(GRP), colnames(EX)))
    set.seed(7L)
    shuffled <- GRP[sample(length(GRP))]
    expect_identical(
        analyzeData(EX, .gp(), groups = shuffled),
        analyzeData(EX, .gp(), groups = GRP)
    )
})

test_that("an x the featStats verb cannot handle is refused directly", {
    # The whole method is derived from analyzeData(x, featStatsParam), so the
    # error names that rather than failing somewhere inside the statistic.
    expect_error(
        analyzeData(data.frame(a = 1), .gp(), groups = GRP),
        "featStatsParam"
    )
})

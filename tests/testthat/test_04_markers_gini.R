
# DATA TO USE
g <- test_data$vis
CLUS <- "leiden_clus"

.gini <- function(...) {
    findGiniMarkers(
        gobject = g,
        expression_values = "normalized",
        cluster_column = CLUS,
        min_feats = 5,
        ...
    )
}

# gate naming ####

test_that("min_expression and min_detection gate the value columns", {
    # The floors apply to `expression` and `detection`, not to the gini
    # coefficients -- which is what the old parameter names implied and never
    # did. Independent ground truth: filter the unfiltered result by hand.
    all_rows <- .gini(min_expression = -Inf, min_detection = -Inf)
    gated <- .gini(min_expression = 1, min_detection = 0.5)

    manual <- all_rows[
        comb_rank <= 5 | (expression > 1 & detection > 0.5)
    ]
    expect_equal(nrow(gated), nrow(manual))
    expect_equal(gated$feats, manual$feats)

    # and they are genuinely restrictive on this data
    expect_lt(nrow(gated), nrow(all_rows))
})

test_that("min_expression and min_detection are not interchangeable", {
    # A feature expressed intensely in a few cells has high `expression` and
    # low `detection`; gating one is not the same as gating the other.
    e_only <- .gini(min_expression = 1, min_detection = -Inf)
    d_only <- .gini(min_expression = -Inf, min_detection = 0.5)
    expect_false(identical(nrow(e_only), nrow(d_only)))
})


# gini gates ####

test_that("min_expression_gini and min_detection_gini are inert by default", {
    expect_equal(.gini(), .gini(
        min_expression_gini = -Inf, min_detection_gini = -Inf
    ))
})

test_that("the gini gates apply to the gini columns", {
    all_rows <- .gini(min_expression = -Inf, min_detection = -Inf)
    gated <- .gini(
        min_expression = -Inf, min_detection = -Inf,
        min_expression_gini = 0.1, min_detection_gini = 0.1
    )

    manual <- all_rows[
        comb_rank <= 5 | (expression_gini > 0.1 & detection_gini > 0.1)
    ]
    expect_equal(nrow(gated), nrow(manual))
    expect_equal(gated$feats, manual$feats)
})

test_that("the gini gates are monotonic and bottom out at the min_feats floor", {
    counts <- vapply(c(0.1, 0.2, 0.3), function(thr) {
        nrow(.gini(min_expression_gini = thr, min_detection_gini = thr))
    }, numeric(1L))

    expect_true(all(diff(counts) < 0))

    # `min_feats` rescues the top 5 per cluster no matter how strict the gate,
    # so the count can never fall below min_feats x n_clusters.
    n_clus <- length(unique(getCellMetadata(g, output = "data.table")[[CLUS]]))
    expect_gte(min(counts), 5L * n_clus)
})


# deprecation ####

test_that("deprecated gate names forward to the renamed arguments", {
    rlang::local_options(lifecycle_verbosity = "quiet")

    expect_equal(
        .gini(min_expr_gini_score = 1, min_det_gini_score = 0.5),
        .gini(min_expression = 1, min_detection = 0.5)
    )
})

test_that("deprecated gate names warn", {
    rlang::local_options(lifecycle_verbosity = "warning")

    expect_warning(
        .gini(min_expr_gini_score = 0.2), class = "lifecycle_warning_deprecated"
    )
    expect_warning(
        .gini(min_det_gini_score = 0.2), class = "lifecycle_warning_deprecated"
    )
})


# dispatchers forward the gates ####

test_that("findMarkers forwards the gini gates", {
    # `findMarkers()` defaults the floors to 0.5 where `findGiniMarkers()` uses
    # 0.2, so set them explicitly rather than relying on the defaults agreeing.
    direct <- .gini(
        min_expression = 0.2, min_detection = 0.2,
        min_expression_gini = 0.15, min_detection_gini = 0.15
    )
    viad <- findMarkers(g,
        method = "gini", expression_values = "normalized",
        cluster_column = CLUS, min_feats = 5,
        min_expression = 0.2, min_detection = 0.2,
        min_expression_gini = 0.15, min_detection_gini = 0.15
    )
    expect_equal(direct, viad)
})

test_that("findMarkers_one_vs_all forwards the gini gates", {
    loose <- findMarkers_one_vs_all(g,
        method = "gini", expression_values = "normalized",
        cluster_column = CLUS, min_feats = 5, verbose = FALSE
    )
    tight <- findMarkers_one_vs_all(g,
        method = "gini", expression_values = "normalized",
        cluster_column = CLUS, min_feats = 5, verbose = FALSE,
        min_expression_gini = 0.2, min_detection_gini = 0.2
    )
    expect_lt(nrow(tight), nrow(loose))
})


# output contract ####

test_that("the returned columns are unchanged", {
    expect_named(.gini(), c(
        "feats", "cluster", "expression", "expression_gini",
        "detection", "detection_gini", "expression_rank",
        "detection_rank", "comb_score", "comb_rank"
    ))
})

test_that("findMarkers_one_vs_all keeps the same column contract", {
    # Mirrors the assertion in test_99_merFISH.R, which only runs when the
    # merFISH dataset can be downloaded. Running it here keeps the contract
    # guarded on every ordinary test run.
    expect_named(
        findMarkers_one_vs_all(g,
            method = "gini", expression_values = "normalized",
            cluster_column = CLUS, min_feats = 1, rank_score = 2,
            verbose = FALSE
        ),
        c(
            "feats", "cluster", "expression", "expression_gini",
            "detection", "detection_gini", "expression_rank",
            "detection_rank", "comb_score", "comb_rank"
        )
    )
})


# regressions ####

test_that("the default gates do not collapse the result to the min_feats floor", {
    # Regression for #1238: pointing the default floors at the gini columns
    # instead of the value columns left nothing able to pass them, because the
    # 0.5 default sits at or above the gini ceiling. The filter went vacuous
    # and the output fell to exactly min_feats per cluster. Anything that
    # re-aims these gates will reproduce that, so assert the result is
    # meaningfully larger than the rescue floor.
    res <- findMarkers_one_vs_all(g,
        method = "gini", expression_values = "normalized",
        cluster_column = CLUS, min_feats = 5, verbose = FALSE
    )
    n_clus <- length(unique(getCellMetadata(g, output = "data.table")[[CLUS]]))
    expect_gt(nrow(res), 5L * n_clus)
})

test_that("a gini coefficient over G values cannot exceed (G - 1) / G", {
    # The documented ceiling. It is what makes a gini threshold non-
    # transferable between runs with different cluster counts, and why the
    # 0.5 default could never pass in a two-group comparison.
    for (G in c(2L, 3L, 5L, 7L, 20L)) {
        expect_equal(
            mygini_fun(c(1, rep(0, G - 1L))), (G - 1L) / G,
            tolerance = 1e-9
        )
    }

    # It must also hold for negative values, since `expression_values =
    # "scaled"` is supported and scaled data is signed.
    set.seed(42)
    for (G in c(2L, 5L, 7L)) {
        vals <- vapply(seq_len(300L), function(i) {
            mygini_fun(stats::rnorm(G, mean = 0, sd = 10))
        }, numeric(1L))
        vals <- vals[is.finite(vals)]
        expect_lte(max(vals), (G - 1L) / G + 1e-9)
    }
})

test_that("one_vs_all gini coefficients are capped at 0.5", {
    # Each iteration compares one cluster against the rest pooled, so G = 2.
    res <- findMarkers_one_vs_all(g,
        method = "gini", expression_values = "normalized",
        cluster_column = CLUS, min_feats = 5, verbose = FALSE
    )
    expect_lte(max(res$expression_gini), 0.5 + 1e-9)
    expect_lte(max(res$detection_gini), 0.5 + 1e-9)
})

test_that("no combination of gates can push the result below min_feats", {
    # All four gates are OR'd with the `min_feats` rescue, so the result
    # shrinks towards min_feats per cluster and never below it.
    starved <- .gini(
        min_expression = Inf, min_detection = Inf,
        min_expression_gini = Inf, min_detection_gini = Inf
    )
    n_clus <- length(unique(getCellMetadata(g, output = "data.table")[[CLUS]]))
    expect_equal(nrow(starved), 5L * n_clus)
    expect_true(all(starved$comb_rank <= 5L))
})

test_that("the gate comparisons are strict", {
    # A value exactly equal to its threshold does not pass. Detection values
    # are fractions over a cluster's cell count, so exact hits are common.
    all_rows <- .gini(min_expression = -Inf, min_detection = -Inf)
    beyond <- all_rows[comb_rank > 5L] # outside the min_feats rescue
    thr <- beyond$detection[[1L]]

    gated <- .gini(min_expression = -Inf, min_detection = thr)
    expect_false(any(gated[comb_rank > 5L]$detection == thr))
    expect_true(any(beyond$detection == thr))
})

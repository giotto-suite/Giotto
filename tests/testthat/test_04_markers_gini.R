
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


# min_length ####

test_that("min_length is inert at its 0 default", {
    expect_equal(.gini(), .gini(min_length = 0))
})

test_that("min_length lifts the (G - 1) / G ceiling", {
    # Padding the per-cluster vector to n entries raises the attainable
    # coefficient to (n - 1) / n regardless of the cluster count, which is what
    # makes gini values comparable between runs.
    expect_equal(
        mygini_fun(c(1, 0), min_length = 16), 15 / 16,
        tolerance = 1e-9
    )
    expect_equal(
        mygini_fun(c(1, rep(0, 4)), min_length = 16), 15 / 16,
        tolerance = 1e-9
    )
    # never shortens: a vector already longer than min_length is untouched
    v <- c(1, rep(0, 19))
    expect_equal(mygini_fun(v, min_length = 16), mygini_fun(v))
})

test_that("min_length reaches the gini scores through every entry point", {
    # Regression guard of the D1 kind: a parameter that exists but is never
    # forwarded. one_vs_all matters most here, since it always compares two
    # groups and is therefore capped at 0.5 without padding.
    base_pw <- .gini()
    pad_pw <- .gini(min_length = 16)
    expect_gt(max(pad_pw$expression_gini), max(base_pw$expression_gini))

    ova <- function(...) {
        findMarkers_one_vs_all(g,
            method = "gini", expression_values = "normalized",
            cluster_column = CLUS, min_feats = 5, verbose = FALSE, ...
        )
    }
    expect_lte(max(ova()$expression_gini), 0.5 + 1e-9)
    expect_gt(max(ova(min_length = 16)$expression_gini), 0.5)

    direct <- findGiniMarkers_one_vs_all(g,
        cluster_column = CLUS, expression_values = "normalized",
        min_feats = 5, verbose = FALSE, min_length = 16
    )
    expect_gt(max(direct$expression_gini), 0.5)
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

test_that("rank_score is inert at its Inf default", {
    expect_equal(.gini(), .gini(rank_score = Inf))
})

test_that("rank_score filters on the raw rank, not the rescaled weight", {
    # Regression: the ranks used to be overwritten in place by their own
    # rescaling into [1, 0.1], so `expression_rank <= rank_score` could never
    # fail for any rank_score >= 1 and the argument did nothing at its own
    # default. Anything that reintroduces that will make these equal.
    expect_lt(nrow(.gini(rank_score = 1)), nrow(.gini(rank_score = Inf)))
    expect_lt(nrow(.gini(rank_score = 1)), nrow(.gini(rank_score = 2)))

    # and the returned ranks are positions, not weights
    r <- .gini()
    expect_true(all(r$expression_rank >= 1))
    expect_true(all(r$expression_rank == as.integer(r$expression_rank)))
})

test_that("rank_score = 1 keeps only clusters topping the feature", {
    kept <- .gini(rank_score = 1, min_expression = -Inf, min_detection = -Inf)
    beyond <- kept[comb_rank > 5L] # outside the min_feats rescue
    expect_true(all(beyond$expression_rank == 1L))
    expect_true(all(beyond$detection_rank == 1L))
})

test_that("clusters tied at the top all count as rank 1", {
    # With ties.method = "average" a two-way tie gives both clusters 1.5, so
    # `<= 1` rejected the feature entirely instead of crediting each winner.
    # Detection is a fraction over cluster size, so top ties are common.
    r <- .gini(min_expression = -Inf, min_detection = -Inf)
    tied <- r[, .(n_top = sum(detection_rank == 1L)), by = feats]
    expect_true(any(tied$n_top > 1L))
    # no feature may end up with no rank-1 cluster at all
    expect_true(all(tied$n_top >= 1L))
})

test_that("in one_vs_all only rank_score = 1 can bind", {
    # Each iteration compares two groups, so ranks are only ever 1 or 2 and
    # any threshold above 1 admits everything. Documented in @details.
    f <- function(rs) {
        findGiniMarkers_one_vs_all(g,
            cluster_column = CLUS, expression_values = "normalized",
            min_feats = 1, verbose = FALSE, rank_score = rs
        )
    }
    loose <- f(Inf)
    expect_equal(f(2), loose)
    expect_lt(nrow(f(1)), nrow(loose))
    expect_lte(max(loose$expression_rank), 2L)
    expect_lte(max(loose$detection_rank), 2L)
})

test_that("findMarkers_one_vs_all forwards rank_score", {
    # Regression: the gini branch declared rank_score but never passed it on.
    loose <- findMarkers_one_vs_all(g,
        method = "gini", expression_values = "normalized",
        cluster_column = CLUS, min_feats = 1, verbose = FALSE,
        rank_score = Inf
    )
    tight <- findMarkers_one_vs_all(g,
        method = "gini", expression_values = "normalized",
        cluster_column = CLUS, min_feats = 1, verbose = FALSE,
        rank_score = 1
    )
    expect_lt(nrow(tight), nrow(loose))
})

test_that("comb_rank still follows comb_score within each cluster", {
    # `comb_score` is built from the rescaled weights, which kept the original
    # default tie handling when the rank was split out; only the returned rank
    # columns changed. So the score -> rank relationship must be intact.
    r <- .gini(min_expression = -Inf, min_detection = -Inf)
    chk <- r[, .(ok = !is.unsorted(comb_rank[order(-comb_score)])),
        by = cluster
    ]
    expect_true(all(chk$ok))
    expect_false(any(is.na(r$comb_score)))
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

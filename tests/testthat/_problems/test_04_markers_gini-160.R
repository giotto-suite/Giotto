# Extracted from test_04_markers_gini.R:160

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
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

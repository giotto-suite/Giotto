# Extracted from test_04_markers_gini.R:142

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
expect_named(.gini(), c(
        "feats", "cluster", "expression", "expression_gini",
        "detection", "detection_gini", "expression_rank",
        "detection_rank", "comb_score", "comb_rank"
    ))

# Extracted from test_04_feat_stats_grouped.R:90

# prequel ----------------------------------------------------------------------
g <- test_data$vis
EX <- getExpression(g, values = "normalized", output = "matrix")
CL <- getCellMetadata(g, output = "data.table")$leiden_clus
LV <- levels(droplevels(factor(CL)))
.fs <- function(...) {
    analyzeData(EX, analyzeParam("feat_stats"), ...)
}

# test -------------------------------------------------------------------------
expect_error(.fs(groups = CL[-1L]), "one entry per cell")

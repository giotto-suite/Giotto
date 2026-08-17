# Extracted from test_04_feat_stats_grouped.R:75

# prequel ----------------------------------------------------------------------
g <- test_data$vis
EX <- getExpression(g, values = "normalized", output = "matrix")
CL <- getCellMetadata(g, output = "data.table")$leiden_clus
LV <- levels(droplevels(factor(CL)))
.fs <- function(...) {
    analyzeData(EX, analyzeParam("feat_stats"), ...)
}

# test -------------------------------------------------------------------------
lo <- .fs(groups = CL, stats = c("sum", "nnz"))
hi <- analyzeData(EX, analyzeParam("feat_stats", detection_threshold = 1),
        groups = CL, stats = c("sum", "nnz"))
expect_equal(lo$total_expr, hi$total_expr)
expect_true(sum(hi$nr_cells) < sum(lo$nr_cells))

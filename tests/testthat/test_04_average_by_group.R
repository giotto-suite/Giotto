
# .average_by_group() is a temporary local copy of the GiottoClass helpers,
# carried until analyzeData(x, featStatsParam) reaches this branch. These
# tests pin the property that made the originals wrong: group membership must
# be keyed on cell_ID, not on position.

g <- test_data$vis
EX <- getExpression(g, values = "normalized", output = "matrix")
MD <- getCellMetadata(g, output = "data.table")
CLUS <- "leiden_clus"

.ref_by_id <- function(fun) {
    lv <- unique(MD[[CLUS]])
    out <- vapply(lv, function(k) {
        ids <- MD[get(CLUS) == k][["cell_ID"]]
        fun(EX[, colnames(EX) %in% ids, drop = FALSE])
    }, numeric(nrow(EX)))
    colnames(out) <- paste0("cluster_", lv)
    out
}

test_that("group means select cells by identifier, not position", {
    got <- as.matrix(Giotto:::.average_by_group(
        g, NULL, NULL, CLUS, "normalized"
    ))
    expect_equal(got, .ref_by_id(Matrix::rowMeans), ignore_attr = TRUE)
})

test_that("detection fractions select cells by identifier", {
    got <- as.matrix(Giotto:::.average_by_group(
        g, NULL, NULL, CLUS, "normalized",
        detection_threshold = 0
    ))
    ref <- .ref_by_id(function(m) Matrix::rowSums(m > 0) / ncol(m))
    expect_equal(got, ref, ignore_attr = TRUE)
})

test_that("results are invariant to cell metadata row order", {
    # holds regardless of any fixture's incidental ordering
    cm <- getCellMetadata(g, output = "cellMetaObj")
    set.seed(9)
    cm[] <- cm[][sample(nrow(cm[]))]
    g2 <- setGiotto(g, cm, verbose = FALSE)

    a <- Giotto:::.average_by_group(g, NULL, NULL, CLUS, "normalized")
    b <- Giotto:::.average_by_group(g2, NULL, NULL, CLUS, "normalized")
    expect_equal(a[, sort(colnames(a))], b[, sort(colnames(b))])
})

test_that("the column contract downstream code relies on is preserved", {
    got <- Giotto:::.average_by_group(g, NULL, NULL, CLUS, "normalized")
    expect_true(all(grepl("^cluster_", colnames(got))))
    expect_equal(rownames(got), rownames(EX))
    # column order follows first appearance in the metadata, as before
    expect_equal(colnames(got), paste0("cluster_", unique(MD[[CLUS]])))
})

test_that("a missing metadata column errors", {
    expect_error(
        Giotto:::.average_by_group(g, NULL, NULL, "not_a_column", "normalized"),
        "not found"
    )
})

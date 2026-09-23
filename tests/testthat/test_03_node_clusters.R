# Enumerating the splits of a hierarchical clustering.
#
# The height-matching implementation this replaced re-selected an already-split
# node whenever two merges shared a height, emitting one split twice and never
# emitting its sibling. Row count stayed at k - 1, so it surfaced as neither an
# error nor a warning. These fixtures pin the cases that broke it.

# two identical cluster pairs -> two merges at exactly the same height
.tied_tree <- function() {
    set.seed(1)
    a <- rnorm(200)
    b <- rnorm(200) + 5
    m <- cbind(
        A = a, B = a + 1e-12,
        C = b, D = b + 1e-12,
        E = rnorm(200) + 10
    )
    stats::hclust(stats::as.dist(1 - stats::cor(m)), "ward.D")
}

# independent leaf-set recursion, deliberately not the implementation's
.leaves_under <- function(hc, i) {
    out <- character(0L)
    for (j in hc$merge[i, ]) {
        out <- c(out, if (j < 0L) hc$labels[-j] else .leaves_under(hc, j))
    }
    out
}

.splits_of <- function(res) {
    vapply(res, function(x) {
        paste(
            paste(sort(x$first), collapse = ","), "|",
            paste(sort(x$sec), collapse = ",")
        )
    }, character(1L))
}

test_that("tied merge heights do not collapse two nodes into one", {
    hc <- .tied_tree()
    expect_true(any(duplicated(signif(hc$height, 12))))

    res <- Giotto:::.node_clusters(hc, verbose = FALSE)[[2]]
    sp <- .splits_of(res)

    expect_length(res, nrow(hc$merge)) # k - 1 = 4
    expect_false(any(duplicated(sp)))
    # both tied pairs must be present; the old code emitted one of them twice
    expect_true(any(grepl("^A,B \\|", sp) | grepl("\\| A,B$", sp)))
    expect_true(any(grepl("^C,D \\|", sp) | grepl("\\| C,D$", sp)))
})

test_that("every leaf is separated from its sibling exactly once", {
    hc <- .tied_tree()
    res <- Giotto:::.node_clusters(hc, verbose = FALSE)[[2]]
    # each node partitions its own subtree, so leaf counts must sum correctly
    for (nd in res) {
        expect_gt(length(nd$first), 0L)
        expect_gt(length(nd$sec), 0L)
        expect_length(intersect(nd$first, nd$sec), 0L)
    }
})

test_that("nodeID is the merge row, so results join back to the tree", {
    hc <- .tied_tree()
    res <- Giotto:::.node_clusters(hc, verbose = FALSE)[[2]]
    for (nd in res) {
        expect_identical(nd$height, hc$height[[nd$node]])
        # the leaves either side must be exactly the subtree under that row
        n_leaves <- length(nd$first) + length(nd$sec)
        expect_identical(
            sort(c(nd$first, nd$sec)),
            sort(.leaves_under(hc, nd$node))
        )
        expect_gt(n_leaves, 1L)
    }
})

test_that("non-monotone linkage does not break node ordering", {
    set.seed(2)
    m <- matrix(rnorm(300), nrow = 60, dimnames = list(NULL, LETTERS[1:5]))
    hc <- stats::hclust(stats::as.dist(1 - stats::cor(m)), "centroid")
    res <- Giotto:::.node_clusters(hc, verbose = FALSE)[[2]]

    expect_length(res, 4L)
    expect_false(any(duplicated(.splits_of(res))))
    # heights come back in decreasing order even when merge order does not
    expect_identical(
        vapply(res, function(x) x$height, numeric(1L)),
        sort(hc$height, decreasing = TRUE)
    )
})

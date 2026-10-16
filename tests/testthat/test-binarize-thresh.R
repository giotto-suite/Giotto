# Per-feature thresholds for thresholdParam("binarize"). Ground truth: a
# brute-force 2-means over every split point, and the rank rule
# `rank(-x, ties = "average") <= N * p / 100`.

.brute_2means <- function(x) {
    u <- sort(unique(x))
    sse <- vapply(u[-length(u)], function(t) {
        lo <- x[x <= t]; hi <- x[x > t]
        sum((lo - mean(lo))^2) + sum((hi - mean(hi))^2)
    }, numeric(1))
    u[which.min(sse)]
}
.rank_on <- function(x, p) rank(-x, ties.method = "average") <= length(x) / 100 * p

.bin_mat <- function(G = 12L, N = 60L, seed = 4L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(G, N, density = 0.45,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    dimnames(m) <- list(paste0("g", seq_len(G)), paste0("c", seq_len(N)))
    m
}

test_that("histThreshold kmeans is the brute-force 2-means split", {
    set.seed(1)
    for (i in 1:20) {
        x <- c(rep(0, sample(5:30, 1)), round(rexp(sample(5:40, 1), 0.3), 1))
        r <- rle(sort(x))
        expect_equal(histThreshold(r$values, r$lengths, method = "kmeans"),
                     .brute_2means(x))
    }
    v <- c(0, 1, 2, 5, 6, 9)
    n <- c(20, 4, 3, 2, 5, 1)
    expect_equal(histThreshold(v, n, upper = v + 0.5),
                 histThreshold(v, n) + 0.5)
    expect_identical(histThreshold(5, 3), Inf)
})

test_that("threshold = \"kmeans\" splits each feature at its 2-means optimum", {
    m <- .bin_mat()
    dm <- as.matrix(m)
    ref <- t(apply(dm, 1, function(x) x > .brute_2means(x))) * 1
    p <- thresholdParam("binarize", threshold = "kmeans")
    expect_equal(as.matrix(processData(m, p)), ref, ignore_attr = TRUE)
    expect_equal(processData(dm, p), ref, ignore_attr = TRUE)
})

test_that("threshold = \"rank\" reproduces the average-rank rule", {
    m <- .bin_mat()
    dm <- as.matrix(m)
    for (pr in c(10, 30)) {
        p <- thresholdParam("binarize", threshold = "rank", percentage_rank = pr)
        ref <- t(apply(dm, 1, .rank_on, p = pr)) * 1
        expect_equal(as.matrix(processData(m, p)), ref, ignore_attr = TRUE)
    }
})

test_that("a numeric threshold may be one per feature", {
    m <- .bin_mat()
    dm <- as.matrix(m)
    thr <- seq(2, 7, length.out = nrow(m))
    p <- thresholdParam("binarize", threshold = thr)
    expect_equal(as.matrix(processData(m, p)), (dm > thr) * 1, ignore_attr = TRUE)
    expect_equal(processData(dm, p), (dm > thr) * 1, ignore_attr = TRUE)
    expect_error(processData(m, thresholdParam("binarize", threshold = 1:2)),
                 "one per row")
    expect_error(thresholdParam("binarize", threshold = "median"))
})

test_that("a raster takes \"kmeans\" as its otsu default and refuses rank", {
    set.seed(2)
    r <- terra::rast(nrows = 30, ncols = 30,
                     vals = c(rnorm(450, 5), rnorm(450, 15)))
    expect_equal(terra::values(processData(r, thresholdParam("binarize"))),
                 terra::values(processData(r,
                     thresholdParam("binarize", threshold = "kmeans"))))
    expect_error(processData(r, thresholdParam("binarize", threshold = "rank")),
                 "raster")
})

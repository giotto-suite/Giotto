# .gini_rows() takes the gini coefficient of every row of a matrix at once,
# replacing a per-feature mygini_fun() call in the marker path.
#
# mygini_fun() is the reference. The requirement is bit-identity, not
# approximate agreement: the coefficients feed rank ties and strict `>`
# thresholds downstream, so a last-bit difference can change which features are
# returned. `expect_identical` throughout is deliberate.

.gr <- get(".gini_rows", envir = asNamespace("Giotto"))

# per-row reference, the implementation being replaced
.ref <- function(m, min_length = 0) {
    vapply(seq_len(nrow(m)),
        function(i) mygini_fun(m[i, ], min_length = min_length), numeric(1L))
}


test_that("bit-identical to mygini_fun across value shapes", {
    set.seed(1L)
    cases <- list(
        nonneg = matrix(runif(200 * 7, 0, 10), 200),
        wide = matrix(runif(200 * 20, 0, 10), 200),
        negatives = matrix(rnorm(200 * 7), 200),
        sparse_counts = matrix(rpois(200 * 7, 0.5), 200),
        constant = matrix(rep(3, 20 * 7), 20),
        heavy_ties = matrix(sample(c(0, 0, 0, 1, 2), 200 * 7, TRUE), 200),
        two_groups = matrix(runif(200 * 2, 0, 10), 200)
    )
    for (nm in names(cases)) {
        expect_identical(.gr(cases[[nm]]), .ref(cases[[nm]]), info = nm)
    }
})

test_that("bit-identical under min_length padding", {
    set.seed(2L)
    m <- matrix(runif(200 * 7, 0, 10), 200)
    neg <- matrix(rnorm(200 * 7), 200)
    ties <- matrix(sample(c(0, 1), 200 * 7, TRUE), 200)

    expect_identical(.gr(m, min_length = 20), .ref(m, min_length = 20))
    expect_identical(.gr(neg, min_length = 20), .ref(neg, min_length = 20))
    expect_identical(.gr(ties, min_length = 20), .ref(ties, min_length = 20))

    # min_length below the actual width never pads, so it is a no-op
    wide <- matrix(runif(50 * 20, 0, 10), 50)
    expect_identical(.gr(wide, min_length = 5), .gr(wide))
})

test_that("padding raises the ceiling, which is why it exists", {
    # A feature top-of-two cannot exceed (2-1)/2; padded to 20 the same values
    # are measured against (20-1)/20, so scores compare across group counts.
    # Here the padding is with the minimum, 0, so the padded vector is maximally
    # unequal and lands exactly on the higher ceiling.
    m <- matrix(c(10, 0), nrow = 1)
    expect_equal(.gr(m), (2 - 1) / 2)
    expect_equal(.gr(m, min_length = 20), (20 - 1) / 20)

    # Padding raises the CEILING, not the score, and is not monotone: a row
    # whose minimum is well above zero becomes MORE equal when padded with it,
    # so its coefficient falls. This is why `min_length` is documented as making
    # runs comparable rather than as a sensitivity knob, and why the default is
    # 0 -- it reorders features within a run.
    m2 <- matrix(c(10, 4), nrow = 1)
    expect_lt(.gr(m2, min_length = 20), .gr(m2))
})

test_that("an all-zero row is NaN, as mygini_fun returns", {
    m <- matrix(0, 3, 7)
    expect_identical(.gr(m), .ref(m))
    expect_true(all(is.nan(.gr(m))))
})

test_that("empty input returns an empty result rather than failing", {
    expect_identical(.gr(matrix(numeric(0), 0, 7)), numeric(0))
    expect_identical(.gr(matrix(numeric(0), 5, 0)), numeric(5))
})


# integration with the scoring step ####

test_that(".gini_score_dt agrees with the per-feature path it replaced", {
    gsd <- get(".gini_score_dt", envir = asNamespace("Giotto"))
    set.seed(3L)
    nf <- 60L; ng <- 5L
    dt <- data.table::data.table(
        feats = rep(paste0("g", seq_len(nf)), times = ng),
        cluster = rep(paste0("c", seq_len(ng)), each = nf),
        expression = runif(nf * ng, 0, 10),
        detection = runif(nf * ng)
    )
    got <- gsd(data.table::copy(dt), min_length = 0, min_expression = 0.2,
        min_detection = 0.2, min_expression_gini = -Inf,
        min_detection_gini = -Inf, rank_score = Inf, min_feats = 5)

    # reference coefficients, computed the old way
    ref_e <- dt[, .(g = mygini_fun(expression)), by = feats]
    ref_d <- dt[, .(g = mygini_fun(detection)), by = feats]
    expect_identical(
        got$expression_gini, ref_e$g[match(got$feats, ref_e$feats)]
    )
    expect_identical(
        got$detection_gini, ref_d$g[match(got$feats, ref_d$feats)]
    )
})

test_that("a non-rectangular table falls back and still scores", {
    # The grouped statistic emits the full cross product, so this should not
    # arise -- but the reshape is only valid when it holds, and the fallback is
    # what makes that assumption safe rather than load-bearing.
    gsd <- get(".gini_score_dt", envir = asNamespace("Giotto"))
    set.seed(4L)
    nf <- 30L; ng <- 4L
    dt <- data.table::data.table(
        feats = rep(paste0("g", seq_len(nf)), times = ng),
        cluster = rep(paste0("c", seq_len(ng)), each = nf),
        expression = runif(nf * ng, 0, 10),
        detection = runif(nf * ng)
    )
    holed <- dt[-c(3L, 40L, 77L)]   # a few (feature, cluster) cells missing
    got <- gsd(holed, min_length = 0, min_expression = 0, min_detection = 0,
        min_expression_gini = -Inf, min_detection_gini = -Inf,
        rank_score = Inf, min_feats = 5)
    expect_true(nrow(got) > 0L)
    expect_false(anyNA(got$expression_gini))
})

test_that("fewer than two groups is refused with a real message", {
    gsd <- get(".gini_score_dt", envir = asNamespace("Giotto"))
    one <- data.table::data.table(
        feats = paste0("g", 1:10), cluster = "c1",
        expression = runif(10), detection = runif(10)
    )
    expect_error(
        gsd(one, min_length = 0, min_expression = 0, min_detection = 0,
            min_expression_gini = -Inf, min_detection_gini = -Inf,
            rank_score = Inf, min_feats = 5),
        "at least 2"
    )
})

# ============================================================================
# Spatial autocorrelation statistics from the edge-statistics contract.
#
# Every function here takes the table `analyzeData(x, edgeStatsParam(),
# network =)` returns -- per-feature power sums, degree-weighted sums and the
# edge product sum, plus the "graph" attribute -- and nothing else. That table
# is the same whichever backend computed it (an in-memory matrix, a GiottoDisk
# store), so the statistics are written once.
#
# Significance is the randomisation null: the exact mean and variance of the
# statistic over all permutations of the values across cells (Cliff & Ord
# 1981), the moments shuffling would estimate, in closed form.
#
# Returned columns: feats, statistic, expected, variance, z, p.value. The
# alternative is positive autocorrelation (clustering).
# ============================================================================

# S0, S1, S2 of the symmetric weight matrix the network stands for, each edge
# entering it as w_ij = w_ji = w_e.
.autocor_graph_constants <- function(st) {
    g <- attr(st, "graph")
    if (is.null(g)) {
        stop("[spatial autocorrelation] the input lacks the \"graph\" ",
            "attribute of an edgeStatsParam result.", call. = FALSE)
    }
    list(
        N = g$n_cells,
        S0 = 2 * g$sum_w,
        S1 = 4 * g$sum_w2,
        S2 = 4 * g$sum_deg2,
        sum_w = g$sum_w
    )
}

# Central moments from power sums: sum of squared and fourth-power deviations,
# and the kurtosis b2 the randomisation variance needs.
.autocor_central <- function(st, N) {
    m <- st$sum / N
    z2 <- st$sumsq - N * m^2
    z4 <- st$sum4 - 4 * m * st$sum3 + 6 * m^2 * st$sumsq - 3 * N * m^4
    list(mean = m, z2 = z2, b2 = N * z4 / z2^2)
}

.autocor_result <- function(st, statistic, expected, variance, z) {
    data.table::data.table(
        feats = st$feats,
        statistic = statistic,
        expected = expected,
        variance = variance,
        z = z,
        p.value = stats::pnorm(z, lower.tail = FALSE)
    )
}

# Moran's I. The edge sum of centred values follows from the raw one:
# sum_e w (x_i - m)(x_j - m) = edge_prod - m * deg_sum + m^2 * sum_w.
.moran_test <- function(st) {
    g <- .autocor_graph_constants(st)
    N <- g$N; S0 <- g$S0; S1 <- g$S1; S2 <- g$S2
    cz <- .autocor_central(st, N)
    cross <- st$edge_prod - cz$mean * st$deg_sum + cz$mean^2 * g$sum_w
    I <- (N / S0) * 2 * cross / cz$z2
    EI <- -1 / (N - 1)
    b2 <- cz$b2
    EI2 <- (N * ((N^2 - 3 * N + 3) * S1 - N * S2 + 3 * S0^2) -
        b2 * ((N^2 - N) * S1 - 2 * N * S2 + 6 * S0^2)) /
        ((N - 1) * (N - 2) * (N - 3) * S0^2)
    V <- EI2 - EI^2
    .autocor_result(st, I, EI, V, (I - EI) / sqrt(V))
}

# Geary's C. sum_e w (x_i - x_j)^2 = deg_sumsq - 2 * edge_prod. Clustering
# lowers C, so z is taken as (1 - C) / sd to keep "greater" one-sided.
.geary_test <- function(st) {
    g <- .autocor_graph_constants(st)
    N <- g$N; S0 <- g$S0; S1 <- g$S1; S2 <- g$S2
    cz <- .autocor_central(st, N)
    C <- (N - 1) * 2 * (st$deg_sumsq - 2 * st$edge_prod) / (2 * S0 * cz$z2)
    b2 <- cz$b2
    V <- ((N - 1) * S1 * (N^2 - 3 * N + 3 - (N - 1) * b2) -
        0.25 * (N - 1) * S2 * (N^2 + 3 * N - 6 - (N^2 - N + 2) * b2) +
        S0^2 * (N^2 - 3 - (N - 1)^2 * b2)) /
        (N * (N - 2) * (N - 3) * S0^2)
    .autocor_result(st, C, 1, V, (1 - C) / sqrt(V))
}

# Black-black join count, on a 0/1 input (a binarized feature): the weighted
# number of edges with both ends on. Randomisation moments with the number of
# on cells fixed. Also reports the on count, the 2x2 edge table's odds ratio
# and, where the expected count is small, a Poisson tail instead of the
# normal one.
.joincount_test <- function(st, small = 10) {
    g <- .autocor_graph_constants(st)
    N <- g$N; S0 <- g$S0; S1 <- g$S1; S2 <- g$S2
    n1 <- st$sum
    ff <- function(k) {
        r <- rep(1, length(n1))
        for (j in 0:(k - 1)) r <- r * (n1 - j) / (N - j)
        r
    }
    p2 <- ff(2); p3 <- ff(3); p4 <- ff(4)
    BB <- st$edge_prod
    EB <- 0.5 * S0 * p2
    V <- 0.25 * (S1 * p2 + (S2 - 2 * S1) * p3 + (S0^2 + S1 - S2) * p4) - EB^2
    out <- .autocor_result(st, BB, EB, V, (BB - EB) / sqrt(V))
    use_pois <- EB < small & BB == round(BB)
    if (any(use_pois)) {
        out$p.value[use_pois] <- stats::ppois(BB[use_pois] - 1,
            EB[use_pois], lower.tail = FALSE)
    }
    disc <- st$deg_sum - 2 * BB          # n10 + n01
    n00 <- g$sum_w - BB - disc
    out[, "n_on" := n1]
    out[, "odds_ratio" := (BB * n00) / (disc / 2)^2]
    out[]
}

# Cauchy combination of p-values from dependent tests (ACAT): stays valid
# under arbitrary dependence, unlike Fisher's method. One p per row of `p`.
.acat <- function(p) {
    p <- as.matrix(p)
    p <- pmin(pmax(p, 1e-300), 1 - 1e-16)
    t <- rowMeans(tan((0.5 - p) * pi))
    0.5 - atan(t) / pi
}

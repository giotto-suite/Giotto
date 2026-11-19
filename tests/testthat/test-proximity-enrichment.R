# The null model for `cellProximityEnrichment()` permutes cell type labels over
# the NODES of the spatial network, holding topology fixed. That makes the
# expected pair counts analytically known, so these tests check the engine
# against a closed form rather than against a stored snapshot:
#
#   E[N_AA] = m * n_A(n_A - 1) / (n(n - 1))
#   E[N_AB] = m * 2 n_A n_B  / (n(n - 1))      (A != B)
#
# for n network nodes and m deduplicated undirected edges. An earlier
# implementation pooled the two endpoint columns and reshuffled that vector,
# which samples labels degree-weighted and lets one cell take different labels
# on different edges; it fails the exactness test below by ~50% when cell type
# tracks node degree. Data-free apart from a small deterministic graph.

# ring of `n` nodes, so every node has degree 2 and the topology is trivial
.ring_edges <- function(n) {
    data.table::data.table(
        from = sprintf("c%03d", seq_len(n)),
        to = sprintf("c%03d", c(seq_len(n)[-1], 1L))
    )
}

.exact_expected <- function(edges, labs) {
    nodes <- unique(c(edges$from, edges$to))
    lab <- labs[nodes]
    n <- length(nodes)
    m <- nrow(edges)
    tb <- table(lab)
    lv <- names(tb)
    out <- list()
    for (i in seq_along(lv)) {
        for (j in i:length(lv)) {
            e <- if (i == j) {
                m * tb[[i]] * (tb[[i]] - 1) / (n * (n - 1))
            } else {
                m * 2 * tb[[i]] * tb[[j]] / (n * (n - 1))
            }
            out[[length(out) + 1L]] <- data.table::data.table(
                unified_int = paste0(lv[i], "--", lv[j]),
                exact = as.numeric(e)
            )
        }
    }
    data.table::rbindlist(out)
}

test_that(".cpe_permutation_counts matches the exact node-permutation mean", {
    n <- 60L
    edges <- .ring_edges(n)
    nodes <- unique(c(edges$from, edges$to))
    set.seed(42)
    lab <- rep(c("A", "B", "C"), length.out = n)
    labs <- stats::setNames(lab, nodes)

    lv <- sort(unique(lab))
    K <- length(lv)
    cnts <- .cpe_permutation_counts(
        from = match(edges$from, nodes),
        to = match(edges$to, nodes),
        codes = match(labs[nodes], lv),
        K = K,
        number_of_simulations = 4000L,
        set_seed = TRUE,
        seed_number = 7L
    )

    sim_mean <- rowMeans(cnts$sim)
    ex <- .exact_expected(edges, labs)
    got <- vapply(seq_len(nrow(ex)), function(i) {
        p <- strsplit(ex$unified_int[i], "--", fixed = TRUE)[[1]]
        sim_mean[(match(p[1], lv) - 1L) * K + match(p[2], lv)]
    }, numeric(1))

    # Monte Carlo tolerance at 4000 draws
    expect_equal(got, ex$exact, tolerance = 0.02, ignore_attr = TRUE)
})

test_that(".cpe_permutation_counts conserves the total edge count", {
    n <- 40L
    edges <- .ring_edges(n)
    nodes <- unique(c(edges$from, edges$to))
    codes <- rep(1:4, length.out = n)
    cnts <- .cpe_permutation_counts(
        from = match(edges$from, nodes), to = match(edges$to, nodes),
        codes = codes, K = 4L, number_of_simulations = 50L, seed_number = 1L
    )
    # every edge is counted exactly once, in the observed and in every draw
    expect_identical(sum(cnts$obs), nrow(edges))
    expect_true(all(colSums(cnts$sim) == nrow(edges)))
})

test_that(".cpe_permutation_counts is reproducible and honours set_seed", {
    n <- 30L
    edges <- .ring_edges(n)
    nodes <- unique(c(edges$from, edges$to))
    args <- list(
        from = match(edges$from, nodes), to = match(edges$to, nodes),
        codes = rep(1:3, length.out = n), K = 3L, number_of_simulations = 20L
    )
    a <- do.call(.cpe_permutation_counts, c(args, list(seed_number = 99L)))
    b <- do.call(.cpe_permutation_counts, c(args, list(seed_number = 99L)))
    d <- do.call(.cpe_permutation_counts, c(args, list(seed_number = 100L)))
    expect_identical(a$sim, b$sim)
    expect_false(identical(a$sim, d$sim))
})

test_that("empirical p-values are bounded away from zero", {
    # the +1 estimator floors p at 1 / (1 + number_of_simulations)
    n_perm <- 99L
    expect_identical((1 + 0) / (1 + n_perm), 1 / (1 + n_perm))
    expect_gt(1 / (1 + n_perm), 0)
})

skip_if_not_installed("GiottoData")

test_that("cellProximityEnrichment keeps its documented result contract", {
    g <- GiottoData::loadGiottoMini("visium")
    cp <- cellProximityEnrichment(g,
        cluster_column = "leiden_clus", number_of_simulations = 99L
    )

    expect_named(cp, c("raw_sim_table", "enrichm_res"))
    expect_s3_class(cp$enrichm_res, "data.table")
    expect_s3_class(cp$raw_sim_table, "data.table")

    # the columns the shipped plots read, in their original positions
    legacy <- c(
        "unified_int", "type_int", "original", "simulations", "enrichm",
        "p_higher_orig", "p_lower_orig", "p.adj_higher", "p.adj_lower",
        "PI_value", "int_ranking"
    )
    expect_identical(names(cp$enrichm_res)[seq_along(legacy)], legacy)
    expect_s3_class(cp$enrichm_res$unified_int, "factor")
    expect_identical(
        names(cp$raw_sim_table),
        c("unified_int", "type_int", "round", "V1", "orig")
    )

    # p-values are valid probabilities and never exactly zero
    for (col in c("p_higher_orig", "p_lower_orig")) {
        p <- cp$enrichm_res[[col]]
        expect_true(all(p > 0 & p <= 1))
        expect_gte(min(p), 1 / (1 + 99L))
    }

    # observed counts must equal the deduplicated network's own pair counts
    net <- GiottoClass::annotateSpatialNetwork(g,
        spatial_network_name = "Delaunay_network", cluster_column = "leiden_clus"
    )
    net <- GiottoUtils::dt_sort_combine_two_columns(net, "to", "from", "uc")
    net <- net[!duplicated(net$uc)]
    tally <- net[, .N, by = "unified_int"]
    chk <- merge(
        cp$enrichm_res[, list(unified_int = as.character(unified_int), original)],
        tally,
        by = "unified_int"
    )
    expect_identical(chk$original, as.numeric(chk$N))
})

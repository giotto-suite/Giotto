# Tests for the bag model behind runNucleusDeconv(): a cell-level classifier
# trained from group-level targets (R/bag-model.R). The DWLS expression loss
# itself is tested in spatialDWLS; here: the gradient of the assembled loss,
# reproducibility, that training recovers planted cell labels, and the edge
# cases of the groups.

# Cells of three types with type-specific features, pooled into groups of
# 3-10 cells; the targets are the groups' true proportions and expression.
.bag_fixture <- function(n_groups = 120L, d = 20L, seed = 1L) {
    set.seed(seed)
    types <- c("A", "B", "C")
    n <- sample(3:10, n_groups, replace = TRUE)
    lab <- sample(types, sum(n), replace = TRUE, prob = c(0.5, 0.3, 0.2))
    ids <- paste0("c", seq_along(lab))
    centre <- matrix(stats::rnorm(3 * d, sd = 1.5), 3, d,
        dimnames = list(types, NULL))
    X <- centre[lab, ] + matrix(stats::rnorm(length(lab) * d), length(lab), d)
    rownames(X) <- ids
    gid <- paste0("g", rep(seq_len(n_groups), n))
    groups <- data.table::data.table(group = gid, cell_ID = ids)
    prop <- as.matrix(table(factor(gid, unique(gid)), factor(lab, types)))
    prop <- prop / rowSums(prop)
    prop <- matrix(prop, nrow(prop), dimnames = list(unique(gid), types))
    genes <- paste0("f", 1:30)
    sig <- matrix(1, 30, 3, dimnames = list(genes, types))
    sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
    counts <- sig %*% t(prop * n)
    counts[] <- stats::rpois(length(counts), counts)
    list(X = X, lab = stats::setNames(lab, ids), groups = groups,
        prop = prop, types = types,
        expr = log2(counts + 1),
        clusters = stats::setNames(rep("all", n_groups), unique(gid)))
}

sig_of <- function(f) {
    genes <- rownames(f$expr)
    sig <- matrix(1, length(genes), 3, dimnames = list(genes, f$types))
    sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
    sig
}

test_that("the assembled loss has the gradient of its finite differences", {
    f <- .bag_fixture(n_groups = 12L, d = 5L)
    learner <- .bag_learner("linear")
    set.seed(2)
    theta <- learner$init(ncol(f$X), f$types)
    cells <- .bag_cells(f$groups, rownames(f$X))
    ids <- names(cells)[1:8]
    losses <- list(
        .bag_loss("proportion", target = f$prop, weight = 0.99),
        .bag_loss("confidence", weight = 0.01)
    )
    if (requireNamespace("spatialDWLS", quietly = TRUE) &&
            exists("dwls_design", asNamespace("spatialDWLS"))) {
        design <- spatialDWLS::dwls_design(f$expr, sig_of(f), f$clusters)
        # J = 1 caps every weight at 1: the weights are then fixed and the
        # gradient is exact, not only with the weights held
        design$J[] <- 1
        losses <- c(losses, list(.bag_loss("dwls", design = design)))
    }
    g <- .bag_batch(theta, learner, f$X, cells, ids, losses)$grad
    h <- 1e-6
    for (k in seq_len(6)) {
        i <- sample(nrow(theta$W), 1)
        j <- sample(ncol(theta$W), 1)
        up <- dn <- theta
        up$W[i, j] <- up$W[i, j] + h
        dn$W[i, j] <- dn$W[i, j] - h
        num <- (.bag_batch(up, learner, f$X, cells, ids, losses, FALSE) -
            .bag_batch(dn, learner, f$X, cells, ids, losses, FALSE)) / (2 * h)
        expect_equal(unname(g$W[i, j]), num, tolerance = 1e-5)
    }
    j <- 2
    up <- dn <- theta
    up$b[j] <- up$b[j] + h
    dn$b[j] <- dn$b[j] - h
    expect_equal(unname(g$b[j]),
        (.bag_batch(up, learner, f$X, cells, ids, losses, FALSE) -
            .bag_batch(dn, learner, f$X, cells, ids, losses, FALSE)) / (2 * h),
        tolerance = 1e-5)
})

test_that("training from group proportions recovers the cell labels", {
    f <- .bag_fixture()
    m <- .bag_train(f$X, f$groups,
        losses = list(.bag_loss("proportion", target = f$prop, weight = 0.99),
            .bag_loss("confidence", weight = 0.01)),
        types = f$types, epochs = 60L, lr = 1e-2, batch_size = 16L, seed = 1L)
    P <- .bag_predict(m, f$X)
    expect_identical(dimnames(P), list(rownames(f$X), f$types))
    expect_equal(unname(rowSums(P)), rep(1, nrow(P)))
    expect_gt(mean(colnames(P)[max.col(P)] == f$lab), 0.9)
    expect_true(is.finite(m$best_epoch))
    expect_length(m$groups$train, 96L)
})

test_that("training from group expression alone recovers the cell labels", {
    skip_if_not_installed("spatialDWLS")
    skip_if_not(exists("dwls_design", asNamespace("spatialDWLS")))
    f <- .bag_fixture()
    design <- spatialDWLS::dwls_design(f$expr, sig_of(f), f$clusters)
    m <- .bag_train(f$X, f$groups,
        losses = list(.bag_loss("dwls", design = design),
            .bag_loss("confidence", weight = 0.01)),
        types = f$types, epochs = 60L, lr = 1e-2, batch_size = 16L, seed = 1L)
    P <- .bag_predict(m, f$X)
    expect_gt(mean(colnames(P)[max.col(P)] == f$lab), 0.9)
})

test_that("the same seed gives the same model and leaves the session RNG alone", {
    f <- .bag_fixture(n_groups = 40L)
    L <- list(.bag_loss("proportion", target = f$prop))
    set.seed(10)
    before <- .Random.seed
    m1 <- .bag_train(f$X, f$groups, L, types = f$types, epochs = 5L, seed = 3L)
    expect_identical(.Random.seed, before)
    m2 <- .bag_train(f$X, f$groups, L, types = f$types, epochs = 5L, seed = 3L)
    expect_identical(m1$theta, m2$theta)
    m3 <- .bag_train(f$X, f$groups, L, types = f$types, epochs = 5L, seed = 4L)
    expect_false(identical(m1$theta, m3$theta))
})

test_that("a warm start switches losses and picks the epoch in the main phase", {
    f <- .bag_fixture(n_groups = 40L)
    m <- .bag_train(f$X, f$groups,
        losses = list(.bag_loss("proportion", target = f$prop, phase = "warm"),
            .bag_loss("confidence", weight = 0.01)),
        types = f$types, epochs = 8L, warm = 5L, seed = 1L)
    expect_gt(m$best_epoch, 5L)
})

test_that("groups without cells or targets are dropped; too few is an error", {
    f <- .bag_fixture(n_groups = 40L)
    extra <- data.table::data.table(
        group = c("empty", "no_target"), cell_ID = c("not_a_cell", "c1"))
    cells <- .bag_cells(rbind(f$groups, extra), rownames(f$X))
    expect_false("empty" %in% names(cells))
    m <- .bag_train(f$X, rbind(f$groups, extra),
        list(.bag_loss("proportion", target = f$prop)),
        types = f$types, epochs = 2L, seed = 1L)
    expect_false("no_target" %in% unlist(m$groups))
    expect_error(
        .bag_train(f$X, f$groups[group %in% c("g1", "g2")],
            list(.bag_loss("proportion", target = f$prop)),
            types = f$types, epochs = 2L),
        "fewer than 10 groups")
    expect_error(.bag_loss("proportion", target = unname(f$prop)), "row names")
    expect_error(.bag_train(unname(f$X), f$groups, list(), types = f$types),
        "row names")
})

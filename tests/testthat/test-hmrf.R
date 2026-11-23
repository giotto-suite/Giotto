options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
    giotto.no_python_warn = TRUE)

sim_tissue <- function(n = 300, seed = 1) {
    set.seed(seed)
    xy <- matrix(runif(2 * n), ncol = 2) * 10
    dom <- 1L + (xy[, 1] > 5)
    y <- matrix(rnorm(n * 6), n) + c(-1.5, 1.5)[dom]
    ids <- paste0("c", seq_len(n))
    rownames(y) <- ids
    knn <- dbscan::kNN(xy, k = 5)$id
    e <- unique(t(apply(cbind(rep(seq_len(n), 5), as.vector(knn)), 1, sort)))
    list(y = y, dom = dom,
        net = data.frame(from = ids[e[, 1]], to = ids[e[, 2]]))
}

test_that("clusterData(hmrfParam) returns one domain per cell", {
    s <- sim_tissue()
    res <- clusterData(s$y, hmrfParam(k = 2, beta = 2, nstart = 10,
        zscore = "none"), network = s$net)
    expect_named(res, c("cell_ID", "cluster"))
    expect_equal(res$cell_ID, rownames(s$y))
    expect_equal(dim(attr(res, "prob")), c(nrow(s$y), 2L))
    expect_equal(unname(rowSums(attr(res, "prob"))), rep(1, nrow(s$y)))
    expect_gt(aricode::ARI(res$cluster, s$dom), 0.9)
})

test_that("blocks are independent sets covering every cell once", {
    s <- sim_tissue()
    adj <- Giotto:::.hmrf_adjacency(s$net, rownames(s$y))
    b <- Giotto:::.hmrf_blocks(adj)
    expect_setequal(unlist(b), seq_len(nrow(s$y)))
    for (bl in b) expect_equal(sum(adj[bl, bl]), 0)
})

test_that("cells without neighbours are handled", {
    s <- sim_tissue()
    net <- s$net[s$net$from != "c1" & s$net$to != "c1", ]
    res <- clusterData(s$y, hmrfParam(k = 2, beta = 2, nstart = 5,
        zscore = "none"), network = net)
    expect_false(anyNA(res$cluster))
})

test_that("matches smfishHmrf on identical inputs", {
    skip_if_not_installed("smfishHmrf")
    s <- sim_tissue()
    ids <- rownames(s$y)
    adj <- Giotto:::.hmrf_adjacency(s$net, ids)
    blocks <- Giotto:::.hmrf_blocks(adj)
    p <- hmrfParam(k = 2, beta = 2, nstart = 10, zscore = "none")@param
    init <- Giotto:::.hmrf_init(s$y, p)
    nb <- split(c(s$net$to, s$net$from), c(s$net$from, s$net$to))[ids]
    nl <- lapply(nb, function(v) sort(match(v, ids)))
    numnei <- lengths(nl)
    nei <- matrix(-1L, length(ids), max(numnei))
    for (i in seq_along(ids)) nei[i, seq_len(numnei[i])] <- nl[[i]]
    ref <- smfishHmrf::smfishHmrf.hmrfem.multi(y = s$y, neighbors = nei,
        beta = 2, numnei = as.integer(numnei), blocks = blocks,
        mu = init$mu, sigma = init$sigma, verbose = FALSE, err = 1e-7,
        maxit = 50, dampFactor = init$damp)
    nat <- Giotto:::.hmrf_fit(s$y, adj, blocks, init, p)
    expect_identical(as.integer(nat$class), as.integer(ref$class))
    expect_equal(nat$prob, ref$prob, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(nat$likelihood, ref$likelihood, tolerance = 1e-10)
})

test_that("runHMRF stores domains, probabilities and its param", {
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    g <- runHMRF(g, feats_to_use = "hvf", k = 4, betas = c(0, 2),
        nstart = 5, verbose = FALSE)
    md <- pDataDT(g)
    expect_true(all(c("hmrf_k4_b0", "hmrf_k4_b2") %in% names(md)))
    enr <- getSpatialEnrichment(g, name = "hmrf_k4_b2", output = "data.table")
    expect_equal(nrow(enr), nrow(md))
    h <- g@parameters[[length(g@parameters)]]
    expect_s4_class(attr(h, "attachments")$param, "hmrfParam")
})

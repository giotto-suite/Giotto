.wnn_fixture <- function(n = 450, protein = c("structured", "noise"),
    dup = 0L, dup_far = FALSE, seed = 7, knn_dims = 1:8, knn_k = 20,
    backend = NULL) {
    protein <- match.arg(protein)
    withr::local_seed(seed)
    cells <- sprintf("c_%04d", seq_len(n))
    cl <- rep_len(1:3, n)
    xr <- matrix(stats::rnorm(n * 20, sd = 0.5), n, 20)
    xr[, 1:3] <- xr[, 1:3] + diag(6, 3)[cl, ]
    xp <- matrix(stats::rnorm(n * 10, sd = 0.5), n, 10)
    if (protein == "structured") {
        xp[, 1:2] <- xp[, 1:2] + diag(6, 2)[pmin(cl, 2L), ]
    }
    if (dup > 0L) {
        xr[seq_len(dup), ] <- if (dup_far) 100 else
            matrix(xr[1, ], dup, 20, byrow = TRUE)
        xp[seq_len(dup), ] <- if (dup_far) 100 else
            matrix(xp[1, ], dup, 10, byrow = TRUE)
    }
    dimnames(xr) <- list(cells, paste0("Dim.", 1:20))
    dimnames(xp) <- list(cells, paste0("Dim.", 1:10))
    mr <- Matrix::Matrix(matrix(stats::rpois(10 * n, 5), 10, n,
        dimnames = list(paste0("g", 1:10), cells)), sparse = TRUE)
    mp <- Matrix::Matrix(matrix(stats::rpois(5 * n, 5), 5, n,
        dimnames = list(paste0("p", 1:5), cells)), sparse = TRUE)
    g <- createGiottoObject(
        expression = list(cell = list(rna = list(raw = mr),
            protein = list(raw = mp))),
        backend = backend, verbose = FALSE
    )
    g <- setGiotto(g, createDimObj(coordinates = xr, name = "pca",
        spat_unit = "cell", feat_type = "rna", method = "pca"),
        verbose = FALSE)
    g <- setGiotto(g, createDimObj(coordinates = xp, name = "protein.pca",
        spat_unit = "cell", feat_type = "protein", method = "pca"),
        verbose = FALSE)
    for (ft in c("rna", "protein")) {
        g <- createNearestNetwork(g, feat_type = ft, type = "kNN", k = knn_k,
            dimensions_to_use = knn_dims, verbose = FALSE)
    }
    list(g = g, cluster = cl)
}

.wnn_paper <- function(g, integrated_snn = "giotto", add_to_metadata = TRUE,
    ...) {
    runWNN(g, graph = "knn", integrated_snn = integrated_snn,
        add_to_metadata = add_to_metadata, ...)
}

.wnn_knn <- function(g, name = "integrated_kNN") {
    .wnn_stored_integrated_knn(g, "cell", "rna", network_name = name)
}

.wnn_weights <- function(g, w_names = c("w_rna", "w_protein")) {
    out <- lapply(w_names, function(n) {
        x <- getMultiomics(g, feat_type = "rna_protein", result_name = n)
        x[[1]]
    })
    names(out) <- w_names
    out
}

.snn_sigma_literal <- function(X, nn_self, nearest_dist, n = 20L) {
    N <- nrow(X)
    k <- ncol(nn_self)
    A <- Matrix::sparseMatrix(i = rep(seq_len(N), k), j = as.vector(nn_self),
        x = 1, dims = c(N, N))
    S <- methods::as(A %*% Matrix::t(A), "CsparseMatrix")
    S@x <- S@x / (2 * k - S@x)
    out <- numeric(N)
    for (i in seq_len(N)) {
        idx <- S@i[(S@p[i] + 1):S@p[i + 1]] + 1L
        val <- S@x[(S@p[i] + 1):S@p[i + 1]]
        o <- order(val, method = "radix")
        val <- val[o]
        idx <- idx[o]
        n_i <- min(n, length(val))
        dists <- numeric(0)
        for (j in seq_along(val)) {
            if (length(dists) < n_i || val[j] == val[n_i]) {
                res <- sqrt(sum((X[idx[j], ] - X[i, ])^2))
                if (nearest_dist[i] > 0) res <- max(res - nearest_dist[i], 0)
                dists <- c(dists, res)
            } else {
                break
            }
        }
        out[i] <- if (length(dists) > n_i) {
            mean(sort(dists, decreasing = TRUE)[seq_len(n_i)])
        } else {
            mean(dists)
        }
    }
    out
}

.snn_weights_literal <- function(emb, k_nn = 20L) {
    l2 <- lapply(emb, function(x) {
        y <- x / sqrt(rowSums(x^2))
        y[!is.finite(y)] <- 0
        y
    })
    nn <- lapply(l2, function(x) {
        cbind(seq_len(nrow(x)), dbscan::kNN(x, k = k_nn - 1L, sort = TRUE)$id)
    })
    near <- lapply(names(l2), function(r) {
        sqrt(rowSums((l2[[r]] - l2[[r]][nn[[r]][, 2], ])^2))
    })
    names(near) <- names(l2)
    pred <- function(x, idx) {
        acc <- 0
        for (j in 2:ncol(idx)) acc <- acc + x[idx[, j], ]
        acc / (ncol(idx) - 1)
    }
    score <- list()
    for (r in names(l2)) {
        r2 <- setdiff(names(l2), r)
        d_w <- pmax(sqrt(rowSums((l2[[r]] - pred(l2[[r]], nn[[r]]))^2)) -
            near[[r]], 0)
        d_c <- pmax(sqrt(rowSums((l2[[r]] - pred(l2[[r]], nn[[r2]]))^2)) -
            near[[r]], 0)
        sig <- .snn_sigma_literal(l2[[r]], nn[[r]], near[[r]], n = k_nn)
        s <- exp(-d_w / sig) / (exp(-d_c / sig) + 1e-4)
        score[[r]] <- pmin(pmax(s, 0), 200)
    }
    tot <- Reduce(`+`, lapply(score, exp))
    lapply(score, function(s) exp(s) / tot)
}

test_that("graph = \"dense\" (the default) is Giotto 4.3.6", {
    for (protein in c("structured", "noise")) {
        fx <- .wnn_fixture(protein = protein)
        ref <- suppressMessages(.wnn_dense_reference(fx$g,
            reduction_names = c("pca", "protein.pca"),
            w_names = c("w_rna", "w_protein")))
        g <- suppressMessages(runWNN(fx$g))
        expect_identical(g@multiomics, ref@multiomics)
        if (protein == "noise") {
            expect_gt(sum(is.nan(.wnn_weights(g)$w_rna)), 0)
        }
        ref <- suppressMessages(.wnn_umap_dense_reference(ref))
        g <- suppressMessages(runIntegratedUMAP(g))
        expect_identical(g@multiomics, ref@multiomics)
        for (ft in c("rna", "protein")) {
            expect_identical(
                getDimReduction(g, feat_type = ft, reduction_method = "umap",
                    name = "integrated.umap", output = "matrix"),
                getDimReduction(ref, feat_type = ft, reduction_method = "umap",
                    name = "integrated.umap", output = "matrix"))
        }
        expect_identical(
            getNearestNetwork(g, feat_type = "rna", nn_type = "kNN",
                name = "integrated_kNN", output = "data.table"),
            getNearestNetwork(ref, feat_type = "rna", nn_type = "kNN",
                name = "integrated_kNN", output = "data.table"))
    }
})

test_that("graph = \"dense\" warns about the other graph's arguments", {
    fx <- .wnn_fixture()
    expect_warning(suppressMessages(runWNN(fx$g, engine = "hnsw",
        knn_range = 50)), "ignored with graph = \"dense\": knn_range, engine")
    g <- suppressMessages(runWNN(fx$g))
    expect_warning(suppressMessages(runIntegratedUMAP(g, name = "u")),
        "ignored with graph = \"dense\": name")
    gk <- .wnn_paper(fx$g)
    expect_warning(runIntegratedUMAP(gk, graph = "knn", force = TRUE,
        matrix_result_name = "m"),
        "ignored with graph = \"knn\": matrix_result_name, force")
    expect_no_warning(runIntegratedUMAP(gk, graph = "knn"))
})

test_that("graph = \"dense\" needs two feature types and notes its size", {
    withr::local_options(giotto.warn_wnn_dense = TRUE)
    expect_no_message(.wnn_dense_notice(1000L))
    expect_message(.wnn_dense_notice(60000L), "26.8 GB each")
    expect_message(.wnn_dense_notice(60000L), "graph = \"knn\"")
    withr::local_options(giotto.warn_wnn_dense = FALSE)
    expect_no_message(.wnn_dense_notice(60000L))

    fx <- .wnn_fixture()
    g <- fx$g
    cells <- spatIDs(g)
    xa <- matrix(stats::rnorm(length(cells) * 5), ncol = 5,
        dimnames = list(cells, paste0("Dim.", 1:5)))
    g <- setGiotto(g, createExprObj(Matrix::Matrix(matrix(1, 3,
        length(cells), dimnames = list(paste0("a", 1:3), cells)),
        sparse = TRUE), name = "raw", feat_type = "atac"), verbose = FALSE)
    g <- setGiotto(g, createDimObj(coordinates = xa, name = "atac.pca",
        spat_unit = "cell", feat_type = "atac", method = "pca"),
        verbose = FALSE)
    expect_error(runWNN(g, feat_types = c("rna", "protein", "atac")),
        "integrates two feature types")
})

test_that("reduction_names defaults to the stored names", {
    fx <- .wnn_fixture()
    ft <- c("rna", "protein")
    expect_identical(.wnn_reduction_names(fx$g, "cell", ft, c("pca", "pca")),
        c("pca", "protein.pca"))
    xr <- getDimReduction(fx$g, feat_type = "rna", output = "matrix")
    g <- setGiotto(fx$g, createDimObj(coordinates = xr, name = "rna.pca",
        spat_unit = "cell", feat_type = "rna", method = "pca"),
        verbose = FALSE)
    expect_identical(.wnn_reduction_names(g, "cell", ft, c("pca", "pca")),
        c("rna.pca", "protein.pca"))
    expect_identical(.wnn_reduction_names(g, "cell", ft, c("lsi", "lsi")),
        c("lsi", "protein.lsi"))
})

test_that("weights are never NaN", {
    cases <- list(
        list(protein = "noise", dup = 0L, far = FALSE),
        list(protein = "structured", dup = 40L, far = FALSE),
        list(protein = "structured", dup = 40L, far = TRUE)
    )
    for (cs in cases) {
        fx <- .wnn_fixture(protein = cs$protein, dup = cs$dup,
            dup_far = cs$far)
        for (l2 in c(TRUE, FALSE)) {
            g <- .wnn_paper(fx$g, l2_norm = l2, knn_range = 60)
            w <- .wnn_weights(g)
            expect_true(all(is.finite(w$w_rna)))
            expect_true(all(w$w_rna >= 0 & w$w_rna <= 1))
            expect_equal(unname(w$w_rna + w$w_protein),
                rep(1, length(w$w_rna)), tolerance = 1e-12)
            expect_true(all(is.finite(.wnn_knn(g)$dist)))
        }
    }
})

test_that("a zero bandwidth gives finite kernels", {
    expect_identical(.wnn_kernel(c(0, 1, 0, 2), c(0, 0, 1, 1)),
        c(1, 0, 1, exp(-2)))
    expect_identical(.wnn_kernel(matrix(c(0, 1, 0, 3), 2), c(0, 1)),
        matrix(c(1, exp(-1), 1, exp(-3)), 2))
    expect_identical(.wnn_l2(matrix(0, 2, 3)), matrix(0, 2, 3))

    fx <- .wnn_fixture(dup = 40L, dup_far = TRUE)
    x <- getDimReduction(fx$g, feat_type = "protein", name = "protein.pca",
        output = "matrix")
    nb <- dbscan::kNN(x, k = 20)$id
    rho <- .wnn_row_min(.wnn_nb_dist(x, nb))
    s <- .wnn_sigma_snn(x, nb, rho)
    expect_true(all(is.finite(s)))
    expect_true(all(s[1:40] == 0))
    expect_true(all(s[-(1:40)] > 0))
    w <- .wnn_weights(.wnn_paper(fx$g, dimensions_to_use = list(1:20, 1:10)))
    expect_equal(unname(w$w_rna[1:40]), rep(0.5, 40))
})

test_that("the snn bandwidth matches a literal port of the WNN reference code", {
    withr::local_seed(3)
    n <- 600
    x <- matrix(stats::rnorm(n * 10), n, 10) +
        matrix(stats::rnorm(80, sd = 3), 8, 10)[sample.int(8, n, TRUE), ]
    x[1:30, ] <- matrix(x[1, ], 30, 10, byrow = TRUE)
    kn <- dbscan::kNN(x, k = 19, sort = TRUE)
    rho <- kn$dist[, 1]
    expect_equal(.wnn_sigma_snn(x, kn$id, rho, chunk = 137L),
        .snn_sigma_literal(x, cbind(seq_len(n), kn$id), rho),
        tolerance = 1e-12)
})

test_that("snn weights equal the WNN reference weights at k = 19", {
    fx <- .wnn_fixture()
    g <- .wnn_paper(fx$g, k = 19, dimensions_to_use = list(1:20, 1:10))
    w <- .wnn_weights(g)
    emb <- list(
        rna = getDimReduction(fx$g, feat_type = "rna", output = "matrix"),
        protein = getDimReduction(fx$g, feat_type = "protein",
            name = "protein.pca", output = "matrix")
    )
    ref <- .snn_weights_literal(emb, k_nn = 20L)
    expect_equal(unname(w$w_rna), unname(ref$rna), tolerance = 1e-10)
    expect_equal(unname(w$w_protein), unname(ref$protein), tolerance = 1e-10)
})

test_that("the sparse integrated kNN matches the dense first-order top k", {
    fx <- .wnn_fixture(n = 600)
    g <- .wnn_paper(fx$g, knn_range = 200)
    inp <- .wnn_inputs(fx$g, "cell", c("rna", "protein"),
        c("pca", "pca"), c("pca", "protein.pca"), NULL, NULL,
        l2_norm = TRUE)
    nb <- lapply(inp$X, function(x) dbscan::kNN(x, k = 20)$id)
    mw <- .wnn_modality_weights(inp$X, nb)
    w_stored <- .wnn_weights(g)
    expect_equal(unname(mw$w$rna), unname(w_stored$w_rna),
        tolerance = 1e-14)
    n <- length(inp$cells)
    tt <- 0
    for (m in names(inp$X)) {
        dd <- as.matrix(stats::dist(inp$X[[m]]))
        kern <- .wnn_kernel(pmax(dd - mw$rho[[m]], 0), mw$bw[[m]])
        tt <- tt + mw$w[[m]] * kern
    }
    diag(tt) <- -Inf
    dense <- t(apply(tt, 1, function(r) order(-r)[1:20]))
    id <- .wnn_knn(g)$id
    overlap <- vapply(seq_len(n), function(i) {
        length(intersect(id[i, ], dense[i, ])) / 20
    }, numeric(1))
    expect_gte(mean(overlap), 0.99)
    net <- getNearestNetwork(g, feat_type = "rna", nn_type = "kNN",
        name = "integrated_kNN", output = "data.table")
    ref <- tt[cbind(match(net$from, inp$cells), match(net$to, inp$cells))]
    expect_equal(net$affinity, ref, tolerance = 1e-12)
})

test_that("results are stored where downstream functions find them", {
    fx <- .wnn_fixture()
    g <- .wnn_paper(fx$g, matrix_result_name = "theta_weighted_matrix")
    n <- length(spatIDs(g))
    mo <- g@multiomics$cell$rna_protein$WNN
    expect_setequal(names(mo), c("w_rna", "w_protein",
        "theta_weighted_matrix"))
    expect_identical(names(getMultiomics(g, feat_type = "rna_protein",
        result_name = "w_protein")), "protein")
    nn <- .wnn_knn(g)
    id <- nn$id
    d <- nn$dist
    expect_identical(dim(id), c(n, 20L))
    expect_identical(rownames(id), spatIDs(g))
    expect_false(any(id == seq_len(n)))
    expect_true(all(d[, -1] >= d[, -20]))
    tw <- mo$theta_weighted_matrix
    expect_s4_class(tw, "dgCMatrix")
    expect_identical(Matrix::nnzero(tw), n * 20L)
    expect_identical(rownames(tw), spatIDs(g))

    net <- getNearestNetwork(g, feat_type = "rna", nn_type = "kNN",
        name = "integrated_kNN")
    expect_s4_class(net, "nnNetObj")
    expect_identical(net@misc$nn_input$method, "WNN")
    edt <- getNearestNetwork(g, feat_type = "rna", nn_type = "kNN",
        name = "integrated_kNN", output = "data.table")
    expect_true(all(c("from", "to", "weight", "distance", "affinity") %in%
        names(edt)))
    expect_equal(edt$weight, 1 / (1 + edt$distance))
    cells <- spatIDs(g)
    expect_equal(tw[cbind(match(edt$from, cells), match(edt$to, cells))],
        edt$affinity)
    expect_null(g@nn_network$cell$rna$sNN$integrated_kNN)
    snet <- getNearestNetwork(g, feat_type = "rna", nn_type = "sNN",
        name = "integrated_sNN")
    expect_identical(snet@misc$nn_input$type, "sNN")
    expect_identical(snet@misc$nn_input$minimum_shared, 5L)
    sdt <- getNearestNetwork(g, feat_type = "rna", nn_type = "sNN",
        name = "integrated_sNN", output = "data.table")
    pair <- function(a, b) paste(pmin(a, b), pmax(a, b))
    expect_true(all(pair(sdt$from, sdt$to) %in% pair(edt$from, edt$to)))
    expect_identical(anyDuplicated(pair(sdt$from, sdt$to)), 0L)
    expect_true(all(sdt$shared >= 1))
    expect_equal(sdt$weight, 1 / (1 + sdt$distance))
    knn <- structure(list(dist = unname(d), id = unname(id), k = 20L,
        sort = TRUE), class = c("kNN", "NN"))
    ref <- GiottoClass::createNetwork(matrix(0, n, 1),
        GiottoClass::sNNNetworkParam(k = 20L, output = "data.table"),
        node_ids = cells, nn_network = knn, verbose = FALSE)
    expect_setequal(pair(sdt$from, sdt$to), pair(ref$from, ref$to))

    for (ft in c("rna", "protein")) {
        meta <- pDataDT(g, feat_type = ft)
        expect_equal(meta$w_rna + meta$w_protein, rep(1, n))
    }

    g <- runIntegratedUMAP(g, graph = "knn")
    for (ft in c("rna", "protein")) {
        u <- getDimReduction(g, feat_type = ft, reduction_method = "umap",
            name = "integrated.umap", output = "matrix")
        expect_identical(dim(u), c(n, 2L))
    }
    g <- doLeidenCluster(g, feat_type = "rna", nn_network_to_use = "kNN",
        network_name = "integrated_kNN", name = "integrated_leiden",
        resolution = 0.5)
    cl <- pDataDT(g, feat_type = "rna")$integrated_leiden
    expect_identical(length(cl), n)
    g <- doLeidenCluster(g, feat_type = "rna", nn_network_to_use = "sNN",
        network_name = "integrated_sNN", name = "integrated_leiden_snn",
        resolution = 0.5)
    expect_true(length(unique(pDataDT(g)$integrated_leiden_snn)) > 1L)
    p <- plotUMAP(g, dim_reduction_name = "integrated.umap",
        cell_color = "integrated_leiden", show_plot = FALSE,
        return_plot = TRUE, save_plot = FALSE)
    expect_s3_class(p, "ggplot")
    expect_false(any(rapply(g@multiomics, function(x) {
        length(dim(x)) == 2L && nrow(x) == n && ncol(x) == n && is.matrix(x)
    }, how = "unlist")))
})

test_that("custom names, and runIntegratedUMAP embeds the stored graph", {
    fx <- .wnn_fixture()
    g <- .wnn_paper(fx$g, w_names = c("wr", "wp"), integrated_k = 15,
        integrated_network_name = "wnn15", integrated_feat_type = "mo",
        add_to_metadata = FALSE, integrated_snn = "none")
    expect_null(g@nn_network$cell$rna$sNN)
    expect_setequal(names(g@multiomics$cell$mo$WNN), c("wr", "wp"))
    expect_false("wr" %in% names(pDataDT(g, feat_type = "rna")))
    expect_error(runIntegratedUMAP(g, graph = "knn", network_name = "wnn15"),
        "larger integrated_k")
    g <- runIntegratedUMAP(g, graph = "knn", network_name = "wnn15", k = 15,
        name = "mo.umap")
    u <- getDimReduction(g, feat_type = "rna", reduction_method = "umap",
        name = "mo.umap", output = "matrix")
    nn <- .wnn_knn(g, "wnn15")
    GiottoUtils::local_seed(1234)
    ref <- uwot::umap(X = NULL, n_neighbors = 15,
        nn_method = GiottoClass::nnToUwot(list(id = nn$id, dist = nn$dist)),
        spread = 5, min_dist = 0.01)
    expect_identical(unname(u), unname(ref))
    g2 <- runIntegratedUMAP(g, graph = "knn", network_name = "wnn15", k = 15,
        name = "mo.umap")
    expect_identical(getDimReduction(g2, feat_type = "rna",
        reduction_method = "umap", name = "mo.umap", output = "matrix"), u)

    g3 <- runIntegratedUMAP(g, graph = "knn", network_name = "wnn15", k = 10,
        name = "mo.umap10", seed = 5, n_epochs = 50)
    GiottoUtils::local_seed(5)
    ref10 <- uwot::umap(X = NULL, n_neighbors = 10,
        nn_method = GiottoClass::nnToUwot(list(id = nn$id[, 1:10],
            dist = nn$dist[, 1:10])), spread = 5, min_dist = 0.01,
        n_epochs = 50)
    expect_identical(unname(getDimReduction(g3, feat_type = "rna",
        reduction_method = "umap", name = "mo.umap10", output = "matrix")),
        unname(ref10))

    expect_warning(runIntegratedUMAP(g, graph = "knn", network_name = "wnn15",
        k = 15, integrated_feat_type = "mo"), "integrated_feat_type")
    expect_no_warning(runIntegratedUMAP(g, graph = "knn",
        network_name = "wnn15", k = 15))
    expect_error(runIntegratedUMAP(fx$g, graph = "knn"), "run runWNN")
})

test_that("the jaccard sNN equals a dense Jaccard graph of the kNN", {
    fx <- .wnn_fixture()
    g <- .wnn_paper(fx$g, integrated_snn = "jaccard", integrated_k = 15)
    cells <- spatIDs(g)
    n <- length(cells)
    id <- .wnn_knn(g)$id
    A <- matrix(0, n, n)
    A[cbind(rep(seq_len(n), 15), as.vector(id))] <- 1
    S <- A %*% t(A)
    J <- S / (2 * 15 - S)
    keep <- upper.tri(J) & J >= 1 / 15
    ref <- data.table::data.table(a = cells[row(J)[keep]],
        b = cells[col(J)[keep]], w = J[keep])
    sdt <- getNearestNetwork(g, feat_type = "rna", nn_type = "sNN",
        name = "integrated_sNN", output = "data.table")
    got <- data.table::data.table(a = pmin(sdt$from, sdt$to),
        b = pmax(sdt$from, sdt$to), w = sdt$weight)
    ref[, `:=`(a2 = pmin(a, b), b2 = pmax(a, b))]
    m <- merge(got, ref[, .(a = a2, b = b2, w_ref = w)], by = c("a", "b"))
    expect_identical(nrow(m), nrow(ref))
    expect_identical(nrow(got), nrow(ref))
    expect_equal(m$w, m$w_ref, tolerance = 1e-14)
    expect_equal(sdt$distance, 1 - sdt$weight)
    expect_identical(getNearestNetwork(g, feat_type = "rna", nn_type = "sNN",
        name = "integrated_sNN")@misc$nn_input$prune, 1 / 15)

    g2 <- .wnn_paper(fx$g, integrated_snn = "jaccard", integrated_k = 15,
        snn_prune = 0.3, integrated_network_name = "wnn",
        integrated_snn_name = "wsnn")
    s2 <- getNearestNetwork(g2, feat_type = "rna", nn_type = "sNN",
        name = "wsnn", output = "data.table")
    expect_identical(nrow(s2), sum(keep & J >= 0.3))
    expect_true(all(s2$weight >= 0.3))
    g2 <- doLeidenCluster(g2, feat_type = "rna", nn_network_to_use = "sNN",
        network_name = "wsnn", name = "leiden_wsnn", resolution = 0.5)
    expect_true(length(unique(pDataDT(g2)$leiden_wsnn)) > 1L)
})

test_that("score_cap bounds the scores", {
    fx <- .wnn_fixture()
    w0 <- .wnn_weights(.wnn_paper(fx$g, score_cap = 0))
    expect_equal(unname(w0$w_rna), rep(0.5, length(w0$w_rna)))
    wi <- .wnn_weights(.wnn_paper(fx$g, score_cap = Inf))
    expect_true(all(is.finite(wi$w_rna)))
    expect_true(all(wi$w_rna >= 0 & wi$w_rna <= 1))
})

test_that("bad embeddings and inputs are refused", {
    fx <- .wnn_fixture()
    xp <- getDimReduction(fx$g, feat_type = "protein", name = "protein.pca",
        output = "matrix")
    xp[5, 3] <- NA
    xp[9, 1] <- Inf
    g <- setGiotto(fx$g, createDimObj(coordinates = xp,
        name = "protein.pca", spat_unit = "cell", feat_type = "protein",
        method = "pca"), verbose = FALSE)
    expect_error(.wnn_paper(g, dimensions_to_use = 1:8),
        "protein has NA or infinite values for 2 cells")
    expect_error(.wnn_paper(fx$g, dimensions_to_use = 10),
        "l2_norm needs at least 2 dimensions of rna \\(got 1\\)")
    expect_no_error(.wnn_paper(fx$g, dimensions_to_use = 10, l2_norm = FALSE))
    expect_error(.wnn_paper(fx$g, k = 1), "k must be at least 2")
    expect_error(.wnn_paper(fx$g, feat_types = "rna"))
    expect_error(.wnn_paper(fx$g, integrated_snn = "snn"))
    small <- .wnn_fixture(n = 15, knn_k = 5)
    expect_error(.wnn_paper(small$g, k = 15, knn_range = 20),
        "must be smaller than the number of cells")
})

test_that("reruns replace the results", {
    fx <- .wnn_fixture()
    a <- .wnn_paper(fx$g)
    b <- .wnn_paper(a)
    expect_identical(.wnn_weights(a), .wnn_weights(b))
    expect_identical(.wnn_knn(a), .wnn_knn(b))
    expect_identical(pDataDT(a, feat_type = "protein"),
        pDataDT(b, feat_type = "protein"))
})

test_that("a non-PCA reduction takes its default name", {
    fx <- .wnn_fixture()
    xp <- getDimReduction(fx$g, feat_type = "protein", name = "protein.pca",
        output = "matrix")
    g <- setGiotto(fx$g, createDimObj(coordinates = xp, name = "protein.lsi",
        spat_unit = "cell", feat_type = "protein", method = "lsi"),
        verbose = FALSE)
    a <- .wnn_paper(g, reduction_methods = c("pca", "lsi"),
        dimensions_to_use = list(1:20, 1:10))
    b <- .wnn_paper(fx$g, dimensions_to_use = list(1:20, 1:10))
    expect_identical(.wnn_weights(a), .wnn_weights(b))
    rec <- a@nn_network$cell$rna$kNN$integrated_kNN@misc$nn_input
    expect_identical(rec$reduction_names, c("pca", "protein.lsi"))
})

test_that("modalities are aligned by cell ID", {
    fx <- .wnn_fixture()
    ref <- .wnn_weights(.wnn_paper(fx$g))
    g <- fx$g
    xp <- getDimReduction(g, feat_type = "protein", name = "protein.pca",
        output = "matrix")
    o <- rev(seq_len(nrow(xp)))
    g <- setGiotto(g, createDimObj(coordinates = xp[o, ],
        name = "protein.pca", spat_unit = "cell", feat_type = "protein",
        method = "pca"), verbose = FALSE)
    g <- createNearestNetwork(g, feat_type = "protein", type = "kNN",
        k = 20, dimensions_to_use = 1:8, verbose = FALSE)
    w <- .wnn_weights(.wnn_paper(g))
    expect_equal(w$w_rna, ref$w_rna, tolerance = 1e-12)

    g <- setGiotto(fx$g, createDimObj(coordinates = xp[-(1:5), ],
        name = "protein.pca", spat_unit = "cell", feat_type = "protein",
        method = "pca"), verbose = FALSE)
    expect_error(.wnn_paper(g, dimensions_to_use = 1:8), "5 missing")
})

test_that("dimensions, networks and k are resolved and checked", {
    fx <- .wnn_fixture(knn_dims = 1:6)
    g <- .wnn_paper(fx$g)
    rec <- g@nn_network$cell$rna$kNN$integrated_kNN@misc$nn_input
    expect_identical(rec$dimensions_to_use, list(rna = 1:6, protein = 1:6))
    g <- .wnn_paper(fx$g, dimensions_to_use = list(1:12, 1:40))
    rec <- g@nn_network$cell$rna$kNN$integrated_kNN@misc$nn_input
    expect_identical(rec$dimensions_to_use, list(rna = 1:12, protein = 1:10))

    g <- .wnn_paper(fx$g, network_names = c("kNN.pca", "kNN.pca"), k = 10)
    expect_true(all(is.finite(.wnn_weights(g)$w_rna)))
    expect_error(.wnn_paper(fx$g, network_names = c("kNN.pca", "kNN.pca"),
        k = 25), "at least that many")
    expect_error(.wnn_paper(fx$g, reduction_names = c("rna.pca", "protein.pca")))
    expect_error(.wnn_paper(fx$g, integrated_k = 30, knn_range = 20),
        "cannot exceed")
    g <- .wnn_paper(fx$g, dimensions_to_use = list(protein = 1:3, rna = 1:5))
    rec <- g@nn_network$cell$rna$kNN$integrated_kNN@misc$nn_input
    expect_identical(rec$dimensions_to_use, list(rna = 1:5, protein = 1:3))
    expect_error(.wnn_paper(fx$g, dimensions_to_use = list(RNA = 1:3, protein = 1:3)),
        "\\[runWNN\\] the names of `dimensions_to_use`")
    expect_error(.wnn_paper(fx$g, dimensions_to_use = list(1:3)),
        "\\[runWNN\\] `dimensions_to_use` needs one entry")
})

test_that("runIntegratedUMAP refuses an integrated kNN made before subsetting", {
    fx <- .wnn_fixture()
    g <- .wnn_paper(fx$g)
    keep <- spatIDs(g)[seq(1, length(spatIDs(g)), by = 2)]
    gs <- subsetGiotto(g, feat_type = "rna", cell_ids = keep)
    expect_error(runIntegratedUMAP(gs, graph = "knn"),
        "rerun runWNN\\(\\) or\\s+runBlockPCA\\(\\) after subsetting")
    gs <- .wnn_paper(gs)
    gs <- runIntegratedUMAP(gs, graph = "knn")
    u <- getDimReduction(gs, feat_type = "rna", reduction_method = "umap",
        name = "integrated.umap", output = "matrix")
    expect_identical(rownames(u), spatIDs(gs))
})

test_that("the hnsw engine gives the same result twice", {
    skip_if_not_installed("RcppHNSW")
    fx <- .wnn_fixture()
    a <- .wnn_paper(fx$g, engine = "hnsw", n_threads = 2, knn_range = 100)
    b <- .wnn_paper(fx$g, engine = "hnsw", n_threads = 2, knn_range = 100)
    expect_identical(a@multiomics, b@multiomics)
    expect_identical(.wnn_knn(a), .wnn_knn(b))
    ref <- .wnn_paper(fx$g, knn_range = 100)
    expect_equal(.wnn_weights(a)$w_rna, .wnn_weights(ref)$w_rna,
        tolerance = 1e-6)
})

test_that("more than two modalities are integrated", {
    fx <- .wnn_fixture()
    g <- fx$g
    cells <- spatIDs(g)
    withr::local_seed(11)
    xa <- matrix(stats::rnorm(length(cells) * 5), ncol = 5,
        dimnames = list(cells, paste0("Dim.", 1:5)))
    ma <- Matrix::Matrix(matrix(stats::rpois(3 * length(cells), 5), 3,
        dimnames = list(paste0("a", 1:3), cells)), sparse = TRUE)
    g <- setGiotto(g, createExprObj(ma, name = "raw", spat_unit = "cell",
        feat_type = "atac"), verbose = FALSE)
    g <- setGiotto(g, createDimObj(coordinates = xa, name = "atac.pca",
        spat_unit = "cell", feat_type = "atac", method = "pca"),
        verbose = FALSE)
    g <- .wnn_paper(g, feat_types = c("rna", "protein", "atac"))
    mo <- g@multiomics$cell$rna_protein_atac$WNN
    w <- cbind(mo$w_rna[[1]], mo$w_protein[[1]], mo$w_atac[[1]])
    expect_true(all(is.finite(w)))
    expect_equal(unname(rowSums(w)), rep(1, nrow(w)))
    expect_lt(median(w[, 3]), 0.05)
})

test_that("runWNN works on a GiottoDisk-backed object", {
    skip_if_not_installed("GiottoDisk")
    vault <- withr::local_tempdir()
    fx_mem <- .wnn_fixture()
    fx <- .wnn_fixture(backend = file.path(vault, "v"))
    expect_false(is.null(fx$g@source))
    knn <- fx$g@nn_network$cell$protein$kNN$kNN.pca[]
    expect_true(inherits(knn, "dataStore"))
    g <- .wnn_paper(fx$g)
    net <- g@nn_network$cell$rna$kNN$integrated_kNN[]
    expect_true(inherits(net, "dataStore"))
    expect_true(inherits(g@nn_network$cell$rna$sNN$integrated_sNN[],
        "dataStore"))
    g_mem <- .wnn_paper(fx_mem$g)
    expect_equal(.wnn_weights(g)$w_rna, .wnn_weights(g_mem)$w_rna,
        tolerance = 1e-12)
    expect_identical(.wnn_knn(g)$id, .wnn_knn(g_mem)$id)
    g <- .wnn_paper(fx$g, integrated_snn = "jaccard")
    expect_true(inherits(g@nn_network$cell$rna$sNN$integrated_sNN[],
        "dataStore"))
    g <- runIntegratedUMAP(g, graph = "knn")
    g <- doLeidenCluster(g, feat_type = "rna", nn_network_to_use = "kNN",
        network_name = "integrated_kNN", name = "integrated_leiden",
        resolution = 0.5)
    expect_identical(length(unique(pDataDT(g)$integrated_leiden)) > 1L, TRUE)
})

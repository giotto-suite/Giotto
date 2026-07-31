# Tests for pcaParam("auto") + autoPcaParam routing. allMatrix defaults
# to IRLBA; substrates in other packages register their own method.
# dry_run = TRUE returns the resolved concrete pcaParam.

.tiny_dgc <- function(n_genes = 40L, n_cells = 200L, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.4,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", sprintf("%03d", seq_len(n_genes)))
    colnames(m) <- paste0("c", sprintf("%04d", seq_len(n_cells)))
    m
}


# Factory ####

test_that("pcaParam('auto') returns autoPcaParam with knobs populated", {
    p <- pcaParam("auto", ncp = 30, center = TRUE, scale = FALSE,
        set_seed = TRUE, seed_number = 42L)
    expect_s4_class(p, "autoPcaParam")
    expect_true(is(p, "pcaParam"))
    expect_equal(p$method, "auto")
    expect_equal(p$ncp, 30L)
    expect_equal(p$center, TRUE)
    expect_equal(p$scale, FALSE)
    expect_equal(p$seed_number, 42L)
    expect_equal(p$dry_run, FALSE)   # default
})

test_that("pcaParam('auto', dry_run = TRUE) sets the dry_run flag", {
    p <- pcaParam("auto", dry_run = TRUE)
    expect_equal(p$dry_run, TRUE)
})

test_that("pcaParam defaults to 'auto' (first match.arg choice)", {
    p <- pcaParam()
    expect_s4_class(p, "autoPcaParam")
})


# dry_run: substrate resolution without running PCA ####

test_that("reduceData(dgc, auto + dry_run) returns irlbaPcaParam", {
    m <- .tiny_dgc()
    resolved <- reduceData(m,
        pcaParam("auto", ncp = 5, dry_run = TRUE))
    expect_s4_class(resolved, "irlbaPcaParam")
    expect_false(is(resolved, "autoPcaParam"))
})

test_that("dry_run preserves user knobs through the rebuild", {
    m <- .tiny_dgc()
    resolved <- reduceData(m,
        pcaParam("auto", ncp = 15, center = FALSE, scale = FALSE,
            set_seed = TRUE, seed_number = 7L,
            n_oversamples = 25L, n_power_iter = 4L,
            dry_run = TRUE))
    expect_equal(resolved$ncp, 15L)
    expect_equal(resolved$center, FALSE)
    expect_equal(resolved$scale, FALSE)
    expect_equal(resolved$seed_number, 7L)
    expect_equal(resolved$n_oversamples, 25L)
    expect_equal(resolved$n_power_iter, 4L)
    # dry_run slot is stripped from the concrete param (not a valid arg)
    expect_null(resolved$dry_run)
})


# End-to-end: reduceData routes to concrete flavor ####

test_that("reduceData(dgc, autoPcaParam) matches reduceData(dgc, irlbaPcaParam)", {
    m <- .tiny_dgc(seed = 2)
    auto_res <- reduceData(m,
        pcaParam("auto", ncp = 5, center = TRUE, scale = FALSE,
            set_seed = TRUE, seed_number = 42L))
    irlba_res <- reduceData(m,
        pcaParam("irlba", ncp = 5, center = TRUE, scale = FALSE,
            set_seed = TRUE, seed_number = 42L))
    expect_equal(auto_res$d, irlba_res$d, tolerance = 1e-10)
    expect_equal(dim(auto_res$u), dim(irlba_res$u))
    expect_equal(dim(auto_res$v), dim(irlba_res$v))
})

test_that("reduceData(dense_matrix, autoPcaParam) works via allMatrix inheritance", {
    m <- as.matrix(.tiny_dgc(n_genes = 30, n_cells = 50, seed = 4))
    res <- reduceData(m, pcaParam("auto", ncp = 5,
        center = TRUE, scale = FALSE,
        set_seed = TRUE, seed_number = 42L))
    expect_named(res, c("u", "d", "v", "sdev", "eigenvalues"),
        ignore.order = TRUE)
    expect_equal(length(res$d), 5L)
})


# ANY fallback: unregistered substrates warn + route to IRLBA.

test_that("reduceData(ANY, autoPcaParam) warns and falls back to irlba", {
    setClass("noAutoRouting", representation(x = "matrix"))
    on.exit(removeClass("noAutoRouting"), add = TRUE)
    x <- new("noAutoRouting",
        x = as.matrix(.tiny_dgc(n_genes = 20, n_cells = 40, seed = 5)))

    expect_warning(
        res <- reduceData(x, pcaParam("auto", ncp = 5,
            center = TRUE, scale = FALSE,
            set_seed = TRUE, seed_number = 42L,
            dry_run = TRUE)),
        "no substrate-specific auto routing registered"
    )
    expect_s4_class(res, "irlbaPcaParam")
})

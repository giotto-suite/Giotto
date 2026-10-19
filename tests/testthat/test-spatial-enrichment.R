# Characterization tests for the sign-matrix enrichment methods.
#
# These pin what PAGE, rank and hypergeometric produce *today*, before the
# analyzeParam refactor moves the arithmetic behind a verb generic. Nothing
# here asserts that the numbers are right -- only that they do not change. The
# reference values were produced by the pre-refactor implementations on the
# fixture below; a diff in any of them means the refactor changed behaviour.
#
# Where a test pins a bug rather than a behaviour it says so, and is updated
# by the commit that fixes it.

skip_if_no_mini <- function() {
    skip_if_not_installed("GiottoData")
}

# 634 genes x 624 spots. `python_path = NA` skips the conda probe, which this
# object does not need -- nothing below touches a python module.
.enrich_fixture <- function(n_markers = 150L, seed = 1L) {
    gwith_options(list(giotto.use_conda = FALSE, giotto.verbose = FALSE), {
        g <- suppressMessages(
            GiottoData::loadGiottoMini("visium", python_path = NA)
        )
        genes <- rownames(getExpression(g, values = "normalized",
                                        output = "matrix"))
        set.seed(seed)
        types <- c("typeA", "typeB", "typeC")
        sm <- matrix(0L, nrow = length(genes), ncol = length(types),
                     dimnames = list(genes, types))
        for (j in seq_along(types)) {
            sm[sample(length(genes), n_markers), j] <- 1L
        }
        list(g = g, sm = sm, genes = genes)
    })
}


# Written-out references for PAGE and hypergeometric: the statistics as the
# methods define them, computed directly on a dense matrix, one cell at a
# time. They are the ground truth the implementations are checked against.

# fold change of each feature over its reference level: the mean of the
# de-logged values, logged back in the same base
.ref_fc <- function(x, base = 2) {
    x <- as.matrix(x)
    x - log(rowMeans(base^x - 1) + 1, base = base)
}

.ref_page <- function(x, sm, min_overlap = 5) {
    fc <- .ref_fc(x)
    types <- colnames(sm)[colSums(sm[intersect(rownames(sm), rownames(fc)), ,
        drop = FALSE] > 0) > min_overlap]
    sapply(types, function(t) {
        g <- intersect(rownames(sm)[sm[, t] > 0], rownames(fc))
        apply(fc, 2, function(v) (mean(v[g]) - mean(v)) * sqrt(length(g)) / sd(v))
    })
}

.ref_hyper <- function(x, sm, top = 5, fc = .ref_fc(x)) {
    u <- intersect(rownames(sm), rownames(fc))
    sapply(colnames(sm), function(t) {
        g <- intersect(rownames(sm)[sm[, t] > 0], u)
        apply(fc, 2, function(v) {
            on <- names(v)[v > stats::quantile(v, 1 - top / 100)]
            on <- intersect(on, u)
            -log10(stats::phyper(length(intersect(on, g)) - 1, length(g),
                length(u) - length(g), length(on), lower.tail = FALSE))
        })
    })
}

# score columns of a result table, in type order
.scores <- function(dt) {
    as.matrix(dt[, setdiff(names(dt), "cell_ID"), with = FALSE])
}

.ref_hyper_plain <- function(x, sm, top) {
    .ref_hyper(x, sm, top, fc = x - rowMeans(x))
}

test_that("runPAGEEnrich matches the written-out PAGE z-score", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    res <- runPAGEEnrich(f$g, sign_matrix = f$sm,
                         return_gobject = FALSE, verbose = FALSE)
    expect_s4_class(res, "spatEnrObj")
    dt <- res[]
    x <- getExpression(f$g, values = "normalized", output = "matrix")
    expect_identical(dt$cell_ID, colnames(x))
    expect_setequal(names(dt), c("cell_ID", "typeA", "typeB", "typeC"))
    expect_equal(.scores(dt), .ref_page(x, f$sm),
        tolerance = 1e-10, ignore_attr = TRUE)

    expect_identical(res@method, "PAGE")
    expect_identical(res@name, "PAGE")
    expect_setequal(names(res@misc), c(
        "expr_values_used", "reverse_log_scale", "logbase",
        "p_values_calculated", "output_enrichment_scores",
        "include_depletion"
    ))
})

test_that("PAGE p-values are the permutation null in closed form", {
    skip_if_no_mini()
    f <- .enrich_fixture()
    x <- getExpression(f$g, values = "normalized", output = "matrix")
    got <- analyzeData(x, enrichParam("PAGE", p_value = TRUE,
        verbose = FALSE), sign_matrix = f$sm)
    z <- .ref_page(x, f$sm)
    G <- nrow(x)
    k <- colSums(f$sm > 0)
    ref <- stats::pnorm(t(t(z) / sqrt((G - k) / G)), lower.tail = FALSE)
    expect_equal(.scores(got), ref, tolerance = 1e-10, ignore_attr = TRUE)

    # against brute-force permutations on a few cells: random sets of the
    # type's size, their z's mean and sd per cell, and the normal tail
    fc <- .ref_fc(x)[, 1:5]
    set.seed(3)
    zp <- replicate(20000, {
        g <- sample(nrow(fc), k[["typeA"]])
        (colMeans(fc[g, , drop = FALSE]) - colMeans(fc)) *
            sqrt(k[["typeA"]]) / apply(fc, 2, sd)
    })
    perm <- stats::pnorm(z[1:5, "typeA"], rowMeans(zp), apply(zp, 1, sd),
        lower.tail = FALSE)
    expect_equal(got$typeA[1:5], unname(perm), tolerance = 0.02)
})


test_that("runRankEnrich output is unchanged", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # `expression_values` must be given explicitly: the default vector cannot
    # survive this function's own match.arg call (see the bug test below).
    res <- runRankEnrich(f$g, sign_matrix = f$sm,
                         expression_values = "normalized",
                         return_gobject = FALSE)

    expect_s4_class(res, "spatEnrObj")
    dt <- res[]
    expect_identical(nrow(dt), 624L)
    expect_setequal(names(dt), c("cell_ID", "typeA", "typeB", "typeC"))
    expect_equal(head(dt$typeA, 3),
                 c(0.887582, 0.880433, 0.877802), tolerance = 1e-5)
    expect_equal(head(dt$typeB, 3),
                 c(0.882421, 0.883287, 0.875983), tolerance = 1e-5)
    expect_equal(head(dt$typeC, 3),
                 c(0.883911, 0.888095, 0.882159), tolerance = 1e-5)

    expect_identical(res@method, "rank")
    expect_identical(res@name, "rank")
})


test_that("runHyperGeometricEnrich matches the written-out test", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    res <- runHyperGeometricEnrich(f$g, sign_matrix = f$sm,
                                   return_gobject = FALSE)
    expect_s4_class(res, "spatEnrObj")
    dt <- res[]
    x <- getExpression(f$g, values = "normalized", output = "matrix")
    expect_identical(dt$cell_ID, colnames(x))
    expect_equal(.scores(dt), .ref_hyper(x, f$sm),
        tolerance = 1e-10, ignore_attr = TRUE)

    p <- runHyperGeometricEnrich(f$g, sign_matrix = f$sm, p_value = TRUE,
        return_gobject = FALSE)[]
    expect_equal(.scores(p), 10^-.ref_hyper(x, f$sm),
        tolerance = 1e-10, ignore_attr = TRUE)

    expect_identical(res@method, "hypergeometric")
    expect_identical(res@name, "hypergeometric")
})

test_that("each cell is cut at its own top percentage", {
    # the cut used to be recycled down the rows, comparing cells against
    # other cells' thresholds
    x <- matrix(c(1:20, 20:1, rep(5, 19), 30), nrow = 20,
        dimnames = list(paste0("g", 1:20), c("a", "b", "c")))
    sm <- matrix(0, 20, 2, dimnames = list(rownames(x), c("hi", "lo")))
    sm[16:20, "hi"] <- 1
    sm[1:5, "lo"] <- 1
    got <- analyzeData(x, enrichParam("hypergeometric", top_percentage = 25,
        reverse_log_scale = FALSE), sign_matrix = sm)
    expect_equal(.scores(got), .ref_hyper_plain(x, sm, 25),
        ignore_attr = TRUE)
    # cell a's top quarter is g16..g20, cell b's g1..g5
    expect_gt(got$hi[1], got$lo[1])
    expect_gt(got$lo[2], got$hi[2])
})


test_that("the enrichment result lands in the gobject", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    g2 <- runPAGEEnrich(f$g, sign_matrix = f$sm, verbose = FALSE)
    expect_s4_class(g2, "giotto")
    expect_true("PAGE" %in% list_spatial_enrichments_names(
        g2, spat_unit = "cell", feat_type = "rna"
    ))

    enr <- getSpatialEnrichment(g2, name = "PAGE", output = "spatEnrObj")
    direct <- runPAGEEnrich(f$g, sign_matrix = f$sm,
                            return_gobject = FALSE, verbose = FALSE)
    expect_equal(enr[], direct[])

    # every run appends one parameter-history entry
    expect_gt(length(g2@parameters), length(f$g@parameters))
    last <- g2@parameters[[length(g2@parameters)]]
    expect_identical(unname(last[["method used"]]), "PAGE")
})


test_that("runSpatialEnrich routes to the same result as the direct call", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    for (m in c("PAGE", "rank", "hypergeometric")) {
        via <- runSpatialEnrich(f$g, enrich_method = m, sign_matrix = f$sm,
            expression_values = "normalized",
            return_gobject = FALSE, verbose = FALSE)
        direct <- switch(m,
            PAGE = runPAGEEnrich(f$g, sign_matrix = f$sm,
                expression_values = "normalized",
                return_gobject = FALSE, verbose = FALSE),
            rank = runRankEnrich(f$g, sign_matrix = f$sm,
                expression_values = "normalized", return_gobject = FALSE),
            hypergeometric = runHyperGeometricEnrich(f$g, sign_matrix = f$sm,
                expression_values = "normalized", return_gobject = FALSE)
        )
        expect_equal(via[], direct[], info = m)
    }
})


# --- bug fixes ---------------------------------------------------------------

test_that("runRankEnrich works on its own default expression_values", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # `expression_values` defaults to c("normalized", "raw", "scaled",
    # "custom") but the choices were built as unique(c("normalized", "scaled",
    # "custom", expression_values)) -- the same four in a different order. Not
    # identical to the arg, so match.arg refused a length-4 arg and the
    # function could not be called without naming a value explicitly.
    res <- runRankEnrich(f$g, sign_matrix = f$sm, return_gobject = FALSE)
    expect_s4_class(res, "spatEnrObj")

    explicit <- runRankEnrich(f$g, sign_matrix = f$sm,
        expression_values = "normalized", return_gobject = FALSE)
    expect_equal(res[], explicit[])
})


test_that("runPAGEEnrich honours output_enrichment", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # The wrapper passed the literal c("original", "zscore") down to
    # .page_dt_method(), which match.arg'd it back to "original". The user's
    # choice was discarded, so PAGE always returned unscaled scores.
    orig <- runPAGEEnrich(f$g, sign_matrix = f$sm, output_enrichment = "original",
        return_gobject = FALSE, verbose = FALSE)[]
    zsc <- runPAGEEnrich(f$g, sign_matrix = f$sm, output_enrichment = "zscore",
        return_gobject = FALSE, verbose = FALSE)[]

    expect_false(isTRUE(all.equal(orig$typeA, zsc$typeA)))
    # "zscore" standardizes within cell type
    expect_equal(mean(zsc$typeA), 0, tolerance = 1e-8)
    expect_equal(stats::sd(zsc$typeA), 1, tolerance = 1e-8)
    # the default is still "original"
    expect_equal(
        runPAGEEnrich(f$g, sign_matrix = f$sm,
            return_gobject = FALSE, verbose = FALSE)[],
        orig
    )
})


test_that("runRankEnrich(p_value = TRUE) returns p-values", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # The permutation branch recursed into runRankEnrich() without
    # return_gobject = FALSE, so it got a giotto object back and then subset it
    # as a table. p_value = TRUE errored in fitdistrplus every time.
    res <- runRankEnrich(f$g, sign_matrix = f$sm, p_value = TRUE,
                         n_times = 20, return_gobject = FALSE)
    expect_s4_class(res, "spatEnrObj")

    dt <- res[]
    score_cols <- c("typeA", "typeB", "typeC")
    for (cl in score_cols) {
        expect_true(all(dt[[cl]] >= 0 & dt[[cl]] <= 1), info = cl)
    }
    # cell_ID must survive as an ID, not be swept into the gamma transform
    expect_type(dt$cell_ID, "character")
    expect_identical(nrow(dt), 624L)

    # and p_value = FALSE still gives the scores, unchanged
    plain <- runRankEnrich(f$g, sign_matrix = f$sm, return_gobject = FALSE)
    expect_equal(head(plain[]$typeA, 3),
                 c(0.887582, 0.880433, 0.877802), tolerance = 1e-5)
})


test_that("runSpatialEnrich forwards the PAGE-only arguments it accepts", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # min_overlap_genes, max_block and verbose are in runSpatialEnrich()'s
    # signature but were never passed on to runPAGEEnrich(). Setting them did
    # nothing, with no warning -- the router is a hand-written switch and
    # forgot three of its own formals.
    sm <- cbind(f$sm, typeD = 0L)
    set.seed(9)
    sm[sample(nrow(sm), 20), "typeD"] <- 1L

    via <- runSpatialEnrich(f$g, enrich_method = "PAGE", sign_matrix = sm,
        min_overlap_genes = 30, return_gobject = FALSE, verbose = FALSE)[]
    direct <- runPAGEEnrich(f$g, sign_matrix = sm,
        min_overlap_genes = 30, return_gobject = FALSE, verbose = FALSE)[]

    # typeD has 20 markers and is dropped by both
    expect_false("typeD" %in% names(via))
    expect_setequal(names(via), names(direct))
    expect_equal(via, direct)

    # and at the default threshold it is kept, so the test above is
    # discriminating rather than vacuous
    kept <- runSpatialEnrich(f$g, enrich_method = "PAGE", sign_matrix = sm,
        return_gobject = FALSE, verbose = FALSE)[]
    expect_true("typeD" %in% names(kept))
})


test_that("every runSpatialEnrich formal reaches a method", {
    # The router used to accept min_overlap_genes, max_block and verbose and
    # forward none of them. This asserts the property rather than the three
    # instances: every formal of the router, apart from the ones it consumes
    # itself, must be a formal of the method it dispatches to.
    router <- names(formals(runSpatialEnrich))
    own <- c("gobject", "enrich_method")

    targets <- list(
        PAGE = runPAGEEnrich,
        rank = runRankEnrich,
        hypergeometric = runHyperGeometricEnrich
    )
    # arguments that belong to exactly one method
    method_only <- c(
        min_overlap_genes = "PAGE", include_depletion = "PAGE",
        max_block = "PAGE", verbose = "PAGE",
        ties_method = "rank", rbp_p = "rank", num_agg = "rank",
        n_times = NA, top_percentage = "hypergeometric"
    )

    for (arg in setdiff(router, own)) {
        owner <- if (arg %in% names(method_only)) method_only[[arg]] else NULL
        if (is.null(owner)) {
            # shared: must be a formal of all three
            for (m in names(targets)) {
                expect_true(arg %in% names(formals(targets[[m]])),
                            info = paste(arg, "->", m))
            }
        } else if (!is.na(owner)) {
            expect_true(arg %in% names(formals(targets[[owner]])),
                        info = paste(arg, "->", owner))
        }
    }

    # and the router body actually names each one
    body_txt <- paste(deparse(body(runSpatialEnrich)), collapse = " ")
    for (arg in setdiff(router, own)) {
        expect_match(body_txt, arg, fixed = TRUE, info = arg)
    }
})


# --- the param family --------------------------------------------------------

test_that("enrichParam builds the right class with the right defaults", {
    expect_true(isVirtualClass("enrichParam"))
    for (cl in c("pageEnrichParam", "rankEnrichParam", "hyperEnrichParam")) {
        expect_true(extends(cl, "enrichParam"), info = cl)
        expect_true(extends(cl, "analyzeParam"), info = cl)
    }

    p <- enrichParam("PAGE")
    expect_s4_class(p, "pageEnrichParam")
    expect_identical(p$min_overlap_genes, 5)
    expect_false(p$include_depletion)
    expect_identical(p$output_enrichment, "original")

    r <- enrichParam("rank")
    expect_s4_class(r, "rankEnrichParam")
    expect_identical(r$ties_method, "average")
    expect_identical(r$rbp_p, 0.99)

    h <- enrichParam("hypergeometric")
    expect_s4_class(h, "hyperEnrichParam")
    expect_identical(h$top_percentage, 5)

    # method name is case-insensitive, and unknown methods are refused
    expect_s4_class(enrichParam("page"), "pageEnrichParam")
    expect_error(enrichParam("gsva"))
    # an inapplicable value is refused rather than ignored
    expect_error(enrichParam("PAGE", output_enrichment = "quantile"))
})


test_that("the result contract is enforced at the seam", {
    good <- data.table::data.table(cell_ID = c("a", "b"), tA = c(1, 2))
    expect_identical(.enrich_check_contract(good), good)

    expect_error(.enrich_check_contract(as.data.frame(good)), "data.table")
    expect_error(
        .enrich_check_contract(data.table::data.table(id = "a", tA = 1)),
        "cell_ID"
    )
    expect_error(
        .enrich_check_contract(data.table::data.table(cell_ID = 1L, tA = 1)),
        "must be character"
    )
    expect_error(
        .enrich_check_contract(data.table::data.table(cell_ID = "a")),
        "no cell-type score columns"
    )
    expect_error(
        .enrich_check_contract(
            data.table::data.table(cell_ID = "a", tA = "x")
        ),
        "must be numeric"
    )
})


test_that("another package could contribute an enrichment method", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # Stand-in for what GiottoDisk (or any extension) does: define a subclass
    # of the VIRTUAL parent outside this package and attach one method. If
    # this passes, the seam is real -- nothing between the wrapper and the
    # arithmetic names a specific engine.
    setClass("fakeEnrichParam", contains = "enrichParam")
    on.exit(removeClass("fakeEnrichParam"), add = TRUE)
    setMethod("analyzeData", signature(x = "ANY", param = "fakeEnrichParam"),
        function(x, param, ..., sign_matrix) {
            data.table::data.table(
                cell_ID = colnames(x),
                constant = rep(param$fill, ncol(x))
            )
        }
    )
    on.exit(removeMethod("analyzeData", c("ANY", "fakeEnrichParam")),
            add = TRUE)

    param <- new("fakeEnrichParam", param = list(fill = 7))
    out <- .enrich_run(
        gobject = f$g, param = param, sign_matrix = f$sm,
        method = "fake", name = "fake", return_gobject = FALSE
    )
    expect_s4_class(out$enrObj, "spatEnrObj")
    expect_identical(out$enrObj@method, "fake")
    expect_identical(unique(out$enrObj[]$constant), 7)

    # and it goes into the gobject through the same path
    g2 <- .enrich_run(
        gobject = f$g, param = param, sign_matrix = f$sm,
        method = "fake", name = "fake", return_gobject = TRUE
    )$gobject
    expect_true("fake" %in% list_spatial_enrichments_names(
        g2, spat_unit = "cell", feat_type = "rna"
    ))
})


test_that("the verb runs on a bare matrix, with no gobject", {
    skip_if_no_mini()
    f <- .enrich_fixture()
    m <- getExpression(f$g, values = "normalized", output = "matrix")

    # The point of moving the arithmetic onto analyzeData: it is testable
    # without an object, and the numbers match what the wrapper produces.
    for (spec in list(
        list(p = enrichParam("rank"), w = function() runRankEnrich(
            f$g, sign_matrix = f$sm, return_gobject = FALSE)[]),
        list(p = enrichParam("hypergeometric"), w = function()
            runHyperGeometricEnrich(
                f$g, sign_matrix = f$sm, return_gobject = FALSE)[])
    )) {
        direct <- analyzeData(m, spec$p, sign_matrix = f$sm)
        via <- data.table::as.data.table(spec$w())
        data.table::setcolorder(via, names(direct))
        expect_equal(as.data.frame(direct), as.data.frame(via),
                     info = class(spec$p))
    }
})


test_that("cells scored in parts against a shared reference match all at once", {
    skip_if_no_mini()
    f <- .enrich_fixture()
    m <- getExpression(f$g, values = "normalized", output = "matrix")
    ref <- log(Matrix::rowMeans(2^m - 1) + 1, base = 2)
    parts <- split(seq_len(ncol(m)), seq_len(ncol(m)) > 300)

    # what a disk-backed engine does: score each window of cells against the
    # reference of all of them
    for (p in list(
        enrichParam("PAGE", verbose = FALSE),
        enrichParam("PAGE", p_value = TRUE, verbose = FALSE),
        enrichParam("hypergeometric"),
        enrichParam("hypergeometric", p_value = TRUE)
    )) {
        whole <- analyzeData(m, p, sign_matrix = f$sm)
        pieces <- data.table::rbindlist(lapply(parts, function(j) {
            analyzeData(m[, j], p, sign_matrix = f$sm, reference = ref)
        }))
        expect_equal(pieces, whole, tolerance = 1e-12, info = class(p))
    }
    expect_error(
        analyzeData(m, enrichParam("PAGE"), sign_matrix = f$sm,
                    reference = unname(ref)),
        "named by feature"
    )
})

test_that("rank's reverse_log_scale and logbase are inert, and say so", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # Silent on the defaults -- an existing call is not doing anything wrong.
    expect_no_warning(
        base <- runRankEnrich(f$g, sign_matrix = f$sm, return_gobject = FALSE)
    )
    # Deprecated when passed explicitly.
    expect_warning(
        runRankEnrich(f$g, sign_matrix = f$sm, logbase = 2,
                      return_gobject = FALSE),
        "logbase"
    )
    expect_warning(
        runRankEnrich(f$g, sign_matrix = f$sm, reverse_log_scale = FALSE,
                      return_gobject = FALSE),
        "reverse_log_scale"
    )

    # And genuinely inert: the scores do not move.
    for (args in list(list(logbase = 10), list(reverse_log_scale = FALSE),
                      list(reverse_log_scale = TRUE, logbase = 10))) {
        got <- suppressWarnings(do.call(runRankEnrich, c(
            list(f$g, sign_matrix = f$sm, return_gobject = FALSE), args)))
        expect_equal(got[], base[], info = paste(names(args), collapse = ","))
    }
})

test_that("ranking is invariant to the transform those arguments would apply", {
    # The reason the arguments cannot work, asserted rather than claimed:
    # rank is invariant to any monotonic per-gene transform, so neither the
    # reverse-log step nor per-gene centring can move a single rank.
    skip_if_not_installed("sparseMatrixStats")
    set.seed(1)
    x <- matrix(rpois(600, 3) + runif(600), 20, 30)

    r_raw <- sparseMatrixStats::rowRanks(x, ties.method = "average")
    expect_identical(
        r_raw, sparseMatrixStats::rowRanks(2^x - 1, ties.method = "average"))
    expect_identical(
        r_raw, sparseMatrixStats::rowRanks(10^x - 1, ties.method = "average"))
    expect_identical(
        r_raw,
        sparseMatrixStats::rowRanks(x - rowMeans(x), ties.method = "average"))

    # and the second ranking, which is what the method actually reports
    expect_identical(
        sparseMatrixStats::colRanks(-r_raw),
        sparseMatrixStats::colRanks(
            -sparseMatrixStats::rowRanks(2^x - 1, ties.method = "average")))
})


test_that("the router does not fire rank's deprecation on its own behalf", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # runSpatialEnrich() has reverse_log_scale/logbase in its own signature and
    # used to forward them to every method, so routing to rank warned even when
    # the caller never mentioned them.
    expect_no_warning(
        runSpatialEnrich(f$g, enrich_method = "rank", sign_matrix = f$sm,
                         return_gobject = FALSE)
    )
    # but asking for them explicitly is still told
    expect_warning(
        runSpatialEnrich(f$g, enrich_method = "rank", sign_matrix = f$sm,
                         logbase = 10, return_gobject = FALSE),
        "logbase"
    )
    # and the other two methods, where they do work, stay quiet
    for (m in c("PAGE", "hypergeometric")) {
        expect_no_warning(
            runSpatialEnrich(f$g, enrich_method = m, sign_matrix = f$sm,
                             logbase = 2, return_gobject = FALSE,
                             verbose = FALSE),
            message = m
        )
    }
})

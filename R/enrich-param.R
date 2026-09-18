# enrichParam ####
#
# Sign-matrix enrichment as an `analyzeData` verb, following the `markersParam`
# arrangement in differential_expression.R: a VIRTUAL parent under
# `analyzeParam`, one concrete class per method, a factory per flavour, and
# `analyzeData` methods that dispatch on the *expression object* rather than on
# the giotto object. Dispatching there is what lets a disk-backed store arrive
# intact and route to its own streaming method, in GiottoDisk, with no change
# here.
#
# The three methods share a return contract -- see `.enrich_check_contract()`
# -- so `runSpatialEnrich()` no longer needs to know which one it called, and
# a fourth engine can be contributed from outside this package.

#' @title Sign-matrix enrichment params
#' @name enrich_param
#' @description
#' Parameter objects for the sign-matrix spatial enrichment methods. Pass one
#' to [analyzeData()] together with a sign matrix, or use the
#' [runPAGEEnrich()] / [runRankEnrich()] / [runHyperGeometricEnrich()]
#' wrappers, which build the param for you.
#'
#' Each method has its own factory because their parameters do not overlap:
#' `min_overlap_genes`, `include_depletion` and `max_block` are PAGE-only,
#' `ties_method`, `rbp_p` and `num_agg` are rank-only, and `top_percentage` is
#' hypergeometric-only. A single flat signature would accept all of them for
#' every method and silently ignore the inapplicable ones, which is the
#' behaviour this replaces.
#'
#' @param method character. `"PAGE"`, `"rank"` or `"hypergeometric"`.
#' @param ... method-specific parameters; see Details.
#' @details
#' Shared by all three:
#' \describe{
#'   \item{`reverse_log_scale`}{logical. undo a log transform before averaging
#'     (default `TRUE`).}
#'   \item{`logbase`}{numeric. base for `reverse_log_scale` (default 2).}
#'   \item{`output_enrichment`}{`"original"` (default) or `"zscore"`.}
#'   \item{`p_value`}{logical. return p-values instead of scores.}
#' }
#' PAGE: `min_overlap_genes` (5), `include_depletion` (`FALSE`), `n_times`
#' (1000), `max_block` (2e7), `verbose` (`TRUE`).
#'
#' rank: `ties_method` (`"average"`), `n_times` (1000), `rbp_p` (0.99),
#' `num_agg` (100).
#'
#' hypergeometric: `top_percentage` (5).
#' @md
#' @returns an [enrichParam-class]-inheriting object
#' @examples
#' p <- enrichParam("PAGE", min_overlap_genes = 10)
#' p$min_overlap_genes
#' @family feature set enrichment
#' @seealso [analyze_param], [runSpatialEnrich()]
NULL

#' @rdname enrich_param
#' @exportClass enrichParam
setClass("enrichParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname enrich_param
#' @exportClass pageEnrichParam
setClass("pageEnrichParam", contains = "enrichParam")

#' @rdname enrich_param
#' @exportClass rankEnrichParam
setClass("rankEnrichParam", contains = "enrichParam")

#' @rdname enrich_param
#' @exportClass hyperEnrichParam
setClass("hyperEnrichParam", contains = "enrichParam")


# param factories ####

#' @rdname enrich_param
#' @export
enrichParam <- function(method = "PAGE", ...) {
    method <- match.arg(tolower(method), c("page", "rank", "hypergeometric"))
    switch(method,
        "page" = .enrich_param_page(...),
        "rank" = .enrich_param_rank(...),
        "hypergeometric" = .enrich_param_hyper(...)
    )
}

# Fields shared by all three methods.
#' @keywords internal
#' @noRd
.enrich_param_common <- function(p) {
    p$reverse_log_scale <- p$reverse_log_scale %null% TRUE
    p$logbase <- as.numeric(p$logbase %null% 2)
    p$output_enrichment <- match.arg(
        p$output_enrichment %null% "original", c("original", "zscore")
    )
    p$p_value <- isTRUE(p$p_value)
    p
}

#' @keywords internal
#' @noRd
.enrich_param_page <- function(...) {
    p <- .enrich_param_common(new("pageEnrichParam", param = list(...)))
    p$min_overlap_genes <- as.numeric(p$min_overlap_genes %null% 5)
    p$include_depletion <- isTRUE(p$include_depletion)
    p$n_times <- as.numeric(p$n_times %null% 1000)
    p$max_block <- as.numeric(p$max_block %null% 20e6)
    p$verbose <- p$verbose %null% TRUE
    p
}

#' @keywords internal
#' @noRd
.enrich_param_rank <- function(...) {
    p <- .enrich_param_common(new("rankEnrichParam", param = list(...)))
    p$ties_method <- match.arg(
        p$ties_method %null% "average", c("average", "max")
    )
    p$n_times <- as.numeric(p$n_times %null% 1000)
    p$rbp_p <- as.numeric(p$rbp_p %null% 0.99)
    p$num_agg <- as.numeric(p$num_agg %null% 100)
    p
}

#' @keywords internal
#' @noRd
.enrich_param_hyper <- function(...) {
    p <- .enrich_param_common(new("hyperEnrichParam", param = list(...)))
    p$top_percentage <- as.numeric(p$top_percentage %null% 5)
    p
}


# result contract ####

# What every enrichment engine must return, and the only thing the wrappers
# and `createSpatEnrObj()` rely on: a data.table with a character `cell_ID`
# column and one numeric column per cell type. Anything a method wants to
# report beyond that travels as attr(x, "detail") and no shared code reads it.
#
# Checked rather than assumed so that a contributed engine fails at the seam,
# naming what it got wrong, instead of somewhere inside createSpatEnrObj().
#' @keywords internal
#' @noRd
.enrich_check_contract <- function(x, what = "enrichment method") {
    if (!data.table::is.data.table(x)) {
        stop("[", what, "] must return a data.table, not ",
            paste(class(x), collapse = "/"), call. = FALSE)
    }
    if (!"cell_ID" %in% names(x)) {
        stop("[", what, "] result has no `cell_ID` column.", call. = FALSE)
    }
    if (!is.character(x$cell_ID)) {
        stop("[", what, "] `cell_ID` must be character, not ",
            class(x$cell_ID)[[1L]], ".", call. = FALSE)
    }
    score_cols <- setdiff(names(x), "cell_ID")
    if (length(score_cols) == 0L) {
        stop("[", what, "] result has no cell-type score columns.",
            call. = FALSE)
    }
    not_num <- score_cols[
        !vapply(x[, score_cols, with = FALSE], is.numeric, logical(1L))
    ]
    if (length(not_num) > 0L) {
        stop("[", what, "] score columns must be numeric; these are not: ",
            paste(not_num, collapse = ", "), call. = FALSE)
    }
    x
}


# Topic definitions. These own the `@name`; the setMethod blocks below use
# `@rdname` so roxygen still emits each method's \alias{} -- a `@name` on a
# setMethod block suppresses it, and R CMD check then reports the method as
# undocumented. Same arrangement as markersParam in differential_expression.R.

#' @name enrich_page
#' @title PAGE enrichment
#' @md
#' @returns a `data.table` of `cell_ID` and one column per cell type
#' @family feature set enrichment
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' expr <- GiottoClass::getExpression(g, values = "normalized",
#'                                    output = "matrix")
#' feats <- rownames(expr)
#' sm <- makeSignMatrixPAGE(
#'     sign_names = c("typeA", "typeB"),
#'     sign_list = list(feats[1:150], feats[151:300])
#' )
#'
#' # the verb runs on a bare matrix, with no giotto object involved
#' res <- analyzeData(expr, enrichParam("PAGE"), sign_matrix = sm)
#' head(res)
NULL

#' @name enrich_rank
#' @title Rank enrichment
#' @md
#' @returns a `data.table` of `cell_ID` and one column per cell type
#' @family feature set enrichment
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' expr <- GiottoClass::getExpression(g, values = "normalized",
#'                                    output = "matrix")
#' feats <- rownames(expr)
#' sm <- makeSignMatrixPAGE(
#'     sign_names = c("typeA", "typeB"),
#'     sign_list = list(feats[1:150], feats[151:300])
#' )
#'
#' # the verb runs on a bare matrix, with no giotto object involved
#' res <- analyzeData(expr, enrichParam("rank"), sign_matrix = sm)
#' head(res)
NULL

#' @name enrich_hyper
#' @title Hypergeometric enrichment
#' @md
#' @returns a `data.table` of `cell_ID` and one column per cell type
#' @family feature set enrichment
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' expr <- GiottoClass::getExpression(g, values = "normalized",
#'                                    output = "matrix")
#' feats <- rownames(expr)
#' sm <- makeSignMatrixPAGE(
#'     sign_names = c("typeA", "typeB"),
#'     sign_list = list(feats[1:150], feats[151:300])
#' )
#'
#' # the verb runs on a bare matrix, with no giotto object involved
#' res <- analyzeData(expr, enrichParam("hypergeometric"), sign_matrix = sm)
#' head(res)
NULL


# PAGE ####

#' @title PAGE enrichment
#' @rdname enrich_page
#' @param x expression values. A `matrix`, a `Matrix`, or anything an
#'   engine registered against [pageEnrichParam-class] accepts.
#' @param param a [pageEnrichParam-class].
#' @inheritParams enrichment_params
#' @param ... unused. Declared because [analyzeData()] carries it, so an
#'   engine registered from another package can accept its own arguments.
#' @md
#' @returns a `data.table` of `cell_ID` and one column per cell type
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "pageEnrichParam"),
    function(x, param, ..., sign_matrix) {
        res <- .page_dt_method(
            sign_matrix = sign_matrix,
            expr_values = as.matrix(x),
            min_overlap_genes = param$min_overlap_genes,
            logbase = param$logbase,
            reverse_log_scale = param$reverse_log_scale,
            output_enrichment = param$output_enrichment,
            p_value = param$p_value,
            include_depletion = param$include_depletion,
            n_times = param$n_times,
            max_block = param$max_block,
            verbose = param$verbose
        )
        out <- data.table::as.data.table(res[["matrix"]])
        # PAGE's long-form table is the one method-specific extra any caller
        # has ever wanted; it rides along rather than widening the contract.
        data.table::setattr(out, "detail", res[["DT"]])
        .enrich_check_contract(out, "analyzeData(pageEnrichParam)")
    }
)


# rank ####

#' @title Rank enrichment
#' @rdname enrich_rank
#' @inheritParams enrich_page
#' @param param a [rankEnrichParam-class].
#' @md
#' @returns a `data.table` of `cell_ID` and one column per cell type
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "rankEnrichParam"),
    function(x, param, ..., sign_matrix) {
        .enrich_check_contract(
            .rank_dt_method(sign_matrix = sign_matrix, expr_values = x,
                            param = param),
            "analyzeData(rankEnrichParam)"
        )
    }
)

# The rank-biased-precision sweep, lifted verbatim out of runRankEnrich() so
# that it can be called on a matrix. The permutation branch now recurses into
# this function rather than back through the wrapper -- which is what made
# `p_value = TRUE` unreachable before.
#' @keywords internal
#' @noRd
.rank_dt_method <- function(sign_matrix, expr_values, param) {
    interGene <- intersect(rownames(sign_matrix), rownames(expr_values))
    if (length(interGene) < 100) {
        stop("Please check the gene numbers or names of scRNA-seq. The names
            of scRNA-seq should be consistent with spatial data.")
    }

    enrichment <- matrix(
        data = NA,
        nrow = dim(sign_matrix)[2],
        ncol = dim(expr_values)[2]
    )

    # No per-gene reference level is computed here, and `reverse_log_scale` /
    # `logbase` are deliberately unused.
    #
    # This method's statistic is a rank of a rank: genes are ranked across cells
    # and those ranks are then ranked within each cell. Ranking is invariant to
    # any monotonic per-gene transform, so `logbase^x - 1` and subtracting a
    # per-gene mean both leave every rank exactly where it was -- verified in
    # test-spatial-enrichment.R against `sparseMatrixStats::rowRanks()`.
    #
    # The previous code computed `log(rowMeans(logbase^expr - 1) + 1)` and never
    # read it. Beyond being dead, `logbase^expr` on a dgCMatrix returns dense:
    # roughly 460 MB allocated and discarded on a 377 x 151,782 store.

    # fold change and ranking
    ties_1 <- ties_2 <- param$ties_method
    if (param$ties_method == "max") {
        ties_1 <- "min"
        ties_2 <- "max"
    }
    # else ties_1=ties_2 is equal to random
    geneFold <- expr_values
    geneFold <- sparseMatrixStats::rowRanks(geneFold, ties.method = ties_1)
    rankFold <- t(sparseMatrixStats::colRanks(-geneFold, ties.method = ties_2))

    rownames(rankFold) <- rownames(expr_values)
    colnames(rankFold) <- colnames(expr_values)

    for (i in seq_len(dim(sign_matrix)[2])) {
        signames <- rownames(sign_matrix)[which(sign_matrix[, i] > 0)]
        interGene <- intersect(signames, rownames(rankFold))
        filterSig <- sign_matrix[interGene, ]
        filterRankFold <- rankFold[interGene, ]

        multiplyRank <- (filterRankFold * filterSig[, i])^(1 / 2)
        rpb <- (1.0 - param$rbp_p) * (param$rbp_p^(multiplyRank - 1))

        vectorX <- rep(NA, dim(filterRankFold)[2])
        for (j in seq_len(dim(filterRankFold)[2])) {
            toprpb <- sort(rpb[, j], decreasing = TRUE)
            vectorX[j] <- sum(toprpb[seq_len(param$num_agg)])
        }
        enrichment[i, ] <- vectorX
    }

    rownames(enrichment) <- colnames(sign_matrix)
    colnames(enrichment) <- colnames(rankFold)
    enrichment <- t(enrichment)

    if (param$output_enrichment == "zscore") {
        enrichment <- scale(enrichment)
    }

    enrichmentDT <- data.table::data.table(cell_ID = rownames(enrichment))
    enrichmentDT <- cbind(enrichmentDT, data.table::as.data.table(enrichment))

    if (isTRUE(param$p_value)) {
        # Suggests-only, and reached solely by this branch -- checked here
        # rather than at the top so the common path does not require it.
        # It went unguarded while this branch errored regardless; now that
        # it works, a missing install is the only thing that stops it.
        package_check(pkg_name = "fitdistrplus", repository = "CRAN")

        random_rank <- .do_rank_permutation(
            sc_gene = rownames(sign_matrix), n = param$n_times
        )
        no_p <- param
        no_p$p_value <- FALSE
        random_enr <- .rank_dt_method(
            sign_matrix = random_rank, expr_values = expr_values, param = no_p
        )
        # by name, not position: cell_ID is character and unlisting it into
        # the gamma fit is how this used to fail
        score_cols <- names(random_enr)[
            vapply(random_enr, is.numeric, logical(1L))
        ]
        background <- unlist(
            random_enr[, score_cols, with = FALSE], use.names = FALSE
        )
        fit.gamma <- fitdistrplus::fitdist(
            background, distr = "gamma", method = "mle"
        )
        own_cols <- names(enrichmentDT)[
            vapply(enrichmentDT, is.numeric, logical(1L))
        ]
        enrichmentDT[, (own_cols) := lapply(.SD, function(v) {
            stats::pgamma(
                v, fit.gamma$estimate[1], rate = fit.gamma$estimate[2],
                lower.tail = FALSE, log.p = FALSE
            )
        }), .SDcols = own_cols]
    }

    enrichmentDT[]
}


# hypergeometric ####

#' @title Hypergeometric enrichment
#' @rdname enrich_hyper
#' @inheritParams enrich_page
#' @param param a [hyperEnrichParam-class].
#' @md
#' @returns a `data.table` of `cell_ID` and one column per cell type
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "hyperEnrichParam"),
    function(x, param, ..., sign_matrix) {
        .enrich_check_contract(
            .hyper_dt_method(sign_matrix = sign_matrix, expr_values = x,
                             param = param),
            "analyzeData(hyperEnrichParam)"
        )
    }
)

# Lifted verbatim out of runHyperGeometricEnrich().
#' @keywords internal
#' @noRd
.hyper_dt_method <- function(sign_matrix, expr_values, param) {
    if (param$reverse_log_scale == TRUE) {
        expr_values <- param$logbase^expr_values - 1
    }

    interGene <- intersect(rownames(expr_values), rownames(sign_matrix))
    inter_sign_matrix <- sign_matrix[interGene, ]

    aveExp <- log2(2 * (Matrix::rowMeans(2^(expr_values - 1), dims = 1)) + 1)
    foldChange <- expr_values - aveExp

    top_q <- 1 - param$top_percentage / 100
    quantilecut <- apply(
        foldChange, 2, stats::quantile, probs = top_q, na.rm = TRUE
    )
    expbinary <- t_flex(1 * t_flex(foldChange > quantilecut))

    markerGenes <- rownames(inter_sign_matrix)
    expbinaryOverlap <- expbinary[markerGenes, ]
    total <- length(markerGenes)
    enrichment <- matrix(
        data = NA,
        nrow = dim(inter_sign_matrix)[2],
        ncol = dim(expbinaryOverlap)[2]
    )

    for (i in seq_len(dim(inter_sign_matrix)[2])) {
        signames <- rownames(inter_sign_matrix)[
            which(inter_sign_matrix[, i] == 1)
        ]
        vectorX <- NULL
        for (j in seq_len(dim(expbinaryOverlap)[2])) {
            cellsiggene <- names(expbinaryOverlap[
                which(expbinaryOverlap[, j] == 1), j
            ])
            x <- length(intersect(cellsiggene, signames))
            m <- length(signames)
            n <- total - m
            k <- length(intersect(cellsiggene, markerGenes))
            vectorX <- append(vectorX, 0 - log10(stats::phyper(
                x, m, n, k, log.p = FALSE, lower.tail = FALSE
            )))
        }
        enrichment[i, ] <- vectorX
    }

    rownames(enrichment) <- colnames(inter_sign_matrix)
    colnames(enrichment) <- colnames(expbinaryOverlap)
    enrichment <- t(enrichment)

    if (param$output_enrichment == "zscore") {
        enrichment <- scale(enrichment)
    }

    enrichmentDT <- data.table::data.table(cell_ID = rownames(enrichment))
    enrichmentDT <- cbind(enrichmentDT, data.table::as.data.table(enrichment))

    if (isTRUE(param$p_value)) {
        score_cols <- names(enrichmentDT)[
            vapply(enrichmentDT, is.numeric, logical(1L))
        ]
        enrichmentDT[, (score_cols) := lapply(.SD, function(v) 10^(-v)),
            .SDcols = score_cols]
    }

    enrichmentDT[]
}


# shared gobject plumbing ####

# Everything the four enrichment wrappers used to do around the arithmetic:
# resolve the defaults, fetch the expression object, run the verb, build the
# spatEnrObj, append the parameter-history entry and set it back. It existed
# four times, character for character apart from the method name and the misc
# fields, which is why the wrappers drifted apart (see the four fixes that
# precede this commit).
#
# `expression_values` is passed through to @misc and to the history entry
# unresolved, as it always was, while `values` is what was actually fetched.
#' @keywords internal
#' @noRd
.enrich_run <- function(gobject,
    param,
    sign_matrix,
    method,
    name,
    spat_unit = NULL,
    feat_type = NULL,
    values = "normalized",
    expression_values = values,
    misc = list(),
    history = character(),
    densify = FALSE,
    return_gobject = TRUE) {
    spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    expr_values <- getExpression(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = values, output = "exprObj"
    )
    mat <- expr_values[]
    # `densify` is about a sparse Matrix that a method cannot index into. A
    # disk-backed store is exempt: materializing it here is exactly what the
    # streaming method exists to avoid, and it would turn an out-of-core run
    # into an out-of-memory one.
    if (isTRUE(densify) && !inherits(mat, "dataStore")) {
        mat <- Matrix::as.matrix(mat)
    }

    # Dispatched on the expression object, not called directly: getExpression()
    # returns whatever the slot holds, so a disk-backed store arrives here
    # intact and routes to its own method.
    res <- analyzeData(x = mat, param = param, sign_matrix = sign_matrix)

    enrObj <- createSpatEnrObj(
        name = name,
        method = method,
        enrichment_data = res,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = expr_values@provenance,
        misc = misc
    )

    if (!isTRUE(return_gobject)) {
        return(list(enrObj = enrObj, detail = attr(res, "detail")))
    }

    spenr_names <- list_spatial_enrichments_names(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )
    if (name %in% spenr_names) {
        cat(name, " has already been used, will be overwritten")
    }

    parameters_list <- gobject@parameters
    update_name <- paste0(length(parameters_list), "_spatial_enrichment")
    parameters_list[[update_name]] <- history
    gobject@parameters <- parameters_list

    list(gobject = setGiotto(gobject, enrObj), detail = attr(res, "detail"))
}

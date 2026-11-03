# nucleus deconvolution ####
#
# Types every cell of a segmentation (e.g. H&E nuclei) from the spots it lies
# in: a classifier on per-cell features (e.g. image embeddings stored as
# expression of the cell spatial unit) is trained so that the mean prediction
# of each spot's cells explains the spot's expression through spatialDWLS's
# dampened WLS, starting from the spot's deconvolved proportions. The pieces
# are the bag model (bag-model.R) and spatialDWLS::dwls_design(); this file is
# the gobject plumbing.


#' @title Type cells from the spots they lie in
#' @name runNucleusDeconv
#' @description
#' Predict the cell type of every cell of a segmentation, such as nuclei on
#' an H&E image, from spot-level data: per-cell features, the spots' expression
#' and a cell type signature. A classifier on the features is trained so that
#' the average prediction over the cells in a spot explains the spot's
#' expression the way spatialDWLS fits it (dampened weighted least squares),
#' and stays close to the spot's deconvolved proportions. It then types every
#' cell, also those between spots.
#' @inheritParams signature_analysis_params
#' @inheritParams deconvolution_params
#' @param spat_unit spatial unit (e.g. "cell"). Here the cells to type, e.g.
#'   `"nucleus"`. They need spatial locations (e.g. polygon centroids from
#'   [addSpatialCentroidLocations()]) and their features as an expression
#'   matrix, features by cells.
#' @param feat_type feature type (e.g. "rna", "dna", "protein"). Here the
#'   feature type of the cell features, e.g. `"img"`; `NULL` takes the one
#'   they are stored under.
#' @param expression_values character. Which expression values to use, e.g.
#'   "normalized". A method's own default is shown in its Usage section. Here
#'   the cell features; `NULL` takes the first of `spat_unit` and `feat_type`.
#' @param spot_unit spatial unit of the spots. It needs polygons of the same
#'   name, the expression, the cluster column and the deconvolution.
#' @param spot_feat_type feature type of the spot expression
#' @param spot_values expression values of the spots, log-scaled
#' @param logbase numeric. Log base of the expression values (default = 2).
#'   Here of `spot_values`.
#' @param cluster_column character. Cell metadata column holding the cluster
#'   assignment used to group cells. Here the spot clusters of `spot_unit`;
#'   spatialDWLS chooses its dampening constant per cluster.
#' @param deconv_name spatial enrichment of `spot_unit` with the starting cell
#'   type proportions (e.g. from [runDWLSDeconv()], or any method imported
#'   with [createSpatEnrObj()]). Its columns must include the cell types of
#'   `sign_matrix`. `NULL` trains from the expression alone, which can
#'   collapse onto a few cell types on real data.
#' @param loss_weights named numeric: weights of the `dwls` (expression),
#'   `proportion` and `confidence` losses.
#' @param warm integer. epochs that use only the proportion loss before
#'   switching to the expression loss; 0 (default) uses both throughout.
#' @param learner character. the classifier; only `"linear"` (softmax
#'   regression) for now.
#' @param epochs,lr,batch_size,weight_decay training: passes over the
#'   training spots, Adam step size, spots per batch, L2 penalty.
#' @param seed integer. seed for the spot split, initialisation and batch
#'   order; the session's random state is left as it was.
#' @param name character. Name to store the result under in the giotto
#'   object's spatial enrichment slot. `NULL` (default) uses the method's own
#'   name -- see the Usage section. It is also the prefix of the metadata
#'   columns (`"nucleusDeconv"` when `NULL`).
#' @param verbose be verbose
#' @returns
#' With `return_gobject = TRUE`, the giotto object with:
#' * a spatial enrichment `name` on `spat_unit`/`feat_type`: the probability
#'   of each cell type for every cell. Its `misc` holds the fitted model
#'   (`misc$model`), so it can be applied to the features of other cells.
#' * cell metadata columns `<name>_type` (most probable type),
#'   `<name>_confidence` (its probability) and `<name>_spot` (the spot the
#'   cell lies in, `NA` between spots).
#'
#' Otherwise the spatial enrichment object.
#' @details
#' Spots are split 80/10/10 into training, validation and test spots; the
#' epoch with the lowest validation loss is kept. Only the cells inside
#' spots are used to train.
#'
#' The losses: `dwls` is [spatialDWLS::dwls_expression_loss()] on the mean
#' prediction of each spot's cells; `proportion` the mean squared difference
#' between that mean and the starting proportions; `confidence` the mean of
#' `-log(max p)` over the cells, which makes each cell choose a type.
#' @family spatial deconvolution
#' @seealso [hierarchicalCellTypes()] to report a lineage and a confident
#'   subtype; [calculateLabelProportions()] for the spot composition from the
#'   typed cells.
#' @md
#' @export
runNucleusDeconv <- function(gobject,
        spat_unit,
        feat_type = NULL,
        expression_values = NULL,
        spot_unit = NULL,
        spot_feat_type = NULL,
        spot_values = "normalized",
        logbase = 2,
        cluster_column = "leiden_clus",
        sign_matrix,
        deconv_name = "DWLS",
        loss_weights = c(dwls = 1, proportion = 0.99, confidence = 0.01),
        warm = 0L,
        learner = "linear",
        epochs = 100L,
        lr = 1e-3,
        batch_size = 64L,
        weight_decay = 0,
        seed = 1L,
        name = NULL,
        return_gobject = TRUE,
        verbose = NULL) {
    # NSE vars
    cell_ID <- spot <- NULL

    fname <- "[runNucleusDeconv]"
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_character(spat_unit, len = 1L)
    if (is.null(name)) name <- "nucleusDeconv"
    checkmate::assert_character(name, len = 1L)
    checkmate::assert_numeric(loss_weights, lower = 0, names = "named")
    if (!all(names(loss_weights) %in% c("dwls", "proportion", "confidence"))) {
        stop(fname, " `loss_weights` names must be among dwls, proportion ",
            "and confidence.", call. = FALSE)
    }
    w <- c(dwls = 0, proportion = 0, confidence = 0)
    w[names(loss_weights)] <- loss_weights
    if (w[["proportion"]] > 0 && is.null(deconv_name)) w[["proportion"]] <- 0
    if (warm > 0L && (w[["proportion"]] == 0 || w[["dwls"]] == 0)) {
        stop(fname, " `warm` needs both the proportion and the dwls loss.",
            call. = FALSE)
    }
    if (w[["dwls"]] == 0 && w[["proportion"]] == 0) {
        stop(fname, " no loss on the spots: give `deconv_name` or a ",
            "`dwls` weight.", call. = FALSE)
    }

    # the cells' feature type: the one their features are stored under, not
    # the object's default (which is usually the spots' "rna")
    feat_type <- feat_type %null% .nd_first_expr(gobject, spat_unit)$feat_type
    spot_unit <- set_default_spat_unit(gobject, spat_unit = spot_unit)
    spot_feat_type <- set_default_feat_type(gobject,
        spat_unit = spot_unit, feat_type = spot_feat_type)
    if (identical(spot_unit, spat_unit)) {
        stop(fname, " `spot_unit` and `spat_unit` are the same spatial unit.",
            call. = FALSE)
    }
    types <- colnames(sign_matrix)
    if (is.null(types) || anyDuplicated(types)) {
        stop(fname, " `sign_matrix` needs unique cell types as column names.",
            call. = FALSE)
    }

    # groups: the spot each cell's location falls in
    groups <- .nd_groups(gobject, spot_unit = spot_unit, spat_unit = spat_unit)
    vmsg(.v = verbose, fname, sprintf("%d cells in %d spots",
        data.table::uniqueN(groups$cell_ID), data.table::uniqueN(groups$group)))

    # losses
    losses <- list()
    if (w[["dwls"]] > 0) {
        cm <- pDataDT(gobject, spat_unit = spot_unit, feat_type = spot_feat_type)
        if (!cluster_column %in% colnames(cm)) {
            stop(fname, " cluster column '", cluster_column, "' not found in ",
                "the cell metadata of '", spot_unit, "'.", call. = FALSE)
        }
        spot_expr <- getExpression(gobject, spat_unit = spot_unit,
            feat_type = spot_feat_type, values = spot_values,
            output = "matrix")
        .dwls_engine_check()
        design <- spatialDWLS::dwls_design(spot_expr, sign_matrix,
            stats::setNames(as.character(cm[[cluster_column]]), cm$cell_ID),
            logbase = logbase)
        losses$dwls <- .bag_loss("dwls", design = design,
            weight = w[["dwls"]], phase = if (warm > 0L) "main" else "all")
    }
    if (w[["proportion"]] > 0) {
        losses$proportion <- .bag_loss("proportion",
            target = .nd_start_props(gobject, spot_unit, spot_feat_type,
                deconv_name, types),
            weight = w[["proportion"]],
            phase = if (warm > 0L) "warm" else "all")
    }
    if (w[["confidence"]] > 0) {
        losses$confidence <- .bag_loss("confidence",
            weight = w[["confidence"]])
    }

    # features of the cells in spots: train
    ex <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values %null%
            .nd_first_expr(gobject, spat_unit, feat_type)$name,
        output = "exprObj")
    feats <- ex[]
    in_spot <- intersect(colnames(feats), groups$cell_ID)
    X <- t(as.matrix(feats[, in_spot, drop = FALSE]))
    model <- .bag_train(X, groups, losses,
        learner = .bag_learner(learner), types = types,
        epochs = epochs, lr = lr, batch_size = batch_size, wd = weight_decay,
        warm = warm, seed = seed, verbose = verbose)
    rm(X)

    # all cells
    P <- .nd_predict(model, feats)
    res <- data.table::data.table(cell_ID = rownames(P), P)
    enr <- createSpatEnrObj(res,
        name = name,
        method = "nucleusDeconv",
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = prov(ex),
        misc = list(
            model = model,
            spot_unit = spot_unit,
            deconv_name = deconv_name,
            loss_weights = w,
            warm = warm,
            learner = learner,
            epochs = epochs,
            lr = lr,
            seed = seed
        ),
        verbose = FALSE
    )
    if (!isTRUE(return_gobject)) return(enr)

    gobject <- setGiotto(gobject, enr, verbose = FALSE)
    spot_of <- groups[!duplicated(cell_ID)]
    meta <- data.table::data.table(
        cell_ID = rownames(P),
        type = types[max.col(P, ties.method = "first")],
        confidence = apply(P, 1L, max)
    )
    meta[, spot := spot_of$group[match(cell_ID, spot_of$cell_ID)]]
    data.table::setnames(meta, c("type", "confidence", "spot"),
        paste0(name, c("_type", "_confidence", "_spot")))
    gobject <- addCellMetadata(gobject,
        spat_unit = spat_unit, feat_type = feat_type,
        new_metadata = meta, by_column = TRUE, column_cell_ID = "cell_ID")

    parameters_list <- gobject@parameters
    update_name <- paste0(length(parameters_list), "_nucleus_deconvolution")
    parameters_list[[update_name]] <- c(
        "spatial unit" = spat_unit,
        "feature type" = feat_type,
        "spot unit" = spot_unit,
        "starting proportions" = deconv_name %null% "none",
        "loss weights" = paste(names(w), w, sep = "=", collapse = ", "),
        "learner" = learner,
        "best epoch" = model$best_epoch,
        "name" = name
    )
    gobject@parameters <- parameters_list
    gobject
}


#' @title Lineage and confident subtype from cell type probabilities
#' @name hierarchicalCellTypes
#' @description
#' Read out cell type probabilities, from [runNucleusDeconv()] or any spatial
#' enrichment with one probability column per cell type, at two levels: the
#' lineage is the one with the highest summed probability of its types, and
#' the subtype is reported only when it holds at least `tau` of that
#' lineage's probability. Cells whose subtype is not confident keep their
#' lineage.
#' @inheritParams signature_analysis_params
#' @param feat_type feature type (e.g. "rna", "dna", "protein"). `NULL` takes
#'   the one the enrichment is stored under.
#' @param enrichment_name name of the spatial enrichment
#' @param lineage_map named character: the lineage of each cell type, named by
#'   cell type. Types missing from it are their own lineage.
#' @param tau numeric in \[0, 1\]. share of the lineage's probability the top
#'   subtype needs.
#' @param column_name character. the metadata column to write; also the
#'   prefix of the others.
#' @returns the giotto object with cell metadata columns `<column_name>`
#'   (the subtype when confident, the lineage otherwise),
#'   `<column_name>_lineage`, `<column_name>_subtype` (`NA` when not
#'   confident) and `<column_name>_share`; or that table.
#' @family spatial deconvolution
#' @md
#' @export
hierarchicalCellTypes <- function(gobject,
        spat_unit = NULL,
        feat_type = NULL,
        enrichment_name = "nucleusDeconv",
        lineage_map,
        tau = 0.7,
        column_name = "cell_type_hier",
        return_gobject = TRUE) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_character(lineage_map, names = "unique")
    checkmate::assert_number(tau, lower = 0, upper = 1)
    spat_unit <- set_default_spat_unit(gobject, spat_unit = spat_unit)
    if (is.null(feat_type)) {
        # the feature type the enrichment is stored under on this unit
        le <- list_spatial_enrichments(gobject, spat_unit = spat_unit)
        le <- le[le$name == enrichment_name]
        feat_type <- if (nrow(le)) le$feat_type[[1]] else
            set_default_feat_type(gobject, spat_unit = spat_unit)
    }
    enr <- getSpatialEnrichment(gobject, spat_unit = spat_unit,
        feat_type = feat_type, name = enrichment_name, output = "data.table")
    P <- as.matrix(enr[, setdiff(colnames(enr), "cell_ID"), with = FALSE])
    res <- .hier_read_out(P, lineage_map, tau)
    res <- data.table::data.table(cell_ID = enr$cell_ID, res)
    data.table::setnames(res, c("label", "lineage", "subtype", "share"),
        c(column_name, paste0(column_name, c("_lineage", "_subtype", "_share"))))
    if (!isTRUE(return_gobject)) return(res)
    addCellMetadata(gobject, spat_unit = spat_unit, feat_type = feat_type,
        new_metadata = res, by_column = TRUE, column_cell_ID = "cell_ID")
}

# label / lineage / subtype / share for each row of P (cells x types)
.hier_read_out <- function(P, lineage_map, tau) {
    types <- colnames(P)
    lin <- lineage_map[types]
    lin[is.na(lin)] <- types[is.na(lin)]
    lev <- unique(lin)
    LS <- vapply(lev, function(l) rowSums(P[, lin == l, drop = FALSE]),
        numeric(nrow(P)))
    LS <- matrix(LS, nrow(P), dimnames = list(NULL, lev))
    lpred <- lev[max.col(LS, ties.method = "first")]
    within <- P / pmax(LS[, lin, drop = FALSE], 1e-12)
    within[lin[col(P)] != lpred[row(P)]] <- -1
    top <- max.col(within, ties.method = "first")
    share <- within[cbind(seq_len(nrow(P)), top)]
    subtype <- ifelse(share >= tau, types[top], NA_character_)
    data.frame(
        label = ifelse(is.na(subtype), lpred, subtype),
        lineage = lpred,
        subtype = subtype,
        share = share,
        stringsAsFactors = FALSE
    )
}

# Cells of `spat_unit` (as points at their spatial locations) in the polygons
# of `spot_unit`: a table of group (spot ID) and cell_ID. The same relation
# calculateLabelProportions(group_method = "polygon") uses.
.nd_groups <- function(gobject, spot_unit, spat_unit) {
    # NSE vars
    group <- cell_ID <- NULL

    spots <- getPolygonInfo(gobject, name = spot_unit,
        return_giottoPolygon = TRUE, verbose = FALSE)
    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        output = "spatLocsObj", verbose = FALSE)
    rel <- relate(spots, createGiottoPoints(as.points(sl)),
        relation = "intersects", pairs = TRUE, output = "data.table",
        use_names = TRUE)
    data.table::setnames(rel, c("x", "y"), c("group", "cell_ID"))
    rel[, group := as.character(group)]
    rel[, cell_ID := as.character(cell_ID)]
    rel
}

# Starting proportions (spots x types): the deconvolution's cell types of
# `types`, spots without any proportion dropped, rows renormalised.
.nd_start_props <- function(gobject, spot_unit, spot_feat_type, deconv_name,
        types) {
    enr <- getSpatialEnrichment(gobject, spat_unit = spot_unit,
        feat_type = spot_feat_type, name = deconv_name, output = "data.table")
    miss <- setdiff(types, colnames(enr))
    if (length(miss)) {
        stop("[runNucleusDeconv] the deconvolution '", deconv_name,
            "' has no column for cell types ", paste(miss, collapse = ", "),
            ".", call. = FALSE)
    }
    P <- as.matrix(enr[, types, with = FALSE])
    rownames(P) <- as.character(enr$cell_ID)
    P[is.na(P)] <- 0
    P <- P[rowSums(P) > 0, , drop = FALSE]
    P / rowSums(P)
}

# The first expression of `spat_unit` (and `feat_type`, when given): a list of
# its feat_type and name.
.nd_first_expr <- function(gobject, spat_unit, feat_type = NULL) {
    ex <- list_expression(gobject, spat_unit = spat_unit)
    if (!is.null(ex) && !is.null(feat_type)) ex <- ex[ex$feat_type == feat_type]
    if (is.null(ex) || nrow(ex) == 0L) {
        stop("[runNucleusDeconv] no features (expression) for spatial unit '",
            spat_unit, "'", if (!is.null(feat_type)) paste0(
                " and feature type '", feat_type, "'"), ".", call. = FALSE)
    }
    list(feat_type = ex$feat_type[[1]], name = ex$name[[1]])
}

# Predictions for every cell of a features x cells matrix, densified in
# chunks of cells.
.nd_predict <- function(model, feats, chunk_size = 20000L) {
    n <- ncol(feats)
    starts <- seq(1L, max(n, 1L), by = chunk_size)
    do.call(rbind, lapply(starts, function(s) {
        i <- s:min(s + chunk_size - 1L, n)
        .bag_predict(model, t(as.matrix(feats[, i, drop = FALSE])))
    }))
}

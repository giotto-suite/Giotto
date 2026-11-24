# Plots for an edgeicg icfObject (findInteractionChangedFeats(method = "edgeicg")).
# The classic plotICF() family draws the classic test's columns; these read the
# edgeicg fit kept in icfObject$edgeicg: scales, hit classes, the estimated leak,
# per-cell exposures and the niche attribution.

# hit classes in a fixed order and colour (see edgeicg::edgeicg_classify())
.icf_class_colors <- c(
    clean = "#1b1b1b", sender_marker = "#d95f02",
    other_marker = "#7570b3", own_down = "#1b9e77"
)

.icf_edgeicg_fit <- function(icfObject) {
    if (!inherits(icfObject, "icfObject") || is.null(icfObject$edgeicg)) {
        .gstop("icfObject must come from",
            "findInteractionChangedFeats(method = \"edgeicg\")")
    }
    icfObject$edgeicg
}

.icf_type_order <- function(fit, cell_order) {
    lv <- fit$types[fit$fitted]
    if (is.null(cell_order)) return(lv)
    c(intersect(cell_order, lv), setdiff(lv, cell_order))
}

.icf_theme <- function() {
    ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title.position = "plot"
        )
}

#' @title Plots for edgeicg interaction-changed features
#' @name plotICF_edgeicg
#' @description Views of an `icfObject` from
#' [findInteractionChangedFeats()] with `method = "edgeicg"`:
#' * `plotICFOverview()`: hits per receiver and sender type at each scale
#'   (dot size: number of hits; colour: share going up).
#' * `plotICFLeak()`: the estimated transcript leak between touching cells
#'   (`lambda`: the share of a touching sender's expression found in the
#'   receiver), which edgeicg removes before testing.
#' * `plotICFPair()`: all features tested for one receiver and sender at one
#'   scale, effect size against significance, coloured by hit class.
#' * `plotICFExposure()`: receiver expression (share of the cell's raw counts)
#'   against exposure to the sender, in exposure bins; the first bin holds the
#'   unexposed receivers.
#' * `plotICFAttribution()`: for niche hits, how much of the evidence each
#'   sender type explains (Wakefield attribution with a 95% credible set).
#' * `plotICFRegion()`: for a fit with `region_column`, the effect of each
#'   feature within each region (log2 FC with a 95% interval; filled =
#'   significant in that region), titled with the between-region pattern.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param icfObject output of `findInteractionChangedFeats(method = "edgeicg")`
#' @param scale one or more of `"contact"`, `"paracrine"`, `"niche"`
#' @param q_thr,min_lfc significance and effect-size thresholds for a hit
#'   (default: the `min_lfc` of the icfObject)
#' @param classes hit classes to count (`"clean"`, `"sender_marker"`,
#'   `"other_marker"`, `"own_down"`); default `"clean"`
#' @param include_regional count niche hits flagged as regional
#' @param cell_order optional order of cell types on the axes
#' @param receiver,sender receiver and sender cell types
#' @param n_labels features labelled per class
#' @param feats features to show
#' @param n_bins number of exposure bins for the exposed receivers
#' @param regions optional regions to show, in this order
#' @returns ggplot
#' @seealso [findInteractionChangedFeats()], [icfParam()]
#' @md
NULL

#' @rdname plotICF_edgeicg
#' @export
plotICFOverview <- function(gobject, icfObject,
    scale = c("contact", "paracrine", "niche"), q_thr = 0.05, min_lfc = NULL,
    classes = "clean", include_regional = FALSE, cell_order = NULL,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotICFOverview") {
    fit <- .icf_edgeicg_fit(icfObject)
    scale <- match.arg(scale, several.ok = TRUE)
    min_lfc <- min_lfc %null% icfObject$test_info$min_lfc %null% 0.25
    if (q_thr > icfObject$test_info$report_q) {
        warning("ICFscores only hold entries with q < ",
            icfObject$test_info$report_q, call. = FALSE)
    }
    # NSE variables
    p.adj <- log2fc <- class <- regional <- cell_type <- int_cell_type <-
        n <- up <- NULL
    scales <- scale
    sc <- data.table::as.data.table(icfObject$ICFscores)
    sc <- sc[sc$scale %in% scales & p.adj < q_thr & abs(log2fc) >= min_lfc &
        class %in% classes]
    if (!include_regional) sc <- sc[!sc$regional]
    lv <- .icf_type_order(fit, cell_order)
    ov <- sc[, list(n = .N, up = mean(log2fc > 0)),
        by = list(receiver = cell_type, sender = int_cell_type, scale)]
    ov[, `:=`(
        receiver = factor(receiver, rev(lv)), sender = factor(sender, lv),
        scale = factor(scale, c("contact", "paracrine", "niche"))
    )]
    pl <- ggplot2::ggplot(ov, ggplot2::aes(sender, receiver)) +
        ggplot2::geom_point(ggplot2::aes(size = n, colour = up)) +
        ggplot2::facet_wrap(~scale, drop = FALSE) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        ggplot2::scale_size_area(max_size = 7, name = "hits") +
        ggplot2::scale_colour_gradient2(
            low = "#2166ac", mid = "#bdbdbd", high = "#b2182b",
            midpoint = 0.5, limits = c(0, 1), name = "share up"
        ) +
        ggplot2::labs(
            x = "sender (neighbour type)",
            y = "receiver (type whose expression changes)",
            title = "Interaction-changed features per cell-type pair",
            subtitle = sprintf("q < %g, |log2 FC| >= %g, class: %s%s",
                q_thr, min_lfc, paste(classes, collapse = ", "),
                if (include_regional) "" else ", not regional")
        ) +
        .icf_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 45, hjust = 1))
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

#' @rdname plotICF_edgeicg
#' @export
plotICFLeak <- function(gobject, icfObject, cell_order = NULL,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotICFLeak") {
    fit <- .icf_edgeicg_fit(icfObject)
    lambda <- NULL
    lv <- .icf_type_order(fit, cell_order)
    pr <- data.table::data.table(
        receiver = factor(fit$types[fit$pairs$receiver], rev(lv)),
        sender = factor(fit$types[fit$pairs$sender], lv),
        lambda = fit$pairs$lambda
    )[!is.na(lambda)]
    pl <- ggplot2::ggplot(pr, ggplot2::aes(sender, receiver, fill = lambda)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
        ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b",
            name = "leak\n(lambda)") +
        ggplot2::labs(x = "sender (touching neighbour)", y = "receiver",
            title = "Estimated transcript leak between touching cells") +
        .icf_theme() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank()
        )
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

#' @rdname plotICF_edgeicg
#' @export
plotICFPair <- function(gobject, icfObject, receiver, sender,
    scale = c("contact", "paracrine", "niche"), q_thr = 0.05, n_labels = 4L,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotICFPair") {
    fit <- .icf_edgeicg_fit(icfObject)
    scale <- match.arg(scale)
    q <- lfc <- class <- gene <- NULL
    # every tested feature of this pair, not only the reported hits
    r <- edgeicg::edgeicg_results(fit, scale, q_thr = 1, min_lfc = 0)
    rcv <- receiver
    snd <- sender # the table has receiver and sender columns
    r <- r[r$receiver == rcv & r$sender == snd]
    if (!nrow(r)) {
        .gstop("no tests for receiver", receiver, "and sender", sender,
            "at", scale, "scale")
    }
    r <- edgeicg::edgeicg_classify(fit, r)
    r[, class := factor(class, names(.icf_class_colors))]
    lab <- r[q < q_thr][order(q)][, utils::head(.SD, n_labels), by = class]
    text_geom <- if (requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene),
            size = 3, min.segment.length = 0, seed = 1, show.legend = FALSE)
    } else {
        ggplot2::geom_text(data = lab, ggplot2::aes(label = gene), size = 3,
            vjust = -0.7, show.legend = FALSE)
    }
    pl <- ggplot2::ggplot(r, ggplot2::aes(lfc, -log10(q), colour = class)) +
        ggplot2::geom_point(size = 1.4, alpha = 0.8) +
        text_geom +
        ggplot2::geom_hline(yintercept = -log10(q_thr), linetype = 2,
            linewidth = 0.3) +
        ggplot2::scale_colour_manual(values = .icf_class_colors, drop = FALSE) +
        ggplot2::labs(
            x = "log2 FC (receivers at high vs low exposure)",
            y = "-log10 q",
            title = sprintf("%s next to %s (%s scale)", receiver, sender,
                scale)
        ) +
        .icf_theme()
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

#' @rdname plotICF_edgeicg
#' @export
plotICFExposure <- function(gobject, icfObject, receiver, sender, feats,
    scale = c("niche", "paracrine", "contact"), n_bins = 8L,
    spat_unit = NULL, feat_type = NULL,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotICFExposure") {
    fit <- .icf_edgeicg_fit(icfObject)
    scale <- match.arg(scale)
    spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )
    x <- y <- b <- se <- feat <- NULL
    cluster <- icfObject$Giotto_info$cluster
    values <- icfObject$Giotto_info$values
    expr <- getExpression(gobject, spat_unit = spat_unit,
        feat_type = feat_type, values = values, output = "exprObj")[]
    ids <- colnames(expr)
    meta <- getCellMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, output = "data.table", copy_obj = TRUE)
    cell_type <- as.character(meta[[cluster]][match(ids, meta$cell_ID)])
    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        output = "data.table", copy_obj = TRUE)
    coords <- as.matrix(sl[match(ids, sl$cell_ID), c("sdimx", "sdimy")])
    # library size over all features, as in the fit
    if ("total_expr" %in% names(meta)) {
        lib <- meta$total_expr[match(ids, meta$cell_ID)]
    } else {
        cs <- analyzeData(expr, analyzeParam("cell_stats"))
        lib <- cs[["total_expr"]][match(ids, cs[["cells"]])]
    }
    # exposures with the settings of the fit
    es <- fit$settings$exposures
    tc <- match(cell_type, fit$types)
    ex <- edgeicg::edgeicg_exposures(coords, tc, length(fit$types),
        k = es$k, sigma_p = es$sigma_p, sigma_n = es$sigma_n, bin = es$bin)
    ri <- which(cell_type == receiver)
    xe <- switch(scale, contact = ex$C, paracrine = ex$P, niche = ex$N)[
        ri, match(sender, fit$types)]
    feats <- intersect(feats, rownames(expr))
    if (!length(feats)) .gstop("none of the feats are in the expression")
    sub <- expr[feats, ids[ri]]
    if (!inherits(sub, c("matrix", "Matrix"))) {
        sub <- GiottoDisk::storeRead(sub, output = "dgcmatrix",
            max_rows = Inf, max_cols = Inf)[feats, ids[ri], drop = FALSE]
    }
    sub <- as.matrix(sub)
    bins <- integer(length(ri))
    exposed <- xe > 1e-6
    if (any(exposed)) {
        bins[exposed] <- cut(rank(xe[exposed], ties.method = "first"),
            min(n_bins, sum(exposed)), labels = FALSE)
    }
    d <- data.table::rbindlist(lapply(feats, function(f) {
        data.table::data.table(x = xe, y = sub[f, ] / lib[ri] * 100, b = bins,
            feat = f)[, list(x = mean(x), y = mean(y),
                se = stats::sd(y) / sqrt(.N)), by = list(feat, b)]
    }))
    d[, feat := factor(feat, feats)]
    pl <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = y - se, ymax = y + se),
            width = 0, linewidth = 0.4, colour = "grey40") +
        ggplot2::geom_line(linewidth = 0.5, colour = "grey40") +
        ggplot2::geom_point(size = 2, colour = "#1b1b1b") +
        ggplot2::facet_wrap(~feat, scales = "free_y") +
        ggplot2::labs(
            x = sprintf("%s exposure to %s (share of the neighbourhood)",
                scale, sender),
            y = "expression (% of the cell's counts)",
            title = sprintf("%s: expression against exposure to %s",
                receiver, sender),
            subtitle = "bin means +/- SE; first point: unexposed receivers"
        ) +
        .icf_theme()
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

#' @rdname plotICF_edgeicg
#' @export
plotICFAttribution <- function(gobject, icfObject, receiver, feats,
    cell_order = NULL,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotICFAttribution") {
    fit <- .icf_edgeicg_fit(icfObject)
    gene <- sender <- alpha <- in_cs <- NULL
    at <- data.table::as.data.table(fit$attribution)
    if (!nrow(at)) .gstop("the fit has no niche attribution")
    ri <- match(receiver, fit$types) # the table has a receiver column
    gs <- feats
    at <- at[at$receiver == ri & at$gene %in% gs]
    if (!nrow(at)) {
        .gstop("no attributed niche hits for these feats in", receiver)
    }
    lv <- .icf_type_order(fit, cell_order)
    at[, `:=`(
        sender = factor(fit$types[sender], rev(lv)),
        gene = factor(gene, intersect(feats, gene)),
        in_cs = factor(ifelse(in_cs, "yes", "no"), c("yes", "no"))
    )]
    pl <- ggplot2::ggplot(at, ggplot2::aes(alpha, sender, fill = in_cs)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::facet_wrap(~gene, nrow = 1) +
        ggplot2::scale_fill_manual(values = c(yes = "#1b1b1b", no = "#bdbdbd"),
            name = "in 95%\ncredible set", drop = FALSE) +
        ggplot2::scale_x_continuous(limits = c(0, 1),
            breaks = c(0, 0.5, 1)) +
        ggplot2::labs(x = "attribution weight", y = "sender",
            title = sprintf("Which neighbour type explains the niche effect in %s?",
                receiver)) +
        .icf_theme()
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

#' @rdname plotICF_edgeicg
#' @export
plotICFRegion <- function(gobject, icfObject, receiver, sender, feats,
    scale = c("contact", "paracrine", "niche"), regions = NULL,
    q_thr = 0.05, min_lfc = 0.25,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotICFRegion") {
    fit <- .icf_edgeicg_fit(icfObject)
    if (is.null(fit$region)) {
        .gstop("the fit has no region effects; rerun",
            "findInteractionChangedFeats() with region_column")
    }
    scale <- match.arg(scale)
    # local copies first: the table has columns named like these arguments
    rcv <- receiver
    snd <- sender
    fs <- feats
    sc <- scale
    region <- lfc <- se <- effect <- q <- z <- pattern <- panel <- NULL
    lo <- hi <- sig <- NULL
    r <- .icf_region_table(fit, q_thr = q_thr, min_lfc = min_lfc, all = TRUE)
    r <- r[r$cell_type == rcv & r$int_cell_type == snd & r$scale == sc &
        r$feats %in% fs]
    if (!nrow(r)) {
        .gstop("no region tests for", receiver, "next to", sender, "at",
            scale, "scale for these features")
    }
    lv <- if (is.null(regions)) fit$region$labels else regions
    r <- r[r$region %in% lv]
    r[, region := factor(region, rev(lv))]
    # log2 FC and its interval on the same scale (lfc is linear in the effect
    # for small effects)
    k <- ifelse(abs(r$effect) > 0, r$lfc / r$effect, NA_real_)
    r[, `:=`(lo = lfc - 1.96 * se * abs(k), hi = lfc + 1.96 * se * abs(k),
        sig = !is.na(q) & q < q_thr & abs(lfc) >= min_lfc)]
    r[, panel := sprintf("%s (%s)", feats, pattern)]
    r[, panel := factor(panel, unique(panel[match(fs, feats)]))]
    pl <- ggplot2::ggplot(r, ggplot2::aes(lfc, region)) +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi),
            height = 0, linewidth = 0.5, colour = "grey40") +
        ggplot2::geom_point(ggplot2::aes(shape = sig), size = 2.4,
            colour = "#1b1b1b", fill = "#d95f02") +
        ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
            labels = c(`TRUE` = "significant", `FALSE` = "not significant"),
            name = NULL) +
        ggplot2::facet_wrap(~panel, scales = "free_x") +
        ggplot2::labs(
            x = "log2 FC within the region (high vs low exposure)", y = NULL,
            title = sprintf("%s next to %s by region (%s scale)", receiver,
                sender, scale)
        ) +
        .icf_theme()
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

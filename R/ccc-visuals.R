# Plots for an lrmap result (spatCellCellcom(method = "lrmap")). The classic
# plotCCcom*() family draws one scale of it; these read the lrmap fit kept in
# attr(result, "lrmap"): calls per cell-type pair and scale, one pair on the
# tissue (per-cell local score, lrmap::lrmap_local()), and expression-only
# versus spatial calls.

.ccc_lrmap_fit <- function(cccResult) {
    fit <- attr(cccResult, "lrmap")
    if (is.null(fit)) {
        .gstop("cccResult must be the result of",
            "spatCellCellcom(method = \"lrmap\")")
    }
    fit
}

#' @title Plots for lrmap ligand-receptor co-location
#' @name plotLRmap
#' @description Views of the result of [spatCellCellcom()] with
#' `method = "lrmap"`:
#' * `plotLRmapOverview()`: called ligand-receptor pairs per sender and
#'   receiver cell type at each scale.
#' * `plotLRmapPair()`: one pair on the tissue. `view = "tissue"` colours
#'   every receiver cell by its local score (receptor x ligand it is exposed
#'   to from sender cells at the scale); `view = "zoom"` shows a window with
#'   the ligand in sender cells (circles, orange) and the receptor in receiver
#'   cells (squares, blue). Expression before leak removal.
#' * `plotLRmapCompare()`: for a result with the `"expression"` scale, how
#'   the expression-only calls (mean ligand x mean receptor, no spatial
#'   information) split into spatially supported, co-located cell types
#'   without pair support, and cell types that never co-locate.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param cccResult output of `spatCellCellcom(method = "lrmap")`
#' @param scale scales to show; for `plotLRmapPair()` one scale (default: the
#'   spatial scale where the pair is most significant)
#' @param test which call rule (see [cccParam()]): `"type"` (default),
#'   `"shuffle"` or `"colocation"`
#' @param q_thr,min_fold q-value and fold thresholds for a call
#' @param cell_order optional order of cell types on the axes
#' @param sender,receiver sender (ligand) and receiver (receptor) cell types
#' @param ligand,receptor ligand and receptor names as in the result
#' @param view `"tissue"` or `"zoom"`
#' @param zoom_size side of the zoom window (coordinate units, microns)
#' @param zoom_center optional centre (x, y) of the zoom window; default: the
#'   window with the highest summed local score
#' @param point_size point size (default from the number of cells)
#' @returns ggplot
#' @seealso [spatCellCellcom()], [cccParam()]
#' @md
NULL

#' @rdname plotLRmap
#' @export
plotLRmapOverview <- function(gobject, cccResult,
    scale = c("contact", "paracrine", "niche"),
    test = c("type", "shuffle", "colocation"), q_thr = 0.05,
    min_fold = 1.25, cell_order = NULL,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotLRmapOverview") {
    fit <- .ccc_lrmap_fit(cccResult)
    scale <- match.arg(scale, several.ok = TRUE)
    test <- match.arg(test)
    scale <- intersect(scale, fit$settings$scales)
    # NSE variables
    sender <- receiver <- n <- NULL
    r <- data.table::as.data.table(lrmap::lrmap_results(fit, scale,
        test = test, q_thr = q_thr, min_fold = min_fold))
    ov <- r[, list(n = .N), by = list(sender, receiver, scale)]
    lv <- fit$types
    if (!is.null(cell_order)) lv <- c(intersect(cell_order, lv), setdiff(lv, cell_order))
    ov[, `:=`(
        sender = factor(sender, rev(lv)), receiver = factor(receiver, lv),
        scale = factor(scale, c("contact", "paracrine", "niche"))
    )]
    pl <- ggplot2::ggplot(ov, ggplot2::aes(receiver, sender, fill = n)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = n), size = 2.8) +
        ggplot2::facet_wrap(~scale, drop = FALSE) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        ggplot2::scale_fill_gradient(low = "#d6e4f5", high = "#2a6fbb",
            name = "pairs") +
        ggplot2::labs(
            x = "receiver (receptor)", y = "sender (ligand)",
            title = "Ligand-receptor pairs called per cell-type pair",
            subtitle = sprintf("test: %s, q < %g, fold >= %g", test, q_thr,
                min_fold)
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

#' @rdname plotLRmap
#' @export
plotLRmapPair <- function(gobject, cccResult, sender, receiver, ligand,
    receptor, scale = NULL, view = c("tissue", "zoom"), zoom_size = 500,
    zoom_center = NULL, point_size = NULL, spat_unit = NULL, feat_type = NULL,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotLRmapPair") {
    fit <- .ccc_lrmap_fit(cccResult)
    view <- match.arg(view)
    info <- attr(cccResult, "ccc_info")
    spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )
    tb <- data.table::as.data.table(fit$table)
    rows <- tb[tb$sender == sender & tb$receiver == receiver &
        tb$ligand == ligand & tb$receptor == receptor &
        tb$scale != "expression"]
    if (!nrow(rows)) {
        .gstop("this sender, receiver and pair were not tested")
    }
    if (is.null(scale)) scale <- rows$scale[which.min(rows$q_type)]
    row <- rows[rows$scale == scale]
    if (!nrow(row)) .gstop("scale", scale, "is not in the result for this pair")
    # subunit genes from the ligand-receptor table used for the fit
    lr <- as.data.frame(attr(cccResult, "lr"))
    genes_of <- function(name, col, sub) {
        i <- match(name, lr[[col]])
        v <- if (sub %in% names(lr) && !is.na(i)) lr[[sub]][i] else NA
        if (is.na(v) || !nzchar(v)) name else strsplit(v, "_", fixed = TRUE)[[1]]
    }
    lg <- genes_of(ligand, "ligand", "ligand_subunits")
    rg <- genes_of(receptor, "receptor", "receptor_subunits")
    inp <- .ccc_inputs(gobject, spat_unit, feat_type, info$cluster,
        info$values)
    ids <- inp$ids
    st <- fit$settings
    loc <- lrmap::lrmap_local(.ccc_counts(inp$expr, ids),
        inp$coords[ids, , drop = FALSE], inp$cell_type[ids], sender, receiver,
        lg, rg, scale, lib = inp$lib[ids], k = st$k, sigma_p = st$sigma_p,
        sigma_n = st$sigma_n, bin_p = st$bin_p, bin_n = st$bin_n)
    d <- data.table::data.table(x = inp$coords[ids, 1], y = inp$coords[ids, 2],
        role = ifelse(loc$type == sender, "sender",
            ifelse(loc$type == receiver, "receiver", "other")),
        ligand = loc$ligand, receptor = loc$receptor, local = loc$local)
    # NSE variables
    x <- y <- role <- local <- col <- s <- NULL
    ttl <- sprintf("%s → %s:  %s → %s (%s)", sender, receiver,
        ligand, receptor, scale)
    sub <- sprintf("beyond cell types: %.2f×, q = %s   co-location: %.2f×, q = %s",
        row$fold_type, format(signif(row$q_type, 2)), row$fold,
        format(signif(row$q, 2)))
    ps <- point_size %null% if (view == "tissue") {
        max(0.05, min(1, 3000 / nrow(d)))
    } else 1.6
    if (view == "tissue") {
        dr <- d[role == "receiver"][order(local)]
        pl <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
            ggplot2::geom_point(data = d[role == "other"], colour = "#d9d8d3",
                size = ps * 0.6) +
            ggplot2::geom_point(data = d[role == "sender"], colour = "#f3b99c",
                size = ps * 0.8) +
            ggplot2::geom_point(data = dr, ggplot2::aes(colour = local),
                size = ps) +
            ggplot2::scale_colour_gradient(low = "#d6e4f5", high = "#0d3b75",
                name = "local\nscore") +
            ggplot2::coord_fixed() +
            ggplot2::labs(title = ttl, subtitle = paste0(sub,
                "\nreceivers coloured by receptor × ligand exposure; ",
                "senders light orange"), x = NULL, y = NULL) +
            .icf_theme()
    } else {
        if (is.null(zoom_center)) {
            dr <- d[role == "receiver" & !is.na(local)]
            k <- dr[, list(s = sum(local)), by = list(
                gx = floor(x / zoom_size), gy = floor(y / zoom_size))][
                which.max(s)]
            zx <- c(k$gx, k$gx + 1) * zoom_size; zy <- c(k$gy, k$gy + 1) * zoom_size
        } else {
            zx <- zoom_center[1] + c(-0.5, 0.5) * zoom_size
            zy <- zoom_center[2] + c(-0.5, 0.5) * zoom_size
        }
        z <- d[x >= zx[1] & x < zx[2] & y >= zy[1] & y < zy[2]]
        ramp <- function(v, lo, hi) {
            v <- pmin(v, stats::quantile(v[v > 0], 0.95, na.rm = TRUE))
            m <- max(v, 1e-12, na.rm = TRUE)
            grDevices::rgb(grDevices::colorRamp(c(lo, hi))(pmin(v / m, 1)),
                maxColorValue = 255)
        }
        z[role == "sender", col := ramp(ligand, "#fbe0d3", "#a8380d")]
        z[role == "receiver", col := ramp(receptor, "#d6e4f5", "#0d3b75")]
        pl <- ggplot2::ggplot(z, ggplot2::aes(x, y)) +
            ggplot2::geom_point(data = z[role == "other"], colour = "#d9d8d3",
                size = ps * 0.4) +
            ggplot2::geom_point(data = z[role == "sender"],
                ggplot2::aes(fill = col), shape = 21, colour = "white",
                stroke = 0.2, size = ps) +
            ggplot2::geom_point(data = z[role == "receiver"],
                ggplot2::aes(fill = col), shape = 22, colour = "white",
                stroke = 0.2, size = ps) +
            ggplot2::scale_fill_identity() +
            ggplot2::coord_fixed(xlim = zx, ylim = zy, expand = FALSE) +
            ggplot2::labs(title = ttl, subtitle = sprintf(paste0("%s\n",
                "orange circles: %s in %s; blue squares: %s in %s ",
                "(light = 0, dark = high)"), sub, ligand, sender, receptor,
                receiver), x = NULL, y = NULL) +
            .icf_theme()
    }
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

#' @rdname plotLRmap
#' @export
plotLRmapCompare <- function(gobject, cccResult,
    test = c("type", "shuffle", "colocation"), q_thr = 0.05, min_fold = 1.25,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(), default_save_name = "plotLRmapCompare") {
    fit <- .ccc_lrmap_fit(cccResult)
    test <- match.arg(test)
    if (!"expression" %in% fit$settings$scales) {
        .gstop("the result has no expression scale; run spatCellCellcom()",
            "with scales = c(\"contact\", \"paracrine\", \"niche\",",
            "\"expression\")")
    }
    cm <- data.table::as.data.table(lrmap::lrmap_compare(fit, test = test,
        q_thr = q_thr, min_fold = min_fold))
    # NSE variables
    expr_call <- spatial_call <- types_colocated <- group <- n <- share <-
        NULL
    lv <- c("spatially supported", "cell types co-locate, pair not beyond",
        "cell types never co-locate")
    ex <- cm[expr_call == TRUE]
    ex[, group := ifelse(spatial_call, lv[1],
        ifelse(types_colocated, lv[2], lv[3]))]
    ov <- ex[, list(n = .N), by = group][, group := factor(group, rev(lv))]
    ov[, share := n / sum(n)]
    pl <- ggplot2::ggplot(ov, ggplot2::aes(n, group, fill = group)) +
        ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%d (%.0f%%)", n,
            100 * share)), hjust = -0.1, size = 3.2) +
        ggplot2::scale_fill_manual(values = stats::setNames(
            c("#2a6fbb", "#9cc0ec", "#c8541f"), lv)) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0.25))) +
        ggplot2::labs(x = "expression-only calls", y = NULL,
            title = "Expression-only calls and spatial support",
            subtitle = sprintf(paste0("%d of %d tests called on expression ",
                "alone; spatial test: %s"), nrow(ex), nrow(cm), test)) +
        .icf_theme()
    GiottoVisuals::plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}

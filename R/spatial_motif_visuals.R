# Plots for multicellular motif results.
#
# All of these consume the motif result contract -- topology, color_tuple,
# observed, z, p_adj -- and nothing engine-specific. In particular none of them
# parses `motif_id`: its grammar belongs to whichever engine produced it, so a
# plot that took it apart would quietly tie these figures to smotif. Topology
# and colours are read from their own columns.


# helpers ####

# GiottoVisuals' .motif_theme() is internal to that package, so rather than reach
# in with ::: this mirrors its sizing (8pt legend / axis text, white panel) as
# a local theme. Keeps these figures looking like the rest of the suite without
# depending on another package's private API.
#' @keywords internal
#' @noRd
.motif_theme <- function(...) {
    base <- ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = 8),
            legend.title = ggplot2::element_text(size = 8),
            axis.title = ggplot2::element_text(size = 8),
            axis.text = ggplot2::element_text(size = 8),
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 8, hjust = 0),
            strip.background = ggplot2::element_rect(
                fill = "grey95", colour = NA
            ),
            panel.grid.minor = ggplot2::element_blank()
        )
    # overrides go in a second theme() layer; folding them into the first would
    # pass the same formal twice when a caller blanks something set above
    extra <- list(...)
    if (length(extra)) base <- base + do.call(ggplot2::theme, extra)
    base
}

#' @keywords internal
#' @noRd
.motif_require_contract <- function(x, what = "result") {
    if (!inherits(x, "data.frame")) {
        .gstop(sprintf("`%s` must be a motif result data.table", what))
    }
    need <- c("motif_id", "topology", "color_tuple", "observed", "z", "p_adj")
    miss <- setdiff(need, colnames(x))
    if (length(miss)) {
        .gstop(sprintf(
            "`%s` is missing column(s): %s. Expected the output of %s.",
            what, paste(miss, collapse = ", "), "cellProximityMotifs()"
        ))
    }
    invisible(TRUE)
}

#' @keywords internal
#' @noRd
.motif_label <- function(res) {
    vapply(res$color_tuple, paste, character(1L), collapse = "-")
}

# Vertex coordinates and edges for each topology, drawn in canonical slot
# order so position 1 is the distinguished vertex (a wedge or claw centre, a
# paw's degree-3 vertex) -- the same convention the colour tuple uses.
#' @keywords internal
#' @noRd
.motif_glyph_layout <- function(topology) {
    switch(topology,
        "edge" = list(
            xy = cbind(c(-0.5, 0.5), c(0, 0)),
            e = cbind(1, 2)
        ),
        "open" = list(
            xy = cbind(c(0, -0.8, 0.8), c(0.5, -0.4, -0.4)),
            e = rbind(c(1, 2), c(1, 3))
        ),
        "closed" = list(
            xy = cbind(c(0, -0.8, 0.8), c(0.7, -0.5, -0.5)),
            e = rbind(c(1, 2), c(1, 3), c(2, 3))
        ),
        "claw" = list(
            xy = cbind(c(0, 0, -0.85, 0.85), c(0.1, 1, -0.7, -0.7)),
            e = rbind(c(1, 2), c(1, 3), c(1, 4))
        ),
        "path" = list(
            xy = cbind(c(-1.1, -0.4, 0.4, 1.1), c(0, 0.45, -0.45, 0)),
            e = rbind(c(1, 2), c(2, 3), c(3, 4))
        ),
        "cycle" = list(
            xy = cbind(c(-0.7, 0.7, 0.7, -0.7), c(0.7, 0.7, -0.7, -0.7)),
            e = rbind(c(1, 2), c(2, 3), c(3, 4), c(1, 4))
        ),
        "paw" = list(
            # triangle on the left, pendant on the right, so the triangle is
            # visually obvious; slot 1 is the degree-3 vertex joining them
            xy = cbind(c(0.05, -0.75, -0.75, 1.05), c(0, 0.65, -0.65, 0)),
            e = rbind(c(1, 2), c(1, 3), c(2, 3), c(1, 4))
        ),
        "diamond" = list(
            xy = cbind(c(-0.6, 0.6, 0, 0), c(0, 0, 0.85, -0.85)),
            e = rbind(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4))
        ),
        "K4" = list(
            xy = cbind(c(-0.7, 0.7, 0.7, -0.7), c(0.7, 0.7, -0.7, -0.7)),
            e = rbind(
                c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4)
            )
        ),
        NULL
    )
}


# plotMotifEnrichment ####

#' @title plotMotifEnrichment
#' @name plotMotifEnrichment
#' @description
#' Triage plot for deciding which motifs are worth pursuing out of the
#' thousands a size-4 run produces. Enriched and depleted motifs appear in one
#' frame, split by topology.
#'
#' The default `style = "effect"` plots effect size against how often the motif
#' actually occurs, with significance as colour. That is deliberate: permutation
#' p-values are floored at `1 / (n_perm + 1)`, so on a real run most significant
#' classes sit exactly on that floor and a volcano's y-axis collapses to a
#' single line carrying no information. Effect size and occurrence count both
#' vary, and together they are what separates a real finding from a large
#' z-score built on three occurrences.
#'
#' `style = "volcano"` gives the conventional view for when p-values are not
#' saturated -- raise `n_perm` to lower the floor.
#' @inheritParams plot_output_params
#' @param gobject giotto object
#' @param motif_result output of [cellProximityMotifs()]
#' @param style `"effect"` plots effect size against occurrence count with
#'   significance as colour; `"volcano"` plots effect size against
#'   `-log10(p_adj)`. See Description for why `"effect"` is the default.
#' @param p_thresh adjusted p-value threshold for calling a motif significant
#' @param min_observed drop motif classes seen fewer than this many times.
#'   Rare classes carry unstable effect sizes and dominate the axes.
#' @param top number of most significant motifs to label
#' @param facet_topology draw one panel per topology
#' @returns ggplot
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' m <- cellProximityMotifs(g, cluster_column = "leiden_clus", size = 3L)
#' plotMotifEnrichment(g, m)
#' }
#' @export
plotMotifEnrichment <- function(gobject,
    motif_result,
    style = c("effect", "volcano"),
    p_thresh = 0.05,
    min_observed = 5,
    top = 15L,
    facet_topology = TRUE,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "plotMotifEnrichment") {
    .motif_require_contract(motif_result, "motif_result")
    style <- match.arg(style)
    d <- data.table::as.data.table(motif_result)
    d <- d[d$observed >= min_observed, ]
    if (!nrow(d)) {
        .gstop(sprintf(
            "no motif class was observed at least %s times; lower min_observed",
            min_observed
        ))
    }
    d$label <- .motif_label(d)
    d$neglogp <- -log10(pmax(d$p_adj, .Machine$double.xmin))
    d$direction <- data.table::fifelse(
        d$p_adj > p_thresh, "n.s.",
        data.table::fifelse(d$z > 0, "enriched", "depleted")
    )
    d$direction <- factor(d$direction, c("enriched", "depleted", "n.s."))

    lab <- d[d$p_adj <= p_thresh, ]
    if (nrow(lab)) {
        # setorder() takes column names, not expressions
        lab <- utils::head(lab[order(-abs(lab$z)), ], top)
    }

    pal <- c(enriched = "#B2182B", depleted = "#2166AC", n.s. = "grey75")

    if (style == "volcano") {
        pl <- ggplot2::ggplot(
            d, ggplot2::aes(x = .data$z, y = .data$neglogp)
        ) +
            ggplot2::geom_hline(
                yintercept = -log10(p_thresh), linetype = "dashed",
                linewidth = 0.3, colour = "grey60"
            ) +
            ggplot2::geom_point(
                ggplot2::aes(size = .data$observed, colour = .data$direction),
                alpha = 0.75
            ) +
            ggplot2::scale_size_continuous(
                transform = "log10", name = "occurrences", range = c(0.8, 5)
            ) +
            ggplot2::labs(
                x = "standardized effect (z)",
                y = expression(-log[10] ~ "adjusted p")
            )
    } else {
        pl <- ggplot2::ggplot(
            d, ggplot2::aes(x = .data$z, y = .data$observed)
        ) +
            ggplot2::geom_point(
                ggplot2::aes(colour = .data$direction), alpha = 0.8, size = 1.8
            ) +
            ggplot2::scale_y_continuous(
                transform = "log10", name = "occurrences"
            ) +
            ggplot2::labs(x = "standardized effect (z)")
    }

    pl <- pl +
        ggplot2::geom_vline(
            xintercept = 0, linewidth = 0.3, colour = "grey70"
        ) +
        ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
        .motif_theme()

    if (nrow(lab) && requireNamespace("ggrepel", quietly = TRUE)) {
        pl <- pl + ggrepel::geom_text_repel(
            data = lab,
            ggplot2::aes(label = .data$label),
            size = 2.4, max.overlaps = 20, min.segment.length = 0
        )
    } else if (nrow(lab)) {
        pl <- pl + ggplot2::geom_text(
            data = lab, ggplot2::aes(label = .data$label),
            size = 2.4, vjust = -0.8, check_overlap = TRUE
        )
    }
    if (isTRUE(facet_topology) && data.table::uniqueN(d$topology) > 1L) {
        pl <- pl + ggplot2::facet_wrap(~topology, scales = "free_x")
    }

    return(plot_output_handler(
        gobject = gobject, plot_object = pl,
        save_plot = save_plot, return_plot = return_plot,
        show_plot = show_plot, default_save_name = default_save_name,
        save_param = save_param
    ))
}


# plotMotifGlyphs ####

#' @title plotMotifGlyphs
#' @name plotMotifGlyphs
#' @description
#' Draws each motif as the little graph it actually is -- nodes coloured by
#' cell type, edges as drawn -- in a grid ranked by effect size.
#'
#' A `motif_id` like `size4_paw_A-B-B-C` is not readable at a glance; the same
#' motif as four coloured dots joined the way the tissue joins them is.
#' @inheritParams plot_output_params
#' @param gobject giotto object
#' @param motif_result output of [cellProximityMotifs()]
#' @param top number of motifs to draw
#' @param order_by `"z"`, `"p_adj"` or `"observed"`
#' @param direction which motifs to draw: `"enriched"`, `"depleted"` or
#'   `"both"`
#' @param min_observed drop motif classes seen fewer than this many times
#' @param ncol columns in the grid
#' @param cell_color_code optional named vector of cell type colours. Defaults
#'   to the object's discrete palette, so glyphs match the spatial plots.
#' @returns ggplot
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' m <- cellProximityMotifs(g, cluster_column = "leiden_clus", size = 3L)
#' plotMotifGlyphs(g, m, top = 12)
#' }
#' @export
plotMotifGlyphs <- function(gobject,
    motif_result,
    top = 12L,
    order_by = c("z", "p_adj", "observed"),
    direction = c("enriched", "depleted", "both"),
    min_observed = 5,
    ncol = 4L,
    cell_color_code = NULL,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "plotMotifGlyphs") {
    .motif_require_contract(motif_result, "motif_result")
    order_by <- match.arg(order_by)
    direction <- match.arg(direction)

    d <- data.table::as.data.table(motif_result)
    d <- d[d$observed >= min_observed, ]
    d <- switch(direction,
        enriched = d[d$z > 0, ],
        depleted = d[d$z < 0, ],
        both = d
    )
    if (!nrow(d)) .gstop("no motif class passed the filters")
    data.table::setorderv(
        d, order_by,
        order = if (order_by == "p_adj") 1L else -1L
    )
    if (direction == "both" && order_by == "z") {
        data.table::setorderv(d, "z", order = -1L)
    }
    d <- utils::head(d, top)
    d$label <- .motif_label(d)

    types <- sort(unique(unlist(d$color_tuple)))
    cols <- cell_color_code
    if (is.null(cols)) {
        cols <- set_default_color_discrete_cell(
            instrs = instructions(gobject)
        )(length(types))
        names(cols) <- types
    }

    nodes <- list()
    edges <- list()
    for (i in seq_len(nrow(d))) {
        lay <- .motif_glyph_layout(d$topology[i])
        if (is.null(lay)) next
        ct <- d$color_tuple[[i]]
        k <- length(ct)
        if (nrow(lay$xy) != k) next
        # One cell type per line. Pasted onto a single line these overflow the
        # facet strip and are silently clipped at both ends -- "Secretory
        # epithelial" three times is well past what a 2-inch panel can show.
        panel <- sprintf(
            "%s  (n=%s, z=%.1f)\n%s",
            d$topology[i], d$observed[i], d$z[i],
            paste(d$color_tuple[[i]], collapse = "\n")
        )
        nodes[[length(nodes) + 1L]] <- data.table::data.table(
            panel = panel, ord = i,
            x = lay$xy[, 1L], y = lay$xy[, 2L], cell_type = ct
        )
        edges[[length(edges) + 1L]] <- data.table::data.table(
            panel = panel, ord = i,
            x = lay$xy[lay$e[, 1L], 1L], y = lay$xy[lay$e[, 1L], 2L],
            xend = lay$xy[lay$e[, 2L], 1L], yend = lay$xy[lay$e[, 2L], 2L]
        )
    }
    if (!length(nodes)) .gstop("no motif had a known topology layout")
    nd <- data.table::rbindlist(nodes)
    ed <- data.table::rbindlist(edges)
    lev <- unique(nd$panel[order(nd$ord)])
    nd$panel <- factor(nd$panel, lev)
    ed$panel <- factor(ed$panel, lev)

    pl <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = ed,
            ggplot2::aes(
                x = .data$x, y = .data$y,
                xend = .data$xend, yend = .data$yend
            ),
            colour = "grey45", linewidth = 0.5
        ) +
        ggplot2::geom_point(
            data = nd,
            ggplot2::aes(x = .data$x, y = .data$y, fill = .data$cell_type),
            shape = 21, size = 5, stroke = 0.4, colour = "grey20"
        ) +
        ggplot2::scale_fill_manual(values = cols, name = "cell type") +
        ggplot2::facet_wrap(~panel, ncol = ncol) +
        ggplot2::coord_fixed(xlim = c(-1.35, 1.35), ylim = c(-1.2, 1.3)) +
        .motif_theme() +
        ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 6.5, lineheight = 0.95),
            panel.spacing = ggplot2::unit(0.35, "lines")
        )

    return(plot_output_handler(
        gobject = gobject, plot_object = pl,
        save_plot = save_plot, return_plot = return_plot,
        show_plot = show_plot, default_save_name = default_save_name,
        save_param = save_param
    ))
}


# spatMotifPlot ####

#' @title spatMotifPlot
#' @name spatMotifPlot
#' @description
#' Shows where a motif actually sits in the tissue: the cells taking part in
#' each occurrence, over a faint background of all other cells.
#'
#' Above `density_threshold` cells the occurrences are drawn as a
#' two-dimensional density instead of individual points, so the figure stays
#' readable on a million-cell section rather than becoming a solid block.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param motif_id one or more `motif_id` values to show
#' @param cluster_column cell metadata column with the cell type labels the
#'   motifs were computed from
#' @param spatial_network_name spatial network the motifs were computed on
#' @param size motif size the ids came from
#' @param point_size,background_size point sizes for participating and
#'   non-participating cells
#' @param density_threshold above this many participating cells, draw a density
#'   surface rather than points
#' @param max_per_class cap on occurrences fetched per motif class
#' @returns ggplot
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' m <- cellProximityMotifs(g, cluster_column = "leiden_clus", size = 3L)
#' spatMotifPlot(g, m$motif_id[1], cluster_column = "leiden_clus")
#' }
#' @export
spatMotifPlot <- function(gobject,
    motif_id,
    cluster_column,
    spat_unit = NULL,
    feat_type = NULL,
    spatial_network_name = "Delaunay_network",
    size = 3L,
    point_size = 1.2,
    background_size = 0.4,
    density_threshold = 20000L,
    max_per_class = 20000L,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "spatMotifPlot") {
    package_check("smotif", repository = "github:drieslab/smotif")
    spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    ig <- .motif_network_as_igraph(gobject, spat_unit, spatial_network_name)
    vids <- igraph::V(ig)$name
    meta <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type, output = "data.table", copy_obj = TRUE
    )
    lab <- .motif_align(meta, vids, cluster_column)
    el <- igraph::as_edgelist(ig, names = FALSE)

    inst <- smotif::motif_instances(
        smotif::SpatialGraph(
            nodes = data.table::data.table(
                cell_id = vids, x = 0, y = 0,
                cell_type = lab, sample_id = "s"
            ),
            edges = data.table::data.table(
                source = vids[el[, 1L]], target = vids[el[, 2L]],
                sample_id = "s"
            )
        ),
        motif_ids = motif_id, size = size, max_per_class = max_per_class
    )
    if (!nrow(inst)) {
        .gstop(sprintf(
            "no occurrences of %s were found",
            paste(motif_id, collapse = ", ")
        ))
    }

    locs <- getSpatialLocations(gobject,
        spat_unit = spat_unit,
        output = "data.table"
    )
    hit <- unique(inst$cell_id)
    bg <- locs[!locs$cell_ID %in% hit, ]
    fg <- merge(
        locs, unique(inst[, c("motif_id", "cell_id")]),
        by.x = "cell_ID", by.y = "cell_id"
    )

    pl <- ggplot2::ggplot() +
        ggplot2::geom_point(
            data = bg,
            ggplot2::aes(x = .data$sdimx, y = .data$sdimy),
            colour = "grey88", size = background_size
        )
    if (nrow(fg) > density_threshold) {
        pl <- pl + ggplot2::stat_density_2d(
            data = fg,
            ggplot2::aes(
                x = .data$sdimx, y = .data$sdimy,
                fill = ggplot2::after_stat(.data$level)
            ),
            geom = "polygon", alpha = 0.6, bins = 12
        ) +
            ggplot2::scale_fill_viridis_c(name = "occurrence\ndensity")
    } else {
        pl <- pl + ggplot2::geom_point(
            data = fg,
            ggplot2::aes(x = .data$sdimx, y = .data$sdimy,
                colour = .data$motif_id),
            size = point_size
        ) +
            ggplot2::scale_colour_manual(
                values = getColors("viridis", data.table::uniqueN(fg$motif_id)),
                name = NULL
            )
    }
    pl <- pl +
        ggplot2::coord_fixed() +
        ggplot2::labs(
            x = NULL, y = NULL,
            subtitle = sprintf(
                "%s occurrence(s), %s participating cells",
                data.table::uniqueN(inst$instance), length(hit)
            )
        ) +
        .motif_theme()

    return(plot_output_handler(
        gobject = gobject, plot_object = pl,
        save_plot = save_plot, return_plot = return_plot,
        show_plot = show_plot, default_save_name = default_save_name,
        save_param = save_param
    ))
}

#' @title findCellTypesFromEnrichment
#' @name findCellTypesFromEnrichment
#' @param gobject Giotto Object
#' @param spat_unit spatial unit in which the enrichment information is stored
#' @param feat_type feature type for which the enrichment information was
#' calculated
#' @param enrichment_name name of the spatial enrichment, as listed by
#' `list_spatial_enrichments_names()` (default "PAGE")
#' @param return_frequency_table see details. Default FALSE
#' @param min_score leave a cell unassigned (`NA`) when its top score is below
#' this (default NULL: no threshold). Not used for an enrichment of p-values.
#' @param max_padj leave a cell unassigned (`NA`) when the adjusted p-value of
#' its top type, from the `<enrichment_name>_padj` enrichment that
#' `p_value = TRUE` stores, is above this (default NULL: no threshold)
#' @param return_gobject add the assignment to the cell metadata as
#' `<enrichment_name>_type` and return the giotto object (default FALSE)
#' @returns a data.table of `cell_ID` and `probable_cell_type`, a frequency
#' table, or the giotto object (`return_gobject = TRUE`)
#' @details
#' Each cell is assigned the cell type with the highest score in the
#' enrichment or, for an enrichment of p-values, the lowest p-value. Ties go
#' to the first type. `min_score` and `max_padj` leave cells without a
#' confident type unassigned.
#'
#' By setting return_frequency_table to TRUE, a table of the number of cells
#' per cell type is returned instead, unassigned cells counted as `NA`.
#' @md
#' @export
findCellTypesFromEnrichment <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        enrichment_name = "PAGE",
        return_frequency_table = FALSE,
        min_score = NULL,
        max_padj = NULL,
        return_gobject = FALSE) {
    if (!inherits(gobject, "giotto")) {
        stop("gobject needs to be a giotto object")
    }
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    get_enr <- function(name) {
        getSpatialEnrichment(gobject,
            spat_unit = spat_unit, feat_type = feat_type, name = name,
            output = "spatEnrObj"
        )
    }

    enr <- get_enr(enrichment_name)
    is_p <- isTRUE(enr@misc$p_values_calculated)
    if (is_p && !is.null(min_score)) {
        stop("`min_score` applies to scores; \"", enrichment_name,
            "\" holds p-values.", call. = FALSE)
    }
    dt <- enr[]
    types <- setdiff(names(dt), "cell_ID")
    vals <- as.matrix(dt[, types, with = FALSE])
    top <- max.col(if (is_p) -vals else vals, ties.method = "first")
    type <- types[top]

    if (!is.null(min_score)) {
        type[vals[cbind(seq_along(top), top)] < min_score] <- NA_character_
    }
    if (!is.null(max_padj)) {
        padj_name <- paste0(enrichment_name, "_padj")
        padj <- get_enr(padj_name)[]
        padj <- as.matrix(padj[match(dt$cell_ID, padj$cell_ID),
            types, with = FALSE])
        type[padj[cbind(seq_along(top), top)] > max_padj] <- NA_character_
    }

    out <- data.table::data.table(cell_ID = dt$cell_ID,
                                  probable_cell_type = type)
    if (isTRUE(return_gobject)) {
        meta <- data.table::setnames(data.table::copy(out),
            "probable_cell_type", paste0(enrichment_name, "_type"))
        return(addCellMetadata(gobject,
            spat_unit = spat_unit, feat_type = feat_type,
            new_metadata = meta, by_column = TRUE,
            column_cell_ID = "cell_ID"
        ))
    }
    if (return_frequency_table) {
        return(table(out$probable_cell_type, useNA = "ifany"))
    }
    out
}

#' @title plotCellTypesFromEnrichment
#' @name plotCellTypesFromEnrichment
#' @param gobject Giotto Object
#' @param spat_unit spatial unit in which the enrichment information is stored
#' @param feat_type feature type for which the enrichment information was
#' calculated
#' @inheritParams findCellTypesFromEnrichment
#' @param title Title of the generated plot.
#'  Default `paste0(spat_unit,"cell types (maximum", enrichment_name, ")")`
#' @inheritParams plot_output_params
#' @returns plot
#' @details
#'
#' This function generates a bar plot of cell types vs the frequency
#' of that cell type in the data. These cell type results are
#' based on the provided `enrichment_name`, and will be determined
#' by the maximum value of the z-score or p-value for a given cell or
#' annotation.
#'
#' @export
plotCellTypesFromEnrichment <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        enrichment_name = "PAGE",
        min_score = NULL,
        max_padj = NULL,
        title = NULL,
        save_param = list(),
        default_save_name = "cell_types_from_enrichment",
        save_plot = NULL,
        show_plot = NULL,
        return_plot = NULL) {
    # guard clauses handled at first step downstream
    # therefore, omitting here.
    id_and_types <- findCellTypesFromEnrichment(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        enrichment_name = enrichment_name,
        min_score = min_score,
        max_padj = max_padj,
        return_frequency_table = FALSE
    )

    # data.table column
    probable_cell_type <- NULL

    if (is.null(title)) {
        title <- paste0(
            spat_unit, "cell types (maximum", enrichment_name, ")"
        )
    }

    pl <- ggplot2::ggplot(id_and_types, aes(x = probable_cell_type)) +
        ggplot2::geom_bar() +
        ggplot2::theme(
            axis.text.x = element_text(angle = 45),
            axis.ticks.length.x = unit(1.5, "cm")
        ) +
        ggplot2::labs(
            title = title,
            x = "Cell Type",
            y = "Frequency"
        )

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}

#' @title pieCellTypesFromEnrichment
#' @name pieCellTypesFromEnrichment
#' @param gobject Giotto Object
#' @param spat_unit spatial unit in which the enrichment information is stored
#' @param feat_type feature type for which the enrichment information was
#' calculated
#' @inheritParams findCellTypesFromEnrichment
#' @param title Title of the generated plot.
#'  Default `paste0(spat_unit,"cell types (maximum", enrichment_name, ")")`
#' @inheritParams plot_output_params
#' @returns pie plot
#' @details
#'
#' This function generates a pie chart of cell types by frequency.
#' These cell type resutls are based on the provided `enrichment_name`,
#' and will be determined by the maximum value of the z-score
#' or p-value for a given cell or annotation.
#' @export
pieCellTypesFromEnrichment <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        enrichment_name = "PAGE",
        min_score = NULL,
        max_padj = NULL,
        title = NULL,
        save_param = list(),
        default_save_name = "cell_types_from_enrichment_pie",
        save_plot = NULL,
        show_plot = NULL,
        return_plot = NULL) {
    # guard clauses handled one step downstream

    freq_table <- findCellTypesFromEnrichment(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        enrichment_name = enrichment_name,
        min_score = min_score,
        max_padj = max_padj,
        return_frequency_table = TRUE
    )

    freq_dt <- data.table::data.table(freq_table)

    colnames(freq_dt) <- c("cell_type", "num_cells")
    # data.table vars
    cell_type <- num_cells <- perc <- NULL

    cell_types <- unique(freq_dt$cell_type)
    total_cells <- sum(freq_dt$num_cells)

    for (i in cell_types) {
        # hackish, admittedly
        nullvar <- freq_dt[cell_type == i, perc := num_cells / sum(
            freq_dt$num_cells
        ) * 100]
    }
    rm(nullvar) # saves memory

    pl <- ggplot2::ggplot(
        as.data.frame(freq_dt),
        aes(
            x = "",
            y = perc,
            fill = cell_type
        )
    ) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = getDistinctColors(length(cell_types))) +
        theme_void() +
        labs(title = paste(
            spat_unit,
            " Cell Types (",
            as.character(total_cells),
            " Cells)"
        ))

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}

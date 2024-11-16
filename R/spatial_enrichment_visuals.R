#' @title findCellTypesFromEnrichment
#' @name findCellTypesFromEnrichment
#' @param gobject Giotto Object
#' @param spat_unit spatial unit in which the enrichment information is stored
#' @param feat_type feature type for which the enrichment information was
#' calculated
#' @param enrichment_name name of the spatial enrichment
#' i.e. output from GiottoClass::list_spatial_enrichment_names()
#' Default value is "PAGE_Z_score"
#' @param return_frequency_table see details. Default FALSE
#' @returns table
#' @details
#' This function returns a two-column matrix, one column
#' will contain cell IDs and the other will contain the
#' most likely cell type based on the provided enrichment.
#'
#' By setting return_frequency_table to TRUE, a table
#' will instead be returned. This table will detail
#' the number of occurences for each cell type, based on
#' the provided enrichment.
#'
#' The cell types are assigned by determining the maximum
#' value of the z-score or -log10(p-value) for a cell and
#' the associated cell types from the enrichment.
#'
#' @export
findCellTypesFromEnrichment <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        enrichment_name = "PAGE_z_score",
        return_frequency_table = FALSE) {
    # guard clauses

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

    # data.table variables
    probable_cell_type <- cell_ID <- NULL

    # extract p-value or z-socre from provided enrichment
    pz_enrich <- getSpatialEnrichment(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        name = enrichment_name
    )

    enrich_is_p_value <- pz_enrich@misc$p_values_calculated

    if (colnames(pz_enrich)[[1]] != "cell_ID") {
        selected_cols <- colnames(pz_enrich)[colnames(pz_enrich) != "cell_ID"]
        data.table::setcolorder(pz_enrich, c("cell_ID", selected_cols))
    }


    n_c <- ncol(pz_enrich)

    # Find the cell type column that corresponds to the
    # maximum value within a row and assign it into a
    # new column, mapping a cell to it's most likely type
    if (enrich_is_p_value) {
        pz_enrich[, probable_cell_type := names(
            .SD
        )[max.col(-.SD)], .SDcols = 2:n_c]
    } else {
        pz_enrich[, probable_cell_type := names(
            .SD
        )[max.col(.SD)], .SDcols = 2:n_c]
    }

    cell_ID_and_types_pz_enrich <- pz_enrich[, .(cell_ID, probable_cell_type)]

    if (return_frequency_table) {
        pz_enrich_cell_type_frequencies <- table(
            cell_ID_and_types_pz_enrich$probable_cell_type
        )
        return(pz_enrich_cell_type_frequencies)
    }

    return(cell_ID_and_types_pz_enrich)
}

#' @title plotCellTypesFromEnrichment
#' @name plotCellTypesFromEnrichment
#' @param gobject Giotto Object
#' @param spat_unit spatial unit in which the enrichment information is stored
#' @param feat_type feature type for which the enrichment information was
#' calculated
#' @param enrichment_name name of the spatial enrichment
#'  i.e. output from GiottoClass::list_spatial_enrichment_names()
#'  Default value is "PAGE_Z_score"
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
        enrichment_name = "PAGE_z_score",
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
#' @param enrichment_name name of the spatial enrichment
#'  i.e. output from GiottoClass::list_spatial_enrichment_names()
#'  Default value is "PAGE_Z_score"
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
        enrichment_name = "PAGE_z_score",
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

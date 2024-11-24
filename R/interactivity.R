# ** interactive plotting ####

#' Select image regions by plotting interactive polygons
#'
#' @description Plot interactive polygons on an image and retrieve the polygons
#' coordinates.
#' @param x A `ggplot` or `rast` plot object to draw polygons on
#' @param width,height An integer, defining the width/height in pixels.
#' @param ... Graphical parameters passed on to `polygon` or `geom_point`.
#'
#' @returns A `data.table` containing x,y coordinates from the plotted polygons.
#'
#' @export
plotInteractivePolygons <- function(
        x,
        width = "auto",
        height = "auto",
        ...) {
    package_check(pkg_name = "miniUI", repository = "CRAN")
    package_check(pkg_name = "shiny", repository = "CRAN")

    # data.table variables
    y <- name <- NULL

    if (is.null(x)) stop("plot object is empty")

    ## find min and max values for spatRaster image
    if ("SpatRaster" %in% class(x)) {
        ui <- miniUI::miniPage(
            miniUI::gadgetTitleBar("Plot Interactive Polygons"),
            miniUI::miniContentPanel(
                shiny::textInput(
                    "polygon_name",
                    label = "Polygon name",
                    value = "polygon 1"
                ),
                shiny::sliderInput("xrange",
                    label = "x coordinates",
                    min = min(terra::ext(x))[1],
                    max = max(terra::ext(x))[1],
                    value = c(
                        min(terra::ext(x))[1],
                        max(terra::ext(x))[1]
                    )
                ),
                shiny::sliderInput("yrange",
                    label = "y coordinates",
                    min = min(terra::ext(x))[2],
                    max = max(terra::ext(x))[2],
                    value = c(
                        min(terra::ext(x))[2],
                        max(terra::ext(x))[2]
                    )
                ),
                shiny::plotOutput("plot", click = "plot_click")
            )
        )
    } else { ## find min and max values for non-spatRaster image
        ui <- miniUI::miniPage(
            miniUI::gadgetTitleBar("Plot Interactive Polygons"),
            miniUI::miniContentPanel(
                shiny::textInput(
                    "polygon_name",
                    label = "Polygon name",
                    value = "polygon 1"
                ),
                shiny::sliderInput("xrange",
                    label = "x coordinates",
                    min = min(x[["layers"]][[1]]$data$sdimx),
                    max = max(x[["layers"]][[1]]$data$sdimx),
                    value = c(
                        min(x[["layers"]][[1]]$data$sdimx),
                        max(x[["layers"]][[1]]$data$sdimx)
                    )
                ),
                shiny::sliderInput("yrange",
                    label = "y coordinates",
                    min = min(x[["layers"]][[1]]$data$sdimy),
                    max = max(x[["layers"]][[1]]$data$sdimy),
                    value = c(
                        min(x[["layers"]][[1]]$data$sdimy),
                        max(x[["layers"]][[1]]$data$sdimy)
                    )
                ),
                shiny::plotOutput("plot", click = "plot_click")
            )
        )
    }

    server <- function(input, output, session) {
        output$plot <- shiny::renderPlot(
            {
                if ("ggplot" %in% class(x)) {
                    x$coordinates$default <- TRUE
                    x +
                        geom_polygon(
                            data = clicklist(),
                            aes(x, y, color = name),
                            alpha = 0,
                            show.legend = FALSE,
                            ...
                        ) +
                        coord_fixed(
                            xlim = c(input$xrange[1], input$xrange[2]),
                            ylim = c(input$yrange[1], input$yrange[2])
                        ) +
                        theme(legend.position = "none")
                } else {
                    terra::plot(x)
                    lapply(
                        split(clicklist(), by = "name"),
                        function(x) graphics::polygon(x$x, x$y, ...)
                    )
                }
            },
            res = 96,
            width = width,
            height = height
        )

        clicklist <- shiny::reactiveVal(data.table::data.table(
            x = numeric(), y = numeric(), name = character()
        )) # empty table
        shiny::observeEvent(input$plot_click, {
            click_x <- input$plot_click$x
            click_y <- input$plot_click$y
            polygon_name <- input$polygon_name
            temp <- clicklist() # get the table of past clicks
            temp <- rbind(temp, data.table::data.table(
                x = click_x, y = click_y, name = polygon_name
            ))
            clicklist(temp)
        })


        output$info <- shiny::renderTable(clicklist())

        shiny::observeEvent(input$done, {
            returnValue <- clicklist()
            shiny::stopApp(returnValue)
        })
    }
    shiny::runGadget(ui, server)
}


#' Get cells located within the polygons area
#'
#' @param gobject A Giotto object
#' @param polygon_name name of polygon selections
#' @param spat_unit spatial unit, default = 'cell'
#' @param spat_loc_name name of spatial locations to use, default = 'raw'
#' @param polygons character. A vector with polygon names to extract cells
#' from. If NULL, cells from all polygons are retrieved
#'
#' @returns A terra 'SpatVector' with cell ID, x y coordinates, and polygon ID
#' where each cell is located in.
#' @examples
#' ## Plot interactive polygons
#' g <- GiottoData::loadGiottoMini("visium")
#' my_polygon_coords <- data.frame(
#'     poly_ID = rep("polygon1", 3),
#'     sdimx = c(5477, 5959, 4720), sdimy = c(-4125, -2808, -5202)
#' )
#'
#' ## Add polygon coordinates to Giotto object
#' my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords,
#'     name = "selections"
#' )
#' g <- addGiottoPolygons(
#'     gobject = g,
#'     gpolygons = list(my_giotto_polygons)
#' )
#'
#' ## Add polygon IDs to cell metadata
#' addPolygonCells(g)
#'
#' ## Get only cells from polygon 1
#' getCellsFromPolygon(g)
#'
#' @export
getCellsFromPolygon <- function(
        gobject,
        polygon_name = "selections",
        spat_unit = "cell",
        spat_loc_name = "raw",
        polygons = NULL) {
    if (!inherits(gobject, "giotto")) {
        stop("gobject needs to be a giotto object")
    }

    ## get polygons spatial info
    polygon_spatVector <- getPolygonInfo(
        gobject = gobject,
        polygon_name = polygon_name,
        return_giottoPolygon = FALSE
    )

    ## get cell spatial locations
    spatial_locs <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table",
        copy_obj = TRUE
    )

    ## convert cell spatial locations to spatVector
    cells_spatVector <- terra::vect(
        as.matrix(spatial_locs[, c("sdimx", "sdimy")]),
        type = "points",
        atts = spatial_locs
    )

    polygonCells <- terra::intersect(cells_spatVector, polygon_spatVector)

    if (!is.null(polygons)) {
        polygonCells <- terra::subset(
            polygonCells, polygonCells$poly_ID %in% polygons
        )
    }

    return(polygonCells)
}


#' Add corresponding polygon IDs to cell metadata
#'
#' @param gobject A Giotto object
#' @param polygon_name name of polygon selections
#' @param feat_type feature name where metadata will be added
#' @param spat_unit spatial unit
#' @param spat_loc_name name of spatial locations to use
#' @param polygons polygon names to plot (e.g. 'polygon_1'). If NULL, plots
#' all available polygons
#' @param na.label polygon label for cells located outside of polygons area.
#' Default = "no_polygon"
#'
#' @returns A Giotto object with a modified cell_metadata slot that includes the
#' polygon name where each cell is located or no_polygon label if the cell is
#' not located within a polygon area
#' @examples
#' ## Plot interactive polygons
#' g <- GiottoData::loadGiottoMini("visium")
#' my_polygon_coords <- data.frame(
#'     poly_ID = rep("polygon1", 3),
#'     sdimx = c(5477, 5959, 4720),
#'     sdimy = c(-4125, -2808, -5202)
#' )
#'
#' ## Add polygon coordinates to Giotto object
#' my_giotto_polygons <- createGiottoPolygon(
#'     my_polygon_coords,
#'     name = "selections"
#' )
#'
#' g <- addGiottoPolygons(
#'     gobject = g,
#'     gpolygons = list(my_giotto_polygons)
#' )
#'
#' ## Add polygon IDs to cell metadata
#' g <- addPolygonCells(g)
#' pDataDT(g)
#' @export
addPolygonCells <- function(
        gobject,
        polygon_name = "selections",
        spat_unit = "cell",
        spat_loc_name = "raw",
        feat_type = "rna",
        polygons = NULL,
        na.label = "no_polygon") {
    ## verify gobject
    if (!inherits(gobject, "giotto")) {
        stop("gobject needs to be a giotto object")
    }

    ## get cells within each polygon
    polygon_cells <- getCellsFromPolygon(
        gobject = gobject,
        polygon_name = polygon_name,
        spat_unit = spat_unit,
        spat_loc_name = spat_loc_name,
        polygons = polygons
    )
    polygon_cells <- data.table::as.data.table(polygon_cells)
    data.table::setnames(polygon_cells, old = "poly_ID", new = polygon_name)

    ## get original cell metadata
    cell_metadata <- getCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )


    ## add polygon ID to cell metadata
    polygon_cells <- polygon_cells[, c("cell_ID", polygon_name), with = FALSE]
    new_cell_metadata <- data.table::merge.data.table(
        x = cell_metadata,
        y = polygon_cells,
        by = "cell_ID", all.x = TRUE
    )

    ## assign a default ID to cells outside of polygons
    selection_values <- new_cell_metadata[[polygon_name]]
    selection_values <- ifelse(
        is.na(selection_values), na.label, selection_values
    )
    new_cell_metadata[, c(polygon_name) := selection_values]

    ## keep original order of cells
    new_cell_metadata <- new_cell_metadata[match(
        cell_metadata$cell_ID,
        new_cell_metadata$cell_ID
    ), ]

    gobject <- addCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        new_metadata = new_cell_metadata[, -1]
    )

    return(gobject)
}


#' Compare gene expression between polygon areas
#'
#' @param gobject A Giotto object
#' @param polygon_name name of polygon selections
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param selected_feats vector of selected features to plot
#' @param expression_values gene expression values to use
#' ("normalized", "scaled", "custom")
#' @param method method to use to detect differentially expressed feats
#' ("scran", "gini", "mast")
#' @param \dots Arguments passed to \link[ComplexHeatmap]{Heatmap}
#'
#' @returns A ComplexHeatmap::Heatmap object
#' @examples
#' ## Plot interactive polygons
#' g <- GiottoData::loadGiottoMini("visium")
#' my_polygon_coords <- data.frame(
#'     poly_ID = rep("polygon1", 3),
#'     sdimx = c(5477, 5959, 4720), sdimy = c(-4125, -2808, -5202)
#' )
#'
#' ## Add polygon coordinates to Giotto object
#' my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords,
#'     name = "selections"
#' )
#' g <- addGiottoPolygons(
#'     gobject = g,
#'     gpolygons = list(my_giotto_polygons)
#' )
#'
#' ## Add polygon cells
#' g <- addPolygonCells(g)
#'
#' comparePolygonExpression(g)
#' @export
comparePolygonExpression <- function(
        gobject,
        polygon_name = "selections",
        spat_unit = "cell",
        feat_type = "rna",
        selected_feats = "top_genes",
        expression_values = "normalized",
        method = "scran",
        ...) {
    # verify gobject
    if (!inherits(gobject, "giotto")) {
        stop("gobject needs to be a giotto object")
    }

    # get expression
    my_expression <- getExpression(gobject,
        values = expression_values,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "matrix"
    )

    # get cell_ID and poly_ID from metadata
    my_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )

    my_metadata <- my_metadata[, c("cell_ID", polygon_name), with = FALSE]

    if (length(selected_feats) == 1 && selected_feats == "top_genes") {
        # find top features
        scran_results <- findMarkers_one_vs_all(gobject,
            spat_unit = "cell",
            feat_type = "rna",
            method = method,
            expression_values = "normalized",
            cluster_column = polygon_name,
            min_feats = 10
        )

        selected_feats <- scran_results[, head(.SD, 2), by = "cluster"]$feats
    }
    # select features
    my_expression <- my_expression[selected_feats, ]

    # convert to data frame
    my_rownames <- rownames(my_expression)

    # calculate zscore

    my_zscores <- my_expression

    for (gene in my_rownames) {
        mean_expression_gene <- mean(my_expression[gene, ])
        sd_expression_gene <- stats::sd(my_expression[gene, ])
        for (cell in colnames(my_expression)) {
            my_zscores[gene, cell] <- (
                my_expression[gene, cell] - mean_expression_gene) /
                sd_expression_gene
        }
    }

    # calculate mean zscore per polygon
    my_zscores_mean <- data.table::data.table(feat_ID = my_rownames)

    for (i in unique(my_metadata[[polygon_name]])) {
        my_cells <- my_metadata[my_metadata[[polygon_name]] == i, "cell_ID"]
        my_sub_zscores <- my_zscores[, my_cells$cell_ID]
        mean_zscores <- Matrix::rowMeans(my_sub_zscores)
        my_zscores_mean <- cbind(my_zscores_mean, mean_zscores)
    }

    my_zscores_mean <- as.matrix(my_zscores_mean[, -1])
    colnames(my_zscores_mean) <- unique(my_metadata[[polygon_name]])
    rownames(my_zscores_mean) <- my_rownames

    # plot heatmap
    my_heatmap <- ComplexHeatmap::Heatmap(my_zscores_mean,
        heatmap_legend_param = list(title = "Normalized mean z score"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_order = mixedsort(colnames(my_zscores_mean)),
        ...
    )
    return(my_heatmap)
}

#' Compare cell types percent per polygon
#'
#' @param gobject A Giotto object
#' @param polygon_name name of polygon selections
#' @param spat_unit spatial unit. Default = "cell"
#' @param feat_type feature type. Default =  "rna"
#' @param cell_type_column column name within the cell metadata table to use
#' @param ... Additional parameters passed to ComplexHeatmap::Heatmap
#'
#' @returns A ComplexHeatmap::Heatmap
#' @examples
#' ## Plot interactive polygons
#' g <- GiottoData::loadGiottoMini("visium")
#' my_polygon_coords <- data.frame(
#'     poly_ID = rep("polygon1", 3),
#'     sdimx = c(5477, 5959, 4720), sdimy = c(-4125, -2808, -5202)
#' )
#'
#' ## Add polygon coordinates to Giotto object
#' my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords,
#'     name = "selections"
#' )
#' g <- addGiottoPolygons(
#'     gobject = g,
#'     gpolygons = list(my_giotto_polygons)
#' )
#'
#' ## Add polygon cells
#' g <- addPolygonCells(g)
#'
#' compareCellAbundance(g)
#' @export
compareCellAbundance <- function(
        gobject,
        polygon_name = "selections",
        spat_unit = "cell",
        feat_type = "rna",
        cell_type_column = "leiden_clus",
        ...) {
    # verify gobject
    if (!inherits(gobject, "giotto")) {
        stop("gobject needs to be a giotto object")
    }

    # get poly_ID and cell_type from metadata
    my_metadata <- getCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    columns_to_select <- c(polygon_name, cell_type_column)
    my_metadata <- my_metadata[, columns_to_select, with = FALSE]

    # count cell_type per polygon
    my_cell_counts <- table(my_metadata)

    my_cell_percent <- 100 * my_cell_counts / rowSums(my_cell_counts)

    # convert to matrix
    my_matrix <- Matrix::as.matrix(my_cell_percent)

    rownames(my_matrix) <- rownames(my_cell_percent)
    colnames(my_matrix) <- colnames(my_cell_percent)

    # plot heatmap
    my_heatmap <- ComplexHeatmap::Heatmap(t_flex(my_matrix),
        heatmap_legend_param = list(title = "Cell type percent\nper polygon"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        ...
    )
    return(my_heatmap)
}


#' Plot stored polygons
#'
#' @param gobject A Giotto object with polygon coordinates
#' @param polygon_name name of polygon selections
#' @param x A ggplot2, spatPlot or terra::rast object
#' @param spat_unit spatial unit
#' @param polygons character. Vector of polygon names to plot. If NULL, all
#' polygons are plotted
#' @param ... Additional parameters passed to ggplot2::geom_polygon() or
#' graphics::polygon
#'
#' @returns A ggplot2 image
#' @examples
#' ## Plot interactive polygons
#' g <- GiottoData::loadGiottoMini("visium")
#' my_polygon_coords <- data.frame(
#'     poly_ID = rep("polygon1", 3),
#'     sdimx = c(5477, 5959, 4720), sdimy = c(-4125, -2808, -5202)
#' )
#'
#' ## Add polygon coordinates to Giotto object
#' my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords,
#'     name = "selections"
#' )
#' g <- addGiottoPolygons(
#'     gobject = g,
#'     gpolygons = list(my_giotto_polygons)
#' )
#'
#' ## Add polygon cells
#' g <- addPolygonCells(g)
#'
#' ## Create spatplot
#' x <- spatPlot2D(g, return_plot = TRUE)
#'
#' plotPolygons(g, x = x)
#' @export
plotPolygons <- function(
        gobject,
        polygon_name = "selections",
        x,
        spat_unit = "cell",
        polygons = NULL,
        ...) {
    ## verify gobject
    if (!inherits(gobject, "giotto")) {
        stop("gobject must be a Giotto object")
    }

    y <- geom <- NULL

    ## verify plot exists
    if (is.null(x)) stop("A plot object must be provided")

    ## get polygons spatial info
    polygon_spatVector <- getPolygonInfo(
        gobject = gobject,
        polygon_name = polygon_name,
        return_giottoPolygon = FALSE
    )

    coordinates <- terra::geom(polygon_spatVector, df = TRUE)

    if (!is.null(polygons)) {
        ## replace polygon names
        for (i in seq_along(unlist(polygon_spatVector[["poly_ID"]]))) {
            coordinates$geom <- replace(
                coordinates$geom,
                coordinates$geom == i,
                unlist(polygon_spatVector[["poly_ID"]])[i]
            )
        }

        coordinates <- coordinates[coordinates$geom %in% polygons, ]
    }

    ## plot over ggplot or spatPlot
    if (inherits(x, "ggplot")) {
        x +
            geom_polygon(
                data = coordinates,
                aes(x, y, colour = factor(geom), group = geom),
                alpha = 0, show.legend = FALSE, ...
            ) +
            theme(legend.position = "none")
    } else {
        terra::plot(x)
        lapply(
            split(coordinates, by = "name"),
            function(x) graphics::polygon(x$x, x$y, ...)
        )
    }
}

# ** 3D interactive plotting ####

#' Plot interactive 3D spatial plot
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param cell_color character. What to color cells by
#' (e.g. metadata col or spatial enrichment col)
#' @param cell_color_code character. discrete colors to use. Palette to use or
#' named vector of colors
#' @param point_size size of point (cell)
#' @param width plot width
#' @param height plot height
#'
#' @returns data.table with selected cell_IDs, spatial coordinates, and
#' cluster_ID.
#' @export
plotInteractive3D <- function(
        gobject, spat_unit = "cell", feat_type = "rna",
        cell_color = "leiden_clus",
        cell_color_code = NULL, point_size = 0.5,
        width = "100%", height = "400px") {
    package_check(
        c("plotly", "miniUI", "shiny"),
        repository = c("CRAN:plotly", "CRAN:miniUI", "CRAN:shiny")
    )

    # NSE vars
    sdimx <- sdimy <- sdimz <- cell_ID <- NULL

    cell_metadata_table <- pDataDT(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    spatial_coordinates <- getSpatialLocations(gobject,
        spat_unit = spat_unit,
        output = "data.table"
    )

    data <- merge(cell_metadata_table, spatial_coordinates)
    extent <- c(
        xmin = min(data$sdimx), xmax = max(data$sdimx),
        ymin = min(data$sdimy), ymax = max(data$sdimy),
        zmin = min(data$sdimz), zmax = max(data$sdimz)
    )
    sorted_colors <- mixedsort(unique(data[[cell_color]]))

    ui <- miniUI::miniPage(
        miniUI::gadgetTitleBar("Slide to select axis range"),
        miniUI::miniContentPanel(
            shiny::fluidRow(
                shiny::column(4,
                    offset = 1,
                    # Move the slide bar to select z-axis ranges
                    shiny::sliderInput("xrange",
                        label = "x-axis",
                        min = extent[["xmin"]],
                        max = extent[["xmax"]],
                        value = c(extent[["xmin"]], extent[["xmax"]])
                    ),
                    shiny::sliderInput("yrange",
                        label = "y-axis",
                        min = extent[["ymin"]],
                        max = extent[["ymax"]],
                        value = c(extent[["ymin"]], extent[["ymax"]])
                    ),
                    shiny::sliderInput("zrange",
                        label = "z-axis",
                        min = extent[["zmin"]],
                        max = extent[["zmax"]],
                        value = c(extent[["zmin"]], extent[["zmax"]])
                    )
                ),
                shiny::column(4,
                    offset = 2,
                    shiny::checkboxGroupInput("clusters",
                        label = "clusters",
                        choices = sorted_colors,
                        selected = sorted_colors
                    ),
                )
            ),

            # Plot
            plotly::plotlyOutput("plot1", width = width, height = height)
        ),
    )

    server <- function(input, output, session) {
        selected_data <- shiny::reactive({
            data[data[[cell_color]] %in% input$clusters, ] %>%
                plotly::filter(
                    sdimx >= input$xrange[1] & sdimx <= input$xrange[2] &
                        sdimy >= input$yrange[1] & sdimy <= input$yrange[2] &
                        sdimz >= input$zrange[1] & sdimz <= input$zrange[2]
                ) %>%
                plotly::select(cell_ID, sdimx, sdimy, sdimz, cell_color)
        })

        # Render the plot
        output$plot1 <- plotly::renderPlotly({
            selected_data <- selected_data()

            # Plot the data
            plotly::plot_ly(selected_data,
                x = ~sdimx,
                y = ~sdimy,
                z = ~sdimz,
                color = as.factor(selected_data[[5]]),
                colors = cell_color_code,
                size = point_size
            )
        })

        # Handle the Done button being pressed.

        output$info <- shiny::renderTable(selected_data())

        shiny::observeEvent(input$done, {
            returnValue <- selected_data()
            shiny::stopApp(returnValue)
        })
    }

    shiny::runGadget(ui, server)
}

#' Create a local anndata zarr folder
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression expression values to extract (e.g. "raw", "normalized",
#' "scaled")
#' @param output_path path to create and save the anndata zarr folder
#'
#' @return local anndata zarr folder
#' @export
#'
#' @examples
#' # using the mini visium object
#' giotto_object <- GiottoData::loadGiottoMini("visium")
#'
#' giottoToAnndataZarr(giotto_object,
#'     expression = "raw",
#'     output_path = tempdir()
#' )
#'
#' # using the mini vizgen object
#' giotto_object <- GiottoData::loadGiottoMini("vizgen")
#'
#' giottoToAnndataZarr(giotto_object,
#'     spat_unit = "aggregate",
#'     expression = "scaled",
#'     output_path = tempdir()
#' )
giottoToAnndataZarr <- function(
        gobject, spat_unit = NULL,
        feat_type = NULL, expression = "raw",
        output_path) {
    proc <- basilisk::basiliskStart(GiottoClass::instructions(
        gobject = gobject,
        param = "python_path"
    ))
    on.exit(basilisk::basiliskStop(proc))

    success <- basilisk::basiliskRun(
        proc,
        function(
        gobject,
        output_path,
        expression) {
            anndata <- reticulate::import("anndata")
            zarr <- reticulate::import("zarr")

            # extract expression matrix
            X <- t(as.matrix(GiottoClass::getExpression(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                values = expression,
                output = "matrix"
            )))

            # extract cell metadata
            obs <- as.data.frame(GiottoClass::getCellMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "data.table"
            ))

            rownames(obs) <- obs$cell_ID

            # extract feature metadata
            var <- as.data.frame(GiottoClass::getFeatureMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                output = "data.table"
            ))

            obsm <- list()

            # extract spatial locations
            spatial_locs <- as.data.frame(GiottoClass::getSpatialLocations(
                gobject = gobject,
                spat_unit = spat_unit,
                output = "data.table"
            ))

            if (!is.null(spatial_locs)) {
                rownames(spatial_locs) <- spatial_locs$cell_ID
                spatial_locs <- spatial_locs[obs$cell_ID, ]
                spatial_locs_matrix <- as.matrix(spatial_locs[, 1:2])

                obsm[["spatial"]] <- spatial_locs_matrix
            }

            # extract pca
            dim_reducs_pca <- GiottoClass::getDimReduction(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                reduction_method = "pca",
                output = "matrix"
            )

            if (!is.null(dim_reducs_pca)) {
                obsm[["pca"]] <- dim_reducs_pca[obs$cell_ID, ]
            }

            # extract umap
            dim_reducs_umap <- GiottoClass::getDimReduction(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                reduction_method = "umap",
                name = "umap",
                output = "matrix"
            )

            if (!is.null(dim_reducs_umap)) {
                obsm[["umap"]] <- dim_reducs_umap[obs$cell_ID, ]
            }

            # extract tSNE
            dim_reducs_tsne <- GiottoClass::getDimReduction(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                reduction_method = "tsne",
                name = "tsne",
                output = "matrix"
            )

            if (!is.null(dim_reducs_tsne)) {
                obsm[["tsne"]] <- dim_reducs_tsne[obs$cell_ID, ]
            }

            adata <- anndata$AnnData(X = X, obs = obs, var = var, obsm = obsm)

            adata$write_zarr(output_path)
            return(TRUE)
        },
        gobject = gobject, output_path = output_path, expression = expression
    )

    return(success)
}

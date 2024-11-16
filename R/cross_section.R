# TODO Create S4 crossSectionObj Class


# cross section helper functions ####

#' @title create_crossSection_object
#' @name create_crossSection_object
#' @description create a crossSection object
#' @param name name of cross section object. (default = cross_section)
#' @param method method to define the cross section plane.
#' @param thickness_unit unit of the virtual section thickness. If "cell",
#' average size of the observed cells is used as length unit. If "natural",
#' the unit of cell location coordinates is used.(default = cell)
#' @param slice_thickness thickness of slice
#' @param cell_distance_estimate_method method to estimate average distance
#' between neighboring cells. (default = mean)
#' @param extend_ratio deciding the span of the cross section meshgrid, as a
#' ratio of extension compared to the borders of the virtual tissue section.
#' (default = 0.2)
#' @param plane_equation a numerical vector of length 4, in the form of
#' c(A,B,C,D), which defines plane Ax+By+Cz=D.
#' @param mesh_grid_n number of meshgrid lines to generate along both
#' directions for the cross section plane.
#' @param mesh_obj object that stores the cross section meshgrid information.
#' @param cell_subset cells selected by the cross section
#' @param cell_subset_spatial_locations locations of cells selected by the
#' cross section
#' @param cell_subset_projection_locations 3D projection coordinates of
#' selected cells onto the cross section plane
#' @param cell_subset_projection_PCA pca of projection coordinates
#' @param cell_subset_projection_coords 2D PCA coordinates of selected cells
#' in the cross section plane
#' @returns crossSection object
create_crossSection_object <- function(
        name = NULL,
        method = NULL,
        thickness_unit = NULL,
        slice_thickness = NULL,
        cell_distance_estimate_method = NULL,
        extend_ratio = NULL,
        plane_equation = NULL,
        mesh_grid_n = NULL,
        mesh_obj = NULL,
        cell_subset = NULL,
        cell_subset_spatial_locations = NULL,
        cell_subset_projection_locations = NULL,
        cell_subset_projection_PCA = NULL,
        cell_subset_projection_coords = NULL) {
    crossSection_obj <- list(
        "method" = method,
        "thickness_unit" = thickness_unit,
        "slice_thickness" = slice_thickness,
        "plane_equation" = plane_equation,
        "mesh_grid_n" = mesh_grid_n,
        "mesh_obj" = mesh_obj,
        "cell_subset" = cell_subset,
        "cell_subset_spatial_locations" = cell_subset_spatial_locations,
        "cell_subset_projection_locations" = cell_subset_projection_locations,
        "cell_subset_projection_PCA" = cell_subset_projection_PCA,
        "cell_subset_projection_coords" = cell_subset_projection_coords
    )
}

#' @title read_crossSection
#' @name read_crossSection
#' @description read a cross section object from a giotto object
#' @param gobject gobject
#' @param spat_unit spatial unit
#' @param name name
#' @param spatial_network_name spatial_network_name
#' @returns crossSectionObjects
#' @keywords internal
read_crossSection <- function(
        gobject,
        spat_unit = NULL,
        name = NULL,
        spatial_network_name = NULL) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )

    if (is.null(spatial_network_name)) {
        stop("spatial_network_name is not specified.")
    }

    sn <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        set_defaults = FALSE,
        copy_obj = FALSE,
        verbose = FALSE,
        output = "spatialNetworkObj"
    )

    cs_list <- slot(sn, "crossSectionObjects")

    if (length(cs_list) == 0L) {
        stop("No cross section object has been created.")
    }

    if (is.null(name)) {
        name <- names(cs_list)[length(cs_list)]

        default_name_msg <- sprintf(
            "cross section object is not specified, \n%s \n'%s'",
            "reading the last one from the existing list:", name
        )
    }

    if (!name %in% names(cs_list)) {
        stop(sprintf(
            "crossSectionObject '%s' has not been created.",
            name
        ))
    }

    crossSection_obj <- cs_list[[name]]

    return(crossSection_obj)
}



#' @title estimateCellCellDistance
#' @name estimateCellCellDistance
#' @description estimate average distance between neighboring cells
#' @param gobject gobject
#' @param spat_unit spatial unit
#' @param spatial_network_name spatial_network_name
#' @param method method
#' @returns matrix
#' @keywords internal
estimateCellCellDistance <- function(
        gobject,
        spat_unit = NULL,
        spatial_network_name = "Delaunay_network",
        method = c("mean", "median")) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )

    net <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "networkDT"
    )

    CellCellDistance <- get_distance(
        networkDT = net,
        method = method
    )

    return(CellCellDistance)
}
#' @title get_sectionThickness
#' @name get_sectionThickness
#' @description get section thickness
#' @param gobject gobject
#' @param spat_unit spatial unit
#' @param thickness_unit thickness_unit
#' @param spatial_network_name spatial_network_name
#' @param cell_distance_estimate_method cell_distance_estimate_method
#' @param plane_equation plane_equation
#' @returns numeric
#' @keywords internal
get_sectionThickness <- function(
        gobject,
        spat_unit = NULL,
        thickness_unit = c("cell", "natural"),
        slice_thickness = 2,
        spatial_network_name = "Delaunay_network",
        cell_distance_estimate_method = c("mean", "median"),
        plane_equation = NULL) {
    thickness_unit <- match.arg(thickness_unit, c("cell", "natural"))

    section_thickness <- switch(thickness_unit,
        "cell" = {
            CellCellDistance <- estimateCellCellDistance(
                gobject = gobject,
                spat_unit = spat_unit,
                method = cell_distance_estimate_method,
                spatial_network_name = spatial_network_name
            )
            CellCellDistance * slice_thickness
        },
        "natural" = slice_thickness
    )

    return(section_thickness)
}

#' @title projection_fun
#' @name projection_fun
#' @description project a point onto a plane
#' @param point_to_project point_to_project
#' @param plane_point plane_point
#' @param plane_norm plane_norm
#' @returns numeric
#' @keywords internal
projection_fun <- function(point_to_project, plane_point, plane_norm) {
    a <- plane_norm[1]
    b <- plane_norm[2]
    c <- plane_norm[3]
    x <- point_to_project[1]
    y <- point_to_project[2]
    z <- point_to_project[3]
    d <- plane_point[1]
    e <- plane_point[2]
    f <- plane_point[3]
    t <- (a * d - a * x + b * e - b * y + c * f - c * z) / (a^2 + b^2 + c^2)
    xp <- x + t * a
    yp <- y + t * b
    zp <- z + t * c
    projection <- c(xp, yp, zp)
    return(projection)
}

#' @title adapt_aspect_ratio
#' @name adapt_aspect_ratio
#' @description adapt the aspact ratio after inserting cross section mesh grid
#' lines
#' @param current_ratio current_ratio
#' @param cell_locations cell_locations
#' @param sdimx sdimx
#' @param sdimy sdimy
#' @param sdimz sdimz
#' @param mesh_obj mesh_obj
#' @returns numeric
#' @keywords internal
adapt_aspect_ratio <- function(
        current_ratio, cell_locations,
        sdimx = NULL, sdimy = NULL, sdimz = NULL,
        mesh_obj = NULL) {
    x_range <- max(cell_locations[[sdimx]]) - min(cell_locations[[sdimx]])
    y_range <- max(cell_locations[[sdimy]]) - min(cell_locations[[sdimy]])
    z_range <- max(cell_locations[[sdimz]]) - min(cell_locations[[sdimz]])

    x_mesh_range <- max(mesh_obj$mesh_grid_lines$mesh_grid_lines_X) - min(
        mesh_obj$mesh_grid_lines$mesh_grid_lines_X
    )
    y_mesh_range <- max(mesh_obj$mesh_grid_lines$mesh_grid_lines_Y) - min(
        mesh_obj$mesh_grid_lines$mesh_grid_lines_Y
    )
    z_mesh_range <- max(mesh_obj$mesh_grid_lines$mesh_grid_lines_Z) - min(
        mesh_obj$mesh_grid_lines$mesh_grid_lines_Z
    )

    if (x_mesh_range > x_range) {
        x_adapt <- x_mesh_range / x_range
    } else {
        x_adapt <- 1
    }
    if (y_mesh_range > y_range) {
        y_adapt <- y_mesh_range / y_range
    } else {
        y_adapt <- 1
    }
    if (z_mesh_range > z_range) {
        z_adapt <- z_mesh_range / z_range
    } else {
        z_adapt <- 1
    }

    new_ratio <- as.numeric(current_ratio) * c(
        as.numeric(x_adapt), as.numeric(y_adapt), as.numeric(z_adapt)
    )
    new_ratio <- new_ratio / min(new_ratio)
    return(new_ratio)
}

# mesh grid line helper functions ####

#' @title extend_vector
#' @name extend_vector
#' @description extend the range of a vector by a given ratio
#' @param x x
#' @param extend_ratio extend_ratio
#' @returns numeric
#' @keywords internal
extend_vector <- function(x, extend_ratio) {
    x_center <- (max(x) + min(x)) / 2
    y <- (x - x_center) * (extend_ratio + 1) + x_center

    return(y)
}

#' @title find_x_y_ranges
#' @name find_x_y_ranges
#' @description get the extended ranges of x and y
#' @param data data
#' @param extend_ratio extend_ratio
#' @returns list
#' @keywords internal
find_x_y_ranges <- function(data, extend_ratio) {
    x_extend <- extend_vector(data[, 1], extend_ratio)
    y_extend <- extend_vector(data[, 2], extend_ratio)

    x_min <- min(x_extend)
    x_max <- max(x_extend)
    y_min <- min(y_extend)
    y_max <- max(y_extend)

    out <- list(
        "x_min" = x_min,
        "x_max" = x_max,
        "y_min" = y_min,
        "y_max" = y_max
    )
}

#' @title create_2d_mesh_grid_line_obj
#' @name create_2d_mesh_grid_line_obj
#' @description create 2d mesh grid line object
#' @param x_min x_min
#' @param x_max x_max
#' @param y_min y_min
#' @param y_max y_max
#' @param mesh_grid_n mesh_grid_n
#' @returns 2d mesh grid line object
#' @keywords internal
create_2d_mesh_grid_line_obj <- function(x_min, x_max, y_min,
    y_max, mesh_grid_n) {
    x_grid <- seq(x_min, x_max, length.out = mesh_grid_n)
    y_grid <- seq(y_min, y_max, length.out = mesh_grid_n)

    mesh_grid_lines_X <- cbind(
        matrix(rep(x_grid, mesh_grid_n), nrow = mesh_grid_n, byrow = TRUE),
        matrix(rep(x_grid, mesh_grid_n), nrow = mesh_grid_n, byrow = FALSE)
    )

    mesh_grid_lines_Y <- cbind(
        matrix(rep(y_grid, mesh_grid_n), nrow = mesh_grid_n, byrow = FALSE),
        matrix(rep(y_grid, mesh_grid_n), nrow = mesh_grid_n, byrow = TRUE)
    )


    mesh_grid_line_obj_2d <- list(
        "mesh_grid_lines_X" = mesh_grid_lines_X,
        "mesh_grid_lines_Y" = mesh_grid_lines_Y
    )
    return(mesh_grid_line_obj_2d)
}

#' @title reshape_to_data_point
#' @name reshape_to_data_point
#' @description reshape a mesh grid line object to data point matrix
#' @param mesh_grid_obj mesh_grid_obj
#' @returns matrix
#' @keywords internal
reshape_to_data_point <- function(mesh_grid_obj) {
    if (length(mesh_grid_obj) == 3) {
        data_points <- cbind(
            as.vector(mesh_grid_obj[[1]]),
            as.vector(mesh_grid_obj[[2]]),
            as.vector(mesh_grid_obj[[3]])
        )
    } else if (length(mesh_grid_obj) == 2) {
        data_points <- cbind(
            as.vector(mesh_grid_obj[[1]]),
            as.vector(mesh_grid_obj[[2]])
        )
    }
    return(data_points)
}

#' @title reshape_to_mesh_grid_obj
#' @name reshape_to_mesh_grid_obj
#' @description reshape a data point matrix to a mesh grid line object
#' @param data_points data_points
#' @param mesh_grid_n mesh_grid_n
#' @returns list
#' @keywords internal
reshape_to_mesh_grid_obj <- function(data_points, mesh_grid_n) {
    if (dim(data_points)[2] == 2) {
        mesh_grid_lines_X <- matrix(
            data_points[, 1],
            nrow = mesh_grid_n, byrow = FALSE
        )
        mesh_grid_lines_Y <- matrix(
            data_points[, 2],
            nrow = mesh_grid_n, byrow = FALSE
        )

        mesh_grid_obj <- list(
            "mesh_grid_lines_X" = mesh_grid_lines_X,
            "mesh_grid_lines_Y" = mesh_grid_lines_Y
        )
    } else if (dim(data_points)[2] == 3) {
        mesh_grid_lines_X <- matrix(
            data_points[, 1],
            nrow = mesh_grid_n, byrow = FALSE
        )
        mesh_grid_lines_Y <- matrix(
            data_points[, 2],
            nrow = mesh_grid_n, byrow = FALSE
        )
        mesh_grid_lines_Z <- matrix(
            data_points[, 3],
            nrow = mesh_grid_n, byrow = FALSE
        )
        mesh_grid_obj <- list(
            "mesh_grid_lines_X" = mesh_grid_lines_X,
            "mesh_grid_lines_Y" = mesh_grid_lines_Y,
            "mesh_grid_lines_Z" = mesh_grid_lines_Z
        )
    }
    return(mesh_grid_obj)
}


#' @title transform_2d_mesh_to_3d_mesh
#' @name transform_2d_mesh_to_3d_mesh
#' @description transform 2d mesh to 3d mesh by reversing PCA
#' @param mesh_line_obj_2d mesh_line_obj_2d
#' @param pca_out pca_out
#' @param center_vec center_vec
#' @param mesh_grid_n mesh_grid_n
#' @returns 3d mesh
#' @keywords internal
transform_2d_mesh_to_3d_mesh <- function(mesh_line_obj_2d, pca_out,
    center_vec, mesh_grid_n) {
    data_point_2d <- reshape_to_data_point(mesh_line_obj_2d)
    center_mat <- matrix(
        rep(center_vec, dim(data_point_2d)[1]),
        nrow = dim(data_point_2d)[1], byrow = TRUE
    )
    data_point_3d <- cbind(
        data_point_2d,
        rep(0, dim(data_point_2d)[1])
    ) %*% t((pca_out$rotation)) + center_mat
    mesh_grid_line_obj_3d <- reshape_to_mesh_grid_obj(
        data_point_3d, mesh_grid_n
    )

    return(mesh_grid_line_obj_3d)
}

#' @title get_cross_section_coordinates
#' @name get_cross_section_coordinates
#' @description get local coordinates within cross section plane
#' @param cell_subset_projection_locations cell_subset_projection_locations
#' @returns data.table
#' @keywords internal
get_cross_section_coordinates <- function(cell_subset_projection_locations) {
    cell_subset_projection_PCA <- stats::prcomp(
        cell_subset_projection_locations
    )

    cell_subset_projection_coords <- cell_subset_projection_PCA$x[
        , c("PC1", "PC2")
    ]

    return(cell_subset_projection_coords)
}

#' @title create_mesh_grid_lines
#' @name create_mesh_grid_lines
#' @description create mesh grid lines for cross section
#' @param cell_subset_projection_locations cell_subset_projection_locations
#' @param extend_ratio extend_ratio
#' @param mesh_grid_n mesh_grid_n
#' @returns mesh grid lines
#' @keywords internal
create_mesh_grid_lines <- function(cell_subset_projection_locations,
    extend_ratio, mesh_grid_n) {
    cell_subset_projection_PCA <- stats::prcomp(
        cell_subset_projection_locations
    )

    cell_subset_projection_coords <- cell_subset_projection_PCA$x[
        , c("PC1", "PC2")
    ]

    x_y_ranges <- find_x_y_ranges(cell_subset_projection_coords, extend_ratio)

    mesh_line_obj_2d <- create_2d_mesh_grid_line_obj(
        x_y_ranges$x_min,
        x_y_ranges$x_max,
        x_y_ranges$y_min,
        x_y_ranges$y_max,
        mesh_grid_n
    )
    center_vec <- apply(
        cell_subset_projection_locations, 2, function(x) mean(x)
    )
    mesh_grid_line_obj_3d <- transform_2d_mesh_to_3d_mesh(
        mesh_line_obj_2d,
        cell_subset_projection_PCA,
        center_vec,
        mesh_grid_n
    )
    return(mesh_grid_line_obj_3d)
}


# cross section creation function ####

#' @title createCrossSection
#' @description Create a virtual 2D cross section.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param spat_loc_name name of spatial locations
#' @param name name of cress section object. (default = cross_sectino)
#' @param spatial_network_name name of spatial network object.
#' (default = Delaunay_network)
#' @param thickness_unit unit of the virtual section thickness. If "cell",
#' average size of the observed cells is used as length unit. If "natural",
#' the unit of cell location coordinates is used. (default = cell)
#' @param slice_thickness thickness of slice. default = 2
#' @param cell_distance_estimate_method method to estimate average distance
#' between neighobring cells. (default = mean)
#' @param extend_ratio deciding the span of the cross section meshgrid, as a
#' ratio of extension compared to the borders of the vitural tissue section.
#' (default = 0.2)
#' @param method method to define the cross section plane.
#' If equation, the plane is defined by a four element numerical vector
#' (equation) in the form of c(A,B,C,D), corresponding to a plane with
#' equation Ax+By+Cz=D.
#' If 3 points, the plane is define by the coordinates of 3 points, as given by
#' point1, point2, and point3.
#' If point and norm vector, the plane is defined by the coordinates of one
#' point (point1) in the plane and the coordinates of one norm vector
#' (normVector) to the plane.
#' If point and two plane vector, the plane is defined by the coordinates of
#' one point (point1) in the plane and the coordinates of two vectors
#' (planeVector1, planeVector2) in the plane.
#' (default = equation)
#' @param equation equation required by method "equation".equations needs to be
#' a numerical vector of length 4, in the form of c(A,B,C,D), which defines
#' plane Ax+By+Cz=D.
#' @param point1 coordinates of the first point required by method
#' "3 points","point and norm vector", and "point and two plane vectors".
#' @param point2 coordinates of the second point required by method "3 points"
#' @param point3 coordinates of the third point required by method "3 points"
#' @param normVector coordinates of the norm vector required by method
#' "point and norm vector"
#' @param planeVector1 coordinates of the first plane vector required by
#' method "point and two plane vectors"
#' @param planeVector2 coordinates of the second plane vector required by
#' method "point and two plane vectors"
#' @param mesh_grid_n numer of meshgrid lines to generate along both directions
#' for the cross section plane.
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @returns giotto object with updated spatial network slot
#' @details Creates a virtual 2D cross section object for a given spatial
#' network object. The users need to provide the definition of the cross
#' section plane (see method).
#' @examples
#' g <- GiottoData::loadGiottoMini("starmap")
#'
#' g <- createCrossSection(
#'     gobject = g,
#'     method = "equation",
#'     equation = c(0, 1, 0, 600),
#'     extend_ratio = 0.6,
#'     name = "new_cs",
#'     return_gobject = TRUE
#' )
#'
#' crossSectionPlot(g, name = "new_cs")
#' @export
createCrossSection <- function(
        gobject,
        spat_unit = NULL,
        spat_loc_name = "raw",
        name = "cross_section",
        spatial_network_name = "Delaunay_network",
        thickness_unit = c("cell", "natural"),
        slice_thickness = 2,
        cell_distance_estimate_method = "mean",
        extend_ratio = 0.2,
        method = c(
            "equation", "3 points", "point and norm vector",
            "point and two plane vectors"
        ),
        equation = NULL,
        point1 = NULL, point2 = NULL, point3 = NULL,
        normVector = NULL,
        planeVector1 = NULL, planeVector2 = NULL,
        mesh_grid_n = 20,
        return_gobject = TRUE,
        verbose = NULL) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )

    # read spatial locations
    spatial_locations <- getSpatialLocations(
        gobject = gobject, spat_unit = spat_unit, name = spat_loc_name,
        set_defaults = FALSE, verbose = FALSE, output = "spatLocsObj"
    )

    spatial_locations <- as.matrix(spatial_locations, id_rownames = TRUE)
    cell_ID_vec <- seq_len(nrow(spatial_locations))
    names(cell_ID_vec) <- rownames(spatial_locations)

    # generate section plane equation

    method <- match.arg(
        method,
        c(
            "equation", "3 points", "point and norm vector",
            "point and two plane vectors"
        )
    )

    switch(method,
        "equation" = {
            if (is.null(equation)) {
                message("equation was not provided.")
            } else {
                plane_equation <- equation
                plane_equation[4] <- -equation[4]
            }
        },
        "point and norm vector" = {
            if (is.null(point1) || is.null(normVector)) {
                message("either point or norm vector was not provided.")
            } else {
                plane_equation <- c()
                plane_equation[seq_len(3)] <- normVector
                plane_equation[4] <- -point1 %*% normVector
            }
        },
        "point and two plane vectors" = {
            if (is.null(point1) ||
                is.null(planeVector1) ||
                is.null(planeVector2)) {
                message("either point or any of the two plane vectors was not
                    provided.")
            } else {
                normVector <- crossprod(planeVector1, planeVector2)
                plane_equation[seq_len(3)] <- normVector
                plane_equation[4] <- -point1 %*% normVector
            }
        },
        "3 points" = {
            if (is.null(point1) || is.null(point2) || is.null(point3)) {
                message("not all three points were provided.")
            } else {
                planeVector1 <- point2 - point1
                planeVector2 <- point3 - point1
                normVector <- crossprod(planeVector1, planeVector2)
                plane_equation[seq_len(3)] <- normVector
                plane_equation[4] <- -point1 %*% normVector
            }
        }
    )

    names(plane_equation) <- c("A", "B", "C", "D")

    # determine section thickness
    thickness_unit <- match.arg(thickness_unit, c("cell", "natural"))
    sectionThickness <- get_sectionThickness(gobject,
        spat_unit = spat_unit,
        thickness_unit = thickness_unit,
        slice_thickness = slice_thickness,
        spatial_network_name = spatial_network_name,
        cell_distance_estimate_method = cell_distance_estimate_method,
        plane_equation = plane_equation
    )

    max_distance_to_section_plane <- sectionThickness / 2

    # calculate distances to cross section
    spatial_locations_mat <- cbind(
        spatial_locations, as.matrix(rep(1, dim(spatial_locations)[1]))
    )
    norm_vec <- function(x) sqrt(sum(x^2))
    distance_to_plane_vector <- abs(spatial_locations_mat %*% as.matrix(
        plane_equation
    ) / norm_vec(plane_equation[1:3]))

    # select cells within section ###
    cell_subset <- distance_to_plane_vector <= max_distance_to_section_plane

    # project the selected cells onto the section plane ###
    cell_subset_spatial_locations <- spatial_locations[cell_subset, ]

    ## find a point on the section plane ##
    if (plane_equation["A"] != 0) {
        plane_point <- c(-plane_equation["D"] / plane_equation["A"], 0, 0)
    } else if (plane_equation["B"] != 0) {
        plane_point <- c(0, -plane_equation["D"] / plane_equation["B"], 0)
    } else if (plane_equation["C"] != 0) {
        plane_point <- c(0, 0, -plane_equation["D"] / plane_equation["C"])
    }
    ## find the projection Xp,Yp,Zp coordinates ##
    cell_subset_projection_locations <- t(apply(
        cell_subset_spatial_locations, 1,
        function(x) {
            projection_fun(x,
                plane_point = plane_point,
                plane_norm = plane_equation[1:3]
            )
        }
    ))

    # get the local coordinates of selected cells on the section plane
    cell_subset_projection_PCA <- stats::prcomp(
        cell_subset_projection_locations
    )
    cell_subset_projection_coords <- get_cross_section_coordinates(
        cell_subset_projection_locations
    )

    # create mesh grid lines for the cross section ###
    mesh_grid_lines <- create_mesh_grid_lines(
        cell_subset_projection_locations, extend_ratio, mesh_grid_n
    )
    mesh_obj <- list("mesh_grid_lines" = mesh_grid_lines)

    ### save and update the spatial object ###

    crossSection_obj <- create_crossSection_object(
        method = method,
        thickness_unit = thickness_unit,
        slice_thickness = slice_thickness,
        cell_distance_estimate_method = cell_distance_estimate_method,
        extend_ratio = extend_ratio,
        plane_equation = plane_equation, mesh_grid_n = mesh_grid_n,
        mesh_obj = mesh_obj, cell_subset = cell_subset,
        cell_subset_spatial_locations = cell_subset_spatial_locations,
        cell_subset_projection_locations = cell_subset_projection_locations,
        cell_subset_projection_PCA = cell_subset_projection_PCA,
        cell_subset_projection_coords = cell_subset_projection_coords
    )


    if (return_gobject) {
        sn <- getSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            name = spatial_network_name,
            copy_obj = FALSE,
            set_defaults = FALSE,
            verbose = FALSE,
            output = "spatialNetworkObj"
        )

        cs_names <- names(sn@crossSectionObjects)
        if (name %in% cs_names) {
            vmsg(.v = verbose, sprintf(
                "name '%s' has already been used, will be overwritten",
                name
            ))
        }

        sn@crossSectionObjects[[name]] <- crossSection_obj

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, sn, verbose = FALSE)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        return(gobject)
    } else {
        return(crossSection_obj)
    }
}


# cross section visual functions ####

####
#' @title crossSectionFeatPlot
#' @name crossSectionFeatPlot
#' @description Visualize cells and feature expression in a virtual cross
#' section according to spatial coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param spat_loc_name name of spatial locations
#' @param crossSection_obj crossSection object
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @param ... parameters for spatFeatPlot2D
#' @returns ggplot
#' @details Description of parameters.
#' @md
#' @seealso [GiottoVisuals::spatGenePlot3D] and [GiottoVisuals::spatFeatPlot2D]
#' @export
crossSectionFeatPlot <- function(gobject = NULL,
    spat_unit = NULL,
    feat_type = NULL,
    spat_loc_name = "raw",
    crossSection_obj = NULL,
    name = NULL,
    spatial_network_name = "Delaunay_network",
    default_save_name = "crossSectionGenePlot",
    ...) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    # load cross section object
    if (is.null(crossSection_obj)) {
        crossSection_obj <- read_crossSection(
            gobject,
            spat_unit = spat_unit,
            name = name,
            spatial_network_name = spatial_network_name
        )
    }

    cell_subset <- crossSection_obj$cell_subset
    cell_subset <- rownames(cell_subset)[which(cell_subset)]
    cell_subset_projection_coords <-
        crossSection_obj$cell_subset_projection_coords

    # modify gobject based on crossSection object
    gobj_sids <- spatIDs(gobject, spat_unit = spat_unit)
    subset_cell_ids <- gobj_sids[gobj_sids %in% cell_subset]
    temp_gobject <- subsetGiotto(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        cell_ids = subset_cell_ids
    )

    sl <- getSpatialLocations(
        gobject = temp_gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "spatLocsObj",
        copy_obj = TRUE,
        verbose = FALSE
    )

    sl[]$sdimx <- cell_subset_projection_coords[, 1]
    sl[]$sdimy <- cell_subset_projection_coords[, 2]
    sl[]$sdimz <- rep(0, dim(cell_subset_projection_coords)[1])

    temp_gobject <- setGiotto(temp_gobject, sl, verbose = FALSE)

    # call spatFeatPlot2D to generate the plots
    GiottoVisuals::spatFeatPlot2D(
        gobject = temp_gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spatial_network_name = spatial_network_name,
        default_save_name = default_save_name,
        ...
    )
}
####

#' @title crossSectionPlot
#' @name crossSectionPlot
#' @description Visualize cells in a virtual cross section according to
#' spatial coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param spat_loc_name name of spatial locations
#' @param crossSection_obj cross section object as alternative input.
#' default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @param ... parameters for spatPlot2D
#' @returns ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{crossSectionPlot}}
crossSectionPlot <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spat_loc_name = "raw",
        crossSection_obj = NULL,
        name = NULL,
        spatial_network_name = "Delaunay_network",
        default_save_name = "crossSectionPlot",
        ...) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    # load cross section object
    if (is.null(crossSection_obj)) {
        crossSection_obj <- read_crossSection(
            gobject = gobject,
            spat_unit = spat_unit,
            name = name,
            spatial_network_name = spatial_network_name
        )
    }

    cell_subset <- crossSection_obj$cell_subset
    cell_subset <- rownames(cell_subset)[which(cell_subset)]
    cell_subset_projection_coords <-
        crossSection_obj$cell_subset_projection_coords

    # modify gobject based on crossSection object
    gobj_sids <- spatIDs(gobject, spat_unit = spat_unit)
    subset_cell_ids <- gobj_sids[gobj_sids %in% cell_subset]
    temp_gobject <- subsetGiotto(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        cell_ids = subset_cell_ids
    )

    sl <- getSpatialLocations(
        gobject = temp_gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "spatLocsObj",
        copy_obj = TRUE,
        verbose = FALSE
    )

    sl[]$sdimx <- cell_subset_projection_coords[, 1]
    sl[]$sdimy <- cell_subset_projection_coords[, 2]
    sl[]$sdimz <- rep(0, dim(cell_subset_projection_coords)[1])

    temp_gobject <- setGiotto(temp_gobject, sl, verbose = FALSE)

    # call spatFeatPlot2D to generate the plots
    spatPlot2D(
        gobject = temp_gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spatial_network_name = spatial_network_name,
        default_save_name = default_save_name,
        ...
    )
}

####
#' @title crossSectionFeatPlot3D
#' @name crossSectionFeatPlot3D
#' @description Visualize cells and feature expression in a virtual cross
#' section according to spatial coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param crossSection_obj cross section object as alternative input.
#' default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param show_other_cells logical. Default = TRUE
#' @param other_cell_color color of cells outside the cross section.
#' default = transparent.
#' @param default_save_name default save name for saving, don't change, change
#' save_name in save_param
#' @param ... parameters for spatGenePlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
crossSectionFeatPlot3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        crossSection_obj = NULL,
        name = NULL,
        spatial_network_name = "Delaunay_network",
        show_other_cells = TRUE,
        other_cell_color = alpha("lightgrey", 0),
        default_save_name = "crossSectionGenePlot3D",
        ...) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    # load cross section object
    if (is.null(crossSection_obj)) {
        crossSection_obj <- read_crossSection(
            gobject = gobject,
            spat_unit = spat_unit,
            name = name,
            spatial_network_name = spatial_network_name
        )
    }

    cell_subset <- crossSection_obj$cell_subset
    cell_subset <- rownames(cell_subset)[which(cell_subset)]
    cell_subset_projection_coords <-
        crossSection_obj$cell_subset_projection_coords

    # modify gobject based on crossSection object
    gobj_sids <- spatIDs(gobject, spat_unit = spat_unit)
    subset_cell_ids <- gobj_sids[gobj_sids %in% cell_subset]

    # call spatGenePlot3D to generate the plots
    spatFeatPlot3D(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        select_cells = subset_cell_ids,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        default_save_name = default_save_name,
        ...
    )
}
####
#' @title crossSectionPlot3D
#' @name crossSectionPlot3D
#' @description Visualize cells in a virtual cross section according to spatial
#' coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param crossSection_obj cross section object as alternative input.
#' default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of cells outside the cross section.
#' default = transparent.
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @param ... parameters for spatPlot3D
#' @returns ggplot
#' @details Description of parameters.
#' @export
crossSectionPlot3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        crossSection_obj = NULL,
        name = NULL,
        spatial_network_name = "Delaunay_network",
        show_other_cells = TRUE,
        other_cell_color = alpha("lightgrey", 0),
        default_save_name = "crossSection3D",
        ...) {
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    # load cross section object
    if (is.null(crossSection_obj)) {
        crossSection_obj <- read_crossSection(
            gobject,
            spat_unit = spat_unit,
            name = name,
            spatial_network_name = spatial_network_name
        )
    }

    cell_subset <- crossSection_obj$cell_subset
    cell_subset <- rownames(cell_subset)[which(cell_subset)]
    cell_subset_projection_coords <-
        crossSection_obj$cell_subset_projection_coords

    # modify gobject based on crossSection object
    gobj_sids <- spatIDs(gobject, spat_unit = spat_unit)
    subset_cell_ids <- gobj_sids[gobj_sids %in% cell_subset]

    # call spatPlot3D to generate the plots
    spatPlot3D(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        select_cells = subset_cell_ids,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        default_save_name = default_save_name,
        ...
    )
}


####
#' @title insertCrossSectionSpatPlot3D
#' @name insertCrossSectionSpatPlot3D
#' @description Visualize the meshgrid lines of cross section together with
#' cells
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param spat_loc_name name of spatial locations
#' @param crossSection_obj cross section object as alternative input.
#' default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param mesh_grid_color color for the meshgrid lines
#' @param mesh_grid_width width for the meshgrid lines
#' @param mesh_grid_style style for the meshgrid lines
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param show_other_cells display not selected cells
#' @param axis_scale axis_scale
#' @param custom_ratio custom_ratio
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @param ... parameters for spatPlot3D
#' @returns ggplot
#' @details Description of parameters.
#' @export
insertCrossSectionSpatPlot3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spat_loc_name = "raw",
        crossSection_obj = NULL,
        name = NULL,
        spatial_network_name = "Delaunay_network",
        mesh_grid_color = "#1f77b4",
        mesh_grid_width = 3,
        mesh_grid_style = "dot",
        sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
        show_other_cells = FALSE,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        default_save_name = "spat3D_with_cross_section",
        ...) {
    package_check("plotly", repository = "CRAN:plotly")

    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    # load cross section object
    if (!is.null(crossSection_obj)) {
        crossSection_obj <- crossSection_obj
    } else {
        crossSection_obj <- read_crossSection(
            gobject,
            spat_unit = spat_unit,
            name = name,
            spatial_network_name = spatial_network_name
        )
    }


    pl <- spatPlot3D(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        sdimx = sdimx,
        sdimy = sdimy,
        sdimz = sdimz,
        show_other_cells = show_other_cells,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE,
        default_save_name = default_save_name,
        ...
    )

    for (i in seq_len(dim(
        crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X
    )[2])) {
        pl <- pl %>% plotly::add_trace(
            x = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X[, i],
            y = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Y[, i],
            z = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Z[, i],
            mode = "lines", type = "scatter3d",
            line = list(
                color = mesh_grid_color,
                width = mesh_grid_width, dash = mesh_grid_style
            )
        )
    }

    sl <- getSpatialLocations(
        gobject = gobject, spat_unit = spat_unit, name = spat_loc_name,
        output = "data.table", copy_obj = TRUE, verbose = FALSE,
        set_defaults = TRUE
    )

    current_ratio <- plotly_axis_scale_3D(
        cell_locations = sl,
        sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
        mode = axis_scale, custom_ratio = custom_ratio
    )

    new_ratio <- adapt_aspect_ratio(
        current_ratio = current_ratio,
        cell_locations = sl,
        sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
        mesh_obj = crossSection_obj$mesh_obj
    )

    pl <- pl %>% plotly::layout(
        showlegend = FALSE,
        scene = list(
            aspectmode = "manual",
            aspectratio = list(
                x = new_ratio[[1]],
                y = new_ratio[[2]],
                z = new_ratio[[3]]
            )
        )
    )

    return(pl)
}
####
#' @title insertCrossSectionFeatPlot3D
#' @name insertCrossSectionFeatPlot3D
#' @description Visualize cells and gene expression in a virtual cross section
#' according to spatial coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param spat_loc_name name of spatial locations
#' @param crossSection_obj cross section object as alternative input.
#' default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param mesh_grid_color color for the meshgrid lines
#' @param mesh_grid_width width for the meshgrid lines
#' @param mesh_grid_style style for the meshgrid lines
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param show_other_cells display not selected cells
#' @param axis_scale axis_scale
#' @param custom_ratio custom_ratio
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot logical. directly save the plot
#' @param save_param list of saving parameters from
#' [GiottoVisuals::all_plots_save_function]
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @param ... parameters for spatGenePlot3D
#' @returns ggplot
#' @details Description of parameters.
#' @md
#' @export
insertCrossSectionFeatPlot3D <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    spat_loc_name = "raw",
    crossSection_obj = NULL,
    name = NULL,
    spatial_network_name = "Delaunay_network",
    mesh_grid_color = "#1f77b4",
    mesh_grid_width = 3,
    mesh_grid_style = "dot",
    sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
    show_other_cells = FALSE,
    axis_scale = c("cube", "real", "custom"),
    custom_ratio = NULL,
    show_plot = NULL, return_plot = NULL, save_plot = NULL,
    save_param = list(),
    default_save_name = "spatGenePlot3D_with_cross_section",
    ...) {
    package_check("plotly", repository = "CRAN:plotly")

    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    # load cross section object
    if (is.null(crossSection_obj)) {
        crossSection_obj <- read_crossSection(
            gobject,
            spat_unit = spat_unit,
            name = name,
            spatial_network_name = spatial_network_name
        )
    }

    pl <- spatFeatPlot3D(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        show_other_cells = FALSE,
        axis_scale = axis_scale,
        custom_ratio = custom_ratio,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE,
        default_save_name = default_save_name,
        ...
    )

    for (i in seq_len(dim(
        crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X
    )[2])) {
        pl <- pl %>% plotly::add_trace(
            x = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X[, i],
            y = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Y[, i],
            z = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Z[, i],
            mode = "lines+markers", type = "scatter3d", color = mesh_grid_color,
            marker = list(color = alpha(mesh_grid_color, 0)),
            line = list(
                color = mesh_grid_color,
                width = mesh_grid_width, dash = mesh_grid_style
            )
        )
    }

    sl <- getSpatialLocations(
        gobject = gobject, spat_unit = spat_unit, name = spat_loc_name,
        output = "data.table", copy_obj = TRUE, verbose = FALSE,
        set_defaults = TRUE
    )


    current_ratio <- plotly_axis_scale_3D(
        cell_locations = sl,
        sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
        mode = axis_scale, custom_ratio = custom_ratio
    )

    new_ratio <- adapt_aspect_ratio(
        current_ratio = current_ratio,
        cell_locations = sl,
        sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
        mesh_obj = crossSection_obj$mesh_obj
    )

    pl <- pl %>% plotly::layout(
        showlegend = FALSE,
        scene = list(
            aspectmode = "manual",
            aspectratio = list(
                x = new_ratio[[1]],
                y = new_ratio[[2]],
                z = new_ratio[[3]]
            )
        )
    )

    cowplot <- pl

    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = cowplot,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}

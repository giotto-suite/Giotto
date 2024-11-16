# interpolateFeature ####
setGeneric(
    "interpolateFeature",
    function(x, y, ...) standardGeneric("interpolateFeature")
)

#' @name interpolateFeature
#' @title Spatial feature interpolation
#' @param x object containing coordinates to use interpolation with
#' @param y data.frame-like. Values for interpolation. Must also have a
#' `cell_ID` column and that matches with `x`.
#' @param verbose be verbose
#' @param ext `SpatExtent`. (optional) extent across which to apply the
#' interpolation. If not provided, will default to the extent of the spatLocsObj
#' expanded by the value of `buffer`. It can be helpful to set this as the
#' extent of any polygons that will be used in aggregation.
#' @param buffer numeric. (optional) default buffer to expand derived extent by
#'  if `ext` is not provided.
#' @param ... additional params to pass downstream methods
NULL






#' @rdname interpolateFeature
#' @param spat_unit (optional) spatial unit to use
#' @param feat_type (optional) feature type to use
#' @param feats character vector. Features to interpolate from the `giotto`
#' object
#' @param spatvalues_params list. Additional list of parameters to pass to
#' [spatValues()] to help with data retrieval from `giotto` object
#' @param spat_loc_name character. Name of spatial locations to use. Values to
#' be interpolated are spatially mapped to these locations by cell_ID.
#' @param name_fmt character. sprintf fmt to apply to `feats` when naming the
#' resulting interpolation `giottoLargeImage` objects. Default is no change.
#' @param savedir character. Output directory. Default is a new `interp_rasters`
#' folder in working directory.
#' @returns
#' `giotto` method returns a `giotto` object with newly made appended
#' feature interpolation rasters as `giottoLargeImages`\cr
#' @export
setMethod(
    "interpolateFeature", signature(x = "giotto", y = "missing"),
    function(
        x,
        spat_unit = NULL,
        feat_type = NULL,
        feats,
        spatvalues_params = list(),
        spat_loc_name = "raw",
        ext = NULL,
        buffer = 50,
        name_fmt = "%s",
        savedir = file.path(getwd(), "interp_rasters"),
        overwrite = FALSE,
        verbose = NULL,
        ...) {
        sl <- NULL

        # This method prepares the data from the giotto object to pass
        # downstream where the actual interpolation happens

        # prep params to pass to lower level methods
        a <- list(...)
        a$name_fmt <- name_fmt
        a$savedir <- savedir
        a$overwrite <- overwrite
        a$verbose <- verbose


        # defaults
        spat_unit <- set_default_spat_unit(
            gobject = x,
            spat_unit = spat_unit
        )
        feat_type <- set_default_feat_type(
            gobject = x,
            spat_unit = spat_unit,
            feat_type = feat_type
        )

        # get spatlocs
        a$x <- getSpatialLocations(
            gobject = x,
            spat_unit = spat_unit,
            name = spat_loc_name,
            output = "spatLocsObj",
            copy_obj = FALSE,
            verbose = FALSE,
            set_defaults = FALSE
        )

        # extent to interpolate across
        if (is.null(ext)) {
            e <- ext(sl)
            ext <- ext(e[] + rep(c(-buffer, buffer), 2L))
        } else {
            checkmate::assert_class(ext, "SpatExtent")
        }
        a$ext <- ext # append to argslist

        # output directory
        if (!dir.exists(savedir)) dir.create(savedir, recursive = TRUE)


        # params and call for data retrieval
        spatvalues_params$gobject <- x
        spatvalues_params$spat_unit <- spat_unit
        spatvalues_params$feat_type <- feat_type
        spatvalues_params$feats <- feats
        spatvalues_params$verbose <- verbose
        a$y <- do.call(spatValues, args = spatvalues_params)


        # call data.frame method
        # Should return a giottoLargeImage of the interpolated raster
        # Note that this object will need to be reconnected.
        interp_img_list <- do.call(interpolateFeature, args = a)

        for (i in seq_along(interp_img_list)) {
            x <- setGiotto(
                gobject = x,
                x = interp_img_list[[i]]
            )
        }
        return(x)
    }
)


#' @rdname interpolateFeature
#' @param rastersize numeric. Length of major axis in px of interpolation
#' raster to create.
#' @param overwrite logical. Whether raster outputs should be overwritten if
#' the same `filename` is provided.
#' @details
#' The data.frame method returns a `giottoLargeImage` linked to an interpolated
#' raster that is written to disk as GeoTIFF.
#' @export
setMethod(
    "interpolateFeature",
    signature(x = "spatLocsObj", y = "data.frame"),
    function(
        x, y,
        ext = NULL,
        buffer = 50,
        rastersize = 500,
        name_fmt = "%s",
        savedir = file.path(getwd(), "interp_rasters"),
        overwrite = FALSE,
        # cores = GiottoUtils::determine_cores(),
        ...) {
        checkmate::assert_character(savedir)
        checkmate::assert_character(name_fmt)
        checkmate::assert_logical(overwrite)
        # checkmate::assert_numeric(cores)
        package_check("gstat", repository = "CRAN:gstat")

        # output directory
        if (!dir.exists(savedir)) dir.create(savedir, recursive = TRUE)

        # get features to iterate across
        y <- data.table::as.data.table(y)
        val_cols <- colnames(y)
        feats <- val_cols[!val_cols %in% c("cell_ID")]

        # combine data to ensure that spatlocs and values are matched
        # This is the data that will be fed to terra::interpolate()
        annotatedlocs <- data.table::merge.data.table(x[], y, by = "cell_ID")

        # extent to interpolate across
        if (is.null(ext)) {
            e <- ext(x)
            ext <- ext(e[] + rep(c(-buffer, buffer), 2L))
        } else {
            checkmate::assert_class(ext, "SpatExtent")
        }
        # create numeric representation since SpatExtent is not passable to
        # workers
        e_numeric <- ext[]


        with_pbar({
            pb <- pbar(along = feats)

            interp_img_list <- lapply_flex(
                feats,
                function(feat) {
                    name <- sprintf(name_fmt, feat)
                    filename <- file.path(savedir, paste0(name, ".tif"))

                    # create subset table with only relevant data
                    data <- annotatedlocs[,
                        c("cell_ID", feat, "sdimx", "sdimy"),
                        with = FALSE
                    ]
                    data.table::setnames(data, old = feat, new = "count")

                    # model to use
                    model <- gstat::gstat(
                        id = feat,
                        formula = count ~ 1,
                        locations = ~ sdimx + sdimy,
                        data = data,
                        nmax = 7,
                        set = list(
                            idp = 0.5
                        )
                    )

                    # regenerate extent
                    interp_e <- terra::ext(e_numeric)

                    # generate raster to interpolate on
                    res <- max(range(interp_e) / rastersize[1L])
                    r <- terra::rast(x = interp_e, res = res)

                    # perform interpolation and save to disk
                    terra::interpolate(
                        object = r,
                        model = model,
                        xyNames = c("sdimx", "sdimy"),
                        debug.level = 0,
                        index = 1L,
                        filename = filename,
                        overwrite = overwrite,
                        # cores = cores,
                        # cpkgs = c("data.table", "gstat"),
                        wopt = list(
                            progress = FALSE,
                            filetype = "COG"
                        )
                    )

                    gimg <- createGiottoLargeImage(
                        raster_object = filename,
                        name = name,
                        use_rast_ext = TRUE,
                        verbose = FALSE
                    )

                    # report progress
                    pb(message = c("interpolating: ", feat))

                    # since this call may be parallelized, the returned
                    # image object is expected to be disconnected
                    return(gimg)
                },
                future.seed = TRUE,
                future.globals = list(
                    e_numeric = e_numeric,
                    annotatedlocs = annotatedlocs,
                    name_fmt = name_fmt,
                    savedir = savedir,
                    # cores = cores,
                    overwrite = overwrite
                ),
                future.packages = c(
                    "terra", "gstat", "data.table"
                )
            )
        })

        # reconnect images
        interp_img_list <- lapply(interp_img_list, GiottoClass::reconnect)
        return(interp_img_list)
    }
)

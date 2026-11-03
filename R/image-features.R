# image features ####
#
# Per-polygon features from an image: a square tile of the image centred on
# each polygon (at the image's own resolution) and the polygon's mask in that
# tile go to a provider, which returns one feature vector per polygon. The
# features are stored as the expression of the polygons' spatial unit, so
# PCA, clustering and runNucleusDeconv() can use them.
#
# Providers are kept ad hoc for now and shaped like the param objects
# (`method` + `param` + `fun`): `fun(tiles, masks, um_per_px, polys)` takes
#   tiles      integer array, n x tile_px x tile_px x 3 (RGB, 0-255; outside
#              the image is white)
#   masks      logical array, n x tile_px x tile_px: the polygon's pixels
#   um_per_px  the image's resolution in the polygons' units
#   polys      the n polygons (SpatVector, polygon units)
# and returns an n x features matrix. A provider resizes and normalises the
# tiles itself, since that is part of what its model expects.


#' @title Image features of polygons
#' @name calculateImageFeatures
#' @description
#' Compute features for every polygon of a spatial unit, e.g. nuclei, from an
#' image, e.g. the H&E stain: a square tile of the image centred on each
#' polygon, with the polygon's mask, goes to a feature provider. The features
#' are stored as the expression of the spatial unit (features by polygons),
#' where [runPCA()], clustering or [runNucleusDeconv()] can use them.
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell"); its polygons are featurised
#' @param image_name name of the image (`giottoLargeImage` or
#'   `giottoAffineImage`), RGB or single channel
#' @param method the feature provider:
#'   * `"uni2"`: UNI2-h, a pathology foundation model (Chen et al. 2024),
#'     through Python. Each tile is resized to 224 px (bicubic), and the
#'     model's last-layer patch tokens are averaged over the polygon's mask,
#'     weighted by the mask's area in each patch. 1536 features. Needs
#'     `torch`, `timm` and `Pillow`, and access to the gated
#'     `MahmoodLab/UNI2-h` model on Hugging Face.
#'   * `"morphology"`: R only. Shape (area, perimeter, circularity,
#'     elongation) and stain (mean and sd of each channel inside the polygon,
#'     mean optical density) features.
#'   * a function `function(tiles, masks, um_per_px, polys)` returning a
#'     polygons by features matrix (see Details).
#' @param tile_um numeric. side of the tile, in the polygons' units
#'   (microns)
#' @param feat_type feature type to store the features under
#' @param name name of the features. `NULL` (default) uses the method's
#'   name (`"custom"` for a function).
#' @param python_env Python environment for `"uni2"`, as for
#'   [doStardistSegmentation()]
#' @param block_px integer. tiles are cut from blocks of about this many
#'   image pixels square, so the image is never read whole.
#' @param batch_size integer. tiles per call of the provider
#' @param verbose be verbose
#' @returns the giotto object with expression `name` of `spat_unit` and
#'   `feat_type`
#' @details
#' A function given as `method` is called for each batch of polygons with
#' `tiles` (integer array, polygons x pixels x pixels x 3, RGB 0-255, white
#' outside the image), `masks` (logical array, polygons x pixels x pixels),
#' `um_per_px` (the image's resolution) and `polys` (the polygons, a
#' `SpatVector`). Tiles are at the image's own resolution, centred on the
#' pixel that holds the polygon's centroid, with an odd side so that pixel is
#' the centre.
#' @seealso [runNucleusDeconv()]
#' @md
#' @export
calculateImageFeatures <- function(gobject,
        spat_unit,
        image_name,
        method = c("uni2", "morphology"),
        tile_um = 112,
        feat_type = "img",
        name = NULL,
        python_env = NULL,
        block_px = 2048L,
        batch_size = 64L,
        verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_character(spat_unit, len = 1L)
    checkmate::assert_number(tile_um, lower = 0)
    provider <- .img_provider(method, python_env = python_env)
    name <- name %null% provider$method

    polys <- getPolygonInfo(gobject, name = spat_unit,
        return_giottoPolygon = TRUE, verbose = FALSE)
    img <- getGiottoImage(gobject, name = image_name)
    feats <- .img_features(polys, img, provider, tile_um = tile_um,
        block_px = block_px, batch_size = batch_size, verbose = verbose)

    ex <- createExprObj(t(feats),
        name = name,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = spat_unit,
        misc = list(method = provider$method, image = image_name,
            tile_um = tile_um)
    )
    setGiotto(gobject, ex, verbose = verbose)
}


# providers ####

.img_provider <- function(method = c("uni2", "morphology"), python_env = NULL) {
    if (is.function(method)) {
        return(list(method = "custom", param = list(), fun = method))
    }
    method <- match.arg(method)
    switch(method,
        "uni2" = .img_provider_uni2(python_env = python_env),
        "morphology" = .img_provider_morphology()
    )
}

# UNI2-h through reticulate. The model is loaded once per session, on the
# first batch.
.img_provider_uni2 <- function(python_env = NULL) {
    mod <- NULL
    list(
        method = "uni2",
        param = list(python_env = python_env),
        fun = function(tiles, masks, um_per_px, polys) {
            if (is.null(mod)) {
                set_giotto_python_path(python_env)
                package_check(c("torch", "timm", "PIL"),
                    repository = c("pip:torch", "pip:timm", "pip:Pillow"))
                mod <<- reticulate::import_from_path("image_features_uni2",
                    path = system.file("python", package = "Giotto"))
                mod$load()
            }
            storage.mode(tiles) <- "integer"
            out <- mod$embed(
                reticulate::r_to_py(tiles)$astype("uint8"),
                reticulate::r_to_py(masks)
            )
            matrix(as.numeric(out), nrow = dim(tiles)[1L])
        }
    )
}

# Shape from the polygons, stain from the pixels under the mask.
.img_provider_morphology <- function() {
    list(
        method = "morphology",
        param = list(),
        fun = function(tiles, masks, um_per_px, polys) {
            polys <- .img_planar(polys)
            area <- terra::expanse(polys, transform = FALSE)
            perim <- terra::perim(polys)
            elong <- vapply(seq_len(nrow(polys)), function(i) {
                # vertices once each (a ring repeats its first at the end)
                xy <- unique(terra::geom(polys[i])[, c("x", "y"), drop = FALSE])
                ev <- eigen(stats::cov(xy), symmetric = TRUE,
                    only.values = TRUE)$values
                sqrt(max(ev[1L], 0) / max(ev[2L], 1e-12))
            }, numeric(1L))
            n <- dim(tiles)[1L]
            stain <- t(vapply(seq_len(n), function(i) {
                m <- masks[i, , ]
                if (!any(m)) m[] <- TRUE
                px <- vapply(1:3, function(k) tiles[i, , , k][m],
                    numeric(sum(m)))
                px <- matrix(px, ncol = 3L)
                od <- -log10(pmax(px, 1) / 255)
                c(colMeans(px), apply(px, 2L, stats::sd),
                    mean(rowSums(od)))
            }, numeric(7L)))
            stain[is.na(stain)] <- 0
            out <- cbind(area, perim, 4 * pi * area / pmax(perim, 1e-12)^2,
                elong, stain)
            colnames(out) <- c("area", "perimeter", "circularity",
                "elongation", "mean_r", "mean_g", "mean_b", "sd_r", "sd_g",
                "sd_b", "mean_od")
            out
        }
    )
}


# tiles ####

# Features of every polygon of `polys` (giottoPolygon) from `img`: a
# polygons x features matrix, rows named by polygon ID.
.img_features <- function(polys, img, provider, tile_um, block_px = 2048L,
        batch_size = 64L, verbose = NULL) {
    fr <- .img_frame(polys, img)
    r <- fr$raster
    pr <- fr$polys
    um_per_px <- fr$um_per_raster_unit * terra::res(r)[1L]
    tn <- as.integer(round(tile_um / um_per_px))
    if (tn %% 2L == 0L) tn <- tn + 1L
    half <- tn %/% 2L
    ids <- polys[]$poly_ID
    vmsg(.v = verbose, sprintf(
        "[calculateImageFeatures] %d polygons, tiles of %d px (%.3g per px)",
        length(ids), tn, um_per_px))

    xy <- terra::crds(terra::centroids(pr))
    row <- terra::rowFromY(r, xy[, 2L])
    col <- terra::colFromX(r, xy[, 1L])
    # a centroid outside the image: the nearest edge pixel, the tile is
    # mostly white
    row[is.na(row)] <- ifelse(xy[is.na(row), 2L] > terra::ymax(r), 1L,
        terra::nrow(r))
    col[is.na(col)] <- ifelse(xy[is.na(col), 1L] < terra::xmin(r), 1L,
        terra::ncol(r))
    pr$.idx <- seq_along(ids)

    block <- paste((row - 1L) %/% block_px, (col - 1L) %/% block_px)
    out <- NULL
    for (b in split(seq_along(ids), block)) {
        cv <- .img_canvas(r, pr, rows = range(row[b]) + c(-half, half),
            cols = range(col[b]) + c(-half, half))
        for (s0 in seq(1L, length(b), by = batch_size)) {
            k <- b[s0:min(s0 + batch_size - 1L, length(b))]
            tiles <- array(0L, c(length(k), tn, tn, 3L))
            masks <- array(FALSE, c(length(k), tn, tn))
            for (j in seq_along(k)) {
                rr <- row[k[j]] - cv$row0 + seq(-half, half)
                cc <- col[k[j]] - cv$col0 + seq(-half, half)
                tiles[j, , , ] <- cv$rgb[rr, cc, , drop = FALSE]
                masks[j, , ] <- cv$lab[rr, cc] == k[j]
            }
            f <- provider$fun(tiles, masks, um_per_px, polys[][k])
            if (is.null(out)) {
                out <- matrix(NA_real_, length(ids), ncol(f),
                    dimnames = list(ids, colnames(f) %null%
                        paste0("f", seq_len(ncol(f)))))
            }
            out[k, ] <- f
        }
    }
    out
}

# The image's raster and the polygons in the raster's coordinates, and how
# many polygon units one raster unit is.
.img_frame <- function(polys, img) {
    if (inherits(img, "giottoAffineImage")) {
        pr <- affine(polys, img@affine, inv = TRUE)[]
        scale <- sqrt(sum(terra::expanse(.img_planar(polys[]),
            transform = FALSE)) / sum(terra::expanse(.img_planar(pr),
            transform = FALSE)))
        return(list(raster = img@raster_object, polys = pr,
            um_per_raster_unit = scale))
    }
    if (inherits(img, "giottoLargeImage")) {
        return(list(raster = img@raster_object, polys = polys[],
            um_per_raster_unit = 1))
    }
    stop("[calculateImageFeatures] the image must be a giottoLargeImage or ",
        "a giottoAffineImage.", call. = FALSE)
}

# Polygon coordinates are planar (microns or pixels): say so, so terra
# measures them as such and does not warn about a missing CRS.
.img_planar <- function(sv) {
    if (!nzchar(terra::crs(sv))) terra::crs(sv) <- "local"
    sv
}

# Pixels rows[1]:rows[2] x cols[1]:cols[2] of the image (white outside it)
# and the polygons rasterised on the same grid (label = .idx, 0 = none).
.img_canvas <- function(r, pr, rows, cols) {
    nr <- rows[2L] - rows[1L] + 1L
    nc <- cols[2L] - cols[1L] + 1L
    rgb <- array(255L, c(nr, nc, 3L))
    lab <- matrix(0L, nr, nc)
    ri <- max(rows[1L], 1L):min(rows[2L], terra::nrow(r))
    ci <- max(cols[1L], 1L):min(cols[2L], terra::ncol(r))
    if (length(ri) && length(ci) && rows[2L] >= 1L && cols[2L] >= 1L &&
            rows[1L] <= terra::nrow(r) && cols[1L] <= terra::ncol(r)) {
        e <- terra::ext(r, cells = terra::cellFromRowCol(r,
            c(ri[1L], ri[length(ri)]), c(ci[1L], ci[length(ci)])))
        sub <- terra::crop(r, e)
        a <- terra::as.array(sub)
        a[is.na(a)] <- 255
        a <- switch(as.character(dim(a)[3L]),
            "1" = a[, , c(1L, 1L, 1L), drop = FALSE],
            "3" = a,
            a[, , 1:3, drop = FALSE])
        rgb[ri - rows[1L] + 1L, ci - cols[1L] + 1L, ] <- as.integer(round(a))
        hit <- terra::relate(pr, terra::as.polygons(e), "intersects")[, 1L]
        if (any(hit)) {
            lr <- terra::rasterize(pr[hit], sub[[1L]], field = ".idx",
                background = 0)
            lab[ri - rows[1L] + 1L, ci - cols[1L] + 1L] <-
                terra::as.matrix(lr, wide = TRUE)
        }
    }
    list(rgb = rgb, lab = lab, row0 = rows[1L] - 1L, col0 = cols[1L] - 1L)
}

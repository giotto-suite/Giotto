# CosMx RNA + protein (multiomic) exports: the in-memory reader's protein
# option, FOV filtering while reading, and the pass-through to GiottoDisk.
#
# Fixture: a two-FOV AtoMx-like export written by the test. RNA and protein
# folders use the same file names, as in the vendor export.
#   cells     c_1_1_1..4 (fov 1), c_1_2_1..3 (fov 2)
#   RNA       every cell but c_1_2_3, plus a background row (cell_ID 0)
#   protein   every cell, rows in reverse order; decimal values, a Channel-*
#             column between the targets, two isotypes, a target with a space
#   polygons  every cell but c_1_1_4, ordered by cell number then fov (as
#             AtoMx writes them); 1,200 vertex rows, and the vertex rows
#             after row 1,000 carry `1e+05` (UPSTREAM_ISSUES #27)
# so keep_cells = "intersect" keeps 5 cells and drops c_1_2_3 (protein only)
# and c_1_1_4 (no polygon).

.cosmx_mo_fixture <- function() {
    root <- tempfile("cosmx_mo_")
    rna_dir <- file.path(root, "rna")
    prot_dir <- file.path(root, "protein")
    dir.create(rna_dir, recursive = TRUE)
    dir.create(prot_dir)
    dt <- data.table::data.table
    fw <- function(x, dir, name, ...) {
        data.table::fwrite(x, file.path(dir, name), ...)
    }

    cells <- dt(fov = c(1L, 1L, 1L, 1L, 2L, 2L, 2L),
                cell_ID = c(1:4, 1:3))
    cells[, cell := sprintf("c_1_%d_%d", fov, cell_ID)]
    rna_cells <- cells[cell != "c_1_2_3"]
    poly_cells <- cells[cell != "c_1_1_4"]

    rna <- dt(
        fov = c(1L, rna_cells$fov), cell_ID = c(0L, rna_cells$cell_ID),
        CD68 = c(9L, 0L, 3L, 0L, 0L, 5L, 0L),
        GeneA = c(4L, 2L, 0L, 0L, 1L, 0L, 0L),
        GeneB = c(0L, 0L, 0L, 7L, 0L, 0L, 1L),
        Negative1 = c(1L, 0L, 0L, 1L, 0L, 0L, 0L),
        SystemControl1 = c(0L, 0L, 1L, 0L, 0L, 0L, 0L)
    )
    fw(rna, rna_dir, "Run_exprMat_file.csv.gz")

    prot <- dt(
        fov = cells$fov, cell_ID = cells$cell_ID,
        CD68 = c(12.5, 230.04, 7, 15.25, 8.5, 99.99, 5112.92),
        `Channel-DNA` = c(500.1, 610.2, 720.3, 830.4, 940.5, 150.6, 260.7),
        EpCAM = c(3.33, 4.44, 1200.75, 6.66, 7.77, 8.88, 9.99),
        `NF-kB p65` = c(101.01, 102.02, 103.03, 104.04, 105, 106.06, 107.07),
        `Ms IgG1` = c(6.1, 7.2, 8.3, 9.4, 5.5, 6.6, 7.7),
        `Rb IgG` = c(5.9, 8.1, 7.2, 6.3, 9.4, 4.5, 3.6),
        check.names = FALSE
    )
    fw(prot[rev(seq_len(nrow(prot)))], prot_dir, "Run_exprMat_file.csv.gz")
    fw(dt(DisplayName = c("CD68", "EpCAM", "NF-kB p65", "Ms IgG1", "Rb IgG"),
          Barcode = NA,
          CodeClass = c("Endogenous", "Endogenous", "Endogenous",
                        "Negative", "Negative"),
          ProbeID = sprintf("CPROT%05d", 1:5)),
       prot_dir, "plex-test.txt", sep = "\t")

    meta <- function(x, protein) {
        m <- dt(fov = x$fov, Area = 100L + seq_len(nrow(x)),
                CenterX_global_px = 1000 + 200 * seq_len(nrow(x)),
                CenterY_global_px = 99800, cell_ID = x$cell_ID,
                Panel = if (protein) "64-plex" else "WTX",
                assay_type = if (protein) "protein" else "RNA")
        if (protein) {
            p <- prot[match(x$cell, cells$cell)]
            m[, nCount_RNA := p$CD68 + p$`Channel-DNA` + p$EpCAM +
                  p$`NF-kB p65`]
            m[, nFeature_RNA := 4L]
            m[, nCount_negprobes := p$`Ms IgG1` + p$`Rb IgG`]
        } else {
            r <- rna[cell_ID != 0L][match(x$cell, rna_cells$cell)]
            m[, nCount_RNA := r$CD68 + r$GeneA + r$GeneB]
            m[, nFeature_RNA := 3L]
            m[, nCount_negprobes := r$Negative1]
        }
        m[, Area.um2 := Area * 0.12028^2]
        m[, cell := x$cell]
        m
    }
    fw(meta(rna_cells, FALSE), rna_dir, "Run_metadata_file.csv.gz")
    fw(meta(cells, TRUE), prot_dir, "Run_metadata_file.csv.gz")

    fovpos <- dt(FOV = 1:2, x_global_px = c(1000, 1600),
                 y_global_px = c(100100, 100100),
                 x_global_mm = c(0.12, 0.19), y_global_mm = c(12, 12))
    fw(fovpos, rna_dir, "Run_fov_positions_file.csv.gz")
    fw(fovpos, prot_dir, "Run_fov_positions_file.csv.gz")

    # 200-vertex circles; the last cell's top vertex is y = 100000 exactly,
    # written by fwrite as `1e+05` (as in the vendor export)
    ang <- seq(0, 2 * pi, length.out = 201L)[-201L]
    polys <- data.table::rbindlist(lapply(seq_len(nrow(poly_cells)), function(k) {
        cy <- if (k == nrow(poly_cells)) 99940 else 99800
        dt(fov = poly_cells$fov[k], cellID = poly_cells$cell_ID[k],
           cell = poly_cells$cell[k],
           x_local_px = round(200 * k + 60 * cos(ang)),
           y_local_px = round(100 + 60 * sin(ang)),
           x_global_px = round(1000 + 200 * k + 60 * cos(ang)),
           y_global_px = round(cy + 60 * sin(ang)))
    }))
    polys[, y_global_px := as.double(y_global_px)]
    # vendor order: by cell number, then fov
    polys <- polys[order(cellID, fov)]
    fw(polys, rna_dir, "Run-polygons.csv.gz")
    fw(polys, prot_dir, "Run-polygons.csv.gz")

    list(root = root, rna_dir = rna_dir, prot_dir = prot_dir, cells = cells,
         rna = rna, prot = prot, polys = polys,
         keep = c("c_1_1_1", "c_1_1_2", "c_1_1_3", "c_1_2_1", "c_1_2_2"))
}

# value and every warning of `expr`
.with_warnings <- function(expr) {
    w <- character()
    val <- withCallingHandlers(expr, warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
    })
    list(value = val, warnings = w)
}

.mo_reader <- function(fx, ...) {
    suppressWarnings(importCosMx(fx$rna_dir, protein_dir = fx$prot_dir, ...))
}

.mo_gobject <- function(fx, ...) {
    .with_warnings(createGiottoCosMxObject(
        fx$rna_dir, protein_dir = fx$prot_dir,
        feat_type = c("rna", "negprobes", "falsecode"),
        split_keyword = list("^Negative", "^SystemControl"),
        load_transcripts = FALSE, load_images = NULL, cores = 1L,
        verbose = FALSE, ...
    ))
}

# expected protein matrix (features x cells) from the fixture table
.prot_mat <- function(fx, feats, ids) {
    p <- fx$prot[match(ids, fx$cells$cell)]
    m <- t(as.matrix(p[, feats, with = FALSE]))
    colnames(m) <- ids
    m
}

fx <- .cosmx_mo_fixture()


test_that("load_protein() splits targets and isotypes and drops Channel-*", {
    x <- .mo_reader(fx)
    expect_identical(x$protein_dir, fx$prot_dir)
    expect_true(is.function(x$load_protein))
    expect_output(show(x), "protein +:")

    pl <- x$load_protein()
    expect_identical(vapply(pl, featType, character(1L)),
                     c("protein", "protein_neg"))
    expect_identical(rownames(pl[[1L]][]), c("CD68", "EpCAM", "NF-kB p65"))
    expect_identical(rownames(pl[[2L]][]), c("Ms IgG1", "Rb IgG"))
    expect_false(any(grepl("^Channel-", unlist(lapply(pl, function(e)
        rownames(e[]))))))

    ids <- colnames(pl[[1L]][])
    expect_setequal(ids, fx$cells$cell)
    expect_identical(as.matrix(pl[[1L]][]),
        .prot_mat(fx, c("CD68", "EpCAM", "NF-kB p65"), ids))
    expect_identical(as.matrix(pl[[2L]][]),
        .prot_mat(fx, c("Ms IgG1", "Rb IgG"), ids))
})

test_that("a protein matrix that will not split as expected warns", {
    x <- .mo_reader(fx)
    w <- .with_warnings(x$load_protein(split_keyword = list("^NoSuch$")))$warnings
    expect_true(any(grepl("no isotype control", w)))
    w <- .with_warnings(x$load_protein(
        path = file.path(fx$rna_dir, "Run_exprMat_file.csv.gz")))$warnings
    expect_true(any(grepl("is `protein_dir` the RNA export", w)))
    expect_length(.with_warnings(x$load_protein())$warnings, 0L)
})

test_that("protein_dir can be cleared with NULL", {
    x <- .mo_reader(fx)
    x$protein_dir <- NULL
    expect_length(x$protein_dir, 0L)
    expect_error(x$load_protein(), "no protein exprMat_file found")
})

test_that("FOV rows of a csv are read in file order, batch by batch", {
    f <- tempfile(fileext = ".csv.gz")
    dt <- data.table::data.table(fov = rep(c(1L, 2L, 3L), each = 5L),
        cellID = 1:15, x_global_px = c(1e5, seq_len(14) + 0.5),
        cell = letters[1:15])
    data.table::fwrite(dt, f)
    cls <- c("integer", "integer", "double", "character")
    got <- .cosmx_read_fovs(f, fovs = c(3L, 1L), drop = "cell",
        col_classes = cls, cores = 1L, batch_lines = 4L)
    expect_identical(got, dt[fov %in% c(1L, 3L), !"cell"])
    expect_identical(nrow(.cosmx_read_fovs(f, fovs = 9L, col_classes = cls,
        cores = 1L)), 0L)
})

test_that("a non-numeric value in the exprMat is a clear error", {
    f <- tempfile(fileext = ".csv.gz")
    con <- gzfile(f, "w")
    writeLines(c("fov,cell_ID,GeneA,GeneB", "1,1,3,0", "1,2,x,1"), con)
    close(con)
    expect_error(.cosmx_read_exprmat(f, c("fov", "cell_ID", "GeneA", "GeneB")),
        "non-numeric values in column GeneA")
})

test_that("a cell with NA values is left out as missing from that matrix", {
    f <- tempfile(fileext = ".csv.gz")
    con <- gzfile(f, "w")
    writeLines(c("fov,cell_ID,GeneA,GeneB", "1,1,3,0", "1,2,NA,1", "1,3,,2"),
        con)
    close(con)
    res <- .with_warnings(
        .cosmx_read_exprmat(f, c("fov", "cell_ID", "GeneA", "GeneB")))
    expect_identical(colnames(res$value), "c_1_1_1")
    expect_match(gsub("\\s+", " ", res$warnings), paste("2 cells in .* have",
        "NA values and are left out, as cells missing from it \\(e.g.",
        "c_1_1_2, c_1_1_3\\)"))

    root <- tempfile("cosmx_na_")
    dir.create(root)
    file.copy(c(fx$rna_dir, fx$prot_dir), root, recursive = TRUE)
    fx2 <- fx
    fx2$rna_dir <- file.path(root, "rna")
    fx2$prot_dir <- file.path(root, "protein")
    p <- data.table::copy(fx$prot)
    p$EpCAM[match("c_1_1_2", fx$cells$cell)] <- NA
    data.table::fwrite(p[rev(seq_len(nrow(p)))],
        file.path(fx2$prot_dir, "Run_exprMat_file.csv.gz"))
    res <- .mo_gobject(fx2, load_expression = TRUE)
    w <- gsub("\\s+", " ", res$warnings)
    expect_true(any(grepl("1 cell in Run_exprMat_file.csv.gz has NA values", w)))
    expect_true(any(grepl(
        "1 with RNA but no protein expression \\(e.g. c_1_1_2\\)", w)))
    expect_identical(spatIDs(res$value), setdiff(fx$keep, "c_1_1_2"))
    expect_identical(
        as.matrix(getExpression(res$value, feat_type = "protein",
            output = "matrix")),
        .prot_mat(fx, c("CD68", "EpCAM", "NF-kB p65"),
            setdiff(fx$keep, "c_1_1_2")))
})

test_that("nothing to load is a clear error", {
    expect_error(createGiottoCosMxObject(fx$rna_dir, load_expression = FALSE,
        load_transcripts = FALSE, load_images = NULL, cores = 1L, verbose = FALSE),
        "nothing to build the object")
})

test_that("an empty split group is dropped by name, not kept empty", {
    x <- .mo_reader(fx)
    pl <- x$load_protein(feat_type = c("protein", "protein_neg", "protein_morph"),
        split_keyword = list("^NoSuch", "^Channel-"), drop_cols = NULL)
    expect_identical(vapply(pl, featType, character(1L)),
        c("protein", "protein_morph"))
    expect_true(all(grepl("^Channel-", rownames(pl[[2L]][]))))
    expect_error(x$load_protein(feat_type = c("protein", "protein_neg", "x")),
        "3 feat_type\\(s\\) given for 2")
})

test_that("isotypes come from the plex file, else from names with IgG", {
    plex <- list.files(fx$prot_dir, pattern = "^plex-", full.names = TRUE)
    expr <- file.path(fx$prot_dir, "Run_exprMat_file.csv.gz")
    kw <- .cosmx_protein_negatives(plex, expr)
    expect_identical(kw, list("^(Ms IgG1|Rb IgG)$"))
    expect_identical(.cosmx_protein_negatives(NULL, expr), list("IgG"))

    # names are regex-escaped
    d <- tempfile("plex_")
    dir.create(d)
    data.table::fwrite(data.table::data.table(
        DisplayName = c("A (ctrl)", "B.1"), CodeClass = "Negative"),
        file.path(d, "plex-x.txt"), sep = "\t")
    kw <- .cosmx_protein_negatives(file.path(d, "plex-x.txt"))
    expect_true(grepl(kw[[1L]], "A (ctrl)"))
    expect_false(grepl(kw[[1L]], "BB1"))

    # no plex file: same split from the "IgG" fallback
    d2 <- tempfile("prot_noplex_")
    dir.create(d2)
    file.copy(setdiff(list.files(fx$prot_dir, full.names = TRUE), plex), d2)
    x <- suppressWarnings(importCosMx(fx$rna_dir, protein_dir = d2))
    pl <- x$load_protein()
    expect_identical(rownames(pl[[2L]][]), c("Ms IgG1", "Rb IgG"))
    expect_identical(nrow(pl[[1L]][]), 3L)
})

test_that("one protein feat_type keeps the isotypes with the targets", {
    x <- .mo_reader(fx)
    pl <- x$load_protein(feat_type = "protein")
    expect_length(pl, 1L)
    expect_identical(featType(pl[[1L]]), "protein")
    expect_setequal(rownames(pl[[1L]][]),
        c("CD68", "EpCAM", "NF-kB p65", "Ms IgG1", "Rb IgG"))
    expect_error(x$load_protein(feat_type = "protein",
        split_keyword = list("IgG")), "1 feat_type\\(s\\) given for 2")
})

test_that("keep_cells = 'intersect' keeps one cell set, in one order", {
    res <- .mo_gobject(fx, load_expression = TRUE)
    g <- res$value
    w <- grep("not loaded", res$warnings, value = TRUE)
    expect_length(w, 1L)
    expect_match(gsub("\\s+", " ", w), paste("2 cells are not loaded: 1 with",
        "protein but no RNA expression \\(e.g. c_1_2_3\\); 1 with expression",
        "but no polygon \\(e.g. c_1_1_4\\)"))

    expect_identical(spatIDs(g), fx$keep)
    for (ft in c("rna", "negprobes", "falsecode", "protein", "protein_neg")) {
        m <- getExpression(g, feat_type = ft, output = "matrix")
        expect_identical(colnames(m), spatIDs(g), info = ft)
    }
    expect_identical(sort(spatIDs(g)),
        sort(getPolygonInfo(g, return_giottoPolygon = TRUE)$poly_ID))

    # protein file rows are reversed: values must follow the ids
    expect_identical(
        as.matrix(getExpression(g, feat_type = "protein", output = "matrix")),
        .prot_mat(fx, c("CD68", "EpCAM", "NF-kB p65"), fx$keep))
    rna <- getExpression(g, feat_type = "rna", output = "matrix")
    r <- fx$rna[cell_ID != 0L]
    r_ids <- sprintf("c_1_%d_%d", r$fov, r$cell_ID)
    expect_equal(unname(rna["CD68", ]),
                 as.numeric(r$CD68[match(fx$keep, r_ids)]))
    expect_identical(rownames(getExpression(g, feat_type = "negprobes",
        output = "matrix")), "Negative1")
})

test_that("protein metadata goes to feat_type protein under protein names", {
    g <- .mo_gobject(fx, load_expression = TRUE)$value
    rm_ <- pDataDT(g, feat_type = "rna")
    pm <- pDataDT(g, feat_type = "protein")
    expect_identical(rm_$cell_ID, spatIDs(g))
    expect_identical(pm$cell_ID, spatIDs(g))

    pmeta <- data.table::fread(file.path(fx$prot_dir,
        "Run_metadata_file.csv.gz"))
    rmeta <- data.table::fread(file.path(fx$rna_dir,
        "Run_metadata_file.csv.gz"))
    expect_identical(rm_$nCount_RNA, rmeta$nCount_RNA[match(fx$keep,
        rmeta$cell)])
    expect_identical(pm$nCount_protein, pmeta$nCount_RNA[match(fx$keep,
        pmeta$cell)])
    expect_identical(pm$nCount_protein_neg,
        pmeta$nCount_negprobes[match(fx$keep, pmeta$cell)])
    expect_true(all(c("nFeature_protein", "Panel_protein",
        "assay_type_protein", "Area.um2", "fov") %in% colnames(pm)))
    expect_false(any(c("nCount_RNA", "nFeature_RNA", "Panel", "Area")
        %in% colnames(pm)))
    expect_identical(pm$fov, rm_$fov)
    expect_identical(pm$fov, c(1L, 1L, 1L, 2L, 2L))

    expect_identical(
        .cosmx_protein_meta_names(c("nCount_RNA", "nCount_negprobes",
            "Panel")),
        c(nCount_RNA = "nCount_protein", nCount_negprobes =
            "nCount_protein_neg", Panel = "Panel_protein"))
})

test_that("protein without RNA expression keeps protein cells with polygons", {
    res <- .mo_gobject(fx, load_expression = FALSE)
    g <- res$value
    # cells in the order of the protein file (written in reverse)
    expect_identical(spatIDs(g), rev(setdiff(fx$cells$cell, "c_1_1_4")))
    expect_match(grep("not loaded", res$warnings, value = TRUE),
                 "1 with expression but no polygon")
    expect_identical(list_expression(g)$feat_type,
                     c("protein", "protein_neg"))
})

test_that("without protein, RNA cells are matched to the polygons", {
    ra <- .with_warnings(createGiottoCosMxObject(fx$rna_dir,
        load_expression = TRUE, load_transcripts = FALSE, load_images = NULL,
        cores = 1L, verbose = FALSE))
    a <- ra$value
    b <- .mo_gobject(fx, load_expression = TRUE, load_protein = FALSE)$value
    expect_false(any(c("protein", "protein_neg") %in%
        list_expression(a)$feat_type))
    expect_false(any(c("protein", "protein_neg") %in%
        list_expression(b)$feat_type))
    ids <- intersect(fx$cells$cell[fx$cells$cell != "c_1_2_3"],
        fx$cells$cell[fx$cells$cell != "c_1_1_4"])
    expect_setequal(spatIDs(a), ids)
    expect_identical(colnames(getExpression(a, output = "matrix")), spatIDs(a))
    expect_identical(pDataDT(a)$cell_ID, spatIDs(a))
    expect_true(any(grepl("1 with expression but no polygon", ra$warnings)))
    expect_identical(spatIDs(a), spatIDs(b))
})

test_that("meta_cols is ignored with a warning without a backend", {
    w <- .with_warnings(createGiottoCosMxObject(fx$rna_dir,
        meta_cols = "fov", load_expression = TRUE, load_transcripts = FALSE,
        load_images = NULL, cores = 1L, verbose = FALSE))$warnings
    expect_true(any(grepl("`meta_cols` is ignored", w)))
})

test_that("the exprMat reader filters FOVs while reading, as the old path", {
    # the old path: whole-file fread, filter, ids, dt_to_matrix
    old_read <- function(path, fovs = NULL, drop = NULL) {
        d <- data.table::fread(path, drop = drop)
        if (!is.null(fovs)) d <- d[fov %in% fovs]
        d <- d[cell_ID != 0L]
        d[, cell_ID := sprintf("c_%d_%d_%d", 1, fov, cell_ID)]
        d[, fov := NULL]
        t_flex(dt_to_matrix(d))
    }
    rpath <- file.path(fx$rna_dir, "Run_exprMat_file.csv.gz")
    ppath <- file.path(fx$prot_dir, "Run_exprMat_file.csv.gz")
    rhdr <- .cosmx_csv_header(rpath)
    phdr <- .cosmx_csv_header(ppath)
    for (fv in list(NULL, 1L, 2L, c(1L, 2L))) {
        for (bv in c(2e7, 3 * length(rhdr))) {
            new <- .cosmx_read_exprmat(rpath, header = rhdr, fovs = fv,
                batch_values = bv, cores = 1L)
            expect_identical(new, old_read(rpath, fv))
        }
        new <- .cosmx_read_exprmat(ppath, header = phdr, fovs = fv,
            drop = "Channel-DNA", batch_values = 2 * length(phdr), cores = 1L)
        expect_identical(new, old_read(ppath, fv, drop = "Channel-DNA"))
    }
    expect_s4_class(.cosmx_read_exprmat(rpath, rhdr), "dgCMatrix")
    expect_s4_class(.cosmx_read_exprmat(ppath, phdr), "dgeMatrix")
    expect_false(any(grepl("_0$", colnames(
        .cosmx_read_exprmat(rpath, rhdr)))))
    expect_error(.cosmx_read_exprmat(rpath, rhdr, fovs = 9L), "no cells found")
})

test_that("FOVs are applied to expression, protein, metadata and polygons", {
    x <- .mo_reader(fx, fovs = 2)
    ex <- x$load_expression()
    expect_identical(colnames(ex[[1L]][]), c("c_1_2_1", "c_1_2_2"))
    pl <- x$load_protein()
    expect_setequal(colnames(pl[[1L]][]), c("c_1_2_1", "c_1_2_2", "c_1_2_3"))
    expect_identical(x$load_cellmeta()[]$cell_ID, c("c_1_2_1", "c_1_2_2"))

    res <- .mo_gobject(fx, load_expression = TRUE, FOVs = 1)
    expect_identical(spatIDs(res$value), c("c_1_1_1", "c_1_1_2", "c_1_1_3"))
    expect_match(grep("not loaded", res$warnings, value = TRUE),
                 "1 with expression but no polygon")
})

test_that("csv polygons with 1e+05 coordinates read with FOVs (#27)", {
    ppath <- file.path(fx$rna_dir, "Run-polygons.csv.gz")
    expect_true(any(grepl("1e+05", readLines(ppath), fixed = TRUE)))
    # the fixture reproduces the issue: int32 inferred from the first rows
    expect_error(GiottoUtils::read_colmatch(ppath, col = "fov",
        values_to_match = 2L), "1e\\+05")

    x <- suppressWarnings(importCosMx(fx$rna_dir, fovs = 2,
        poly_pref = "csv"))
    p <- x$load_polys()
    expect_identical(sort(p$poly_ID), c("c_1_2_1", "c_1_2_2", "c_1_2_3"))
    expect_identical(as.vector(terra::ext(p[]))[["ymax"]], 1e5)

    p <- x$load_polys(keep_ids = c("c_1_2_1", "c_1_2_3"))
    expect_identical(sort(p$poly_ID), c("c_1_2_1", "c_1_2_3"))
})

test_that("a split that matches one feature stays a matrix", {
    x <- suppressWarnings(importCosMx(fx$rna_dir))
    ex <- x$load_expression(feat_type = c("rna", "negprobes"),
        split_keyword = list("^Negative"))
    expect_identical(dim(ex[[2L]][]), c(1L, 6L))
    expect_identical(rownames(ex[[2L]][]), "Negative1")
})

test_that("the reader's protein_dir can be set after creation", {
    x <- suppressWarnings(importCosMx(fx$rna_dir))
    expect_null(x@paths$protein_expr_path)
    expect_error(x$load_protein(), "protein exprMat_file")
    x <- suppressWarnings({
        x$protein_dir <- fx$prot_dir
        x
    })
    expect_identical(basename(x@paths$protein_expr_path),
                     "Run_exprMat_file.csv.gz")
    expect_identical(nrow(x$load_protein()[[1L]][]), 3L)
})


# GiottoDisk pass-through ####
# importCosMx(backend =) returns GiottoDisk::importCosMxDisk(), whose reader's
# create_gobject() takes the same arguments. A mock records what reaches it.

.mock_disk <- function(env, protein = TRUE) {
    fake_reader <- list(create_gobject = function(...) {
        env$create <- list(...)
        "disk gobject"
    })
    f <- function(cosmx_dir = NULL, backend, slide = 1, fovs = NULL,
                  version = "default", micron = FALSE, px2um = 0.12028,
                  poly_pref = c("csv", "mask"), protein_dir = NULL) {
        env$import_vals <- as.list(environment())
        env$import <- as.list(match.call())[-1L]
        fake_reader
    }
    if (!protein) formals(f)$protein_dir <- NULL
    f
}

test_that("protein arguments reach GiottoDisk with a backend", {
    skip_if_not_installed("GiottoDisk")
    rec <- new.env()
    local_mocked_bindings(importCosMxDisk = .mock_disk(rec),
                          .package = "GiottoDisk")
    out <- createGiottoCosMxObject(fx$rna_dir, backend = "a_backend",
        protein_dir = fx$prot_dir, load_expression = TRUE,
        load_transcripts = FALSE, load_images = NULL, FOVs = 2,
        protein_feat_type = c("adt", "adt_neg"),
        protein_metadata_path = "meta.csv", meta_cols = c("fov", "Area"))
    expect_identical(out, "disk gobject")
    expect_identical(rec$create$meta_cols, c("fov", "Area"))
    expect_identical(rec$import_vals$protein_dir, fx$prot_dir)
    expect_identical(rec$import_vals$backend, "a_backend")
    expect_identical(rec$import_vals$fovs, 2)
    expect_identical(rec$create$protein_feat_type, c("adt", "adt_neg"))
    expect_identical(rec$create$keep_cells, "intersect")
    expect_identical(rec$create$protein_metadata_path, "meta.csv")
    expect_false("load_protein" %in% names(rec$create))
    expect_false("protein_expression_path" %in% names(rec$create))

    createGiottoCosMxObject(fx$rna_dir, backend = "a_backend",
        protein_dir = fx$prot_dir, load_protein = FALSE,
        protein_expression_path = "expr.csv")
    expect_identical(rec$create$load_protein, FALSE)
    expect_identical(rec$create$protein_expression_path, "expr.csv")

    r <- importCosMx(fx$rna_dir, backend = "a_backend",
        protein_dir = fx$prot_dir)
    expect_identical(rec$import_vals$protein_dir, fx$prot_dir)
})

test_that("RNA-only calls pass GiottoDisk the same arguments as before", {
    skip_if_not_installed("GiottoDisk")
    rec <- new.env()
    local_mocked_bindings(importCosMxDisk = .mock_disk(rec, protein = FALSE),
                          .package = "GiottoDisk")
    createGiottoCosMxObject(fx$rna_dir, backend = "a_backend",
        load_expression = TRUE)
    expect_identical(names(rec$import), c("cosmx_dir", "backend", "slide",
        "fovs", "version", "micron", "px2um", "poly_pref"))
    expect_identical(sort(names(rec$create)), sort(c("feat_type",
        "split_keyword", "load_images", "image_negative_y",
        "load_expression", "load_cellmeta", "load_transcripts",
        "instructions", "cores", "verbose")))

    # a GiottoDisk without the protein option: a clear error
    err <- tryCatch(createGiottoCosMxObject(fx$rna_dir, backend = "a_backend",
        protein_dir = fx$prot_dir), error = conditionMessage)
    expect_match(gsub("\\s+", " ", err), "cannot read CosMx protein data")
})

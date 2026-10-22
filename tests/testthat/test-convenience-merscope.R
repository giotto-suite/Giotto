# synthetic MERSCOPE region directory ####

.mk_merscope_dir <- function() {
    d <- tempfile("merscope_")
    dir.create(d)
    n_cells <- 6L
    sq <- function(cx, cy, r = 2) {
        terra::geom(terra::vect(sprintf(
            "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
            cx - r, cy - r, cx + r, cy - r, cx + r, cy + r,
            cx - r, cy + r, cx - r, cy - r)), wkb = TRUE)[[1L]]
    }
    cx <- 10 * seq_len(n_cells)
    cy <- rep(c(10, 30), length.out = n_cells)
    tb <- arrow::arrow_table(
        EntityID = seq_len(n_cells),
        ZIndex = rep(0L, n_cells),
        Type = rep("cell", n_cells)
    )
    tb$Geometry <- arrow::Array$create(
        lapply(seq_len(n_cells), function(i) sq(cx[i], cy[i])),
        type = arrow::binary()
    )
    arrow::write_parquet(tb, file.path(d, "cell_boundaries.parquet"))

    genes <- c("A", "B", "Blank-1")
    tx <- data.table::data.table(
        V1 = 0:(3L * n_cells - 1L),
        global_x = rep(cx, each = 3L),
        global_y = rep(cy, each = 3L),
        global_z = 0,
        fov = 0L,
        gene = rep(genes, times = n_cells),
        transcript_id = "t",
        cell_id = rep(seq_len(n_cells), each = 3L)
    )
    data.table::setnames(tx, "V1", "")
    data.table::fwrite(tx, file.path(d, "detected_transcripts.csv"))

    cbg <- data.table::data.table(
        cell = seq_len(n_cells), A = 1:n_cells, B = 2L, `Blank-1` = 1L
    )
    data.table::fwrite(cbg, file.path(d, "cell_by_gene.csv"))
    d
}


# importMerscope ####

test_that("importMerscope returns a MerscopeReader", {
    skip_if_not_installed("arrow")
    r <- importMerscope(.mk_merscope_dir())
    expect_s4_class(r, "MerscopeReader")
})


# createGiottoMerscopeObject ####

test_that("expression is read from cell_by_gene.csv", {
    skip_if_not_installed("arrow")
    d <- .mk_merscope_dir()
    g <- createGiottoMerscopeObject(d, verbose = FALSE)
    m <- as.matrix(getExpression(g, values = "raw", feat_type = "rna",
                                 output = "matrix"))
    expect_equal(unname(m["A", as.character(1:6)]), as.numeric(1:6))
    expect_equal(unname(m["B", as.character(1:6)]), rep(2, 6))
    expect_true("neg_probe" %in% list_expression(g)$feat_type)
})

test_that("y coordinates are negated", {
    skip_if_not_installed("arrow")
    g <- createGiottoMerscopeObject(.mk_merscope_dir(), verbose = FALSE)
    e <- as.vector(ext(getPolygonInfo(g, return_giottoPolygon = TRUE)))
    expect_lt(e[["ymax"]], 0)
})

test_that("calculate_overlap is deprecated", {
    skip_if_not_installed("arrow")
    lifecycle::expect_deprecated(
        createGiottoMerscopeObject(.mk_merscope_dir(),
                                   calculate_overlap = TRUE, verbose = FALSE)
    )
})

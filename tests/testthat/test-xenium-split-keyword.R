# Atera WTA panels name their genomic control probes `Intergenic_Region_*`
# instead of `GenomicControl*`. The default split_keyword has to catch both, or
# the probes load as `rna` and aggregation adds them to the gene expression.

.mk_tx_parquet <- function(dir) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    feats <- c("GENE_A", "GENE_B", "NegControlProbe_00002",
        "UnassignedCodeword_0001", "NegControlCodeword_0500",
        "GenomicControl_0001", "Intergenic_Region_100",
        "DeprecatedCodeword_11440")
    dt <- data.table::data.table(
        transcript_id = seq_along(feats),
        cell_id = "a",
        feature_name = feats,
        x_location = seq_along(feats),
        y_location = seq_along(feats),
        qv = 40
    )
    p <- file.path(dir, "transcripts.parquet")
    arrow::write_parquet(dt, p)
    p
}

test_that("default split puts Intergenic_Region probes in GenomicControl", {
    skip_if_not_installed("arrow")
    d <- file.path(tempdir(), "xen_split_1")
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    gp <- Giotto:::.xenium_transcript(.mk_tx_parquet(d), verbose = FALSE)
    names(gp) <- vapply(gp, objName, character(1L))

    expect_setequal(featIDs(gp$rna), c("GENE_A", "GENE_B"))
    expect_setequal(featIDs(gp$GenomicControl),
        c("GenomicControl_0001", "Intergenic_Region_100"))
    expect_identical(featIDs(gp$DeprecatedCodeword), "DeprecatedCodeword_11440")
})

test_that("the 10x reader defaults agree on feat_type and split_keyword", {
    # guards the copies of the default against drifting apart
    ref <- formals(Giotto:::.xenium_transcript)
    for (f in list(createGiottoXeniumObject, createGiottoAteraObject)) {
        expect_identical(formals(f)$feat_type, ref$feat_type)
        expect_identical(formals(f)$split_keyword, ref$split_keyword)
    }
    # one feat_type per keyword set, plus the default `rna`
    expect_length(eval(ref$feat_type), length(eval(ref$split_keyword)) + 1L)
})

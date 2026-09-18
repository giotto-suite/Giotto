# Documentation consistency for the two signature-based families.
#
# These assert the *property*, derived from the installed help database, rather
# than a fixed list of topics -- so a function added to either family later is
# held to the same rules without anyone remembering to update a test.
#
# What went wrong before these existed: all nine titles were the function name
# repeated, `@param name` was missing its first word in four places so
# `?runRankEnrich` rendered "name: to give to spatial enrichment results",
# `runSpatialEnrich` documented a default it does not have, and thirteen shared
# parameters carried between two and five descriptions each.

# Works both ways: against the installed package's Rd database when there is
# one (R CMD check), and against the source tree's man/ when there is not
# (devtools::load_all + test_file, where nothing has been installed).
.rd_db <- function() {
    db <- tryCatch(tools::Rd_db("Giotto"), error = function(e) list())
    if (length(db) > 0L) return(db)
    dir <- NULL
    for (cand in c("../../man", "../man", "man")) {
        if (dir.exists(cand)) { dir <- cand; break }
    }
    skip_if(is.null(dir), "no Rd database and no man/ directory")
    files <- list.files(dir, pattern = "\\.Rd$", full.names = TRUE)
    skip_if(length(files) == 0L, "man/ is empty")
    db <- lapply(files, tools::parse_Rd)
    names(db) <- basename(files)
    db
}

.rd_tag <- function(rd, tag) {
    out <- character()
    walk <- function(x) {
        if (identical(attr(x, "Rd_tag"), tag)) {
            out <<- c(out, paste(rapply(x, as.character, how = "unlist"),
                                 collapse = ""))
        }
        if (is.list(x)) lapply(x, walk)
    }
    walk(rd)
    trimws(out)
}

# \item{name}{body} pairs out of the \arguments section
.rd_params <- function(rd) {
    args <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\arguments"), rd)
    if (!length(args)) return(list())
    out <- list()
    for (a in args[[1]]) {
        if (!identical(attr(a, "Rd_tag"), "\\item") || length(a) < 2L) next
        nm <- trimws(paste(rapply(a[[1]], as.character, how = "unlist"),
                           collapse = ""))
        bd <- trimws(paste(rapply(a[[2]], as.character, how = "unlist"),
                           collapse = ""))
        bd <- gsub("\\s+", " ", bd)
        for (n in strsplit(nm, ",\\s*")[[1]]) out[[n]] <- bd
    }
    out
}

.family_members <- function(db, family) {
    keep <- vapply(db, function(rd) {
        any(grepl(sprintf("Other %s", family), .rd_tag(rd, "\\seealso"),
                  fixed = TRUE))
    }, logical(1L))
    db[keep]
}

FAMILIES <- c("feature set enrichment", "spatial deconvolution")


test_that("both families have members, found through their family tag", {
    db <- .rd_db()
    for (f in FAMILIES) {
        m <- .family_members(db, f)
        expect_gt(length(m), 3L)
    }
    # and they do not overlap -- deconvolution is a separate family on purpose,
    # because it returns proportions rather than scores
    a <- names(.family_members(db, FAMILIES[[1]]))
    b <- names(.family_members(db, FAMILIES[[2]]))
    expect_length(intersect(a, b), 0L)
})


test_that("no topic's title is just its own function name", {
    db <- .rd_db()
    for (f in FAMILIES) {
        for (nm in names(.family_members(db, f))) {
            ttl <- .rd_tag(.family_members(db, f)[[nm]], "\\title")
            fn <- sub("\\.Rd$", "", nm)
            expect_false(identical(ttl, fn),
                         info = sprintf("%s title is just the function name", nm))
            expect_gt(nchar(paste(ttl, collapse = "")), nchar(fn))
        }
    }
})


test_that("no parameter description starts mid-sentence", {
    # `@param name to give to spatial enrichment results` -- the word the
    # reader needs is the one that is missing.
    db <- .rd_db()
    for (f in FAMILIES) {
        for (nm in names(.family_members(db, f))) {
            for (p in names(.rd_params(.family_members(db, f)[[nm]]))) {
                body <- .rd_params(.family_members(db, f)[[nm]])[[p]]
                if (!nzchar(body)) next
                expect_false(grepl("^(to|of|for|with|in|and|or)\\b", body),
                    info = sprintf("%s: `%s` reads '%s'", nm, p,
                                   substr(body, 1, 50)))
            }
        }
    }
})


test_that("a shared parameter is described one way within each family", {
    db <- .rd_db()
    # Three kinds of legitimate difference, encoded as rules rather than as a
    # list of blessed topics, so a new offender is still caught:
    #
    #  1. a router tags an argument with the method it applies to -- "(rank)
    #     ..." -- which the single-method pages should not carry;
    #  2. one description *extends* another, e.g. makeSignMatrixRank's gobject
    #     adds a caveat the others do not need. A strict superset is a
    #     specialization, not a conflict;
    #  3. `param` names the param class the method dispatches on, so it is
    #     necessarily different on each; `...` carries the method-specific
    #     parameters on the enrichParam() factory but is unused on the
    #     analyzeData() methods, which document it only because it is in
    #     their usage; and runRankEnrich deprecates reverse_log_scale/logbase,
    #     so its text differs on purpose.
    per_topic <- c("param", "...", "reverse_log_scale", "logbase")
    for (f in FAMILIES) {
        members <- .family_members(db, f)
        seen <- list()
        for (nm in names(members)) {
            for (p in names(.rd_params(members[[nm]]))) {
                seen[[p]] <- union(seen[[p]], .rd_params(members[[nm]])[[p]])
            }
        }
        for (p in names(seen)) {
            v <- seen[[p]]
            if (length(v) == 1L) next
            if (p %in% per_topic) next
            # rule 1: strip a leading "(method) " tag
            base <- unique(sub("^\\([^)]+\\)\\s*", "", v))
            # rule 2: collapse variants that are extensions of a shorter one
            base <- base[!vapply(base, function(x) {
                any(vapply(setdiff(base, x), function(y) startsWith(x, y),
                           logical(1L)))
            }, logical(1L))]
            expect_true(length(unique(base)) == 1L,
                info = sprintf("%s in [%s] has %d unrelated descriptions:\n%s",
                               p, f, length(unique(base)),
                               paste("   -", v, collapse = "\n")))
        }
    }
})


test_that("markdown is processed, not printed", {
    # Only runSpatialEnrich had @md, so backticks and [fn()] links rendered
    # literally everywhere else -- including a lifecycle badge that came out as
    # raw inline-R source.
    db <- .rd_db()
    for (f in FAMILIES) {
        members <- .family_members(db, f)
        for (nm in names(members)) {
            txt <- paste(c(.rd_tag(members[[nm]], "\\description"),
                           unlist(.rd_params(members[[nm]]))), collapse = " ")
            expect_false(grepl("`", txt, fixed = TRUE),
                         info = sprintf("%s has a literal backtick", nm))
            expect_false(grepl("\\[[a-zA-Z_.]+\\(\\)\\]", txt),
                         info = sprintf("%s has an unrendered md link", nm))
        }
    }
})

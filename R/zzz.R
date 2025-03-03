# Run on library loading


.onAttach <- function(libname, pkgname) {
    ## print version number ##
    packageStartupMessage("Giotto Suite ", packageVersion("Giotto"))

    check_ver <- getOption("giotto.check_version", TRUE)
    if (isTRUE(check_ver)) {
        check_github_suite_ver()
        options("giotto.check_version" = FALSE)
    }


    ## cores detection ##
    check_core <- getOption("giotto.check_core", TRUE)
    if (isTRUE(check_core)) {
        cores <- determine_cores(cores = NA)
        data.table::setDTthreads(threads = cores)
        options("giotto.check_core" = FALSE)
    }

    # options #

    # GiottoVisuals #
    # ------------- #
    # colors continuous
    init_option("giotto.color_cd_pal", c("blue", "white", "red"))
    init_option("giotto.color_cs_pal", "viridis")
    init_option("giotto.color_c_rev", FALSE)

    # colors discrete
    init_option("giotto.color_d_pal", "distinct")
    init_option("giotto.color_d_rev", FALSE)
    init_option("giotto.color_d_strategy", "interpolate")

    # image resampling
    init_option("giotto.plot_img_max_sample", 5e5)
    init_option("giotto.plot_img_max_crop", 1e8)
    init_option("giotto.plot_img_max_resample_scale", 100)

    # GiottoUtils #
    # ----------- #
    init_option("giotto.verbose", TRUE)
}

.onLoad <- function(libname, pkgname) {
    all_matrix <- c("matrix", "Matrix")
    update_matrix_sig <- FALSE
    
    # suppress known warnings about onload class and method extensions
    suppressWarnings({
        # extensible classunions --------------------------------------------#
        if (requireNamespace("DelayedArray", quietly = TRUE)) {
            getClass("DelayedArray")
            all_matrix <- c(all_matrix, "DelayedArray")
            update_matrix_sig <- TRUE
        }
        if (requireNamespace("dbMatrix", quietly = TRUE)) {
            getClass("dbMatrix")
            all_matrix <- c(all_matrix, "dbMatrix")
            update_matrix_sig <- TRUE
        }
        # signature update
        if (isTRUE(update_matrix_sig)) {
            setClassUnion("allMatrix", members = all_matrix)
        }
        # methods extensions ------------------------------------------------#
        
        if (requireNamespace("dbMatrix", quietly = TRUE)) {
            setMethod("processData",
                signature(x = "dbMatrix", param = "logNormParam"),
                function(x, param) {
                    x[] <- dplyr::mutate(x[], x = x + param$offset)
                    log(x) / log(param$base)
                }
            )
        }
    })
}

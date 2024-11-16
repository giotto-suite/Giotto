#' @title doFeatureSetEnrichment
#' @name doFeatureSetEnrichment
#' @description Preform Gene Set Enrichment Analysis using marker genes
#' @param dryrun do a dry run, default TRUE.
#' @param path_to_GSEA path to GSEA command line executable, e.g. gsea-XXX.jar.
#' See details (1.) for more information.
#' @param GSEA_dataset path to a Human/Mouse collection from GSEA, e.g.
#' Hallmarks C1. See details (2.) for more information.
#' @param GSEA_ranked_file path to .rnk file for GSEA. See details (3.) for
#' more information
#' @param output_folder path to which the GSEA results will be saved. Default
#' is current working directory.
#' @param name_analysis_folder default output subdirectory prefix to which
#' results are saved. Will live within output_folder; equivalent of
#' "Analysis Name" in GSEA Application.
#' @param collapse only 'false' is supported. This will use your dataset as-is,
#' in the original format.
#' @param mode option selected in Advanced Field "Collapsing Mode for
#' Probe Sets => 1 gene"
#' @param norm normalization mode; only meandiv is supported.
#' @param nperm number of permutations, default 1000
#' @param scoring_scheme Default "weighted", equivalent of
#' "enrichment statistic" in GSEA Application
#' @param plot_top_x Default 20, number of enrichment plots to produce.
#' @param set_max default 500, equivalent to "max size; exclude larger sets"
#' in Basic Fields in GSEA Application
#' @param set_min default 15, equivalent to "min size; exclude smaller sets"
#' in Basic Fields in GSEA Application
#' @returns data.table
#' @details
#' NECESSARY PREREQUISITES
#' 1. download and install the COMMAND line (all platforms) gsea-XXX.jar
#' https://www.gsea-msigdb.org/gsea/downloads.jsp
#' 1.1. download zip file
#' 1.2. unzip and move to known location
#' (e.g. in path/to/your/applications/gsea/GSEA_4.3.2)
#'
#' 2. download the Human and Mouse collections
#' https://www.gsea-msigdb.org/gsea/msigdb/index.jsp or zipped folder
#' https://www.gsea-msigdb.org/gsea/downloads.jsp (all downloaded)
#'
#' 3. create ranked gene lists
#' format: data.table or data.frame with 2 columns
#' column 1 = gene symbols, column 2 = weight or rank metric
#' for more details, see:
#' https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29
#'
#' For more information on parameter conventions,
#' please reference GSEA's documentation here:
#' https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Syntax
#' @export
doFeatureSetEnrichment <- function(
        dryrun = TRUE,
        path_to_GSEA = NULL,
        GSEA_dataset = NULL,
        GSEA_ranked_file = NULL,
        output_folder = NULL,
        name_analysis_folder = "my_GSEA_analysis",
        collapse = "false",
        mode = c(
            "Abs_max_of_probes",
            "Max_probe",
            "Median_of_probes",
            "Mean_of_probes",
            "Sum_of_probes"
        ),
        norm = "meandiv",
        nperm = 1000,
        scoring_scheme = "weighted",
        plot_top_x = 20,
        set_max = 500,
        set_min = 15) {
    # set don't run to false as a start
    dont_run <- FALSE

    # SYSTEM CHECK FOR JAVA
    java_not_installed <- as.logical(system("java -version"))
    # returns 0 if java is installed (i.e., command runs successfully),
    # 1 otherwise
    if (java_not_installed) {
        stop(wrap_txt("Java must be installed for doFeatureSetEnrichment() to
                run. Please install Java: https://www.java.com/en/download/",
            errWidth = TRUE
        ))
    }


    mode <- match.arg(mode, choices = c(
        "Abs_max_of_probes",
        "Max_probe",
        "Median_of_probes",
        "Mean_of_probes",
        "Sum_of_probes"
    ))

    if (is.null(output_folder)) {
        output_folder <- paste0(
            getwd(), "/Feature_set_enrichment_results/"
        )
    }

    if (!dir.exists(output_folder)) {
        wrap_msg(paste0(
            "Directory does not yet exist. Creating directory at:",
            output_folder
        ))
        dir.create(output_folder)
    }

    # check for path to GSEA tool
    if (is.null(path_to_GSEA)) {
        stop("Path to the GSEA directory needs to be provided")
    }
    if (!file.exists(path_to_GSEA)) {
        stop("Path to the GSEA directory does not exist")
    }

    path_to_GSEA <- paste0('"', path_to_GSEA, '"')

    # check for path to GSEA dataset .gmt
    if (is.null(GSEA_dataset)) {
        warning("Path to a GSEA dataset needs to be provided, only dryrun will
                work for testing")
        dont_run <- TRUE
        GSEA_dataset <- "test.gmt"
    }

    GSEA_dataset <- paste0('"', GSEA_dataset, '"')

    # check for GSRA ranked file (path or data.frame)
    if (is.null(GSEA_ranked_file)) {
        warning("A ranked gene file needs to be provided, only dryrun will work
                for testing")
        dont_run <- TRUE
        GSEA_ranked_file <- "my_ranked_file.rnk"
    } else if (inherits(GSEA_ranked_file, "character")) {
        message("The ranked list looks like a path to a file")
        if (!file.exists(GSEA_ranked_file)) {
            stop("Path to the ranked file does not exist")
        }
    } else if (inherits(GSEA_ranked_file, "data.frame")) {
        message("The ranked list looks like a data.frame")

        # write data.frame to temporary folder and use that path
        temp_location_path <- paste0(tempdir(), "/", "temp_rankfile.rnk")
        data.table::fwrite(
            x = GSEA_ranked_file,
            file = temp_location_path,
            sep = "\t", row.names = FALSE, col.names = FALSE
        )

        GSEA_ranked_file <- temp_location_path
    }


    # make sure all paths and files can be read by java
    output_folder <- paste0('"', output_folder, '"')
    name_analysis_folder <- paste0('"', name_analysis_folder, '"')

    # 1. identify operating system
    my_os <- get_os()

    # 2. create execution path
    operation <- "GSEAPreranked"

    if (my_os == "windows") {
        execution_path <- paste0(
            path_to_GSEA, "/", "gsea-cli.bat", " ", operation
        )
    } else {
        execution_path <- paste0(
            path_to_GSEA, "/", "gsea-cli.sh", " ", operation
        )
    }

    created_command <- sprintf(
        "%s \\
            -gmx %s \\
            -collapse %s \\
            -mode %s \\
            -norm %s \\
            -nperm %s \\
            -rnk %s \\
            -scoring_scheme %s \\
            -include_only_symbols true \\
            -make_sets true \\
            -plot_top_x %s \\
            -rnd_seed timestamp \\
            -set_max %s \\
            -set_min %s \\
            -zip_report false \\
            -rpt_label %s \\
            -out %s",
        execution_path,
        GSEA_dataset,
        collapse,
        mode,
        norm,
        nperm,
        GSEA_ranked_file,
        scoring_scheme,
        plot_top_x,
        set_max,
        set_min,
        name_analysis_folder,
        output_folder
    )


    if (isTRUE(dryrun) | isTRUE(dont_run)) {
        message("DRYRUN VERSION")
        print(created_command)
    } else {
        message("START GSEA RUN")

        system(created_command)
    }
}

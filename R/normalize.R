# Data Processing and Pipelining ------------------------------------------ #
# ADDING NEW METHODS
# - 1. Define the operation as a function. 
#      (See # internals / implementations section below)
#   ## if adding a new operation type ##
#   - 1.1  Define an extending param class
#   - 1.2 Create an internal function to generate the class with default params 
#         (See # internals / params setup)
#   - 1.3 Update the relevant factory function (normParam, scaleParam, 
#         adjustParam). Add method to to the `match.arg` and register it in the
#         switch statement, pointing at the new default params internal 
#         function.
# - 2. Add the S4 method for the matrix/param combination desired
# - 3. Update the method documentation (if a method was added). Both the
#      operation-specific documentation and under `process_param`, linking
#      to it.

# Documentation ####

#' @name processExpression
#' @title Expression Data Processing
#' @description
#' Perform data transformations, or set up chains of transformations and
#' operations to be applied to expression type data in the `giotto` object.
#' @param gobject `giotto` object
#' @inheritParams processData
#' @param expression_values character. Name of matrix to use
#' @param spat_unit character (optional). spatial unit to use
#' @param feat_type character (optional). feature type to use
#' @param return_gobject logical (optional). Whether to return the `gobject`.
#' When FALSE, the `exprObj` is returned instead.
#' @returns A `giotto` object when `return_gobject = TRUE`. Otherwise, an
#' `exprObj`
#' @seealso [process_param] for processing operations that can be performed
#' 
#' [processData()] for the lower level generic handling these operations
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' # single operation
#' processExpression(g, normParam("library"), name = "library")
#' 
#' # single operation with changed parameter
#' lib <- normParam("library")
#' lib$scalefactor = 1000
#' processExpression(g, lib, name = "library2")
#' 
#' # return the exprObj instead
#' processExpression(g, lib, name = "library2", return_gobject = FALSE)
#' 
#' # chained operation (this is the Giotto standard normalization)
#' processExpression(g,
#'     list(
#'         normParam("library"),
#'         normParam("log"),
#'         scaleParam("zscore", MARGIN = 1),
#'         scaleParam("zscore", MARGIN = 2)
#'     ),
#'     name = "scaled2"
#' )
#' @md
NULL

#' @name processData
#' @title Composable Data Processing
#' @description
#' Perform data transformations, or set up chains of transformations and
#' operations to be applied to matrix type data. `processData()` is a generic
#' for which methods can be defined off both `x` (the data to transform),
#' and `param` (the transform operation).
#' @param x data to transform
#' @param param S4 parameter class defining the transform operation and
#' params affecting it. Can also be a list of several of these objects, acting
#' as a pipeline.
#' @param name character (optional). [Object name][GiottoClass::giotto_schema]
#' to assign to the output. Default `name` changes based on `param` input:
#' * when `param` is `list` or `scaleParam`: `name = "scaled"`
#' * when `param` is `normParam`: `name = "normalized"`
#' * when `param` is `adjustParam`: `name = "custom"`
#' * when `param` is `osmFISHNormParam`: `name = "custom"`
#' * when `param` is `pearsonResidNormParam`: `name = "scaled"`
#' @param \dots additional params to pass
#' @examples
#' m <- matrix(c(0, 0, 3, 2, 0, 5, 4, 0, 0, 1, 12, 0), nrow = 3)
#' 
#' # single operation
#' lib_norm <- normParam("library")
#' lib_norm$scalefactor <- 5000 # alter a default param of library norm
#' processData(m, lib_norm)
#' 
#' # chained operations
#' log_norm <- normParam("log")
#' zscore_rows <- scaleParam("zscore", MARGIN = 1)
#' zscore_cols <- scaleParam("zscore")
#' # this is essentially the same as the default giotto normalization
#' # only difference is the library norm scalefactor change.
#' processData(m, list(lib_norm, log_norm, zscore_rows, zscore_cols))
#' @seealso [process_param] for processing operations that can be performed
#' through `processData()`
#' @seealso [processExpression()] for the way to use this framework with the 
#' `giotto` object
#' @returns The same class as `x`
#' @md
NULL

#' @name process_param
#' @title Data Processing Parameter Classes
#' @description Data processing operations in Giotto Suite can be divided into
#' normalization, scaling, and adjustments. These operations can be selected
#' via the factory functions `normParam()`, `scaleParam()`, and `adjustParam()`, 
#' respectively.
#' 
#' Requested operations are generated as method-specific param classes that
#' contain all the parameters needed to perform them, editable through `$<-`.
#' @param method character. Name of method to use. See details.
#' @param \dots (optional) Additional named parameters relevant to the param 
#' class.
#' @section normParam methods: 
#' 
#' * [`"default"`][norm_default] - default Giotto normalizations steps 
#' (library + log norms)
#' * [`"library"`][norm_library] - library normalization
#' * [`"log"`][norm_log] - log normalization
#' * [`"osmfish"`][norm_osmfish] - osmfish normalization method
#' * [`"pearson"`][norm_pearson] - Lause/Kobak 2020 pearson residuals
#' normalization
#' * [`"quantile"`][norm_quantile] - quantile normalization
#' * [`"tf-idf"`][norm_tfidf] - Term Frequency-Inverse Document Frequency
#' * [`"l2"`][norm_l2] - L2 normalization (also known as Euclidean
#' normalization)
#' * [`"arcsinh"`][norm_arcsinh] - arcsinh transformation
#' 
#' @section scaleParam methods: 
#' 
#' * [`"default"`][scale_default] - default Giotto scaling steps (scale along
#' features then cells)
#' * [`"zscore"`][scale_zscore] - essentially the same as `base::scale()`, but
#' with a `MARGIN` param allowing scaling long either cols or rows
#' 
#' @section adjustParam methods:
#' 
#' * [`"limma"`][adjust_limma] - limma batch correction
#' 
#' @section thresholdParam methods:
#' 
#' * [`"binarize"`][threshold_binarize] - data binarization (matrices and 
#' rasters)
#' * [`"minmax"`][threshold_minmax] - value restriction/clamping (matrices 
#' and rasters)
#' @details
#' Generated params are S4 objects inheriting from `processParam` and one of 
#' `normParam`, `scaleParam`, and `adjustParam`.
#' @seealso [processData()] for the generic used to apply these params
#' @seealso [processExpression()] for the way to use this framework with the 
#' `giotto` object
#' @md
NULL

#' @name norm_default
#' @title Default Giotto Normalization
#' @description
#' Expression matrix normalization method.
#' 
#' Steps:
#' 
#' 1. [Total library size][norm_library] normalization and scaling by 
#' a custom scale-factor.
#' 2. [Log][norm_log] transformation of data.
#' 
#' @section params: 
#' 
#' \tabular{ll}{
#'   `library_size_norm` \tab logical (default = `TRUE`). whether to perform
#'   library size normalization \cr
#'   `scalefactor` \tab numeric (default = 6000). Scalefactor to use after
#'   library size normalization. (skipped if `library_size_norm = FALSE`) \cr
#'   `log_norm` \tab logical (default = `TRUE`). Whether to transform values to
#'   log-scale. \cr
#'   `log_offset` \tab numeric (default = 1). If `log_norm = TRUE`, offset
#'   value to add to expression values to avoid `log(0)`. Note: for sparse-like
#'   matrices (e.g. `dgCMatrix`, 'dbSparseMatrix'), only `log_offset = 1` is
#'   supported (log1p) \cr
#'   `logbase` \tab numeric (default = 2). If `log_norm = TRUE`, log base to
#'   use to log normalize expression values
#' }
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
#' @md
NULL

#' @name norm_library
#' @title Library Size Normalization
#' @description
#' Normalize expression matrix for total library size and then scale by
#' a custom scalefactor.
#' 
#' This method does not work well when any cells/samples
#' have a library size of 0, so filtering prior to this is recommended.
#' 
#' \deqn{\LARGE
#' x'_{i,j} = \frac{x_{i,j}}{\sum_{i} x_{i,j}} \times k
#' }
#' Where:
#' 
#' * (\eqn{x_{i,j}}) is the raw count for feature \eqn{i} in sample \eqn{j}
#' * (\eqn{x'_{i,j}}) is the library normalized and scaled expression value for
#' feature \eqn{i} in sample \eqn{j}
#' * (k) is a scalefactor applied after normalization
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `scalefactor` \tab numeric (default = 6000). Scalefactor to use after 
#'   library size normalization. Expressed as ***k*** in the above equation
#' }
#' @md
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
NULL

#' @name norm_log
#' @title Log Normalization
#' @description
#' Apply a log normalization
#' 
#' \deqn{\LARGE
#' x'_{i,j} = \frac{\log(x_{i,j} + b)}{\log(a)}
#' }
#' Where:
#' 
#' * (\eqn{x_{i,j}}) is the raw count for feature \eqn{i} in sample \eqn{j}
#' * (\eqn{x'_{i,j}}) is the log normalized expression value for feature 
#' \eqn{i} in sample \eqn{j}
#' * (\eqn{a}) is the log base
#' * (\eqn{b}) is an offset value
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `base` \tab numeric (default = 2) log base to use. Expressed as \eqn{a} in
#'   the above equation. \cr
#'   `offset` \tab numeric (default = 1). Offset to add to expression values to
#'   avoid \eqn{\log(0)}. Expressed as \eqn{b} in the above equation.
#'   For sparse-like matrices and `dbMatrix`, only `offset = 1` is supported;
#'   other values would change implicit zeros and require densification.
#' }
#' @md
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
NULL

#' @name norm_osmfish
#' @title osmFISH Normalization
#' @description
#' Normalization method as provided by the osmFISH paper
#' 
#' Steps:
#' 
#' 1. First normalize genes, for each gene divide the counts by the total gene 
#' count and multiply by the total number of genes.
#' 2. Next normalize cells, for each cell divide the normalized gene counts by
#' the total counts per cell and multiply by the total number of cells.
#' 
#' \deqn{\LARGE
#' x'_{i,j} = \frac{x_{i,j}}{\sum_j x_{i,j}} \times n_{\text{features}}
#' }
#'
#' \deqn{\LARGE
#' x''_{i,j} = \frac{x'_{i,j}}{\sum_i x'_{i,j}} \times n_{\text{samples}}
#' }
#' 
#' Where:
#' 
#' * (\eqn{x_{i,j}}) is the raw count for feature \eqn{i} in sample \eqn{j}
#' * (\eqn{x'_{i,j}}) is the feature normalized expression value
#' * (\eqn{x''_{i,j}}) is the final normalized expression value after both
#' feature and cell normalization
#' * (\eqn{n_{\text{samples}}}) is the total number of cells
#' (columns in matrix)
#' * (\eqn{n_{\text{features}}}) is the total number of cells
#' (rows in matrix)
#' 
#' @section params:
#' None
#' @md
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
NULL

#' @name norm_pearson
#' @title Lause/Kobak Pearson Residuals Normalization
#' @description
#' Calculate Pearson residuals with a dispersion adjustment, to identify cells
#' that deviate significantly from what would be expected under independence. 
#' The normalization divides by the standard deviation of the difference, which
#' is adjusted by the dispersion parameter θ.
#' 
#' This normalization is designed for detection of highly variable features and
#' dimension reduction and clustering.
#' 
#' \deqn{\LARGE
#' z_{i,j} = \frac{x_{i,j} - \mu_{i,j}}{\sqrt{\mu_{i,j} + \mu_{i,j}^2 / \theta}}
#' }
#'
#' \deqn{\LARGE
#' \mu_{i,j} = \frac{r_i \cdot c_j}{N}
#' }
#' 
#' Where:
#' * (\eqn{x_{i,j}}) is the raw count for feature \eqn{i} in sample \eqn{j}
#' * (\eqn{\mu_{i,j}}) is the expected value under the model
#' * (\eqn{r_i}) is \eqn{\sum_j x_{i,j}}
#' * (\eqn{c_j}) is \eqn{\sum_i x_{i,j}}
#' * (\eqn{N}) is \eqn{\sum_{i,j} x_{i,j}}
#' * (\eqn{\theta}) is a dispersion parameter
#' * (\eqn{z_{i,j}}) is the Pearson residual clipped to the range 
#' \eqn{[-\sqrt{n}, \sqrt{n}]} where \eqn{n} is the number of columns. This is 
#' done to prevent extreme values from dominating the analysis.
#' 
#' # Note
#' Scaling is not recommended after this normalization since it is already
#' transforming the data to z-score-like values with a dispersion adjustment.
#' It is also not recommended to use this with DGE analysis.
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `theta` \tab dispersion parameter expressed as \eqn{\theta} in the above
#'   formula
#' }
#' 
#' @references Lause, J., Berens, P. & Kobak, D. Analytic Pearson residuals for
#' normalization of single-cell RNA-seq UMI data. Genome Biol 22, 258 (2021).
#' https://doi.org/10.1186/s13059-021-02451-7
#' @md
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
NULL

#' @name norm_quantile
#' @title Quantile Normalization
#' @description
#' Quantile normalization makes the statistical distribution of values in each
#' column identical by replacing the original values with the mean of the
#' values at the same rank across all columns. This removes technical variation
#' while preserving relative differences between features.
#'
#' Steps:
#' 1. Rank the values within each column (average taken in case of ties)
#' 2. Calculate the mean of values at the same rank across all columns
#' 3. Replace each value with the mean value corresponding to its rank
#'
#' \deqn{\LARGE
#' q_{i,j} = \bar{x}_{rank(i,j)}
#' }
#'
#' Where:
#' * (\eqn{rank(i,j)}) is the rank of feature \eqn{i} within column \eqn{j}
#' * (\eqn{\bar{x}_{r}}) where \eqn{r = rank(i,j)} is the mean of values with
#' rank \eqn{r} across all columns
#' * (\eqn{q_{i,j}}) is the quantile-normalized value
#' 
#' # Note
#' Library normalization and log normalization is recommended prior to this
#' normalization.
#' 
#' @section params:
#' None
#'
#' @references Bolstad, B.M., Irizarry, R.A., Astrand, M. et al. A comparison of
#' normalization methods for high density oligonucleotide array data based on
#' variance and bias. Bioinformatics 19, 185–193 (2003).
#' https://doi.org/10.1093/bioinformatics/19.2.185
#' @md
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
NULL

#' @name norm_tfidf
#' @title TF-IDF Normalization
#' @description
#' TF-IDF (Term Frequency-Inverse Document Frequency) normalization is borrowed 
#' from natural language processing to identify features that are highly
#' expressed in specific samples but not widely expressed across the entire 
#' dataset.
#' 
#' There are several different implementations that apply log or binarization
#' to different terms. `sub_method = c(1:3)` and `dgCMatrix` optimizations are
#' based on the ArchR implementations.
#' 
#' \deqn{\LARGE
#' TF_{i,j} = \frac{x_{i,j}}{\sum_{i} x_{i,j}}
#' }
#' 
#' \deqn{\LARGE
#' IDF_{i} = \frac{n_{samples}}{\sum_{j} x_{i,j}}
#' }
#' 
#' \deqn{\LARGE
#' IBDF_{i} = \frac{n_{samples}}{1 + n_{samples \: where \: feature \: i > 0}}
#' }
#' 
#' **Implementations** (`sub_method`):
#' 
#' \deqn{\large
#' (default) \quad TFIDF_{i,j} = TF_{i,j} \times \log(IBDF_{i} + 1)
#' }
#' 
#' \deqn{\large
#' (1) \quad TFIDF_{i,j} = TF_{i,j} \times \log(IDF_{i} + 1)
#' }
#' 
#' \deqn{\large
#' (2) \quad TFIDF_{i,j} = \log(TF_{i,j} \times IDF_{i} \times S + 1) \quad 
#' }
#' 
#' \deqn{\large
#' (3) \quad TFIDF_{i,j} = \log(TF_{i,j} + 1) \times \log(IDF_{i} + 1)
#' }
#' 
#' Where:
#' * (\eqn{x_{i,j}}) is the raw count for feature \eqn{i} in sample \eqn{j}
#' * (\eqn{TF_{i,j}}) is the term frequency of feature \eqn{i} in sample \eqn{j}
#' * (\eqn{IDF_{i}}) is the inverse document frequency of feature \eqn{i}
#' * (\eqn{IBDF_{i}}) is the inverse binarized document frequency of feature \eqn{i}
#' * (\eqn{TFIDF_{i,j}}) is the final TF-IDF normalized value
#' * (\eqn{S}) is a scalefactor (default = 10000)
#' 
#' # Note
#' [L2][norm_l2] normalization is commonly performed after TF-IDF normalization
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `sub_method` \tab Either numeric 1, 2, or 3 or "default". Determines which
#'   set of defaults to use during the TF-IDF calculation. Methods 1-3 map to
#'   the same `LSIMethod` settings in {ArchR}. See sub_method section below.\cr
#'   `log_tf` \tab logical (overrides `sub_method` defaults). Whether to log
#'   transform TF values (includes a +1 offset).\cr
#'   `log_idf` \tab logical (overrides `sub_method` defaults). Whether to log
#'   transform IDF values (includes a +1 offset).\cr
#'   `log_tf_idf` \tab logical (overrides `sub_method` defaults). Whether to log
#'   transform the TF-IDF value (also applies a scalefactor (\eqn{s}) and +1
#'   offset before the log operation).\cr
#'   `binarized_rowsums` \tab logical (overrides `sub_method` defaults).
#'   Whether to calculate IBDF instead of IDF, where the calculation is based
#'   on the presence of a feature as opposed to its count.\cr
#'   `scalefactor` \tab numeric (default = 10000). A scalefactor used when
#'   `log_tf_idf = TRUE`.
#' }
#' 
#' @section `sub_method`:
#' `sub_method` can be one of `"default"` or any of the other implementations
#' from 1 to 3. These apply some defaults to the way that TF-IDF is calculated.
#' The individual `log_` params will override these defaults.
#' 
#' * `"default"` - default Giotto implementation:
#'   * `log_idf = TRUE`
#'   * `binarized_rowsums = TRUE`
#' * `1` - Method introduced in Cusanovich et al. 2018.
#'   * `log_idf = TRUE`
#' * `2` - Method introduced in Stuart et al. 2021.
#'   * `log_tf_idf = TRUE`
#' * `3` - Method 3 in {ArchR} `iterativeLSI()`
#'   * `log_tf = TRUE`
#'   * `log_idf = TRUE`
#' 
#' @md
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
#' @references
#' Cusanovich, D., Reddington, J., Garfield, D. et al. The cis-regulatory 
#' dynamics of embryonic development at single-cell resolution.
#' Nature 555, 538–542 (2018). https://doi.org/10.1038/nature25981
#' 
#' Stuart T, Srivastava A, Madad S, Lareau CA, Satija R. Single-cell chromatin
#' state analysis with Signac. Nat Methods. 2021 Nov;18(11):1333-1341.
#' doi: 10.1038/s41592-021-01282-5.
#' 
#' Granja JM, Corces MR, Pierce SE, Bagdatli ST, Choudhry H, Chang HY,
#' Greenleaf WJ. ArchR is a scalable software package for integrative
#' single-cell chromatin accessibility analysis.
#' Nat Genet. 2021 Mar;53(3):403-411. doi: 10.1038/s41588-021-00790-6.
#' @examples
#' e <- GiottoData::loadSubObjectMini("exprObj")
#' processData(e, normParam("tf-idf"))
#' processData(e, normParam("tf-idf", sub_method = 1))
#' processData(e, normParam("tf-idf", sub_method = 2))
#' processData(e, normParam("tf-idf", sub_method = 3))
NULL

#' @name norm_l2
#' @title L2 Normalization
#' @description
#' L2 normalization (also known as Euclidean normalization) scales each column
#' (sample) in the expression matrix to have unit Euclidean length. This
#' process makes samples with different sequencing depths more comparable and
#' improves the performance of distance-based analyses.
#' 
#' \deqn{\LARGE
#' x'_{i,j} = \frac{x_{i,j}}{\sqrt{\sum_{i} x_{i,j}^2}}
#' }
#' 
#' Where:
#' * (\eqn{x_{i,j}}) is the expression value for feature \eqn{i} in sample \eqn{j}
#' * (\eqn{x'_{i,j}}) is the L2-normalized expression value
#' 
#' # Note
#' L2 normalization can be applied to raw data, but is most commonly used after 
#' other normalization methods such as TF-IDF or log normalization to standardize
#' sample-to-sample comparisons.
#' 
#' @section params:
#' None
#' 
#' @family normalization parameters
#' @seealso [process_param]
#' @returns normalized object
#' @md
NULL

#' @name norm_arcsinh
#' @title Arcsinh Normalization
#' @description
#' A normalization commonly used with intensity-based data (CODEX, CyCIF, IMC).
#' It effectively handles a wide dynamic range and zero/near-zero values while
#' preserving the relative differences between signals of different intensities.
#' 
#' \deqn{\LARGE
#' x'_{i,j} = \operatorname{arcsinh}\left({\frac{x_{i,j}}{c}}\right)
#' }
#' 
#' Where:
#' * (\eqn{x_{i,j}}) is the raw intensity for feature \eqn{i} in sample \eqn{j}
#' * (\eqn{x'_{i,j}}) is the normalized intensity for feature \eqn{i} in
#' sample \eqn{j}
#' * (\eqn{c}) is a cofactor that determines the degree of transformation
#' 
#' # Note
#' The cofactor \eqn{c} prevents over-amplification of small values and allows
#' better differentiation of signals at different intensities.
#' 
#' Common values to use are:
#' * `c = 5` for fluorescence imaging-based proteomics (CODEX, CyCIF) 
#' * `c = 1` or `5` for mass-cytometry-based imaging (IMC).
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `c` \tab numeric (default = 5). Expressed as \eqn{c} in the above
#'   formula.
#' }
#' 
#' @family normalization parameters
#' @seealso [process_param()]
#' @returns normalized object
#' @md
NULL

#' @name scale_default
#' @title Default Giotto Scaling
#' @description
#' 2 step [z-scoring][scale_zscore] along features and samples
#' @section params: 
#' 
#' \tabular{ll}{
#'   `scale_feats` \tab logical (default = `TRUE`) Whether to scale across
#'   features \cr
#'   `scale_cells` \tab logical (default = `TRUE`) Whether to scale across
#'   cells/samples \cr
#'   `scale_order` \tab character. One of either `"first_feats"` or 
#'   `"first_cells"`. When both `scale_feats` and `scale_cells` are `TRUE`,
#'   determines the order in which the 2 scaling operations are performed. \cr
#'   `verbose` \tab logical (default = `TRUE`) Whether to be verbose
#' }
#' 
#' @md
#' @family scaling parameters
#' @seealso [process_param]
#' @returns scaled object
NULL

#' @name scale_zscore
#' @title Z Score Scaling
#' @description
#' Wrapper around `base::scale()` to make it compatible with the
#' [processData()] framework. Additionally provides a `MARGIN` param.
#' 
#' \deqn{\LARGE
#' z_{i,j} = \frac{x_{i,j} - \mu_i}{\sigma_i}
#' }
#'
#' Where:
#' * \eqn{x_{i,j}} is the original value for feature \eqn{i} in sample \eqn{j}
#' * \eqn{\mu_i} is the mean of feature \eqn{i} across all samples
#' * \eqn{\sigma_i} is the standard deviation of feature \eqn{i} across all 
#' samples
#' * \eqn{z_{i,j}} is the resulting scaled value
#' 
#' @section params: 
#' 
#' \tabular{ll}{
#'   `scale` \tab logical (default = `TRUE`) Whether to scale values \cr
#'   `center` \tab logical (default = `TRUE`) Whether to center values\cr
#'   `MARGIN` \tab numeric. Either 1 (rows) or 2 (cols). Direction along which
#'   to perform the operation.
#' }
#' @md
#' @family scaling parameters
#' @returns scaled object
#' @seealso [process_param]
NULL

#' @name adjust_limma
#' @title Limma Batch Correction
#' @description
#' Batch effect removal via [limma::removeBatchEffect()]
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `batch_columns` \tab [svkey][GiottoClass::svkey()] (optional) Up to two
#'   columns of information from a Giotto object with information indicating
#'   batches to remove the effects of. \cr
#'   `covariate_columns` \tab [svkey][GiottoClass::svkey()] (optional) Columns
#'   of information from a Giotto object with information indicating covariates
#'   to regress out.
#' }
#' @returns limmaAdjustParam
#' @examples
#' limma <- adjustParam("limma")
#' limma$covariate_columns <- svkey(feats = c("nr_feats", "total_expr"))
#' 
#' g <- GiottoData::loadGiottoMini("visium")
#' processExpression(g, limma, name = "limma")
#' @family adjustment parameters
#' @seealso [process_param]
#' @md
NULL

#' @name threshold_binarize
#' @title Data Binarization
#' @description
#' Binarize values to 0 and 1 based on a minimal value. For matrices, the
#' default threshold is 0. For rasters, the default is a value determined
#' through sampled (5e5 pixels) otsu.
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `threshold` \tab numeric (optional) Values above or equal to the threshold
#'   will be set to 1. Below will be set to 0. If not provided, defaults to 0
#'   for matrices and a value determined by otsu for rasters. \cr
#'   `drop0` \tab logical (only for `dgCMatrix`, default = FALSE) Whether to
#'   run [Matrix::drop0] after binarization 
#' }
#' @returns binarizeThreshParam
#' @examples
#' e <- GiottoData::loadSubObjectMini("exprObj")
#' # also works with matrix classes
#' bin_e <- processData(e, thresholdParam("binarize"))
#' force(bin_e)
#' 
#' gimg <- GiottoData::loadSubObjectMini("giottoLargeImage")
#' # also works with SpatRasters
#' bin_img <- processData(gimg, thresholdParam("binarize"))
#' plot(bin_img)
#' @family threshold parameters
#' @seealso [process_param]
#' @md
NULL

#' @name threshold_minmax
#' @title Value MinMax Restriction/Clamping
#' @description
#' Set an `upper` and `lower` bound for the data. Values above `upper` will be
#' set to the `upper` value. Values below `lower` will be set to the `lower`
#' value.
#' 
#' @section params:
#' 
#' \tabular{ll}{
#'   `upper` \tab numeric (default = Inf) highest acceptable value. Values
#'   above this will be set to the same as `upper`. \cr
#'   `lower` \tab numeric (default = -Inf) lowest acceptable value. Values
#'   below this will be set to the same as `lower`.\cr
#'   `values` \tab logical (`SpatRaster` only) If `FALSE` values outside the
#'   clamping range become `NA`, if `TRUE`, they get the extreme values
#' }
#' @returns minmaxThreshParam
#' @examples
#' e <- GiottoData::loadSubObjectMini("exprObj")
#' # also works with matrix classes
#' max_e <- processData(e, thresholdParam("minmax", upper = 6))
#' force(max_e)
#' 
#' gimg <- GiottoData::loadSubObjectMini("giottoLargeImage")
#' # also works with SpatRasters
#' mm_img <- processData(gimg,
#'     thresholdParam("minmax", lower = 30, upper = 100)
#' )
#' plot(mm_img)
#' @family threshold parameters
#' @seealso [process_param]
#' @md
NULL


# VIRTUAL classes ####
setClass("normParam", contains = c("VIRTUAL", "processParam"))
setClass("scaleParam", contains = c("VIRTUAL", "processParam"))
setClass("adjustParam", contains = c("VIRTUAL", "processParam"))
setClass("threshParam", contains = c("VIRTUAL", "processParam"))

# access ####
#' @export
.DollarNames.scaleParam <- function(x, pattern) {
    names(x@param)
}
#' @export
.DollarNames.normParam <- function(x, pattern) {
    names(x@param)
}
#' @export
.DollarNames.adjustParam <- function(x, pattern) {
    names(x@param)
}
#' @export
.DollarNames.thresholdParam <- function(x, pattern) {
    names(x@param)
}

# extending method classes ####

#' @rdname process_param
#' @exportClass defaultNormParam
setClass("defaultNormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass libraryNormParam
setClass("libraryNormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass logNormParam
setClass("logNormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass osmFISHNormParam
setClass("osmFISHNormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass pearsonResidNormParam
setClass("pearsonResidNormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass quantileNormParam
setClass("quantileNormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass tfidfNormParam
setClass("tfidfNormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass l2NormParam
setClass("l2NormParam", contains = "normParam")
#' @rdname process_param
#' @exportClass arcsinhNormParam
setClass("arcsinhNormParam", contains = "normParam")

#' @rdname process_param
#' @exportClass defaultScaleParam
setClass("defaultScaleParam", contains = "scaleParam")
#' @rdname process_param
#' @exportClass zscoreScaleParam
setClass("zscoreScaleParam", contains = "scaleParam")

#' @rdname process_param
setClass("limmaAdjustParam", contains = "adjustParam")

#' @rdname process_param
setClass("binarizeThreshParam", contains = "threshParam")
#' @rdname process_param
setClass("minmaxThreshParam", contains = "threshParam")

# allMatrix signature ####
setClassUnion("allMatrix", members = c("matrix", "Matrix"))




# param factories ####

#' @rdname process_param
#' @export
normParam <- function(method = "default", ...) {
    method <- match.arg(tolower(method),
        c("default", "library", "log", "osmfish", "pearson", "quantile", 
        "tf-idf", "l2", "arcsinh")
    )
    switch(method,
        "default" = .norm_param_default(...),
        "library" = .norm_param_lib(...),
        "log" = .norm_param_log(...),
        "osmfish" = .norm_param_osmfish(...),
        "pearson" = .norm_param_pears_resid(...),
        "quantile" = .norm_param_quantile(...),
        "tf-idf" = .norm_param_tfidf(...),
        "l2" = .norm_param_l2(...),
        "arcsinh" = .norm_param_arcsinh(...)
    )
}

#' @rdname process_param
#' @export
scaleParam <- function(method = "default", ...) {
    method <- match.arg(tolower(method),
        c("default", "zscore")
    )
    switch(method,
        "default" = .scale_param_default(...),
        "zscore" = .scale_param_zscore(...)
    )
}

#' @rdname process_param
#' @export
adjustParam <- function(method = "limma", ...) {
    method <- match.arg(tolower(method),
        c("limma")
    )
    switch(method,
        "limma" = .adjust_param_limma(...)
    )
}

#' @rdname process_param
#' @export
thresholdParam <- function(method = "binarize", ...) {
    method <- match.arg(tolower(method),
        c("binarize", "minmax")
    )
    switch(method,
        "binarize" = .thresh_param_binarize(...),
        "minmax" = .thresh_param_minmax(...)
    )
}



# methods ####

# * ANY ####

setMethod("processData",
signature(x = "ANY", param = "ANY"), function(x, param, ...) {
    stop(wrap_txtf(
        "param of class '%s' is not recognized for use with '%s'", 
        class(param), class(x)),
        call. = FALSE)
})

# * giottoLargeImage ####

# TODO make these delayed operations

#' @rdname processData
setMethod("processData",
    signature(x = "giottoLargeImage", param = "list"),
    function(x, param, name = NULL, ...) {
        x[] <- processData(x[], param, ...)
        if (!is.null(name)) {
            objName(x) <- name
        }
        x
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "giottoLargeImage", param = "processParam"),
    function(x, param, name = NULL, ...) {
        x[] <- processData(x[], param, ...)
        if (!is.null(name)) {
            objName(x) <- name
        }
        x
    }
)

# * exprObj ####

#' @rdname processData
setMethod("processData",
    signature(x = "exprObj", param = "list"),
    function(x, param, name = "scaled", ...) {
        x[] <- processData(x[], param, ...)
        objName(x) <- name
        return(x)
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "exprObj", param = "normParam"), 
    function(x, param, name = "normalized", ...) {
        x[] <- processData(x[], param, ...)
        objName(x) <- name
        return(x)
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "exprObj", param = "scaleParam"),
    function(x, param, name = "scaled", ...) {
        x[] <- processData(x[], param, ...)
        objName(x) <- name
        return(x)
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "exprObj", param = "adjustParam"),
    function(x, param, name = "custom", ...) {
        x[] <- processData(x[], param, ...)
        objName(x) <- name
        return(x)
    }
)

setMethod("processData",
    signature(x = "exprObj", param = "processParam"),
    function(x, param, name = "custom", ...) {
        x[] <- processData(x[], param, ...)
        objName(x) <- name
        return(x)
    }
)

# specialized handling for osmfish
setMethod("processData",
    signature(x = "exprObj", param = "osmFISHNormParam"), 
    function(x, param, name = "custom", ...) {
        if (!featType(x) %in% c("rna", "RNA")) {
            warning("Caution: osmFISH normalization was developed for RNA 
                    in situ data",
                    call. = FALSE)
        }
        x[] <- processData(x[], param, ...)
        objName(x) <- name
        return(x)
    }
)

# specialized handling for pearson residual
setMethod("processData",
    signature(x = "exprObj", param = "pearsonResidNormParam"), 
    function(x, param, name = "scaled", ...) {
        if (!featType(x) %in% c("rna", "RNA")) {
            warning("Caution: pearson residual normalization was developed 
                    for RNA count normalization",
                    call. = FALSE)
        }
        x[] <- processData(x[], param, ...)
        objName(x) <- name
        return(x)
    }
)


# * matrix ####

# ** param list ####

#' @rdname processData
setMethod("processData",
    signature(x = "SpatRaster", param = "list"),
    function(x, param, ...) {
        for (p in param) {
            x <- processData(x, p, ...)
        }
        return(x)
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "allMatrix", param = "list"),
    function(x, param, ...) {
        for (p in param) {
            x <- processData(x, p, ...)
        }
        return(x)
    }
)

# ** norm ------------------ ####
# *** library norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "libraryNormParam"),
    function(x, param, ...) {
        libsizes <- colSums_flex(x)
        .libzero_warn(libsizes = libsizes)
        t_flex(t_flex(x) / libsizes) * param$scalefactor
    }
)
setMethod("processData",
    signature(x = "dgCMatrix", param = "libraryNormParam"),
    function(x, param, ...) {
        libsizes <- colSums_flex(x)
        .libzero_warn(libsizes = libsizes)
        # x / colSums(x) equivalent for dgc
        x@x <- .dgc_div_csum_sparse_vector(x, colsums = libsizes) * 
            param$scalefactor
        x
    }
)
# *** log norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "logNormParam"),
    function(x, param, ...) {
        # offset != 1 is invalid for sparse: log(0 + offset) != 0
        if (.is_sparse_like(x) && !isTRUE(param$offset == 1)) {
            stop("`offset != 1` is not supported for sparse matrices. ",
                 "Use offset = 1 (log1p) to preserve sparsity.")
        }
        if (.is_sparse_like(x)) {
            return(log1p(x) / log(param$base)) # sparsity preserved
        }
        return(log(x + param$offset) / log(param$base))
    }
)
# *** osmFISH norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "osmFISHNormParam"),
    function(x, param, ...) {
        # 1. normalize raw expr per gene with scale-factor equal to number of genes
        norm_feats <- (x / rowSums_flex(x)) * nrow(x)
        # 2. normalize per cells with scale-factor equal to number of cells
        t_flex((t_flex(norm_feats) / colSums_flex(norm_feats)) * ncol(x))
    }
)
# *** pearson norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "pearsonResidNormParam"),
    function(x, param, ...) {
        .pears_resid_citation(verbose = param$verbose)
        if (methods::is(x, "HDF5Matrix")) { # avoid importing this class
            .csums <- .csum_nodrop.HDF5Matrix
            .rsums <- .rsum_nodrop.HDF5Matrix
        } else {
            .csums <- .csum_nodrop.Matrix
            .rsums <- .rsum_nodrop.Matrix
        }
        .prnorm(
            x = x, 
            theta = param$theta, 
            .csums = .csums,
            .rsums = .rsums
        )
    }
)
# *** quantile norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "quantileNormParam"),
    function(x, param, ...) {
        .qnorm(x)
    }
)
# *** tf-idf norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "tfidfNormParam"),
    function(x, param, ...) {
        p <- param
        p <- .norm_param_tfidf_method_defaults(p)
        # rely on defaults in .tf_idf_norm() for tf_fun and log_scale_fun
        
        # finalize args list
        p$x <- x
        
        do.call(.tf_idf_norm, args = p@param)
    }
)
setMethod("processData",
    signature(x = "dgCMatrix", param = "tfidfNormParam"),
    function(x, param, ...) {
        p <- param
        p <- .norm_param_tfidf_method_defaults(p)
        # optimizations for dgCMatrix
        tf_fun <- function(mat) {
            # x / colSums(x) equivalent
            mat@x <- .dgc_div_csum_sparse_vector(mat)
            mat
        }
        log_scale_fun <- function(mat, scalef) {
            mat@x <- log(mat@x * scalef + 1)
            mat
        }
        mat_log_fun <- function(mat) {
            mat@x <- log(mat@x + 1)
            mat
        }
        
        # finalize args list
        p$x <- x
        p$tf_fun <- tf_fun
        p$mat_log_fun <- mat_log_fun
        p$log_scale_fun <- log_scale_fun
        p$sub_method <- NULL # cleanup

        do.call(.tf_idf_norm, args = p@param)
    }
)
# *** default norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "defaultNormParam"),
    function(x, param, ...) {
        plist <- list()
        # 1. library size normalization
        if (isTRUE(param$library_size_norm)) {
            plist <- c(plist, normParam("library", 
                scalefactor = param$scalefactor))
        }
        # 2. log normalize
        if (isTRUE(param$log_norm)) {
            plist <- c(plist, normParam("log", 
                logbase = param$logbase, 
                log_offset = param$log_offset)
            )
        }
        processData(x, plist, ...)
    }
)
# *** L2 norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "l2NormParam"),
    function(x, param, ...) {
        .l2_norm(x)
    }
)
# *** arcsinh norm ####
setMethod("processData",
    signature(x = "allMatrix", param = "arcsinhNormParam"),
    function(x, param, ...) {
        .arcsinh_norm(x, c = param$c)
    }
)
setMethod("processData",
    signature(x = "SpatRaster", param = "arcsinhNormParam"),
    function(x, param, ...) {
        .arcsinh_norm(x, c = param$c)
    }
)


# ** threshold ------------ #### 

# *** binarization thresh ####
setMethod("processData",
    signature(x = "allMatrix", param = "binarizeThreshParam"),
    function(x, param, ...) {
        threshold <- param$threshold %null% 0
        bool_mat <- x > threshold
        x[TRUE] <- 0L
        x[bool_mat] <- 1L
        x
    }
)
setMethod("processData",
    signature(x = "dgCMatrix", param = "binarizeThreshParam"),
    function(x, param, ...) {
        threshold <- param$threshold %null% 0
        bool <- x@x > threshold
        drop0 <- param$drop0 %null% FALSE
        # integers not supported for dgCMatrix @x slot
        x@x <- rep.int(0, length(x@x))
        x@x[bool] <- 1
        if (drop0) x <- Matrix::drop0(x)
        x
    }
)
setMethod("processData",
    signature(x = "SpatRaster", param = "binarizeThreshParam"),
    function(x, param, ...) {
        lyrnames <- names(x)
        if (is.null(param$threshold)) {
            threshold <- .otsu(x)
        } else {
            threshold <- param$threshold
        }
        x <- terra::app(x, function(val) {
            ifelse(val >= threshold, 1, 0)
        })
        names(x) <- lyrnames
        x
    }
)

# *** minmax thresh ####
setMethod("processData",
    signature(x = "allMatrix", param = "minmaxThreshParam"),
    function(x, param, ...) {
        upper <- param$upper
        lower <- param$lower
        if (!is.infinite(upper)) {
            x[x > upper] <- upper
        }
        if (!is.infinite(lower)) {
            x[x < lower] <- lower
        }
        x
    }
)
setMethod("processData",
    signature(x = "SpatRaster", param = "minmaxThreshParam"),
    function(x, param, ...) {
        p <- param@param # pull param list
        p$values <- p$values %null% TRUE
        p$x <- x
        x <- do.call(terra::clamp, args = p)
        x
    }
)

# ** scale ----------------- ####
# *** zscore scale ####
setMethod("processData",
    signature("allMatrix", param = "zscoreScaleParam"), 
    function(x, param, ...) {
        if (!param$MARGIN %in% c(1, 2)) {
            stop(
            "processData zscore: 'MARGIN' must be either 1 (rows) or 2 (cols)", 
            call. = FALSE)
        }
        if (param$MARGIN == 1) x <- t_flex(x)
        x <- standardise_flex(x, center = param$center, scale = param$scale)
        if (param$MARGIN == 1) x <- t_flex(x)
        return(x)
    })
# *** default scale ####
setMethod("processData",
    signature(x = "allMatrix", param = "defaultScaleParam"),
    function(x, param, ...) {
        plist <- list()
        s1 <-scaleParam("zscore", center = TRUE, scale = TRUE, MARGIN = 1)
        s2 <-scaleParam("zscore", center = TRUE, scale = TRUE, MARGIN = 2)
        if (isTRUE(param$scale_feats) && isTRUE(param$scale_cells)) {
            scale_order <- match.arg(
                param$scale_order,
                choices = c("first_feats", "first_cells")
            )
            if (scale_order == "first_feats") {
                vmsg(.v = param$verbose, "first scale feats and then cells")
                plist <- c(plist, s1, s2)
            } else if (scale_order == "first_cells") {
                vmsg(.v = param$verbose, "first scale cells and then feats")
                plist <- c(plist, s2, s1)
            } else {
                stop("processData defaultNormParam: scale order must be given", 
                    call. = FALSE)
            }
        } else if (isTRUE(param$scale_feats)) {
            plist <- c(plist, s1)
        } else if (isTRUE(param$scale_cells)) {
            plist <- c(plist, s2)
        }
        processData(x, plist)
    }
)


# ** adjust ####

# *** limma ####

setMethod("processData",
    signature(x = "allMatrix", param = "limmaAdjustParam"),
    function(x, param, context = NULL, ...) {
        package_check("limma")
        if (is.null(context)) {
            c(
                "limma adjustment: `context` arg should be a gobject",
                "containing the columns to use for batches and/or covariates",
                "information."
            ) %>%
                wrap_txt(errWidth = TRUE) %>%
                stop(call. = FALSE)
        }
        batches <- param$batch_columns
        covariates <- param$covariate_columns
        if (is.null(batches) && is.null(covariates)) {
            "limma adjustment: At least one of `batch_columns` or 
            `covariate_columns` must be provided." %>%
                wrap_txt() %>%
                stop(call. = FALSE)
        }

        sample_order <- colnames(x)
        limma_args <- list(x = x, ...)
        # batches
        if (!is.null(batches)) {
            b_dt <- .get_svkey(batches, context, sample_order = sample_order)
            if (ncol(b_dt > 2)) {
                "max of 2 columns are allowed for 'batch_columns'" %>%
                    stop(call. = FALSE)
            } else {
                limma_args$batch <- b_dt[[1]]
                if (ncol(b_dt == 2)) {
                    limma_args$batch2 <- b_dt[[2]]
                }
            }
        }
        # covariates
        if (!is.null(covariates)) {
            c_dt <- .get_svkey(
                covariates, context, sample_order = sample_order)
            limma_args$covariates <- as.matrix(c_dt)
        }
        do.call(limma::removeBatchEffect, args = limma_args) %>%
            as("Matrix")
    })


#' @rdname processExpression
#' @export
processExpression <- function(gobject, param,
    name = NULL,
    expression_values = "raw",
    spat_unit = NULL, 
    feat_type = NULL, 
    return_gobject = TRUE,
    ...) {
    ex <- getExpression(gobject,
        values = expression_values,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "exprObj",
        set_defaults = TRUE
    )
    process_args <- list(
        x = ex,
        param = param,
        ...
    )
    if (!is.null(name)) {
        process_args$name <- name
    }

    # detect svkeys
    if (!is.list(param)) param <- list(param)
    param_dump <- lapply(param, function(p) {
        p@param
    })
    has_svk <- .check_svkey(unlist(param_dump), type = "any")
    
    if (has_svk) process_args$context <- gobject
    
    res <- do.call(processData, args = process_args)
    if(!isTRUE(return_gobject)) return(res)
    setGiotto(gobject, res)
}





#' @title normalizeGiotto
#' @name normalizeGiotto
#' @description fast normalize and/or scale expression values of Giotto object
#' @param gobject `giotto` object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param norm_methods normalization method to use
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param log_norm transform values to log-scale
#' @param log_offset offset value to add to expression matrix, default = 1
#' @param logbase log base to use to log normalize expression values
#' @param scale_feats z-score genes over all cells
#' @param scale_cells z-score cells over all genes
#' @param scale_order order to scale feats and cells
#' @param theta theta parameter for the pearson residual normalization step
#' @param name character. name to use for normalization results
#' @param verbose be verbose
#' @param scale_genes deprecated, use scale_feats
#' @param update_slot deprecated. Use `name` param instead
#' @md
#' @returns `giotto` object
#' @details Currently there are two 'methods' to normalize your raw counts data.
#'
#' A. The standard method follows the standard protocol which can be adjusted
#' using the provided parameters and follows the following order: \cr
#' \itemize{
#'   \item{1. Data normalization for total library size and scaling by a custom
#'   scale-factor.}
#'   \item{2. Log transformation of data.}
#'   \item{3. Z-scoring of data by genes and/or cells.}
#' }
#' B. The normalization method as provided by the osmFISH paper is also
#' implemented: \cr
#' \itemize{
#'   \item{1. First normalize genes, for each gene divide the counts by the
#'   total gene count and multiply by the total number of genes.}
#'   \item{2. Next normalize cells, for each cell divide the normalized gene
#'   counts by the total counts per cell and multiply by the total number of
#'   cells.}
#' }
#' C. The normalization method as provided by Lause/Kobak et al is also
#' implemented: \cr
#' \itemize{
#'   \item{1. First calculate expected values based on Pearson correlations.}
#'   \item{2. Next calculate z-scores based on observed and expected values.}
#' }
#' D. Quantile normalization across features
#' \itemize{
#'   \item{1. Rank feature expression}
#'   \item{2. Define a common distribution by sorting expression values per
#'   feature then finding the mean across all features per index}
#'   \item{3. Apply common distribution to expression information by using
#'   the ranks from step 1 as indices}
#' }
#' By default the latter two results will be saved in the Giotto slot for
#' scaled expression, this can be changed by changing the update_slot parameters
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' normalizeGiotto(g) # default is method A
#' @export
normalizeGiotto <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = "raw",
    norm_methods = c("standard", "pearson_resid", "osmFISH", "quantile"),
    library_size_norm = TRUE,
    scalefactor = 6e3,
    log_norm = TRUE,
    log_offset = 1,
    logbase = 2,
    scale_feats = TRUE,
    scale_genes = deprecated(),
    scale_cells = TRUE,
    scale_order = c("first_feats", "first_cells"),
    theta = 100,
    name = "scaled",
    update_slot = deprecated(),
    verbose = TRUE) {
    ## deprecated arguments
    scale_feats <- deprecate_param(
        scale_genes, scale_feats,
        fun = "normalizeGiotto",
        when = "3.0.0"
    )
    name <- deprecate_param(
        update_slot, name,
        fun = "normalizeGiotto",
        when = "4.1.3"
    )

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    ## default is to start from raw data
    values <- match.arg(expression_values, unique(c("raw", expression_values)))
    raw_expr <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )

    norm_methods <- match.arg(
        arg = norm_methods, choices = c(
            "standard", "pearson_resid", "osmFISH", "quantile"
        )
    )

    # normalization according to standard methods
    gobject <- switch(norm_methods,
        "standard" = .rna_standard_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            library_size_norm = library_size_norm,
            scalefactor = scalefactor,
            log_norm = log_norm,
            log_offset = log_offset,
            logbase = logbase,
            scale_feats = scale_feats,
            scale_cells = scale_cells,
            scale_order = scale_order,
            verbose = verbose
        ),
        "osmFISH" = .rna_osmfish_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            name = name,
            verbose = verbose
        ),
        "pearson_resid" = .rna_pears_resid_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            theta = theta,
            name = name,
            verbose = verbose
        ),
        "quantile" = .quantile_norm(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            name = name,
            verbose = verbose
        )
    )

    ## update parameters used ##

    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try(
        {
            upstream_func <- sys.call(-2)
            fname <- as.character(upstream_func[[1]])
            if (fname == "processGiotto") {
                return(gobject)
            }
        },
        silent = TRUE
    )


    # If this function call is not downstream of processGiotto, update normally
    gobject <- update_giotto_params(gobject, description = "_normalize")

    return(gobject)
}















# internals ####

# * params setup ####
.norm_param_arcsinh <- function(...) {
    p <- new("arcsinhNormParam", ...)
    p$c <- p$c %null% 5
    p
}

.norm_param_lib <- function(...) {
    p <- new("libraryNormParam", param = list(...))
    p$scalefactor <- p$scalefactor %null% 6e3
    p
}
.norm_param_log <- function(...) {
    p <- new("logNormParam", param = list(...))
    p$base <- p$base %null% 2
    p$offset <- p$offset %null% 1
    p
}
.norm_param_osmfish <- function(...) {
    new("osmFISHNormParam", param  = list(...))
}
.norm_param_pears_resid <- function(...) {
    p <- new("pearsonResidNormParam", param = list(...))
    p$theta <- p$theta %null% 100
    p
}
.norm_param_quantile <- function(...) {
    new("quantileNormParam", param = list(...))
}
.norm_param_default <- function(...) {
    p <- new("defaultNormParam", param = list(...))
    p$library_size_norm <- p$library_size_norm %null% TRUE
    p$scalefactor <- p$scalefactor %null% 6e3
    p$log_norm <- p$log_norm %null% TRUE
    p$log_offset <- p$log_offset %null% 1
    p$logbase <- p$logbase %null% 2
    p
}
.norm_param_l2 <- function(...) {
    new("l2NormParam", param = list(...))
}
.norm_param_tfidf <- function(...) {
    p <- new("tfidfNormParam", param = list(...))
    p$scalefactor <- p$scalefactor %null% 1e4
    # `sub_method` defaults happen at time of `processData()` call
    # allows opportunity for the user to alter `sub_method` setting via `$<-`
    
    # `tf_fun`, `log_scale_fun`, `mat_log_fun` params are assigned in specific
    # S4 method calls
    p
}
.norm_param_tfidf_method_defaults <- function(p) {
    checkmate::assert_class(p, "tfidfNormParam")
    p$sub_method <- p$sub_method %null% "default"
    if (!identical(p$sub_method, "default") && !p$sub_method %in% seq_len(3)) {
        stop("[tfidfNormParam] `$sub_method` must be either \"default\" or one of 1, 2, or 3 ",
             call. = FALSE)
    } 
    switch(p$sub_method,
        "1" = {
            # Casanovich et al.
            p$log_tf <- p$log_tf %null% FALSE
            p$log_idf <- p$log_idf %null% TRUE
            p$log_tf_idf <- p$log_tf_idf %null% FALSE
            p$binarized_rowsums <- p$binarized_rowsums %null% FALSE
        },
        "2" = {
            # Stuart et al.
            p$log_tf <- p$log_tf %null% FALSE
            p$log_idf <- p$log_idf %null% FALSE
            p$log_tf_idf <- p$log_tf_idf %null% TRUE
            p$binarized_rowsums <- p$binarized_rowsums %null% FALSE
        },
        "3" = {
            # ArchR method 3
            p$log_tf <- p$log_tf %null% TRUE
            p$log_idf <- p$log_idf %null% TRUE
            p$log_tf_idf <- p$log_tf_idf %null% FALSE
            p$binarized_rowsums <- p$binarized_rowsums %null% FALSE
        },
        {
            # Giotto default
            p$log_tf <- p$log_tf %null% FALSE
            p$log_idf <- p$log_idf %null% TRUE
            p$log_tf_idf <- p$log_tf_idf %null% FALSE
            p$binarized_rowsums <- p$binarized_rowsums %null% TRUE
        }
    )
    p
}

.thresh_param_binarize <- function(...) {
    p <- new("binarizeThreshParam", param = list(...))
    p$threshold <- p$threshold %null% NULL
    p
}
.thresh_param_minmax <- function(...) {
    p <- new("minmaxThreshParam", param = list(...))
    p$lower <- p$lower %null% -Inf
    p$upper <- p$upper %null% Inf
    p
}

.scale_param_zscore <- function(...) {
    p <- new("zscoreScaleParam", param = list(...))
    p$scale <- p$scale %null% TRUE
    p$center <- p$center %null% TRUE
    p$MARGIN <- p$MARGIN %null% 2
    p
}
.scale_param_default <- function(...) {
    p <- new("defaultScaleParam", param = list(...))
    p$scale_feats <- p$scale_feats %null% TRUE
    p$scale_cells <- p$scale_cells %null% TRUE
    p$scale_order <- p$scale_order %null% c("first_feats", "first_cells")
    p$verbose <- p$verbose %null% TRUE
    p
}


.adjust_param_limma <- function(...) {
    p <- new("limmaAdjustParam", param = list(...))
    p@param <- if (is.null(p@param$batch_columns)) {
        c(p@param, list(batch_columns = NULL))
    }
    p@param <- if (is.null(p@param$covariate_columns)) {
        c(p@param, list(covariate_columns = NULL))
    }
    p
}


# * implementations ####

.check_svkey <- function(x, type = c("all", "any")) {
    type <- match.arg(type, choices = c("all", "any"))
    if (!inherits(x, "list")) x <- list(x)
    res <- vapply(x, FUN = inherits, FUN.VALUE = logical(1L), "svkey")
    switch (type,
        "any" = any(res),
        "all" = all(res)
    )
}

# get from gobject and ensure order is correct.
# return without cell_IDs col
.get_svkey <- function(x, gobject, sample_order = NULL) {
    if (!inherits(x, "list")) x <- list(x)
    reslist <- lapply(x, function(key) {
        data <- key@get(gobject)
        if (!is.null(sample_order)) {
            data <- data[match(cell_ID, sample_order)]
        }
        return(data[, -"cell_ID"])
    })
    Reduce(cbind, reslist)
}

.arcsinh_norm <- function(x, c) {
    asinh(x / c)
}

.l2_norm <- function(x) {
    # Calculate column norms (Euclidean length of each column)
    col_norms <- sqrt(colSums_flex(x^2))
    # Avoid division by zero
    col_norms[col_norms == 0] <- 1
    # Normalize each column
    t_flex(t_flex(x) / col_norms)
}

.pears_resid_citation <- function(verbose = NULL) {
    vmsg(.v = verbose, "using 'Lause/Kobak' method to normalize count matrix.
    If used in published research, please cite:
    Jan Lause, Philipp Berens, Dmitry Kobak (2020).
    'Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data'")
}

.libzero_warn <- function(libsizes) {
    if (any(0 == libsizes)) {
        warning(wrap_txt("Total library size or counts for individual spat
            units are 0.
            This will likely result in normalization problems.
            filter (filterGiotto) or impute (imputeGiotto) spatial
            units.")
        )
    }
}

# equivalent to calculating x / colSums(x) for dgCMatrix
# returns just the sparse vector of values that need to be assigned to `@x`
.dgc_div_csum_sparse_vector <- function(x, colsums = NULL) {
    if (is.null(colsums)) colsums <- colSums_flex(x)
    x@x / rep.int(colsums, Matrix::diff(x@p))
}

#' @title Normalize expression matrix for library size
#' @param mymatrix matrix object
#' @param scalefactor scalefactor
#' @returns matrix
#' @keywords internal
#' @noRd
.lib_norm_giotto <- function(mymatrix, scalefactor) {
    libsizes <- colSums_flex(mymatrix)
    .libzero_warn(libsizes = libsizes)

    if (inherits(mymatrix, "dgCMatrix")) {
        norm_expr <- mymatrix
        norm_expr@x <- .dgc_div_csum_sparse_vector(norm_expr, 
            colsums = libsizes) *
            scalefactor
    } else {
        norm_expr <- t_flex(t_flex(mymatrix) / libsizes) * scalefactor
    }
    
    return(norm_expr)
}

#' @title Check if matrix is sparse-like
#' @param x matrix object
#' @returns logical
#' @keywords internal
#' @noRd
# Cannot create class union because DelayedArray does not have sparse class
.is_sparse_like <- function(x) {
    inherits(x, c("sparseMatrix", "dbSparseMatrix")) ||
        (inherits(x, "DelayedArray") && DelayedArray::is_sparse(x))
}

#' @title Log normalize expression matrix
#' @returns matrix
#' @keywords internal
#' @noRd
.log_norm_giotto <- function(mymatrix, base, offset) {
    if (.is_sparse_like(mymatrix) && !isTRUE(offset == 1)) {
        stop("`offset != 1` is not supported for sparse matrices ")
    }

    if (.is_sparse_like(mymatrix)) {
        mymatrix <- log1p(mymatrix) / log(base) # sparsity preserved
    } else {
        mymatrix <- log(as.matrix(mymatrix) + offset) / log(base)
    }
    return(mymatrix)
}


#' @title compute_dbMatrix
#' @description saves dbMatrix to db if global option is set
#' @details
#' Set \code{options(giotto.dbmatrix_compute = FALSE)} if saving dbMatrix
#' after each step of normalization workflow is not desired.
#' @keywords internal
#' @noRd
.compute_dbMatrix <- function(dbMatrix, name, verbose = TRUE) {
    # input validation
    if (!inherits(dbMatrix, "dbMatrix")) {
        stop("dbMatrix must be of class dbMatrix")
    }

    if (!is.character(name)) {
        stop("name must be a character")
    }

    # TODO: update with dbData generic
    con <- dbMatrix:::get_con(dbMatrix)

    # overwrite table by default
    if (name %in% DBI::dbListTables(con)) {
        DBI::dbRemoveTable(con, name)
    }

    if (verbose) {
        msg <- glue::glue("Computing {name} expression matrix on disk...")
        cat(msg)
    }

    dbMatrix[] |>
        dplyr::compute(temporary = FALSE, name = name)

    # TODO: update below with proper setters from dbMatrix
    dbMatrix[] <- dplyr::tbl(con, name) # reassign to computed mat
    dbMatrix@name <- name

    if (verbose) cat("done \n")

    return(dbMatrix)
}

#' @title RNA standard normalization
#' @name .rna_standard_normalization
#' @description standard function for RNA normalization
#' @returns giotto object
#' @keywords internal
#' @noRd
.rna_standard_normalization <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    library_size_norm = TRUE,
    scalefactor = 6e3,
    log_norm = TRUE,
    log_offset = 1,
    logbase = 2,
    scale_feats = TRUE,
    scale_cells = TRUE,
    scale_order = c("first_feats", "first_cells"),
    verbose = TRUE) {
    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: Standard normalization was developed for RNA data \n")
    }

    # evaluate provenance before modifying raw_expr in case h5_file exists
    if (isS4(raw_expr)) {
        provenance <- raw_expr@provenance
    } else {
        provenance <- NULL
    }


    feat_names <- rownames(raw_expr[])
    col_names <- colnames(raw_expr[])

    ## 1. library size normalize
    if (isTRUE(library_size_norm)) {
        norm_expr <- .lib_norm_giotto(
            mymatrix = raw_expr[],
            scalefactor = scalefactor
        )
    } else {
        norm_expr <- raw_expr[]
    }

    ## 2. log normalize
    if (isTRUE(log_norm)) {
        norm_expr <- .log_norm_giotto(
            mymatrix = norm_expr,
            base = logbase,
            offset = log_offset
        )
    }

    ## 3. scale
    if (isTRUE(scale_feats) && isTRUE(scale_cells)) {
        scale_order <- match.arg(
            arg = scale_order, choices = c("first_feats", "first_cells")
        )

        if (scale_order == "first_feats") {
            if (isTRUE(verbose)) {
                vmsg(.v = verbose, "first scale feats and then cells")
            }

            norm_scaled_expr <- t_flex(standardise_flex(
                x = t_flex(norm_expr), center = TRUE, scale = TRUE
            ))
            norm_scaled_expr <- standardise_flex(
                x = norm_scaled_expr, center = TRUE, scale = TRUE
            )
        } else if (scale_order == "first_cells") {
            if (isTRUE(verbose)) {
                vmsg(.v = verbose, "first scale cells and then feats")
            }

            norm_scaled_expr <- standardise_flex(
                x = norm_expr, center = TRUE, scale = TRUE
            )
            norm_scaled_expr <- t_flex(standardise_flex(
                x = t_flex(norm_scaled_expr), center = TRUE, scale = TRUE
            ))
        } else {
            stop("\n scale order must be given \n")
        }
    } else if (isTRUE(scale_feats)) {
        norm_scaled_expr <- t_flex(standardise_flex(
            x = t_flex(norm_expr), center = TRUE, scale = TRUE
        ))
    } else if (isTRUE(scale_cells)) {
        norm_scaled_expr <- standardise_flex(
            x = norm_expr, center = TRUE, scale = TRUE
        )
    } else {
        norm_scaled_expr <- NULL
    }


    ## 4. add cell and gene names back
    if (!is.null(norm_expr)) {
        rownames(norm_expr) <- feat_names
        colnames(norm_expr) <- col_names
    }
    if (!is.null(norm_scaled_expr)) {
        rownames(norm_scaled_expr) <- feat_names
        colnames(norm_scaled_expr) <- col_names
    }

    ## 5. create and set exprObj
    # Save dbMatrix to db
    compute_mat <- getOption("giotto.dbmatrix_compute", default = FALSE)
    if (compute_mat && !is.null(norm_expr) && inherits(norm_expr, "dbMatrix")) {
        norm_expr <- .compute_dbMatrix(
            dbMatrix = norm_expr,
            name = "normalized",
            verbose = verbose
        )
    }

    norm_expr <- create_expr_obj(
        name = "normalized",
        exprMat = norm_expr,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        misc = NULL
    )

    # Save dbMatrix to db
    if (compute_mat && !is.null(norm_scaled_expr) && 
            inherits(norm_scaled_expr, "dbMatrix")) {
        norm_scaled_expr <- .compute_dbMatrix(
            dbMatrix = norm_scaled_expr,
            name = "scaled",
            verbose = verbose
        )
    }

    norm_scaled_expr <- create_expr_obj(
        name = "scaled",
        exprMat = norm_scaled_expr,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        misc = NULL
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(
        gobject, norm_expr, verbose = verbose, initialize = FALSE)
    gobject <- setGiotto(
        gobject, norm_scaled_expr, verbose = verbose, initialize = FALSE)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    ## 6. return Giotto object
    return(initialize(gobject))
}



#' @title RNA osmfish normalization
#' @name .rna_osmfish_normalization
#' @description function for RNA normalization according to osmFISH paper
#' @returns giotto object
#' @keywords internal
#' @noRd
.rna_osmfish_normalization <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    name = "custom",
    verbose = TRUE) {
    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: osmFISH normalization was developed for RNA in situ
                data \n")
    }

    # 1. normalize per gene with scale-factor equal to number of genes
    norm_feats <- (raw_expr[] / rowSums_flex(raw_expr[])) * nrow(raw_expr[])
    # 2. normalize per cells with scale-factor equal to number of cells
    norm_feats_cells <- t_flex((t_flex(norm_feats) /
        colSums_flex(norm_feats)) * ncol(raw_expr[]))

    norm_feats_cells <- createExprObj(
        expression_data = norm_feats_cells,
        name = name,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = raw_expr@provenance
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(gobject, norm_feats_cells, verbose = verbose)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
}


#' @title RNA pearson residuals normalization
#' @name rna_pears_resid_normalization
#' @description function for RNA normalization according to Lause/Kobak et al
#' paper
#' Adapted from https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c#file-pearson_residuals-r
#' @returns giotto object
#' @keywords internal
#' @noRd
.rna_pears_resid_normalization <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    theta = 100,
    name = "scaled",
    verbose = TRUE) {
    # print message with information #
    if (verbose) {
        message("using 'Lause/Kobak' method to normalize count matrix If used in
        published research, please cite:
        Jan Lause, Philipp Berens, Dmitry Kobak (2020).
        'Analytic Pearson residuals for normalization of single-cell RNA-seq UMI
        data' ")
    }

    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: pearson residual normalization was developed for RNA
                count normalization \n")
    }

    if (methods::is(raw_expr[], "HDF5Matrix")) {
        .csums <- .csum_nodrop.HDF5Matrix
        .rsums <- .rsum_nodrop.HDF5Matrix
    } else {
        .csums <- .csum_nodrop.Matrix
        .rsums <- .rsum_nodrop.Matrix
    }

    z <- .prnorm(x = raw_expr[], theta, .csums = .csums, .rsums = .rsums)
    z <- create_expr_obj(
        name = name,
        exprMat = z,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = prov(raw_expr)
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(gobject, z, verbose = verbose)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
}

.quantile_norm <- function(gobject,
    raw_expr,
    feat_type,
    spat_unit,
    name = "quantile",
    verbose = TRUE) {
    z <- .qnorm(x = raw_expr[])
    z <- create_expr_obj(
        name = name,
        exprMat = z,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = prov(raw_expr)
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(gobject, z, verbose = verbose)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
}

# pearson residuals normalization
# x      : raw expression matrix
# .csums : function for colSums that does not drop to vector
# .rsums : function for rowSums that does not drop to vector
.prnorm <- function(x,
    theta = 100,
    .csums = .csum_nodrop.Matrix,
    .rsums = .rsum_nodrop.Matrix) {
    # find 1. colsums, 2. rowsums, 3. matrix sum
    counts_sum0 <- .csums(x)
    counts_sum1 <- .rsums(x)
    counts_sum <- sum(x)

    # get residuals
    mu <- (counts_sum1 %*% counts_sum0) / counts_sum
    z <- (x - mu) / sqrt(mu + mu^2 / theta)

    # clip to be within the range [-sqrt(n), sqrt(n)]
    # This is done to prevent extreme values from dominating the analysis.
    n <- ncol(x)
    z[z > sqrt(n)] <- sqrt(n)
    z[z < -sqrt(n)] <- -sqrt(n)
    return(z)
}



# quantile normalization
.qnorm <- function(x) {
    # apply on features by default
    x <- t_flex(x)
    # Rank the values within each column
    ranked_data <- t_flex(MatrixGenerics::colRanks(x, ties.method = "average"))

    # Calculate the mean of sorted values across all columns
    rank_means <- rowMeans(apply(x, 2, sort))

    # Replace the original values with the rank means
    # TODO revisit for large matrices
    normalized_data <- apply(ranked_data, 2, function(idx) {
        .qnorm_vector(idx, rank_means)
    }) |>
        methods::as("Matrix")

    # Retain the original column names
    colnames(normalized_data) <- colnames(x)
    normalized_data <- t_flex(normalized_data)
    return(normalized_data)
}

# create lookup value vector for quantile norm.
# .5 indices should pull the mean of the adjacent values
# indices: index values with some values being .5, designating ranking ties
# values: values to pull from with the indices
.qnorm_vector <- function(indices, values) {
    sorted_values <- sort(values)
    lower_indices <- floor(indices)
    upper_indices <- ceiling(indices)
    lower_values <- sorted_values[lower_indices]
    upper_values <- sorted_values[upper_indices]
    weights <- indices - lower_indices
    result <- (1 - weights) * lower_values + weights * upper_values
    return(result)
}
# based on ArchR implementations
.tf_idf_norm <- function(x,
    log_tf = FALSE,
    log_idf = TRUE,
    log_tf_idf = FALSE, 
    binarized_rowsums = TRUE,
    scalefactor = 1e4,
    tf_fun = function(mat) {
        mat / colSums_flex(mat)
    },
    mat_log_fun = NULL, # default. uses .logvals() below
    log_scale_fun = function(mat, scalef) {
        log(mat * scalef + 1)
    },
    ...
) {
    # for optimizations on TF calc.
    checkmate::assert_function(tf_fun)
    # for optimizations on logging scaled values
    checkmate::assert_function(log_scale_fun)
    # for optimizations on logging of TF values
    checkmate::assert_function(mat_log_fun, null.ok = TRUE)
    .logvals <- function(vals) {
        log(1 + vals)
    }
    if (is.null(mat_log_fun)) mat_log_fun <- .logvals
    
    # rowsums calc
    if (isTRUE(binarized_rowsums)) {
        # +1 since this will be the denominator
        rsums <- 1 + rowSums_flex(x > 0)
    } else {
        rsums <- rowSums_flex(x)
    }

    # tf calc
    tf <- tf_fun(x)
    if (isTRUE(log_tf)) tf <- mat_log_fun(tf)

    # idf calc
    idf <- ncol(x) / rsums
    if (isTRUE(log_idf)) idf <- .logvals(idf)

    tf_idf <- tf * idf
    if (isTRUE(log_tf_idf)) {
        tf_idf <- log_scale_fun(tf_idf, scalef = scalefactor)
    }
    tf_idf
}

# find optimal otsu threshold
# x should be a SpatRaster
# nbins should usuall be the bitdepth
.otsu <- function(x, nbins = NULL) {
    checkmate::assert_numeric(nbins, null.ok = TRUE)
    vals <- terra::spatSample(x, 
        method = "regular", size = 5e5, na.rm = TRUE, as.df = TRUE
    )[[1]]
    if (is.null(nbins)) {
        nbins <- .bitdepth(vals = vals, return_max = TRUE)
    }
    nbins <- as.integer(nbins)
    
    # create histogram
    h <- hist(vals, breaks = seq(min(vals), max(vals), length.out = nbins + 1),
              plot = FALSE)
    counts <- h$counts
    breaks <- h$breaks

    # Calculate the normalized histogram (probability distribution)
    p <- counts / sum(counts)
    
    # initialize vars
    max_variance <- 0
    optimal_idx <- 0
    # Cumulative sum of probabilities
    cum_sum <- cumsum(p)
    # Cumulative mean
    cum_mean <- cumsum(p * seq_along(p) * (breaks[2] - breaks[1])) / cum_sum
    # Global mean
    global_mean <- sum(p * seq_along(p) * (breaks[2] - breaks[1]))
    
    # For each possible threshold, calculate between-class variance
    for (i in 1:(length(p) - 1)) {
        # Weight of background class
        w0 <- cum_sum[i]
        # Weight of foreground class
        w1 <- 1 - w0
        # If one of the classes is empty, skip this threshold
        if (w0 <= 0 || w1 <= 0) {
            next
        }
        # Mean of background class
        mean0 <- cum_mean[i]
        # Mean of foreground class
        mean1 <- (global_mean - w0 * mean0) / w1
        # Calculate between-class variance
        between_variance <- w0 * w1 * (mean0 - mean1)^2
        # Update optimal threshold if current variance is greater
        if (between_variance > max_variance) {
            max_variance <- between_variance
            optimal_idx <- i
        }
    }
    breaks[optimal_idx + 1]
}

# detect image bitdepth
# x should be a SpatRaster if vals are still needed
.bitdepth <- function(x, vals = NULL, return_max = FALSE) {
    checkmate::assert_numeric(vals, null.ok = TRUE)
    if (is.null(vals)) {
        vals <- terra::spatSample(x, 
            method = "regular", size = 5e5, na.rm = TRUE, as.df = TRUE
        )[[1]]
    }
    res <- ceiling(log(max(vals), base = 2L)) # power of 2 needed to represent
    # value(s)
    res <- 2^ceiling(log(res, base = 2L)) # actual bitdepth
    
    if (isTRUE(return_max)) {
        res <- 2^res - 1
    }
    res
}

.csum_nodrop.Matrix <- function(x) {
    x |>
        Matrix::colSums() |>
        matrix(nrow = 1L) |>
        methods::as("Matrix")
}
.rsum_nodrop.Matrix <- function(x) {
    x |>
        Matrix::rowSums() |>
        matrix(ncol = 1L) |>
        methods::as("Matrix")
}
.csum_nodrop.HDF5Matrix <- function(x) {
    x |>
        MatrixGenerics::colSums2() |>
        matrix(nrow = 1L) |>
        methods::as("HDF5Matrix")
}
.rsum_nodrop.HDF5Matrix <- function(x) {
    x |>
        MatrixGenerics::rowSums2() |>
        matrix(ncol = 1L) |>
        methods::as("HDF5Matrix")
}

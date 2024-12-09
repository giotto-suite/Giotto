
# Giotto 4.1.6 (2024/12/09)

## Bug fixes
* `doScrubletDetect()` seed setting

## Enhancements
* `labelTransfer()` now has `integration_method = "harmony"` for label transferring with an integration pipeline. See ?labelTransfer and the `integration_method` section.
* `importXenium()` `load_transcripts()` can now return a `data.table` rather than the `giottoPoints` representation

## New
* `doMesmerSegmentation()` and `doStardistSegmentation()` segmentation wrappers
* `.varexp()` internal for calculating SVD variance determined with support for partial SVDs
* `.cumvar()` internal for calculating cumulative variance explained
* re-export of `dotPlot()` from GiottoVisuals

## Changes
* GiottoClass req raised to 0.4.5
* GiottoUtils req raised to 0.2.2
* GiottoVisuals req raised to 0.2.10


# Giotto 4.1.5 (2024/11/08)

## Enhancements
* `createGiottoXeniumObject()` auto loading for morphology focus images, image directory loading, auto centroid calculation, allow skipping transcript loading

## Website changes
* New Analysis mini tutorials for showing common processing functions independently of the spatial technology.
* New Slide-seq and OpenST examples.
* New Contributing tab with guidelines for contributing to the package and the website.
* New Visualizations tutorials.
* New Giotto workflow and Core Functions tutorials under Get Started tab.
* New Create and change Giotto instructions tutorial.
* New Spatial Patterns tutorials section.
* New tutorials under Interactivity for regions selection with vitessceR.
* New Multi-samples tutorials section.
* Updated technologies examples.
* Updated tutorials for using Docker and Singularity Giotto containers.
* Homogenized variable names across examples and tutorials.

# Giotto 4.1.4 (2024/10/30)

## Changes
* `createGiottoVisiumObject()` apply a modifier of 0.8461538 to visium spot diameter to reflect actual spot size
* `doLeidenClusterIgraph()` deprecate param `resolution_parameter` in favor of `resolution`

## Enhancements
* `createGiottoVisiumObject()` append multiplicative scalefactor to get micron 
  values from the current coordinate units during Visium object creation. 
  Accessible through `instructions(gobject, "micron_scale")`

# Giotto 4.1.3 (2024/10/27)

## New
* Add `giottoToAnndataZarr()` to create a local anndata zarr folder and interact with the vitessceR package.
* `reduceDims()` API function for dimension reductions
* `runNMF()` implementation that works via RcppML

## Changes
* `runWNN()` and `runIntegratedUMAP()` arguments were updated to make the function flexible to handle any number of modalities.
* update `jackstrawPlot()` to make more flexible and efficient. Changed default params for `scaling`, `centering`, and `feats_to_use` to match `runPCA()`
* change warning when reduction "feats" is selected in `runtSNE()` to error to avoid accidentally wiping the `giotto` object.


# Giotto 4.1.2

## Breaking changes
* remove deprecated `PAGEEnrich()`. Use `runPAGEEnrich()` instead
* remove deprecated `rankEnrich()`. Use `runRankEnrich()` instead
* remove deprecated `hyperGeometricEnrich()`. Use `runHyperGeometricEnrich()` instead
* remove deprecated `createSpatialEnrich()`. Use `runSpatialEnrich()` instead
* remove deprecated `heatmSpatialCorGenes()`. Use `heatmSpatialCorFeats()` instead
* remove deprecated `runPAGEEnrich_OLD()`. Use `runPAGEEnrich()` instead
* remove `do_pca`, `expression_values`, `feats_to_use` args from `runGiottoHarmony()`. Running PCA during the `harmony::RunHarmony()` call is deprecated.

## Enhancements
* add `'quantile'` normalization method to `normalizeGiotto()`

## Changes
* `limma`, `plotly`, and `Rtsne` moved to Suggests
* move `progressr` and `jsonlite` dependencies to GiottoUtils v0.1.12
* remove `reshape2` dependency.

## Bug fixes
* `processGiotto()` can now skip adjust step by default

## New
* `identifyTMAcores()` for assigning IDs to tissue microarray spatial data.
* `labelTransfer()` for transferring labels between giotto objects or subsets thereof. Supercedes `doClusterProjection()`

# Giotto 4.1.1

## Bug fixes
* Allow `giottoInstructions` passing for Xenium convenience functions

## Changes
* Deprecate `screePlot()` `name` in favor of `dim_reduction_name` param


# Giotto 4.1.0 (2024/07/31)

## Breaking changes
* Deprecated `detectSpatialCorGenes()` removed. Use `detectSpatialCorFeats()` instead
* Deprecated `findInteractionChangedGenes()` removed. Use `findInteractionChangedFeats()` instead
* Deprecated `findCellProximityGenes()` removed. Use `findInteractionChangedFeats()` instead
* `createGiottoXeniumObject()` has been overhauled and parameters have changed.

## Bug fixes
* Fix error in `plotInteractivePolygons()` when providing a spatial plot with a continuous scale [#964](https://github.com/drieslab/Giotto/issues/964) by jweis3
* Fix error in DWLS `find_dampening_constant()` when `S[subset, ]` produces only 1 gene.
* Fix error in `interpolateFeatures()` where feature names with `-` or starting with numbers did not work
* Add catch in `runPCAprojectionBatch()` for when ncp requested exceeds number of feats used
* Make `spatCellCellcom()` respect `verbose` flag [#949](https://github.com/drieslab/Giotto/issues/949) by rbutleriii

## New
* Dataset affine registration via interactive shiny app and automated SIFT detection
* Cell segmentation via Cellpose
* `read10xAffineImage()` for reading 10x affine transformed images
* Several modular importer and convenience functions
* ONTraC implementation

## Enhancements
* `print()` methods for `icfObject` and `combIcfObject`

## Changes
* require GiottoUtils (>= 0.1.10)
* require GiottoClass (>= 0.3.3)

# Giotto 4.0.8 (2024/05/22)

## Breaking changes
* `crossSectionGenePlot()` removed. Use `crossSectionFeatPlot()` instead
* `crossSectionGenePlot3D()` removed. Use `crossSectionFeatPlot3D()` instead
* `insertCrossSectionGenePlot3D()` removed Use `insertCrossSectionFeatPlot3D()` instead

## Bug fixes
* `binSpect()` param passing error introduced in _v4.0.6_
* updated `viewHMRFresults3D()` and `viewHMRFresults2D()`
* updated `createCrossSections()`, `insertCrossSectionSpatPlot3D()`, `crossSectionPlot()`, `crossSectionFeatPlot3D()`, `insertCrossSectionFeatPlot3D()`, `crossSectionPlot3D()`, `crossSectionFeatPlot()`

## Changes
* GiottoVisuals (>= 0.2.2), GiottoClass (>= 0.3.1), GiottoUtils (>= 0.1.8) are now required.

# Giotto 4.0.6 (2024/05/13)

## Enhancements
* New `interpolateFeature()` for kriging interpolation of values

## Changes
* GiottoVisuals (>= 0.2.0) and GiottoClass (>= 0.3.0) are now required.



# Giotto 4.0.5 (2024/03/12)

## Bug fixes
* Fix Error "cannot coerce class ‘structure("spatLocsObj", package = "Giotto")’ to a data.frame" in `spatialDE()`

## Enhancements
* `readPolygonVizgenParquet()` now has `calc_centroids = TRUE` by default

# Giotto 4.0.4 (2024/02/28)

## Breaking changes
* Remove `do_manual_adj` and image adjustment params from `createGiottoVisiumObject()`
* `createGiottoVisiumObject()` now creates `giottoLargeImage` for spatial images.
* `exprCellCellcom()` deprecated `gene_set_*` params removed

## Bug fixes
* Fix metadata appending/sorting issues introduced by *GiottoClass v0.1.3* (2024/01/12). Affected functions: `addHMRF()`, `addFeatsPerc()`, `doScrubletDetect()`
* `findNetworkNeighbors()` now has default `spatial_network_name` value of `NULL`
* `get10Xmatrix()` now obeys `split_by_type = FALSE`

## Changes
* Deprecate `set.seed` in favor of `seed` param for `binSpect()`
* `binSpect()` now sets a seed by default for reproducibility
* pkgdown files moved to separate [repo](https://github.com/drieslab/Giotto_website)

## Enhancements
* Use `mixedsort()` for unique clusters metadata info
* Remove unnecessary matrix densification and expose `seed` param in `doScrubletDetect()`
* Remove unnecessary matrix densification in `makeSignMatrixRank()`


# Giotto 4.0.3 (2024/02/20)


## Bug fixes
* Remove old argument `type_default = list(pal = c('blue', 'yellow', 'red'))` in plotRankSpatvsExpr()


# Giotto 4.0.2 (2023/12/21)

## Bug Fixes
* fix bug in `doHclust()`

## Changes
* Move *GiottoClass* back to depends to fix access to some generics


# Giotto 4.0.1 (2023/12/16)

## Breaking changes
* Remove `cell_ids` param for `calculateHVF()` in favor of simpler `random_subset`
* Move *GiottoUtils*, *GiottoClass*, and *GiottoVisuals* to imports

## Added
* Add `parse_affine()` for interpreting affine transform matrices
* Add seed setting to `runGiottoHarmony()`
* Add parallelized calculation for `calculateHVF()` when a *future* plan is set

## Changes
* Fix `createGiottoXeniumObject()` feature metadata reading for `.json` file
* Update `runGiottoHarmony()` to call `harmony::RunHarmony()` 
* Update *Matrix* requirement to >= 1.6.2 (a re-install of *irlba* may resolve issues with Matrix incompatibility.)




# Giotto 4.0.0 (2023/11/29)

## Breaking changes
* Update to modular package organization with the main packages being `GiottoUtils`, `GiottoClass`, `GiottoVisuals`, and `Giotto` as the analytical umbrella package.

## Added

* New File `spatial_enrichment_visuals.R`
* New function `findCellTypesFromEnrichment()` within `spatial_enrichment_visuals.R` to show most probable cell types based on a provided enrichment
* New function `plotCellTypesFromEnrichment()` within `spatial_enrichment_visuals.R` that generates a bar plot of cell types vs frequency based on a provided enrichment
* New function `pieCellTypesFromEnrichment()` within `spatial_enrichment_visuals.R` that generates a pie chart of cell types based on a provided enrichment
* New function `addVisiumPolygons()` within `convenience.R` (along with its requisite internal functions) that adds  circular polygons centered at the spatial locations of a Giotto Object made with Visium data. Takes a Giotto Object and a path to the Visium output file `scalefactors_json.json` as input arguments.
* Added `addVisiumPolygons()` to `createGiottoVisiumObject()` workflow.
* Add `cell_ids` param to `calculateHVF()` to allow calculation of HVFs on a subset of cells
* Add seed setting to `runGiottoHarmony()`
* Update `runGiottoHarmony()` to call `harmony::RunHarmony()`


## Changes
* Update *Matrix* requirement to >= 1.6.3



# Giotto 3.3.1 (2023/08/02)

## Breaking changes

* Change `checkGiottoEnvironment()`. Downgrade from error to message and return FALSE when a provided directory does not exist

## Added

* New file `poly_influence.R`
* New function `showPolygonSizeInfluence()` within `poly_influence.R` to show if cells switch clusters when across resized polygon annotations 
* New function `showCellProportionSwitchedPie()` within `poly_influence.R` to visualize results from `showPolygonSizeInfluence()` in a pie chart
* New function `showCellProportionSwitchedSankey()` within `poly_influence.R` to visualize results from `showPolygonSizeInfluence()` in a Sankey diagram
* New function `makePseudoVisium()` within `giotto_structure.R` to generate a pseudo visium grid of circular spots 
* New function `tessellate()` within `giotto_structure.R` to generate a grid of hexagons or squares for spatial binning
* New file `feature_set_enrichment.R`
* New function `doFeatureSetEnrichment()` within `feature_set_enrichment.R` for GSEA analysis
* New function `doGiottoClustree()` within `clustering.R` for visualizations of leiden clusters at varying resolutions
* New `createArchRProj()` and `CreateGiottoObjectFromArchR()` functions to create a `giotto` object with ATAC or epigenetic features using the *ArchR* pipeline.
* New `giottoMasterToSuite()` function to convert a `giotto` object created with the master version to a Giotto suite object.
* New `readPolygonVizgenParquet()` for updated parquet outputs
* Add *checkmate* to Imports for assertions checking
* Add exported `create` function for `exprObj` creation
* New file `spatial_manipulation.R`
* Add `ext()` methods for `giottoPolygon`, `giottoPoints`, `spatialNetworkObj`, `spatLocsObj`, `giottoLargeImage`
* Add `flip()` methods for `giottoPolygon`, `giottoPoints`, `spatialNetworkObj`, `spatLocsObj`, `SpatExtent`, `giottoLargeImage`
* Add access to terra plotting params for `giottoLargeImage` `plot()` method.

## Changes

* Fix bug in `combine_matrices()`
* Fix bug in `createGiottoObject()` that will not allow object creation without supplied expression information
* Updated `polyStamp()` to replace an apply function with a crossjoin for better performance.
* Updated `spatInSituPlotPoints()` with `plot_last` parameter. Default output now plots polygons above points for better visibility.
* Add check for spatLocsObj for spatlocs in polyStamp()
* Removed various print() and cat() statements throughout. 
* Changed default verbose argument to FALSE for createGiottoObject
* Changed default verbose argument to FALSE for joinGiottoObjects
* Changed default verbose argument to FALSE for createGiottoObjectSubcellular
* Default verbose = FALSE argument added to cellProximityEnrichmentSpots
* Default verbose = FALSE argument added to .specific_CCCScores_spots
* Default verbose = FALSE argument added to runWNN
* Default verbose = FALSE argument added to subset_giotto_points_object 
* Default verbose = FALSE argument added to subset_feature_info_data 
* Default verbose = FALSE argument added to subsetGiotto
* Default verbose = FALSE argument added to subsetGiottoLocsSubcellular 
* Default verbose = FALSE argument added to .createGiottoXeniumObject_subcellular
* Update `readPolygonFilesVizgenHDF5()` add option to return as `data.table` and skip `giottoPolygon` creation. Downstream `giottoPolygon` creation refactored as new internal function
* Update cell segmentation workflow to check for *deepcell* and *PIL* python packages
* Update cell segmentation workflow to return grayscale mask images instead of RGB
* Update `createGiottoVisiumObject()` image h5 scalefactors reading to use partial matching for whether hi or lowres image is supplied
* Update `giottoLargeImage` `plot()` method to use `terra::plot()` instead of `terra::plotRGB()` for grayscale images
* Remove unnecessary prints from `subsetGiotto()`
* Fix bug in `readCellMetadata()` and `readFeatMetadata()`


# Giotto 3.3.0 (2023/04/18)

## Breaking changes

* Set Suite as default branch
* Removed all deprecated accessors from `accessors.R`
* `set_default_feat_type()` error downgraded to warning when no `feat_type`s exist for given `spat_unit`
* update `loadGiotto()` and `saveGiotto()` to allow using long strings as column names in the spatVector objects
* `'active_spat_unit'` and `'active_feat_type'` params that can be set through `instructions()` are now used instead of 'giotto.spat_unit' and 'giotto.feat_type' global options
* removed duplicate `create_dimObject()` internal function. Keeping `create_dim_obj()`



## Added

* New implementations of `anndataToGiotto()` and `giottoToAnnData()` for Nearest Neighbor and Spatial Networks
* New `check_py_for_scanpy()` function, shifting code around in `anndataToGiotto()`
* Add `initialize()` method for `giotto`
* Add exported `create` constructor functions for Giotto S4 subobjects
* Add `activeSpatUnit()` and `activeFeatType()` for getting and setting active defaults on gobject
* New `get_*_list()` internal functions for retrieving list of all objects of a particular class for a spat_unit and feat_type
* Add `instructions()` generic for `giotto` to access and edit `giottoInstructions`
* Add `centroids()` method for `giottoPolygon` to get centroid info
* Add `overlaps()` generic for accessing `overlaps` slot of `giottoPolygon`
* Add `[` and `[<*` (empty) access generics to get the data from main slots of `giottoPolygon` and `giottoPoints`
* Add cores detection to run on package attach. (`getOption('giotto.cores')`)
* Add option to return as `giottoPoints` from `getFeatureInfo` (default is still `SpatVector`)
* Add `spatVector_to_dt2` internal as a barebones alternative to `spatVector_to_dt()`
* Add `getRainbowColors()` color palette
* New `assign_objnames_2_list()` and `assign_listnames_2_obj()` internals for passing list names to object `@name` slots and vice versa
* New test_that `test_createObject.R` script for `read` functions/S4 subobject creation
* New test_that `test_accessors.R` script for `accessor` functions
* New test_that `test_gobject.R` script for gobject consistency checks



## Changes

* Update `installGiottoEnvironment()` and downstream internal functions to allow custom python installation with a new argument, `mini_install_path`.
* Update `checkGiottoEnvironment()` to account for custom python installations with a new argument, `mini_install_path`.
* Update `removeGiottoEnvironment()` to account for custom python installations with a new argument, `mini_path`.
* Update `createGiottoObject()` with new data ingestion pipeline
* Modify `cell_ID`, `feat_ID`, `cell_metadata`, `feat_metadata` slot initialization
* Update `read_expression_data()` and `evaluate_expr_matrix()` to be compatible with `exprObj`
* Change `changeGiottoInstructions()` to allow addition of new params and enforce logical class of known params
* Update and fix bugs in `createGiottoCosMxObject()` associated with polygon placement and generation
* Update `plot()` for `giottoPoints` with faster rasterized plotting option. (Now used by default)
* Fix bug in `doLouvainCluster()` (sub)functions and made them compatible with new Giotto Suite framework.
* Fix bug in `gefToGiotto()` bin_size arguments.
* Update `loadGiotto()` and `saveGiotto()` with path.expand to expand provided file/directory paths
* Organize new and refactored slot `check` functions in `giotto.R` for checking gobject consistency during `initialize()`
* Organize new and refactored `evaluate` functions in `data_evaluation.R` for data wrangling of external data
* Organize new and refactored `read` functions in `data_input.R` for ingesting data and converting to list of Giotto native S4 subobjects
* Organize dummy documentation in `dd.R` for inheriting commonly used documentation
* Moved `create_featureNetwork_object()`, `create_giotto_points_object()`, `create_giotto_polygon_object()` to classes.R
* Moved `depth()` from giotto.R to utilities.R


# Giotto 3.2.0 (2023/02/02)

## Breaking changes

* Removed support for deprecated nesting in `@nn_network` slot
* `createSpatialNetwork()` will now output a `spatialNetworkObj` by default when `return_gobject = FALSE`. It is possible to change this back to the data.table output by setting `output = 'data.table'`
* Set incomplete classes in classes.R as virtual to prevent their instantiation
* Removed `createGiottoCosMxObject()` `aggregate` and `all` workflows until they are updated

## Added

* New `gefToGiotto()` interoperability function to convert gef object from Stereo*seq to giotto
* New `giottoToAnnData()` interoperability function to convert giotto object to squidpy flavor anndata .h5ad file(s)
* New `giottoToSpatialExperiment()` and `spatialExperimentToGiotto()` to convert between Giotto and SpatialExperiment
* New `spatialAutoCorLocal()` and `spatialAutoCorGlobal()` functions to find spatial autocorrelations from expression and cell metadata information
* New `createSpatialWeightMatrix()` function to generate spatial weight matrix from spatial networks for autocorrelation
* Add spatial_interaction_spot.R with functions adapted from master branch for working with the Giotto suite object.
* New exported accessors for slots (experimental)
* New `multiomics` slot in `giotto`
* Add `coord_fix_ratio` param to `spatFeatPlot2D()` and `spatFeatPlot2D_single()`
* Add `order` parameter to `dimFeatPlot2D` and `spatDimFeatPlot2d` to plot and order cells according to the levels of the selected feature ([#477](https://github.com/drieslab/Giotto/issues/477))
* Add `plot()` method for `spatialNetworkObj`
* New `set_row_order_dt()` internal for setting `data.table` to a specific row order by reference
* New `fread_colmatch()` internal for fread loading a subset of rows based on matches in a specified column
* Add missing `create_nn_net_obj()` internal constructor function for S4 `nnNetObj`
* Add `id_col`, `x_col`, `y_col` params to `polyStamp()` to make stamp location input more flexible
* Add `optional` and `custom_msg` params to `package_check()`
* New `wrap()` and `vect()` generics for `giotto`, `giottoPoints`, and `giottoPolygons`
* New `rotate()`, `t()`, and `spatShift` generics for giotto subobject spatial manipulation
* New `spatIDs()` and `featIDs()` generics
* New `objName()` and `objName` generics for setting the names of relevant S4 subobjects
* New `rbind()` generic to append `giottoPolygon` objects
* Add packages `exactextractr` and `sf` to "suggests" packages
* Add package `progressr` to "imports" packages

## Changes

* Move giotto object method*specific creation functions from `giotto.R` to `convenience.R`
* Update `addFeatMetadata()` to handle replacement of existing columns
* Update `show()` method for `giotto`
* Update `show()` method for `spatEnrObj`
* Deprecate older snake_case accessors
* Deprecate `polygon_feat_names` param in favor of `z_indices` in `readPolygonFilesVizgenHDF5()`
* Deprecate `xy_translate_spatial_locations()` in favor of `shift_spatial_locations()`
* Optimize `readPolygonFilesVizgen()`
* Fix bug in `replaceGiottoInstructions()` where instructions with more slots than previous are not allowed
* Fix bug in `loadGiotto()` that prevents proper parsing of filenames when spat_unit or feat_type contains '_' characters
* Fix `loadGiotto()` loss of over*allocation for data.tables*based objects after loading from disk


# Giotto 3.1.0 (2022/12/01)  
  

## Added
* New `initialize()` generic that calls `setalloccol()` for data.table*based S4 subobjects to allow setting by reference
* New `spatUnit`, `spatUnit<*`, `featType`, and `featType<*` feat type generics for S4 subobjects for setting the relevant slots
* Add `hexVertices()` to polygon shape array generation functionality

## Changes  

* Update `createGiottoCosMxObject()` for 3.0 and modularization of functions. 'subcellular' workflow has been tested to work along with an updated tutorial.
* Update grid plotting behavior to set a default number columns to use based on number of elements to plot. Can be overridden by explicitly providing input to `cow_n_col` param
* Fix bug in `annotateGiotto()` after 3.0 update ([#433](https://github.com/drieslab/Giotto/issues/433#issuecomment*1324211224))
* Fix bug in `joinGiottoObjects()` metadata processing
* Update seed setting behavior in [dimension_reduction.R](https://github.com/drieslab/Giotto/blob/suite/R/dimension_reduction.R) and [clustering.R](https://github.com/drieslab/Giotto/blob/suite/R/clustering.R)





# Giotto 3.0.1 (2022/11/20)

## Added

* New system color support detection (based on crayon package logic)
* Add ability to turn off colored text in `show` functions with `options("giotto.color_show" = FALSE)`

## Changes

* Fix bug in `extract_polygon_list()` ([#433](https://github.com/drieslab/Giotto/issues/433#issuecomment*1321221382))
* Update Unicode character printing with `show` functions for Latin1 systems





# Giotto 3.0.0 (2022/11/18)

## Breaking changes

* S4 subobjects framework will require giotto objects to be remade

## Added

* New `createGiottoXeniumObject()` for loading 10x Xenium data
* New S4 subobjects. Details can be found in [classes.R](https://github.com/drieslab/Giotto/blob/suite/R/classes.R)
* New basic generics for S4 subobjects. Mainly the use of `[]` and `[]<*` to get or set information into the main data slot
* New `@provenance` slot in S4 subobjects to track provenance of aggregated information (z_layers used for example)
* New `calculateOverlapPolygonImages()` for calculating overlapped intensities from image*based information (e.g. IMC, IF, MIBI, ...) and polygon data (e.g. cell)
* New `overlapImagesToMatrix()` converts intensity*polygon overlap info into an expression matrix (e.g. cell by protein)
* New `aggregateStacks()` set of functions work with multiple subcellular layers when generating aggregated expression matrices

## Changes

* Update `setter` functions to read the `@spat_unit` and `@feat_type` slots of subobjects to determine nesting
* Update of `show` functions to display color coded nesting names and tree structure





# Giotto 2.1.0 (2022/11/09)

## Breaking changes

* Update of python version to **3.10.2** [details](https://giottosuite.readthedocs.io/en/latest/additionalinformation.html#giotto*suite*2*1*0*2202*11*09)

## Added

* New `anndataToGiotto()` to convert scanpy anndata to Giotto






# Giotto 2.0.0.998

## Added

* New `GiottoData` package to work with spatial datasets associated with Giotto
	* Stores the minidatasets: preprocessed giotto objects that are ready to be used in any function
	* Moved: `getSpatialDataset()` and `loadGiottoMini()` functions to this package
* New `saveGiotto()` and `loadGiotto()` for preserving memory*pointer based objects. In [general_help.R](https://github.com/drieslab/Giotto/blob/suite/R/general_help.R)
  * It saves a Giotto object into a folder using a specific structure. Essentially a wrapper around `saveRDS()` that also works with spatVector and spatRaster pointers.
* New `plotInteractivePolygon()` for plot*interactive polygonal selection of points.
* New polygon shape array creation through `polyStamp()`, `circleVertices`, `rectVertices`. In [giotto_structures.R](https://github.com/drieslab/Giotto/blob/suite/R/giotto_structures.R)
* Add accessor functions `get_CellMetadata` (alias of `pDataDT()`), `set_CellMetadata`, `get_FeatMetadata` (alias of `fDataDT()`), `set_FeatMetadata`. See [accessors.R](https://github.com/drieslab/Giotto/blob/suite/R/accessors.R)
* New `filterDistributions()` to generate histogram plots from expression statistics

## Changes

* Deprecate `plotInteractionChangedGenes()` ,`plotICG()`, `plotCPG()` in favor of `plotInteractionChangedFeatures()` and `plotICF()` and `plotCPF()`
* Deprecate `plotCellProximityGenes()`, in favor of `plotCellProximityFeatures()`
* Deprecate `plotCombineInteractionChangedGenes()`, `plotCombineICG()`, `plotCombineCPG()` in favor of `plotCombineInteractionChangedFeatures()` and `plotCombineICF()`
* Deprecate `findInteractionChangedGenes()`, `findICG()`, `findCPG()` in favor of `findInteractionChangedFeats()` and `findICF`
* Deprecate `filterInteractionChangedGenes()`, `filterICG()`, `filterCPG()` in favor of `filterInteractionChangedFeats()` and `filterICF()`
* Deprecate `combineInteractionChangedGenes()`, `combineICG()`, `combineCPG()` in favor of `combineInteractionChangedFeats()` and `combineICF()`
* Deprecate `combineCellProximityGenes_per_interaction()` in favor of `combineCellProximityFeatures_per_interaction()`

## Breaking changes

* ICF output internal object structure names have changed to use feats instead of genes













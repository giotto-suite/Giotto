# GiottoUtils ####
#' @returns pipe operator
#' @export
GiottoUtils::`%>%`
#' @returns vector of colors
#' @export
GiottoUtils::getRainbowColors
#' @returns vector of colors
#' @export
GiottoUtils::getDistinctColors



# GiottoClass ####

# GiottoClass defined generics only
#' @export
GiottoClass::activeSpatUnit
#' @export
GiottoClass::`activeSpatUnit<-`
#' @export
GiottoClass::activeFeatType
#' @export
GiottoClass::`activeFeatType<-`
#' @export
GiottoClass::ext
#' @export
GiottoClass::`ext<-`
#' @export
GiottoClass::featType
#' @export
GiottoClass::`featType<-`
#' @export
GiottoClass::spatUnit
#' @export
GiottoClass::`spatUnit<-`
#' @export
GiottoClass::instructions
#' @export
GiottoClass::`instructions<-`
#' @export
GiottoClass::prov
#' @export
GiottoClass::`prov<-`
#' @export
GiottoClass::objName
#' @export
GiottoClass::`objName<-`
#' @export
GiottoClass::as.points
#' @export
GiottoClass::as.polygons
#' @export
GiottoClass::as.terra
#' @export
GiottoClass::as.stars
#' @export
GiottoClass::as.sf
#' @export
GiottoClass::as.sp
#' @export
GiottoClass::overlaps
#' @export
GiottoClass::calculateOverlap
#' @export
GiottoClass::overlapToMatrix
#' @export
GiottoClass::centroids
#' @export
GiottoClass::copy
#' @export
GiottoClass::crop
#' @export
GiottoClass::density
#' @export
GiottoClass::flip
#' @export
GiottoClass::spin
#' @export
GiottoClass::spatShift
#' @export
GiottoClass::featIDs
#' @export
GiottoClass::spatIDs
#' @export
GiottoClass::setGiotto
#' @export
GiottoClass::rescale
#' @export
GiottoClass::vect
#' @export
GiottoClass::wrap
#' @export
GiottoClass::createGiottoPoints
#' @export
GiottoClass::createGiottoPolygon

# All other exported functions
#' @export
GiottoClass::addCellMetadata
#' @export
GiottoClass::addFeatMetadata
#' @export
GiottoClass::addGiottoImage
#' @export
GiottoClass::addGiottoImageMG
#' @export
GiottoClass::addGiottoLargeImage
#' @export
GiottoClass::addGiottoPoints
#' @export
GiottoClass::addGiottoPoints3D
#' @export
GiottoClass::addGiottoPolygons
#' @export
GiottoClass::addNetworkLayout
#' @export
GiottoClass::addSpatialCentroidLocations
#' @export
GiottoClass::addSpatialCentroidLocationsLayer
#' @export
GiottoClass::aggregateStacks
#' @export
GiottoClass::aggregateStacksExpression
#' @export
GiottoClass::aggregateStacksLocations
#' @export
GiottoClass::aggregateStacksPolygonOverlaps
#' @export
GiottoClass::aggregateStacksPolygons
#' @export
GiottoClass::anndataToGiotto
#' @export
GiottoClass::annotateGiotto
#' @export
GiottoClass::annotateSpatialGrid
#' @export
GiottoClass::annotateSpatialNetwork
#' @export
GiottoClass::calculateMetaTable
#' @export
GiottoClass::calculateMetaTableCells
#' @export
GiottoClass::calculateOverlapParallel
#' @export
GiottoClass::calculateOverlapPolygonImages
#' @export
GiottoClass::calculateOverlapRaster
#' @export
GiottoClass::calculateOverlapSerial
#' @export
GiottoClass::calculateSpatCellMetadataProportions
#' @export
GiottoClass::changeGiottoInstructions
#' @export
GiottoClass::changeImageBg
#' @export
GiottoClass::checkGiottoEnvironment
#' @export
GiottoClass::circleVertices
#' @export
GiottoClass::combineCellData
#' @export
GiottoClass::combineFeatureData
#' @export
GiottoClass::combineFeatureOverlapData
#' @export
GiottoClass::combineMetadata
#' @export
GiottoClass::combineSpatialCellMetadataInfo
#' @export
GiottoClass::combineToMultiPolygon
#' @export
GiottoClass::convertGiottoLargeImageToMG
#' @export
GiottoClass::createBentoAdata
#' @export
GiottoClass::createCellMetaObj
#' @export
GiottoClass::createDimObj
#' @export
GiottoClass::createExprObj
#' @export
GiottoClass::createFeatMetaObj
#' @export
GiottoClass::createGiottoImage
#' @export
GiottoClass::createGiottoInstructions
#' @export
GiottoClass::createGiottoLargeImage
#' @export
GiottoClass::createGiottoLargeImageList
#' @export
GiottoClass::createGiottoObject
#' @export
GiottoClass::createGiottoObjectSubcellular
#' @export
GiottoClass::createGiottoPolygonsFromDfr
#' @export
GiottoClass::createGiottoPolygonsFromGeoJSON
#' @export
GiottoClass::createGiottoPolygonsFromMask
#' @export
GiottoClass::createMetafeats
#' @export
GiottoClass::createNearestNetObj
#' @export
GiottoClass::createSpatEnrObj
#' @export
GiottoClass::createSpatLocsObj
#' @export
GiottoClass::createSpatNetObj
#' @export
GiottoClass::createSpatialDefaultGrid
#' @export
GiottoClass::createSpatialDelaunayNetwork
#' @export
GiottoClass::createSpatialFeaturesKNNnetwork
#' @export
GiottoClass::createSpatialGrid
#' @export
GiottoClass::createSpatialKNNnetwork
#' @export
GiottoClass::createSpatialNetwork
#' @export
GiottoClass::createSpatialWeightMatrix
#' @export
GiottoClass::createNearestNetwork
#' @export
GiottoClass::cropGiottoLargeImage
#' @export
GiottoClass::distGiottoImage
#' @export
GiottoClass::estimateImageBg
#' @export
GiottoClass::fDataDT
#' @export
GiottoClass::featureNetwork
#' @export
GiottoClass::gefToGiotto
#' @export
GiottoClass::getCellMetadata
#' @export
GiottoClass::getDimReduction
#' @export
GiottoClass::getExpression
#' @export
GiottoClass::getFeatureInfo
#' @export
GiottoClass::getFeatureMetadata
#' @export
GiottoClass::getGiottoImage
#' @export
GiottoClass::getMultiomics
#' @export
GiottoClass::getNearestNetwork
#' @export
GiottoClass::getPolygonInfo
#' @export
GiottoClass::getSpatialEnrichment
#' @export
GiottoClass::getSpatialGrid
#' @export
GiottoClass::getSpatialLocations
#' @export
GiottoClass::getSpatialNetwork
#' @export
GiottoClass::getSpatialEnrichment
#' @export
GiottoClass::getSpatialGrid
#' @export
GiottoClass::getSpatialLocations
#' @export
GiottoClass::getSpatialNetwork
#' @export
GiottoClass::giotto
#' @export
GiottoClass::giottoImage
#' @export
GiottoClass::giottoLargeImage
#' @export
GiottoClass::giottoMasterToSuite
#' @export
GiottoClass::giottoPoints
#' @export
GiottoClass::giottoPolygon
#' @export
GiottoClass::giottoToAnnData
#' @export
GiottoClass::giottoToSeurat
#' @export
GiottoClass::giottoToSeuratV4
#' @export
GiottoClass::giottoToSeuratV5
#' @export
GiottoClass::giottoToSpatialExperiment
#' @export
GiottoClass::hexVertices
#' @export
GiottoClass::hist
#' @export
GiottoClass::installGiottoEnvironment
#' @export
GiottoClass::joinGiottoObjects
#' @export
GiottoClass::loadGiotto
#' @export
GiottoClass::makePseudoVisium
#' @export
GiottoClass::objHistory
#' @export
GiottoClass::orthoGrid
#' @export
GiottoClass::overlapImagesToMatrix
#' @export
GiottoClass::overlapToMatrixMultiPoly
#' @export
GiottoClass::pDataDT
#' @export
GiottoClass::plotGiottoImage
#' @export
GiottoClass::polyStamp
#' @export
GiottoClass::readCellMetadata
#' @export
GiottoClass::readDimReducData
#' @export
GiottoClass::readExprData
#' @export
GiottoClass::readExprMatrix
#' @export
GiottoClass::readFeatData
#' @export
GiottoClass::readFeatMetadata
#' @export
GiottoClass::readGiottoInstructions
#' @export
GiottoClass::readNearestNetData
#' @export
GiottoClass::readPolygonData
#' @export
GiottoClass::readSpatEnrichData
#' @export
GiottoClass::readSpatLocsData
#' @export
GiottoClass::readSpatNetData
#' @export
GiottoClass::reconnectGiottoImage
#' @export
GiottoClass::rectVertices
#' @export
GiottoClass::removeCellAnnotation
#' @export
GiottoClass::removeFeatAnnotation
#' @export
GiottoClass::removeGiottoEnvironment
#' @export
GiottoClass::replaceGiottoInstructions
#' @export
GiottoClass::rescalePolygons
#' @export
GiottoClass::saveGiotto
#' @export
GiottoClass::setCellMetadata
#' @export
GiottoClass::setDimReduction
#' @export
GiottoClass::setExpression
#' @export
GiottoClass::setFeatureInfo
#' @export
GiottoClass::setFeatureMetadata
#' @export
GiottoClass::setGiottoImage
#' @export
GiottoClass::setMultiomics
#' @export
GiottoClass::setNearestNetwork
#' @export
GiottoClass::setPolygonInfo
#' @export
GiottoClass::setSpatialEnrichment
#' @export
GiottoClass::setSpatialGrid
#' @export
GiottoClass::setSpatialLocations
#' @export
GiottoClass::setSpatialNetwork
#' @export
GiottoClass::seuratToGiotto
#' @export
GiottoClass::seuratToGiottoV4
#' @export
GiottoClass::seuratToGiottoV5
#' @export
GiottoClass::showGiottoCellMetadata
#' @export
GiottoClass::showGiottoDimRed
#' @export
GiottoClass::showGiottoExpression
#' @export
GiottoClass::showGiottoFeatInfo
#' @export
GiottoClass::showGiottoFeatMetadata
#' @export
GiottoClass::showGiottoImageNames
#' @export
GiottoClass::showGiottoInstructions
#' @export
GiottoClass::showGiottoNearestNetworks
#' @export
GiottoClass::showGiottoSpatEnrichments
#' @export
GiottoClass::showGiottoSpatGrids
#' @export
GiottoClass::showGiottoSpatLocs
#' @export
GiottoClass::showGiottoSpatNetworks
#' @export
GiottoClass::showGiottoSpatialInfo
#' @export
GiottoClass::showProcessingSteps
#' @export
GiottoClass::smoothGiottoPolygons
#' @export
GiottoClass::spatQueryGiottoPolygons
#' @export
GiottoClass::spatialExperimentToGiotto
#' @export
GiottoClass::stitchFieldCoordinates
#' @export
GiottoClass::stitchGiottoLargeImage
#' @export
GiottoClass::stitchFieldCoordinates
#' @export
GiottoClass::subsetGiotto
#' @export
GiottoClass::subsetGiottoLocs
#' @export
GiottoClass::subsetGiottoLocsMulti
#' @export
GiottoClass::subsetGiottoLocsSubcellular
#' @export
GiottoClass::tessellate
#' @export
GiottoClass::triGrid
#' @export
GiottoClass::updateGiottoImage
#' @export
GiottoClass::updateGiottoImageMG
#' @export
GiottoClass::updateGiottoLargeImage
#' @export
GiottoClass::updateGiottoObject
#' @export
GiottoClass::updateGiottoPointsObject
#' @export
GiottoClass::updateGiottoPolygonObject
#' @export
GiottoClass::writeGiottoLargeImage








# GiottoVisuals ####

#' @export
GiottoVisuals::addGiottoImageToSpatPlot
#' @export
GiottoVisuals::dimCellPlot
#' @export
GiottoVisuals::dimCellPlot2D
#' @export
GiottoVisuals::dimFeatPlot2D
#' @export
GiottoVisuals::dimGenePlot3D # TODO remove in next version
#' @export
GiottoVisuals::dimFeatPlot3D
#' @export
GiottoVisuals::dimPlot
#' @export
GiottoVisuals::dimPlot2D
#' @export
GiottoVisuals::dimPlot3D
#' @export
GiottoVisuals::getColors
#' @export
GiottoVisuals::giottoSankeyPlan
#' @export
GiottoVisuals::plotHeatmap
#' @export
GiottoVisuals::plotMetaDataCellsHeatmap
#' @export
GiottoVisuals::plotMetaDataHeatmap
#' @export
GiottoVisuals::plotPCA
#' @export
GiottoVisuals::plotPCA_2D
#' @export
GiottoVisuals::plotPCA_3D
#' @export
GiottoVisuals::plotStatDelaunayNetwork
#' @export
GiottoVisuals::plotTSNE
#' @export
GiottoVisuals::plotTSNE_2D
#' @export
GiottoVisuals::plotTSNE_3D
#' @export
GiottoVisuals::plotUMAP
#' @export
GiottoVisuals::plotUMAP_2D
#' @export
GiottoVisuals::plotUMAP_3D
#' @export
GiottoVisuals::sankeyLabel
#' @export
GiottoVisuals::`sankeyLabel<-`
#' @export
GiottoVisuals::sankeyPlot
#' @export
GiottoVisuals::sankeyRelate
#' @export
GiottoVisuals::`sankeyRelate<-`
#' @export
GiottoVisuals::sankeySet
#' @export
GiottoVisuals::sankeySetAddresses
#' @export
GiottoVisuals::showColorInstructions
#' @export
GiottoVisuals::showClusterDendrogram
#' @export
GiottoVisuals::showClusterHeatmap
#' @export
GiottoVisuals::showSaveParameters
#' @export
GiottoVisuals::spatCellPlot
#' @export
GiottoVisuals::spatCellPlot2D
#' @export
GiottoVisuals::spatDeconvPlot
#' @export
GiottoVisuals::spatDimCellPlot
#' @export
GiottoVisuals::spatDimCellPlot2D
#' @export
GiottoVisuals::spatDimFeatPlot2D
#' @export
GiottoVisuals::spatDimGenePlot3D # TODO remove in next version
#' @export
GiottoVisuals::spatDimFeatPlot3D
#' @export
GiottoVisuals::spatDimPlot
#' @export
GiottoVisuals::spatDimPlot2D
#' @export
GiottoVisuals::spatDimPlot3D
#' @export
GiottoVisuals::spatFeatPlot2D
#' @export
GiottoVisuals::spatFeatPlot2D_single
#' @export
GiottoVisuals::spatGenePlot3D # TODO remove in next version
#' @export
GiottoVisuals::spatFeatPlot3D
#' @export
GiottoVisuals::spatPlot
#' @export
GiottoVisuals::spatPlot2D
#' @export
GiottoVisuals::spatPlot3D
#' @export
GiottoVisuals::spatInSituPlotDensity
#' @export
GiottoVisuals::spatInSituPlotHex
#' @export
GiottoVisuals::spatInSituPlotPoints
#' @export
GiottoVisuals::spatNetwDistributions
#' @export
GiottoVisuals::spatNetwDistributionsDistance
#' @export
GiottoVisuals::spatNetwDistributionsKneighbors
#' @export
GiottoVisuals::subsetSankeySet
#' @export
GiottoVisuals::violinPlot

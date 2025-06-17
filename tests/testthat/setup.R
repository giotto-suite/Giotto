# Global options for testing
options("lifecycle_verbosity" = "quiet")

# Test Data Setup
# ---------------
setup_test_data <- function() {
    # Load test objects
    vis <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    
    # load example data
    ex <- GiottoData::loadSubObjectMini("exprObj")
    sl <- GiottoData::loadSubObjectMini("spatLocsObj")
    cm <- GiottoData::loadSubObjectMini("cellMetaObj")
    fm <- GiottoData::loadSubObjectMini("featMetaObj")
    sn <- GiottoData::loadSubObjectMini("spatialNetworkObj")
    enr <- GiottoData::loadSubObjectMini("spatEnrObj")
    nn <- GiottoData::loadSubObjectMini("nnNetObj")
    dr <- GiottoData::loadSubObjectMini("dimObj")
    gpoly <- GiottoData::loadSubObjectMini("giottoPolygon")
    gpoints <- GiottoData::loadSubObjectMini("giottoPoints")
    
    list(
        ex = ex, sl = sl, cm = cm, fm = fm, sn = sn,
        enr = enr, nn = nn, dr = dr, gpoly = gpoly,
        gpoints = gpoints, vis = vis
    )
}

test_data <- setup_test_data()
set.seed(1234) # reproducibility seed

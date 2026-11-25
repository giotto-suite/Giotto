test_that("Giotto's own methods declare the paper they implement", {
    expect_match(format(methodCitation(hmrfParam(k = 3)), style = "text"),
        "10.1038/nbt.4260", fixed = TRUE)
    expect_match(format(methodCitation(neighborhoodPCAParam()), style = "text"),
        "10.1038/s41588-024-01664-3", fixed = TRUE)
})

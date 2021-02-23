
## following the walkthrough tutorial, walkthrough.md

library(conos)
library(dplyr)
library(pagoda2)

con <- Conos$new(small_panel.preprocessed, n.cores=1)

test_that("check Conos object, samples", {
	expect_equal(length(con$samples), 2)
})

## warnings due to small panel
suppressWarnings(con$buildGraph(ncomps=25))
suppressWarnings(con$findCommunities(method = igraph::walktrap.community, steps=7))
suppressWarnings(con$embedGraph(alpha=0.001, sgd_batched=1e8))


test_that("check Conos object, number of clusters", {
	expect_equal(length(con$clusters$walktrap$groups), 59)
	expect_equal(length(con$clusters$walktrap$result), 2)
})

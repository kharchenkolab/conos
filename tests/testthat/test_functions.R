
## following the walkthrough tutorial
## https://github.com/kharchenkolab/conos/blob/master/vignettes/walkthrough.md

library(conos)
library(dplyr)
library(pagoda2)

## warnings due to small panel

options(warn=-1)

con <- Conos$new(small_panel.preprocessed, n.cores=1)

test_that("check Conos object, samples", {
	expect_equal(length(con$samples), 2)
})

con$buildGraph(ncomps=25)
con$findCommunities(method = igraph::walktrap.community, steps=7)
con$embedGraph(alpha=0.001, sgd_batched=1e8)  


test_that("check Conos object, number of clusters", {
	expect_equal(length(con$clusters$walktrap$groups), 59)
	expect_equal(length(con$clusters$walktrap$result), 2)
})


## following the walkthrough tutorial
## https://github.com/kharchenkolab/conos/blob/master/vignettes/walkthrough.md

library(conos)
library(dplyr)
library(pagoda2)

## load the data
panel <- conosPanel::panel

## pre-processing
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=1, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

con <- Conos$new(panel.preprocessed, n.cores=1)

test_that("check Conos object, samples", {
	expect_equal(length(con$samples), 4)
})

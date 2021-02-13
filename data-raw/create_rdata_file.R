
## create Rdata file used in tests and examples

library(conos)
library(dplyr)
library(pagoda2)

small_panel = readRDS("panel.rds")  ## located in older versions of Conos

## simply subsetted panel.rds

small_panel[[1]] = small_panel[[1]][c(1:1000), c(1:100)]
small_panel[[2]] = small_panel[[2]][c(1:1000), c(1:100)]

## dim(small_panel[[1]]), 1000, 100
## dim(small_panel[[2]]), 1000, 100
## length(small_panel), 2

small_panel.preprocessed = lapply(small_panel, basicP2proc, n.cores=1, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
save(small_panel.preprocessed, file="small_panel.preprocessed.rda", compress="xz", compression_level=9) 
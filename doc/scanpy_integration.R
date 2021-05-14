## ---- message=FALSE, warning=FALSE--------------------------------------------
install.packages('conosPanel', repos='https://kharchenkolab.github.io/drat/', type='source')

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(conos)
panel <- conosPanel::panel

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=1, min.cells.per.gene=0, n.odgenes=2e3, 
                             get.largevis=FALSE, make.geneknn=FALSE)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con <- Conos$new(panel.preprocessed, n.cores=1)
con$buildGraph(k=15, k.self=5, space='PCA', ncomps=30)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$findCommunities()
con$embedGraph(method="UMAP")

## ---- message=FALSE, warning=FALSE--------------------------------------------
metadata <- data.frame(Cluster=con$clusters$leiden$groups)

## ---- message=FALSE, warning=FALSE--------------------------------------------
## use current directory
exchange_dir <- "."
hdf5file = "example.h5"
saveConosForScanPy(con, output.path=exchange_dir, hdf5_filename=hdf5file, verbose=TRUE)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(rhdf5)
metadata = h5read(paste0(exchange_dir, "/example.h5"), 'metadata/metadata.df')
head(metadata, 4)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(rhdf5)
library(Matrix)
rawcountMat = h5read(paste0(exchange_dir, "/example.h5"), 'raw_count_matrix')
raw_count_matrix = sparseMatrix(x = as.numeric(rawcountMat$data),  
    dims = as.numeric(c(rawcountMat$shape[1], rawcountMat$shape[2])), 
    p = as.numeric(rawcountMat$indptr), 
    i = rawcountMat$indices, index1=FALSE)


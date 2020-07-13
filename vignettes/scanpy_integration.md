Convert Conos object to ScanPy
==============================

```r
library(conos)
```

Load data:


```r
panel <- readRDS(file.path(find.package('conos'), 'extdata', 'panel.rds'))
```

Run Pagoda 2:


```r
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, 
                             get.largevis=FALSE, make.geneknn=FALSE)
```

```
## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
## calculating variance fit ... using gam 172 overdispersed genes ... 172persisting ... done.
## running PCA using 2000 OD genes .... done
## running tSNE using 4 cores:
## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
## calculating variance fit ... using gam 159 overdispersed genes ... 159persisting ... done.
## running PCA using 2000 OD genes .... done
## running tSNE using 4 cores:
## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
## calculating variance fit ... using gam 251 overdispersed genes ... 251persisting ... done.
## running PCA using 2000 OD genes .... done
## running tSNE using 4 cores:
## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
## calculating variance fit ... using gam 168 overdispersed genes ... 168persisting ... done.
## running PCA using 2000 OD genes .... done
## running tSNE using 4 cores:
```

Align datasets:


```r
con <- Conos$new(panel.preprocessed, n.cores=4)
con$buildGraph(k=15, k.self=5, space='PCA', ncomps=30)
```

```
## found 0 out of 6 cached PCA  space pairs ... running 6 additional PCA  space pairs  done
## inter-sample links using  mNN   done
## local pairs local pairs  done
## building graph ..done
```

Find clusters and create embedding:


```r
con$findCommunities()
con$embedGraph(method="UMAP")
```

```
## Convert graph to adjacency list...
## Done
## Estimate nearest neighbors and commute times...
## Estimating hitting distances: 01:43:26.
## Done.
## Estimating commute distances: 01:43:40.
## Hashing adjacency list: 01:43:40.
## Done.
## Estimating distances: 01:43:40.
## Done
## Done.
## All done!: 01:43:44.
## Done
## Estimate UMAP embedding...
```

```
## 01:43:44 UMAP embedding parameters a = 0.0267 b = 0.7906
```

```
## 01:43:44 Read 12000 rows and found 1 numeric columns
```

```
## 01:43:44 Commencing smooth kNN distance calibration using 4 threads
```

```
## 01:43:45 Initializing from normalized Laplacian + noise
```

```
## 01:43:45 Commencing optimization for 1000 epochs, with 332028 positive edges using 4 threads
```

```
## 01:43:54 Optimization finished
```

```
## Done
```

Prepare metadata (can be any type of clustering of all the cells):


```r
metadata <- data.frame(Cluster=con$clusters$leiden$groups)
```

Save data (set `exchange_dir` to your path):


```r
exchange_dir <- "~/scanpy_demo"
hdf5file = "example.h5"
dir.create(exchange_dir)
```



```r
saveConosForScanPy(con, output.path=exchange_dir, hdf5_filename=hdf5file, verbose=TRUE)
```

```
## Merge raw count matrices...	Done.
## Save the embedding...		Done.
## Save graph matrices...		Done.
## Write data to disk...		
```

```
## You created a large dataset with compression and chunking.
## The chunk size is equal to the dataset dimensions.
## If you want to read subsets of the dataset, you should testsmaller chunk sizes to improve read times.
## You created a large dataset with compression and chunking.
## The chunk size is equal to the dataset dimensions.
## If you want to read subsets of the dataset, you should testsmaller chunk sizes to improve read times.
```

```
## All Done!
```

Users can then access the data saved to the HDF5 file, e.g. to access metadata, run:



```r
library(rhdf5)
metadata = h5read(paste0(exchange_dir, "/example.h5"), 'metadata/metadata.df')
head(metadata, 4)
```

```
##                                 CellId           Dataset
## 1 MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1 MantonBM1_HiSeq_1
## 2 MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1 MantonBM1_HiSeq_1
## 3 MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1 MantonBM1_HiSeq_1
## 4 MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1 MantonBM1_HiSeq_1
```

All possible fields included in the output HDF5 file are:

* `raw_count_matrix`: the sparse `dgCMatrix` of raw counts 
	* `raw_count_matrix/data`: the matrix entries
	* `raw_count_matrix/shape`: the matrix dimensions
	* `raw_count_matrix/indices`: 0-based vector of non-zero matrix entries
	* `raw_count_matrix/indptr`: vector of pointers, one for each column (or row), to the initial (zero-based) index of elements in the column (or row)
* `metadata/metadata.df`: the `data.frame` of metadata values
* `genes/genes.df`: the `data.frame` of genes
* `count_matrix`: the sparse `dgCMatrix` of normalized counts, if `cm.norm=TRUE`
	* `count_matrix/data`: the matrix entries
	* `count_matrix/shape`: the matrix dimensions
	* `count_matrix/indices`: 0-based vector of non-zero matrix entries
	* `count_matrix/indptr`: vector of pointers, one for each column (or row), to the initial (zero-based) index of elements in the column (or row)
* `embedding/embedding.df`:  the `data.frame` of the conos embedding, if `embedding=TRUE`
* `pseudopca/pseudopca.df`:  the `data.frame` of an emulated PCA by embedding the graph to a space with `n.dims` dimensions and save it as a pseudoPCA, if `pseudopca=TRUE`
* `pca/pca.df`: the `data.frame` of the PCA of all the samples (not batch corrected), if `pca=TRUE`
* `graph_connectivities`: the sparse `dgCMatrix` of graph connectivites, if `alignment.graph=TRUE`
	* `graph_connectivities/data`: the matrix entries
	* `graph_connectivities/shape`: the matrix dimensions
	* `graph_connectivities/indices`: 0-based vector of non-zero matrix entries
	* `graph_connectivities/indptr`: vector of pointers, one for each column (or row), to the initial (zero-based) index of elements in the column (or row)
* `graph_distances`: the sparse `dgCMatrix` of graph distances, if `alignment.graph=TRUE`
	* `graph_distances/data`: the matrix entries
	* `graph_distances/shape`: the matrix dimensions
	* `graph_distances/indices`: 0-based vector of non-zero matrix entries
	* `graph_distances/indptr`: vector of pointers, one for each column (or row), to the initial (zero-based) index of elements in the column (or row)



In order to read in the `dcGMatrix` again, simply use the `Matrix` package as follows:



```r
library(rhdf5)
library(Matrix)
rawcountMat = h5read(paste0(exchange_dir, "/example.h5"), 'raw_count_matrix')
raw_count_matrix = sparseMatrix(x = as.numeric(rawcountMat$data),  
    dims = as.numeric(c(rawcountMat$shape[1], rawcountMat$shape[2])), 
    p = as.numeric(rawcountMat$indptr), 
    i = rawcountMat$indices, index1=FALSE)
```

**Note:** Please set `index1=FALSE` as the index vectors are 0-based. For more details, see the documentation for [Matrix::sparseMatrix()](https://www.rdocumentation.org/packages/Matrix/versions/1.2-18/topics/sparseMatrix)


Convert Conos Object to ScanPy
==============================

First load the Conos library and the corresponding data:


```r
library(conos)
```

```
## Loading required package: Matrix
```

```
## Loading required package: igraph
```

```
## 
## Attaching package: 'igraph'
```

```
## The following objects are masked from 'package:stats':
## 
##     decompose, spectrum
```

```
## The following object is masked from 'package:base':
## 
##     union
```

```r
panel <- readRDS(file.path(find.package('conos'), 'extdata', 'panel.rds'))
```

Run Pagoda2:


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

Align the datasets:


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

Save the data (set `exchange_dir` to your path):


```r
exchange_dir <- "~/scanpy_demo"
dir.create(exchange_dir)
```

```
## Warning in dir.create(exchange_dir): '/Users/evanbiederstedt/scanpy_demo'
## already exists
```

```r
saveConosForScanPy(con, output.path=exchange_dir, verbose=T)
```

```
## Merge raw count matrices...	Done.
## Save the embedding...		
```

```
## Warning in saveConosForScanPy(con, output.path = exchange_dir, verbose = T): 
##  No embedding found in the conos object. Skipping...
```

```
## Save graph matrices...		Done.
## Write data to disk...		Done.
## All Done!
```


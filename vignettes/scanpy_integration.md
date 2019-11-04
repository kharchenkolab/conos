Convert Conos object to ScanPy
================

``` r
library(conos)
```

Load data:

``` r
panel <- readRDS(file.path(find.package('conos'), 'extdata', 'panel.rds'))
```

Run Pagoda 2:

``` r
library(pagoda2)
```

    ## 

``` r
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, 
                             get.largevis=FALSE, make.geneknn=FALSE)
```

    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 171 overdispersed genes ... 171 persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:
    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 159 overdispersed genes ... 159 persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:
    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 248 overdispersed genes ... 248 persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:
    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 166 overdispersed genes ... 166 persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:

Align datasets:

``` r
con <- Conos$new(panel.preprocessed, n.cores=4)
con$buildGraph(k=15, k.self=5, space='PCA', ncomps=30)
```

    ## found 0 out of 6 cached PCA  space pairs ... running 6 additional PCA  space pairs  done
    ## inter-sample links using  mNN   done
    ## local pairs local pairs  done
    ## building graph ..done

Find clusters and create embedding:

``` r
con$findCommunities()
con$embedGraph(method="UMAP")
```

    ## Convert graph to adjacency list...
    ## Done
    ## Estimate nearest neighbors and commute times...
    ## Estimating hitting distances: 11:23:45.
    ## Done.
    ## Estimating commute distances: 11:24:01.
    ## Hashing adjacency list: 11:24:01.
    ## Done.
    ## Estimating distances: 11:24:02.
    ## Done
    ## Done.
    ## All done!: 11:24:05.
    ## Done
    ## Estimate UMAP embedding...

    ## 11:24:05 UMAP embedding parameters a = 0.0267 b = 0.7906

    ## 11:24:05 Read 12000 rows and found 1 numeric columns

    ## 11:24:05 Commencing smooth kNN distance calibration using 4 threads

    ## 11:24:06 Initializing from normalized Laplacian + noise

    ## 11:24:07 Commencing optimization for 1000 epochs, with 330992 positive edges using 4 threads

    ## 11:24:23 Optimization finished

    ## Done

Prepare metadata (can be any type of clustering of all the cells):

``` r
metadata <- data.frame(Cluster=con$clusters$leiden$groups)
```

Save data (set `exchange_dir` to your path):

``` r
exchange_dir <- "."
dir.create(exchange_dir)
saveConosForScanPy(con, output.path=exchange_dir, metadata.df=metadata, 
                   norm=TRUE, embed=TRUE, pseudo.pca=TRUE, pca=TRUE, 
                   connect=TRUE, verbose=TRUE)
```

    ## Merge raw count matrices...  Done.
    ## Merge count matrices...      Done.
    ## Save the embedding...        Done.
    ## Create psudo-PCA space...    Done.
    ## Save PCA space...            Done.
    ## Save graph matrices...       Done.
    ## Write data to disk...        Done.
    ## All Done!

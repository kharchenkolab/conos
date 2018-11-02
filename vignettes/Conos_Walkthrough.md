Conos Walkthrough
================

In this tutorial we will go over the analysis of a panel of samples using conos. Conos objects can be used to identify clusters of corresponding cells across panels of samples from similar or dissimilar sources, with different degrees of cell type overlap. Here we will identify corresponding clusters accorss a panel of bone marrow (BM) and cord blood (CB) by generating a joint graph with the cells from all the samples. We will use the graph to propagate labels from a single labelled sample to other samples and finally perform differential expression between BM and CM samples.

Preliminary
===========

First of all we need to load the libraries we will use, we will need pagoda2, conos and some helper libraries.

``` r
library(pagoda2)
```

    ## 

    ## Warning: replacing previous import 'igraph::%>%' by 'magrittr::%>%' when
    ## loading 'pagoda2'

``` r
library(conos)
```

    ## Loading required package: Matrix

    ## Warning: replacing previous import 'igraph::%>%' by 'dplyr::%>%' when
    ## loading 'conos'

    ## 
    ## Attaching package: 'conos'

    ## The following objects are masked from 'package:pagoda2':
    ## 
    ##     buildWijMatrix, projectKNNs

``` r
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(nbHelpers)
```

    ## 
    ## Attaching package: 'nbHelpers'

    ## The following object is masked from 'package:pagoda2':
    ## 
    ##     namedNames

Loading the data
================

Next we will load a previously prepared panel of samples. This panel is made up of 16 cord blood and bone marrow samples. It has already been quality controlled so we don't need to worry about that. All samples have been subset to exactly 3000 cells.

``` r
panel <- readRDS(file.path(find.package('conos'),'extdata','panel.rds'))
```

Let's take a look at the panel. The panel is a named list of sparse matrices (type dgCMatrix).

``` r
str(panel,1)
```

    ## List of 8
    ##  $ MantonBM1_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonBM2_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonBM3_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonBM4_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonCB1_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonCB2_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonCB3_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonCB4_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots

Before we continue it is very important to make sure that cells in our panel are uniquely named. No two cells (even in different samples) should be named identically. In this case the cells have been prefixed by sample id, so there will not be any collisions. However in most cases you will have to prefix the cells before continuing.

``` r
head(colnames(panel[[1]]))
```

    ## [1] "MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1"
    ## [2] "MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1"
    ## [3] "MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1"
    ## [4] "MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1"
    ## [5] "MantonBM1_HiSeq_1-TATTACCCAAAGGAAG-1"
    ## [6] "MantonBM1_HiSeq_1-CGCCAAGCATCTGGTA-1"

Next we will subset the panel to just 4 samples. Conos has been tested with 10s of samples, however larger panel sizes increase running time considerably. For brevity we will only use 2 CB and 2 BM samples. After subsetting we will generate pagoda2 apps for these four samples using the basicP2proc function.

``` r
panel <- panel[c('MantonBM1_HiSeq_1','MantonBM2_HiSeq_1','MantonCB1_HiSeq_1','MantonCB2_HiSeq_1')]
panel.p2 <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0)
```

    ## 3000 cells, 18535 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 171 overdispersed genes ... 171 persisting ... done.
    ## running PCA using 3000 OD genes .... done
    ## Estimating embeddings.
    ## running tSNE using 4 cores:
    ## 3000 cells, 17720 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 159 overdispersed genes ... 159 persisting ... done.
    ## running PCA using 3000 OD genes .... done
    ## Estimating embeddings.
    ## running tSNE using 4 cores:
    ## 3000 cells, 18672 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 248 overdispersed genes ... 248 persisting ... done.
    ## running PCA using 3000 OD genes .... done
    ## Estimating embeddings.
    ## running tSNE using 4 cores:
    ## 3000 cells, 17530 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 166 overdispersed genes ... 166 persisting ... done.
    ## running PCA using 3000 OD genes .... done
    ## Estimating embeddings.
    ## running tSNE using 4 cores:

Let's look at the output of our processing. We now have a named list of pagoda2 objects, which is the starting point for the analysis with Conos.

``` r
str(panel.p2,1)
```

    ## List of 4
    ##  $ MantonBM1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant
    ##  $ MantonBM2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant
    ##  $ MantonCB1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant
    ##  $ MantonCB2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant

Processing with Conos
=====================

We will now construct a Conos object for this panel of samples. At this point we haven't calculated anything. We have just generated an object that contains the samples. Note that at this step we also set the n.cores parameter. The graph generation with Conos can take advantage of parallel processing, so use as many physical cores as you have available here.

``` r
con <- Conos$new(panel.p2,n.cores=4)
```

Our original pagoda2 apps are now saved in the conos object (if you are short of memory you can go ahead and delete the originals).

``` r
str(con$samples,1)
```

    ## List of 4
    ##  $ MantonBM1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant
    ##  $ MantonBM2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant
    ##  $ MantonCB1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant
    ##  $ MantonCB2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
    ##   ..and 34 methods, of which 20 are  possibly relevant

We can now plot a panel of these samples using the clusters we identified by examining each sample on its own. We note that each sample has an independent set of clusters that bears no relationship to clusters in other sample (for example note cluster 9).

``` r
con$plotPanel(clustering="multilevel", use.local.clusters=T, title.size=6)
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-9-1.png)

Next we will build the graph emcompasses all the samples. We do that by pairwise projecting samples onto a common space and establishing KNN of mNN pairs between the samples. We then append iwthin sample KNN neighbours to the graph to ensure that all the cell are included in the graph.

``` r
con$buildGraph(k=30, k.self=10, k.self.weight=0.1, space='PCA', matching.method='mNN', metric='angular', data.type='counts', l2.sigma=1e5, var.scale =TRUE, ncomps=50, n.odgenes=1000, return.details=T,neighborhood.average=FALSE,neighborhood.average.k=10, exclude.pairs=NULL, exclude.samples=NULL, common.centering=TRUE , verbose=TRUE)
```

    ## found 0 out of 6 cached PCA  space pairs ... running 6 additional PCA  space pairs  done
    ## inter-sample links using  mNN   done
    ## local pairs  done

We next use the graph we identified to get global clusters. Here we use mutlievel to obtain clusters.

``` r
con$findCommunities(method=multilevel.community, min.group.size=0)
```

We can now plot the clusters we obtained. Note that the cluster numbers between different samples now correspond to the same cell type. Also not the presence of cluster 5 in BM samples only, but not in CB.

``` r
con$plotPanel(font.size=4)
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-12-1.png)

Check an expression pattern of a specific gene across all the individual embeddings.

``` r
con$plotPanel(gene = 'GZMK')
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-13-1.png)

Next we embed and visualize the complete joint graph:

``` r
con$plotGraph()
```

    ## Estimating embeddings.

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-14-1.png)

We note that the graph captures the population structure irrespecively of the sample of origin of each cell.

``` r
con$plotGraph(color.by='sample',mark.groups=F,alpha=0.1,show.legend=T)
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-15-1.png)

Other community detection methods can provide a more sensitive and hierarchical view of the subpopulation structure. Here we run walktrap community detection method on the same joint graph:

``` r
con$findCommunities(method = igraph::walktrap.community, steps=4)
```

Visualize new clusters:

``` r
con$plotPanel(clustering='walktrap',font.size=4)
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-17-1.png)

New clustering, as viewed on a joint graph:

``` r
con$plotGraph(clustering='walktrap')
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-18-1.png)

Label propagation
=================

One of the uses of this graph is to propagate labels. For example in some cases we will only have information about the cell types in one of the samples and we want to automatically label the other samples.

Load cell annotation for one of the samples. Here we are loading annotations we generated with the pagoda2 web tool. We first load the annotations, but because we made the with the lasso tool there is no guarantee that some cells are not included in multiple selctions. For this reason we use removeSelectionOverlaps before conversing the selections to a factor that we can use downstream.

``` r
cellannot <- readPagoda2SelectionFile(file.path(find.package('conos'),'extdata','selections.txt'))
cellannot <- removeSelectionOverlaps(cellannot)
cellannot <- factorFromP2Selection(cellannot)
```

Next we plot our panel with the annotations we made. This is to verify that the annotated cells are indeed in only one sample and that the other samples are unlabelled.

``` r
con$plotPanel(groups = cellannot)
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-20-1.png)

Next let's propagaes the labels from the one annotated sample to the other samples.

``` r
tmp.names <- names(cellannot)
cellannot.c <- as.character(cellannot)
names(cellannot.c) <- tmp.names
new.label.probabilities <- con$propagateLabels(labels = cellannot,verbose = T)
```

This function returns probabilites of the cell belonging to each group, we can assign each cell to the the the cell type with the highest probability.

``` r
new.annot <- colnames(new.label.probabilities)[apply(new.label.probabilities,1,which.max)]
names(new.annot) <- rownames(new.label.probabilities)
head(new.annot)
```

    ## MantonBM1_HiSeq_1-AGGGATGCAGGTGCCT-1 MantonBM2_HiSeq_1-CTGATAGAGCGTTCCG-1 
    ##                              "Tcyto"                              "Tcyto" 
    ## MantonBM1_HiSeq_1-GTAACTGCAGATCCAT-1 MantonBM1_HiSeq_1-CAAGATCTCTGAGTGT-1 
    ##                              "Tcyto"                              "Tcyto" 
    ## MantonBM2_HiSeq_1-GGAAAGCCAGACGCCT-1 MantonBM1_HiSeq_1-ACACCGGAGTAACCCT-1 
    ##                              "Tcyto"                              "Tcyto"

We not see that all our samples have been labelled automagically!

``` r
con$plotPanel(groups = new.annot)
```

![](Conos_Walkthrough_files/figure-markdown_github/unnamed-chunk-23-1.png)

Differential expression
=======================

Once we have identified a joint clustering of cells that captures cell relationships between samples, we want to ask what is different between the cells of these populations between specific samples types, in this case CB and BM samples. Conos provides routines to be able to do that.

First we need to define our sample groups

``` r
samplegroups <- list(
  bm = c("MantonBM1_HiSeq_1","MantonBM2_HiSeq_1"),
  cb = c("MantonCB1_HiSeq_1","MantonCB2_HiSeq_1")
)
```

We can then run differential expression between cells in these groups

``` r
de.multilevel <- getPerCellTypeDE(con, groups=as.factor(new.annot),sample.groups = samplegroups, ref.level='bm', n.cores=4)
```

...and examine the output

``` r
str(de.multilevel[1:3], 2)
```

    ## List of 3
    ##  $ Bcells:List of 3
    ##   ..$ res         :'data.frame': 15032 obs. of  6 variables:
    ##   ..$ cm          :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   ..$ sample.groups:List of 2
    ##  $ Mono  :List of 3
    ##   ..$ res         :'data.frame': 15032 obs. of  6 variables:
    ##   ..$ cm          :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   ..$ sample.groups:List of 2
    ##  $ Tcyto :List of 3
    ##   ..$ res         :'data.frame': 15032 obs. of  6 variables:
    ##   ..$ cm          :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   ..$ sample.groups:List of 2

Let's look at the results for the Bcells

``` r
res <- de.multilevel[['Bcells']]$res
head(res[order(res$padj,decreasing = FALSE),])
```

    ##                baseMean log2FoldChange     lfcSE       stat       pvalue
    ## JCHAIN         767.7050      -4.663820 0.4321216 -10.792842 3.721025e-27
    ## IGKC          6221.7068      -3.565991 0.4312053  -8.269820 1.341612e-16
    ## IGHA1         1797.2574     -13.915045 1.9554493  -7.116034 1.110765e-12
    ## HBG2          4343.2021      11.461933 1.6753089   6.841683 7.826815e-12
    ## IGHG1          328.2118     -10.051699 1.4653478  -6.859600 6.905386e-12
    ## RP11-386I14.4  464.4805       2.670045 0.4179233   6.388839 1.671497e-10
    ##                       padj
    ## JCHAIN        5.535025e-23
    ## IGKC          9.978242e-13
    ## IGHA1         5.507545e-09
    ## HBG2          2.328477e-08
    ## IGHG1         2.328477e-08
    ## RP11-386I14.4 4.143920e-07

With correction
---------------

In certain cases we observe that differential expression will result in the similar genes between multiple cell types. This may be due to genuine biological reasons (similar response), due to background, or due to other effects. Conos can calculate a mean expression vector between the two conditions and subtract this from all the comparisons, so observer the cell-type specific effect.

``` r
## Calculate correction
fc.correction <- getCorrectionVector(con, groups=as.factor(new.annot),sample.groups = samplegroups, ref.level='bm', n.cores=4,correction.method='varianceweighted')
fc.correction[is.na(fc.correction)] <- 0

## Use corrected version
de.multilevel.corrected <- getPerCellTypeDECorrected(con, groups=as.factor(new.annot),sample.groups = samplegroups, ref.level='bm', n.cores=4, correction = fc.correction)
```

``` r
res <- as.data.frame(de.multilevel.corrected[['Mono']]$res)
head(res[order(res$padj,decreasing = FALSE),])
```

    ##        baseMean log2FoldChange     lfcSE      stat       pvalue
    ## APOC1  47.74995      -8.773700 1.5279227 -5.742240 9.343201e-09
    ## MPO    25.56018      -4.927059 1.0865772 -4.534477 5.774648e-06
    ## CDC20  17.45361      -5.509215 1.2398522 -4.443445 8.852984e-06
    ## CCNA2  16.47920      -4.542436 1.1499688 -3.950051 7.813448e-05
    ## CDK1   15.39750      -4.370825 1.1978772 -3.648809 2.634588e-04
    ## AKR1C3 24.39723      -3.503196 0.9838868 -3.560568 3.700540e-04
    ##               padj
    ## APOC1  0.000138214
    ## MPO    0.042712181
    ## CDC20  0.043654066
    ## CCNA2  0.288960828
    ## CDK1   0.779469268
    ## AKR1C3 0.782029897

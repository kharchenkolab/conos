---
title: "Conos Walkthrough"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{"Raster geoms"}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


- [Loading the Data](#loading-the-data)
  * [Pre-processing with Pagoda2](#pre-processing-with-pagoda2)
  * [Pre-processing with Seurat](#pre-processing-with-seurat)
- [Integrating Datasets with Conos](#integrating-datasets-with-conos)
  * [Visualization](#visualization)
  * [Changing embedding parameters](#changing-embedding-parameters)
    + [largeVis](#largevis)
    + [UMAP](#umap)
- [Exploring Hierarchical Community Structure](#exploring-hierarchical-community-structure)
  * [Using code](#using-code)
  * [Using Shiny Application](#using-shiny-application)
- [Label Propagation](#label-propagation)
  * [General workflow](#general-workflow)
- [Differential Expression](#differential-expression)
  * [Cluster markers](#cluster-markers)
  * [Differential expression between sample groups](#de-between-sample-groups)
    + [Simple run](#simple-run)
- [Forcing Better Alignment](#forcing-better-alignment)

In this tutorial, we will go over the analysis of a panel of samples
using Conos. Conos objects can be used to identify clusters of
corresponding cells across panels of samples from similar or dissimilar
sources, with different degrees of cell type overlap. Here we will
identify the clusters of corresponding cells across a panel of bone marrow (BM) and
cord blood (CB) by generating a joint graph with the cells from all the
samples. We will then use this graph to propagate labels from a single
labelled sample to other samples, and finally perform differential
expression between the BM and CB samples.

First, let’s load Conos library:


```r
library(conos)
library(dplyr)
```

# Loading the Data

Next we will load a previously prepared panel of samples. This panel was
made up of 16 cord blood and bone marrow samples, but for convenience, we will here focus on a smaller subset of just 4 samples. All samples have been subset to a size of exactly 3000 cells. 

**Note:** When starting with your own panel, we
recommend filtering out low-count/poor-quality/dying cells, as is standard for quality control.


```r
panel <- readRDS(file.path(find.package('conos'),'extdata','panel.rds'))
```

Let’s take a look at the panel. The panel is a named list of sparse
matrices (type `"dgCMatrix"`).


```r
str(panel, 1)
```

```
## List of 4
##  $ MantonBM1_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
##  $ MantonBM2_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
##  $ MantonCB1_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
##  $ MantonCB2_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
```

Before we continue, it is very important to make sure that cells in our
panel are uniquely named. No two cells (even in different samples)
should be named identically. In this case, the cells have been prefixed
by sample id, so there will not be any collisions. However, in most cases
you will have to prefix the cells before continuing.


```r
head(colnames(panel[[1]]))
```

```
## [1] "MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1"
## [2] "MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1"
## [3] "MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1"
## [4] "MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1"
## [5] "MantonBM1_HiSeq_1-TATTACCCAAAGGAAG-1"
## [6] "MantonBM1_HiSeq_1-CGCCAAGCATCTGGTA-1"
```


To quickly check that the cell names are unique, we can run:


```r
any(duplicated(unlist(lapply(panel,colnames))))
```

```
## [1] FALSE
```


Conos is focused on integration, and relies on either
[pagoda2](https://github.com/hms-dbmi/pagoda2) or
[Seurat](https://satijalab.org/seurat/) to perform dataset
pre-processing.

## Pre-processing with Pagoda2

We will generate pagoda2 objects for poorly-expressed genes from each
individual sample using the `basicP2proc` helper function for quick
processing. As the datasets will be compared to each other, we will turn
off automated dropping of low-expressed genes (using `min.cells.per.gene=0`),
and lower the numbers of local principal components (PCs) estimated for faster processing.

(**Note:** You could run the outer loop in parallel using `mclapply`, however
if executed within RStudio this sometimes causes multithreading problems.
Also, multiprocessing must be disabled in order to obtain exactly the
same individual sample embeddings from one run to another: this can be
done by using `set.seed(1)` and specifying `n.cores=1` in the command
below.)


```r
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
```

```
## 3000cells,33694genes; normalizing ...
```

```
## using plain model
```

```
## winsorizing ...
```

```
## log scale ...
```

```
## done.
```

```
## calculating variance fit ...
```

```
##  using gam
```

```
## 172overdispersed genes ...172
```

```
## persisting ...
```

```
## done.
```

```
## running PCA using 2000 OD genes .
```

```
## .
## .
## .
```

```
##  done
```

```
## creating space of type angular done
## adding data ... done
## building index ... done
## querying ... done
```

```
## running tSNE using4cores:
```

```
## 3000cells,33694genes; normalizing ...
```

```
## using plain model
```

```
## winsorizing ...
```

```
## log scale ...
```

```
## done.
```

```
## calculating variance fit ...
```

```
##  using gam
```

```
## 159overdispersed genes ...159
```

```
## persisting ...
```

```
## done.
```

```
## running PCA using 2000 OD genes .
```

```
## .
## .
## .
```

```
##  done
```

```
## creating space of type angular done
## adding data ... done
## building index ... done
## querying ... done
```

```
## running tSNE using4cores:
```

```
## 3000cells,33694genes; normalizing ...
```

```
## using plain model
```

```
## winsorizing ...
```

```
## log scale ...
```

```
## done.
```

```
## calculating variance fit ...
```

```
##  using gam
```

```
## 251overdispersed genes ...251
```

```
## persisting ...
```

```
## done.
```

```
## running PCA using 2000 OD genes .
```

```
## .
## .
## .
```

```
##  done
```

```
## creating space of type angular done
## adding data ... done
## building index ... done
## querying ... done
```

```
## running tSNE using4cores:
```

```
## 3000cells,33694genes; normalizing ...
```

```
## using plain model
```

```
## winsorizing ...
```

```
## log scale ...
```

```
## done.
```

```
## calculating variance fit ...
```

```
##  using gam
```

```
## 168overdispersed genes ...168
```

```
## persisting ...
```

```
## done.
```

```
## running PCA using 2000 OD genes .
```

```
## .
## .
## .
```

```
##  done
```

```
## creating space of type angular done
## adding data ... done
## building index ... done
## querying ... done
```

```
## running tSNE using4cores:
```

Let’s look at the output of our processing: we now have a named list of
pagoda2 objects, which is the starting point for the analysis with
Conos.


```r
str(panel.preprocessed, 1)
```

```
## List of 4
##  $ MantonBM1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
##  $ MantonBM2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
##  $ MantonCB1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
##  $ MantonCB2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
```

## Pre-processing with Seurat

Alternatively with Seurat, pre-processing can be done in a similar way
using an analogous `basicSeuratProc` helper function. If
you already have a set of Seurat objects (one per dataset), you can just
skip this step and feed them directly to `Conos$new()` as shown below.


```r
library(Seurat)
panel.preprocessed <- lapply(panel, basicSeuratProc)
```

We note that sample pre-processing steps can be used to filter/adjust
the data in custom ways. For instance, one can reduce the impact of the
cell cycle contributions by omitting cycle-annotated genes from the
matrices prior to the pre-processing. Similarly, if it is deemed
appropriate, one can regress out certain signatures using [standard
techniques](https://satijalab.org/seurat/v3.0/cell_cycle_vignette.html#regress-out-cell-cycle-scores-during-data-scaling). Please 
see the Seurat documentation for more [details](https://satijalab.org/seurat/).

# Integrating Datasets with Conos

We will now construct a Conos object for this panel of samples. At this
point we haven’t calculated anything: we have just generated an object
that contains the samples. At this step, we also set the
`n.cores` parameter. Because the graph generation with Conos can take advantage of
parallel processing, feel free to use as many physical cores as you have available
here.


```r
con <- Conos$new(panel.preprocessed, n.cores=4)
```

Our original pagoda2 (or Seurat) objects are now saved in the Conos
object (if you are short of memory you can go ahead and delete the
originals).


```r
str(con$samples,1)
```

```
## List of 4
##  $ MantonBM1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
##  $ MantonBM2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
##  $ MantonCB1_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
##  $ MantonCB2_HiSeq_1:Reference class 'Pagoda2' [package "pagoda2"] with 16 fields
##   ..and 35 methods, of which 21 are  possibly relevant
```

We can now plot a panel of these samples using the clusters we have
identified by examining each sample on its own. Please note that each sample
has an independent set of clusters that bears no relation to
clusters in other samples. For example, notice the presence (and lack thereof) of cluster
9.


```r
con$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)


Next we will build the joint graph that encompasses all the samples. We
do this by pairwise projecting samples onto a common space and
establishing the k-nearest neighbors (kNN) of mutual nearest neighbor (mNN) pairs between the samples. We then append within-sample k-nearest neighbors to the graph to ensure that all of the cells
are included in the graph:

  - We use ‘PCA’ space here which is very fast and will yield good
    integration in most cases.
  - CPCA space should provide more accurate alignment under greater
    dataset-specific distortions.
  - CCA space optimizes conservation of correlation between datasets and
    can give yield very good alignments in low-similarity cases
    (e.g. large evolutionary distances).
  - If your datasets were all measured on the same platform you may also
    want to consider “genes” space which can give better resolution in
    such (simpler) cases.

The other parameters passed to the `buildGraph()` function below are all
default values, but are included for clarity:


```r
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)
```

```
## found 0 out of 6 cached PCA space pairs ...
```

```
## running 6 additional PCA space pairs
```

```
##  done
```

```
## inter-sample links using mNN
```

```
##  done
```

```
## local pairs
```

```
##  done
```

```
## building graph .
```

```
## .
```

```
## done
```


**Note:** As pairwise comparisons may take a while, Conos will cache results
for each space. If you wish to recalculate PCA (as an example) using pairings
with different set of parameters (e.g. more components, different number
of starting over-dispersed genes, etc.), clear the cache first by doing
`con$pairs$PCA <- NULL`.

In the `$buildGraph()` invocation above, we specified
`score.component.variance=TRUE` which estimates the amount of variance
explained by successive PCs (by default this option is off to save
time). We can visualize the results using:


```r
plotComponentVariance(con, space='PCA')
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)


When using the ‘angular’ distance measure (default), it is NOT recommended to
reduce the number of components to a bare minimum indicated by the
“elbow” inflection point----rather, please include 10-20 more (typically 30 components
work well). For the ‘L2’ distance, using fewer components (i.e. at ‘elbow’
value) is sometimes better. (**NOTE:** Remember that if you want to
recalculate projections, clear the cache for that space as detailed above, i.e.
`con$pairs$PCA <- NULL`.)

We next use the graph we identified to get the global clusters. Here we use the
Leiden community detection method to obtain clusters. Increasing the
value of the resolution parameter will result in more fine-grained
clusters, while decreasing it will return coarser clustering.


```r
con$findCommunities(method=leiden.community, resolution=1)
```

## Visualization

We can now plot the clusters we obtained. Note that the number of clusters
between different samples now correspond to the same cell type. 


```r
con$plotPanel(font.size=4)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)


The convenience function `plotClusterBarplots` can be used to examine the composition of the
clusters in terms of samples (top), sample entropy (middle), and cluster size
(bottom):


```r
plotClusterBarplots(con, legend.height = 0.1)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)


Next we can check the expression pattern of a specific gene across all the individual
embeddings. In this case, we investigate the expression pattern of [GZMK](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GZMK):


```r
con$plotPanel(gene = 'GZMK')
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)


Next we embed and visualize the complete joint graph:

**Note:** An embedding estimation will run the first time around. Please see the
`$embedGraph` function for additional embedding options.

Also, both functions `$plotGraph` and `$plotPanel` are constructed off of the
main function `sccore::embeddingPlot` and will pass all visualization parameters
to this main function. So, to get full list of the possible parameters please refer to
`?sccore::embeddingPlot` and the examples below.


```r
con$plotGraph(alpha=0.1)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)


Observe that the graph captures the population structure irrespective
of the sample of origin for each cell:


```r
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)



We can also visualize gene expression on this joint graph embedding, again using "GMZK" as an example:


```r
con$plotGraph(gene='GZMK', title='GZMK expression')
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)



Other community detection methods can provide a more sensitive and
hierarchical view of the subpopulation structure. Here we run the [igraph walktrap
community](https://www.rdocumentation.org/packages/igraph/versions/0.5.1/topics/walktrap.community) detection method on the same joint graph:


```r
con$findCommunities(method = igraph::walktrap.community, steps=7)
```

**Note:** We recommend using a higher number of steps (e.g. 8-10,
though these calculations take much longer). Here we’ll get a lot of smaller clusters. 

**Note:** Different clustering results are kept as a simple list under `con$clusters`.

Now let's visualize these new clusters:


```r
con$plotPanel(clustering='walktrap', font.size=4)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

And here is the new clustering, as viewed on a joint graph:


```r
con$plotGraph(clustering='walktrap')
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)



## Changing embedding parameters

Conos is currently able to use two methods of graph embedding:
[largeVis](https://github.com/lferry007/LargeVis) (default) and
[UMAP](https://github.com/jlmelville/uwot). The UMAP embedding takes a bit longer
to estimate, but will generally give a better quality of the embedding, i.e.
sometimes UMAP will distinguish the slightest difference (which is not detected by
either largeVis or even clustering algorithms). It is best to examine both types of embeddings.

### largeVis

For the description of largeVis parameters, please look at the
`conos::projectKNNs` function. The most influential are `alpha` and
`sgd_batches`. Decreasing alpha results in less compressed clusters, and
increasing `sgd_batches` often helps to avoid cluster intersections and the
spreading out of clusters. Here we take `alpha` to a very low value, for the
sake of example:


```r
con$embedGraph(alpha=0.001, sgd_batched=1e8)  
```


```r
con$plotGraph(clustering='walktrap', size=0.1)
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)

### UMAP

The UMAP embedding supports all parameters, as described in the
[uwot](https://github.com/jlmelville/uwot) package. The two most important
ones are `spread` and `min.dist`, which together control how tight the
clusters are. According to the [python
manual](https://umap-learn.readthedocs.io/en/latest/api.html):

>   - **min.dist:** The effective minimum distance between embedded
>     points. Smaller values will result in a more clustered/clumped
>     embedding where nearby points on the manifold are drawn closer
>     together, while larger values will result on a more even dispersal
>     of points. The value should be set relative to the spread value,
>     which determines the scale at which embedded points will be spread
>     out.
>   - **spread:** The effective scale of embedded points. In combination
>     with min\_dist this determines how clustered/clumped the embedded
>     points are.
There is also a parameter responsible for the trade-off between performance
and accuracy: 

> - **min.prob.lower:** minimal probability of hitting a neighbor, after which the random walk stops. Default: 1e-7.

```r
con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=4, min.prob.lower=1e-3)
```

```
## Estimating hitting distances: 10:48:16.
## Done.
## Estimating commute distances: 10:48:25.
## Hashing adjacency list: 10:48:25.
## Done.
## Estimating distances: 10:48:27.
## Done
## Done.
## All done!: 10:48:33.
```



```r
con$plotGraph(clustering='walktrap', size=0.1)
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27-1.png)


In the example above, the UMAP layout distinguishes many of the very small
subpopulations called by walktrap apparent.

### plotPanel with common embedding

Now we can use this common embedding in `plotPanel` as well:


```r
con$plotPanel(clustering='walktrap', size=0.1, use.common.embedding=TRUE)
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28-1.png)


# Exploring Hierarchical Community Structure

## Using code

Walktrap clustering generates a hierarchical community structure. Let's being by taking a cut of the top dendrogram and visualizing it. Here we’ll take the 40 top clusters.


```r
fc <- greedyModularityCut(con$clusters$walktrap$result, 40)
```

The cut determines a finer clustering (likely overclustering) of the
dataset on its leafs:


```r
con$plotGraph(groups=fc$groups, size=0.1)
```

![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-30-1.png)


Let’s look at the hierarchical structure of these
clusters:


```r
# fc$hc is an hclust structure ... here we will convert it to a dendrogram
dend <- as.dendrogram(fc$hc)
plot(dend)
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31-1.png)



We can modify the dendrogram to show various properties. For instance,
we can alter the width of the edges to reflect how many samples are
contributing to it (normalized entropy). To do so, let’s first define a
factor specifying which samples different samples came from:


```r
samf <- con$getDatasetPerCell()
str(samf)
```

```
##  Factor w/ 4 levels "MantonBM1_HiSeq_1",..: 1 1 1 1 1 1 1 1 1 1 ...
##  - attr(*, "names")= chr [1:12000] "MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1" "MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1" "MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1" "MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1" ...
```

Now we’ll use `dendSetWidthByBreadth()` function to calculate the
entropies of each edge and set the width
accordingly:


```r
dend <- dendSetWidthByBreadth(dend, samf, fc$leafContent, min.width=1, max.width=4)
plot(dend)
```

![plot of chunk unnamed-chunk-33](figure/unnamed-chunk-33-1.png)


Similarly, we can find a factor that labels the cells by the respective tissue from which they originate (in this case BM or CB). To define a factor for this simple
dataset, we’ll simply parse the cell
names:


```r
tissue.factor <- as.factor(setNames(ifelse(grepl('BM',names(samf)),'BM','CB'), names(samf)))
str(tissue.factor)
```

```
##  Factor w/ 2 levels "BM","CB": 1 1 1 1 1 1 1 1 1 1 ...
##  - attr(*, "names")= chr [1:12000] "MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1" "MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1" "MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1" "MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1" ...
```

Now, we can color the dendrogram edges according to the tissue mixture, resulting in a more informative plot:


```r
dend <- dendSetColorByMixture(dend, tissue.factor, fc$leafContent)
plot(dend)
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35-1.png)

## Using Shiny Application

An alternative way to explore this the hierarchical community structure is by
using an interactive app, currently found at https://github.com/kharchenkolab/conosViz. The app also allows users to visualize tissue
composition and sample similarities:



```r
library(conosViz)
conosViz::conosShinyApp(con, N=30)
```


# Label Propagation

One of the uses of this graph is to propagate labels. For example, in
some cases we will only have information about the cell types in one of
the samples and we will want to automatically label the other samples.

We’ll load the annotation from a simple text file (first column giving the cell
name, second giving the cell type), and make a named factor out of
it:


```r
cellannot <- read.table(file.path(find.package('conos'),'extdata','cellannot.txt'),header=FALSE,sep='\t')
cellannot <- setNames(cellannot[,2], cellannot[,1])
```

Next we plot our panel with the annotations we made. This is to verify
that the annotated cells are indeed in only one sample and that the
other samples are unlabelled.



```r
con$plotPanel(groups = cellannot)
```

![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-38-1.png)

Next let’s propagate the labels from the one annotated sample to the
other samples.


```r
new.label.info <- con$propagateLabels(labels = cellannot, verbose=TRUE)
```

```
## Stop after 23 iterations. Norm: 0.0240661
## Min weight: 1.67017e-05, max weight: 0.367879, fading: (10, 0.1)
```

This function returns probabilities, uncertainty scores, and final labels
in the dataset of each cell belonging to each
group:


```r
con$plotPanel(colors=new.label.info$uncertainty, show.legend=TRUE, legend.title="Uncertainty", legend.pos=c(1, 0))
```

![plot of chunk unnamed-chunk-40](figure/unnamed-chunk-40-1.png)

```r
con$plotPanel(groups=new.label.info$labels, show.legend=FALSE)
```

![plot of chunk unnamed-chunk-40](figure/unnamed-chunk-40-2.png)




```r
head(new.label.info$label.distribution)
```

```
##                                        T CD4-CD8-  progenitors    B cells
## MantonBM1_HiSeq_1-GGAACTTCACTGTCGG-1 0.000000e+00 0.000000e+00 0.0000e+00
## MantonBM2_HiSeq_1-CTGATAGAGCGTTCCG-1 2.537335e-06 3.334794e-09 8.6071e-11
## MantonBM1_HiSeq_1-ACTGATGGTGGTGTAG-1 0.000000e+00 0.000000e+00 0.0000e+00
## MantonBM1_HiSeq_1-GGACATTTCCAAACTG-1 0.000000e+00 0.000000e+00 0.0000e+00
## MantonBM1_HiSeq_1-TCATTACAGACAAAGG-1 0.000000e+00 0.000000e+00 0.0000e+00
## MantonBM1_HiSeq_1-GATCGCGGTTGATTCG-1 0.000000e+00 0.000000e+00 0.0000e+00
##                                                NK    T cyto    monocytes
## MantonBM1_HiSeq_1-GGAACTTCACTGTCGG-1 0.0000000000 1.0000000 0.000000e+00
## MantonBM2_HiSeq_1-CTGATAGAGCGTTCCG-1 0.0003790042 0.9996182 6.741504e-14
## MantonBM1_HiSeq_1-ACTGATGGTGGTGTAG-1 0.0000000000 1.0000000 0.000000e+00
## MantonBM1_HiSeq_1-GGACATTTCCAAACTG-1 1.0000000000 0.0000000 0.000000e+00
## MantonBM1_HiSeq_1-TCATTACAGACAAAGG-1 0.0000000000 1.0000000 0.000000e+00
## MantonBM1_HiSeq_1-GATCGCGGTTGATTCG-1 0.0000000000 1.0000000 0.000000e+00
##                                      monomyelocytes plasma cells  dying cells
## MantonBM1_HiSeq_1-GGAACTTCACTGTCGG-1   0.000000e+00 0.000000e+00 0.000000e+00
## MantonBM2_HiSeq_1-CTGATAGAGCGTTCCG-1   4.608114e-08 2.171312e-11 1.738327e-07
## MantonBM1_HiSeq_1-ACTGATGGTGGTGTAG-1   0.000000e+00 0.000000e+00 0.000000e+00
## MantonBM1_HiSeq_1-GGACATTTCCAAACTG-1   0.000000e+00 0.000000e+00 0.000000e+00
## MantonBM1_HiSeq_1-TCATTACAGACAAAGG-1   0.000000e+00 0.000000e+00 0.000000e+00
## MantonBM1_HiSeq_1-GATCGCGGTTGATTCG-1   0.000000e+00 0.000000e+00 0.000000e+00
##                                         erythroid         HSC          pDC
## MantonBM1_HiSeq_1-GGAACTTCACTGTCGG-1 0.000000e+00 0.00000e+00 0.000000e+00
## MantonBM2_HiSeq_1-CTGATAGAGCGTTCCG-1 5.069933e-10 1.12714e-10 1.092946e-12
## MantonBM1_HiSeq_1-ACTGATGGTGGTGTAG-1 0.000000e+00 0.00000e+00 0.000000e+00
## MantonBM1_HiSeq_1-GGACATTTCCAAACTG-1 0.000000e+00 0.00000e+00 0.000000e+00
## MantonBM1_HiSeq_1-TCATTACAGACAAAGG-1 0.000000e+00 0.00000e+00 0.000000e+00
## MantonBM1_HiSeq_1-GATCGCGGTTGATTCG-1 0.000000e+00 0.00000e+00 0.000000e+00
##                                                DC
## MantonBM1_HiSeq_1-GGAACTTCACTGTCGG-1 0.000000e+00
## MantonBM2_HiSeq_1-CTGATAGAGCGTTCCG-1 1.865171e-14
## MantonBM1_HiSeq_1-ACTGATGGTGGTGTAG-1 0.000000e+00
## MantonBM1_HiSeq_1-GGACATTTCCAAACTG-1 0.000000e+00
## MantonBM1_HiSeq_1-TCATTACAGACAAAGG-1 0.000000e+00
## MantonBM1_HiSeq_1-GATCGCGGTTGATTCG-1 0.000000e+00
```

# Differential Expression

## Cluster markers

The first step we can do to understand meaning of the dataset is to look
at the cluster cell markers:


```r
new.annot <- new.label.info$labels
de.info <- con$getDifferentialGenes(groups=new.annot, n.cores=4, append.auc=TRUE)
```

```
## Estimating marker genes per sample
```

```
## Aggregating marker genes
```

```
## Estimating specificity metrics
```

```
## All done!
```

```r
head(de.info$`B cells`)
```

```
##              Gene        M        Z        PValue          PAdj       AUC
## CD74         CD74 1.789081 30.80356 3.620609e-207 1.219928e-202 0.7299556
## HLA-DRA   HLA-DRA 1.979156 28.95478 3.541451e-183 1.193221e-178 0.8681325
## HLA-DPA1 HLA-DPA1 2.170884 27.27370 1.188189e-162 4.003247e-158 0.8710053
## CD79A       CD79A 2.368081 27.24539 2.570550e-162 8.660439e-158 0.9068872
## HLA-DPB1 HLA-DPB1 2.110922 27.03211 8.391264e-160 2.827017e-155 0.8721696
## HLA-DRB1 HLA-DRB1 1.948719 26.18603 5.026161e-150 1.693263e-145 0.8557757
##          Specificity Precision ExpressionFraction
## CD74       0.4629829 0.2442885          0.9971910
## HLA-DRA    0.7625913 0.4158677          0.9747191
## HLA-DPA1   0.8344894 0.4853206          0.9101124
## CD79A      0.9169311 0.6492083          0.8983146
## HLA-DPB1   0.8230977 0.4729574          0.9235955
## HLA-DRB1   0.8016386 0.4412164          0.9129213
```



```r
cowplot::plot_grid(con$plotGraph(groups=new.annot), con$plotGraph(gene="CD74"))
```

![plot of chunk unnamed-chunk-43](figure/unnamed-chunk-43-1.png)



In addition, `getDifferentialGenes` estimates
[specificity](https://en.wikipedia.org/wiki/Sensitivity_and_specificity),
[precision](https://en.wikipedia.org/wiki/Precision_and_recall) and
expression fraction (sum expression of the gene within the cluster
divided by the total expression of this gene). If the `append.auc` flag is
set, it can estimate [ROC
AUC](https://en.wikipedia.org/wiki/Receiver_operating_characteristic#Area_under_the_curve),
but it can take some time. To find the most meaningful markers, it’s
recommended to filter the data by some lower value for the AUC and then
order the results by Z-score or
precision.


```r
de.info$monocytes %>% filter(AUC > 0.75) %>% arrange(-Precision) %>% head()
```

```
##                        Gene        M        Z        PValue          PAdj
## CD14                   CD14 3.226999 15.42978  7.996969e-53  2.685462e-48
## SERPINA1           SERPINA1 3.218886 21.05319  2.255639e-97  7.590001e-93
## RAB31                 RAB31 3.074515 13.79722  1.836815e-42  6.163985e-38
## CSTA                   CSTA 3.115659 23.84884 1.243703e-124 4.188294e-120
## FCN1                   FCN1 3.164729 26.65963 1.847622e-155 6.223716e-151
## RP11-1143G9.4 RP11-1143G9.4 2.621482 18.49146  2.243326e-75  7.541612e-71
##                     AUC Specificity Precision ExpressionFraction
## CD14          0.7691389   0.9915132 0.9076621          0.5473934
## SERPINA1      0.8829085   0.9832366 0.8807462          0.7831754
## RAB31         0.7834741   0.9844757 0.8520761          0.5835308
## CSTA          0.9103738   0.9759765 0.8488995          0.8453791
## FCN1          0.9508512   0.9707518 0.8375067          0.9312796
## RP11-1143G9.4 0.8535018   0.9764050 0.8294157          0.7316351
```


```r
con$plotGraph(gene="CD14")
```

![plot of chunk unnamed-chunk-45](figure/unnamed-chunk-45-1.png)


Or we can plot a heatmap of the top genes (top by AUC, by
default)


```r
plotDEheatmap(con,as.factor(new.annot),de.info, n.genes.per.cluster = 5, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7)
```

![plot of chunk unnamed-chunk-46](figure/unnamed-chunk-46-1.png)



Here w make a smaller heatmap, selecting a subset of cell types and showing
only a hand-picked set of
genes:


```r
gns <- c("GZMB","IL32","CD3E","LYZ","HLA-DRA","IGHD","GNLY","IGHM","GZMK")
plotDEheatmap(con,new.annot,de.info[-c(3,10)], n.genes.per.cluster = 30, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7, labeled.gene.subset = gns)
```

![plot of chunk unnamed-chunk-47](figure/unnamed-chunk-47-1.png)



## Differential expression between sample groups

Next, given a joint clustering of cells that captures the cell relationships
between samples, we can want to ask what is different between the cells
of these populations between specific samples types (in this case, between CB and
BM samples). Conos provides routines for users to do that.

The general approach we suggest for differential expression analysis is
to first pool all the data associated with each cluster (forming a
meta-cell that is analogous bulk RNA-seq measurement of the cells within
each cluster), and then use standard differential expression packages (such as DESeq2 or limma) to compare these “bulk-like” meta-cell samples,
using appropriate design models. In this section we show a convenience
routine called `getPerCellTypeDE` that enables one type of comparison (same
cluster, between sample groups); if however more advanced models are desired
(e.g. additional model variables, etc.), the `getClusterCountMatrices`
command can be used to obtain the meta-cell counts:


```r
str(con$getClusterCountMatrices(), 1)
```

```
## List of 4
##  $ MantonBM1_HiSeq_1: num [1:33694, 1:12] 0 0 0 1 0 0 0 0 46 6 ...
##   ..- attr(*, "dimnames")=List of 2
##  $ MantonBM2_HiSeq_1: num [1:33694, 1:12] 0 0 0 0 0 0 0 0 62 4 ...
##   ..- attr(*, "dimnames")=List of 2
##  $ MantonCB1_HiSeq_1: num [1:33694, 1:12] 0 0 0 0 0 0 0 0 82 7 ...
##   ..- attr(*, "dimnames")=List of 2
##  $ MantonCB2_HiSeq_1: num [1:33694, 1:12] 0 0 0 0 0 0 0 0 150 18 ...
##   ..- attr(*, "dimnames")=List of 2
```

The list above returns a pooled count matrix for each sample, where the
rows are genes and the columns are clusters. A different value for the `groups` parameter can
be supplied.

Back to DE analysis of the cluster states between groups of samples:
First we need to define our sample groups


```r
samplegroups <- list(
  bm = c("MantonBM1_HiSeq_1","MantonBM2_HiSeq_1"),
  cb = c("MantonCB1_HiSeq_1","MantonCB2_HiSeq_1")
)
```

### Simple run

We can then run differential expression between cells in these sample
groups:


```r
de.info <- getPerCellTypeDE(con, groups=as.factor(new.annot), sample.groups = samplegroups, ref.level='bm', n.cores=4)
```

```
## Error in validatePerCellTypeParams(con.obj, groups, sample.groups, ref.level, : You have to install DESeq2 package to use differential expression
```

…and examine the output:


```r
str(de.info[1:3], 2)
```

```
## List of 3
##  $ B cells    :'data.frame':	416 obs. of  9 variables:
##   ..$ Gene              : chr [1:416] "CD74" "HLA-DRA" "HLA-DPA1" "CD79A" ...
##   ..$ M                 : num [1:416] 1.79 1.98 2.17 2.37 2.11 ...
##   ..$ Z                 : num [1:416] 30.8 29 27.3 27.2 27 ...
##   ..$ PValue            : num [1:416] 3.62e-207 3.54e-183 1.19e-162 2.57e-162 8.39e-160 ...
##   ..$ PAdj              : num [1:416] 1.22e-202 1.19e-178 4.00e-158 8.66e-158 2.83e-155 ...
##   ..$ AUC               : num [1:416] 0.73 0.868 0.871 0.907 0.872 ...
##   ..$ Specificity       : num [1:416] 0.463 0.763 0.834 0.917 0.823 ...
##   ..$ Precision         : num [1:416] 0.244 0.416 0.485 0.649 0.473 ...
##   ..$ ExpressionFraction: num [1:416] 0.997 0.975 0.91 0.898 0.924 ...
##  $ DC         :'data.frame':	14 obs. of  9 variables:
##   ..$ Gene              : chr [1:14] "CST3" "HLA-DRB1" "HLA-DRA" "HLA-DPB1" ...
##   ..$ M                 : num [1:14] 3.13 2.4 1.99 2.34 2.41 ...
##   ..$ Z                 : num [1:14] 4.15 3.87 3.69 3.67 3.65 ...
##   ..$ PValue            : num [1:14] 7.37e-05 2.23e-04 4.33e-04 4.78e-04 5.06e-04 ...
##   ..$ PAdj              : num [1:14] 1 1 1 1 1 1 1 1 1 1 ...
##   ..$ AUC               : num [1:14] 0.895 0.821 0.81 0.821 0.827 ...
##   ..$ Specificity       : num [1:14] 0.791 0.696 0.655 0.713 0.725 ...
##   ..$ Precision         : num [1:14] 0.0219 0.0144 0.0129 0.015 0.0156 ...
##   ..$ ExpressionFraction: num [1:14] 1 0.946 0.964 0.929 0.929 ...
##  $ dying cells:'data.frame':	331 obs. of  9 variables:
##   ..$ Gene              : chr [1:331] "MALAT1" "MT-ND2" "MT-CYB" "MT-CO3" ...
##   ..$ M                 : num [1:331] 0.367 0.666 0.727 0.541 0.672 ...
##   ..$ Z                 : num [1:331] 15.6 12.7 11.4 11.2 10.7 ...
##   ..$ PValue            : num [1:331] 1.13e-53 4.31e-36 2.20e-29 3.15e-28 5.05e-26 ...
##   ..$ PAdj              : num [1:331] 3.79e-49 1.45e-31 7.40e-25 1.06e-23 1.70e-21 ...
##   ..$ AUC               : num [1:331] 0.501 0.504 0.512 0.505 0.51 ...
##   ..$ Specificity       : num [1:331] 0.00789 0.02093 0.05138 0.01455 0.03695 ...
##   ..$ Precision         : num [1:331] 0.0493 0.0496 0.0504 0.0497 0.0501 ...
##   ..$ ExpressionFraction: num [1:331] 0.995 0.988 0.975 0.997 0.983 ...
```


Let’s look at the results for the B cells:


```r
res <- de.info[['B cells']]$res
head(res[order(res$padj,decreasing = FALSE),])
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'head': argument 1 is not a vector
```

# Forcing Better Alignment

As can be seen from the sample distribution plot, different samples (in
particular, those representing different tissues, i.e. BM or CB in our case) form separate
subclusters within the clusters of major cell types. Conos allows users to
force better alignment through i) adjustment of the `alignment.strength
parameter`, and ii) through rebalancing of edge weights based on a
specific factor (e.g. tissue to which the cell belongs) using the
`balance.edge.weights`
parameter.



```r
con$buildGraph(k=15, k.self=5, alignment.strength=0.3, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)
```

```
## found 6 out of 6 cached PCA space pairs ...
```

```
##  done
```

```
## inter-sample links using mNN
```

```
##  done
```

```
## local pairs
```

```
##  done
```

```
## building graph .
```

```
## .
```

```
## done
```


We can re-generate the embedding and visualize the sample distribution
again:


```r
con$embedGraph()
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)
```

![plot of chunk unnamed-chunk-54](figure/unnamed-chunk-54-1.png)


We can also check the entropy, as described above:


```r
con$findCommunities()
plotClusterBarplots(con, legend.height = 0.1)
```

![plot of chunk unnamed-chunk-55](figure/unnamed-chunk-55-1.png)

# View Conos in Pagoda2 Application

Users may also interactively explore Conos objects in the [Pagoda2](https://github.com/kharchenkolab/pagoda2) application. The process is very similar to the Pagoda2 walkthrough. 

After constructing the `con` object as shown above, users can save to a serialized `*.bin` file and upload into the pagoda application with the `p2app4conos()` function:


```r
library(pagoda2)
p2app = p2app4conos(conos=con) ## users may also save the *.bin
```

```
## 12000cells,21532genes; normalizing ...
```

```
## using plain model
```

```
## winsorizing ...
```

```
## log scale ...
```

```
## done.
```

```
## calculating variance fit ...
```

```
##  using gam
```

```
## 426overdispersed genes ...426
```

```
## persisting ...
```

```
## done.
```

```
## running PCA using 3000 OD genes .
```

```
## .
## .
## .
```

```
##  done
```

```
## creating space of type angular done
## adding data ... done
## building index ... done
## querying ... done
```

```
## running tSNE using2cores:
```

```
## creating space of type angular done
## adding data ... done
## building index ... done
## querying ... done
## Writing... matsparse
## 		p array size: 21533 [First entry value: 0]
## 		i array size: 10727910 [First entry value: 2675]
## 		x array size: 10727910 [First entry value: 0.261147]
## 	Sparse matrix header information
## 		dim1=12000
## 		dim2=21532
## 		pStartOffset=32
## 		iStartOffset=86164
## 		xStartOffset=42997804
## 		dimnames1StartOffset=85909444
## 		dimnames2StartOffset=86377446
## 		dimnames2EndOffset=86592656
## Writing... mataspect
## 		p array size: 17 [First entry value: 0]
## 		i array size: 50136 [First entry value: 0]
## 		x array size: 50136 [First entry value: -1.89612]
## 	Sparse matrix header information
## 		dim1=12000
## 		dim2=16
## 		pStartOffset=32
## 		iStartOffset=100
## 		xStartOffset=200644
## 		dimnames1StartOffset=401188
## 		dimnames2StartOffset=869190
## 		dimnames2EndOffset=869468
## Writing... sparseMatrixTransp
## 		p array size: 12001 [First entry value: 0]
## 		i array size: 10727910 [First entry value: 17]
## 		x array size: 10727910 [First entry value: 0.306041]
## 	Sparse matrix header information
## 		dim1=21532
## 		dim2=12000
## 		pStartOffset=32
## 		iStartOffset=48036
## 		xStartOffset=42959676
## 		dimnames1StartOffset=85871316
## 		dimnames2StartOffset=86086526
## 		dimnames2EndOffset=86554528
## Making File from payload...
## 	File format information
## 		Index entry size is 140 bytes
## 		File header size is 48 bytes
## 	Preparing header...
## 	Total index size is: 2660 bytes
## 	Constructing index...
## 		Writing entry 0 of size 205 blocks (or 6717440 bytes)
## 		Writing entry 1 of size 1 blocks (or 32768 bytes)
## 		Writing entry 2 of size 46 blocks (or 1507328 bytes)
## 		Writing entry 3 of size 15 blocks (or 491520 bytes)
## 		Writing entry 4 of size 49 blocks (or 1605632 bytes)
## 		Writing entry 5 of size 1 blocks (or 32768 bytes)
## 		Writing entry 6 of size 1 blocks (or 32768 bytes)
## 		Writing entry 7 of size 1 blocks (or 32768 bytes)
## 		Writing entry 8 of size 1 blocks (or 32768 bytes)
## 		Writing entry 9 of size 6 blocks (or 196608 bytes)
## 		Writing entry 10 of size 22 blocks (or 720896 bytes)
## 		Writing entry 11 of size 22 blocks (or 720896 bytes)
## 		Writing entry 12 of size 6 blocks (or 196608 bytes)
## 		Writing entry 13 of size 6 blocks (or 196608 bytes)
## 		Writing entry 14 of size 6 blocks (or 196608 bytes)
## 		Writing entry 15 of size 6 blocks (or 196608 bytes)
## 		Writing entry 16 of size 2643 blocks (or 86605824 bytes)
## 		Writing entry 17 of size 27 blocks (or 884736 bytes)
## 		Writing entry 18 of size 2642 blocks (or 86573056 bytes)
## Free up entry payloads
```

```r
## p2app4conos(conos=con, file="conosApp1.bin", save=TRUE)
show.app(app=p2appconos, name='conos_app')
```

```
## starting httpd help server ...
```

```
##  done
```

```
## 
## Server started on host 127.0.0.1 and port 18151 . App urls are:
## 
## 	http://127.0.0.1:18151/custom/RookTest
```

```
## Error in inherits(app, "RhttpdApp"): object 'p2appconos' not found
```












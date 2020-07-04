Conos Walkthrough
================

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

``` r
library(conos)
library(dplyr)
```

# Loading the Data

Next we will load a previously prepared panel of samples. This panel was
made up of 16 cord blood and bone marrow samples, but for convenience, we will here focus on a smaller subset of just 4 samples. All samples have been subset to a size of exactly 3000 cells. 

**Note:** When starting with your own panel, we
recommend filtering out low-count/poor-quality/dying cells, as is standard for quality control.

``` r
panel <- readRDS(file.path(find.package('conos'),'extdata','panel.rds'))
```

Let’s take a look at the panel. The panel is a named list of sparse
matrices (type `"dgCMatrix"`).

``` r
str(panel, 1)
```

    ## List of 4
    ##  $ MantonBM1_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonBM2_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonCB1_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##  $ MantonCB2_HiSeq_1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots

Before we continue, it is very important to make sure that cells in our
panel are uniquely named. No two cells (even in different samples)
should be named identically. In this case, the cells have been prefixed
by sample id, so there will not be any collisions. However, in most cases
you will have to prefix the cells before continuing.

``` r
head(colnames(panel[[1]]))
```

    ## [1] "MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1"
    ## [2] "MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1"
    ## [3] "MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1"
    ## [4] "MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1"
    ## [5] "MantonBM1_HiSeq_1-TATTACCCAAAGGAAG-1"
    ## [6] "MantonBM1_HiSeq_1-CGCCAAGCATCTGGTA-1"

To quickly check that the cell names are unique, we can run:

``` r
any(duplicated(unlist(lapply(panel,colnames))))
```

    ## [1] FALSE

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

``` r
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
```

    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 171 overdispersed genes ... 171persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:
    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 159 overdispersed genes ... 159persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:
    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 248 overdispersed genes ... 248persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:
    ## 3000 cells, 33694 genes; normalizing ... using plain model winsorizing ... log scale ... done.
    ## calculating variance fit ... using gam 166 overdispersed genes ... 166persisting ... done.
    ## running PCA using 2000 OD genes .... done
    ## running tSNE using 4 cores:

Let’s look at the output of our processing: we now have a named list of
pagoda2 objects, which is the starting point for the analysis with
Conos.

``` r
str(panel.preprocessed, 1)
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

## Pre-processing with Seurat

Alternatively with Seurat, pre-processing can be done in a similar way
using an analogous `basicSeuratProc` helper function. If
you already have a set of Seurat objects (one per dataset), you can just
skip this step and feed them directly to `Conos$new()` as shown below.

``` r
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

``` r
con <- Conos$new(panel.preprocessed, n.cores=4)
```

Our original pagoda2 (or Seurat) objects are now saved in the Conos
object (if you are short of memory you can go ahead and delete the
originals).

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

We can now plot a panel of these samples using the clusters we have
identified by examining each sample on its own. Please note that each sample
has an independent set of clusters that bears no relation to
clusters in other samples. For example, notice the presence (and lack thereof) of cluster
9.

``` r
con$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

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
    (e.g. large evolutionary distances).
  - If your datasets were all measured on the same platform you may also
    want to consider “genes” space which can give better resolution in
    such (simpler) cases.

The other parameters passed to the `buildGraph()` function below are all
default values, but are included for clarity:

``` r
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)
```

    ## found 0 out of 6 cached PCA  space pairs ... running 6 additional PCA  space pairs  done
    ## inter-sample links using  mNN   done
    ## local pairs local pairs  done
    ## building graph ..done

**Note:** As pairwise comparisons may take a while, Conos will cache results
for each space. If you wish to recalculate PCA (as an example) using pairings
with different set of parameters (e.g. more components, different number
of starting over-dispersed genes, etc.), clear the cache first by doing
`con$pairs$PCA <- NULL`.

In the `$buildGraph()` invocation above, we specified
`score.component.variance=TRUE` which estimates the amount of variance
explained by successive PCs (by default this option is off to save
time). We can visualize the results using:

``` r
plotComponentVariance(con, space='PCA')
```

![](walkthrough_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

When using the ‘angular’ distance measure (default), it is NOT recommended to
reduce the number of components to a bare minimum indicated by the
“elbow” inflection point----rather, please include 10-20 more (typically 30 components
work well). For the ‘L2’ distance, using fewer components (i.e. at ‘elbow’
value) is sometimes better. (**NOTE:** Remember that if you want to
recalculate projections, clear the cache for that space as detailed above, i.e.
`con$pairs$PCA <- NULL`.)

We next use the graph we identified to get the global clusters. Here we use the
Leiden community detection method to obtain clusters. Increasing the
value of the resolution parameter will result in more fine-grained
clusters, while decreasing it will return coarser clustering.

``` r
con$findCommunities(method=leiden.community, resolution=1)
```

## Visualization

We can now plot the clusters we obtained. Note that the number of clusters
between different samples now correspond to the same cell type. Also not
the presence of cluster 5 in BM samples only, but not in CB.

``` r
con$plotPanel(font.size=4)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

The convenience function `plotClusterBarplots` can be used to examine the composition of the
clusters in terms of samples (top), sample entropy (middle), and cluster size
(bottom):

``` r
plotClusterBarplots(con, legend.height = 0.1)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Next we can check the expression pattern of a specific gene across all the individual
embeddings. In this case, we investigate the expression pattern of [GZMK](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GZMK):

``` r
con$plotPanel(gene = 'GZMK')
```

![](walkthrough_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Next we embed and visualize the complete joint graph:

**Note:** An embedding estimation will run the first time around. Please see the
`$embedGraph` function for additional embedding options.

Also, both functions `$plotGraph` and `$plotPanel` are constructed off of the
main function `conos::embeddingPlot` and will pass all visualization parameters
to this main function. So, to get full list of the possible parameters please refer to
`?conos::embeddingPlot` and the examples below.

``` r
con$plotGraph(alpha=0.1)
```

    ## Estimating embeddings.

![](walkthrough_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Observe that the graph captures the population structure irrespective
of the sample of origin for each cell:

``` r
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

We can also visualize gene expression on this joint graph embedding, again using "GMZK" as an example:

``` r
con$plotGraph(gene='GZMK', title='GZMK expression')
```

![](walkthrough_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Other community detection methods can provide a more sensitive and
hierarchical view of the subpopulation structure. Here we run the [igraph walktrap
community detection method](https://www.rdocumentation.org/packages/igraph/versions/0.5.1/topics/walktrap.community) on the same joint graph:

``` r
con$findCommunities(method = igraph::walktrap.community, steps=7)
```

**Note:** We recommend using a higher number of steps (e.g. 8-10,
though these calculations take much longer). Here we’ll get a lot of smaller clusters. 

**Note:** Different clustering results are kept as a simple list under `con$clusters`.

Now let's visualize these new clusters:

``` r
con$plotPanel(clustering='walktrap',font.size=4)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

And here is the new clustering, as viewed on a joint graph:

``` r
con$plotGraph(clustering='walktrap')
```

![](walkthrough_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

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

``` r
con$embedGraph(alpha=0.001, sgd_batched=1e8)  
```

    ## Estimating embeddings.

``` r
con$plotGraph(clustering='walktrap', size=0.1)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

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

``` r
con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=4, min.prob.lower=1e-3)
```

    ## Convert graph to adjacency list...
    ## Done
    ## Estimate nearest neighbors and commute times...
    ## Estimating hitting distances: 17:56:45.
    ## Done.
    ## Estimating commute distances: 17:56:48.
    ## Hashing adjacency list: 17:56:48.
    ## Done.
    ## Estimating distances: 17:56:49.
    ## Done
    ## Done.
    ## All done!: 17:56:51.
    ## Done
    ## Estimate UMAP embedding...

    ## 17:56:51 Read 12000 rows and found 1 numeric columns

    ## 17:56:52 Commencing smooth kNN distance calibration using 4 threads

    ## 17:56:52 Initializing from normalized Laplacian + noise

    ## 17:56:53 Commencing optimization for 1000 epochs, with 351362 positive edges using 4 threads

    ## 17:57:09 Optimization finished

    ## Done

``` r
con$plotGraph(clustering='walktrap', size=0.1)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

In the example above, the UMAP layout distinguishes many of the very small
subpopulations called by walktrap apparent.

### plotPanel with common embedding

Now we can use this common embedding in `plotPanel` as well:

``` r
con$plotPanel(clustering='walktrap', size=0.1, use.common.embedding=TRUE)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

# Exploring Hierarchical Community Structure

## Using code

Walktrap clustering generates a hierarchical community structure. Let's being by taking a cut of the top dendrogram and visualizing it. Here we’ll take the 40 top clusters.

``` r
fc <- greedyModularityCut(con$clusters$walktrap$result, 40)
```

The cut determines a finer clustering (likely overclustering) of the
dataset on its leafs:

``` r
con$plotGraph(groups=fc$groups, size=0.1)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

Let’s look at the hierarchical structure of these
clusters:

``` r
# fc$hc is an hclust structure ... here we will convert it to a dendrogram
dend <- as.dendrogram(fc$hc)
plot(dend)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

We can modify the dendrogram to show various properties. For instance,
we can alter the width of the edges to reflect how many samples are
contributing to it (normalized entropy). To do so, let’s first define a
factor specifying which samples different samples came from:

``` r
samf <- con$getDatasetPerCell()
str(samf)
```

    ##  Factor w/ 4 levels "MantonBM1_HiSeq_1",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:12000] "MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1" "MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1" "MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1" "MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1" ...

Now we’ll use `dendSetWidthByBreadth()` function to calculate the
entropies of each edge and set the width
accordingly:

``` r
dend <- dendSetWidthByBreadth(dend, samf, fc$leafContent, min.width=1, max.width=4)
plot(dend)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Similarly, we can find a factor that labels the cells by the respective tissue from which they originate (in this case BM or CB). To define a factor for this simple
dataset, we’ll simply parse the cell
names:

``` r
tissue.factor <- as.factor(setNames(ifelse(grepl('BM',names(samf)),'BM','CB'), names(samf)))
str(tissue.factor)
```

    ##  Factor w/ 2 levels "BM","CB": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:12000] "MantonBM1_HiSeq_1-TCTATTGGTCTCTCGT-1" "MantonBM1_HiSeq_1-GAATAAGTCACGCATA-1" "MantonBM1_HiSeq_1-ACACCGGTCTAACTTC-1" "MantonBM1_HiSeq_1-TCATTTGGTACGCTGC-1" ...

Now, we can color the dendrogram edges according to the tissue mixture, resulting in a more informative plot:

``` r
dend <- dendSetColorByMixture(dend, tissue.factor, fc$leafContent)
plot(dend)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

## Using Shiny Application

An alternative way to explore this the hierarchical community structure is by
using an interactive app. The app also allows users to visualize tissue
composition and sample similarities:

``` r
conosShinyApp(con,N=30)
```

# Label Propagation

One of the uses of this graph is to propagate labels. For example, in
some cases we will only have information about the cell types in one of
the samples and we will want to automatically label the other samples.

We’ll load the annotation from a simple text file (first column giving the cell
name, second giving the cell type), and make a named factor out of
it:

``` r
cellannot <- read.table(file.path(find.package('conos'),'extdata','cellannot.txt'),header=FALSE,sep='\t')
cellannot <- setNames(cellannot[,2], cellannot[,1])
```

Next we plot our panel with the annotations we made. This is to verify
that the annotated cells are indeed in only one sample and that the
other samples are unlabelled.

``` r
con$plotPanel(groups = cellannot)
```

    ## Warning: Factor `Group` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`
    
    ## Warning: Factor `Group` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`
    
    ## Warning: Factor `Group` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`
    
    ## Warning: Factor `Group` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`
    
    ## Warning: Factor `Group` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`
    
    ## Warning: Factor `Group` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    ## Warning: Removed 1 rows containing missing values (geom_label_repel).
    
    ## Warning: Removed 1 rows containing missing values (geom_label_repel).
    
    ## Warning: Removed 1 rows containing missing values (geom_label_repel).

![](walkthrough_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

Next let’s propagate the labels from the one annotated sample to the
other samples.

``` r
new.label.info <- con$propagateLabels(labels = cellannot, verbose=TRUE)
```

This function returns probabilities, uncertainty scores, and final labels
in the dataset of each cell belonging to each
group:

``` r
con$plotPanel(colors=new.label.info$uncertainty, show.legend=TRUE, legend.title="Uncertainty", legend.pos=c(1, 0))
```

![](walkthrough_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

``` r
con$plotPanel(groups=new.label.info$labels, show.legend=FALSE)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-40-2.png)<!-- -->

``` r
head(new.label.info$label.distribution)
```

    ##                                        T CD4-CD8-  progenitors
    ## MantonBM1_HiSeq_1-ATCACGAAGGAGTCTG-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-TTTGCGCGTAGCGCTC-1 1.522970e-06 2.367580e-08
    ## MantonBM1_HiSeq_1-GCGCGATAGAAGGACA-1 0.000000e+00 0.000000e+00
    ## MantonBM1_HiSeq_1-ACAGCTATCGATAGAA-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-AGCTCCTCATGTCTCC-1 9.402861e-07 1.940159e-08
    ## MantonBM1_HiSeq_1-CCCAATCCACCATCCT-1 0.000000e+00 0.000000e+00
    ##                                           B cells           NK
    ## MantonBM1_HiSeq_1-ATCACGAAGGAGTCTG-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-TTTGCGCGTAGCGCTC-1 1.708361e-08 2.533198e-09
    ## MantonBM1_HiSeq_1-GCGCGATAGAAGGACA-1 0.000000e+00 0.000000e+00
    ## MantonBM1_HiSeq_1-ACAGCTATCGATAGAA-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-AGCTCCTCATGTCTCC-1 1.548886e-08 1.763021e-09
    ## MantonBM1_HiSeq_1-CCCAATCCACCATCCT-1 0.000000e+00 0.000000e+00
    ##                                            T cyto monocytes monomyelocytes
    ## MantonBM1_HiSeq_1-ATCACGAAGGAGTCTG-1 0.000000e+00 1.0000000   0.0000000000
    ## MantonBM2_HiSeq_1-TTTGCGCGTAGCGCTC-1 5.414639e-10 0.9982103   0.0015851989
    ## MantonBM1_HiSeq_1-GCGCGATAGAAGGACA-1 0.000000e+00 1.0000000   0.0000000000
    ## MantonBM1_HiSeq_1-ACAGCTATCGATAGAA-1 0.000000e+00 1.0000000   0.0000000000
    ## MantonBM2_HiSeq_1-AGCTCCTCATGTCTCC-1 4.431973e-10 0.9990664   0.0007616298
    ## MantonBM1_HiSeq_1-CCCAATCCACCATCCT-1 0.000000e+00 1.0000000   0.0000000000
    ##                                      plasma cells  dying cells
    ## MantonBM1_HiSeq_1-ATCACGAAGGAGTCTG-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-TTTGCGCGTAGCGCTC-1 1.879407e-07 2.435256e-07
    ## MantonBM1_HiSeq_1-GCGCGATAGAAGGACA-1 0.000000e+00 0.000000e+00
    ## MantonBM1_HiSeq_1-ACAGCTATCGATAGAA-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-AGCTCCTCATGTCTCC-1 1.324181e-07 2.126047e-07
    ## MantonBM1_HiSeq_1-CCCAATCCACCATCCT-1 0.000000e+00 0.000000e+00
    ##                                         erythroid          HSC
    ## MantonBM1_HiSeq_1-ATCACGAAGGAGTCTG-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-TTTGCGCGTAGCGCTC-1 6.017598e-09 2.270053e-09
    ## MantonBM1_HiSeq_1-GCGCGATAGAAGGACA-1 0.000000e+00 0.000000e+00
    ## MantonBM1_HiSeq_1-ACAGCTATCGATAGAA-1 0.000000e+00 0.000000e+00
    ## MantonBM2_HiSeq_1-AGCTCCTCATGTCTCC-1 4.369539e-09 1.832747e-09
    ## MantonBM1_HiSeq_1-CCCAATCCACCATCCT-1 0.000000e+00 0.000000e+00
    ##                                               pDC           DC
    ## MantonBM1_HiSeq_1-ATCACGAAGGAGTCTG-1 0.000000e+00 0.0000000000
    ## MantonBM2_HiSeq_1-TTTGCGCGTAGCGCTC-1 1.807266e-06 0.0002007056
    ## MantonBM1_HiSeq_1-GCGCGATAGAAGGACA-1 0.000000e+00 0.0000000000
    ## MantonBM1_HiSeq_1-ACAGCTATCGATAGAA-1 0.000000e+00 0.0000000000
    ## MantonBM2_HiSeq_1-AGCTCCTCATGTCTCC-1 1.606282e-06 0.0001690121
    ## MantonBM1_HiSeq_1-CCCAATCCACCATCCT-1 0.000000e+00 0.0000000000

# Differential Expression

## Cluster markers

The first step we can do to understand meaning of the dataset is to look
at the cluster cell markers:

``` r
new.annot <- new.label.info$labels
de.info <- con$getDifferentialGenes(groups=new.annot, n.cores=4, append.auc=TRUE)
```

    ## Estimating marker genes per sample
    ## Aggregating marker genes
    ## Estimating specificity metrics
    ## All done!

``` r
head(de.info$`B cells`)
```

    ##              Gene        M        Z        PValue          PAdj       AUC
    ## CD74         CD74 1.783430 30.71596 5.358033e-206 1.805336e-201 0.7285243
    ## HLA-DRA   HLA-DRA 1.970110 28.84419 8.653330e-182 2.915566e-177 0.8659041
    ## HLA-DPA1 HLA-DPA1 2.166771 27.26773 1.398515e-162 4.711877e-158 0.8699558
    ## CD79A       CD79A 2.367299 27.25606 1.922121e-162 6.475817e-158 0.9059042
    ## HLA-DPB1 HLA-DPB1 2.104391 26.93159 1.263872e-158 4.257986e-154 0.8697674
    ## HLA-DRB1 HLA-DRB1 1.946088 26.22140 1.989544e-150 6.702575e-146 0.8553614
    ##          Specificity Precision ExpressionFraction
    ## CD74       0.4631641 0.2448390          0.9944103
    ## HLA-DRA    0.7626656 0.4161074          0.9703745
    ## HLA-DPA1   0.8348111 0.4865189          0.9077697
    ## CD79A      0.9172838 0.6508323          0.8960313
    ## HLA-DPB1   0.8230977 0.4729574          0.9189491
    ## HLA-DRB1   0.8021024 0.4428455          0.9116825

``` r
cowplot::plot_grid(con$plotGraph(groups=new.annot), con$plotGraph(gene="CD74"))
```

![](walkthrough_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

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

``` r
de.info$monocytes %>% filter(AUC > 0.75) %>% arrange(-Precision) %>% head()
```

    ##       Gene        M        Z        PValue          PAdj       AUC
    ## 1     CD14 3.226926 15.41709  9.724657e-53  3.265637e-48 0.7692823
    ## 2 SERPINA1 3.220854 21.07233  1.507243e-97  5.071723e-93 0.8832630
    ## 3    RAB31 3.073102 13.76918  2.703236e-42  9.070979e-38 0.7833055
    ## 4     CSTA 3.114940 23.81856 2.559077e-124 8.617947e-120 0.9100943
    ## 5     FCN1 3.164136 26.62166 5.080555e-155 1.711385e-150 0.9506489
    ## 6     G0S2 3.177618 16.25523  1.673163e-58  5.620658e-54 0.7979838
    ##   Specificity Precision ExpressionFraction
    ## 1   0.9914237 0.9066798          0.5477745
    ## 2   0.9831445 0.8800799          0.7839763
    ## 3   0.9842970 0.8503460          0.5833828
    ## 4   0.9756997 0.8471148          0.8451039
    ## 5   0.9704726 0.8359084          0.9311573
    ## 6   0.9805675 0.8298722          0.6166172

``` r
con$plotGraph(gene="CD14")
```

![](walkthrough_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

Or we can plot a heatmap of the top genes (top by AUC, by
default)

``` r
plotDEheatmap(con,as.factor(new.annot),de.info, n.genes.per.cluster = 5, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

Here w make a smaller heatmap, selecting a subset of cell types and showing
only a hand-picked set of
genes:

``` r
gns <- c("GZMB","IL32","CD3E","LYZ","HLA-DRA","IGHD","GNLY","IGHM","GZMK")
plotDEheatmap(con,new.annot,de.info[-c(3,10)], n.genes.per.cluster = 30, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7, labeled.gene.subset = gns)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

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
(e.g. additional model variables, etc.), the `getClusterCountMatrices`
command can be used to obtain the meta-cell counts:

``` r
str( con$getClusterCountMatrices() , 1)
```

    ## List of 4
    ##  $ MantonBM1_HiSeq_1: num [1:33694, 1:12] 0 0 0 1 0 0 0 0 44 5 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##  $ MantonBM2_HiSeq_1: num [1:33694, 1:12] 0 0 0 0 0 0 0 0 66 4 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##  $ MantonCB1_HiSeq_1: num [1:33694, 1:12] 0 0 0 0 0 0 0 0 84 7 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##  $ MantonCB2_HiSeq_1: num [1:33694, 1:12] 0 0 0 0 0 0 0 0 153 20 ...
    ##   ..- attr(*, "dimnames")=List of 2

The list above returns a pooled count matrix for each sample, where the
rows are genes and the columns are clusters. A different value for the `groups` parameter can
be supplied.

Back to DE analysis of the cluster states between groups of samples:
First we need to define our sample groups

``` r
samplegroups <- list(
  bm = c("MantonBM1_HiSeq_1","MantonBM2_HiSeq_1"),
  cb = c("MantonCB1_HiSeq_1","MantonCB2_HiSeq_1")
)
```

### Simple run

We can then run differential expression between cells in these sample
groups:

``` r
de.info <- getPerCellTypeDE(con, groups=as.factor(new.annot), sample.groups = samplegroups, ref.level='bm', n.cores=4)
```

…and examine the output:

``` r
str(de.info[1:3], 2)
```

    ## List of 3
    ##  $ B cells    :List of 3
    ##   ..$ res          :'data.frame':    33694 obs. of  6 variables:
    ##   ..$ cm           : num [1:33694, 1:4] 0 0 0 1 0 0 0 0 22 1 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ sample.groups:List of 2
    ##  $ DC         : logi NA
    ##  $ dying cells:List of 3
    ##   ..$ res          :'data.frame':    33694 obs. of  6 variables:
    ##   ..$ cm           : num [1:33694, 1:4] 0 0 0 0 0 0 0 0 9 1 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ sample.groups:List of 2

Let’s look at the results for the B cells:

``` r
res <- de.info[['B cells']]$res
head(res[order(res$padj,decreasing = FALSE),])
```

    ##                baseMean log2FoldChange     lfcSE      stat       pvalue
    ## RP11-386I14.4 637.71953       3.618939 0.2941141 12.304543 8.562006e-35
    ## HBA2          437.33335       3.386807 0.3697592  9.159495 5.214086e-20
    ## CH17-373J23.1 435.34473       3.383322 0.3844136  8.801255 1.352946e-18
    ## CD69          621.47862       2.479981 0.2993561  8.284387 1.187195e-16
    ## IGHA1          94.53635      -6.202638 0.7607468 -8.153354 3.539679e-16
    ## HMGB2         173.36582      -3.096187 0.3862006 -8.017043 1.083212e-15
    ##                       padj
    ## RP11-386I14.4 1.341666e-30
    ## HBA2          4.085237e-16
    ## CH17-373J23.1 7.066888e-15
    ## CD69          4.650835e-13
    ## IGHA1         1.109336e-12
    ## HMGB2         2.828989e-12

# Forcing Better Alignment

As can be seen from the sample distribution plot, different samples (in
particular, those representing different tissues, i.e. BM or CB in our case) form separate
subclusters within the clusters of major cell types. Conos allows users to
force better alignment through i) adjustment of the `alignment.strength
parameter`, and ii) through rebalancing of edge weights based on a
specific factor (e.g. tissue to which the cell belongs) using the
`balance.edge.weights`
parameter.

``` r
con$buildGraph(k=15, k.self=5, alignment.strength=0.3, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)
```

    ## found 6 out of 6 cached PCA  space pairs ...  done
    ## inter-sample links using  mNN   done
    ## local pairs local pairs  done
    ## building graph ..done

We can re-generate the embedding and visualize the sample distribution
again:

``` r
con$embedGraph()
```

    ## Estimating embeddings.

``` r
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

We can also check the entropy, as described above:

``` r
con$findCommunities()
plotClusterBarplots(con, legend.height = 0.1)
```

![](walkthrough_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

# conos
## Clustering on Network of Samples
* What is Conos? 
It's a package to wire together large collections of single-cell RNA-seq datasets. It focuses on uniform mapping of homologous cell types across heterogeneous sample collections. For instance, a collection of dozens of peripheral blood samples from cancer patients, combined with dozens of controls. And perhaps also including samples of a related tissue, such as lymph nodes.

* How does it work? 
![overview](http://pklab.med.harvard.edu/peterk/conos/Figure1_take3.pk.png)
Conos applies one of many error-prone methods to align each pair of samples in a collection, establishing weighted inter-sample cell-to-cell links, creating a global joint graph. Cells of the same type will tend to map to each other across many such pair-wise comparisons, forming cliques, that can be recognized as clusters (graph communities). 

* What does it produce?
In essense, Conos will take a large, potentially heterogeneous panel of samples and will produce clustering grouping similar cell subpopulations togehter in a way that will be robust to inter-sample variation:
![example](http://pklab.med.harvard.edu/peterk/conos/bm_uniform_labels_trim.png)

* What are the advantages over existing alignment methods? 
Conos is robust to heterogeneity of samples within collection, as well as noise. The ability to resolve finer subpopulation structure improves as the size of the panel increases.

* What do I need to run it?
Conos is an R package. Currently, it supports analysis of sample collections analyzed using [pagoda2](https://github.com/hms-dbmi/pagoda2). We will add support for Seurat collections shortly.

## Installation

```r
devtools::install_github("hms-dbmi/conos")
```

## Usage example
Please see [walkthrough](http://pklab.med.harvard.edu/peterk/conos/walkthrough.nb.html) for an example of Conos analysis.

Given a list of individual processed samples (`pl`), Conos processing can be as simple as this:
```r
# construct conos object, where pl is a list of pagoda2 objects 
con <- Conos$new(pl)

# build graph
con$buildGraph(verbose=T)

# find communities
con$findCommunities()

# plot panel with joint clustering results
con$plotPanel()
```

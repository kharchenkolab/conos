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

## Installing Conos as Docker Container
The fastest and most efficient way to get Conos on a mac or windows system is through a docker container. The docker distribution is current as of October 2018 and also includes the (Pagoda2 package)[https://github.com/hms-dbmi/pagoda2]. To start a docker container, first [install docker](https://docs.docker.com/install/) on your platform and then start the pagoda container with the following command in the shell:

```
docker run -p 8787:8787 docker.io/barkasn/pagoda2
```
The first time you run the command it will download several images so make sure that you have fast internet access setup. You can then point your browser to http://localhost:8787/ to get an Rstudio environment with pagoda2 and conos installed. Explore the docker (--mount option)[https://docs.docker.com/storage/volumes/] to allow access of the docker image to your local files.


### Native

(Install pagoda2 and it's dependencies)[https://github.com/hms-dbmi/pagoda2]
```r
devtools::install_github("hms-dbmi/conos")
```

## Usage example
Please see [walkthrough](http://pklab.med.harvard.edu/peterk/conos/walkthrough.nb.html) for an example of Conos analysis.

A more comprehensive tutorial is available [here](vignettes/Conos_Walkthrough.Rmd)

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

[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/conos.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/conos)
[![CRAN status](https://www.r-pkg.org/badges/version/conos)](https://cran.r-project.org/package=conos)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/conos)](https://cran.r-project.org/package=conos)


<img src="https://github.com/kharchenkolab/conos/blob/dev/inst/conos_logo.png" align="right" height="140">

# conos

- [Introduction](#conos-clustering-on-network-of-samples)
- [Basics of using conos](#basics-of-using-conos)
- [Tutorials](#tutorials)
  * [Conos walkthrough](#conos-walkthrough-start-here)
  * [Advanced Conos workflows](#advanced-conos-workflows)
  * [Data sources (Seurat, anndata, loom, lstar)](#data-sources-seurat-anndata-loom-lstar)
  * [Adjusting alignment strength](#adjusting-alignment-strength)
- [Installation](#installation)
  * [Running conos via Docker](#running-conos-via-docker)
- [References](#references)
  
## Conos: Clustering On Network Of Samples

* **What is conos?**
Conos is an R package to wire together large collections of single-cell RNA-seq datasets, which allows for both the identification of recurrent cell clusters and the propagation of information between datasets in multi-sample or atlas-scale collections. It focuses on the uniform mapping of homologous cell types across heterogeneous sample collections. For instance, users could investigate a collection of dozens of peripheral blood samples from cancer patients combined with dozens of controls, which perhaps includes samples of a related tissue such as lymph nodes.

* **How does it work?**
![overview](http://pklab.med.harvard.edu/peterk/conos/Figure1_take3.pk.png)
Conos applies one of many error-prone methods to align each pair of samples in a collection, establishing weighted inter-sample cell-to-cell links. The resulting joint graph can then be analyzed to identify subpopulations across different samples. Cells of the same type will tend to map to each other across many such pairwise comparisons, forming cliques that can be recognized as clusters (graph communities). 

   Conos processing can be divided into three phases:
    * **Phase 1: Filtering and normalization** Each individual dataset in the sample panel is filtered and normalized using standard packages for single-dataset processing: either `pagoda2` or `Seurat`. Specifically, Conos relies on these methods to perform cell filtering, library size normalization, identification of overdispersed genes and, in the case of pagoda2, variance normalization. (Conos is robust to variations in the normalization procedures, but it is recommended that all of the datasets be processed uniformly.)
    * **Phase 2: Identify multiple plausible inter-sample mappings** Conos performs pairwise comparisons of the datasets in the panel to establish an initial error-prone mapping between cells of different datasets. 
    * **Phase 3: Joint graph construction** These inter-sample edges from Phase 2 are then combined with lower-weight intra-sample edges during the joint graph construction. The joint graph is then used for downstream analysis, including community detection and label propagation. For a comprehensive description of the algorithm, please refer to our [publication](https://doi.org/10.1038/s41592-019-0466-z).

* **What does it produce?**
In essence, conos will take a large, potentially heterogeneous panel of samples and will produce clustering grouping similar cell subpopulations together in a way that will be robust to inter-sample variation:  
![example](http://pklab.med.harvard.edu/peterk/conos/bm_uniform_labels_trim.png)

* **What are the advantages over existing alignment methods?** 
Conos is robust to heterogeneity of samples within a collection, as well as noise. The ability to resolve finer subpopulation structure improves as the size of the panel increases.


## Basics of using conos

Each sample is first processed on its own with [pagoda2](https://github.com/kharchenkolab/pagoda2)
(or Seurat). From a named list of count matrices (`cms`), one call per sample builds a `Pagoda2` object
carrying the PCA reduction Conos aligns on:

```r
library(conos)
library(pagoda2)

samples <- lapply(cms, function(cm) Pagoda2$from(cm)$run(steps = c("variance", "pca")))
```

Then the whole joint analysis is a handful of calls on the Conos object:

```r
con <- Conos$new(samples)   # collect the samples

con$runGraph()              # align every pair of samples into one joint graph
con$runClustering()         # joint clusters (Leiden communities)
con$runEmbedding()          # joint 2-D embedding (UMAP)

con$plotGraph()             # the joint embedding (colour by cluster / sample / gene)
con$plotPanel()             # the same, faceted per sample
con$plotMarkerDotPlot()     # top marker genes per cluster
con$propagateLabels(labels = cellannot)  # transfer annotations from one sample to the rest
```

`runGraph()`, `runClustering()`, `runEmbedding()` and `runMarkers()` are the recommended verbs (the older
`buildGraph()` / `findCommunities()` / `embedGraph()` / `getDifferentialGenes()` names still work but are
deprecated). Samples can be `Pagoda2` or `Seurat` objects, or read from files — see the
[data-sources tutorial](doc/conos-data-sources.ipynb). For full documentation of the class, run `?Conos`.


## Tutorials


Please see the following tutorials for detailed examples of how to use conos.

The tutorials are rendered Jupyter notebooks — GitHub displays them directly (no download needed):

### Conos walkthrough (start here):
The standard workflow on a panel of samples — pre-processing, building the joint graph, clustering,
embedding, marker dot plots, and label transfer.
* [Jupyter notebook](doc/conos-walkthrough.ipynb)

### Advanced Conos workflows:
Choosing the alignment space (reciprocal PCA / CCA), `planIntegration()`, disk-backed memory control
(`pairs.storage`), resolution scanning, label-transfer methods, and multimodal facets.
* [Jupyter notebook](doc/conos-advanced.ipynb)

### Data sources (Seurat, anndata, loom, lstar):
Building a panel from Seurat objects in memory, from files written by other pipelines (anndata `.h5ad`,
Seurat `.h5seurat`, loom, lstar zarr — read via pagoda2's `from*()` constructors), and from a mix of object
types in one panel.
* [Jupyter notebook](doc/conos-data-sources.ipynb)

### Adjusting alignment strength:
Controlling how forcefully samples are pulled together — `alignment.strength`, and "supervised" alignment
that down-weights edges within a chosen factor.
* [Jupyter notebook](doc/conos-alignment-strength.ipynb)

## Installation

To install the stable version from [CRAN](https://cran.r-project.org/package=conos), use:

```r
install.packages('conos')
```

To install the latest version of `conos`, use:

```r
install.packages('devtools')
devtools::install_github('kharchenkolab/conos')
```

Conos pre-processes each sample with [pagoda2](https://github.com/kharchenkolab/pagoda2) (or Seurat). The
example data used in the tutorials is the `conosPanel` package:

```r
install.packages('conosPanel', repos = 'https://kharchenkolab.github.io/drat/', type = 'source')
```

Some optional features need extra packages: `Seurat` (Seurat samples), `SeuratDisk` (`.h5seurat` files),
`hdf5r` (`.h5ad` / loom files), and `lstar` (zarr stores / collection round-trips). Each is required only
when you use the corresponding data source.

#### System dependencies

The dependencies are inherited from [pagoda2](https://github.com/kharchenkolab/pagoda2). Note that this package also has the dependency [igraph](https://igraph.org/r/), which requires various libraries to install correctly. Please see the installation instructions at that page for more details, along with the github README [here](https://github.com/igraph/rigraph).

##### Ubuntu dependencies

To install system dependencies using `apt-get`, use the following:
```sh
sudo apt-get update
sudo apt-get -y install libcurl4-openssl-dev libssl-dev libxml2-dev libgmp-dev libglpk-dev
```

##### Red Hat-based distributions dependencies

For Red Hat distributions using `yum`, use the following command:

```sh
sudo yum update
sudo yum install openssl-devel libcurl-devel libxml2-devel gmp-devel glpk-devel
```

##### Mac OS

Using the Mac OS package manager [Homebrew](https://brew.sh/), try the following command:

```sh
brew update
brew install openssl curl-openssl libxml2 glpk gmp
```
(You may need to run `brew uninstall curl` in order for `brew install curl-openssl` to be successful.)

If you hit issues installing `conos` on Mac OS, see the wiki page for further instructions:
[Installing conos for Mac OS](https://github.com/kharchenkolab/conos/wiki/Installing-conos-for-Mac-OS)


### Running conos via Docker

If your system configuration is making it difficult to install `conos` natively, an alternative way to get `conos` running is through a docker container.

**Note:** On Mac OS X, Docker Machine has Memory and CPU limits. To control it, please check instructions either for [CLI](https://stackoverflow.com/questions/32834082/how-to-increase-docker-machine-memory-mac/32834453#32834453) or for [Docker Desktop](https://docs.docker.com/docker-for-mac/#advanced).

#### Ready-to-run Docker image

The docker distribution has the latest version and also includes the [pagoda2 package](https://github.com/kharchenkolab/pagoda2). To start a docker container, first [install docker](https://docs.docker.com/install/) on your platform and then start the `pagoda2` container with the following command in the shell:

```
docker run -p 8787:8787 -e PASSWORD=pass pkharchenkolab/conos:latest
```

The first time you run this command, it will download several large images so make sure that you have fast internet access setup. You can then point your browser to http://localhost:8787/ to get an Rstudio environment with `pagoda2` and `conos` installed (please log in using credentials username=`rstudio`, password=`pass`). Explore the [docker --mount option](https://docs.docker.com/storage/volumes/) to allow access of the docker image to your local files.

**Note:** If you already downloaded the docker image and want to update it, please pull the latest image with: 
```
docker pull pkharchenkolab/conos:latest
```

#### Building Docker image from the Dockerfile

If you want to build image by your own, download the [Dockerfile](https://github.com/kharchenkolab/conos/blob/main/docker/Dockerfile) (available in this repo under `/docker`) and run to following command to build it:
```
docker build -t conos .
```
This will create a "conos" docker image on your system (please be patient, as the build could take approximately 30-50 minutes to finish).
You can then run it using the following command:
```
docker run -d -p 8787:8787 -e PASSWORD=pass --name conos -it conos
```


## References

If you find this software useful for your research, please cite the corresponding [paper](https://doi.org/10.1038/s41592-019-0466-z):

```
Barkas N., Petukhov V., Nikolaeva D., Lozinsky Y., Demharter S., Khodosevich K., & Kharchenko P.V. 
Joint analysis of heterogeneous single-cell RNA-seq dataset collections. 
Nature Methods, (2019). doi:10.1038/s41592-019-0466-z
```

The R package can be cited as:

```
Viktor Petukhov, Nikolas Barkas, Peter Kharchenko, and Evan
Biederstedt (2021). conos: Clustering on Network of Samples. R
package version 1.5.4.
```

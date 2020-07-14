# Conos

- [Conos: Clustering on Network of Samples](#conos-clustering-on-network-of-samples)
- [Tutorials](#tutorials)
  * [Usage: Alignment of Datasets](#usage-alignment-of-datasets)
  * [Integration with ScanPy](#integration-with-scanpy)
  * [Running RNA velocity on a conos object](#running-rna-velocity-on-a-conos-object)
- [Installation](#installation)
  * [Native installation](#native-installation)
    + [System dependencies](#system-dependencies)
      - [Ubuntu dependencies](#ubuntu-dependencies)
      - [Red Hat-based distributions dependencies](#red-hat-based-distributions-dependencies)
      - [OS X](#os-x)
  * [Installing Conos as a Docker container](#installing-conos-as-docker-container)
    + [Ready-to-run Docker image](#ready-to-run-docker-image)
    + [Building Docker image from the Dockerfile](#building-docker-image-from-the-dockerfile)
- [Reference](#reference)
  
## Conos: Clustering On Network Of Samples

* **What is Conos?**
Conos is a package to wire together large collections of single-cell RNA-seq datasets, which allows for both the identification of recurrent cell clusters and the propagation of information between datasets in multi-sample or atlas-scale collections. It focuses on the uniform mapping of homologous cell types across heterogeneous sample collections. For instance, users could investigate a collection of dozens of peripheral blood samples from cancer patients combined with dozens of controls, which perhaps includes samples of a related tissue such as lymph nodes.

* **How does it work?**
![overview](http://pklab.med.harvard.edu/peterk/conos/Figure1_take3.pk.png)
Conos applies one of many error-prone methods to align each pair of samples in a collection, establishing weighted inter-sample cell-to-cell links. The resulting joint graph can then be analyzed to identify subpopulations across different samples. Cells of the same type will tend to map to each other across many such pairwise comparisons, forming cliques that can be recognized as clusters (graph communities). 

  To elaborate in more detail, Conos processing can be divided into three phases:
    * **Phase 1: Filtering and normalization** Each individual dataset in the sample panel is filtered and normalized using standard packages for single-dataset processing: either `pagoda2` or `Seurat`. Specifically, Conos relies on these methods to perform cell filtering, library size normalization, identification of overdispersed genes and, in the case of pagoda2, variance normalization. (Conos is robust to variations in the normalization procedures, but it is recommended that all of the datasets be processed uniformly.)
    * **Phase 2: Identify multiple plausible inter-sample mappings** Conos performs pairwise comparisons of the datasets in the panel to establish an initial error-prone mapping between cells of different datasets. 
    * **Phase 3: Joint graph construction** These inter-sample edges from Phase 2 are then combined with lower-weight intra-sample edges during the joint graph construction. The joint graph is then used for downstream analysis, including community detection and label propagation. For a comprehensive description of the algorithm, please refer to our [publication](https://doi.org/10.1038/s41592-019-0466-z).

* **What does it produce?**
In essence, Conos will take a large, potentially heterogeneous panel of samples and will produce clustering grouping similar cell subpopulations together in a way that will be robust to inter-sample variation:  
![example](http://pklab.med.harvard.edu/peterk/conos/bm_uniform_labels_trim.png)

* **What are the advantages over existing alignment methods?** 
Conos is robust to heterogeneity of samples within a collection, as well as noise. The ability to resolve finer subpopulation structure improves as the size of the panel increases.

* **What do I need to run it?**
Conos is an R package. Currently, it supports pre-processing (filtering, normalization, etc.) of the individual datasets using [pagoda2](https://github.com/hms-dbmi/pagoda2) or [Seurat](https://satijalab.org/seurat/).


## Tutorials

To see the class documentation, run `?Conos`.

### Usage: Alignment of Datasets

Please see the [Conos tutorial walkthrough](vignettes/walkthrough.md) for a detailed example of how to use Conos. The overall runtime of the tutorial should be approximately 5 minutes.

Additional tutorials for Conos include: 
* [Adjustment of Alignment Strength with Conos](vignettes/adjust_alignment_strength.md)
* [Integrating RNA-seq and ATAC-seq](http://pklab.med.harvard.edu/peterk/conos/atac_rna/example.html).

Given a list of individual processed samples (`pl`), Conos processing can be as simple as this:
```r
# construct Conos object, where pl is a list of pagoda2 objects 
con <- Conos$new(pl)

# build graph
con$buildGraph()

# find communities
con$findCommunities()

# plot joint graph
con$plotGraph()

# plot panel with joint clustering results
con$plotPanel()
```

### Integration with ScanPy

For integration with ScanPy, you need to save Conos files on disk from your R session, and 
then upload these files from Python. See the following tutorials:
- [Save Conos for ScanPy](vignettes/scanpy_integration.md)
- [Load ScanPy from Conos](vignettes/scanpy_integration.ipynb)

### Running RNA velocity on a Conos object

First of all, in order to obtain an RNA velocity plot from a Conos object you have to use the [dropEst](https://github.com/hms-dbmi/dropEst) pipeline to align and annotate your single-cell RNA-seq measurements. You can see [this tutorial](http://pklab.med.harvard.edu/velocyto/notebooks/R/SCG71.nb.html) and [this shell script](http://pklab.med.harvard.edu/velocyto/mouseBM/preprocess.sh) to see how it can be done. In this example we specifically assume that when running dropEst you have used the **-V** option to get estimates of unspliced/spliced counts from the dropEst directly. Secondly, you need the [velocyto.R](http://velocyto.org/) package for the actual velocity estimation and visualisation.

After running dropEst you should have 2 files for each of the samples: 
- `sample.rds` (matrix of counts)
- `sample.matrices.rds` (3 matrices of exons, introns and spanning reads)

The `.matrices.rds` files are the velocity files. Load them into R in a list (same order as you give to Conos). Load, preprocess and integrate with Conos the count matrices (`.rds`) as you normally would. Before running the velocity, you must at least create an embedding and run the leiden clustering. Finally, you can estimate the velocity as follows:  
```r
### Assuming con is your Conos object and cms.list is the list of your velocity files ###

library(velocyto.R)

# Preprocess the velocity files to match the Conos object
vi <- velocityInfoConos(cms.list = cms.list, con = con, 
                        n.odgenes = 2e3, verbose = TRUE)

# Estimate RNA velocity
vel.info <- vi %$%
  gene.relative.velocity.estimates(emat, nmat, cell.dist = cell.dist, 
                                   deltaT = 1, kCells = 25, fit.quantile = 0.05, n.cores = 4)

# Visualise the velocity on your Conos embedding 
# Takes a very long time! 
# Assign to a variable to speed up subsequent recalculations
cc.velo <- show.velocity.on.embedding.cor(vi$emb, vel.info, n = 200, scale = 'sqrt', 
                                          cell.colors = ac(vi$cell.colors, alpha = 0.5), 
                                          cex = 0.8, grid.n = 50, cell.border.alpha = 0,
                                          arrow.scale = 3, arrow.lwd = 0.6, n.cores = 4, 
                                          xlab = "UMAP1", ylab = "UMAP2")

# Use cc=cc.velo$cc when running again (skips the most time consuming delta projections step)
show.velocity.on.embedding.cor(vi$emb, vel.info, cc = cc.velo$cc, n = 200, scale = 'sqrt', 
                               cell.colors = ac(vi$cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 15, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 2,
                               do.par = F, cell.border.alpha = 0.1, n.cores = 4,
                               xlab = "UMAP1", ylab = "UMAP2")

```


## Installation

Native installations have been tested in Linux and Mac OS. Normally, installations should take under 10 minutes.

### Native installation

Please make sure that the `devtools` package is installed (use `install.packages("devtools")` if installation is needed).
Then install [pagoda2](https://github.com/hms-dbmi/pagoda2) (or Seurat), then install `conos`:
```r
devtools::install_github("hms-dbmi/conos")
```

If you have problems with `sccore` package, run `devtools::install_github("hms-dbmi/sccore")` before installing `conos`.

#### System dependencies

The dependencies are inherited from [pagoda2](https://github.com/hms-dbmi/pagoda2):

##### Ubuntu dependencies

To install system dependencies using `apt-get`, use the following:
```sh
sudo apt-get update
sudo apt-get -y install libcurl4-openssl-dev libssl-dev
```

##### Red Hat-based distributions dependencies

For Red Hat distributions using `yum`, use the following command:

```sh
yum install openssl-devel libcurl-devel
```

##### OS X

Using the Mac OS package manager [Homebrew](https://brew.sh/), try the following command:

```sh
brew install openssl curl-openssl
```
(You may need to run `brew uninstall curl` in order for `brew install curl-openssl` to be successful.)

**Note:** It is possible to install `pagoda2` and `conos` on OS X, however some users have reported issues with the OpenMP configuration. For instructions, see the [pagoda2](https://github.com/hms-dbmi/pagoda2#mac-dependencies) README.

### Installing Conos as Docker container

If your system configuration is making it difficult to install `conos` natively, an alternative way to get `conos` running is through a docker container.

**Note:** on OS X, Docker Machine has Memory and CPU limits. To control it, please check instructions either for [CLI](https://stackoverflow.com/questions/32834082/how-to-increase-docker-machine-memory-mac/32834453#32834453) or for [Docker Desktop](https://docs.docker.com/docker-for-mac/#advanced).

#### Ready-to-run Docker image

The docker distribution has the latest version and also includes the [Pagoda2 package](https://github.com/hms-dbmi/pagoda2). To start a docker container, first [install docker](https://docs.docker.com/install/) on your platform and then start the `pagoda2` container with the following command in the shell:

```
docker run -p 8787:8787 -e PASSWORD=pass docker.io/vpetukhov/conos:latest
```

The first time you run this command, it will download several large images so make sure that you have fast internet access setup. You can then point your browser to http://localhost:8787/ to get an Rstudio environment with `pagoda2` and `conos` installed (log in using credentials *rstudio* / *pass*). Explore the docker [--mount option]([https://docs.docker.com/storage/volumes/) to allow access of the docker image to your local files.

**Note:** If you already downloaded the docker image and want to update it, please pull the latest image with: 
```
docker pull vpetukhov/conos:latest
```

#### Building Docker image from the Dockerfile

If you want to build image by your own, download the [Dockerfile](https://github.com/hms-dbmi/conos/blob/master/dockers/Dockerfile) (available in this repo under `/dockers`) and run to following command to build it:
```
docker build -t conos .
```
This will create a "conos" docker image on your system (please be patient, as the build takes approximately 30-50 minutes).
You can then run it using the following command:
```
docker run -d -p 8787:8787 -e PASSWORD=pass --name conos -it conos
```


## Reference

If you find this software useful for your research, please cite the corresponding paper:

Barkas N., Petukhov V., Nikolaeva D., Lozinsky Y., Demharter S., Khodosevich K. & Kharchenko P.V. Joint analysis of heterogeneous single-cell RNA-seq dataset collections. Nat. Methods, (2019). [doi:10.1038/s41592-019-0466-z](https://doi.org/10.1038/s41592-019-0466-z)

# conos

- [Clustering on Network of Samples](#clustering-on-network-of-samples)
- [Installation](#installation)
  * [Native installation](#native-installation)
    + [System dependencies](#system-dependencies)
      - [Ubuntu Dependencies](#ubuntu-dependencies)
      - [Red-Hat-based distributions Dependencies](#red-hat-based-distributions-dependencies)
      - [OS X](#os-x)
- [Installing Conos as Docker Container](#installing-conos-as-docker-container)
  * [Ready-to-run docker image](#ready-to-run-docker-image)
  * [Building docker image on the fly](#building-docker-image-on-the-fly)
- [Usage example](#usage-example)
  * [Alignment of datasets](#alignment-of-datasets)
  * [Integration with ScanPy](#integration-with-scanpy)
- [Reference](#reference)
  
## Clustering on Network of Samples

* What is Conos? 
It's a package to wire together large collections of single-cell RNA-seq datasets. It focuses on uniform mapping of homologous cell types across heterogeneous sample collections. For instance, a collection of dozens of peripheral blood samples from cancer patients, combined with dozens of controls. And perhaps also including samples of a related tissue, such as lymph nodes.

* How does it work? 
![overview](http://pklab.med.harvard.edu/peterk/conos/Figure1_take3.pk.png)
Conos applies one of many error-prone methods to align each pair of samples in a collection, establishing weighted inter-sample cell-to-cell links, creating a global joint graph. Cells of the same type will tend to map to each other across many such pair-wise comparisons, forming cliques, that can be recognized as clusters (graph communities). 

* What does it produce?
In essense, Conos will take a large, potentially heterogeneous panel of samples and will produce clustering grouping similar cell subpopulations together in a way that will be robust to inter-sample variation:  
![example](http://pklab.med.harvard.edu/peterk/conos/bm_uniform_labels_trim.png)

* What are the advantages over existing alignment methods? 
Conos is robust to heterogeneity of samples within collection, as well as noise. The ability to resolve finer subpopulation structure improves as the size of the panel increases.

* What do I need to run it?
Conos is an R package. Currently, it supports pre-processing (filtering, normalization, etc.) of the individual datasets using [pagoda2](https://github.com/hms-dbmi/pagoda2) or [Seurat](https://satijalab.org/seurat/).

## Installation

Native installations have been tested in Linux. Normal installation should take <10min.

### Native installation

Please make sure devtools package is installed (use `install.packages("devtools")` to install it if needed).
Then install [pagoda2](https://github.com/hms-dbmi/pagoda2) (or Seurat), then install conos:
```r
devtools::install_github("hms-dbmi/conos")
```

#### System dependencies

The dependencies are inherited from [pagoda2](https://github.com/hms-dbmi/pagoda2):

##### Ubuntu Dependencies

Install system dependencies, example here provided for Ubuntu
```sh
sudo apt-get update
sudo apt-get -y install libcurl4-openssl-dev libssl-dev
```

##### Red-Hat-based distributions Dependencies

```sh
yum install openssl-devel libcurl-devel
```

##### OS X

It is possible to install pagoda2 and Conos on OS X, however some users have reported issues with OpenMP configuration. For instructions see [pagoda2](https://github.com/hms-dbmi/pagoda2#mac-dependencies) readme.

## Installing Conos as Docker Container

If your system configuration is making it difficult to install Conos natively, an alternative way to get Conos running is through a docker container.

**Note:** on OS X, Docker Machine has Memory and CPU limit. To control it, please check instructions either for [CLI](https://stackoverflow.com/questions/32834082/how-to-increase-docker-machine-memory-mac/32834453#32834453) or for [Docker Desktop](https://docs.docker.com/docker-for-mac/#advanced).

### Ready-to-run docker image

The docker distribution has the latest version and also includes the [Pagoda2 package](https://github.com/hms-dbmi/pagoda2). To start a docker container, first [install docker](https://docs.docker.com/install/) on your platform and then start the pagoda container with the following command in the shell:

```
docker run -p 8787:8787 -e PASSWORD=pass docker.io/vpetukhov/conos:latest
```

The first time you run the command it will download several images so make sure that you have fast internet access setup. You can then point your browser to http://localhost:8787/ to get an Rstudio environment with pagoda2 and conos installed (log in using credentials *rstudio* / *pass*). Explore the docker [--mount option]([https://docs.docker.com/storage/volumes/) to allow access of the docker image to your local files.

**Note:** if you already downloaded the docker image and want to update it, please run 
```
docker pull vpetukhov/conos:latest
```

### Building docker image on the fly

If you want to build image by your own, download the [Dockerfile](https://github.com/hms-dbmi/conos/blob/master/dockers/Dockerfile) (available in this repo under `/dockers`) and run to following command to build it:
```
docker build -t conos .
```
This will create a "conos" docker image on your system (be patient, as the build takes ~30-50min or so).
You can then run it using the following command:
```
docker run -d -p 8787:8787 -e PASSWORD=pass --name conos -it conos
```

## Usage example

### Alignment of datasets

Please see [Conos tutorial](vignettes/walkthrough.md) for detailed usage. The overall runtime of the tutorial should be ~5 minutes.

Additional examples: [forcing better alignment](vignettes/adjust_alignment_strength.md), [integrating RNA-seq and ATAC-seq](http://pklab.med.harvard.edu/peterk/conos/atac_rna/example.html).

Given a list of individual processed samples (`pl`), Conos processing can be as simple as this:
```r
# construct conos object, where pl is a list of pagoda2 objects 
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

For integration with ScanPy you need to save Conos files on disk from R session, 
than upload these files from Python. See the following tutorials:
- [Save Conos for ScanPy](vignettes/scanpy_integration.Rmd)
- [Load ScanPy from Conos](vignettes/scanpy_integration.ipynb)

### Running RNA velocity on a conos object

First of all, in order to obtain an RNA velocity plot from a conos object you have to use the [dropEst](https://github.com/hms-dbmi/dropEst) pipeline to align and annotate your single-cell RNA-seq measurments. You can see [this tutorial](http://pklab.med.harvard.edu/velocyto/notebooks/R/SCG71.nb.html) and [this shell script](http://pklab.med.harvard.edu/velocyto/mouseBM/preprocess.sh) to see how it can be done. In this example we specifically assume that when running dropEst you have used the **-V** option to get estimates of unspliced/spliced counts from the dropEst directly. Secondly, you need the [velocyto.R](http://velocyto.org/) package for the actual velocity estimation and visualisation.

After running dropEst you should have 2 files for each of the samples: 
- `sample.rds` (matrix of counts)
- `sample.matrices.rds` (3 matrices of exons, introns and spanning reads)

The `.matrices.rds` files are the velocity files. Load them into R in a list (same order as you give to conos). Load, preprocess and integrate with conos the count matrices (`.rds`) as you normally would. Before running the velocity you have to had at least created an embedding and run the leiden clustering. Finally, you can esitmate the velocity:  
```r
### Assuming con is your conos object and velo_list is the list of your velocity files ###

library(velocyto.R)

# Preprocess the velocity files to match the conos obejct
vi <- VelocityInfoConos(cms_list = velo_list, 
                        con = con, n.odgenes = 2e3, verbose = TRUE)

# Estimate RNA velocity
vel_info <- vi %$%
  gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells=25, cell.dist=cell.dist, 
                                   fit.quantile=0.05, n.cores=4)

# Visualise the velocity on your conos embedding 
# Takes a very long time! 
# Assign to a variable to speed up subsequent recalculations
cc.velo <- show.velocity.on.embedding.cor(vi$emb, vel_info, n=200, scale='sqrt', 
                                          cell.colors=ac(vi$cell.colors, alpha=0.5), 
                                          cex=0.8, grid.n = 50, cell.border.alpha = 0,
                                          arrow.scale = 3, arrow.lwd = 0.6, n.cores=4, 
                                          xlab = "UMAP1", ylab = "UMAP2")

# Use cc=cc.velo$cc when running again (skips the most time consuming delta projections step)
show.velocity.on.embedding.cor(vi$emb, vel_info, n=200, scale='sqrt', 
                               cell.colors=ac(vi$cell.colors, alpha=0.5), 
                               cex=0.8, arrow.scale=15, show.grid.flow=TRUE, 
                               min.grid.cell.mass=0.5, grid.n=40, arrow.lwd=2,
                               do.par=F, cell.border.alpha = 0.1, n.cores=4,
                               xlab = "UMAP1", ylab = "UMAP2", cc=cc.velo$cc)

```

## Reference

If you find this pipeline useful for your research, please consider citing the paper:

Barkas N., Petukhov V., Nikolaeva D., Lozinsky Y., Demharter S., Khodosevich K. & Kharchenko P.V. Joint analysis of heterogeneous single-cell RNA-seq dataset collections. Nat. Methods, (2019). [doi:10.1038/s41592-019-0466-z](https://doi.org/10.1038/s41592-019-0466-z)

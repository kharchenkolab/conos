## Upcoming

## [1.5.1] - 2024-22-January
- Remove C++11 flag
- Fixed various roxygen2 typos and variable names for docs consistency
- S3 method exports, buildWijMatrix.*

## [1.5.0] - 2022-17-October
- Fixed compilation issue on certain OS configurations


## [1.4.9] - 2022-29-September
- Fix multithreading with largeVis, GitHub issue #129

## [1.4.8] - 2022-25-September
- Add back velocityInfoConos()
- Remove 'Matrix.utils' as dependency


## [1.4.7] - 2022-04-September
- Fixes for Matrix, updated versions 1.4.2, 1.5.0
- Use `inherits()` for conditionals with class()


## [1.4.6] - 2022-30-March

- Fixed bug with `snn` argument in `buildGraph()` due to the parameter `snn.k.self`
- Fixed bug in `getDifferentialGenes()` based on warnings from sccore::plapply()

## [1.4.5] - 2022-20-January

### Changed
- Removed unused R packages from "Suggests" (now in sccore), i.e. 'ggrastr', 'pROC', 'pbapply'
- Modified adjustedRand.c from clues to C++ function, revised R code for internal function adjustedRand() in conclass.R

## [1.4.4] - 2021-08-November

### Changed
- Added `fail.on.error=TRUE` in some plapplys
- Re-added `getGeneExpression()` methods for Seurat (lost in merge 552408f)
- Switched to CircleCI


## [1.4.3] - 2021-02-August

### Changed
- Fix the function `parseCellGroups()`, check if clustering exists


## [1.4.2] - 2021-28-June

### Changed

- Add scaling to `scaledMatricesSeurat()`, `scaledMatricesSeuratV3()`
- Change `sccore::plapply()` in `updatePairs()`

## [1.4.1] - 2021-14-May

### Added

- support Seurat objects in `getOdGenesUniformly` and `con$correctGenes`

### Removed

- functions `collapseCellsByType` and `colSumByFactor` are moved to [sccore](https://github.com/kharchenkolab/sccore/)
- removed strong dependency on drat repositories; only used now for the vignettes in the README

## [1.4.0] - 2021-23-Feb

### Changed

- extensive revisions for CRAN upload, including roxygen2 documentation
- replaced relevant C++ code and Rcpp functions with 
N2R and leidenAlg
- vignettes edits, detailing `p2app4conos()` for rendering Conos to pagoda2 application
- updated Dockerfile
- extensively revised vignettes and moved them, due to on CRAN build + check duration limits
- README revisions for clarity

### Added 

- added `getGeneExpression()` for Seurat v2 and v3 (January 2021)
- add parameter `raster.dpi` in `con$plotEmbedding()` to replace `raster.height` and `raster.width`, given these parameters are defunct with rewrite of `ggrastr` (v0.2.0)[https://github.com/VPetukhov/ggrastr/releases/tag/v0.2.0]
- Rjnmf added as Rcpp function
- auxilliary package conosPanel used


### Removed


## [1.3.1] - 2020-24-09

### Changed

- allow multiple embeddings in conclass (July 2020)
- Improved `plotDEheatmap` function
- Fixed bug with `balancing.factor.per.sample` in `buildGraph`
- Fixed some installation problems
- Improve R6 documentation
- Changed `std::cout` to `Rcpp::Rcout` (July 2020)
- Revised README, vignettes (July 2020)


### Added

- multiple embeddings in Conos object (July 2020)
- Write to HDF5 for `saveConosForScanPy()` (July 2020)
- added `ht_opt$message = FALSE` for ComplexHeatmap (July 2020)
- Added checks for `getPerCellTypeDE()` for errors, removing NAs (July 2020)
- Added `ht_opt$message = FALSE` for ComplexHeatmap (July 2020)
- LICENSE (July 2020)

### Removed

- Removed `getCorrectionVector()` and `getPerCellTypeDECorrected` (2 July 2020)
- Removed all neighborhood averaging via `neighborhood.average` (4 July 2020)
- Removed `raster.height` and `raster.width` from `con$plotEmbedding()`, given these parameters are defunct with rewrite of `ggrastr` (v0.2.0)[https://github.com/VPetukhov/ggrastr/releases/tag/v0.2.0]


## [1.3.0] - 2020-19-3

### Changed

- Moved some code to the new package `sccore`
- Fixed inconsistent use of parameters for different spaces in `buildGraph`
- Various small fixes
- Fixed the number of components calculated for the simple PCA rotation
- Conos is R6 class now (instead of refClass)
	
### Added

- Functionality for PAGA graph collapsing
- Parameters `k.same.factor` and `balancing.factor.per.sample` to `buildGraph`. 
  It can be used to improve alignment between different conditions: with `same.factor.downweight`
  it gives the system similar to `k.self` and `k.self.weight`
- plotDEheatmap() function for viewing marker genes
- Function `convertToPagoda2` to create Pagoda 2 from Conos. Helpful for PagodaWebApp.
	
## [1.2.1] - 2019-12-3

### Changed

- Fixed `getDifferentialGenes`
- Fixed testing clustering stability

## [1.2.0] - 2019-11-27

### Changed

- Added mean M value column to the diff. expression output
- Optimized plotting with coloring by genes
- `getDifferentialGenes` uses first clustering by default
- Fixed bug with `collapseCellsByType`. **Note:** probably will affect DE results.
- Added re-normalization of edge weights to fix problem with negative edge weights during label propagation
- Now in plotting 'groups' aren't ignoted if 'gene' is provided: it's used to subset cells for plotting.
- UMAP now set `n_sgd_threads` from `uwot` to `n.cores` by default. It gives much better parallelization, but kills reproducibility.
  Use `n.sgd.cores=1` to get reproducible embeddings.
- Account for `target.dims` in UMAP embedding
- Fixed estimation of `cor.based` with `alingnment.strength == 0`. It removes edges with negative correlation and reduce down-weight of inter-sample edges, which can change results of the alignment.
- Changed default value of `fixed.initial.labels` in `propagateLabels` from `FALSE` to `TRUE`. Presumably, `FALSE` should never be used.
- New output format for label propagation (list of "labels", "uncertainty" and "label.distribution")
- Numerous small bug fixes and small validations for correct arguments
- ScanPy integration tutorials to refelect the changes in `saveConosForScanPy`

### Added

- Metrics to masure specifisity of cell type markers to DE info in `getDifferentialGenes` (parameters `append.specifisity.metrics` and `append.auc`)
- Implementation of label propagateion based on matrix equations (*occured to be too slow*)
- Function `findSubcommunities` to increase resolution for specific clusters
- Parameter `subgroups` to `embeddingPlot`. It allows to plot only cells, belonging to the specified subgroups
- Parameter `keep.limits` to `embeddingPlot`
- Added metrics to masure specifisity of cell type markers to DE info in `getDifferentialGenes` (parameters `append.specifisity.metrics` and `append.auc`)
- `velocityInfoConos` function for RNA velocity analysis on samples integrated with conos (together with supplementary functions `prepareVelocity` and `pcaFromConos`)
- "Running RNA velocity on a conos object" section in README.md (explains usage of the `velocityInfoConos` function)
- Function `getJointCountMatrix` to conos obect
- New possibilities to customize the output of `saveConosForScanPy` 
- Function `parseCellGroups` to parse properly cell groupings depending on user settings
- Function `estimteWeightEntropyPerCell` to visualize alignment quality per cell

## [1.1.2] - 2019-07-16

### Changed

- Changed a description line in the getClusterPrivacy() doc to fix installation under some R3.6 versions (issue 32)

## [1.1.1] - 2019-07-15

### Changed

- Fixed docker build: use of BiocManager, reference to master instead of dev
- Updated src/Makevars to remove the CXX directive, which trips up older versions of R (3.2.x)

## [1.1.0] - 2019-07-02

### Added

- Support for CCA space

### Changed

- `buildGraph` now use PCA space as the default
- fixed common variance rescaling to use geometric mean of the target

## [1.0.3] - 2019-07-02

### Added

- Support for Seurat v3 objects

## [1.0.2] - 2019-06-26

### Added

- Functions to export Conos object to ScanPy

## [1.0.1] - 2019-05-1

### Fixed

- Default value for `cluster.sep.chr` in DE functions is changed from '+' to '<!!>', 
  as it shouldn't be normally present in cluster names
- Removed Boost dependency
- Fixed version of Seurat and fpc packages in Docker

## [1.0.0] - 2019-04-17

### Fixed

- Fixed ggplot2 namespace for function calls in `plotClusterStability`
- Renamed `stable.tree.clusters`, `get.cluster.graph` and `scan.k.modularity`
- Removed exports of largeVis internals, `get.cluster.graph` and `get_nearest_neighbors`
- `embedding` is stored with samples by rows now (i.e. not transposed anymore)
- Using scale.data instead of data in Seurat if provided

### Deprecated

- multitrap.community and multimulti.community functions

## [0.1.0] - 2019-04-17

### Added

- Pre-release version

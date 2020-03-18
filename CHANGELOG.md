## Upcoming

### Changed

- Moved some code to the new sccore package
- Fixed inconsistent use of parameters for different spaces in `buildGraph`
- Various small fixes
- Fixed the number of components calculated for the simple PCA rotation
	
### Added

- Functionality for PAGA graph collapsing
- Parameters `k.same.factor` and `balancing.factor.per.sample` to `buildGraph`. 
  It can be used to improve alignment between different conditions: with `same.factor.downweight`
  it gives the system similar to `k.self` and `k.self.weight`
- plotDEheatmap() function for viewing marker genes
	
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

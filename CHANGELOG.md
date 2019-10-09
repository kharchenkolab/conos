## [1.1.2] - 2017-07-16
 - Changed a description line in the getClusterPrivacy() doc to fix installation under some R3.6 versions (issue 32)

## [1.1.1] - 2019-07-15

## Changed

- Fixed docker build: use of BiocManager, reference to master instead of dev
- Updated src/Makevars to remove the CXX directive, which trips up older versions of R (3.2.x)

## [1.1.0] - 2019-07-02

## Added

- Support for CCA space

## Changed

- `buildGraph` now use PCA space as the default
- fixed common variance rescaling to use geometric mean of the target

## [1.0.3] - 2019-07-02

## Added

- Support for Seurat v3 objects

## [1.0.2] - 2019-06-26

## Added

- Functions to export Conos object to ScanPy

## [1.0.1] - 2019-05-1

## Fixed

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

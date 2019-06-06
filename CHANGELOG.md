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

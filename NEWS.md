# conos 2.0.0

Major release, coordinated with pagoda2 2.0 and lstar.

## New features

* `runGraph()`, `runClustering()`, `runEmbedding()` and `runMarkers()` are the **preferred verbs**
  (aligned with the pagoda2.1 `runX(method=)` vocabulary). The previous names `buildGraph()`,
  `findCommunities()`, `embedGraph()` and `getDifferentialGenes()` are retained as **deprecated-but-working
  aliases** and now emit a deprecation warning.
* `runClustering()` (and `findCommunities()`) now accepts a **string** method name
  (`"leiden"`, `"walktrap"`, `"louvain"`/`"multilevel"`, `"infomap"`, `"fastgreedy"`, `"labelprop"`,
  `"leadingeigen"`) in addition to a function.
* `planIntegration()` — a preliminary step that polls each sample's molecular modalities (pagoda2.1
  facets / Seurat assays) and reports, per modality, the shared-feature count and per-pair overlap with a
  usable / marginal / not-integrable verdict (surfacing e.g. independently-called scATAC peaks that cannot
  reconcile without consistent pre-processing). Records the resolved default modality on the object.
* `getModalities()` / `getDefaultModality()` accessors (pagoda2.1 facets, Seurat assays, or a single
  `"RNA"` for legacy/single-modality objects).
* `plotMarkerDotPlot()` — top per-cluster markers rendered through `sccore::dotPlot` (the conos counterpart
  of pagoda2.1's `plotMarkerDotPlot`).
* `scanResolution()` — scan a community method's `resolution` over a range, reporting cluster count and
  modularity.
* `buildGraph(pairs.storage = c("keep", "drop", "disk"))` — optionally drop or offload the O(n^2) per-pair
  alignment rotations to free memory on large panels (`"disk"` is restored on the next build).

## Changes

* The default Leiden `n.iterations` is now **5** (was leidenAlg's default of 2) for better convergence on
  large joint graphs.
* The default alignment `space` remains `"PCA"` (reciprocal PCA); `"CPCA"` stays available but non-default,
  and `"CCA"` is opt-in.

## Bug fixes

* `snn = TRUE` graphs build again under R >= 4.2 (an `&&` applied to a length-2 `snn.quantile` previously
  raised "length = 2 in coercion to logical(1)").
* `scanKModularity(plot = TRUE)` now draws its plot (it was built and silently discarded) and its x/y axis
  labels are no longer swapped.
* `plotClusterStability(what = "dend")` returns the dendrogram (invisibly) instead of `NULL`.
* `propagateLabels(method = "solver", solver = "Matrix")` no longer spuriously warns about `rmumps`
  (the warning is now raised only when the `"mumps"` solver is requested but unavailable).
* `propagateLabels(method = "solver")` no longer errors when the labelled (reference) cells miss one or
  more cluster levels (the label-indicator matrix is now sized to the number of levels).

## Performance and memory

* `getJointCountMatrix()` no longer allocates a dense zero-block when sample gene sets differ (the merge
  stays sparse).
* SNN edge weights are computed only at the mutual-nearest-neighbour entries, without densifying the
  n1 x n2 neighbour product.
* The `buildGraph` path is verified to stream sample data (block reads), so lstar-backed pagoda2.1 samples
  keep peak memory flat.

## Documentation / CRAN

* `@return` documentation added to all sample accessor generics and to `saveDEasCSV()`; the
  `getSampleNamePerCell()` example is now valid.
* Robust pagoda2.1 sample accessors (`getPca`/`getEmbedding`/`getClustering`), and symmetric `Conos`
  methods for `getPca()` / `getCountMatrix()` / `getEmbedding()`.

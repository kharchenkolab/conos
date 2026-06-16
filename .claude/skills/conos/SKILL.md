---
name: conos
description: >-
  Use when integrating MULTIPLE single-cell samples into one joint analysis with the conos R package:
  wiring a panel of datasets into a single joint graph, then jointly clustering, embedding, finding
  markers, and transferring (propagating) cell-type labels across samples. Covers building a panel from
  pagoda2 and/or Seurat objects (v3/v4/v5) or from files written by other pipelines (AnnData/.h5ad,
  Seurat/.h5seurat, loom, lstar/zarr — read via pagoda2's from*() constructors), MIXED panels, the
  recommended verbs runGraph/runClustering/runEmbedding/runMarkers (old buildGraph/findCommunities/
  embedGraph/getDifferentialGenes are deprecated aliases), alignment-strength control (alignment.strength,
  supervised balancing.factor.per.cell), plotGraph/plotPanel/plotMarkerDotPlot, propagateLabels, and
  exporting a whole collection (lstar zarr, Seurat v5 split assay, AnnData for scanpy). Conos integrates
  in GRAPH space — there is no batch-corrected expression matrix. Keywords: conos, integration, joint
  graph, alignment, panel, multi-sample, batch, runGraph, runClustering, runEmbedding, runMarkers, Leiden,
  UMAP, propagateLabels, label transfer, plotMarkerDotPlot, alignment.strength, mNN, reciprocal PCA, CPCA,
  CCA, pagoda2, Seurat, scanpy, lstar, collection, atlas.
when_to_use: >-
  Use for two or more single-cell samples that should be analyzed together — recurrent cell types across
  a panel, atlas-scale collections, cross-sample label transfer, or correcting batch/donor/tissue
  separation while preserving biology. Each sample is first pre-processed on its own (pagoda2 or Seurat);
  conos then aligns them. For a SINGLE dataset, use the pagoda2 skill instead.
avoid_when: >-
  Do not use for a single dataset (use pagoda2), for trajectory/RNA-velocity analysis, or for a
  Seurat/scanpy-native integration unless the user explicitly wants conos. Conos does not produce a
  corrected expression matrix; if the user needs integrated/denoised expression values, say so.
keywords: [conos, integration, joint graph, multi-sample, panel, alignment, batch correction, runGraph, runClustering, runEmbedding, runMarkers, Leiden, UMAP, label transfer, propagateLabels, plotMarkerDotPlot, alignment.strength, reciprocal PCA, CPCA, CCA, mNN, pagoda2, Seurat, scanpy, lstar, collection, atlas]
domain: genomics
---

# conos — joint analysis of single-cell sample collections

conos wires together a **panel of single-cell datasets** into one **joint graph**, then analyzes that
graph: shared cell types map to each other across samples (forming cross-sample communities), so you get
one clustering, one embedding, one set of markers, and the ability to **transfer labels** from any sample
to the rest. It is robust to sample heterogeneity, and resolution improves as the panel grows. Repo:
`~/p21/conos` (work on branch `dev`).

## The model in three sentences

- A **panel** is a *named list of per-sample objects* (each a `Pagoda2` or `Seurat` object, already
  normalized + PCA'd on its own); `Conos$new(panel)` collects them.
- conos performs error-prone pairwise alignments between samples, then combines those inter-sample edges
  with intra-sample edges into ONE **joint graph** over the union of all cells; clustering, embedding,
  markers and label transfer all operate on that graph.
- Integration is in **graph space**: there is **no batch-corrected expression matrix**. Per-sample raw
  counts stay as they are; the joint layer is graph + embedding + clustering.

## The recommended API (use these verbs)

`runGraph()`, `runClustering()`, `runEmbedding()`, `runMarkers()` are the preferred verbs (aligned with
pagoda2.1's `runX(method=)` vocabulary). The old names — `buildGraph()`, `findCommunities()`,
`embedGraph()`, `getDifferentialGenes()` — still work but are **deprecated aliases** that emit a warning.

```r
library(conos)
library(pagoda2)

# 1. pre-process each sample on its own (one Pagoda2 object per count matrix; run() adds the PCA)
samples <- lapply(cms, function(cm) Pagoda2$from(cm)$run(steps = c("variance", "pca")))

# 2. the joint analysis is a handful of calls on the Conos object
con <- Conos$new(samples, n.cores = 4)
con$runGraph()          # align every pair of samples into one joint graph
con$runClustering()     # joint Leiden communities
con$runEmbedding()      # joint 2-D embedding (largeVis by default; method = "UMAP" also)

# 3. inspect
con$plotGraph()                         # joint embedding, colour by cluster / sample / gene
con$plotPanel(clustering = "leiden")    # the same, faceted per sample
con$plotMarkerDotPlot()                 # specific per-cluster markers

# 4. transfer annotations from one (or some) labelled sample(s) to the rest
new.labels <- con$propagateLabels(labels = cellannot)$labels
```

Samples may be `Pagoda2` or `Seurat` objects, or read from files — see `reference/data_sources.md`. Build
each Seurat sample the usual way (`NormalizeData`→`FindVariableFeatures`→`ScaleData`→`RunPCA`); conos
aligns on the PCA.

## Main usage patterns

- **Standard integration** — the four-step block above. Defaults: `space = "PCA"` (reciprocal PCA),
  Leiden clustering, largeVis embedding. Full parameters: `reference/workflow.md`.
- **Tune the alignment** — if samples still separate by batch/tissue, raise `alignment.strength` (0→1) or
  use *supervised* alignment (`balancing.factor.per.cell`, `same.factor.downweight`). Caution: over-mixing
  merges real populations. See `reference/alignment.md`.
- **Markers across a panel** — `runMarkers()` / `plotMarkerDotPlot()` work on pagoda2, Seurat, and **mixed**
  panels; markers are computed consistently regardless of sample type. See `reference/markers.md`.
- **Label transfer** — `propagateLabels()` diffuses labels over the joint graph and returns per-cell
  labels + uncertainty. See `reference/label_transfer.md`.
- **Build a panel from any source** — pagoda2/Seurat objects in memory, or `.h5ad`/`.h5seurat`/loom/lstar
  files via `Pagoda2$from*()`, including mixed panels. See `reference/data_sources.md`.
- **Export / interchange a whole collection** — to an lstar zarr store (`write_conos`/`read_conos`), a
  Seurat v5 split-assay object (`write_seurat`), or an AnnData for scanpy (`lstar convert`). See
  `reference/export_interop.md`.

## Key principles (do not violate)

- **Use the recommended verbs.** `runGraph`/`runClustering`/`runEmbedding`/`runMarkers`, not the deprecated
  `buildGraph`/`findCommunities`/`embedGraph`/`getDifferentialGenes`.
- **Process samples uniformly.** conos is robust to normalization differences, but align like with like —
  pre-process every sample the same way; don't mix wildly different pipelines if avoidable.
- **No corrected expression matrix.** Integration is in graph space. Never fabricate a "corrected" or
  "integrated" expression matrix; `getJointCountMatrix()` returns per-sample raw/normalized values merged,
  not a batch-corrected matrix.
- **Judge alignment by the sample-coloured embedding.** Well-aligned samples interleave; reach for
  `alignment.strength` only when residual separation is batch, not biology, and re-check markers after.
- **Stay disk-backed / streaming where the sample backend is.** lstar-backed pagoda2 samples stream;
  don't pull whole matrices into memory to compute things yourself when an accessor or the per-sample
  object can do it.
- **Markers must be consistent across a mixed panel.** pagoda2 samples use pagoda2's own disk-backed DE;
  other samples use the same Wilcoxon core via `sccore::matrixDE()` — so results don't depend on sample
  type.

## Reference files (read the one you need)

- `reference/workflow.md` — the full `Conos$new` → `runGraph` → `runClustering` → `runEmbedding` API:
  every parameter that matters (alignment `space`, `k`/`k.self`, SNN, `pairs.storage`), `scanResolution`,
  `scanKModularity`, `plotGraph`/`plotPanel`, and the deprecated-alias map.
- `reference/data_sources.md` — building a panel from pagoda2/Seurat (v3/v4/v5) objects, from files
  (`Pagoda2$fromAnnData/fromH5Seurat/fromLoom/fromLstar`), mixed panels, `planIntegration()`,
  `getModalities()`/`getDefaultModality()`.
- `reference/markers.md` — `runMarkers()`, `plotMarkerDotPlot()` (balanced selection, `min.auc`),
  class-agnostic + mixed-panel markers (pagoda2 path vs `sccore::matrixDE`), `saveDEasCSV()`.
- `reference/label_transfer.md` — `propagateLabels()` (diffusion vs solver), uncertainty, common patterns.
- `reference/alignment.md` — `alignment.strength`, supervised alignment (`balancing.factor.per.cell`,
  `same.factor.downweight`), and how to choose values without over-mixing.
- `reference/export_interop.md` — exporting a whole collection via lstar (zarr round-trip, Seurat v5,
  AnnData), reading a collection back from a Seurat v5 `.rds`, and the recommended **Python (scanpy)**
  workflow. (No corrected matrix is ever fabricated.)

Tutorials (rendered notebooks) live in `doc/`: `conos-walkthrough`, `conos-advanced`,
`conos-data-sources`, `conos-alignment-strength`.

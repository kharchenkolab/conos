# Data sources: building a panel from anything

A conos panel is just a **named list of per-sample objects**. They need not all be pagoda2 objects built
this session. conos works with pagoda2 objects, Seurat objects (v3/v4/v5) in memory, samples read from
files written by other pipelines, and **mixed** panels.

## pagoda2 objects (the default)

```r
samples <- lapply(cms, function(cm) Pagoda2$from(cm)$run(steps = c("variance", "pca")))
con <- Conos$new(samples)
```

`run(steps = c("variance", "pca"))` is the minimum conos needs (variance normalization + the PCA it aligns
on). A full `p2$run()` also embeds/clusters each sample, which is unnecessary for conos but harmless.

## Seurat objects in memory (v3/v4/v5)

Hand `Seurat` objects straight to `Conos$new()`. Build each the usual way; conos aligns on the PCA:

```r
seurat.samples <- lapply(cms, function(cm) {
  s <- CreateSeuratObject(counts = cm)
  s <- NormalizeData(s, verbose = FALSE)
  s <- FindVariableFeatures(s, nfeatures = 2000, verbose = FALSE)
  s <- ScaleData(s, verbose = FALSE)
  RunPCA(s, npcs = 30, verbose = FALSE)
})
con <- Conos$new(seurat.samples)
con$runGraph(space = "PCA", ncomps = 30, n.odgenes = 2000)
```

Seurat v5 / `SeuratObject >= 5` works: the count-matrix accessors use `layer=` on v5 and `slot=` on v3/v4
(via the internal `getSeuratAssayData()`), so `GetAssayData(slot=)` being defunct in SeuratObject 5 is
handled. `getCountMatrix.Seurat` returns the **log-normalized `data` layer** (non-negative), not
`scale.data`.

## Reading samples from files (no Python required)

pagoda2's `from*()` constructors read each format directly via HDF5/zarr. Each returns a `Pagoda2` object;
add the variance + PCA reduction conos needs:

```r
s1 <- Pagoda2$fromAnnData("sample1.h5ad")$run(steps = c("variance", "pca"))     # anndata / scanpy
s2 <- Pagoda2$fromH5Seurat("sample2.h5seurat")$run(steps = c("variance", "pca")) # Seurat (SeuratDisk)
s3 <- Pagoda2$fromLoom("sample3.loom")$run(steps = c("variance", "pca"))         # loom
s4 <- Pagoda2$fromLstar("sample4.zarr")$run(steps = c("variance", "pca"))        # lstar (zarr) store
con <- Conos$new(list(s1, s2, s3, s4))
```

Notes:
- `fromAnnData` / `fromLstar` preserve gene names exactly (best for round-trips).
- `fromH5Seurat` reads any valid `.h5seurat`; if a file lacks stored gene names (e.g. a Seurat v5 object
  written by a `SeuratObject < 5` build of SeuratDisk, which can't serialize `Assay5` feature names),
  pagoda2 **warns** and falls back to positional names rather than silently mislabelling genes.
- Constructor reader options go in `reader.args = list(layer = "counts", sample.name = ...)`.

## Mixed panels

A panel may combine sources — e.g. one `Pagoda2` and one `Seurat` object:

```r
mixed <- list(
  bm_pagoda2 = Pagoda2$from(cm1)$run(steps = c("variance", "pca")),
  cb_seurat  = RunPCA(ScaleData(FindVariableFeatures(NormalizeData(CreateSeuratObject(counts = cm2)))))
)
con <- Conos$new(mixed)
con$runGraph()   # heterogeneous panel: a type-agnostic scaled-matrix builder, not "assume the first class"
```

`runGraph` detects heterogeneous panels and routes through `scaledMatricesMixed` (per-sample
`getCountMatrix` + per-gene unit-variance scaling), instead of assuming all samples share the first
sample's class. Markers on a mixed panel are computed consistently (see `markers.md`).

## Inspecting integrability and modalities

```r
con$planIntegration(min.common.features = 5)   # per modality: shared-feature count + per-pair overlap +
                                               # usable/marginal/not-integrable verdict; records default modality
getModalities(sample)        # pagoda2.1 facets, Seurat assays, or "RNA" for single-modality
getDefaultModality(sample)
```

`planIntegration()` surfaces, e.g., independently-called scATAC peaks that can't reconcile without
consistent pre-processing.

## Accessors (work across Pagoda2 / Seurat / Conos)

`getPca`, `getCountMatrix`, `getRawCountMatrix`, `getEmbedding`, `getClustering`, `getOverdispersedGenes`,
`getCellNames`, `getGenes`, `getGeneExpression`, `getDatasetPerCell`, `getJointCountMatrix(raw = FALSE)`.
Prefer these to reaching into object internals — they keep disk-backed samples streaming.

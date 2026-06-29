# Data sources: building a panel from anything

A conos panel is just a **named list of per-sample objects**. They need not all be pagoda2 objects built
this session. conos works with pagoda2 objects, Seurat objects (v3/v4/v5) in memory, samples read from
files written by other pipelines, and **mixed** panels.

## pagoda2 objects (the default)

Use the flavor-routed `preprocess_p2()` helper from **SKILL.md Step 1** (handles pagoda2 1.x CRAN
*and* 2.0 devel — never call `Pagoda2$from`/`$run` unguarded, they don't exist on CRAN 1.x):

```r
samples <- lapply(cms, preprocess_p2)   # preprocess_p2 defined in Step 1
con <- Conos$new(samples)
```

It produces the minimum conos needs (variance normalization + the PCA it aligns on). On devel a full
`p2$run()` would also embed/cluster each sample — unnecessary for conos but harmless.

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

**This is flavor-dependent.** pagoda2 **devel (2.0)** has direct `from*()` constructors for every format
(HDF5/zarr); pagoda2 **CRAN (1.x)** reads only **10x** natively (`read10xMatrix()`/`read.10x.matrices()`).
Route through this helper (it reuses `preprocess_p2()` from Step 1):

```r
# Read ONE sample file into a pre-processed Pagoda2 object, across pagoda2 flavors.
read_sample_p2 <- function(path, format = c("10x", "h5ad", "h5seurat", "loom", "lstar"), n.cores = 1) {
  format <- match.arg(format)
  if (pagoda2_is_devel()) {                       # devel: native readers for every format
    reader <- switch(format,
      "10x"      = pagoda2::Pagoda2$from10x,
      "h5ad"     = pagoda2::Pagoda2$fromAnnData,
      "h5seurat" = pagoda2::Pagoda2$fromH5Seurat,
      "loom"     = pagoda2::Pagoda2$fromLoom,
      "lstar"    = pagoda2::Pagoda2$fromLstar)
    reader(path, n.cores = n.cores)$run(steps = c("variance", "pca"), verbose = FALSE)
  } else {                                          # CRAN 1.x: only 10x is native
    if (format != "10x")
      stop("pagoda2 ", utils::packageVersion("pagoda2"), " (CRAN) reads only 10x directly. For '",
           format, "', install pagoda2 devel (>= 2.0), or load the matrix with another package ",
           "(anndata/SeuratDisk/etc.) and pass it to preprocess_p2().")
    cm <- pagoda2::read10xMatrix(path, verbose = FALSE)
    preprocess_p2(cm, n.cores = n.cores)
  }
}

# devel: any of these; CRAN 1.x: the 10x line works, the others stop() with the upgrade hint.
s1 <- read_sample_p2("sample1.h5ad",     "h5ad")       # anndata / scanpy   (devel only)
s2 <- read_sample_p2("sample2.h5seurat", "h5seurat")   # Seurat (SeuratDisk) (devel only)
s3 <- read_sample_p2("sample3.loom",     "loom")       # loom               (devel only)
s4 <- read_sample_p2("sample4.zarr",     "lstar")      # lstar (zarr) store  (devel only)
s5 <- read_sample_p2("sample5_10x_dir",  "10x")        # 10x                 (both flavors)
con <- Conos$new(list(s1, s2, s3, s4, s5))
```

Notes:
- `fromAnnData` / `fromLstar` (devel) preserve gene names exactly (best for round-trips).
- `fromH5Seurat` (devel) reads any valid `.h5seurat`; if a file lacks stored gene names (e.g. a Seurat v5
  object written by a `SeuratObject < 5` build of SeuratDisk, which can't serialize `Assay5` feature
  names), pagoda2 **warns** and falls back to positional names rather than silently mislabelling genes.
- Constructor reader options (devel) go in `reader.args = list(layer = "counts", sample.name = ...)`.
- **CRAN 1.x + non-10x:** convert upstream — e.g. read the `.h5ad` with the `anndata`/`zellkonverter`
  package to a matrix, then `preprocess_p2(cm)` — or just install pagoda2 devel.

## Mixed panels

A panel may combine sources — e.g. one `Pagoda2` and one `Seurat` object:

```r
mixed <- list(
  bm_pagoda2 = preprocess_p2(cm1),    # flavor-routed (Step 1); NOT a bare Pagoda2$from
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

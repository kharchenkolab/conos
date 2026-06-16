# Export & interoperability (via lstar)

A `Conos` object is a **whole collection**: every sample's raw counts *plus* the joint analysis layer —
graph, embedding, clustering — over the union of cells. Because conos integrates in **graph space**, there
is **no batch-corrected expression matrix**. lstar writes the whole collection to native multi-sample
formats *without fabricating* such a matrix: per-sample raw counts are preserved, and the joint graph /
embedding / clusters travel alongside.

lstar provides the conversions (R package `lstar`, plus a `lstar convert` CLI). Conos has no native
exporter — use lstar.

## lstar (zarr) store — native round-trip

```r
library(lstar)
lstar_write(write_conos(con), "panel.lstar.zarr")   # Conos -> L* collection -> zarr on disk
con2 <- read_conos("panel.lstar.zarr")              # reopen as a live Conos:
                                                    # samples + joint graph/clustering/embedding restored
```

`read_conos()` reconstitutes per-sample `Pagoda2` objects (raw counts + PCA) and restores the joint layer,
so you can keep plotting, finding markers and transferring labels without recomputing the integration.
(Re-running `runGraph()` would recompute the per-sample variance model, which is not stored.)

## Seurat v5

```r
so <- write_seurat(write_conos(con))   # -> Seurat v5: a *split* assay (one raw layer per sample) +
                                       #    the joint graph as a Seurat Graph + embedding as a DimReduc +
                                       #    sample / cluster labels in meta.data
saveRDS(so, "panel_seurat_v5.rds")

# read a collection BACK from a Seurat v5 .rds:
con3 <- read_conos(read_seurat(readRDS("panel_seurat_v5.rds")))   # joint graph restored
```

`read_seurat()` reads `SeuratObject::Graphs(so)`, so the joint integration graph survives the round-trip
(not just the per-sample layers). `read_conos()` finds the graph by the field name `graph`, or falls back
to any `cells × cells` relation field, so a generic Seurat-v5 integration object opens as a `Conos` too.

## AnnData (for scanpy)

The stored collection flattens to a single AnnData: `X` = **raw joint counts** (no corrected matrix), the
joint graph in `obsp` (aliased to `connectivities`, with a matching `uns['neighbors']`) so scanpy's graph
tools run on the conos graph directly, the embedding in `obsm`, sample/cluster in `obs`. Convert the store
with the lstar CLI:

```bash
lstar convert panel.lstar.zarr panel.h5ad
```

(The same `lstar convert` also moves individual objects between formats, e.g.
`lstar convert sample.rds sample.h5ad`.)

## Running conos from a Python (scanpy) workflow

conos itself is R-only, but the two directions add up to a clean recommendation for scanpy users: keep
data as `.h5ad` and bridge through one short R step.

```r
# conos_step.R  (run: Rscript conos_step.R)
library(conos); library(pagoda2); library(lstar)
samples <- lapply(c("s1.h5ad", "s2.h5ad", "s3.h5ad"),
                  function(f) Pagoda2$fromAnnData(f)$run(steps = c("variance", "pca")))
con <- Conos$new(samples)
con$runGraph(); con$runClustering(); con$runEmbedding()
lstar_write(write_conos(con), "panel.lstar.zarr")
```

```bash
Rscript conos_step.R
lstar convert panel.lstar.zarr integrated.h5ad
```

```python
import scanpy as sc
adata = sc.read_h5ad("integrated.h5ad")
sc.tl.leiden(adata)   # runs on the conos joint graph in obsp['connectivities']
sc.tl.umap(adata)     # conos's own embedding is already in adata.obsm
```

There is **no pure-Python conos**; the R step *is* the integration. To stay entirely in Python, a
Python-native method (scVI, Harmony, BBKNN) is the alternative — but it won't reproduce conos's
joint-graph alignment.

## Note: lstar lives in `~/p21/lstar`

The conos-related conversions are in lstar's R profiles (`write_conos`/`read_conos`,
`write_seurat`/`read_seurat`) and exercised by `conformance/conos.sh`. lstar is local git only.

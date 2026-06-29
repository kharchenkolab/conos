---
name: conos-scrna-integration
description: Integrate MULTIPLE single-cell samples with the conos R package — pre-process each sample (pagoda2 or Seurat), wire them into one joint graph, then jointly cluster, embed (UMAP), find markers, and propagate (transfer) cell-type labels across samples. Covers panels built from pagoda2 and/or Seurat objects or files (h5ad/h5seurat/loom/lstar), mixed panels, alignment-strength control, and export to lstar/Seurat-v5/AnnData.
when_to_use: Use for two or more single-cell samples that should be analyzed together — recurrent cell types across a panel, atlas-scale collections, cross-sample label transfer, or removing batch/donor/tissue separation while preserving biology. Each sample is pre-processed on its own first; conos aligns them in graph space.
avoid_when: Do not use for a SINGLE dataset (use the pagoda2 single-dataset recipe), for trajectory/RNA-velocity analysis, or for a Seurat/scanpy-native integration unless the user explicitly wants conos. Conos does NOT produce a batch-corrected expression matrix — if the user needs integrated/denoised expression values, say so.
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [R, conos, pagoda2, sccore, leidenAlg, igraph]
keywords: [conos, integration, joint graph, multi-sample, panel, alignment, batch correction, runGraph, runClustering, runEmbedding, runMarkers, Leiden, UMAP, label transfer, propagateLabels, plotMarkerDotPlot, alignment.strength, reciprocal PCA, CPCA, CCA, mNN, pagoda2, Seurat, scanpy, lstar, collection, atlas]
produces: [joint_umap_clusters.png, joint_umap_samples.png, panel_clusters.png, marker_dotplot.png, label_transfer_umap.png, cluster_markers_*.csv, conos_integrated.rds]
domain: genomics
source: "conos 2.0 (GitHub dev) — source-verified R6 Conos methods (R/conclass.R, R/conos.R) + doc/ tutorials."
---

# Multi-sample single-cell integration with conos

Wire a **panel** of single-cell samples into one **joint graph**, then analyze that
graph: cells of the same type map to each other across samples (forming cross-sample
communities), giving one joint clustering, one joint embedding, one set of markers, and
the ability to **transfer cell-type labels** from any annotated sample to the rest.

Each sample is first pre-processed on its own (with pagoda2 or Seurat); conos then aligns
them. Integration happens in **graph space**: per-sample raw counts stay as they are and
there is **no batch-corrected expression matrix** — the joint layer is the graph +
embedding + clustering. For a single dataset, use the pagoda2 single-dataset recipe
instead; for converting/exporting the result between formats, see the lstar recipe.

## Bundled references — load on demand

This SKILL.md is self-contained for the standard integration. Load a reference only for a
variant, a parameter detail, or troubleshooting:

- `references/workflow_and_graph.md` — `Conos$new()`, `runGraph()` (alignment `space`,
  `k`/`k.self`, `ncomps`, `n.odgenes`, SNN, `pairs.storage`), `runClustering()` (string
  vs function methods, `scanResolution`, `scanKModularity`), `runEmbedding()`
  (largeVis vs UMAP), `plotGraph`/`plotPanel`, and the deprecated-alias map.
- `references/data_sources.md` — building a panel from pagoda2 (CRAN 1.x *and* devel 2.0) /
  Seurat (v3/v4/v5) objects, from files (the flavor-aware `read_sample_p2()`; devel reads
  h5ad/h5seurat/loom/lstar, CRAN 1.x reads 10x), **mixed** panels, `planIntegration()`,
  modality accessors.
- `references/alignment_strength.md` — `alignment.strength` and supervised alignment
  (`balancing.factor.per.cell`, `same.factor.downweight`); how to dial mixing without
  over-merging.
- `references/markers.md` — `runMarkers()`, `plotMarkerDotPlot()` (balanced selection,
  `min.auc`), class-agnostic + mixed-panel markers, `saveDEasCSV()`.
- `references/label_transfer.md` — `propagateLabels()` (diffusion vs solver), uncertainty,
  common patterns.
- `references/export_and_interop.md` — export a whole collection via lstar (zarr round-trip,
  Seurat v5, AnnData) and the Python/scanpy recommendation. (No corrected matrix is faked.)

## Install

This recipe requires **conos >= 2.0** (the `run*` verbs + the pagoda2-flavor-agnostic
accessor shim) and **sccore >= 1.1.0** (the native marker-heatmap engine); both are on
GitHub until the next CRAN release. **pagoda2 either flavor works** — CRAN 1.x (`main`) or
devel 2.0 — because conos 2.0 reads counts from both; the recipe keeps whatever pagoda2 you
have installed and only branches its own per-sample calls (see Step 1). The workflow also
needs `leidenAlg` (clustering) and `uwot` (UMAP). `conosPanel` (a `drat` repo) supplies the
example panel used in tutorials.

The version guards below matter: a **stale CRAN conos 1.5.4** must be upgraded (it has no
`runGraph`/`runClustering`), and a plain `requireNamespace` check would silently skip it —
that mismatch (devel-API recipe code on a CRAN-stack install) is the classic "it fails on
both" trap.

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in c("remotes", "ggplot2", "uwot", "leidenAlg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# sccore >= 1.1.0 (native marker-heatmap engine pagoda2/conos call) — upgrade if older
if (!requireNamespace("sccore", quietly = TRUE) ||
    utils::packageVersion("sccore") < "1.1.0") {
  remotes::install_github("kharchenkolab/sccore", ref = "dev", upgrade = "never")
}
# conos >= 2.0 REQUIRED (run* verbs + flavor-agnostic shim) — upgrade a stale CRAN 1.5.4 too
if (!requireNamespace("conos", quietly = TRUE) ||
    utils::packageVersion("conos") < "2.0.0") {
  remotes::install_github("kharchenkolab/conos", ref = "dev", upgrade = "never")
}
# pagoda2: EITHER flavor is fine. Install devel only if pagoda2 is entirely absent; an existing
# CRAN (1.x) install is kept and handled by Step 1's branch. (Need devel for non-10x file readers
# — h5ad/h5seurat/loom/lstar; see references/data_sources.md.)
if (!requireNamespace("pagoda2", quietly = TRUE)) {
  remotes::install_github("kharchenkolab/pagoda2", ref = "devel", upgrade = "never")
}
# example panel (optional; tutorials use it):
if (!requireNamespace("conosPanel", quietly = TRUE)) {
  install.packages("conosPanel", repos = "https://kharchenkolab.github.io/drat/", type = "source")
}

library(conos)
library(pagoda2)
cat("conos", as.character(packageVersion("conos")),
    "| pagoda2", as.character(packageVersion("pagoda2")),
    if (is.function(tryCatch(Pagoda2$from, error = function(e) NULL))) "(devel API)" else "(CRAN 1.x API)", "\n")
```

## Decisions to surface up front

Tell the user these are the integration-defining choices:

1. **Per-sample pre-processing** — every sample is normalized + PCA'd on its own before
   conos. Use the SAME pipeline for all samples (pagoda2 for all, or Seurat for all);
   conos is robust to differences but uniform pre-processing is recommended. Mixed
   pagoda2 + Seurat panels work (see `references/data_sources.md`).
2. **Alignment space** — default `space = "PCA"` (reciprocal PCA: fast, robust). `"CPCA"`
   and `"CCA"` are opt-in. Use the default unless you have a reason.
3. **Alignment strength** — default is deliberately gentle. Only raise `alignment.strength`
   (or use supervised alignment) if the sample-coloured embedding shows batch-like
   separation; over-mixing merges real populations (`references/alignment_strength.md`).
4. **Clustering resolution** — Leiden by default. Check cluster count/sizes; scan
   `resolution` if over/under-split (`references/workflow_and_graph.md`).
5. **No corrected expression** — conos integrates the graph, not the expression matrix.
   Markers and label transfer run on the joint graph; there is no denoised expression
   output to export.
6. **Compute footprint** — `n.cores` parallelizes the (dominant) pairwise alignments; on
   large panels use `pairs.storage = "drop"`/`"disk"` to cap memory.

Show the user these figures as the analysis proceeds:

- `joint_umap_clusters.png`, `joint_umap_samples.png` (the sample-coloured one is how you
  judge alignment)
- `panel_clusters.png`
- `marker_dotplot.png`
- `label_transfer_umap.png` (if transferring labels)

---

## Step 1 — Pre-process each sample into a panel

A panel is a **named list** of per-sample objects, each already normalized with a PCA
reduction. The most common case is one `Pagoda2` object per count matrix.

> **pagoda2 has two API generations — the recipe supports BOTH.** pagoda2 **1.x (CRAN,
> the `main` branch)** and **2.0+ (the `devel` branch)** build a sample differently, and a
> snippet written for one ERRORS on the other (e.g. `Pagoda2$from` doesn't exist on 1.x →
> `"attempt to apply non-function"`). Detect the flavor once and route. conos 2.0 reads
> counts from *either* flavor internally (its accessor shim falls back from
> `getExpressionBlock()` to `$counts`), so **only the per-sample construction/IO needs
> branching — not the conos steps.** Define this helper and reuse it everywhere a sample is
> built from a matrix:

```r
# TRUE on pagoda2 >= 2.0 (devel): the unified Pagoda2$from(...)$run(steps=) API exists.
# FALSE on pagoda2 1.x (CRAN/main): only the classic Pagoda2$new(...) + explicit methods.
pagoda2_is_devel <- function() is.function(tryCatch(pagoda2::Pagoda2$from, error = function(e) NULL))

# Build ONE pre-processed Pagoda2 sample from a raw counts matrix (genes x cells).
# Produces exactly what conos needs: variance normalization + a PCA reduction. Flavor-routed.
preprocess_p2 <- function(cm, n.cores = 1) {
  if (pagoda2_is_devel()) {
    # pagoda2 >= 2.0 (devel): unified constructor + step pipeline.
    pagoda2::Pagoda2$from(cm, n.cores = n.cores, verbose = FALSE)$run(
      steps = c("variance", "pca"), verbose = FALSE)
  } else {
    # pagoda2 1.x (CRAN/main): classic constructor + the two explicit steps conos needs.
    # (On devel these legacy methods still work but warn — that's why we branch, not unify.)
    p <- pagoda2::Pagoda2$new(cm, n.cores = n.cores, log.scale = TRUE, verbose = FALSE)
    p$adjustVariance(plot = FALSE, verbose = FALSE)
    p$calculatePcaReduction(nPcs = 30, n.odgenes = 2000, verbose = FALSE)
    p
  }
}

# `cms` is a named list of raw count matrices (genes x cells), one per sample.
# Names of the list become the sample IDs.
samples <- lapply(cms, preprocess_p2)
```

Seurat objects (v3/v4/v5) work too — build each the usual way
(`NormalizeData`→`FindVariableFeatures`→`ScaleData`→`RunPCA`) and put them in the list;
conos aligns on the PCA (the Seurat path is pagoda2-flavor-independent). Samples can also be
read from files; on **devel** pagoda2 the `from*()` constructors read `.h5ad`/`.h5seurat`/
loom/lstar directly, while **CRAN (1.x)** reads only 10x natively (`read10xMatrix()`) — the
branched `read_sample_p2()` helper in `references/data_sources.md` handles both. A panel may
MIX object types.

**Report:** number of samples, and per-sample cells × genes. To try the recipe without
your own data, the bundled `conos::small_panel.preprocessed` is a ready 2-sample panel
(`data("small_panel.preprocessed", package = "conos")`).

For panel-building variants (Seurat, files, mixed, `planIntegration()`), read
`references/data_sources.md`.

---

## Step 2 — Build the joint graph

`Conos$new()` collects the panel; `runGraph()` aligns every pair of samples and assembles
the joint graph over the union of cells.

```r
con <- Conos$new(samples, n.cores = 4)
con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 30, n.odgenes = 2000, verbose = FALSE)
```

- `k` / `k.self` — inter-sample and within-sample neighbours.
- `space = "PCA"` — reciprocal PCA (default). `alignment.strength = <0..1>` pulls samples
  together more forcefully (default: gentle — leave unset unless needed).
- `runGraph()` is the recommended verb; `buildGraph()` is a deprecated alias.

**Report:** the joint graph is in `con$graph` (an igraph over all cells). Print
`length(igraph::V(con$graph))` (total cells) and confirm `!is.null(con$graph)`.

For all `runGraph` parameters (SNN, `pairs.storage`, CPCA/CCA), read
`references/workflow_and_graph.md`. For alignment tuning, read
`references/alignment_strength.md`.

---

## Step 3 — Joint clustering + embedding, then JUDGE the alignment

Cluster and embed the joint graph, then colour the embedding **by sample** — that is how
you judge alignment: well-aligned samples interleave; residual batch shows as
sample-specific regions.

```r
con$runClustering(verbose = FALSE)                                   # Leiden communities
con$runEmbedding(method = "UMAP", min.visited.verts = 100, verbose = FALSE)  # joint 2-D embedding
```

```r
p_clusters <- con$plotGraph(clustering = "leiden")
ggplot2::ggsave("joint_umap_clusters.png", p_clusters, width = 7, height = 6, dpi = 120, bg = "white")

p_samples <- con$plotGraph(color.by = "sample", mark.groups = FALSE, alpha = 0.1, show.legend = TRUE)
ggplot2::ggsave("joint_umap_samples.png", p_samples, width = 7.5, height = 6, dpi = 120, bg = "white")

# the same joint clustering laid over each sample separately, on the JOINT embedding:
# use.common.embedding = TRUE is required here — conos pre-processing builds only variance + PCA
# (Step 1), so the samples have NO per-sample embedding; FALSE would error ("No 'tSNE' embedding").
p_panel <- con$plotPanel(clustering = "leiden", use.common.embedding = TRUE)
ggplot2::ggsave("panel_clusters.png", p_panel, width = 9, height = 8, dpi = 120, bg = "white")
```

**Report:** cluster count (`nlevels(con$clusters$leiden$groups)`), largest/smallest cluster
sizes, and — from the sample-coloured plot — whether samples interleave or separate. If
they separate by batch, raise `alignment.strength` (`references/alignment_strength.md`)
and re-run from Step 2.

`runClustering` accepts a string method (`"leiden"`, `"walktrap"`, …) or a function;
`runEmbedding` does largeVis (default) or UMAP. Resolution scanning is in
`references/workflow_and_graph.md`.

---

## Step 4 — Joint markers

`runMarkers()` finds genes distinguishing the joint clusters, computed **consistently**
across the panel (pagoda2 samples use pagoda2's disk-backed path; Seurat/other samples use
the same Wilcoxon core via `sccore::matrixDE`). `plotMarkerDotPlot()` shows specific
per-cluster markers.

```r
de <- con$runMarkers(z.threshold = 3.0, upregulated.only = FALSE, verbose = FALSE)

p_dot <- con$plotMarkerDotPlot(n.genes.per.group = 3, min.auc = 0.6)
ggplot2::ggsave("marker_dotplot.png", p_dot, width = 12, height = 9, dpi = 120, bg = "white")

# saveDEasCSV(de.results, saveprefix) writes ONE csv per cluster: paste0(saveprefix, cluster, ".csv").
# There is no `file=` argument; the 2nd arg is a path PREFIX, not a single filename.
saveDEasCSV(de, "cluster_markers_")   # -> cluster_markers_<cluster>.csv, one per joint cluster
```

**Report:** per-cluster top markers; whether they are cluster-specific or dominated by
broad/house-keeping genes (raise `min.auc`/`z.threshold` if so). `runMarkers()` is the
recommended verb (`getDifferentialGenes()` is deprecated). For ranking modes and the
mixed-panel DE details, read `references/markers.md`.

---

## Step 5 — Transfer labels across samples (optional)

If one sample (or a subset of cells) is annotated, `propagateLabels()` diffuses those
labels over the joint graph and reports a per-cell confidence.

```r
# cellannot: a NAMED factor of known labels (names = cell IDs of the annotated sample(s)).
res <- con$propagateLabels(labels = cellannot, method = "diffusion")
# res$labels (per-cell label), res$uncertainty (1 - max posterior), res$label.distribution

p_lab <- con$plotGraph(groups = res$labels)
ggplot2::ggsave("label_transfer_umap.png", p_lab, width = 7, height = 6, dpi = 120, bg = "white")
```

**Report:** fraction of cells confidently labelled (e.g. `mean(res$uncertainty < 0.25)`),
and which clusters got which labels. High-uncertainty regions are where samples lack a
good cross-mapping. For solver vs diffusion and patterns, read
`references/label_transfer.md`.

---

## Step 6 — Save the integrated object

```r
saveRDS(con, "conos_integrated.rds")   # native R object for continuation
```

The whole collection (samples + joint graph/clustering/embedding) also exports to an lstar
zarr store, a Seurat v5 split-assay object, or an AnnData for scanpy — there is no
corrected matrix, the joint graph travels natively. See `references/export_and_interop.md`
(and the lstar recipe).

**Report:** the saved `.rds`, the marker CSV, and the figures produced.

---

## Batch variant

For `args == "batch"` (an orchestrator integrating each of several panels): skip the
per-step figures and "Report" footers, still build the graph/clustering/embedding and save
`conos_integrated.rds` (+ `cluster_markers.csv`), and print ONE summary line
(`samples=… cells=… clusters=… markers=…`). Do not emit `joint_umap_*`/`panel_clusters`/
`marker_dotplot` PNGs in batch mode unless explicitly requested — a multi-panel run would
otherwise flood the chat with figures the user never asked for.

---

## Final response checklist

When the integration is complete, summarize in order:

- Panel: number of samples, per-sample object type (pagoda2 / Seurat / mixed), cells × genes.
- Alignment: `space`, `k`, any `alignment.strength`/supervised setting, and whether the
  sample-coloured embedding shows good mixing.
- Clustering/embedding: method, cluster count, largest/smallest cluster sizes.
- Markers: top per-cluster markers used; whether broad genes dominated.
- Label transfer (if run): reference sample, fraction confidently labelled.
- Saved outputs: `conos_integrated.rds`, `cluster_markers.csv`, and the figures.
- Caveats: residual batch, over-mixing risk, no corrected-expression output, small clusters.

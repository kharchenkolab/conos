# Markers: runMarkers and plotMarkerDotPlot

Find genes that distinguish the joint clusters, computed **consistently across the panel** regardless of
how each sample was pre-processed.

## `runMarkers()`

```r
con$runMarkers(clustering = NULL, groups = NULL, z.threshold = 3.0, upregulated.only = FALSE,
               verbose = TRUE, append.specificity.metrics = TRUE, append.auc = TRUE)
```

- `clustering` — name of a clustering in `con$clusters` (default: the most recent / default clustering).
- `groups` — a custom per-cell factor instead of a stored clustering.
- `z.threshold = 3.0` — minimum rank-based Z to report.
- `upregulated.only = FALSE` — keep down-regulated genes too.
- `append.auc`, `append.specificity.metrics` — add AUC and precision/expression-fraction columns
  (used by the balanced dot-plot selection).

Returns (and stores) a per-cluster data frame: gene, `Z`, `M` (log2 fold-change), `AUC`, specificity
metrics, adjusted p-value. The DE core is a **Wilcoxon / Mann-Whitney rank test**: rank-based Z plus log2
fold change, BH-adjusted in log space.

## How markers are computed (per-sample, then aggregated)

For each sample, conos computes the per-sample Wilcoxon test, then aggregates Z across samples. The
per-sample step is **class-aware**:

- **pagoda2 samples** use pagoda2's own optimized, disk-backed `runMarkers()` path.
- **Seurat / other samples** are scored with the shared `sccore::matrixDE()` — the *same* Wilcoxon Z /
  log-fold-change core (verified to reproduce pagoda2's numbers exactly).

So a **mixed panel** is scored consistently: the result does not depend on whether a sample is a pagoda2
or Seurat object. The input to DE is each sample's normalized (log) layer, accessed via
`getCountMatrix(sample, transposed = TRUE)` (disk-backed where the backend is) — not a whole-matrix grab.

## `plotMarkerDotPlot()` — specific per-cluster markers

```r
con$plotMarkerDotPlot(clustering = NULL, groups = NULL, n.genes.per.group = 5,
                      z.threshold = 1, min.auc = 0.6, cols = c("grey88", "firebrick3"),
                      dot.scale = 6)
```

Genes are selected with pagoda2.1's **"balanced" rule**, not raw `Z`:

- keep genes that are up-regulated **and** discriminative (`Z >= z.threshold` **and** `AUC >= min.auc`),
- rank by the **harmonic mean of precision and expression fraction** (an F1-like score), not raw `Z`,
- assign each gene to its single best cluster (dedup), and order genes by cluster.

This keeps ubiquitous house-keeping / mitochondrial genes from dominating the plot (raw-Z selection picks
them up). Rendered through `sccore::dotPlot`. Works on pagoda2, Seurat, and mixed panels.

If a dot plot looks dominated by broad genes, raise `min.auc` and/or `z.threshold`.

## Saving

```r
de <- con$runMarkers(z.threshold = 3.0, upregulated.only = FALSE)   # runMarkers, not the deprecated getDifferentialGenes
saveDEasCSV(de, "markers_")   # writes ONE csv per cluster: paste0(prefix, cluster, ".csv") -> markers_<cluster>.csv
```

`saveDEasCSV(de.results, saveprefix, gene.metadata = NULL)` — the 2nd arg is a path PREFIX (it appends
`<cluster>.csv`), **not** a `file=` filename, and there is no `file=` argument. Compute with
`runMarkers()`; `getDifferentialGenes()` is the deprecated alias.

## Pseudobulk per cluster

```r
mats <- con$getClusterCountMatrices(clustering = "leiden", common.genes = TRUE, omit.na.cells = TRUE)
# list of per-sample dense matrices (genes x clusters) — pseudo-bulk sums; useful for cross-sample DE.
```

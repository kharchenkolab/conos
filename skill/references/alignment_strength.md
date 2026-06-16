# Adjusting alignment strength

By default conos aligns **gently** — it corrects obvious batch differences while preserving each sample's
structure. When samples from different tissues/conditions still form sample-specific sub-clusters within a
shared cell type, and that residual separation is batch rather than biology, pull them together more
forcefully. Two knobs, both on `runGraph()`.

**Caution:** stronger alignment trades structure for mixing. Push too far and genuinely distinct
populations merge. Always re-check that the populations you care about stay separated (re-run
`runClustering()` + `plotMarkerDotPlot()` after any change).

## How to judge alignment

Colour the joint embedding **by sample**: well-aligned samples interleave; residual batch shows up as
sample-specific regions.

```r
con$plotGraph(color.by = "sample", mark.groups = FALSE, alpha = 0.1, show.legend = TRUE)
```

## 1. `alignment.strength` (global dial, 0 → 1)

0 = the gentle default; 1 = maximum. Higher values pull mutual neighbours across samples together more
aggressively.

```r
con$runGraph(alignment.strength = 0.3)   # moderate; 0.2–0.4 is usually plenty
con$runGraph(alignment.strength = 0.7)   # strong — watch for distinct populations collapsing
```

Build into separately named embeddings to compare:

```r
con$runGraph(alignment.strength = 0.3)
con$runEmbedding(method = "UMAP", embedding.name = "strength_0.3", min.visited.verts = 100)
con$plotGraph(embedding = "strength_0.3", color.by = "sample")
```

## 2. Supervised alignment (target a named factor)

When you know *which* factor the unwanted separation follows (tissue, donor, batch), down-weight edges
*within* that factor so alignment works hardest *across* it — more surgical than turning up global
strength.

```r
# a per-cell factor over the panel's cells (here: tissue from the sample names)
tissue <- ifelse(grepl("BM", con$getDatasetPerCell()), "BM", "CB")
tissue <- setNames(factor(tissue), names(con$getDatasetPerCell()))

con$runGraph(balancing.factor.per.cell = tissue,   # the factor to balance across
             same.factor.downweight = 0.1,          # < 1: scale down within-factor edges
             alignment.strength = 0.3)
```

`balancing.factor.per.cell` names the factor; `same.factor.downweight` (< 1) scales down within-factor
edges, so cross-factor links carry relatively more weight. There is also `balancing.factor.per.sample`
for a sample-level factor, and `balance.edge.weights = TRUE` to equalize per-sample edge mass.

## Choosing a value

- Start from the **default** (no `alignment.strength`). Only reach for these when the sample-coloured plot
  shows batch-like separation the default didn't remove.
- Raise `alignment.strength` gradually (0.2–0.4 usually suffices).
- Prefer **supervised** alignment when you can name the nuisance factor — it corrects that axis without
  over-mixing everything.
- After any change, re-cluster and re-check markers; if previously distinct populations merged, back off.

# Workflow: Conos$new → runGraph → runClustering → runEmbedding

The joint analysis is a sequence of methods on one `Conos` object. Each mutates the object in place and
stores results under canonical names; later steps read what earlier ones produced.

```r
con <- Conos$new(samples, n.cores = 4)   # samples: named list of Pagoda2 / Seurat objects
con$runGraph()
con$runClustering()
con$runEmbedding()
```

## `Conos$new(x, n.cores = ..., verbose = TRUE)`

`x` is a named list of per-sample objects (`Pagoda2` or `Seurat`), each already normalized with a PCA
reduction. Names become the sample IDs (`getDatasetPerCell()`); if unnamed, samples are auto-named.
`n.cores` sets the default thread budget for the pairwise alignments.

## `runGraph()` — build the joint graph

Aligns every pair of samples, then assembles inter- + intra-sample edges into one graph. Key parameters
(defaults shown):

- `k = 15` — inter-sample neighbours per cell (the cross-sample mapping strength).
- `k.self = 10`, `k.self.weight = 0.1` — intra-sample neighbours and their relative weight.
- `space = "PCA"` — alignment space. `"PCA"` = **reciprocal PCA** (default, fast, robust); `"CPCA"` =
  common PCA; `"CCA"` / `"PCA"` variants opt-in; `"genes"` aligns in gene space. Use the default unless
  you have a reason.
- `matching.method = "mNN"` — mutual nearest neighbours; `metric = "angular"`.
- `ncomps = 30`, `n.odgenes = 2000` — components and overdispersed genes used for alignment.
- `alignment.strength = NULL` — 0→1 dial to pull samples together more forcefully (default = gentle). See
  `alignment.md`.
- `balancing.factor.per.cell`, `same.factor.downweight = 1.0` — *supervised* alignment: down-weight edges
  within a named factor (tissue/donor/batch) so alignment works across it. See `alignment.md`.
- `snn = FALSE`, `snn.quantile = 0.9` — optionally build a shared-NN graph.
- `pairs.storage = c("keep","drop","disk")` — keep (default), drop, or offload the O(n²) per-pair
  rotations to disk to save memory on large panels (`"disk"` is restored on the next build).
- `exclude.samples = NULL`, `verbose = TRUE`.

The graph lands in `con$graph` (an igraph object over the union of cells).

## `runClustering()` — joint communities

`runClustering(method = leiden, min.group.size = 0, name = NULL, test.stability = FALSE, resolution = ...,
verbose = TRUE)`.

- `method` accepts a **function** (`leiden.community`, `walktrap.community`, …) **or a string**:
  `"leiden"`, `"walktrap"`, `"louvain"`/`"multilevel"`, `"infomap"`, `"fastgreedy"`, `"labelprop"`,
  `"leadingeigen"`. Default is Leiden (with `n.iterations = 5`).
- `name` — store the result under a custom name (default is the method name, e.g. `leiden`); results live
  in `con$clusters[[name]]$groups` (a named factor over all cells).
- `test.stability = TRUE` — bootstrap cluster stability (`stability.subsamples`,
  `stability.subsampling.fraction`); inspect with `con$plotClusterStability()`.

Scan resolution / k:

```r
scanResolution(con, resolutions = seq(0.1, 2, by = 0.1), method = leiden.community, plot = TRUE)  # clusters & modularity vs resolution
scanKModularity(con, min = 3, max = 50, by = 1, plot = TRUE)                                       # modularity vs graph k
```

## `runEmbedding()` — joint 2-D embedding

`runEmbedding(method = "largeVis", embedding.name = method, M = 1, gamma = 1, alpha = 0.1,
perplexity = NA, sgd_batches = 1e8, seed = 1, target.dims = 2, verbose = TRUE, ...)`.

- `method = "largeVis"` (default, fast) or `"UMAP"`. For UMAP, useful extras: `min.visited.verts`,
  `min.prob`, `n.neighbors` (passed through). UMAP runs on the commute-time neighbours of the joint graph.
- `embedding.name` — store under a name so several embeddings coexist (`con$embeddings[[name]]`);
  plot with `plotGraph(embedding = name)`.
- largeVis embeddings can be slow with default `sgd_batches`; lower for quick previews.

## Plotting

```r
con$plotGraph(color.by = "cluster")                 # color.by: "cluster" | "sample" | gene (with gene=)
con$plotGraph(color.by = "sample", mark.groups = FALSE, alpha = 0.1, show.legend = TRUE)
con$plotGraph(gene = "CD3E")                          # expression of one gene on the joint embedding
con$plotGraph(clustering = "leiden", embedding = "umap")
con$plotPanel(clustering = "leiden")                  # one facet per sample, shared layout
```

`plotGraph` args of note: `color.by`, `clustering`, `embedding`, `groups` (a custom factor), `gene`,
`subset`, `mark.groups`, `alpha`, `show.legend`, `plot.theme`.

## Deprecated aliases (still work, warn)

| deprecated            | use instead       |
|-----------------------|-------------------|
| `buildGraph()`        | `runGraph()`      |
| `findCommunities()`   | `runClustering()` |
| `embedGraph()`        | `runEmbedding()`  |
| `getDifferentialGenes()` | `runMarkers()` |

## Memory / performance

- `pairs.storage = "drop"` or `"disk"` frees the per-pair alignment rotations on large panels.
- `getJointCountMatrix()` merges sample matrices **sparsely** (no dense zero-block when gene sets differ).
- lstar-backed pagoda2 samples stream during `runGraph` (block reads) — peak memory stays flat.
- `n.cores` parallelizes the pairwise alignments (the dominant cost).

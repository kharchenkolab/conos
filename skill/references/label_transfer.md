# Label transfer: propagateLabels

Transfer cell-type (or any categorical) labels from labelled cells to the rest of the panel by diffusing
them over the **joint graph**. This is conos's headline cross-sample capability: annotate one (or a few)
samples, then propagate to all.

## API

```r
res <- con$propagateLabels(labels, method = "diffusion", ...)
# res$labels            : named character vector — the inferred label for EVERY cell in the panel
# res$uncertainty       : named numeric in [0,1] — 1 - max posterior (higher = less confident)
# res$label.distribution: matrix (cells x labels) — full posterior per cell
```

- `labels` — a **named** vector/factor of known labels; names are cell IDs (a subset of the panel's
  cells, e.g. one annotated sample's cells). Unlabelled cells are inferred.
- `method`:
  - `"diffusion"` (default) — label diffusion over the graph (`propagateLabelsDiffusion`). Robust general
    choice.
  - `"solver"` — a linear-solver formulation (`propagateLabelsSolver`); pass `solver = "Matrix"`
    (default) or `solver = "mumps"` (needs `rmumps`). Use when you want the harmonic-function solution.

## Typical pattern

```r
# cellannot: named factor of annotations for the cells of one reference sample
res <- con$propagateLabels(labels = cellannot, method = "diffusion")

con$plotGraph(groups = res$labels)                 # view propagated labels on the joint embedding
con$plotGraph(colors = res$uncertainty)            # where is the transfer least confident?

# keep only confident calls
confident <- res$labels[res$uncertainty < 0.25]
```

## Notes

- Propagation runs on `con$graph`, so `runGraph()` must have been called. Embedding/clustering are not
  required for label transfer (but help you inspect it).
- The reference (labelled) cells may miss some label levels relative to others — the label-indicator
  matrix is sized to the number of levels, so this no longer errors.
- `method = "solver", solver = "Matrix"` does not spuriously warn about `rmumps` (that warning is raised
  only when the `"mumps"` solver is actually requested but unavailable).
- Inspect `res$uncertainty` before trusting calls; high uncertainty regions are where samples don't have a
  good cross-mapping (rare types present in only one sample, or genuine novelty).

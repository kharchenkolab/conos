# Conos upgrade proposal

*Draft 2026-06-14 (rev. 2026-06-15). Target: conos `dev` (currently v1.5.4) → **conos 2.0**, shipped in a single coordinated CRAN wave (lstar first, then pagoda2 2.0 + conos 2.0 together). No intermediate conos 1.6 CRAN release.*

## 0. Goals & guiding principles

Conos is the integration boundary of the kharchenkolab single-cell stack: it turns a
collection of per-sample [pagoda2](https://github.com/kharchenkolab/pagoda2) (or Seurat)
objects into one joint graph, clustering, and embedding, and it is the **only** surface
through which `cacoa` reaches pagoda2. pagoda2 has just been rewritten (2.0.0: R6, multimodal
"facets", generic `runX(method=)` verbs, RcppHNSW kNN, RSpectra SVD, lstar Zarr disk-backing,
`sccore`-based plotting). This proposal brings conos in line and addresses its three biggest
pain points, in priority order: **memory use, speed, and integration quality.**

**Guiding principles**

1. **Memory, speed, integration quality first.** These are the levers that decide whether
   conos is usable at atlas scale and whether the result is trustworthy. A concrete target:
   **the whole conos pipeline should be runnable in a disk-backed, economical-RAM mode** —
   peak memory bounded by the joint graph + embedding + the current pair's working set, not
   by holding every sample's matrices (or all O(n²) rotations) in RAM at once.
2. **Additive-with-deprecation on the public API.** conos 2.0 is a major version, so default
   flips and the verb realignment are fair game — but do them *gracefully*: ship the new `runX`
   verbs as preferred with the **old names kept as deprecated-but-working aliases** (`lifecycle`/
   `.Deprecated()`), and reserve outright *removals* (`.Defunct`) for the next major. cacoa is off
   CRAN and migrates against the deprecated-but-working API at its own pace. So "conservative" is
   now a user-experience / cacoa-migration choice, not a CRAN constraint.
3. **Default flips land in 2.0, clearly noted.** Because this is the major bump, the decided
   default changes (e.g. `space`→CPCA, Leiden iterations) ship here with a prominent NEWS entry —
   no need to defer them further.
4. **One coordinated wave; lstar first.** Submit **lstar** (the leaf dependency) first, then
   **pagoda2 2.0 + conos 2.0 together** (co-dependent, noted in each `cran-comments.md`). The
   pagoda2.1 accessor fixes (§4) are still a hard prerequisite — pagoda2 2.0's revdep check needs a
   conos that works against the new pagoda2. (`sccore` may also need to go in the leaf tier — see §11.)

**Release framing**

A **single coordinated CRAN wave**, not a bridge-then-break sequence:

- **lstar 1.x → CRAN first** — leaf dependency; gates the wave date (currently 0.0.1 with known
  issues, see §11).
- **pagoda2 2.0 + conos 2.0 → CRAN together.** conos 2.0 carries *everything ready and worth
  shipping*: the correctness/accessor fixes, the disk-backed economical mode, speed, graphics
  consistency, the decided default flips, CRAN hygiene, **and** the API realignment done additively
  (runX preferred, old names deprecated-but-working).
- **The non-breaking "1.6" set is NOT a separate CRAN release** — at most a `dev` checkpoint tag
  (not merged to main). It folds into conos 2.0.
- **Deferred to a later conos release** (post-wave, GitHub-first, coordinated with the cacoa
  migration): removing the deprecated aliases, the DE unification (breaking output shape, couples to
  cacoa), and the unbaked integration-quality research (full multimodal/WNN, joint >2-sample
  reductions, scalable solver).

---

## 1. Memory

**The real main path** (per the maintainer) is: *preprocessed collection → cross-pair rotations →
joint graph → cluster → embed → markers → propagate/diffuse labels.* The joint count matrix and
`correctGenes` "adjusted expression" are **off** this path and rarely used (adjusted-expression
values are typically not useful), so the earlier "dense joint matrix" finding is **not** the
main-path RAM driver. On the actual hot path, peak RAM is held by:

- **(A) the sample objects themselves** — every sample's counts + reductions are resident unless the
  sample is disk-backed. For large in-memory collections this is likely the dominant consumer.
- **(B) O(n²) cross-pair rotations** (`self$pairs`), retained for the object's lifetime.
- **(C) per-pair alignment intermediates × fork copies** during `buildGraph` (scaled matrices,
  covariances), plus transient dense adjacency rebuilds.

The unifying fix is the **disk-backed economical mode** of §5 (stream A, disk-back B, free C).
Items ranked for that goal:

| # | Item | Evidence | Fix | Effort |
|---|---|---|---|---|
| **1.1** | **(A) Sample matrices held resident.** conos holds samples by reference (not duplicated — good), but the samples' own count/reduction matrices stay in RAM. The accessor path is already half-wired to stream from pagoda2.1/lstar; the requirement is that conos **never materializes a whole sample matrix** (read by gene/cell block) so disk-backed samples actually cut peak RAM. | `access_wrappers.R:24-73`; densifying consumers `R/conclass.R:871-873,951` | Guarantee block/streamed reads end-to-end; see §5. **The main lever for economical RAM.** | Med |
| **1.2** | **(B) O(n²) pairwise rotations retained.** `self$pairs[[space]]` keeps one rotation per pair (CPCA `CPC`; CCA `u/v/ul/vl`; JNMF `rot1/rot2/z`) for the object's lifetime — sometimes useful for re-runs, but a real drag (~5k pairs at 100 samples). | `R/conclass.R:1062-1092` | **Three-level `pairs.storage = c("keep","disk","drop")`** — keep (current), disk-back the rotations (the sweet spot for large panels), or drop after the graph is built. | Low–med |
| **1.3** | **(C) Per-pair intermediates × fork copies + transient adjacency.** `papply`'s `sapply(res, class)` materializes all results in the parent at once; CPCA stacks `k` dense p×p covariance matrices (`abind` cube); several methods rebuild a full n×n adjacency on the fly. | `R/conos.R:410,186,220-242`; `R/conclass.R:689,362`; `src/spcov.cpp:14` | Free intermediates promptly; lighter `papply` error check / chunking; accumulate CPCA covariance incrementally; avoid repeated adjacency materialization. | Low–med |
| **1.4** | **(latent) Dense `getJointCountMatrix`.** `mergeCountMatrices`/`extendMatrix` are defined **twice** (sparse in `conos.R`, dense in `integrations.R`); with **no `Collate:`** the dense version wins, so the *rarely-used* export / DE-specificity / ScanPy / velocity paths densify. **Off the main path, so not the headline** — but a cheap correctness fix. | `R/integrations.R:6-26` vs `R/conos.R:1029-1049` | Delete the dense dup (or add `Collate:`); keep sparse. | Trivial |

**Net:** the economical-RAM win comes from §5 (stream samples, disk-back/drop rotations), not from
the joint-matrix bug. 1.4 stays as a cheap fix for the export paths that do use it.

---

## 2. Speed

Parallelism in conos is **across sample pairs** — the pairwise loop forks over `combn(samples, 2)`
via `papply`/`mclapply`. So when #pairs ≫ #cores the inner per-pair kNN/SVD should **not** also
thread (that would oversubscribe). The levers, with that in mind:

| # | Item | Evidence | Fix | Effort |
|---|---|---|---|---|
| **2.1** | **kNN backend.** `N2R::Knn`/`crossKnn` are called with `nThreads` hardcoded to `1` (and N2R 1.0.5's threading is a no-op anyway). Because alignment is parallelized **across** pairs, per-pair threading is usually unnecessary — but the **within-sample** kNN (`getLocalNeighbors`, once per sample) and the **few-samples / large-cells** regime *do* leave cores idle. | `R/conos.R:596,852-856` | Pluggable `.conos_knn_sparse()` (RcppHNSW preferred, N2R fallback) — for consistency with pagoda2, better recall, faster single-thread; **thread only where pair-level parallelism doesn't already saturate cores** (within-sample step, few-pair panels). Not a blanket "8× via threading." | Med |
| **2.2** | **SVD is all `irlba`.** | `R/conos.R:185,288,344`; `R/integrations.R:313` | Swap to `RSpectra::svds` (near drop-in, consistent with pagoda2.1) for the per-pair decompositions; offer randomized SVD for the few-pass disk-backed case. | Low |
| **2.3** | **O(n²) pairwise alignment** with no subsampling (the `# TODO: add random subsampling for very large panels` is unimplemented). This — not within-pair threading — is the real scaling wall for large panels. | `R/conclass.R:1045,1058` | Landmark/anchor-sample or subsampled-pair scheme; the pair cache already supports incremental population. | Med–high (later) |
| **2.4** | **Per-pair re-decomposition — mostly intentional.** Recomputing per pair is *correct* in the general case: each pair should be compared the most sensible way for *that* pair (common genes between those two samples — gene sets can differ; a pair-specific CCA; etc.). The only redundancy is the **simple common case**: identical gene space across datasets + PCA/CPCA, where the per-sample PCA could be reused (conos already pulls the dataset's stored PCA for the *within*-sample step). | `R/conos.R:288`, `access_wrappers.R:9` | In the common-gene + PCA/CPCA case, reuse the precomputed per-sample reduction for the inter-sample step too; keep per-pair recomputation as the default general path. | Low–med |

---

## 3. Integration quality (key)

**Pipeline.** For each sample pair, `buildGraph` computes a shared low-dim space
(`space=`), projects both samples in, and links cells by mutual NN (mNN); within-sample kNN
edges are added at `k.self.weight = 0.1` (10× weaker) so cross-sample links drive community
structure. Edges are summed, optionally SNN-reweighted / balanced, then Leiden-clustered and
embedded. The quality levers and their current defaults:

| Lever | Default | Effect / concern |
|---|---|---|
| `space` | **`PCA`** (current) → **`CPCA`** (decided) | `PCA` is *concatenated per-pair PCAs*, the weakest aligner — correspondence rests entirely on mNN. **`CPCA`** (common principal components, `src/cpca.cpp`) is the principled shared-subspace choice and the original-paper recommendation — **the agreed default.** `CCA` aligns *too aggressively* (over-merges distinct populations) and stays **opt-in, not a default**. |
| `ncomps` | **40** | **Per-pair** alignment dimensionality (the shared subspace between *two* samples) — *not* a whole-atlas representation; atlas-scale resolution comes from the *graph* connectivity across all pairs, not from this rank. So 40 is reasonable and raising it is not clearly motivated. (Minor: `quickCPCA`/`quickCCA` default internally to 100 but `buildGraph` truncates to 40 — an inconsistency worth reconciling for clarity, not a quality problem.) |
| `k` / `k.self` / `k.self.weight` | 15 / 10 / 0.1 | Mixing strength vs structure preservation; the key tuning axis. |
| `matching.method` | `mNN` | mNN (conservative) vs NN (permissive). |
| `n.odgenes` | 2000, **intersected** across all samples | Heterogeneous panels shrink the shared feature set (`commonOverdispersedGenes`). |
| Leiden | `resolution=1`, **`n.iterations=2`** | Low iterations; single fixed resolution; modularity objective (resolution-limit problem on large graphs). |
| `snn` / `balance.edge.weights` | FALSE / FALSE | Both measurably improve type-coherence / confounder control but are off and under-documented. |

**Quality-improvement suggestions** (the heart of the upgrade):

- **3.1 — Better defaults (behavior-affecting; stage carefully).** **Switch the default `space`
  `PCA` → `CPCA`** (*decided* — CPCA is the principled shared-subspace aligner and the original-paper
  default; `CCA` over-merges and stays opt-in, never a default). Also bump Leiden `n.iterations`
  2 → ~5–10 and surface `resolution` + a resolution-sweep helper. (`ncomps=40` is **per-pair** and
  fine — see the levers table — so it is *not* on this list.) *These ship in conos 2.0 (the major
  bump) with a prominent NEWS entry flagging the reproducibility change.*
- **3.2 — Reuse the precomputed reduction in the common case (not a blanket change).** Per-pair
  recomputation is the *right* behavior in general — each pair should be aligned the most sensible
  way for that pair (common genes between those two samples, a pair-specific CCA, etc.). The narrow
  win is the **simple common case**: identical gene space across datasets + PCA/CPCA, where the
  per-sample PCA (which conos already pulls for the within-sample step, `access_wrappers.R:9`) can
  also serve the inter-sample step — avoiding the redundant per-pair PCA and making the two edge
  types geometrically consistent. Implementation: let `getPca` accept a reduction name (§4) and add
  a fast path in `getPcaBasedNeighborMatrix` for the precomputed-reduction case; **do not** force it
  on the CCA / differing-gene paths.
- **3.3 — Multimodal / WNN integration (biggest capability gap; deferred research, post-wave).** Conos has **no**
  collection-level multimodal path (`grep` for wnn/facet/atac/adt is empty). pagoda2.1 now produces
  per-sample WNN / `runReduction(facets=...)` joint reductions. Add a path to align samples on a
  joint multimodal reduction (accept per-sample cell-embeddings as alignment coordinates; multimodal
  or per-modality-combined mNN).
- **3.4 — Integration QC, including over-integration.** Conos already has
  `estimateWeightEntropyPerCell` (batch-mixing entropy), `scanKModularity`, and stability tests, but
  they're scattered and opt-in, and only detect *under*-mixing. Package a one-call "integration
  report" (per-cell + per-cluster mixing entropy, cross-sample mNN edge fraction) **and** an
  over-integration diagnostic (local label-purity vs mixing trade-off) so users can tune
  `space`/`k`/`alignment.strength` principledly.
- **3.5 — Scalable high-quality label propagation.** `propagateLabels` has a high-quality solver
  (ZGL harmonic) explicitly flagged "inappropriate for >20k cells" and a `TODO` to replace it; large
  data silently falls to the lower-quality diffusion method. Add a Laplacian-aware iterative solver
  (preconditioned CG) so quality scales.
- **3.6 — Joint >2-sample reductions / global anchoring (deferred research, post-wave).** Move beyond strictly
  pairwise spaces toward a shared latent basis across the whole panel (multi-sample CPCA generalizes
  naturally), reducing reliance on emergent graph connectivity for cross-panel consistency.

---

## 4. pagoda2.1 alignment & the accessor contract (correctness — revdep blocker)

The coupling lives in `R/access_wrappers.R`. It has been **partially** patched (raw counts /
expression / od-genes route through pagoda2.1 accessors via method-probing — good). Still broken
against pagoda2.1:

| conos need | pagoda2.1 accessor | current conos call | status |
|---|---|---|---|
| raw counts | `getRawCounts(...)` | method-probed | OK |
| normalized expr | `getExpressionBlock(..., scale.variance=)` | method-probed (but variance scaling hand-rolled from `misc$varinfo`) | OK / improve |
| od genes | `getOdGenes(n)` | `sample$getOdGenes(n)` | OK |
| reduction | *(no getter — request `getReduction()` upstream)* | `sample$reductions$PCA` (raw, fragile) | **fragile** |
| embedding | `getEmbedding(type, name)` | `sample$embeddings$PCA[[type]]` (wrong key) | **broken** |
| clustering | `getGrouping()` / `cellMeta` | `sample$clusters$PCA[[type]]` (wrong key) | **broken** |
| (legacy) | — | `p2app4conos`/`convertToPagoda2` touch removed `$counts` | **throws** |
| facet/modality | `getFacet`/`listFacets`/`getFacetMembership` | none | **absent** |

**Suggestions:**
- **4.1 (P0)** Fix the broken keys: `getEmbedding` → `sample$getEmbedding(type, name)` /
  `embeddings[["counts"]][[name]]`; `getClustering` → `getGrouping`/`cellMeta`; `getPca` →
  method-probe a reduction getter, fall back to `reductions[[sample$defaults$reduction]]` (not
  hardcoded `$reductions$PCA`). Remove `$counts` reach-ins from `p2app4conos`/`convertToPagoda2`.
  Refs: `access_wrappers.R:9,358,393`; `integrations.R:421-491,334-361`.
- **4.2 (P1)** Replace the hand-rolled variance scaling in `scaledMatricesP2` with
  `getExpressionBlock(scale.variance=TRUE)` — removes the last reach-in into a pagoda2 internal
  (`R/conos.R:18-46`).
- **4.3 (P1)** Make accessors facet-aware: thread a `facet=` (default = sample's `defaultFacet`)
  through the count/expr/gene/cellname accessors, so "align on the ATAC facet" becomes a parameter
  rather than impossible. (Multimodal *integration* is §3.3; this is just the plumbing.)
- **4.4 (P2)** Formalize a single documented `ConosSampleAdapter` capability contract (the table
  above) to replace the scattered per-class S4 reach-ins (resolves the `# TODO: package-independent
  wrapper`, `conclass.R:110`), centralizing orientation/facet handling.

---

## 5. Disk-backed, economical-RAM mode (a first-class goal)

The target: **the whole conos pipeline runnable with peak RAM ≈ joint graph + embedding + the
current pair's working set**, independent of collection size. This is the unifying frame for the
memory work (§1). Four things must hold simultaneously:

1. **Sample access streams from disk.** The accessors are already **half-wired** — the patched
   `getRawCounts`/`getExpressionBlock` route through pagoda2.1's streaming lstar Zarr reads
   (`access_wrappers.R:24-73`). The requirement is that conos **never materializes a whole sample
   matrix** (read by gene/cell block), so an lstar-backed sample genuinely lowers peak RAM rather
   than being pulled into a dense copy. No new conos dependency — lstar stays optional, on the
   pagoda2 side.
2. **Cross-pair rotations can live on disk** — the `pairs.storage = c("keep","disk","drop")` lever
   (§1.2). "disk" is the default for the economical mode: rotations are written to a scratch store
   and memory-mapped/read on demand when re-building the graph.
3. **Per-pair intermediates are freed promptly** (§1.3) and the pairwise fork fan-out is bounded so
   transient working memory ≈ (cores × one pair), not (all pairs at once).
4. **Off-path densifiers are sparse/streamed** when used (§1.4, `getJointCountMatrix`/pseudobulk).

Deliverable: a `con$buildGraph(..., pairs.storage="disk")`-style economical path plus a documented
"large-collection / disk-backed" recipe, validated to hold peak RSS roughly flat as samples are
added. This is the highest-value memory work and the through-line of the §1 items.

---

## 6. Graphics consistency

Embeddings already go through `sccore::embeddingPlot` (good); everything else is inconsistent:

- **Three theme strategies** coexist — a near-duplicate of `themePagoda2` (`adjustTheme`) + bare
  `theme_bw()` + no theme — so one session yields visually different plots
  (`R/conclass.R:990`; `R/plot.R:189,263,329`).
- **Palette drift:** `rainbow()` in the DE heatmap vs `fac2col` elsewhere (same cluster, different
  color across plots); the DE-heatmap expression palette differs from pagoda2.1's
  (`R/plot.R:387,582`).
- **No marker dot plot** (pagoda2.1 ships `plotMarkerDotPlot` on `sccore::dotPlot`) — the biggest
  "plot the same thing" gap.
- **Robustness:** device-dependent panel sizing via `dev.size()` (`R/plot.R:63`); a base-graphics
  side effect inside `plotClusterStability(what='dend')` that returns `NULL`; a `scanKModularity`
  plot that's built then discarded (`R/conos.R:1021`); per-sample legend duplication in `plotPanel`.

**Suggestions (mostly non-breaking):** adopt one shared theme (`themePagoda2`, ideally moved into
`sccore` so both packages import it); route all categorical colors through `sccore::fac2col`/
`fac2palette` and continuous through `val2col`/`val2ggcol`; add a grouping/palette resolver
(mirroring pagoda2.1's `resolveGrouping`/`resolveFactorColors`) so a cluster has one color across
embedding, barplot, and heatmap; add `plotMarkerDotPlot` as a thin `sccore::dotPlot` wrapper with
pagoda2.1's defaults; replace `dev.size()` with explicit/NULL sizing; fix the dend and
`scanKModularity` robustness bugs; correct the `plot.theme` roxygen type.

---

## 7. API changes (priority-ranked)

**Rule:** in conos 2.0, ship the realignment *additively with deprecation* — new form preferred,
old names kept as deprecated-but-working aliases (`lifecycle`/`.Deprecated()`); reserve outright
**removals** (`.Defunct`) for the next major. Conos is already ~60% aligned with pagoda2.1 in spirit
(it has `method=` on `embedGraph`/`propagateLabels` and name-keyed `con$clusters`/`embeddings`/
`pairs`); the gaps are consistency, not architecture.

**Ships in conos 2.0 (the wave), ranked by value:**

1. **`findCommunities(method=)` should also accept a string** (`"leiden"`/`"walktrap"`/...), not only
   a live function object. Validated enum + internal dispatch; passing the function still works.
   *Highest agent-accessibility win, zero breakage.* (`R/conclass.R:431`)
2. **Fix the asymmetric/broken accessors** (also §4.1): `getEmbedding(con, name)` honor
   `con$embeddings[[name]]` (it currently ignores `name`); add `Conos` methods for
   `getPca`/`getCountMatrix` (currently a bare dispatch error). (`access_wrappers.R:358,380`)
3. **Accept aliased argument names** rather than renaming: add `con` everywhere `con.obj`/`conos.obj`/
   `conos` appear, and `n.cores` where `threads`/`nthreads` appear — old names keep working; document
   the canonical one. Also drop the hardcoded `n.cores=30` in `stableTreeClusters` (a bug).
4. **`runGraph`/`runClustering`/`runEmbedding`/`runMarkers` as the preferred verbs**, with
   `buildGraph`/`findCommunities`/`embedGraph`/`getDifferentialGenes` retained as
   **deprecated-but-working aliases** (`.Deprecated()`). Lets cacoa and users migrate before the
   next-major removal.

**Deferred to the next conos release (post-wave, coordinated with cacoa migration):**

5. Remove the deprecated aliases (`.Defunct`).
6. **Unify DE** — `getDifferentialGenes` + `getPerCellTypeDE` + `getBetweenCellTypeDE` →
   one `runMarkers(method = "wilcoxon"|"pseudobulk")` with one result container (name-keyed under
   `con$markers[[name]]`), consumed uniformly by `plotDEheatmap`/`saveDE*`. Collapses the
   near-duplicate functions and fixes the `saveDEasJSON(gene.metadata=NULL)` bug
   (`R/de_functions.R:285,349,417`). *Breaking output shape — couples to cacoa, so it waits.*
7. Canonicalize argument names with removal; restructure `con$graph` → name-keyed `con$graphs[[method]]`
   (so a CPCA and a CCA graph can coexist, matching `embeddings`/`clusters`).

---

## 8. CRAN-readiness

conos is mid-CRAN-prep (v1.5.4 "fix for CRAN"). Remaining gaps:

- **Undocumented `@return`** on all 11 S4 accessor generics (`access_wrappers.R`) and `saveDEasCSV` —
  likely the biggest WARNING source; fix the malformed `getSampleNamePerCell` example.
- **Undeclared `grid`** dependency (used in `plotDEheatmap`); add to Imports.
- **`ComplexHeatmap` is a hard Bioconductor Import** but only used behind a `requireNamespace` guard;
  move to Suggests and qualify the `importFrom`s.
- **Unguarded heavy examples** in `man/Conos.Rd` run `buildGraph`+`findCommunities` at check time;
  wrap in `\donttest{}`.
- **No shipped vignette** → `knitr`/`rmarkdown`/`ggrastr`/`shinycssloaders` are unused Suggests
  (NOTE). Add a real `vignettes/` + `VignetteBuilder: knitr`, or drop the unused Suggests.
- `utils::globalVariables` gaps (`aRI`/`cluster`/`jc`); two bare `T`/`F` literals; ~1 MB of removable
  `inst/` payload (`scanpy_integration.ipynb` 430 KB).
- Thin test suite (5 `test_that` blocks) — expand with `skip_if_not_installed` for Suggests paths.
- **Maintainer change.** `Authors@R` currently lists the former maintainer (E. Biederstedt, who has
  left) as `cre`. Move `cre` to the current maintainer, sync the `Maintainer:`/`Author:` fields (or
  delete them and let R derive from `Authors@R`), and note the handover in `cran-comments.md` (CRAN
  emails the new maintainer to confirm). Same edit applies to pagoda2 and lstar.

---

## 9. Consolidated roadmap

### conos 2.0 — the wave release (with pagoda2 2.0; after lstar)

Everything ready and worth shipping goes here. The "1.6" conservative set is folded in (at most a
`dev` checkpoint tag, not a separate CRAN release).

**Tier A — correctness + the economical-RAM mode (highest value):**
- §4.1 Fix pagoda2.1 sample accessors (`getEmbedding`/`getClustering`/`getPca` keys; `$counts`
  reach-ins) — **the pagoda2 2.0 revdep prerequisite.**
- §5 / §1.1 **Disk-backed sample access**: guarantee block/streamed reads so lstar-backed samples
  lower peak RAM (the main-path memory lever).
- §1.2 **`pairs.storage = keep | disk | drop`**; §1.3 free per-pair intermediates / bound fork
  fan-out — together these deliver the **flat-peak-RAM disk-backed mode** (§5).
- §2.1 Pluggable kNN backend (RcppHNSW), threaded only where pair-level parallelism leaves cores idle.
- §12 **C++/concurrency parity (audit-driven):** fork-safe threading (the pagoda2 N2R lesson),
  threaded/streamed view-aware reducers (the `colSumByFac` disparity), built on `sccore_par.hpp`.
- §8 CRAN basics (`@return`, declare `grid`, `ComplexHeatmap`→Suggests, `\donttest` examples, unused
  Suggests, **maintainer change**).

**Tier B — leanness, consistency, integration-quality plumbing:**
- §2.2 irlba → RSpectra; §2.4/§3.2 common-case reduction reuse; §4.2 `scaledMatricesP2` via
  `getExpressionBlock`; §1.4 sparse `getJointCountMatrix` (cheap, off-path).
- §7.1–7.4 API realignment, **additive + deprecation**: string `method=`; symmetric accessors;
  argument-name aliases; `runX` preferred with old names deprecated-but-working.
- §6 Graphics consistency: shared theme + `sccore` palettes + grouping resolver + `plotMarkerDotPlot`;
  fix robustness bugs (see §11 for the sccore decision).
- §3.4 Integration-QC wrapper (mixing entropy + over-integration); §3.5 surface SNN/edge-balancing.

**Tier C — decided default flips (major-version-appropriate, NEWS-flagged):**
- §3.1 **Default `space` PCA→CPCA** (decided; `CCA` stays opt-in — too aggressive); Leiden
  `n.iterations` 2→~5–10 + a resolution-sweep helper. (`ncomps=40` left as-is — per-pair.)

### Deferred to the next conos release (post-wave, GitHub-first, with cacoa migration)

- §7.5–7.7 Remove deprecated aliases (`.Defunct`); **unify DE** (breaking output shape, couples to
  cacoa); restructure `con$graph` → `con$graphs[[method]]`; canonicalize args with removal.
- §3.3 Multimodal / WNN collection-level integration; §4.3 facet-aware accessors; §4.4 adapter contract.
- §2.3 Landmark/subsampled pairing; §3.5 scalable Laplacian solver; §3.6 joint >2-sample reductions.

---

## 10. Decisions & open questions

**Decided:**
- **Release shape** — no intermediate conos 1.6 CRAN release; **one coordinated wave** (lstar first,
  then pagoda2 2.0 + conos 2.0 together). The non-breaking set folds into conos 2.0 (at most a `dev`
  tag, not merged to main). cacoa stays off CRAN and migrates later on GitHub.
- **Default `space` = `CPCA`** (currently `PCA`). CPCA is the principled shared-subspace aligner;
  `CCA` is too aggressive (over-merges) and stays opt-in, never a default. (Behavior-affecting — land
  with a NEWS note in 2.0.)
- **`ncomps=40` stays** — it's per-pair alignment rank, not whole-atlas representation.
- **API realignment is additive in 2.0** — `runX` preferred, old names deprecated-but-working;
  removals deferred to the next major.

**Open:**
1. **Common-case reduction reuse** — per-pair recomputation stays the default (it's correct for
   differing gene sets / pair-specific CCA). Should the common-case shortcut (identical genes +
   CPCA/PCA → reuse the per-sample PCA) be **auto-detected** or an explicit opt-in flag?
2. **`pairs.storage` default** — for the economical mode, default to `"disk"` when available, or
   keep `"keep"` and let users opt into disk/drop? Same question for whether disk-backed samples
   should auto-select the economical path.
3. **Multimodal scope** — for the *next* (post-wave) release: full collection-level WNN, or start
   with "align on a chosen facet"?
4. **sccore** — keep, slim, or fold in? (See §11.)

---

## 11. Shared package: sccore

`sccore` is the shared **leaf dependency** of both pagoda2 and conos — low-level utilities plus the
plotting primitives both packages already use (`embeddingPlot`, `fac2col`/`val2col`, `dotPlot`,
`plapply`, …). Because both depend on it, sccore sits in the same "goes first" tier as lstar **if it
needs to change**; if the published sccore already suffices and we don't consolidate, leave it
untouched.

Two questions to settle during the conos work:

1. **Necessity (version pin).** Does the *published* CRAN sccore already provide everything the new
   pagoda2 + conos use, or do they rely on `dev`-only sccore (local is `dev` 1.0.7)? If the latter,
   sccore must be released **before/with lstar**, with pagoda2/conos pinning `sccore (>= <new>)`.

2. **Keep / slim / fold — judged per item, not wholesale.** The criterion: **a shared package is
   worth it only for substantial, genuinely-shared code — not for thin pass-throughs.**
   - **Not worth keeping in sccore:** thin one-liner wrappers that just call another package (e.g. a
     `knn_sparse` that only forwards to RcppHNSW, or a one-line SVD wrapper). These add a dependency
     hop and indirection for no real reuse — **inline them in each consumer** (pagoda2 already has
     `.pagoda2_knn_sparse`/`.pagoda2_truncated_svd` locally; conos can do the same). Centralizing a
     trivial wrapper isn't a benefit.
   - **Worth keeping/centralizing:** larger pieces with real shared logic and value — the common
     **plotting functions** (`embeddingPlot`, `dotPlot`, palette/`fac2col`/`val2col`, and a shared
     ggplot theme), parallel helpers, and other non-trivial utilities used by ≥2 packages. These are
     exactly what makes the §6 graphics-consistency goal achievable (one implementation, no drift),
     so the shared *graphics layer* is the strongest reason for sccore to exist.

### Inventory findings (done, 2026-06-15)

**Verdict: keep sccore — it's genuinely substantial and shared across ≥3 packages (incl. cacoa).**
The actual work is the reverse of "fold sccore in": fold conos's *drifted duplicates* back onto sccore.

**Necessity / release.** No consumer needs an *unreleased* sccore feature, but: (a) pagoda2's pin
`sccore (>= 0.1.1)` is **stale and wrong** — it actually requires **`>= 1.0.6`** (uses the
`embeddingPlot` S4 generic with `object=` from 1.0.6 and `dotPlot(scale.center=)` from 1.0.4); fix the
pin. (b) `dev` is 2 commits ahead of the `v1.0.7` tag with a **dotPlot speedup** (vectorized
`colSumByFactor` aggregation + a new exported `dotPlotData`) and a **deprecated-Matrix-coercion fix**
(CRAN-hygiene-relevant). → cut a small **sccore 1.0.8** in the leaf tier (with lstar) and pin both
consumers to it. Minor bump, low risk.

**Keep (the reason sccore exists).** The **C++ backbone** (`colSumByFactor`, `jsDist`,
`propagate_labels`, `smooth_count_matrix`, `get_nearest_neighbors`, and the `sccore_par.hpp` threading
header), the **plotting layer** (`embeddingPlot`, `dotPlot`/`dotPlotData`, `fac2col`, `val2col`), and
the **graph/DE algorithms** (`collapseCellsByType`, `getClusterGraph`, `collapseGraph*`, `multi2dend`,
graph-embedding pipeline, `propagateLabels*`, `smoothSignalOnGraph`, `mergeCountMatrices`,
`appendSpecificityMetricsToDE`, `plapply`).

**Slim (inline the one-liners, drop the export):** `sn` (`names(x)<-x`), `setMinMax` (`pmax/pmin`),
`heatFilter` (1-line kernel → fold into `smoothSignalOnGraph`), `splitVectorByNodes` (→ fold into
`graphToAdjList`).

**The real win — dedupe conos's drifted copies onto sccore.** conos carries local copies that
**shadow** sccore and have drifted apart: `extendMatrix`/`mergeCountMatrices` (TWO copies — and the
`integrations.R` one is the **dense-matrix bug of §1.4**), `propagateLabels{,Diffusion,Solver}`,
`embedGraphUmap`/`embedKnnGraph`/`graphToAdjList`/`splitVectorByNodes`. Reconcile signatures (conos's
versions grew extra args — port those into sccore), then have conos call `sccore::`. This kills the
maintenance drift **and fixes the dense-matrix RAM bug for free** (sccore's `mergeCountMatrices` is
sparse + parallel). Also: pagoda2's `papply` → `sccore::plapply` (already TODO-flagged); `sn` →
`sccore::sn`.

**Cross-package note — divergent `colSumByFactor` (don't unify).** pagoda2 *reworked*
column-sum-by-factor into its own `colSumByFacView` (`pagoda2/src/misc2.cpp`) — view-aware (applies
plain/CLR/TF-IDF inline), threaded, disk-backed — **separate** from sccore's plain in-memory
`colSumByFactor`. These are legitimately distinct and **stay separate**. sccore's `colSumByFactor`
is still load-bearing for sccore-internal (`dotPlotData`, `collapseCellsByType`,
`appendSpecificityMetricsToDE`) **and conos** (`conclass.R:943` pseudobulk); pagoda2 no longer calls
it. The plotting layer (`dotPlot`/`embeddingPlot`) *does* flow into pagoda2. **Binding rule: every
sccore change is rebuilt + tested against both pagoda2 (`devel`) and conos (`dev`) before it lands.**

**cacoa.** Several sccore exports are used by *neither* pagoda2 nor conos directly (`saveDeAsJson`,
`collapseGraph{Paga,Sum}`, `smoothSignalOnGraph`, `val2col`, `multi2dend`, standalone
`propagateLabels`) — historically cacoa's. cacoa is off-CRAN and **will be reworked**, so it no longer
gates sccore changes — we can consolidate/slim freely (verifying pagoda2 + conos per the rule above).

**Dependency notes — replace `pROC` with a built-in AUC (planned, covers both consumers).**
`appendSpecificityMetricsToDE`'s AUC branch is on pagoda2's *default* marker path
(`runMarkers(append.auc=TRUE)`) and is reachable from conos's DE too — so it's genuinely used by both,
and replacing it in sccore fixes both at once. The predictor is binary (binarized counts), so AUC has
a closed vectorized form (Mann-Whitney identity: per gene, `AUC = (a·d + 0.5·(a·c + b·d))/(n_pos·n_neg)`
with `a/c` = #expressed among in-/out-of-cluster cells) — faster than the current per-column
`apply(pROC::auc)` and it **drops `pROC` from Imports** (→ Suggests, used only by the equivalence
test). Plan: implement the vectorized AUC, add a test asserting numeric equivalence to `pROC::auc`,
and verify the AUC column is unchanged in pagoda2 *and* conos markers. Other single-justifier deps
(`irlba` smoothSignalOnGraph, `uwot` graph-UMAP, `pbmcapply` plapply) — keep. No Bioconductor deps
anywhere — all CRAN.

---

## 12. C++ / concurrency parity with pagoda2 (threading, fork-safety, view-aware reducers)

pagoda2 invested heavily in its compiled/concurrency layer during the rewrite; **conos will need the
same treatment** and currently has had little of it. Three related threads:

**(a) View-aware / disk-backed / threaded reducers — the `colSumByFac` disparity.** pagoda2 reworked
column-sum-by-factor into `colSumByFacView` (`pagoda2/src/misc2.cpp`): applies the normalization view
(plain/CLR/TF-IDF) inline, threads under a fork-safe `if(ncores>1)` OpenMP guard, and streams from the
lstar store. conos's pseudobulk (`getClusterCountMatrices` → sccore's plain in-memory `colSumByFactor`)
and its other reducers have none of this. Re-examine conos's reduce/pseudobulk/`correctGenes` paths for
the same wins: threading, streaming from disk-backed samples, and view-aware aggregation where
applicable. (Keep sccore's plain `colSumByFactor` and pagoda2's `colSumByFacView` separate — §11 — but
conos can adopt the threaded/streamed pattern, ideally via a shared primitive built on `sccore_par.hpp`.)

**(b) kNN backend (ties to §2.1).** conos's `N2R::Knn`/`crossKnn` run single-threaded (`nThreads=1`)
and N2R 1.0.5's threading is a no-op. Move to RcppHNSW (§2.1) — and make the threading fork-safe per (c).

**(c) Thread- AND fork-safety — the pagoda2 lesson.** pagoda2 hit a real bug: N2R's parallel query
shared a `visited_list_`, corrupting results under threading; the fix was per-thread searcher pools
(`BatchSearch*`) plus a fork-safe `if(nThreads>1)` serial path. conos parallelizes **across pairs** with
`mclapply` (fork), so any C++ that threads internally (OpenMP) or holds shared mutable state must be
fork-safe when run inside a forked worker. **Audit conos's `src/`** — `cpca.cpp`/`spcov.cpp`,
`Rjnmf.cpp`, `edgeweights.cpp`/`edge_rebalancing.cpp`, `propagate_labels.cpp`, `graph_embedding.cpp`,
`largeVis.cpp` — for: threaded? OpenMP? shared-mutable-state hazards under fork? Apply pagoda2's pattern
(fork-safe `if(ncores>1)` guards, per-thread scratch, no shared buffers). The shared `sccore_par.hpp`
threading header is the natural basis.

### Audit findings (done, 2026-06-15)

OpenMP **is** enabled (`src/Makevars` passes `$(SHLIB_OPENMP_CXXFLAGS)`), but there is **no fork-safe
`if(ncores>1)` guard anywhere** — the core pagoda2 idiom is absent. `sccore` is `Imports`-only, **not
`LinkingTo`**, so `sccore_par.hpp` isn't available to conos C++ — **adding `sccore` to `LinkingTo` is a
prerequisite** for the shared-primitive approach. **No *live* oversubscription bug today:** the
OpenMP-threaded files (`graph_embedding.cpp::get_nearest_neighbors`, `largeVis.cpp::sgd`,
`edgeweights.cpp::referenceWij`) all run **top-level, once per object**; the forked pair-loop
(`updatePairs` `plapply` over pairs) calls only **serial** reducers (`spcov`, `cpcaF`, `RjnmfC`). So
the hazards are latent/correctness + "must-guard-before-threading," not a current corruption.

**P0 — fork-safety / correctness (cheap, mechanical, no behavior change):**
- Add a fork-safe `if(ncores>1)` (or `omp_get_level()==0`) guard to **every** `omp_set_num_threads` /
  `#pragma omp parallel for` site (`graph_embedding.cpp:185,232,198,245,274`; `largeVis.cpp:274,278`;
  `edgeweights.cpp:110,116`; `checkfunctions.cpp::checkCRAN`). Makes `ncores==1` truly serial and is
  the precondition for threading anything that can land under a fork.
- Remove two stray `#pragma omp barrier` outside any parallel region (`graph_embedding.cpp:309`,
  `largeVis.cpp:289`); drop the unnecessary `#pragma omp critical` in `graph_embedding.cpp:209`
  (each write index is distinct); verify the per-id disjointness of `edge_weight[]` writes in
  `edgeweights.cpp::similarityOne`.

**P1 — reducers that should adopt pagoda2's view-aware / threaded / streamed pattern:**
- **`spcov.cpp` is the headline `colSumByFacView` analog.** Dense gene×gene covariance, called per
  sample **inside the forked pair-loop**, then `abind`-stacked into a dense k×p×p cube (`conos.R:185`)
  for `cpcaF` — conos's covariance densification hotspot (ties to §1 memory). Rework into a
  streamed/blocked, fork-safe covariance that feeds `cpcaF` without materializing the cube.
  **Largest, highest-value item; gated on `sccore` in `LinkingTo`.**
- `propagate_labels.cpp::smooth_count_matrix_c` — the per-edge accumulation loop is embarrassingly
  parallel and runs top-level; thread it (per-thread partials + guard) and avoid the full per-iteration
  dense copy.
- `edge_rebalancing.cpp::getSumWeightMatrix` — already a weight-sum-by-factor reducer with inline
  normalization (closest existing match); tiny output, low priority, mostly worth aligning to the
  shared primitive.

**P2 / notes:** `RjnmfC` (and any BLAS) runs inside the forked pair-loop → **BLAS thread
oversubscription**; set BLAS threads=1 in workers (R/env fix, mirroring pagoda2's
`.pagoda2_with_blas_threads`). `adjustedRand.cpp` is O(n²) serial under the stability fork (and uses raw
`malloc`/`free` — make interrupt-safe). **Fine as-is:** `gradients.cpp` (reentrant by design),
`deltacut.cpp`, `edgeFilter.cpp` (sequential greedy). The matching plan was anticipated here; the audit
confirms it.

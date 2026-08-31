# Design choices and rationale

Why the aligned package is built the way it is. Each entry states **the choice**, **why**, the
**evidence** behind it, and the **trade-off** accepted. Companion to `changes.md` (which lists *what*
changed); this file explains *why*. Numbers come from the reports in this folder.

Two goals are in tension throughout, and it is worth naming them up front:

- **Fidelity** — reproduce R `scrattch.bigcat`'s results, so the pipelines are interchangeable.
- **Correctness / efficiency** — do the statistically or computationally better thing.

Where they conflict, the choice is called out explicitly under **Deliberate divergences from R**.

---

## At a glance

| Area | Choice | Primary reason |
|---|---|---|
| Cell-level KNN | approximate (annoy) | exact is infeasible at 10⁵–10⁶ cells |
| `annoy_trees` | **50** (not the `log₂N`≈17 default) | graph stability across seeds; matches R |
| `annoy_seed` | **fixed (=1)**, decoupled from the run seed | identical graph across seeds, ~zero cost |
| Community detection | **Leiden** (`vtraag`) | largest single self-consistency lever; matches R |
| Graph weighting | **SNN Jaccard** (`jaccard_snn`) | matches R `jaccard_big`; edge overlap 0.31 → 0.90 |
| Edge pruning | `0.05`, gated at >50k cells | replicates R's `nbin>1` gating |
| Cluster-level (merge) KNN | **exact `scipy.cdist`** | faster *and* more accurate at ~10²–10³ centroids |
| Merge candidate metric | **Euclidean** (not correlation) | matches R; the dominant merge fix |
| Merge control flow | **`aligned`** (one merge/round + `<th/2` extras) | prevents premature merging; mirrors `merge_cl_big` |
| DE score | capped at **20/gene** | matches R `de_stats_pair`; without it Python over-splits |
| DE thresholds | ln-converted (`×ln2`) | R is log2(CPM+1), Python is ln(CPM+1) |
| `max_cl_size` | **disabled** (`None`) by default | matches the reference runs' `max.cl.size=Inf` |
| `score_thresh` | **`[300, 150]`** (top, recursive) | mirrors R `clust.R`'s two-level `de.score.th` |
| Merge internals | numpy row updates, GC suppressed | 74 % of merge time was pandas overhead |
| Final-merge space | selectable `markers`/`latent`/`pca` | R uses marker genes; tc had only latent/PCA |
| Port style | **additive only** (nothing removed) | non-breaking for existing tc users |

---

## 1. Graph construction

### Cell-level KNN is approximate — by necessity
**Choice:** approximate KNN on the reduced space. **Backend: pynndescent (default)**, annoy available.
**Why:** exact KNN is O(N²). At 194k cells that is 3.8×10¹⁰ distances (~280 GB); at 3M cells ~72 TB.
There is no exact option at this scale.
**Why pynndescent as the default:** measured against annoy — closer agreement with R (ARI 0.840 vs
0.829), exactly reproducible at a fixed seed, more stable across seeds (0.912 vs 0.875), and **8.3×
faster at 1M cells (43.7 s vs 361.7 s) for 1.6× the peak memory (4.6 GB vs 2.9 GB)**. Below ~400k cells
it is also the *lighter* of the two. See [9_knn_backends.md](9_knn_backends.md).
**When to use annoy instead:** byte-level comparability with `scrattch.bigcat`, or data that genuinely
cannot be indexed in memory — far beyond 1M cells, or a high-dimensional input rather than a 32-d latent
(the benchmark covers only the latter).
**Trade-off:** the graph is stochastic in principle, which is why the next two choices exist. Note
`annoy_trees` applies to the annoy backend only; NN-descent has no forest.

### `annoy_trees = 50`, not the stock `log₂N` default
**Choice:** pin to 50.
**Why:** tc's default (`None` → `max(1, int(log₂ N))`) is data-dependent and *low* — only ~17 for 194k
cells. A small forest returns seed-sensitive neighbors. R uses a fixed **50** (`get_knn_batch`), so 50
also aligns the two.
**Evidence:** first-step diff-seed self-consistency 0.51 → 0.57 directly, and — more importantly — it
stabilizes the graph enough that the fixed-seed lever becomes nearly redundant (`2_self_consistency.md`).
**Trade-off:** slightly slower build; negligible.

### `annoy_seed` fixed and decoupled from the run seed
**Choice:** `annoy_seed=1` always, independent of `random_seed`.
**Why:** guarantees a bit-identical KNN graph across runs at any seed. R gets this for free — its Annoy
(`BiocNeighbors::buildAnnoy`) uses a fixed internal seed that ignores `set.seed()`.
**Evidence:** +0.20 self-consistency under Louvain, +0.04 under Leiden. Small once trees + Leiden are in
place, but free, so kept.

### Leiden over Louvain
**Choice:** `louvain_method='vtraag'` (leidenalg).
**Why:** Leiden is robust enough that residual approximate-KNN noise doesn't reshuffle labels; it is also
what the R reference runs use.
**Evidence:** the single largest self-consistency lever — 0.57 → 0.869 at 50 trees with a varying seed.
**Trade-off:** Python `leidenalg` and R igraph `cluster_leiden` are *different implementations*, so
cell-perfect cross-pipeline identity is unattainable regardless.

### SNN Jaccard graph, not KNN-edge Jaccard
**Choice:** `weighting_method='jaccard_snn'` — Jaccard over the shared-neighbor crossproduct (`B·Bᵀ`),
so every pair sharing ≥1 neighbor gets an edge.
**Why:** this is what R `jaccard_big` builds. Stock tc weighted only existing KNN edges — a structurally
different, much sparser graph.
**Evidence:** edge Jaccard vs R **0.31 → 0.90**, weight Pearson **0.23 → 0.988**; first-step cross ARI
0.52 → ~0.88 (with the merge fix). The biggest single fidelity gain in the project.

### Edge pruning gated at >50k cells
**Choice:** drop Jaccard weights ≤ 0.05, but **only when the current (sub)graph has >50,000 cells**.
**Why:** replicates R `jaccard_big`, which prunes inside `if (nbin>1)` where `nbin = ceil(n/50000)` — so
small sub-clusters are returned unpruned. Recomputed per recursion level, exactly as R does.

---

## 2. Merging

### Cluster-level candidate KNN is **exact** (`scipy.cdist`) — a deliberate divergence from R
**Choice:** compute exact Euclidean distances between cluster centroids; no ANN index.
**Why:** merging operates on ~10²–10³ centroids, where exact is *cheaper than approximating*.
**Evidence (measured, d=32, k=4):**

| M centroids | exact `cdist` | annoy (50 trees) | annoy recall @50 | annoy recall @100 (R's setting) |
|---:|---:|---:|---:|---:|
| 450 | **3.8 ms** | 28.7 ms | 0.921 | 0.984 |
| 2,000 | **61.7 ms** | 159.5 ms | 0.752 | 0.916 |
| 3,000 | ~140 ms | — | 0.694 | 0.879 |

So at realistic cluster counts exact is both **faster and more accurate**; annoy's recall degrades as M
grows because `search_k` (= k × trees) covers a shrinking fraction of the space.
**Trade-off:** R uses annoy at *both* levels (cell `ntrees=50`, cluster `ntrees=100`), so this is a small,
intentional source of R↔Python difference — we chose correctness over bit-fidelity. It is a plausible
contributor to the "same clusters, different merge partner" cases in `5_final_merging.md`.
**Scaling caveat:** `calculate_similarity` builds the full M×M matrix every merge round — 72 MB at
M=3,000 (fine), but O(M²) would need chunking beyond ~10k clusters.

### Euclidean candidate metric, not correlation
**Choice:** `cdist_normalized` (Euclidean) for merge candidate selection.
**Why:** matches R `get_knn_pairs` (`method="Annoy.Euclidean"`). Stock tc used correlation distance for
>2 dims, which fed the merge a *different candidate set* than R.
**Evidence:** merge-on-shared-partition ARI **0.68 → 0.95**, cluster count 48 → 41 (R = 42). The dominant
merge fix.

### `merge_mode='aligned'` — one merge per round, with a strict extras gate
**Choice:** each round, re-scan **all** clusters' k-nearest, then merge the single lowest-score pair
unconditionally plus extras **only if `score < score_th/2`**; recompute means and invalidate cached DE
scores for anything touching a merged cluster.
**Why:** prevents **premature merging**. Stock tc merged *every* candidate below threshold from one
scoring pass, so a borderline pair could be committed before a prior merge revealed a better partner.
The aligned rule mirrors R `merge_cl_big`: merge the best pair, recompute, repeat.
**Evidence:** ARI vs R 0.61 → 0.68 (rules #1/#3), then → 0.95 with the Euclidean metric.
**Trade-off:** many more rounds. A `merge_mode='fast'` option is retained (batch-merge, ~0.906 ARI) for
when speed matters more than fidelity.

### The same merge function serves the per-step and final merges
**Choice:** `merge_clusters` → `merge_clusters_by_de` is used both inside each clustering step and by
`final_merge` — mirroring R, where both call `merge_cl_big`.
**Why:** one implementation to align, test and optimize.
**Evidence:** because the per-step merge was already aligned, the final merge required no new work — and
run on the same input, R and Python agree at **ARI 0.997** (`5_final_merging.md`).

### `max_cl_size` disabled by default
**Choice:** `None` → use **all** cells per cluster in the DE test.
**Why:** matches the R reference runs (`max.cl.size=Inf`). It also sidesteps R's `sample_cells`
non-determinism, which only triggers when the cap is finite.
**Evidence:** effect on results is negligible (R self-consistency 0.93 → 0.92, cross 0.52 → 0.52).
**Caveat:** R's `max.cl.size` *subsamples cells*; the Python parameter only caps the count passed to the
t-test. If a finite cap is ever needed, that is a genuine remaining difference.

---

## 3. Differential expression

### DE score capped at 20 per gene
**Choice:** `min(-log10(padj), 20)` before summing.
**Why:** R `de_stats_pair` does `tmp[tmp>20] = 20`. Stock tc summed uncapped (allowing `inf`), so Python
scores were systematically inflated and fewer pairs fell below `score_thresh`.
**Evidence:** with the cap, R `fast_limma` ≡ Python `ebayes` (score r = **1.000**, 100 % merge-decision
agreement). Without it, Python over-splits to **78–84** clusters vs R's 45–49.

### ln-converted thresholds
**Choice:** `lfc_thresh` and `low_thresh` multiplied by ln2 (e.g. `1 → 0.6931472`).
**Why:** R normalizes with **log2(CPM+1)**, Python with **ln(CPM+1)**. The t-statistic, p-values and DE
score are scale-invariant, but thresholds read directly in log units are not.

### Two-level `score_thresh = [300, 150]`
**Choice:** a 2-element list resolved to top-level vs recursive thresholds.
**Why:** mirrors R `clust.R`, which uses `de.score.th=300` for the top level and `150` for refinement —
previously a single value applied at every level in Python.

---

## 4. Performance (behavior-preserving)

Every optimization below was verified **bit-identical** to the path it replaced.

### numpy-level merge updates
**Choice:** hold cluster means/vars/present as numpy arrays with a `label→row` map inside the merge loop;
rebuild small DataFrames only for the per-round DE call.
**Why:** cProfile showed **74 % of merge time** was `numpy.array` construction driven by pandas
`.loc[row]=` / `.drop` on 17,277-column frames (dtype introspection, ~69k times) — not real computation.
**Evidence:** merge 8,918 s → **62 s** (aligned, ~144×); full pipeline 2 h 33 m → **18.9 min**, with an
end-to-end check showing the old and new paths produce **identical** partitions (ARI 1.000).
**Subtlety preserved:** the original variance formula's `(mean2 − mean_comb)²` term was silently zero (a
pandas view was overwritten before use); the numpy version reproduces that exactly rather than "fixing"
it, so results stay identical.

### `no_gc_collect` around `multipletests`
**Choice:** neutralize `gc.collect()` during statsmodels' `multipletests`.
**Why:** it calls `gc.collect()` on *every* call, firing thousands of full collections in the per-pair DE
loop. `gc.disable()` cannot stop an explicit call.
**Evidence:** ~10× on the DE-heavy path.

### numpy `de_pairs_ebayes`
**Choice:** replace the per-pair 8-column DataFrame + `filter_gene_stats` with numpy boolean masks.
**Why:** per-pair pandas construction dominated after the GC fix.
**Note:** `process_pair` (the parallel marker-selection variant) intentionally still uses the pandas path.

---

## 5. Final merge

### Selectable reduced space: `markers` / `latent` / `pca`
**Choice:** `final_merge(..., space=...)` with optional `rm_genes`; default preserves prior behavior
(`latent` if an embedding is configured, else `pca`).
**Why:** R's `merge.R` builds `rd.dat` from **raw marker-gene expression** (no PCA) — its rationale being
that final merging should be judged on the genes that define clusters. tc offered only PCA-of-markers or
the latent, so R's actual behavior was unavailable.
**Trade-off:** the two spaces select *different candidate pairs*, so they are not interchangeable; the
choice is exposed rather than hard-coded.

---

## 6. Packaging

### Additive port — nothing removed
**Choice:** the alignment changes were applied to upstream `transcriptomic_clustering` **without**
removing `cluster_louvain_phenograph`, `normalization.py`, or `pairwise_DEGs.py`, and without renaming
the package.
**Why:** keeps the change non-breaking for existing tc users and keeps the PR reviewable — it adds and
improves, never deletes. Upstream's own newer fixes (`dimension_reduction`, `filter_known_modes`,
`pairwise_DEGs`) were preserved rather than overwritten with older forked copies.

---

## Deliberate divergences from R

| Divergence | Why we chose it |
|---|---|
| **Exact `cdist` for merge candidates** (R uses annoy) | Faster *and* more accurate at 10²–10³ centroids; annoy recall falls to ~0.88 at 3k clusters even at R's 100 trees. |
| **`max_cl_size` caps a count, doesn't subsample** | Default is disabled anyway; subsampling is what makes R's `sample_cells` non-deterministic. |
| **Different Leiden implementation** (`leidenalg` vs igraph) | Unavoidable; both are Leiden with resolution 1. |

These are *known and bounded* — they are why cross-pipeline agreement plateaus at ARI ≈ 0.83 (full
recursive) rather than 1.0, alongside the irreducible different-Annoy-forest effect.

## Deferred / open choices

- **pynndescent as a cell-level backend.** Measured at **1M cells**: **43.7 s / 4.61 GB** vs annoy's
  **361.7 s / 2.87 GB** — ~**8× faster** with higher recall, at ~1.6× the memory (annoy's index is
  memory-mapped and out-of-core, so its memory grows ~1.6 GB/M cells vs pynndescent's ~4 GB/M).
  Not implemented — annoy remains the only backend. See `knn_backends.md`.
- **Chunked / exact-KNN cluster similarity** for >10k clusters, to avoid the O(M²) matrix.
- **Marker-space final merge is available but not the default** — the latent path remains default for
  backward compatibility.

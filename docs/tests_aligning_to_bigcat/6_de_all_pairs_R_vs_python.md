# `de_all_pairs`: R `scrattch.bigcat` vs Python `transcriptomic_clustering`

Head-to-head test of the exhaustive all-pairs DE step, run on the **same clustering** so only the
implementation differs.

**Verdict: the two are statistically equivalent, and the residual difference is mostly R's own
randomness.** R was tested in two configurations -- **R-200** (its default: cluster statistics from a
random 200 cells per cluster) and **R-ALL** (all cells). Across all 91,378 pairs of the 428-cluster
taxonomy:

Every measurement, across the three regimes, **grouped by which dataset it comes from**.

### From the per-pair summary (`de_summary/` — the `num` / `score` columns, all 91,378 pairs)

| measurement | what it tells you | **R-200** (as shipped) | **R-ALL** (all cells) | **R-ALL + Python `present` fix** |
|---|---|---:|---:|---:|
| **`num` Pearson** | do the two agree on *how separable* each pair is | 0.99884 | 0.999922 | **0.999999** |
| **`score` Pearson** | same, for the DE score that drives merging | 0.99862 | 0.999973 | **0.999999** |
| **`num` median abs diff** | typical per-pair disagreement, in genes | 4 genes | **0** | **0** |
| **pairs with the same `num`** | how often the gene *counts* match exactly (not the genes — see the detail block) | 17.8 % | 65.7 % | **88.7 %** |
| **mean `num`** (R / Python) | is either side systematically more permissive | 438.49 / 436.58 | 434.98 / 436.58 | 434.98 / **434.77** |
| **direction** R>Py / R<Py / equal | *which way* the residual leans — a bias, or symmetric noise | 43.7 / 38.5 / 17.8 % | 10.1 / **24.2** / 65.7 % | 11.3 / **0.0** / 88.7 % |
| **merge agreement @ `score<100`** | would the two merge the same pairs at a strict threshold | 99.998 % | 100.0000 % | **100.0000 %** |
| **merge agreement @ `score<150`** | …at R's operating threshold — the number that matters | 99.982 % | 99.9989 % | **100.0000 %** |
| **merge agreement @ `score<300`** | …at a permissive threshold | 99.967 % | 99.9956 % | **100.0000 %** |
| **R reproducible run-to-run?** | can R even reproduce itself — the precondition for any comparison | **no** (r = 0.9980) | **yes** (r = 1.00000) | yes |

### From the per-gene detail (`de_parquet/` — 3,000 sampled pairs, `up` / `down` reported separately)

| measurement | what it tells you | **R-200** | **R-ALL** | **R-ALL + fix** |
|---|---|---:|---:|---:|
| **gene-set Jaccard** | do they pick the *same genes* per pair, not just as many | 0.9385 / 0.9355 | 0.9954 / 0.9950 | **0.9997 / 0.9997** |
| **identical gene sets** | fraction of pairs whose gene lists match **exactly** | 27.8 / 26.6 % | 82.9 / 81.5 % | **96.0 / 95.7 %** |
| **top-10 overlap** | agreement on the genes anyone actually reads | 0.9174 / 0.9218 | 0.9620 / 0.9637 | **0.9625 / 0.9652** |
| **rank Spearman** | do shared genes get ranked the same way | 0.9872 / 0.9869 | 0.9939 / 0.9938 | **0.9939 / 0.9938** |
| **genes kept** (R / Python) | is one side keeping more genes per pair | 141.8 / 141.5 | 141.1 / 141.5 | **141.1 / 141.1** |

### From the marker queries (`select_top_pos_markers_ds`, `n.markers=3`, all 428 clusters)

| measurement | what it tells you | **R-200** | **R-ALL** | **R-ALL + fix** |
|---|---|---:|---:|---:|
| **identical set** | the end product — same 3 markers per cluster | 78.5 % | 94.9 % | **96.5 %** |
| **same order** | …and in the same order | 69.9 % | 92.3 % | **93.9 %** |
| **same top-1** | agreement on the single best marker | 90.0 % | 97.4 % | **97.7 %** |
| **mean shared of 3** | average overlap per cluster | 2.773 | 2.944 | **2.960** |
| **mean Jaccard** | overlap, chance-insensitive | 0.889 | 0.973 | **0.9811** |
| **shared 3 of 3** | clusters where all three markers match (order ignored) | 336 (**78.5 %**) | 406 (**94.9 %**) | **413 (96.5 %)** |
| **shared 2 of 3** | one marker differs | 87 (20.3 %) | 20 (4.7 %) | **13 (3.0 %)** |
| **shared 1 of 3** | two markers differ | 5 (1.2 %) | 2 (0.5 %) | **2 (0.5 %)** |
| **shared 0 of 3** | no marker in common — complete disagreement | **0 (0.0 %)** | **0 (0.0 %)** | **0 (0.0 %)** |

The 3-of-3 row is by definition the same as "identical set" (order-independent overlap of 3 means the
sets match). The distribution adds what the summary statistic hides: **no cluster has ever disagreed
completely**, in any regime, and the residual is one marker on a shrinking minority of clusters —
87 → 20 → 13 clusters. For reference, the two R configurations differ from *each other* by about as much
as R-200 differs from Python: **R-ALL vs R-200 = 352 / 73 / 3 / 0** (82.2 % fully shared).

Reading it left to right: **R-200 → R-ALL** is what raising `max.cl.size` buys (R stops subsampling and
becomes reproducible); **R-ALL → +fix** is what the Python `present` fix buys. The `direction` row is the
clearest single indicator — `R < Python` runs 38.5 % → 24.2 % → **0.0 %**, i.e. the bias is gone, and the
remaining 11.3 % one-directional residual is the DE engine (#3 below).

Two independent defects were found and each is now accounted for: **R subsamples 200 cells per cluster**
(fix: raise `max.cl.size`) and **Python computed detection rates via `np.mean` on a sparse mask**, which
returns one ulp above an exact half and wrongly passes the strict `q1 > 0.5` filter (fixed here). With
both addressed, the two pipelines agree on **100 % of merge decisions** and the only remaining
difference is a small one-directional offset between `fast_limma` and `ebayes`.

The sharpest result is this: **at R's default settings, two R runs of the same input disagree with each
other more than R disagrees with Python.** `get_cl_stats_big` subsamples 200 cells per cluster with a
non-reproducible seed. Raising `max.cl.size` removes nearly all the disagreement and makes R
deterministic -- so most of what looks like a cross-language difference is R disagreeing with itself.

What remains at R-ALL is a small *systematic* difference (Python keeps marginally more DE genes, by well
under one gene per pair), which is the genuine implementation gap and was previously masked by the
subsampling noise.

On runtime, Python is **~2.1x faster in parallel** measured on the same node (98 s vs 209 s on n291) and
~2.5x faster serial (cross-node, indicative only).

Two differences found here have since been **aligned to R**: the detail's `logPval` is no longer capped,
and the detail now encodes direction in `P1`/`P2` the way scrattch.bigcat does. Neither changes any
measurement in this report.

## The two R configurations (read this first)

R was run in **two configurations**, and every result below is labelled with which one it uses. The
difference is a single argument to `get_cl_stats_big`:

| label | `max.cl.size` | cells used for cluster means / present / variance | reproducible? |
|---|---|---|---|
| **R-200** | `200` (**R's default**) | a random **200 cells per cluster**, re-drawn every run | **no** |
| **R-ALL** | `1e9` (above every cluster size) | **all cells** -- `sample_cells` short-circuits and never touches the RNG | **yes** |

Python always uses **all cells** and is always reproducible, so it is a single configuration throughout.

**R-200 is what you get out of the box**, so it is the honest baseline for "how do the two pipelines
compare as shipped". **R-ALL is the like-for-like comparison**, since it feeds both pipelines the same
statistics. Where both were run, both are shown.

| result | R-200 | R-ALL |
|---|:--:|:--:|
| Summary agreement (91,378 pairs) | yes | yes |
| Merge-eligibility | yes | yes |
| R self-consistency (two runs) | yes | yes |
| Per-gene detail (3,000 pairs) | yes | yes |
| Top-3 markers (428 clusters) | yes | yes |
| Runtime, same node (n291) | -- | yes |
| Runtime, cross-node table | yes | yes |
| Serial vs parallel | yes | -- |

Datasets on disk: `summary_R_full_parallel.csv` / `de_parquet_R_full_parallel/` are **R-200**;
`summary_R_allcells.csv` / `de_parquet_R_allcells/` are **R-ALL**; `summary_py_newschema.csv` /
`de_parquet_py_newschema/` are Python.

## Implementation differences at a glance

Every point of comparison examined in this report, with the measured effect. Ordered by impact within
each group. Rows struck through were divergences that have since been **fixed on the Python side**.

### A. Affects results

| # | | R `scrattch.bigcat` | Python `transcriptomic_clustering` | measured impact |
|---|---|---|---|---|
| 1 | **Cells used for cluster moments** (means / present / variance) | **random 200 per cluster** (`get_cl_stats_big(max.cl.size=200)`), reseeded from clock/PID each call | **all cells**, deterministic | **the dominant effect at R's defaults.** Pairs with the same gene count 17.8 % (R-200) vs 65.7 % (R-ALL); top-3 markers 78.5 % vs 94.9 %; and it makes **R itself non-reproducible** — two R runs differ on 77 % of pairs |
| 2 | ~~**`present` (detection rate) arithmetic**~~ | counts in C++, **one division** → exact `k/n` | ~~`np.mean` on a sparse mask → 1 ulp above an exact half~~ → **now counts and divides once** | **RESOLVED — was the entire `R < Python` bias.** 19,308 entries sit at exactly `q1 = 0.5`; 5,971 disagreed, Python passing where R failed in **all** of them. Fixing it took merge agreement to **100.0000 %** and markers to 96.5 %, and is **1.83× faster** |
| 3 | **DE engine** | `fast_limma` | `ebayes` | **the only difference left.** Verified per gene on identical inputs: **0.000 %** of genes match even at 1e-6, Python's `-log10(padj)` uniformly **lower** (median −0.061, p-values ~15 % larger), so Python is more conservative and R keeps more genes in **11.3 %** of pairs, never fewer |
| 4 | ~~**`logPval` in the per-gene detail**~~ | **uncapped** (finite max 318.8, plus `Inf` on underflow) | ~~capped at 20.0~~ → **now uncapped, incl. `Inf`** | **RESOLVED** — detail only. Summary scores were always identical, since both sum values truncated at 20 |
| 5 | **Variance derivation** | reconstructed from `cl.sqr.means` as `E[x²] − E[x]²` (cancellation-prone; rel. error to **1e-03**, the largest of the three inputs) | `np.var(..., ddof=1)` computed directly | **none measurable.** Swapping R's variance into Python's DE shifted the mean gene count by **+0.000** and changed 0.2 % of pairs — the eBayes moderation pools across all 428 clusters and absorbs it |
| 6 | **Log normalization** | `log2(CPM+1)`, `low.th = 1` | `ln(CPM+1)`, `low_thresh = ln(2)` | **none.** Both mean "CPM > 1"; the implied **CPM ratio is exactly 1.000000** (median, mean, IQR) and `lfc_thresh = ln(2)` is exactly log2FC 1 |

### B. Checked and found aligned — no difference to find

| # | | both sides | how it was checked |
|---|---|---|---|
| 7 | **`present` values themselves** | identical to float precision | after the #2 fix, exact `k/n` on both sides; before it, 4.5e-13 apart on 49.6 % of 7,394,556 values — enough to flip the strict `q1 > 0.5` at exact halves |
| 8 | **Cluster means precision** | agree to **1.4e-07** mean relative (max 8.7e-06) | independent CPM accumulation, not a formula difference; **flips 0.00 genes per pair** across the `lfc` threshold |
| 9 | **Multiple-testing correction** | **Holm** | R `p.adjust(pval)` (its default is holm); Python `multipletests(..., method='holm')` |
| 10 | **`cl.size` (the *n* in the test)** | **full** cluster sizes | R `table(cl)`, Python `{c: len(v)}`. Worth stating: R uses full *n* even when its moments came from 200 cells |
| 11 | **Pair enumeration and ordering** | **lexicographic** on the label string; bins by **numeric** rank | `create_pairs` / `tc.create_pairs`; `'19' < '3'`, so `cl.bin` ordering and pair ordering are deliberately different bases on both sides |
| 12 | **`top_n`, DE thresholds, summary fields** | `top_n = 500`; `q1/q2/qdiff/padj/lfc/min_genes/score_thresh` converted exactly; `num`, `up_num`, `down_num`, `score`, `up_score`, `down_score` | the input clustering is the same file for both runs |
| 13 | **Score cap when summing** | truncate each gene at **20** into `up_score` / `down_score` | R `get_de_truncate_score_sum`, Python `SCORE_CAP`. This is why the summary was identical even while the detail differed (#4) |
| 14 | **Detail content volume** | row counts within **0.3 %** (R-200 27,637,938 / R-ALL 27,484,657 / Py 27,565,123) | the 791 MB vs 371 MB on-disk gap is **parquet encoding, not extra data** |

### C. Differs, but not in the numbers

| # | | R | Python | note |
|---|---|---|---|---|
| 15 | ~~**Detail schema / partition layout**~~ | direction in **`P1`/`P2`** — each pair stored under both orientations, so `bin.x`/`bin.y` fill the full 5×5 grid = **25** dirs | ~~direction in `sign`, `bin.x ≤ bin.y` = 15 dirs~~ → **now matches R** | **RESOLVED.** Verified on 7.3 M R rows that `pair` and `sign` are both derivable from `P1`/`P2` (0 violations). Filtering on `sign` is a trap — it is relative to lexicographic order (`'17' < '3'`) |
| 16 | **`cl.bin` default sizing** | `cl.bin.size = min(100, n_clusters/mc.cores)` — **layout depends on core count** | `bin_size = 100` fixed, **decoupled from `n_jobs`** | deliberate: keeps a dataset's on-disk layout reproducible and queryable regardless of how many workers wrote it |
| 17 | **Worker count rule** | `min(mc.cores, ceiling(n_pairs/5000))` — self-caps by workload | `min(n_jobs, n_bin_pair_tasks)` — no workload cap | both landed on **15 workers** here, which is what makes the runtime comparison fair. It is also why the 60-cluster subset was misleading: R capped itself to 1 core |
| 18 | **Parallel backend** | `mclapply` **fork** — frames inherited copy-on-write, free | `Pool` + initializer — frames **pickled** to each worker (~59 MB each) | why R scales **5.9×** and Python **3.2×** on the same 15 workers |
| 19 | **Serial vs parallel equivalence** | **not identical** — 70,321 of 91,378 pairs differ | **bitwise identical** on every field | R's mismatch is not a parallelism bug: each job called `get_cl_stats_big` separately and drew a different subsample (#1). At R-ALL two runs agree at r = 1.00000 |
| 20 | **Runtime** (`de_all_pairs`, 91,378 pairs) | 209 s parallel / 1,545 s serial | **98 s** parallel / 616 s serial | **~2.1× faster in parallel on the same node** (both n291) and ~2.5× serial. Node variation here is ~2×, so cross-node numbers cannot support a speed claim — see the runtime section |
| 21 | **`gc.collect()` in multiple-testing correction** | n/a | neutralized via `no_gc_collect` | statsmodels calls it on every `multipletests`; it was **92 % of runtime**. Pure speed (25×), outputs bitwise identical |

## Changes made to the Python implementation

**Three performance fixes**, all verified to leave results bit-identical:

1. **`no_gc_collect` around `multipletests`** -- statsmodels calls `gc.collect()` on every call, which
   profiling showed was **92 % of runtime** (14.0 s of 15.2 s). This alone gave **23.6x**.
2. **`_de_one_pair` rewritten with numpy masks** instead of a per-pair DataFrame + `filter_gene_stats`.
   Worth ~1.0x on its own -- the pandas overhead was negligible next to the gc problem -- but it keeps
   the QC path consistent with the merge path.
3. **`n_jobs` extended to the in-memory path**, so parallelism no longer requires writing to disk.

Before these, the same 1,770-pair subset took 366 s; after, 14.6 s.

**Two alignments to R**, made after the comparisons below surfaced them. Both change only how the detail
is stored, never which genes are selected or how they are ranked:

4. **Detail `logPval` uncapped** -- R stores the raw `-log10(padj)` and truncates only the copy it sums
   into the score. Python now does the same, including R's `Inf` for `padj` underflow.
   ([details](#logpval-capping-found-here-since-aligned-to-r))
5. **Direction encoded in `P1`/`P2`** -- `P1` is now always the cluster the gene is up in, as
   scrattch.bigcat does, so a pair's genes span both orientations of its bin-pair.
   ([details](#the-detail-schema-direction-lives-in-p1p2-python-now-aligned))

**One bug fix**, which *does* change results -- and fixes a genuine defect rather than matching R for its
own sake:

6. **`present` computed as an exact rational.** `get_cluster_means_inmemory` counted detections with
   `np.mean` on a sparse mask, which returns one ulp above an exact half and wrongly passes the strict
   `q1 > 0.5` filter. It now counts and divides once. This removed the entire `R < Python` bias, took
   merge agreement to 100 %, and is **1.83x faster**.
   ([details](#the-fix-and-what-it-bought))

Every full-scale comparison was re-run against output written *after* the `logPval` and `P1`/`P2` changes (results unchanged, as expected), and again after the
`present` fix (results improved -- see
[The fix](#the-fix-and-what-it-bought)).

## Setup

| | value |
|---|---|
| Input clustering | `final_merge_compare/Py_final_merged_latent.csv` — the final-merged partition, **identical for both** |
| Clusters tested | **all 428** → **91,378 pairs**, 194,221 cells (a 60-cluster / 1,770-pair subset was also run; same conclusions) |
| DE params (R) | `de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=10)` |
| DE params (Python) | same, ln-converted: `lfc_thresh=0.6931472` (= log2 1), `padj_thresh=0.01`, `cluster_size_thresh=10` |
| R stats | `get_cl_stats_big(big.dat, cl, stats=c("means","present","sqr_means"))` |
| Python stats | `tc.get_cluster_means(adata, ...)` → means / present / variances |
| `top_n` | 500 both sides |

Both write the same per-pair summary: `up_num`, `down_num`, `num`, `up_score`, `down_score`, `score`.

## Agreement (full scale: 428 clusters, 91,378 pairs)

**Scope: everything in this section comes from the per-pair SUMMARY** -- R's `de_summary/` (read with
`open_dataset(summary.dir)`) versus Python's `de_summary_py/` (read with `read_de_pairs`), each dumped
to `summary_{R,py}_full_parallel.csv` and joined on the unordered cluster pair. The compared fields are
the six summary statistics `num`, `up_num`, `down_num`, `score`, `up_score`, `down_score`
(R spells them `up.num` / `up.score`, normalized before joining).

The per-gene detail (`de_parquet/`) is compared separately in "Per-gene detail" below.

All 91,378 pairs matched on both sides (0 R-only, 0 Python-only).

**R-200 vs Python** (R as shipped):

| field | R-200 mean | Python mean | Pearson | Spearman | median abs diff | within 10 % |
|---|---:|---:|---:|---:|---:|---:|
| `num` | 438.49 | 436.58 | **0.99884** | 0.99804 | 4 genes | 94.7 % |
| `up_num` | 203.68 | 202.87 | 0.99904 | 0.99819 | 2 genes | 91.1 % |
| `down_num` | 234.81 | 233.71 | 0.99895 | 0.99861 | 2 genes | 91.2 % |
| `score` | 6368.79 | 6320.41 | **0.99862** | 0.99831 | 50.8 | 95.6 % |
| `up_score` | 2989.95 | 2967.87 | 0.99879 | 0.99827 | 21.17 | 92.2 % |
| `down_score` | 3378.84 | 3352.53 | 0.99873 | 0.99864 | 21.65 | 92.4 % |

**R-ALL vs Python** (same statistics on both sides):

| field | R-ALL mean | Python mean | Pearson | Spearman | median abs diff | within 10 % |
|---|---:|---:|---:|---:|---:|---:|
| `num` | 434.98 | 436.58 | **0.99992** | 0.99988 | **0** | **99.9 %** |
| `up_num` | 202.17 | 202.87 | 0.99995 | 0.99987 | **0** | 99.3 % |
| `down_num` | 232.81 | 233.71 | 0.99992 | 0.99989 | **0** | 99.3 % |
| `score` | 6309.12 | 6320.41 | **0.99997** | 0.99994 | 4.55 | **100.0 %** |
| `up_score` | 2962.85 | 2967.87 | 0.99998 | 0.99992 | 0.74 | 99.7 % |
| `down_score` | 3346.27 | 3352.53 | 0.99997 | 0.99992 | 0.91 | 99.6 % |

The median absolute difference in gene count drops from **4 genes to 0** purely by removing R's
subsample.

**Merge-eligibility agreement.** This is a per-pair boolean taken straight from the summary's `score`
column -- `score < threshold` on each side, compared pair by pair over all 91,378 pairs. It is *not* a
run of the iterative merge algorithm (which merges the lowest-scoring pair, recomputes cluster stats and
repeats, so a single early disagreement can cascade). Read it as "how often do the two agree that a pair
is a merge candidate", not "the two produce the same merged taxonomy":

| threshold | R-200 | R-ALL | Python | agreement, **R-200** | agreement, **R-ALL** |
|---|---:|---:|---:|---:|---:|
| `score < 100` | 5 | 5 | 5 | 99.998 % | **100.000 %** |
| `score < 150` | 41 | 38 | 37 | 99.982 % | **99.999 %** |
| `score < 300` | 233 | 247 | 243 | 99.967 % | **99.996 %** |

At R's operating threshold of 150, the two disagree on **16 of 91,378 pairs (0.018 %)**.

## Where the difference comes from: R-200 subsamples 200 cells per cluster

| | R > Python | R < Python | equal |
|---|---:|---:|---:|
| `num` | 43.7 % | 38.5 % | 17.8 % |

Neither implementation looks consistently more permissive at R's default settings -- the difference
appears to be symmetric noise. (It is: with the subsampling removed the symmetry disappears and a small
systematic bias emerges -- see [The decisive test](#the-decisive-test-r-with-all-cells).) Running R
twice identifies the source:

```r
get_cl_stats_big <- function(big.dat, cl, max.cl.size = 200, stats = c("means"), ...) {
    sampled.cells = sample_cells(cl, max.cl.size)   # <- random 200 cells per cluster
    cl = cl[sampled.cells]
```

`get_cl_stats_big` does not use every cell. It draws **200 cells per cluster** through `sample_cells`,
which reseeds from the clock/PID (the `set.seed(NULL)` behaviour documented in
[2_self_consistency.md](2_self_consistency.md)). Clusters here average 454 cells, so **R computes its
cluster means from roughly 44 % of the data, and from a different 44 % on every run.** Python uses all
cells.

The consequence, measured on the same 91,378 pairs:

| comparison | `num` Pearson | median abs diff | same `num` |
|---|---:|---:|---:|
| **R-200 run 1 vs R-200 run 2** | 0.9980 | 4 genes | 23.0 % |
| **R-200 vs Python** | **0.9988** | 4 genes | 17.8 % |
| **R-ALL run 1 vs R-ALL run 2** | **1.00000** | **0** | **100 %** |
| **Python run 1 vs Python run 2** | **1.0000** | 0 | **100 %** |

**At R's default settings, R agrees with Python better than it agrees with itself.** That is exactly what subsampling predicts:
an R-vs-R comparison contains two independent 200-cell draws, an R-vs-Python comparison contains only
one (R's) plus Python's exact full-data computation. Deviation from perfect correlation is `1 - r =
0.0020` for R-vs-R against `0.0012` for R-vs-Python -- a ratio of 1.7, close to the factor of 2 expected
if independent sampling noise dominates.

So the residual disagreement is **not** primarily an implementation gap. The implementation-level effect
that does exist is **R's variance reconstruction** (`sqr_mean - mean^2` vs Python's direct
`np.var(ddof=1)`), measured and quantified in
[What remains](#what-remains-is-small-and-traced-to-its-source). The `log2` vs `ln` difference
contributes nothing -- both sides' thresholds are exact conversions of each other.

This is consistent with the earlier finding that the DE engines themselves are equivalent (R
`fast_limma` = Python `ebayes`, score r = 1.000 on identical inputs).

**Practical note.** If you need a reproducible R run, compute `cl.stats` once and reuse the saved object
(as the test harness here does), or pass `max.cl.size` large enough to disable subsampling. Otherwise
re-running post-clustering QC will shift DE scores by a few genes per pair, which at a threshold of 150
is enough to change a handful of merge decisions.

## The decisive test: R-ALL (R using all cells)

This section is the roll-up of the **R-ALL** column shown throughout the report. Setting `max.cl.size`
above every cluster size makes `get_cl_stats_big` use **all cells** -- `sample_cells` returns the full
set when `sample.size >= length(cells)`, so it never reaches the RNG either. That feeds both pipelines
the same statistics and isolates the true implementation gap.

```bash
sbatch _run_de_all_pairs_R.sh 0 20 _allcells 1000000000   # 4th arg = max.cl.size
```

**1. R becomes exactly reproducible.** Two independent all-cells runs:

| | `num` | `score` |
|---|---|---|
| **R-ALL** run 1 vs run 2 | r = **1.00000**, 100 % identical | r = **1.00000**, 100 % identical |
| **R-200** run 1 vs run 2 (for contrast) | r = 0.9980, 23.0 % identical | r = 0.9974 |

So the subsampling was the *only* source of R's non-determinism. The serial-vs-parallel mismatch
reported earlier disappears once `max.cl.size` is raised.

**2. Agreement with Python improves sharply at every level.**

| | **R-200** vs Python | **R-ALL** vs Python |
|---|---:|---:|
| `num` Pearson | 0.99884 | **0.99992** |
| `num` median abs diff | 4 genes | **0 genes** |
| pairs with the same `num` | 17.8 % | **65.7 %** |
| `score` Pearson | 0.99862 | **0.99997** |
| `score` median abs diff | 50.8 | **4.55** |
| detail gene-set Jaccard (up/down) | 0.9385 / 0.9355 | **0.9954 / 0.9950** (median **1.0000**) |
| detail identical gene sets | 27.8 % / 26.6 % | **82.9 % / 81.5 %** |
| detail top-10 overlap | 0.9174 / 0.9218 | **0.9620 / 0.9637** |
| detail rank Spearman | 0.9872 / 0.9869 | **0.9939 / 0.9938** |
| **top-3 markers: identical set** | 78.5 % | **94.9 %** |
| top-3 markers: same order | 69.9 % | **92.3 %** |
| top-3 markers: same top-1 | 90.0 % | **97.4 %** |
| top-3 markers: mean Jaccard | 0.8890 | **0.9729** |

Merge-eligibility agreement becomes near-perfect: **100.000 %** at `score < 100`, **99.999 %** at 150,
**99.996 %** at 300.

For context, dropping the subsample moves R about as far as the R-Python gap itself
(**R-ALL vs R-200**: `num` r = 0.99891, median diff 3 genes, markers Jaccard 0.9091) -- i.e. **most of
what looked like an R-vs-Python difference was R disagreeing with itself.**

### What remains, and where it actually came from

With the subsampling gone, a real bias surfaced that had been masked:

| | R > Python | R < Python | equal |
|---|---:|---:|---:|
| `num`, **R-200** | 43.7 % | 38.5 % | 17.8 % |
| `num`, **R-ALL** | 10.1 % | **24.2 %** | 65.7 % |

Python kept **+1.60 genes per pair** more than R (436.58 vs 434.98). Tracking that down took four
measurements, and the first three hypotheses were all wrong -- recorded here because the wrong turns are
the instructive part.

**Step 1 -- compare the three DE inputs element-wise** (428 x 17,277 = 7,394,556 cluster-gene values):

| input | R vs Python | first read |
|---|---|---|
| `present` | max diff **4.5e-13** | "identical" |
| `means` | mean rel. 1.4e-07, max 8.7e-06 | negligible |
| `variance` | median rel. 3.3e-07, 99.99th pct **4.5e-04** | **largest error -> prime suspect** |

**Step 2 -- inputs or the DE test?** Feeding R's own statistics through Python's DE
(`_diag_swap.py`, all 428 clusters):

| arm | mean `num` |
|---|---:|
| **A** Python stats -> Python DE | 436.58 |
| **B** R's stats -> Python DE | 434.77 |
| **R** R's stats -> R's DE | 434.98 |

| isolates | Pearson | pairs identical | direction |
|---|---:|---:|---|
| **B vs A** -- the **inputs** | 0.99992 | 73.5 % | B<A 25.1 %, B>A 1.4 % |
| **B vs R** -- the **DE test** | 1.00000 | 88.7 % | B<R 11.3 %, B>R 0.0 % |

**Step 3 -- which input?** Swapping them one at a time (`_diag_which_input.py`):

| arm | mean `num` | shift vs all-Python | pairs differing |
|---|---:|---:|---:|
| R **variance** only | 436.58 | **+0.000** | **0.2 %** |
| R **means + present** only | 434.77 | **-1.809** | **26.3 %** |
| all R stats | 434.77 | -1.809 | 26.5 % |

**The variance -- the largest numerical difference -- contributed nothing.** It is smoothed away by the
eBayes moderation pooled across all 428 clusters. And the means turned out to flip **0.00 genes per
pair** across the `lfc` threshold, with the R/Python CPM ratio exactly **1.000000** (median, mean and
IQR), so the normalisation is identical and the 1e-07 differences are pure summation rounding.

That left `present` -- the input whose difference was *smallest*.

**Step 4 -- the mechanism.** `present` is a detection rate, `k / n` cells above `low_th`, so it is an
**exact rational** -- and for a 454-cell cluster, `227/454` is **exactly 0.5**. The `q1` filter is a
**strict** `q1 > q1_thresh` with `q1_thresh = 0.5`. The two pipelines computed that rational differently:

```r
# R (Rcpp_parallel_transpose_dense.cpp): count into an integer, then ONE division
res(col_i, cluster) += 1.0;
col = col / cluster_id_number[i];        // 227.0 / 454.0 = 0.50000000000000000000
```

```python
# Python (get_cluster_means_inmemory), before the fix:
np.mean((sliced_X > low_th), axis=0)     # -> 0.50000000000000144329
```

`np.mean` on a **sparse** boolean mask dispatches to scipy's sparse mean, which is a matrix-vector
product against a vector of `1/n` -- it adds `1/n` to itself `k` times. Reproduced exactly:

```
sum of 227 copies of 1/454 by repeated addition = 0.50000000000000144329   > 0.5 -> True
227 / 454 by a single division                  = 0.50000000000000000000   > 0.5 -> False
```

One ulp, and a half-detected gene passes the strict filter in Python and fails it in R. Measured over
the full matrix:

| | |
|---|---|
| `present` entries differing at all | **3,669,074 of 7,394,556 (49.6 %)** -- the earlier "0 differ" used a 1e-12 tolerance that hid a 4.5e-13 error |
| entries sitting at **exactly 0.5** | 19,308 |
| entries where `q1 > 0.5` **disagrees** | **5,971 -- every one at an exact tie** |
| direction | **Python passes / R fails: 5,971. The reverse: 0.** |
| reach | 85 of 428 clusters, ~70 genes each, x 427 pairs per cluster |

**The smallest numerical difference of the three was the only one that mattered**, because it is the only
quantity that is hard-thresholded at a value its discrete values land on *exactly*. Ranking the
candidates by error magnitude -- which is what steps 1-2 did -- pointed at the wrong one every time.

### The fix, and what it bought

`get_cluster_means_inmemory` now counts and divides once, matching both R's C++ and this repo's own
`get_cluster_means_backed` (which had always done it correctly -- the two Python paths had drifted):

```python
present_count = np.asarray((sliced_X > low_th).sum(axis=0)).ravel()
present_cluster_means_lst.append(present_count / n_cells_in_cluster)
```

It is also **1.83x faster** (10.82 ms -> 5.92 ms per cluster; 4.63 s -> 2.54 s over 428 clusters),
because scipy's sparse mean does a full matvec where a reduction suffices. There was no speed/accuracy
trade-off -- the old form was slower *and* wrong.

Re-running the full 428-cluster comparison against **R-ALL** (the full three-regime table is at the
[top of this report](#de_all_pairs-r-scrattchbigcat-vs-python-transcriptomic_clustering)):

| R-ALL vs Python | before fix | **after fix** |
|---|---:|---:|
| `num` Pearson | 0.999922 | **0.999999** |
| `score` Pearson | 0.999973 | **0.999999** |
| pairs with the same `num` | 65.7 % | **88.7 %** |
| mean `num` (R 434.98) | 436.58 | **434.77** |
| direction **R < Python** | 24.2 % | **0.0 %** |
| direction R > Python | 10.1 % | 11.3 % |
| merge agreement @ `score<100` | 100.0000 % | **100.0000 %** |
| merge agreement @ `score<150` | 99.9989 % | **100.0000 %** |
| merge agreement @ `score<300` | 99.9956 % | **100.0000 %** |
| **top-3 markers identical** | 94.9 % | **96.5 %** |
| top-3 markers, same top-1 | 97.4 % | **97.7 %** |
| top-3 markers, mean Jaccard | 0.9729 | **0.9811** |
| **detail gene-set Jaccard** (up/down) | 0.9954 / 0.9950 | **0.9997 / 0.9997** |
| detail identical gene sets | 82.9 / 81.5 % | **96.0 / 95.7 %** |
| detail genes kept (R / Python) | 141.1 / 141.5 | **141.1 / 141.1** |

**The `R < Python` direction disappears entirely, merge decisions agree on all 91,378 pairs at every
threshold, and the per-gene detail becomes near-exact** — Jaccard 0.9997 with the median pair
gene-for-gene identical, and genes-kept now matching to the first decimal (141.1 vs 141.1, from
141.1 vs 141.5). 73.7 % of pairs were unchanged by the fix; the rest shifted by the expected
+1.809 genes.

**What is left is exactly the DE engines, and nothing else.** After the fix, R vs Python shows 88.7 %
identical with R>Python 11.3 % and R<Python 0.0 % -- **the same numbers the swap experiment measured for
the `fast_limma` vs `ebayes` arm alone** (88.7 % identical, B<R 11.3 %, B>R 0.0 %).

That agreement of summary statistics is suggestive but not proof, because a log2->ln conversion artifact
would produce the same signature. So the engines were compared **per gene on identical inputs**
(`_diag_engine.py`: Python's DE run on R's own statistics *with detail*, versus R's stored detail; 300
pairs, 90,805 genes present in both):

| per-gene `logPval`, identical inputs | |
|---|---|
| identical to 1e-6 | **0.000 %** |
| identical to 1e-9 | 0.000 % |
| **Python > R** | **0.00 % of genes** |
| median difference | **-0.061** (mean -0.097, max 0.920) |
| correlation | 0.9999999528 |
| genes only Python kept / only R kept | 2 / 26 |

**The engines really do differ, and not at the margins.** *Every* gene's p-value differs, always in the
same direction -- Python's `-log10(padj)` is uniformly lower, i.e. its p-values are ~15 % larger. A
threshold or conversion artifact would flip genes *at* a boundary and leave the rest untouched; this
shifts the entire distribution. The correlation of 0.99999995 shows it is a smooth systematic offset,
not noise.

So Python's `ebayes` is slightly **more conservative** than R's `fast_limma`, which is exactly why R
keeps more genes in 11.3 % of pairs and never fewer. The likely origin is the variance moderation --
`df_prior` / `var_prior` estimation, or R's pooled `get_cl_sigma` versus Python's per-gene moderated
variance. Neither is "wrong"; they are different eBayes formulations. Chasing it further is only
worthwhile if exact per-gene p-value parity is required.

> **Scope warning.** `get_cluster_means` also feeds `merge_clusters` and `final_merge`, so this fix
> changes clustering output, not only post-clustering QC. Before the fix Python admitted ~5,971 marginal
> genes it should not have, making clusters look marginally more separable and therefore **under-calling
> merges**. The effect on the merge path has **not** been measured here.

### Practical recommendation

Two changes, one per pipeline, make R and Python agree on **100 % of merge decisions**:

1. **R:** raise `max.cl.size` in `get_cl_stats_big` (or compute `cl.stats` once and reuse it). Costs
   almost nothing -- `cl.stats` went 64 s -> 69 s, since that step is I/O-bound on the parquet backend
   rather than limited by cell count -- and it buys exact reproducibility plus a jump from 78.5 % to
   94.9 % marker agreement.
2. **Python:** the `present` fix is already applied here. It is strictly a correctness improvement
   (exact rational arithmetic, no accumulation), it is 1.83x faster, and it takes marker agreement to
   96.5 % and merge agreement to 100 %.

Neither is a matter of matching the other pipeline for its own sake: R's subsampling makes its own
output irreproducible, and Python's `np.mean` was returning a detection rate that is provably wrong at
exact halves.

## Runtime (full scale, and the effect of parallelism)

*(`de_all_pairs` runs on the `428 x 17,277` stat matrices, so its cost is the same for R-200 and R-ALL;
only the preceding `get_cl_stats_big` step differs -- 64 s vs 69 s.)*

| run | cores requested | workers actually used | `de_all_pairs` | node | speedup from parallelism |
|---|---:|---:|---:|---|---:|
| **R-200 parallel** | `mc.cores=20` | 15 | **260 s** | n294 | **5.9x** |
| **R-200 serial** | `mc.cores=1` | 1 | **1,545 s** | n297 | -- |
| **R-ALL parallel** | `mc.cores=20` | 15 | **209 s** | n291 | -- |
| **Python parallel** | `n_jobs=16` | 15 | **194 s** | n286 | **3.2x** |
| **Python parallel** (re-run) | `n_jobs=16` | 15 | **98 s** | n291 | -- |
| **Python serial** | `n_jobs=1` | 1 | **616 s** | n294 | -- |

*(R additionally spends 64 s in `get_cl_stats_big`; the Python side spends its setup time loading the
h5ad and computing cluster means. Only `de_all_pairs` is compared above.)*

**The four runs above were on four different nodes, and node variation here is ~2x** -- the identical
Python parallel job took 194 s on n286 and 98 s on n291, same partition, same 16 CPUs. That is larger
than the R-vs-Python gap, so the cross-node table cannot support a speed claim on its own.

**Same-node comparison.** Two full-scale parallel runs landed on **n291**, which gives a fair number.
`de_all_pairs` operates on the `428 x 17,277` stat matrices, so its cost does not depend on
`max.cl.size` and the all-cells R run is directly comparable:

| run (both n291, 15 workers) | `de_all_pairs` |
|---|---:|
| **R-ALL** (job 25235624) | **209 s** |
| **Python** (job 25235552) | **98 s** |

**Python is ~2.1x faster than R in parallel on the same node** -- larger than the earlier cross-node
estimate of 1.34x, and this time measured against a controlled comparison. The serial numbers below
remain cross-node and should be treated as indicative only.

**R scales better across workers than Python does** (5.9x vs 3.2x on the same 15 workers), which narrows
the gap in parallel. The reason is how each ships data to workers: R's `mclapply` **forks**, so the
cluster means/present frames are inherited copy-on-write at zero cost, while Python's `Pool` **pickles**
them to each of the 15 workers through the initializer (~59 MB per frame). Python pays a fixed startup
cost that R does not. Switching the Python pool to a fork start method with module-level globals instead
of `initargs` would recover part of that -- an open optimization, not done here.

### How the worker count is actually decided

Requesting N cores does not mean N workers -- both sides reduce it:

- **R** caps by workload first: `mc.cores = min(mc.cores, ceiling(n_pairs/5000))`. At 91,378 pairs that
  is 19; the work is then split into bin-pair tasks, giving `min(19, 15) = 15`.
- **Python** has no workload cap; the on-disk path uses `min(n_jobs, n_bin_pair_tasks) = min(16, 15) = 15`.

Both therefore ran **15 workers** here, which is what makes the runtime comparison fair.

> This is also why the earlier 60-cluster subset was misleading: it has only **1,770 pairs** and, at
> `bin_size=100`, only **1 bin -> 1 task**. R capped itself to `ceiling(1770/5000) = 1` core and Python
> had a single task, so **both ran serial** despite being asked for 20/16 cores.

### Why 3.2x and not 15x

The bin-pair tasks are uneven: a diagonal bin-pair holds `C(100,2) = 4,950` pairs, an off-diagonal one
`100 x 100 = 10,000`, and the last bin is a 28-cluster stub (smallest task 378 pairs). Makespan is set
by the largest task, capping speedup at `91,378 / 10,000 ~ 9.1x`. The rest of the gap is fixed cost:
each of the 15 workers receives the cluster means/present frames (428 x 17,277 floats, ~59 MB each)
through the pool initializer. Smaller `bin_size` would even out the load at the cost of more partitions.

### Correctness of the parallelization

*(This subsection uses **R-200**, R's default.)*

**Python serial and parallel outputs are bitwise identical across all 91,378 pairs** (every field of the
summary). Parallelism changes only scheduling, never results.

**R-200's serial and parallel outputs are not identical** -- 70,321 of 91,378 pairs differ in `num` (max
614 genes) and 91,061 differ in `score` (max 11,838). This is *not* a bug in R's parallelism: the two
jobs each called `get_cl_stats_big` independently and drew different 200-cell subsamples, so they were
computing DE on different inputs. This is confirmed by R-ALL: with `max.cl.size` raised, two independent runs are bit-for-bit identical
(r = 1.00000), so the serial/parallel mismatch disappears entirely. It
does mean you cannot use an R re-run as a reference for "did my change alter the output" unless the
`cl.stats` object is held fixed.

## Per-gene detail: do they pick the same genes?

Everything above tests the **summary**. This section tests `de_parquet/` -- whether the two pipelines
select the *same genes* and rank them the same way.

Sampled **200 pairs from each of the 12 populated partitions = 3,000 pairs**, comparing R's detail rows
against Python's for the same pair, `up` and `down` directions separately:

| metric | **R-200** up | **R-200** down | **R-ALL** up | **R-ALL** down |
|---|---:|---:|---:|---:|
| gene-set Jaccard, mean | 0.9385 | 0.9355 | **0.9954** | **0.9950** |
| gene-set Jaccard, median | 0.9576 | 0.9581 | **1.0000** | **1.0000** |
| gene sets exactly identical | 27.8 % | 26.6 % | **82.9 %** | **81.5 %** |
| top-10 marker overlap | 0.9174 | 0.9218 | **0.9620** | **0.9637** |
| rank correlation (Spearman) | 0.9872 | 0.9869 | **0.9939** | **0.9938** |
| genes kept, R mean | 141.8 | 170.9 | 141.1 | 170.1 |
| genes kept, Python mean | 141.5 | 170.6 | 141.5 | 170.6 |

**With R as shipped (R-200) the gene lists agree at ~94 %; feeding both sides the same statistics
(R-ALL) takes that to 99.5 %, with the median pair gene-for-gene identical.** Shared genes are ranked
nearly identically either way (rho 0.987 -> 0.994), and the top-10 markers overlap 92 % -> 96 %. Genes-kept counts
match to within 0.3 genes on average, so neither side is systematically keeping more.

The ~6 % that differ are the marginal genes sitting on the filter thresholds, and by the section above
they are largely R's subsampling noise rather than an implementation gap. Note that the *identical-set*
rate is low (27 %) while Jaccard is high (0.94): with ~150 genes per pair, a single gene differing
breaks exact equality, so the exact-match rate is a much harsher metric than it looks.

**Size on disk is not a content difference.** The detail datasets differ (R 791 MB vs Python 371 MB) but
the row counts do not: **R-200 27,637,938 / R-ALL 27,484,657 vs Python 27,565,123 rows**, all within
0.3 %. The gap is parquet
encoding/compression, not extra data.

### `logPval` capping (found here, since aligned to R)

The comparison runs surfaced one genuine column difference: R's detail stored **uncapped** `logPval`
(values of 20.87, 20.41 appear) while Python **capped at 20.0**.

R keeps the raw per-gene value and truncates only the copy it sums into the score:

```r
up.genes = setNames(-log10(df[up,"padj"]), up)   # stored: UNCAPPED
tmp = up.genes; tmp[tmp > 20] = 20               # scored: capped at 20
up.score <- sum(tmp)
```

Python has been changed to match: `_logpval` no longer truncates, and `SCORE_CAP` is applied only when
summing `up_score` / `down_score`. The contract is now the same on both sides:

```
summary.score == sum(min(detail.logPval, 20))     # over the pair's genes
```

**Nothing measured in this report changes.** The score was already computed from capped values on both
sides, so every summary statistic, merge-decision count, and runtime above is unaffected; gene ranking
is by `padj`, so the gene-level Jaccard and rank correlations are unaffected too. Only the value stored
in the detail column changed, and only for genes whose `padj` is below 1e-20.

Regression tests: `test_logpval_is_not_capped` (raw value survives: 1e-25 -> 25.0) and
`test_score_is_sum_of_capped_detail` (the score/detail contract, including the up/down split) in
`tests/test_post_clustering_qc.py`.

## The detail schema: direction lives in `P1`/`P2` (Python now aligned)

This is the one structural difference that matters when reading the per-gene detail, and it explains
the partition counts:

| | partitions | detail on disk | rows |
|---|---:|---:|---:|
| R | **25** (= 5 x 5, the full bin grid) | 791 MB | 27,637,938 |
| Python | **15** (= 5 x 6 / 2, upper triangle) | 371 MB | 27,565,123 |

**R writes each pair under both orientations, one per direction.** For the pair `1` / `10`:

| stored as | rows | `sign` | meaning |
|---|---:|---|---|
| `P1=1, P2=10` | 289 | all `up` | the genes higher in cluster **1** |
| `P1=10, P2=1` | 24 | all `down` | the genes higher in cluster **10** |

Both records carry `pair = "1_10"`. So in R, **`P1` is always the cluster the gene is up in**, and `sign`
is redundant with the `P1`/`P2` orientation. Python stores the pair once with `P1`/`P2` fixed in
lexicographic order and puts the direction in `sign` alone (228 `up` + 27 `down` in one record).

Two consequences:

1. **Partition count.** Because R materializes both orientations, `bin.x`/`bin.y` cover the full
   `5 x 5` grid; Python's normalized `bin.x <= bin.y` covers the `15`-cell upper triangle. Same content,
   different tree -- and **not cross-readable**: each pipeline's reader only finds its own layout.
2. **Row counts stay comparable.** R splits a pair's genes across two records (289 + 24) where Python
   keeps them in one (228 + 27), so the totals are within 0.26 % rather than 2x. The 791 MB vs 371 MB
   gap is parquet encoding/compression, not extra data.

Verified on 7,265,907 R detail rows: **`sign` and `pair` are both fully derivable from `P1`/`P2`**
(`sign == 'up'` exactly when `P1 < P2` lexicographically; `pair` is the two labels in lexicographic
order), with zero violations. They are conveniences, not information.

### Python has been changed to match

Python now uses R's convention: the detail stores `P1` = the cluster the gene is up in, splitting each
pair's genes across the two orientations, with `pair` canonical and `sign` retained (redundant) for
readability. Reading "the genes up in cluster `g`" is now the same one-liner on both sides:

```python
rows where P1 == g          # both pipelines
```

Consequences, all verified by test:
- Python's **detail** now fills the full bin grid (3 bins -> 9 partitions in the unit test, 5 -> 25 at
  full scale), matching R. The **summary** is unchanged: one row per unordered pair, upper triangle only.
- `read_de_pairs` now opens both orientations of a bin-pair and returns exactly the rows the in-memory
  path produces; total row count is unchanged (nothing duplicated or dropped).
- `_dir_genes` (used by `check_triplet`) now selects on `P1` rather than `sign`.
- Tests: `test_detail_encodes_direction_in_p1_p2` plus the updated `test_de_all_pairs_binned_roundtrip`.

One practical warning this change removes: filtering the detail on `sign` is a trap, because `sign` is
relative to lexicographic order. For the pair `('3','17')`, `'17' < '3'`, so `sign == 'up'` means
"higher in **17**", not in 3. `P1 == '3'` is unambiguous.

## Final marker lists: the comparison that matters most

`top_n = 500` genes are stored per pair, but downstream work does not use them directly -- it reduces
them to a handful of markers per cluster. That reduction is scrattch.bigcat's `select_top_pos_markers_ds`:
for each cluster `g`, take the genes up in `g` against every other cluster and score each gene by

```
score(gene) = sum over partners of (1000 - rank)        # R's get_gene_score_ds, max.num = 1000
```

then keep the top `n.markers`. Because `1000 - rank` is between 500 and 999 for every stored gene, the
score is dominated by **how many partners a gene is a top marker against**, with rank as a tiebreak.

Run at `n.markers = 3` over all 428 clusters, on each pipeline's own DE output:

| | result |
|---|---|
| **Validation** -- the same Python scorer run on **R's** DE output | **428/428 clusters identical, same genes in the same order** (Jaccard 1.000) |
| **R vs Python** -- each on its own DE output, for both **R-200** and **R-ALL** | see below |

The validation line matters: it shows the Python port of the scoring is exact, so everything below is a
difference in the **DE data**, not in the marker algorithm.

| metric (428 clusters, 3 markers each) | **R-200** vs Python | **R-ALL** vs Python | **R-ALL** vs Python **+ `present` fix** |
|---|---:|---:|---:|
| identical list (same genes, same order) | 299/428 (69.9 %) | 92.3 % | **93.9 %** |
| identical set (order ignored) | 336/428 (78.5 %) | 94.9 % | **96.5 %** |
| same top-1 marker | 385/428 (90.0 %) | 97.4 % | **97.7 %** |
| mean markers shared, of 3 | 2.773 | 2.944 | **2.960** |
| mean Jaccard | 0.889 | 0.973 | **0.9811** |
| clusters sharing 3 / 2 / 1 / 0 | 336 / 87 / 5 / **0** | -- / -- / -- / **0** | 343 / 80 / 5 / **0** |

(The `present` fix also lifts the **R-200** comparison, from 78.5 % to 80.1 % identical sets -- but R's
subsampling remains the limiting factor there.)

For reference, the two R configurations differ from *each other* by about as much as R-200 differs from
Python: **R-ALL vs R-200 = 82.2 % identical sets, Jaccard 0.909.**

**No cluster disagrees completely** in either configuration. With R as shipped the strongest marker
agrees 90 % of the time, rising to 97.4 % once R uses all cells. Disagreement concentrates in the
*third* marker -- these R-200 examples are typical:

```
cl  2:  R [Lama2, Ndnf, Itga4]        Py [Lama2, Ndnf, Rspo2]
cl  4:  R [Slc6a5, Slc32a1, Gad1]     Py [Slc6a5, Slc32a1, Pax2]
cl 13:  R [Srd5a2, Sox6, Glis3]       Py [Srd5a2, Sox6, Aox1]
```

This is the sharpest illustration of what the 200-cell subsampling costs. Agreement degrades along the
chain as each step makes a harder discrete choice from nearly-tied candidates:

| level | agreement |
|---|---:|
| pair scores (continuous) | r = **0.999** |
| gene sets per pair (~150 genes) | Jaccard **0.94** |
| top-3 markers per cluster | Jaccard **0.89**, 78.5 % identical |

None of this is an implementation gap -- the scorer reproduces R exactly on R's own data. It is R
re-drawing 200 cells per cluster, propagated through to a top-3 cut where a few points of score
difference reorders near-ties. Holding `cl.stats` fixed on the R side would remove most of it.

Scripts: `_markers_R.R` (calls `select_top_pos_markers_ds` directly), `_markers_py.py` (the port, which
runs on either convention), `_compare_markers.py`. Jobs 25235192 (R, 32 s) and 25235201 (Python, 16 s).

## Reproduce

```bash
cd ../de_all_paris        # test_clustering/de_all_paris
# <n_clusters> <cores> <output tag>;  n_clusters=0 means all 428
# R-200 (default max.cl.size=200)
sbatch _run_de_all_pairs_R.sh   0 20 _full_parallel
sbatch _run_de_all_pairs_R.sh   0 1  _full_serial
sbatch _run_de_all_pairs_py.sh  0 16 _full_parallel
sbatch _run_de_all_pairs_py.sh  0 1  _full_serial
bash   _collect_runtimes.sh               # -> runtimes_full.csv
python _compare_full.py                   # -> comparison_metrics_full.csv  (summary agreement)
sbatch _compare_detail.sh                 # -> detail_comparison.csv        (per-gene agreement)

# R-ALL: all cells instead of the default 200-cell subsample (4th arg = max.cl.size)
sbatch _run_de_all_pairs_R.sh 0 20 _allcells  1000000000
sbatch _run_de_all_pairs_R.sh 0 20 _allcells2 1000000000   # replicate, to test determinism
sbatch _compare_allcells.sh               # summary + detail + markers, all three comparisons

# Python re-run after the `present` fix, and its comparison
sbatch _run_de_all_pairs_py.sh 0 16 _fixed
sbatch _compare_fixed.sh                  # R-ALL vs Python, before/after + markers

# where does the residual come from? (needs R's cl.stats exported first)
sbatch _export_clstats_R.sh               # -> Rclstats_{means,present,sqr_means}.parquet
sbatch _diag_residual.sh                  # element-wise means / present / variance comparison
sbatch _diag_swap.sh                      # R's stats through Python's DE: inputs vs test
sbatch _diag_which_input.sh               # swap means/present/variance one at a time
sbatch _diag_means.sh                     # is the normalisation different? (no: CPM ratio 1.000000)
sbatch _diag_ties.sh                      # exact ties at q1 = 0.5, and which side they fall
sbatch _diag_engine.sh                    # per-gene logPval: fast_limma vs ebayes on identical inputs

# final marker lists (n.markers = 3)
sbatch _markers_R.sh  _full_parallel 30 3 # R-200 -> markers_R.csv + marker_genes_universe.csv
sbatch _markers_py.sh                     # -> markers_py.csv + markers_pyimpl_on_Rdata.csv
python _compare_markers.py
```

Both sides read the **same** input clustering (`final_merge_compare/Py_final_merged_latent.csv`) and
build `cl.bin` / `cl_bin` explicitly with `bin_size = 100`, as `post_clustering_qc.R` does.
Jobs: **R-200** 25233660 (parallel) / 25233659 (serial); Python 25233662 (parallel) / 25233661 (serial),
re-run post-schema-change as 25235552; detail comparison 25234160; **R-ALL** 25235624 / 25235663 with
comparison 25235664; residual diagnostics 25235942 (stats export), 25235952 / 25235995
(input comparison), 25236012 (input-vs-test swap), 25236200 (per-input swap), 25236256 (normalisation),
25236292 (exact ties), 25236733 (per-gene engine comparison); Python re-run with the fix
25236405 / 25236406.

The whole-pipeline consequences of the `present` fix are in
[6_full_recursive_comparision_after_fixing_present.md](6_full_recursive_comparision_after_fixing_present.md).

`_compare_detail.py` walks one bin-pair at a time and reads only the matching `(i,j)` and `(j,i)`
partitions from each side, so peak memory stays bounded despite the ~28 M rows/side.
Usage example: `transcriptomic_clustering/examples/de_all_pairs_example.py`.

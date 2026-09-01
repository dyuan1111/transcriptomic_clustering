# Post-clustering QC: replicating the R doublet / low-quality workflow

Can `transcriptomic_clustering` reproduce the reference R workflow
(`clustering_bigcat/post_clustering_qc.R`, the "label doublet and low quality clusters" section)?

**Verdict: yes, now — after fixing three fidelity gaps.** Verified against R on the same dataset:
doublets-by-triplet and low-quality clusters come out **identical**, and the one apparent difference in
marker-based doublets turned out to be **a bug in R**, which detects the same cluster but loses its name
on the way out.

| | R | Python |
|---|---|---|
| doublets by triplet | 286, 300, 312 | **286, 300, 312** |
| low-quality candidates | 154 | **154** (same list) |
| doublets by marker | *(empty — R bug)* | **393** (correct) |

Of the three gaps found on the Python side, one had broken doublet detection completely (it reported
zero doublets on a dataset with three) and one made the workflow unrunnable at scale.

## The R workflow, and what it calls

| R step | R function | Python |
|---|---|---|
| per-cluster QC stats from `qc.csv` | `dplyr` group_by/summarize | plain pandas (in the example) |
| enumerate all cluster pairs | `create_pairs` | `tc.create_pairs` |
| bin map | `cl.bin` data.frame | `tc.make_cl_bin` |
| candidate doublet triplets | `find_triplets_big` | `find_triplets` |
| score the triplets | `find_doublets_all_big` -> `check_triplet_big` | `find_doublets` -> `check_triplet` |
| doublets by marker | `find_doublet_by_marker` | `find_doublet_by_marker` |
| low-quality clusters | `find_low_quality_big` | `find_low_quality` |

## Gaps found

### 1. `check_triplet` capped gene sets R deliberately leaves uncapped

R applies `rank <= top.n` to the four *scored* gene sets but **not** to the two *overlap* sets:

```r
up.genes  = de.df %>% filter(P1==cl1 & P2==cl2 & rank <= top.n)   # capped
tmp.genes = de.df %>% filter(P1==cl  & P2==cl2)                    # NOT capped
olap.up.genes1 = intersect(tmp.genes, names(up.genes))
```

Python capped all six. The overlap sets are what `olap_ratio_up_1` / `olap_ratio_down_1` are computed
from — and those are exactly the quantities the doublet call thresholds on — so capping them shrinks the
overlap and **understates the ratio**, making genuine doublets harder to detect.

Fixed: `_dir_genes` takes `top_n=None` for "no rank cap", and `check_triplet` uses it for those two sets.

**Measured impact — this one broke doublet detection entirely.** Scoring all 3,670 triplets both ways:

| | capped (old) | uncapped (R-faithful) | changed on |
|---|---:|---:|---:|
| `score` | 0.2784 | **0.5139** | 99.2 % of triplets |
| `olap_ratio_up_1` | 0.2290 | **0.4903** | 96.5 % |
| `olap_ratio_down_1` | 0.2290 | **0.4903** | 96.5 % |
| **doublet clusters called** | **none** | **286, 300, 312** | — |

Capping the overlap sets roughly **halved every score**. Since the workflow keeps clusters with
`score > 0.8` and `olap sum > 1.4`, nothing cleared the bar: the old code reported **zero doublets on a
dataset that has three**. It would have failed silently — an empty doublet list looks like a clean
taxonomy, not like a bug.

### 2. `find_doublets` returned one row per candidate; R keeps every tested triplet

This looks harmless until you notice an asymmetry in the R workflow:

| | threshold |
|---|---|
| `find_doublets_all_big` **early-stops** at | `score > 0.8` **and** olap sum > **1.6** |
| the workflow then **selects** at | `score > 0.8` **and** olap sum > **1.4** |

A triplet scoring between 1.4 and 1.6 never triggers the early stop, but *does* pass selection — and R
has it, because it writes every tested triplet to disk. Returning only the best-scoring row per
candidate silently drops those, so clusters R calls doublets can be missed.

Fixed: `find_doublets` returns every tested triplet; the early-stop semantics are unchanged.

### 3. No way to run it without the whole detail in memory

R reads each triplet's rows from the partitions of its three clusters. Python required the entire detail
as a DataFrame — 27.5 M rows for 428 clusters, and worse on a larger taxonomy.

Fixed: `find_doublets(None, triplets, root=..., cl_bin=...)` streams per triplet, as R does. The
in-memory path is retained and is much faster when the detail fits.

## What already matched

- `find_low_quality` ≡ `find_low_quality_big` (same filter, same `cl` / `cl_low` orientation).
- `find_doublet_by_marker` ≡ R's, and it drops markers absent from the matrix automatically — R does
  that by hand (its `rm_` step) and errors otherwise.
- `check_triplet`'s score, `sum(olap.score) / sum(up + down + up2 + down2)`, is identical.
- `find_triplets` was aligned in the previous round: R's `arrange(cl.up, down.num.x + down.num.y)`
  ordering, the `pair1` / `pair2` columns, and the `all_pairs` / `select_cl` restrictions.

## End-to-end run

All four steps run on the full 428-cluster dataset (`summary_R_allcells` + `de_parquet_R_allcells`),
job 25237379, **1 m 36 s total**:

| step | result |
|---|---|
| `find_triplets` | 3,670 rows -> **19 candidate clusters** |
| detail load (148 involved clusters) | 2,266,061 rows, 10 s |
| `find_doublets` | 3,655 triplets scored, 14 s |
| **doublet clusters** (score > 0.8, olap sum > 1.4) | **286, 300, 312** |
| `find_low_quality` | 437 rows, 154 distinct low-quality candidates |
| `find_doublet_by_marker` | 1 cluster flagged (all 11 markers present) |

Top-scoring triplets:

```
 cl cl1 cl2    score  olap_up_1  olap_down_1
300 229 331    0.880       1.00         0.76
312 157 331    0.830       0.98         0.68
286 128 331    0.825       0.86         0.82
```

The marker hit is cluster **393**, expressing `Slc32a1` (3.69, neuronal) together with `Opalin` (4.81)
and `Sox10` (2.80, oligodendrocyte) — two lineages above the 3.5 threshold in one cluster, which is what
the check is for.

**Impact of each gap, measured (`_gap_impact.py`):**

| gap | changes results? | evidence |
|---|---|---|
| **1. capped overlap sets** | **yes, decisively** | scores halved; **0 doublets called instead of 3** |
| 2. best-row-only | not on this dataset | 0 triplets fall in the 1.4-1.6 band, so both give 286/300/312 |
| 3. in-memory only | no — a scalability limit, not a numerical one | streamed and in-memory paths agree |

Gap 2 is worth fixing anyway: R's selection threshold (1.4) really is looser than its early-stop
threshold (1.6), so a dataset with a triplet between them would lose a doublet call. But this run is
**not** evidence that it mattered — here that band is empty.

Note also that only 3,655 of 3,670 triplets were scored — the early stop fires for a few candidates, so
their remaining triplets are never tested. That matches R.

## Verified against R on the same dataset

The three nominated cluster lists were compared directly — R's workflow (`_qc_workflow_R.R`,
job 25237494) and Python's, both reading `de_summary_R_allcells` / `de_parquet_R_allcells` with the
same `cl.bin`, marker list and thresholds:

| step | R | Python | identical |
|---|---|---|---|
| `find_triplets` | 3,670 rows / 19 candidates | 3,670 rows / 19 candidates | **yes** |
| **doublets by triplet** | **286, 300, 312** | **286, 300, 312** | **yes** |
| low-quality (`low.df` rows / distinct `cl.low`) | 437 / **154** | 437 / **154** | **yes** (same cluster list) |
| doublets by marker | *(empty)* | **393** | **no — see below** |

### The marker difference is an R bug, not a disagreement

Both implementations identify the same cluster. Running R's internals directly:

```
clusters with >1 marker group above 3.5: 1
which: 393
    Neuron Olig Astro Endo VLMC Peri SMC Micro
393   3.69 4.81  0.53 0.22    0    0   0  0.22
```

R flags 393 exactly as Python does — `Slc32a1` at 3.69 (neuronal) together with `Opalin` at 4.81
(oligodendrocyte). The loss happens in the **return value**:

```r
return(t(cl.means[unlist(markers), doublet]))
```

When exactly **one** cluster is flagged, `cl.means[genes, doublet]` drops to a vector (R's `drop=TRUE`
default), so `t()` produces a `1 x n_genes` matrix whose `rownames` are `NULL` — the cluster name is
gone, and the colnames are the gene names instead. The workflow's
`cl.doublets_byMarker <- rownames(doublet.means)` therefore returns **nothing**:

```
find_doublet_by_marker returned: dim 1x11
rownames(out): NULL
colnames(out): Slc32a1,Slc17a7,Slc17a6,Sox10,Opalin,...
```

The bug only bites when exactly one cluster is flagged — which is the common case for marker-based
doublets — and it fails silently: an empty list reads as "no marker doublets found". Python's
`find_doublet_by_marker` returns `cluster_means.loc[is_doublet, all_markers]`, which keeps the index
regardless of how many rows survive, so it is correct here. **The deviation is deliberate and not
reproduced.**

### Runtime

| | |
|---|---|
| R `find_doublets_all_big` | **1,943 s** (32 min) |
| Python `find_doublets` (in-memory) | **14 s**, after a 10 s detail load |

Both score the same 3,655 triplets. The difference is I/O strategy: R re-reads the parquet for every
triplet, Python loads the 2.27 M relevant rows once. It shows most on cluster 140, which has 1,442
triplets and never early-stops, so every one is tested. Python's `root=` path exists for when the
detail does not fit in memory and matches R's profile.

## The example

`examples/4_post_clustering_qc_example.py` replicates the workflow block by block, each annotated with the
R lines it mirrors, and writes the same three CSVs (`doubletsClusters`, `doubletsClusters_byMarker`,
`lowQClusters`).

Thresholds are collected at the bottom of the example with their R provenance, because several are
workflow-level rather than function defaults:

| step | value | from |
|---|---|---|
| `find_triplets` | `min_up_num=30, max_down_num=10, min_de_num=50` | R `find_triplets_big` defaults |
| `find_doublets` | `top_n=50, score_th=0.8, olap_th=1.6` | R `find_doublets_all_big` defaults |
| doublet selection | `score > 0.8 & olap sum > 1.4` | R workflow line 84 — **not** 1.6 |
| `find_doublet_by_marker` | `th=3.5` | R default |
| `find_low_quality` | `low_th=2` | R default |
| low-quality selection | `umi_y*1.5 < umi_x & gene_y < 4000`, then `cl_size > 1.5*cl_low_size` | R workflow lines 120 / 129 |

## Known quirks kept on purpose

**Triplets are emitted twice, mirrored.** R's self-join only filters `cl.down.x != cl.down.y` and never
deduplicates the unordered parent pair, so each triplet appears as both `(x, y)` and `(y, x)`.
`check_triplet` returns the same answer either way, so it doubles the work but not the result. Left
matching R rather than silently diverging.

**R's rename typo is not reproduced.** `find_triplets_big` has
`tmp.cols = c("up.num.new.x", "up.num.new.x", "down.num.new.x", "down.num.new.y")` — `up.num.new.x`
twice, `up.num.new.y` missing — so that one column keeps its `.new` name in R's output. The Python
version names all four consistently.

## Reproduce

```bash
cd ../de_all_paris
sbatch _qc_workflow_check.sh     # Python: all four steps on the 428-cluster dataset
sbatch _qc_workflow_R.sh         # R: the same workflow on the same dataset (~33 min)
python _compare_qc_workflow.py   # compares the three nominated cluster lists
python _gap_impact.py            # quantifies what each gap was worth
sbatch _marker_check_R.sh        # isolates the R drop=TRUE bug in find_doublet_by_marker
```

The DE dataset these steps consume is characterised in
[6_de_all_pairs_R_vs_python.md](6_de_all_pairs_R_vs_python.md); the marker-query layer over the same
dataset is in [7_marker_queries_R_vs_python.md](7_marker_queries_R_vs_python.md).

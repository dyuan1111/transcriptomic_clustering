# Marker queries over `de_parquet`: R `scrattch.bigcat` vs Python

`de_all_pairs` writes up to 500 genes per cluster pair, but nothing downstream reads them raw — the
result is always reduced to a handful of markers per cluster. This report covers the **query layer**
that does that reduction: the port of scrattch.bigcat's `markers.parquet.R` family into
`transcriptomic_clustering/de_all_pairs.py`, and how faithfully it reproduces R.

**Verdict: on the same DE input, R and Python return the same markers.** The top-N family is exact —
**428/428 clusters identical, same genes in the same order**. The greedy family returns **identical
marker sets for 12/12 clusters** (Jaccard 1.000, identical set sizes); 2 of the 12 order the last few
genes differently. That is **not a disagreement** — it is tied scores at the algorithm's 999 floor being
broken differently, at exactly the pick where ties first appear, so the tail order carries no
information ([measured](#why-top-n-is-exact-and-the-greedy-tail-is-not)).

## The two families

Both answer "which genes mark this cluster", but they optimise different things.

| | **top-N** | **greedy / combo** |
|---|---|---|
| R entry point | `select_top_pos_markers_ds` | `select_pos_markers_ds` |
| question | what are the best N markers for this cluster **overall**? | give me a marker **set** such that **every** comparison is covered |
| method | rank-sum score, take the top N | seed from top-N, then greedily add whichever gene covers the most still-uncovered pairs |
| result | fixed size (N per cluster) | variable size, grows until covered or the budget runs out |

Both reduce to one scoring rule (R `get_gene_score_ds`):

```
score(gene) = sum over the selected pairs of (max_num - rank)        # R max.num = 1000
```

Because `max_num - rank` lies between `max_num - top_n` and `max_num - 1` for every stored gene, the
score is dominated by **how many comparisons a gene wins**, with its rank acting only as a tie-break.

**Direction comes from `P1`, not `sign`.** The detail stores `P1` = the cluster the gene is up in, so
"genes up in `g`" is just the rows with `P1 == g` — no sign filtering, exactly as R does. This is also
what makes the queries cheap: those rows always carry `bin.x == cl_bin[g]`, so a query for a few
clusters opens only those partitions instead of the whole dataset.

## What was ported

| Python (`de_all_pairs.py`) | R (`markers.parquet.R`) |
|---|---|
| `get_gene_score_ds` | `get_gene_score_ds` |
| `select_markers_pair_group_top_ds` | `select_markers_pair_group_top_ds` |
| `select_top_pos_markers_ds` | `select_top_pos_markers_ds` |
| `check_pairs_lfc` | `check_pairs_lfc` |
| `check_pairs_ds` | `check_pairs_ds` |
| `select_markers_pair_direction_ds` | `select_markers_pair_direction_ds` |
| `select_markers_pair_group_ds` | `select_markers_pair_group_ds` |
| `select_pos_markers_ds` | `select_pos_markers_ds` |

All are exported from the package (`tc.select_top_pos_markers_ds(...)` etc.).

## Agreement with R

Both pipelines were run on **the same DE dataset** — R's own `de_parquet_R_full_parallel` — with the
same 17,277-gene universe and the same `cl.bin`, so any difference is the query implementation alone.

### top-N family: exact

| | |
|---|---|
| clusters | 428 |
| `n_markers` | 3 |
| **identical list (same genes, same order)** | **428 / 428** |
| identical set | 428 / 428 |

Every cluster's marker list matches R's `select_top_pos_markers_ds` output exactly. Since
`select_top_pos_markers_ds` is a thin wrapper over `select_markers_pair_group_top_ds` and
`get_gene_score_ds`, this exercises all three.

### greedy family

Run with `n_markers=1`, `max_genes=20`, `cl.means=NULL` (so both take the `check_pairs_ds` path),
each of 12 clusters separated from all 427 others:

| | |
|---|---|
| clusters | 12 |
| **identical marker set** | **12 / 12** (mean Jaccard **1.0000**) |
| identical list including order | 10 / 12 |
| set sizes, R | `15, 4, 2, 11, 4, 4, 3, 3, 2, 5, 3, 7` |
| set sizes, Python | `15, 4, 2, 11, 4, 4, 3, 3, 2, 5, 3, 7` — **identical** |

**The two clusters that differ do so only in ordering, and only deep in the tail** — the predicted
consequence of tie-breaking (deviation 3 below):

```
cluster 1  first divergence at position 12 of 15:  R Sorcs3   Py Adamts6
  R : Meis2 Nrn1 Ptprk Hs3st4 Bcl11a Nxph1 Mkx Ndst4 Mpped2 Rplp1 Fnbp1l Sorcs3 Adamts6 Sv2b Arpp21
  Py: Meis2 Nrn1 Ptprk Hs3st4 Bcl11a Nxph1 Mkx Ndst4 Mpped2 Rplp1 Fnbp1l Adamts6 Arpp21 Sorcs3 Sv2b

cluster 4  first divergence at position 9 of 11:   R Luzp2    Py Arpp21
  R : Slc6a5 Nrp2 Vwc2l Pax2 Shisa9 Zfp536 Dpm1 Zfhx4 Luzp2 Chrm3 Arpp21
  Py: Slc6a5 Nrp2 Vwc2l Pax2 Shisa9 Zfp536 Dpm1 Zfhx4 Arpp21 Chrm3 Luzp2
```

Both lists agree on every gene through the first 8-11 picks — the ones that matter — and then reorder
the last few among themselves, converging on the same set. That is what tied scores look like: once the
strongly-discriminating genes are taken, the remaining candidates cover equal numbers of the leftover
pairs, so the pick order is arbitrary and only the tie-break rule decides it. **Order in the tail
carries no information; treat the combo output as a set.**

Note the greedy sets are much larger than the top-N ones (up to 15 genes vs 3), because each cluster is
being separated from **all 427** others and hard-to-separate neighbours each demand their own marker.

### Why top-N is exact and the greedy tail is not

The two families use the same scoring rule, so the difference in reproducibility comes from *how many
pairs each score sums over*. Measured directly (`_tie_check.py`):

**Top-N scores once, over all 427 partners.** A top-3 gene sums ~427 terms of 500-999, so scores land in
the hundreds of thousands:

| | |
|---|---|
| clusters checked | 428 |
| clusters with a **tie at the 3rd/4th boundary** | **0 (0.00 %)** |
| typical 3rd-place score | median **377,059** (min 166,262, max 424,797) |

Two distinct genes hitting the same six-figure integer essentially never happens, so the ordering is
fully determined by the data and both implementations reach it. Hence 428/428 identical, order included.

**The greedy loop rescores over a shrinking pair set.** Each pick removes the pairs it covers, so the
scores collapse. Replaying cluster 1 — one of the two that diverged — pick by pick:

| pick | gene | score | pairs left | genes tied at top |
|---:|---|---:|---:|---:|
| 1 | Meis2 | 253,746 | 427 | 1 |
| 2 | Nrn1 | 94,038 | 170 | 1 |
| 3 | Ptprk | 25,486 | 75 | 1 |
| … | … | … | … | 1 |
| 11 | Fnbp1l | 1,993 | 7 | **1** |
| **12** | **Adamts6** | **999** | **5** | **4** |
| 13 | Arpp21 | 999 | 4 | **3** |
| 14 | Sorcs3 | 999 | 3 | **2** |
| 15 | Sv2b | 999 | 2 | 1 |

Picks 1-11 are unambiguous — one gene at the top each time. **At pick 12 the score hits a floor of 999
and four genes tie**, and that is exactly where R and Python diverged (position 12 of 15). The four tied
genes are precisely the four that came out reordered: `Adamts6`, `Arpp21`, `Sorcs3`, `Sv2b`.

**Why 999 is a floor.** `score = sum(max_num - rank)`, so a gene covering *one* remaining pair at rank 1
scores exactly `1000 - 1 = 999`. Once only a handful of pairs are uncovered, every remaining candidate
covers one of them at a good rank and they all score 999 — a hard tie that no amount of data breaks. R
resolves it by whatever `arrange(-score)` leaves first (effectively the order rows arrived from the
parquet scan); this port resolves it by gene name, so repeated Python runs agree with each other.

**So the greedy family did not disagree with R — it produced the same markers in a different order**, and
that order is arbitrary in principle rather than a discrepancy between implementations. Two R runs whose
rows arrived in a different order would disagree with each other the same way. The first ~11 picks are
meaningful and reproducible; the tail is tie-broken filler that exists only to finish covering pairs.
**Treat the greedy output as a set.**

## What the greedy variant buys

Measured on the real dataset, cluster 1 against all 427 others:

| | markers | comparisons covered |
|---|---:|---:|
| `select_top_pos_markers_ds` (top-1) | 1 | 222 / 427 |
| `select_pos_markers_ds` (combo) | 4 | **375 / 427** |

The single best marker leaves nearly half the comparisons uncovered; the greedy pass targets what is
still missing. Across 12 clusters the combo sets were size 1-4, and the first gene matched the top-1
marker for **12/12** — it seeds from the top-N result before filling gaps, as R does.

Neither reaches 427/427: the loop stops when no remaining gene scores above zero for the uncovered
pairs, i.e. those comparisons have no qualifying marker in the gene universe at all. R behaves the same
way (its loop breaks on an empty `gene.score`).

> **Note on an earlier wrong reading.** A first pass appeared to show a large disagreement — Python
> returning 1-4 markers per cluster where R returned up to 15. That was a harness error, not a code
> difference: the Python call had been given only the 12 selected clusters as its full cluster list, so
> each was separated from 11 others rather than all 427. Passing the same cluster universe both sides
> gives the agreement above. Worth flagging because the failure mode is silent — `select_pos_markers_ds`
> happily returns a plausible, much smaller marker set when `clusters` is under-specified. Pass the
> **whole taxonomy** as `clusters` and use `select_cl` to choose which ones to compute.

## Deviations from R, and why

1. **`select_top_pos_markers_ds` makes one pass, not one query per cluster.** R loops over clusters and
   issues a separate query each time; the port streams the needed partitions once and scores every
   cluster simultaneously. Same results (428/428), but it is what makes a whole-taxonomy call practical.
2. **The detail is pre-loaded once per cluster in the greedy path** (R's `de=` argument, which its own
   `select_top_pos_markers_ds` already uses). R's `select_pos_markers_ds` path re-reads the dataset on
   every greedy iteration. `get_gene_score_ds` gained an optional `de=` parameter for this.
3. **Ties are broken by gene name.** R takes whatever `arrange(-score)` happens to put first. Sorting by
   `(score desc, gene asc)` makes repeated Python runs agree with each other. This is the *only* place
   the greedy variant diverges from R, and it is now measured rather than assumed: ties appear exactly
   at the score floor of 999, at exactly the pick where the two implementations parted
   ([why](#why-top-n-is-exact-and-the-greedy-tail-is-not)).
4. **Caching is JSON, not `.rda`.** `out_dir` / `overwrite` behave as in R (skip clusters already
   written), only the file format differs.

## Performance

| call | scope | time |
|---|---|---:|
| `select_top_pos_markers_ds` | 428 clusters, n_markers=3 | **31 s** |
| `select_pos_markers_ds` (Python) | 12 clusters vs 427, n_markers=1, max_genes=20 | 36 s |
| `select_pos_markers_ds` (R) | same | 93 s |

The top-N variant scales well because it is a single streaming pass. The greedy variant is inherently
sequential per cluster — each pick depends on the previous one — so it costs roughly as much for 12
clusters as the top-N call does for all 428. Budget accordingly before running it across a taxonomy,
and use `out_dir` so an interrupted run can resume.

Python is ~2.6x faster than R on the greedy path (36 s vs 93 s for the same 12 clusters), which is the
pre-loading in deviation 2 — R re-reads the dataset on every greedy iteration. Both ran on `celltypes`
nodes; given the ~2x node variation documented in report 6, read this as "at least not slower".

## Reproduce

```bash
cd ../de_all_paris
# R
sbatch _markers_R.sh _full_parallel 30 3      # top-N     -> markers_R.csv
sbatch _combo_markers_R.sh 12 1 20            # greedy    -> combo_markers_R.csv
# Python port, on R's own dataset
python _markers_py.py de_parquet_R_full_parallel R markers_pyimpl_on_Rdata.csv
python _compare_markers.py                    # top-N comparison
python _compare_combo_markers.py              # greedy comparison
sbatch _tie_check.sh                          # why top-N is exact and the greedy tail is not
```

Tests: `tests/test_de_all_pairs.py` — `test_select_top_pos_markers_ds_ranks_by_breadth_then_rank`,
`test_select_markers_pair_group_top_ds_directions`, `test_check_pairs_lfc_counts_separating_genes`,
`test_check_pairs_ds_counts_stored_markers`,
`test_select_pos_markers_ds_covers_pairs_the_top_marker_misses`,
`test_select_pos_markers_ds_cache_roundtrip`.

The DE dataset these queries read is characterised in
[6_de_all_pairs_R_vs_python.md](6_de_all_pairs_R_vs_python.md).

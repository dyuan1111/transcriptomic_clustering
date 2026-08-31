# Clustering comparison reports — Python `transcriptomic_clustering` vs R `scrattch.bigcat`

Findings from aligning the Python single-cell clustering pipeline with the R reference, on the same
194,221 spinal-cord neurons (17,277 genes, 32-d scVI latent) — and, in report 10, on the 285,230-cell
multi-platform WMB2 `TH-EPI-Glut` neighborhood. Numbered files follow the order the work was done;
images live in `images/`.

> **Published copy.** These reports are working documents copied from the Allen cluster
> (`test_clustering/reports/`); they are the validation evidence for the `bigcat-alignment` branch.
> Two things to know when reading off-site: absolute paths (`/allen/...`, `/home/...`) refer to the
> Allen filesystem and are kept as provenance, and the harness/probe scripts the reports cite
> (`_merge_R.R`, `_compare_py.py`, ...) live in that working area, not in this repo.

## What we changed
- **[changes.md](changes.md)** — Every source-code and parameter change, **with the measured effect of
  each** (SNN Jaccard graph, DE-score cap at 20, fixed Annoy seed, aligned merge rules #1/#3/#4, and the
  merge-speed fixes). Start here for *what* was modified and *why*.

- **[design_choices.md](design_choices.md)** — **Why the package is built this way**: every design
  decision (exact vs approximate KNN, Leiden, SNN graph, DE cap, merge control flow, the perf rewrites)
  with its rationale, evidence, and trade-off — plus the deliberate divergences from R.

## The analysis, in order
1. **[1_R_vs_Python_early_alignment.md](1_R_vs_Python_early_alignment.md)** — The **early/partial**
   alignment snapshot (stock Python: python-louvain, KNN-edge graph, per-run Annoy, stock merge → cross
   ARI 0.42). Establishes the key lesson: **make each pipeline self-consistent before comparing them.**
2. **[2_self_consistency.md](2_self_consistency.md)** — What made **each pipeline** stable across seeds,
   organized R part / Python part. (R: the `sample_cells` seeding fix; Python: `annoy_trees=50` + fixed
   Annoy seed + Leiden + SNN, 0.51 → 0.91.) Both now ~0.91.
3. **[3_firststep_comparison.md](3_firststep_comparison.md)** — **Cross-pipeline** agreement at one
   clustering step: early 0.47 → fully aligned **~0.88** after the SNN graph + DE-merge alignment.
4. **[4_full_recursive_comparison.md](4_full_recursive_comparison.md)** — The whole recursive pipeline:
   R vs Python `aligned`/`fast` (ARI/NMI, confusion heatmaps), plus the **runtime analysis** and the merge
   bottleneck fix (74 % pandas `.loc` overhead → numpy, ~8× faster end-to-end).
5. **[5_final_merging.md](5_final_merging.md)** — The **final merge** (R `merge_cl_big` vs tc
   `final_merge`): the reduced-space difference (marker genes vs latent), and the isolation test showing
   the two implementations agree at **ARI 0.997** — so the final merge contributes ~nothing to the gap.

6. **[6_de_all_pairs_R_vs_python.md](6_de_all_pairs_R_vs_python.md)** — The **all-pairs DE** step head to
   head on the same clustering. Found **two defects, one per pipeline**: R subsamples 200 cells per
   cluster with a non-reproducible seed, and Python computed detection rates with `np.mean` on a sparse
   mask (one ulp above an exact half, wrongly passing the strict `q1 > 0.5`). With both addressed the two
   agree on **100 % of merge decisions** and **96.5 %** of top-3 marker lists; the only residue is a
   uniform p-value offset between `fast_limma` and `ebayes`.

6b. **[6_full_recursive_comparision_after_fixing_present.md](6_full_recursive_comparision_after_fixing_present.md)**
   — The **whole pipeline re-run** with the `present` fix, versus R before and after final merging. The
   fix removes ~5,971 spurious marginal genes and merges one extra cluster at each stage, but moves
   agreement with R by only ~0.0004 ARI: end-to-end agreement is limited by compounding KNN/Leiden
   nondeterminism, not by DE genes.

7. **[7_marker_queries_R_vs_python.md](7_marker_queries_R_vs_python.md)** — The **marker query layer**
   over `de_parquet` (`select_top_pos_markers_ds`, `select_pos_markers_ds` and helpers), ported from
   `markers.parquet.R`. On the same DE input the two agree: **428/428** clusters identical for the
   top-N family, **12/12 identical marker sets** for the greedy family.

8. **[8_post_clustering_qc_workflow.md](8_post_clustering_qc_workflow.md)** — Replicating the reference
   R **doublet / low-quality** workflow (`post_clustering_qc.R`). Found and fixed **three fidelity
   gaps** on the Python side (one had suppressed doublet detection entirely). Verified against R on the
   same data: doublets-by-triplet (**286, 300, 312**) and low-quality clusters (**154**) match exactly;
   the one marker-doublet difference is **a bug in R**, which finds the cluster but drops its name.

9. **[9_knn_backends.md](9_knn_backends.md)** — The new **`knn_method='pynndescent'`** option measured
   end to end against annoy and R, before and after final merging. pynndescent agrees with R slightly
   better (**ARI 0.840 vs 0.829**) and under-splits less. Self-consistency across seeds: pynndescent
   **0.912**, annoy **0.875**, **R 0.814** — so **both Python backends agree with R more closely than R
   agrees with itself**, and pynndescent is exactly deterministic at a fixed seed.

10. **[10_harmonize_multimodal.md](10_harmonize_multimodal.md)** — The **batch-aware merging port**
    (`merge_cl_multiple` -> `merge_clusters(..., batch_aware_merging=...)`) and its validation on the
    285,230-cell multi-platform WMB2 `TH-EPI-Glut` neighborhood. Five tests: R's own reproducibility
    (nondeterministic solely from a hidden 300-cell statistics subsample; byte-reproducible once
    statistics are supplied), Python's reproducibility (cell-identical across nodes), the isolated
    merge-rule comparisons (**ARI 1.000000** first-step; **identical partitions** in exhaustive mode,
    0.9999981 at R's production `pairBatch=100`), and the end-to-end pipeline comparison against CK's
    production run (ARI 0.64 -- dominated by the different recursive clusterings, not the merge).

## Supporting
- **[three_way_comparison.md](three_way_comparison.md)** — R vs **stock** `transcriptomic_clustering` vs
  the aligned package at the first step (self-consistency + cross).
  *(The former `knn_backends.md` no longer exists — its content, plus the measured 50k–1M scaling
  benchmark, is in [9_knn_backends.md](9_knn_backends.md).)*

## Headline numbers (full recursive run, pre-final-merge)

| Run | Clusters (R=465) | ARI vs R | NMI vs R | Wall-clock |
|---|---:|---:|---:|---|
| R (scrattch.bigcat) | 465 | — | — | 2 h 37 m |
| Python aligned (numpy merge) | 454 | **0.829** | 0.912 | **18.9 min** |
| Python aligned (before merge fix) | 454 | 0.829 | 0.912 | 2 h 33 m |
| Python fast | 481 | 0.761 | 0.902 | 46 min |
| **Python current code (`present` fix)** | **453** | **0.829** | 0.912 | **17.3 min** |

After final merging (R = 438): Python pre-fix 428 clusters at ARI 0.828, current code **427** at
**0.828**. Full pipeline: R **2 h 39 m** vs Python **19 m 42 s** (~8x).

After the final merge: R 465 → **438**, Python 454 → **428**, cross ARI **0.828** (unchanged) — see
[5_final_merging.md](5_final_merging.md).

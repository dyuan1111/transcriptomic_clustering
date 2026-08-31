# Full recursive pipeline after the `present` fix: Python vs R

Re-run of the **whole end-to-end pipeline** with the corrected detection-rate (`present`) calculation,
compared against R at both stages: **before final merging** and **after final merging**.

**Verdict: the fix does exactly what it should to the DE step, and almost nothing to the final
partition.** Agreement with R improves by ~0.0004 ARI at both stages -- real but negligible -- because
whole-pipeline agreement is limited by a completely different bottleneck (compounding KNN/Leiden
nondeterminism across recursion levels), not by marginal DE genes. Cluster counts drop by exactly one at
both stages, in the direction predicted: Python had been **under-calling merges**.

## What changed

The `present` bug and its fix are documented in
[6_de_all_pairs_R_vs_python.md](6_de_all_pairs_R_vs_python.md#the-fix-and-what-it-bought). In short:
detection rates are exact rationals `k/n`, `q1` is filtered with a **strict** `q1 > 0.5`, and
`np.mean()` on a sparse mask returned one ulp *above* an exact half -- so half-detected genes passed in
Python and failed in R. `get_cluster_means_inmemory` now counts and divides once.

| | value |
|---|---|
| Data | spinal-cord neuron subset, 194,221 cells x 17,277 genes, scVI 32-d latent |
| Clustering params | identical to [4_full_recursive_comparison.md](4_full_recursive_comparison.md): k=15, annoy seed 1, `jaccard_snn`, Leiden res 1.0, split_size 10, score_thresh [300,150], `merge_mode='aligned'` |
| Final-merge params | `score_thresh=100`, `merge_mode='aligned'`, k=4 |
| R reference | `full_R_leiden.csv` (465 clusters) and `R_final_merged_latent.csv` (438 clusters) |
| Python **pre-fix** | `clusters_numpy_454.csv` (454) and `Py_final_merged_latent.csv` (428) |
| Python **present-fix** | `clustering_hicatMPI_presentfix/out/clusters_{before,after}_final_merge.csv` |

> **Scope caveat.** The present-fix run used the current `transcriptomic_clustering` repo, which carries
> *all* the changes made during this work -- the annoy per-worker load fix, the `final_merge` marker-space
> option, the numpy `_de_one_pair`, `no_gc_collect`, and the `present` fix. The pre-fix baseline is the
> older pybigcat run. So the deltas below are **"current code vs the old baseline"**, not the `present`
> fix in isolation. The `present` fix is the only one of those expected to change *results*; the others
> were each verified result-neutral when introduced.

## Before final merging

| comparison | clusters | ARI | NMI |
|---|---|---:|---:|
| **R vs Python pre-fix** | 465 vs 454 | 0.8291 | 0.9121 |
| **R vs Python present-fix** | 465 vs **453** | **0.8294** | **0.9123** |
| Python pre-fix vs present-fix | 454 vs 453 | 0.9995 | 0.9988 |

## After final merging

| comparison | clusters | ARI | NMI |
|---|---|---:|---:|
| **R vs Python pre-fix** | 438 vs 428 | 0.8275 | 0.9115 |
| **R vs Python present-fix** | 438 vs **427** | **0.8279** | **0.9119** |
| Python pre-fix vs present-fix | 428 vs 427 | 0.9993 | 0.9984 |

## Runtime

| stage | job | elapsed |
|---|---|---:|
| **R** clustering | 23280208 | **2 h 37 m 03 s** |
| **R** final merge | 25050682 | 2 m 08 s |
| **R total** | | **2 h 39 m 11 s** |
| **Python pre-fix** clustering | 23303147 | 18 m 54 s |
| **Python pre-fix** final merge | 25053435 | 2 m 32 s |
| **Python pre-fix total** | | **21 m 26 s** |
| **Python present-fix** clustering | 25236606 | **17 m 19 s** |
| **Python present-fix** final merge | 25236607 | 2 m 23 s |
| **Python present-fix total** | | **19 m 42 s** |

**Python runs the full pipeline ~8x faster than R** (19 m 42 s vs 2 h 39 m 11 s). The present-fix run is
1 m 35 s faster than the pre-fix run, which is consistent with the fix being **1.83x faster** on the
`present` step itself -- but these runs were months apart on different nodes, and node variation on this
cluster is ~2x (see the runtime caveat in
[6_de_all_pairs_R_vs_python.md](6_de_all_pairs_R_vs_python.md#runtime-full-scale-and-the-effect-of-parallelism)),
so **the 1 m 35 s should not be read as a real speedup**. The 8x gap over R is far larger than node
noise and is real.

## What this means

**1. The fix moved the partition in the predicted direction, by the predicted (small) amount.**
Cluster count fell by exactly one at both stages (454 -> 453, 428 -> 427). Before the fix Python admitted
~5,971 spurious marginal genes at exactly `q1 = 0.5`, making clusters look slightly *more* separable than
they were and therefore **under-calling merges**. Removing them lets slightly more merging happen. That
is the direction predicted in the DE report, now confirmed end-to-end.

**2. Agreement with R barely moved, and that is expected.** ARI improved 0.8291 -> 0.8294 before merging
and 0.8275 -> 0.8279 after. Compare that with the DE-step results, where the same fix took per-pair merge
agreement from 99.982 % to **100.000 %** and marker agreement from 94.9 % to 96.5 %. The difference in
scale is the point:

| level | what limits agreement |
|---|---|
| a **single** DE comparison on a fixed partition | marginal genes -- the `present` fix removes essentially all of it |
| the **whole recursive pipeline** | per-level KNN / graph / Leiden nondeterminism, compounding across every recursion level |

Whole-pipeline ARI (~0.83) is far below the isolated shared-partition merge agreement (~0.95, report 4),
precisely because each recursion level's nondeterminism accumulates. The `present` fix cannot touch that,
so it cannot move whole-pipeline ARI much. **A correctness fix in the DE step is not a lever on
end-to-end clustering agreement.**

**3. Python is highly stable under the change.** Pre-fix vs present-fix agree at ARI 0.9995 (before) and
0.9993 (after) -- the fix perturbs ~0.05 % of the partition. This is a safe change to adopt: it corrects
a real defect without destabilising existing results.

**4. Final merging slightly reduces agreement with R** (0.8294 -> 0.8279 present-fix; 0.8291 -> 0.8275
pre-fix), consistent with report 5 -- the shared merge routine's approximate-vs-exact neighbour search
introduces its own small divergence on top of the clustering difference.

## Reproduce

```bash
# Python, current code (transcriptomic_clustering, not pybigcat):
#   hicatMPI now selects the package via HICAT_PKG (default transcriptomic_clustering);
#   at the time of this run that lived in a separate hicatMPI_tc/, since consolidated
cd clustering_hicatMPI_presentfix && ./submit_pipeline.sh     # run_final_merging=true
# jobs 25236606 (clustering) -> 25236607 (final merge)

cd ..
python _compare_presentfix_recursive.py     # -> presentfix_recursive_metrics.csv
                                            #    presentfix_recursive_runtime.csv
```

R references are unchanged from report 4 (`_run_full_R_leiden.R`, job 23280208) and report 5
(`final_merge_compare/_final_merge_R.R`, job 25050682).

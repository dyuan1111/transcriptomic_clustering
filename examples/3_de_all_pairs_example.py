"""
Step 3 — markers: exhaustive pairwise DE with de_all_pairs, on the FINAL clusters.

Computes DE for every pair of clusters once, writes the results as binned parquet, then queries
pairs and cluster markers from it -- scrattch.bigcat's de_all_pairs + de_parquet/de_summary
workflow: compute once -> persist -> query later. Run it on the final merged clustering from
step 2 (or the hicatMPI pipeline).

Sizing: the per-gene detail costs ~300 bytes/row (all pairs of 454 clusters ~ 30 GB), so results
stream to parquet partitioned by cluster bin; nothing has to fit in memory. Requires pyarrow.
For a small number of pairs there is also an in-memory form (drop out_dir/summary_dir; returns
(summary, detail) DataFrames).

Output-format and query details (column semantics, orientation conventions, efficient lookups):
docs/tests_aligning_to_bigcat/6_de_all_pairs_R_vs_python.md and the de_all_pairs docstrings.

Next step: 4_post_clustering_qc_example.py (doublet / low-quality calls from this parquet).
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering/')
import transcriptomic_clustering as tc
from transcriptomic_clustering.de_all_pairs import (
    select_top_pos_markers_ds, select_pos_markers_ds
)

# Everything this script writes goes under OUT_DIR (the parquet datasets, the marker cache, the
# marker csv). Step 4 reads the same directories.
OUT_DIR = '/path/to/your/de_out'
DE_PARQUET = os.path.join(OUT_DIR, 'de_parquet')     # per-gene detail   (R out.dir)
DE_SUMMARY = os.path.join(OUT_DIR, 'de_summary')     # per-pair summary  (R summary.dir)
os.makedirs(OUT_DIR, exist_ok=True)

# One threshold dict, mirroring scrattch.bigcat's de_param(). Log-scale thresholds are
# ln-converted: R normalizes with log2(CPM+1), this package with ln(CPM+1), so R's low.th=1 and
# lfc.th=1 (log2) both become ln(2)=0.6931472. de_all_pairs consumes the gene-level filters;
# score_thresh/min_genes are carried for the downstream merge/QC steps.
THRESHOLDS = {
    'low_thresh': 0.6931472,      # R low.th    = 1 (log2)
    'padj_thresh': 0.01,          # R padj.th
    'lfc_thresh': 0.6931472,      # R lfc.th    = 1 (log2)
    'q1_thresh': 0.5,             # R q1.th
    'q2_thresh': None,            # R q2.th
    'qdiff_thresh': 0.7,          # R q.diff.th
    'cluster_size_thresh': 10,    # R min.cells
    'score_thresh': 150,          # R de.score.th
    'min_genes': 5,               # R min.genes
}

adata = sc.read('/path/to/normalized_adata.h5ad')          # adata.X already ln(CPM+1) normalized

# The clustering csv is KEYED BY CELL NAME (sample_id, cl) -- from 2_final_merging_example.py
# (clusters_final.csv) or the hicatMPI pipeline (out/clusters_after_final_merge.csv). Join by name,
# never by row position, and fail loudly if adata and clustering are out of sync.
labels = pd.read_csv('/path/to/your/clusters_final.csv', index_col=0)['cl'].astype(str)
labels.index = labels.index.astype(str)
labels = labels.reindex(adata.obs_names)
if labels.isna().any():
    raise ValueError(f"{int(labels.isna().sum())} cells in the adata have no cluster label "
                     f"-- the adata and the clustering are out of sync")
cluster_assignments = {c: np.flatnonzero((labels == c).to_numpy()) for c in labels.unique()}

# de_all_pairs works from cluster statistics, not the cell matrix
cluster_by_obs = np.empty(adata.n_obs, dtype=object)
for label, idx in cluster_assignments.items():
    cluster_by_obs[idx] = label
means, present, variances = tc.get_cluster_means(
    adata, cluster_assignments, cluster_by_obs, low_th=THRESHOLDS['low_thresh']
)
cl_size = {label: len(idx) for label, idx in cluster_assignments.items()}

# ---- compute all pairs, stream to parquet --------------------------------------------------------
# cl_bin pins the on-disk layout (R cl.bin, cl.bin.size=100): partitions stay identical across
# runs, and de_all_pairs saves it to <dir>/_cl_bin.json so read_de_pairs finds it automatically.
cl_bin = tc.make_cl_bin(cluster_assignments.keys(), bin_size=100)
tc.de_all_pairs(
    means, variances, present, cl_size, THRESHOLDS,
    top_n=500,                   # R de_selected_pairs default
    out_dir=DE_PARQUET,
    summary_dir=DE_SUMMARY,
    cl_bin=cl_bin,
    n_jobs=16,                   # bin-pairs computed in parallel (R mc.cores)
)

# ---- query it -------------------------------------------------------------------------------------
# summary: one row per pair [pair, P1, P2, up_num, down_num, num, up_score, down_score, score];
# small enough to load whole. detail: one row per gene; query specific pairs so only their
# partitions are read (either orientation works).
summary = tc.read_de_pairs(DE_SUMMARY)
print(summary.nsmallest(10, 'score'))                              # most-similar pairs
one_pair = tc.read_de_pairs(DE_PARQUET, pairs=[('3', '17')])     # per-gene detail for one pair

# ---- cluster markers, via the ported R-validated functions ---------------------------------------
# (validated against scrattch.bigcat: 428/428 clusters identical top-N lists, 12/12 identical
# greedy marker sets -- docs/tests_aligning_to_bigcat/7_marker_queries_R_vs_python.md)
all_clusters = list(cluster_assignments.keys())

# top-N positive markers per cluster (R select_top_pos_markers_ds)
top_markers = select_top_pos_markers_ds(DE_PARQUET, clusters=all_clusters, n_markers=3)

# coverage-based positive markers (R select_pos_markers_ds): keeps adding genes until every
# cluster-vs-other comparison is covered n_markers times. n_markers=1, max_genes=20 is the exact
# configuration validated against R (12/12 identical marker sets); out_dir caches per-cluster
# results so a long run can resume. Optionally pass cl_means=means.T to enable R's fold-change
# check (a supported branch, not part of the validated configuration).
pos_markers = select_pos_markers_ds(DE_PARQUET, clusters=all_clusters, n_markers=1,
                                    max_genes=20, out_dir=os.path.join(OUT_DIR, 'marker_cache'))

pd.Series({c: ','.join(g) for c, g in top_markers.items()}, name='markers') \
  .rename_axis('cl').to_csv(os.path.join(OUT_DIR, 'cluster_top_markers.csv'))

# OPTIONAL utility (not part of the numbered workflow) -- pairwise DEGs for a handful of
# clusters, in memory, without the parquet layer. This is the Python analog of scrattch.bigcat's
# select_markers (markers.R:13): DE over every pair of the given clusters, returning the union
# of the top n_markers up + n_markers down genes per pair. For all-pairs DE persisted to disk
# and per-cluster marker queries, use 3_de_all_pairs_example.py instead.
import sys

import numpy as np
import pandas as pd
import scanpy as sc

# Skip this line if transcriptomic_clustering is installed
sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering/')
from transcriptomic_clustering.pairwise_DEGs import pairwise_degs

thresholds = {
    'q1_thresh': 0.5,
    'q2_thresh': None,
    'cluster_size_thresh': 10,    # R min.cells
    'qdiff_thresh': 0.7,
    'padj_thresh': 0.01,
    'lfc_thresh': 0.6931472,      # ln(2) = R lfc.th=1 (log2) -- adata.X is ln(CPM+1)
    'score_thresh': 100,
    'low_thresh': 0.6931472,      # ln(2) = R low.th=1 (log2)
    'min_genes': 5
}

# adata.X must already be ln(CPM+1) normalized
adata = sc.read('/path/to/normalized_adata.h5ad')

# Load the clustering (sample_id -> cl, from 1_clustering_example.py / 2_final_merging_example.py /
# the hicatMPI pipeline) and join it by CELL NAME.
labels = pd.read_csv('/path/to/your/clusters_final.csv', index_col=0)['cl'].astype(str)
labels.index = labels.index.astype(str)
labels = labels.reindex(adata.obs_names)
if labels.isna().any():
    raise ValueError("adata and clustering are out of sync (cells without labels)")

# Pick the clusters to compare and subset the adata to just their cells.
CLUSTERS = ['1', '2']
keep = labels.isin(CLUSTERS).to_numpy()
adata = adata[keep, :].copy()
labels = labels[keep]
obs_by_cluster = {c: np.flatnonzero((labels == c).to_numpy()) for c in CLUSTERS}

# returns a set of markers (20 up regulated and 20 down regulated per pair)
degs = pairwise_degs(
    adata,
    obs_by_cluster,
    thresholds,
    n_markers=20,
    de_method='ebayes',
    n_jobs=30
)
print(degs)

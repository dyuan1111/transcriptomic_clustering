"""
Step 1 — iterative (recursive) clustering with transcriptomic_clustering, single process.

This is the plain-tc version: one Python process runs the whole recursion depth-first
(`iter_clust`). It is the right tool for small-to-medium datasets; for large datasets on an HPC
cluster, use the MPI pipeline instead (hicatMPI: submit_pipeline_pooled.sh /
submit_pipeline_batchaware.sh), which runs this same algorithm manager/worker-parallel and
produces the same partition.

Input : an AnnData (raw counts or ln(CPM+1)-normalized) and, optionally, a latent space (csv).
Output: clusters.csv -- one row per cell (sample_id, cl). Keyed by CELL NAME, not by index, so it
        cannot be silently mis-paired with a reordered/subset adata downstream.

Next step: 2_final_merging_example.py
"""
import sys

import numpy as np
import pandas as pd
import scanpy as sc

# Skip this line if transcriptomic_clustering is installed
sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering/')

import transcriptomic_clustering as tc
from transcriptomic_clustering.iterative_clustering import (
    build_cluster_dict, iter_clust, OnestepKwargs
)

# Load the data. adata.X may hold raw counts or ln(CPM+1)-normalized values.
adata = sc.read('path/to/your/data.h5ad')

# Normalize to ln(CPM+1). Skip if adata.X is already normalized -- but it must be NATURAL log
# (equivalent to sc.pp.normalize_total(target_sum=1e6) + sc.pp.log1p): the thresholds below are in
# ln units (ln(2) = 0.6931472 corresponds to R's 1 in log2), and log2 data would silently shift
# every threshold.
adata = tc.normalize(adata)

# Add the scVI latent space. Skip if adata.obsm['X_scVI'] is already present; leave
# latent_component=None below to cluster on PCA instead of an embedding.
scvi = pd.read_csv('path/to/scvi_latent_space.csv', index_col=0)
adata.obsm['X_scVI'] = np.asarray(scvi.loc[adata.obs_names])

# Set up the clustering parameters
def setup_transcriptomic_clustering():
    means_vars_kwargs = {          # PCA path only (ignored when clustering on a latent)
        'low_thresh': 0.6931472,
        'min_cells': 4
    }
    highly_variable_kwargs = {     # PCA path only
        'max_genes': 4000
    }
    pca_kwargs = {                 # PCA path only
        'cell_select': 30000,
        'n_comps': 50,
        'svd_solver': 'randomized'
    }
    filter_pcs_kwargs = {          # PCA path only
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore',
        'zth': 2,
        'max_pcs': None,
    }
    filter_known_modes_kwargs = {
        # 'known_modes': known_modes_df,
        'similarity_threshold': 0.7
    }
    latent_kwargs = {
        'latent_component': "X_scVI"   # None = PCA path; a str = adata.obsm key to cluster on
    }
    cluster_louvain_kwargs = {
        'k': 15,
        'nn_measure': 'euclidean',
        'knn_method': 'pynndescent',   # default; 'annoy' matches scrattch.bigcat (then also set
                                       # 'annoy_trees': 50 -- R's fixed value)
        'knn_seed': 1,                 # fixed -> reproducible KNN graph, decoupled from random_seed
        'louvain_method': 'vtraag',    # Leiden
        'weighting_method': 'jaccard_snn',
        'n_jobs': 30,
        'resolution': 1.0,
    }
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10,     # R min.cells
            'qdiff_thresh': 0.7,
            'padj_thresh': 0.01,
            'lfc_thresh': 0.6931472,       # ln(2) = R lfc.th=1 (log2)
            'score_thresh': 100,           # R de.score.th
            'low_thresh': 0.6931472,       # ln(2) = R low.th=1 (log2)
            'min_genes': 5
        },
        'k': 4,
        'de_method': 'ebayes',
        # Multi-platform / multi-modality data only. Leave this out (the default) and clusters are
        # merged on pooled statistics, which cannot tell a real difference from a platform artefact.
        # Set it and each batch is tested separately: genes whose fold change flips direction between
        # batches are discarded, and a pair merges only if no batch can still tell the two apart.
        # 'batch_aware_merging': {
        #     'batch_obs': 'platform',        # column of adata.obs naming each cell's batch
        #     'lfc_conservation_th': 0.7,     # fraction of batches that must agree on a gene
        #     'conservation_lfc_th': 0.6931472,  # optional; the fold change a gene must clear to
        #                                        # count as agreeing (default: thresholds['lfc_thresh'])
        #     'thresholds': {                 # optional per-batch overrides; unlisted batches
        #         '10X_nuclei_v3': {'q1_thresh': 0.3},  # inherit the values above (WMB2 production:
        #     },                              # nuclei sample less RNA -> lower detection bar)
        # },
    }
    onestep_kwargs = OnestepKwargs(
        means_vars_kwargs = means_vars_kwargs,
        highly_variable_kwargs = highly_variable_kwargs,
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        latent_kwargs = latent_kwargs,
        cluster_louvain_kwargs = cluster_louvain_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs
    )
    return onestep_kwargs

onestep_kwargs = setup_transcriptomic_clustering()

# Run the iterative clustering. Needs a tmp folder for intermediate results.
# `clusters` is a list of cell-index lists (one per cluster); `markers` is the union of the DE
# genes that drove the splits (a diagnostic -- proper markers are computed on the FINAL clusters
# as a separate step, see 3_de_all_pairs_example.py).
clusters, markers = iter_clust(
    adata,
    min_samples=10,                # min cluster size to keep sub-clustering (R split.size)
    onestep_kwargs=onestep_kwargs,
    random_seed=123,
    tmp_dir="/path/to/your/tmp"
)

# Save the clustering KEYED BY CELL NAME (sample_id -> cl). This csv is the artifact every later
# step consumes (final merge, DE, QC); index-based formats are avoided because they fail silently
# when paired with a reordered or subset adata.
cl = pd.Series(index=adata.obs_names, dtype=object, name='cl')
for i, cells in enumerate(clusters):
    cl.iloc[cells] = i + 1
cl.index.name = 'sample_id'
cl.to_csv('/path/to/your/clusters.csv')

# Optional: the split-time marker union, one gene per line (diagnostic only; see docstring above).
pd.Series(sorted(markers), name='gene').to_csv('/path/to/your/markers_split_time.csv', index=False)

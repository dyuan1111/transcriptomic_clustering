"""
Step 2 — final merging: the global DE-based cleanup after recursive clustering.

All clusters from the whole recursion are re-examined together and near-duplicates are merged
(same merge rule as the per-step merge, typically at a stricter score threshold).

Input : the AnnData (same normalization as step 1) and clusters.csv from 1_clustering_example.py
        -- keyed by CELL NAME, joined to the adata by name, refusing to run if they are out of sync.
Output: clusters_final.csv (sample_id, cl).

Next step: 3_de_all_pairs_example.py (markers on the final clusters).
"""
import sys

import numpy as np
import pandas as pd
import scanpy as sc

# Skip this line if transcriptomic_clustering is installed
sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering/')

from transcriptomic_clustering.final_merging import final_merge, FinalMergeKwargs

# Load the data; adata.X must be ln(CPM+1)-normalized, exactly as in step 1.
adata = sc.read('path/to/your/data.h5ad')

# Add the scVI latent space. Skip if adata.obsm['X_scVI'] is already present; set
# latent_component=None below to shortlist merge candidates on PCA-of-markers instead.
scvi = pd.read_csv('path/to/scvi_latent_space.csv', index_col=0)
adata.obsm['X_scVI'] = np.asarray(scvi.loc[adata.obs_names])

# ---- load the clustering by CELL NAME, never by raw index ---------------------------------------
# A name-keyed join is order-independent and fails loudly if the adata and the clustering have been
# decoupled (reordered, subset, or simply the wrong file).
cl_df = pd.read_csv('/path/to/your/clusters.csv', index_col=0)
cl_df.index = cl_df.index.astype(str)
labels = cl_df['cl'].reindex(adata.obs_names)
if labels.isna().any():
    raise ValueError(f"{int(labels.isna().sum())} cells in the adata have no label in clusters.csv "
                     f"-- the adata and the clustering are out of sync")
if len(cl_df.index.difference(adata.obs_names)):
    raise ValueError("clusters.csv contains cells absent from the adata "
                     "-- the adata and the clustering are out of sync")
codes, uniques = pd.factorize(labels.values)
# lists, not arrays: merge_clusters mutates assignments with list.extend()
clusters = [np.flatnonzero(codes == i).tolist() for i in range(len(uniques))]

# ---- parameters ----------------------------------------------------------------------------------
# The pca/filter/project blocks are used ONLY when latent_component is None (candidate shortlist on
# PCA-of-markers); on the latent path they are ignored and may be omitted.
def setup_merging():
    pca_kwargs = {
        # 'cell_select': 30000,  # do NOT set for final merging -- cells are sampled per cluster
        'n_comps': 50,
        'svd_solver': 'randomized'
    }
    filter_pcs_kwargs = {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore',
        'zth': 2,
        'max_pcs': None}
    filter_known_modes_kwargs = {
        'known_modes': 'log2ngene',
        'similarity_threshold': 0.7}
    project_kwargs = {}
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10,     # R min.cells
            'qdiff_thresh': 0.7,
            'padj_thresh': 0.01,
            'lfc_thresh': 0.6931472,       # ln(2) = R lfc.th=1 (log2)
            'score_thresh': 100,           # R de.score.th for the final merge
            'low_thresh': 0.6931472,       # ln(2) = R low.th=1 (log2)
            'min_genes': 5
        },
        'k': 4,
        'de_method': 'ebayes',
        'n_markers': None,  # None bypasses the marker calculation (the slow step); markers are
                            # computed properly on the final clusters in step 3 anyway
        # Multi-platform / multi-modality data: same key and meaning as in step 1's
        # merge_clusters_kwargs -- makes the FINAL merge batch-aware; the reduced space that
        # shortlists candidate pairs is unchanged.
        # 'batch_aware_merging': {
        #     'batch_obs': 'platform',
        #     'lfc_conservation_th': 0.7,
        #     'thresholds': {'10X_nuclei_v3': {'q1_thresh': 0.3}},
        # },
    }
    latent_kwargs = {
        'latent_component': "X_scVI"   # None = PCA-of-markers; a str = adata.obsm key
    }

    return FinalMergeKwargs(
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        project_kwargs = project_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs,
        latent_kwargs = latent_kwargs
    )

merge_kwargs = setup_merging()

# Run the final merging. marker_genes is required only for the PCA-of-markers shortlist
# (latent_component=None); load it from step 1's markers_split_time.csv in that case:
#   markers = set(pd.read_csv('/path/to/your/markers_split_time.csv')['gene'])
clusters_after_merging, markers_after_merging = final_merge(
    adata,
    clusters,
    n_samples_per_clust=20,
    random_seed=2024,
    n_jobs=30,                     # number of cores to use
    return_markers_df=False,
    final_merge_kwargs=merge_kwargs
)

# Save KEYED BY CELL NAME, as in step 1.
out = pd.Series(index=adata.obs_names, dtype=object, name='cl')
for i, cells in enumerate(clusters_after_merging):
    out.iloc[list(cells)] = i + 1
out.index.name = 'sample_id'
out.to_csv('/path/to/your/clusters_final.csv')

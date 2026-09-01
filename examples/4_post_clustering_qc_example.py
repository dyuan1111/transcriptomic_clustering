"""
Step 4 — post-clustering QC: doublet and low-quality cluster calls.

Replicates the reference R workflow (`clustering_bigcat/post_clustering_qc.R`, "label doublet and
low quality clusters") with the ported find_* functions; the joins/filters between them are the R
workflow's own script-level logic, annotated with the R lines they mirror. Validated against R on
the same data (docs/tests_aligning_to_bigcat/8_post_clustering_qc_workflow.md): doublet and
low-quality cluster calls match exactly.

Inputs : the de_parquet / de_summary dataset from 3_de_all_pairs_example.py, the clustering csv it
         was built from, the adata (for cluster means), and a per-cell QC table
         (sample_id, nGene, nCount, pct_counts_mt -- R's qc.csv).
Outputs: R's final_result/*.csv -- doublets by triplet, doublets by marker, low-quality clusters.
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering/')
import transcriptomic_clustering as tc
from transcriptomic_clustering.post_clustering_qc import (
    find_doublet_by_marker, find_doublets, find_low_quality, find_triplets,
)

DE_DIR     = "/path/to/your/de_out"               # OUT_DIR of 3_de_all_pairs_example.py
DE_PARQUET = os.path.join(DE_DIR, "de_parquet")   # per-gene detail   (R de.dir)
DE_SUMMARY = os.path.join(DE_DIR, "de_summary")   # per-pair summary  (R summary.dir)
CLUSTERS   = "/path/to/your/clusters_final.csv"   # sample_id -> cl, as in step 3
QC_CSV     = "/path/to/your/qc.csv"               # sample_id, nGene, nCount, pct_counts_mt
OUT_DIR    = "/path/to/your/final_result"         # this script's outputs (R final_result/)
os.makedirs(OUT_DIR, exist_ok=True)

# ---- clustering (joined by cell name), cl.bin, all pairs ----------------------------------------
adata = sc.read('/path/to/normalized_adata.h5ad')
cl = pd.read_csv(CLUSTERS, index_col=0)['cl'].astype(str)
cl.index = cl.index.astype(str)
labels = cl.reindex(adata.obs_names)
if labels.isna().any():
    raise ValueError("adata and clustering are out of sync (cells without labels)")
cn = sorted(labels.unique(), key=lambda c: int(c) if c.isdigit() else c)
# cl.bin: LOAD the one saved inside the parquet dataset (step 3 wrote <de_parquet>/_cl_bin.json --
# the Python equivalent of R's cl.bin.rda). Never rebuild it here: the saved map is the dataset's
# actual on-disk layout.
cl_bin = tc.load_cl_bin(DE_PARQUET)
all_pairs = pd.DataFrame(tc.create_pairs(cn), columns=['P1', 'P2'])
all_pairs['pair'] = all_pairs.P1 + '_' + all_pairs.P2
all_pairs['pair_id'] = np.arange(1, len(all_pairs) + 1)

# cluster means (clusters x genes), for the marker-based doublet check
cluster_assignments = {c: np.flatnonzero((labels == c).to_numpy()) for c in cn}
cluster_by_obs = labels.to_numpy()
cluster_means, _, _ = tc.get_cluster_means(adata, cluster_assignments, cluster_by_obs,
                                           low_th=0.6931472)

# per-cluster QC statistics (R: anno.df %>% group_by(cl) %>% summarize(...))
qc = pd.read_csv(QC_CSV)
anno = pd.DataFrame({'sample_id': cl.index, 'cl': cl.values}).merge(qc, on='sample_id', how='left')
cl_qc_stats = anno.groupby('cl').agg(
    gene_counts=('nGene', 'mean'), umi_counts=('nCount', 'mean'),
    mt_pct=('pct_counts_mt', 'mean'), size=('sample_id', 'size'),
).reset_index()

# ---- 1. doublets by triplet ----------------------------------------------------------------------
summary = tc.read_de_pairs(DE_SUMMARY)
summary['P1'] = summary.P1.astype(str); summary['P2'] = summary.P2.astype(str)
if 'pair' not in summary.columns:
    summary['pair'] = summary.P1 + '_' + summary.P2

triplets = find_triplets(summary, all_pairs=all_pairs)          # R find_triplets_big
# root= streams each triplet's rows from its clusters' partitions (R find_doublets_all_big)
doublets = find_doublets(None, triplets, root=DE_PARQUET, cl_bin=cl_bin,
                         top_n=50, score_th=0.8, olap_th=1.6)
doublets = doublets.merge(cl_qc_stats, on='cl', how='left')
# R workflow line 84 selects at olap > 1.4 while find_doublets stops early at 1.6 -- an asymmetry
# faithfully reproduced; it works because find_doublets returns every tested triplet.
cl_doublets = sorted(set(doublets[(doublets.score > 0.8) &
                                  (doublets.olap_ratio_up_1 + doublets.olap_ratio_down_1 > 1.4)]['cl']))
print(f"doublet clusters (by triplet): {len(cl_doublets)}")

# ---- 2. doublets by marker -----------------------------------------------------------------------
markers = {  # broad-class markers (R's list); missing genes are filtered automatically
    'Neuron': ['Slc32a1', 'Slc17a7', 'Slc17a6'], 'Olig': ['Sox10', 'Opalin'], 'Astro': ['Aqp4'],
    'Endo': ['Ly6c1'], 'VLMC': ['Slc6a13'], 'Peri': ['Kcnj8'], 'SMC': ['Acta2'], 'Micro': ['C1qc'],
}
doublet_means = find_doublet_by_marker(cluster_means, markers=markers, th=3.5).round(2)
cl_doublets_by_marker = sorted(doublet_means.index.astype(str))
print(f"doublet clusters (by marker): {len(cl_doublets_by_marker)}")

# ---- 3. low-quality clusters ---------------------------------------------------------------------
low_df = find_low_quality(summary, low_th=2)                    # R find_low_quality_big
low_df = low_df[~low_df.cl.isin(cl_doublets) & ~low_df.cl_low.isin(cl_doublets)]
low_df = (low_df.merge(cl_qc_stats, on='cl', how='left')
                .merge(cl_qc_stats, left_on='cl_low', right_on='cl', how='left',
                       suffixes=('_x', '_y')))
# R workflow lines 120/129: cl_low reads much less RNA than its counterpart, and is much smaller
select_low = low_df[(low_df.umi_counts_y * 1.5 < low_df.umi_counts_x) &
                    (low_df.gene_counts_y < 4000)]
cl_low_df = select_low.groupby('cl_low', as_index=False).agg(size=('size_x', 'sum'))
cl_size = anno.groupby('cl', as_index=False).agg(cl_size=('sample_id', 'size'))
cl_low_df = cl_low_df.merge(cl_size, left_on='cl_low', right_on='cl', how='left') \
                     .rename(columns={'size': 'cl_low_size'})
cl_low = sorted(set(cl_low_df[cl_low_df.cl_size > 1.5 * cl_low_df.cl_low_size]['cl_low']))
print(f"low-quality clusters: {len(cl_low)}")

# ---- export: the same three CSVs the R workflow writes -------------------------------------------
pd.DataFrame({'cl': cl_doublets}).to_csv(f"{OUT_DIR}/doubletsClusters.csv", index=False)
pd.DataFrame({'cl': cl_doublets_by_marker}).to_csv(f"{OUT_DIR}/doubletsClusters_byMarker.csv", index=False)
pd.DataFrame({'cl': cl_low}).to_csv(f"{OUT_DIR}/lowQClusters.csv", index=False)

# Parameter provenance (all matching the R workflow / validated in report 8):
#   find_triplets defaults (min_up_num=30, max_down_num=10, min_de_num=50); find_doublets
#   top_n=50/score_th=0.8/olap_th=1.6; selection score>0.8 & olap sum>1.4; marker th=3.5;
#   find_low_quality low_th=2; low selection umi_y*1.5<umi_x & gene_y<4000, then size>1.5x.

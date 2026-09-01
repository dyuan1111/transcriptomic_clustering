"""
Post-clustering QC: doublet and low-quality cluster calls.

A line-for-line replication of the reference R workflow
(`clustering_bigcat/post_clustering_qc.R`, the "label doublet and low quality clusters" section),
using `transcriptomic_clustering`. Each block is annotated with the R it corresponds to.

Inputs
------
  * a written de_parquet / de_summary dataset (see de_all_pairs_example.py)
  * the cluster assignment used to build it
  * a per-cell QC table (nGene, nCount, pct_counts_mt) -- R's qc.csv
  * cluster means, for the marker-based doublet check

Outputs (matching R's final_result/*.csv)
  * doublet clusters, by triplet
  * doublet clusters, by marker
  * low-quality clusters
"""
import os
import numpy as np
import pandas as pd
import transcriptomic_clustering as tc
from transcriptomic_clustering.post_clustering_qc import (
    find_doublet_by_marker, find_doublets, find_low_quality, find_triplets,
)

# ---- paths ------------------------------------------------------------------------------------
DE_PARQUET = "de_parquet"          # per-gene detail   (R de.dir)
DE_SUMMARY = "de_summary"          # per-pair summary  (R summary.dir)
CLUSTERS   = "clusters.csv"        # sample_id -> cl   (R merge.result$cl)
QC_CSV     = "qc.csv"              # sample_id, nGene, nCount, pct_counts_mt
OUT_DIR    = "final_result"

os.makedirs(OUT_DIR, exist_ok=True)

# ---- cluster assignment + cl.bin ---------------------------------------------------------------
# R:  cl = merge.result$cl
#     cn <- as.character(sort(unique(cl)))
#     cl.bin = data.frame(cl=cn, bin = ceiling((1:length(cn)/cl.bin.size)))
cl = pd.read_csv(CLUSTERS, index_col=0).iloc[:, -1].astype(str)
cn = sorted(cl.unique(), key=lambda c: int(c) if c.isdigit() else c)
cl_bin = tc.make_cl_bin(cn, bin_size=100)

# R:  pairs = create_pairs(cn); all.pairs = pairs
all_pairs = pd.DataFrame(tc.create_pairs(cn), columns=['P1', 'P2'])
all_pairs['pair'] = all_pairs.P1 + '_' + all_pairs.P2
all_pairs['pair_id'] = np.arange(1, len(all_pairs) + 1)

# ---- per-cluster QC statistics -----------------------------------------------------------------
# R:  anno.df = data.frame(sample_id=names(cl), cl=as.character(cl)) %>% left_join(qc)
#     cl.qc.stats = anno.df %>% group_by(cl) %>% summarize(gene.counts=mean(nGene), ...)
qc = pd.read_csv(QC_CSV)
anno = pd.DataFrame({'sample_id': cl.index, 'cl': cl.values}).merge(qc, on='sample_id', how='left')
cl_qc_stats = anno.groupby('cl').agg(
    gene_counts=('nGene', 'mean'),
    umi_counts=('nCount', 'mean'),
    mt_pct=('pct_counts_mt', 'mean'),
    size=('sample_id', 'size'),
).reset_index()

# =================================================================================================
# 1. Doublets by triplet
# =================================================================================================
# R:  triplets = find_triplets_big(de.summary="de_summary", all.pairs=all.pairs)
summary = tc.read_de_pairs(DE_SUMMARY)
summary['P1'] = summary.P1.astype(str)
summary['P2'] = summary.P2.astype(str)
if 'pair' not in summary.columns:
    summary['pair'] = summary.P1 + '_' + summary.P2

triplets = find_triplets(summary, all_pairs=all_pairs)
print(f"triplets: {len(triplets):,} over {triplets.cl_up.nunique()} candidate clusters")

# R:  tmp = find_doublets_all_big(de.dir="de_parquet", triplets=triplets, cl.bin=cl.bin, ...)
#     doublets = open_dataset("doublets_result") %>% collect() %>% arrange(-score)
# Pass root= so each triplet's rows are read from its three clusters' partitions, exactly as R
# does -- the whole detail does not have to fit in memory.
doublets = find_doublets(None, triplets, root=DE_PARQUET, cl_bin=cl_bin,
                         top_n=50, score_th=0.8, olap_th=1.6)
doublets = doublets.sort_values('score', ascending=False)

# R:  doublets %>% left_join(cl.qc.stats, by="cl") %>% left_join(..., by=c("cl1"="cl")) %>% ...
doublets = (doublets
            .merge(cl_qc_stats, left_on='cl', right_on='cl', how='left')
            .merge(cl_qc_stats, left_on='cl1', right_on='cl', how='left', suffixes=('', '_cl1'))
            .merge(cl_qc_stats, left_on='cl2', right_on='cl', how='left', suffixes=('', '_cl2')))

# R:  cl.doublets = doublets %>% filter(score > 0.8 & olap.ratio.up.1 + olap.ratio.down.1 > 1.4)
#                            %>% pull(cl) %>% unique
# NOTE the asymmetry, faithfully reproduced: find_doublets STOPS at olap > 1.6, but the selection
# below uses 1.4 -- so triplets that scored between the two are kept as evidence but never
# triggered an early stop. This only works because find_doublets returns every tested triplet.
cl_doublets = sorted(set(doublets[(doublets.score > 0.8) &
                                  (doublets.olap_ratio_up_1 + doublets.olap_ratio_down_1 > 1.4)]['cl']))
print(f"doublet clusters (by triplet): {len(cl_doublets)}")

# =================================================================================================
# 2. Doublets by marker
# =================================================================================================
# R:  markers = list(Neuron=c("Slc32a1","Slc17a7","Slc17a6"), Olig=c("Sox10","Opalin"), ...)
#     markers_ = lapply(markers, toupper)   # the reference data uses upper-case gene symbols
#     doublet.means = find_doublet_by_marker(cl.means, markers=markers_)
markers = {
    'Neuron': ['Slc32a1', 'Slc17a7', 'Slc17a6'],
    'Olig':   ['Sox10', 'Opalin'],
    'Astro':  ['Aqp4'],
    'Endo':   ['Ly6c1'],
    'VLMC':   ['Slc6a13'],
    'Peri':   ['Kcnj8'],
    'SMC':    ['Acta2'],
    'Micro':  ['C1qc'],
}
# cluster_means is clusters x genes. R drops markers missing from the matrix by hand
# (its `rm_` step); find_doublet_by_marker does that filtering for you.
cluster_means = pd.read_pickle("cluster_means.pkl")     # or recompute with tc.get_cluster_means
doublet_means = find_doublet_by_marker(cluster_means, markers=markers, th=3.5).round(2)
cl_doublets_by_marker = sorted(doublet_means.index.astype(str))
print(f"doublet clusters (by marker): {len(cl_doublets_by_marker)}")

# =================================================================================================
# 3. Low-quality clusters
# =================================================================================================
# R:  low.df = find_low_quality_big(ds.summ, pairs=all.pairs)
#     low.df = low.df %>% filter(!cl %in% cl.doublets & !cl.low %in% cl.doublets)
low_df = find_low_quality(summary, low_th=2)
low_df = low_df[~low_df.cl.isin(cl_doublets) & ~low_df.cl_low.isin(cl_doublets)]

# R:  low.df %>% left_join(cl.qc.stats) %>% left_join(cl.qc.stats, by=c("cl.low"="cl"))
#     -> .x columns describe `cl` (the good cluster), .y columns describe `cl.low`
low_df = (low_df
          .merge(cl_qc_stats, left_on='cl', right_on='cl', how='left')
          .merge(cl_qc_stats, left_on='cl_low', right_on='cl', how='left', suffixes=('_x', '_y')))

# R:  select.low.df = low.df %>% filter((umi.counts.y * 1.5 < umi.counts.x) & (gene.counts.y < 4000))
select_low = low_df[(low_df.umi_counts_y * 1.5 < low_df.umi_counts_x) &
                    (low_df.gene_counts_y < 4000)]

# R:  cl.low.df = select.low.df %>% group_by(cl.low) %>% summarize(size=sum(size.x))
#     cl.low.df = cl.low.df %>% left_join(cl.size, by=c("cl.low"="cl")) %>% rename(cl.low.size=size)
#     cl.low = cl.low.df %>% filter(cl.size > 1.5*cl.low.size) %>% pull(cl.low) %>% unique
cl_low_df = select_low.groupby('cl_low', as_index=False).agg(size=('size_x', 'sum'))
cl_size = anno.groupby('cl', as_index=False).agg(cl_size=('sample_id', 'size'))
cl_low_df = cl_low_df.merge(cl_size, left_on='cl_low', right_on='cl', how='left')
cl_low_df = cl_low_df.rename(columns={'size': 'cl_low_size'})
cl_low = sorted(set(cl_low_df[cl_low_df.cl_size > 1.5 * cl_low_df.cl_low_size]['cl_low']))
print(f"low-quality clusters: {len(cl_low)}")

# =================================================================================================
# 4. Export -- the same three CSVs the R workflow writes
# =================================================================================================
pd.DataFrame({'cl': cl_doublets}).to_csv(f"{OUT_DIR}/doubletsClusters.csv", index=False)
pd.DataFrame({'cl': cl_doublets_by_marker}).to_csv(f"{OUT_DIR}/doubletsClusters_byMarker.csv", index=False)
pd.DataFrame({'cl': cl_low}).to_csv(f"{OUT_DIR}/lowQClusters.csv", index=False)
print(f"wrote {OUT_DIR}/doubletsClusters.csv, doubletsClusters_byMarker.csv, lowQClusters.csv")

# ---- notes --------------------------------------------------------------------------------------
#
# Thresholds, and where they come from (R post_clustering_qc.R):
#   find_triplets    min_up_num=30, max_down_num=10, min_de_num=50   -- R find_triplets_big defaults
#   find_doublets    top_n=50, score_th=0.8, olap_th=1.6             -- R find_doublets_all_big defaults
#   selection        score > 0.8 AND olap_up_1 + olap_down_1 > 1.4   -- R workflow line 84 (NOT 1.6)
#   find_doublet_by_marker  th=3.5                                   -- R default
#   find_low_quality low_th=2                                        -- R default
#   low-quality selection   umi_y*1.5 < umi_x AND gene_y < 4000,
#                           then cl_size > 1.5 * cl_low_size         -- R workflow lines 120/129
#
# Triplets are emitted twice, mirrored (cl_down_x/cl_down_y swapped) -- R's self-join does the same
# and check_triplet gives the same answer either way. It doubles the work but not the result.

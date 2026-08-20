"""
Exhaustive pairwise DE with de_all_pairs.

Computes DE for every pair of clusters once, writes the results as binned parquet, then reads
individual pairs back cheaply. This mirrors scrattch.bigcat's de_all_pairs + de_parquet/de_summary
workflow: compute once -> persist -> query any pair later.

Two ways to run it:
  1. in memory  -- returns (summary, detail) DataFrames. Fine for a modest number of pairs.
  2. to disk    -- pass out_dir/summary_dir. Results stream to parquet partitioned by cluster bin,
                   so an all-pairs run never holds every pair in memory. Requires pyarrow.

Rough sizing: the per-gene detail costs ~300 bytes/row, so all pairs of 454 clusters is ~30 GB.
Use the on-disk form for anything beyond a few hundred pairs.
"""
import pickle
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(1, '/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering_allen/')
import transcriptomic_clustering as tc

# One threshold dict, mirroring scrattch.bigcat's de_param() so the same object can be carried
# through DE, merging and QC exactly as R passes de.param around. de_all_pairs consumes only the
# gene-level filters and ignores the rest.
#
#   R de_param()        Python key             consumed by de_all_pairs?
#   low.th      = 1     low_thresh             no  - used to compute `present` (get_cluster_means)
#   padj.th     = 0.01  padj_thresh            yes
#   lfc.th      = 1     lfc_thresh             yes
#   q1.th       = 0.5   q1_thresh              yes
#   q2.th       = NULL  q2_thresh              yes
#   q.diff.th   = 0.7   qdiff_thresh           yes
#   min.cells   = 10    cluster_size_thresh    yes
#   de.score.th = 150   score_thresh           no  - merge / QC decision threshold
#   min.genes   = 5     min_genes              no  - merge decision
#
# NOTE log-scale thresholds are ln-converted: scrattch.bigcat normalizes with log2(CPM+1) while this
# package uses ln(CPM+1), so low.th = 1 and lfc.th = 1 (log2) both become 0.6931472 (ln).
# This is the equivalent of: de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=10)
THRESHOLDS = {
    'low_thresh': 0.6931472,      # R low.th    = 1
    'padj_thresh': 0.01,          # R padj.th
    'lfc_thresh': 0.6931472,      # R lfc.th    = 1
    'q1_thresh': 0.5,             # R q1.th
    'q2_thresh': None,            # R q2.th
    'qdiff_thresh': 0.7,          # R q.diff.th
    'cluster_size_thresh': 10,    # R min.cells
    'score_thresh': 150,          # R de.score.th  (used by the merge/QC steps below)
    'min_genes': 5,               # R min.genes    (used by the merge/QC steps below)
}

adata = sc.read('/path/to/normalized_adata.h5ad')          # adata.X already ln(CPM+1) normalized

# ---------------------------------------------------------------- input: the clustering to test
# final_merge (like onestep_clust / iter_clust) returns a LIST OF CELL-INDEX LISTS -- one list per
# cluster, holding positional indices into `adata`. It carries no cluster names; the pipeline labels
# clusters by position, 1-based (see hicatMPI/final_merging.py).
#
#   clusters, markers = tc.final_merge(adata, cluster_assignments, ..., final_merge_kwargs=...)
#
# In a normal workflow you reload what the pipeline already wrote:
#   import pickle
#   with open('out/clusters_after_final_merge.pkl', 'rb') as f:
#       clusters = pickle.load(f)          # List[List[int]]
clusters = pickle.load(open('out/clusters_after_final_merge.pkl', 'rb'))

# de_all_pairs wants {cluster_label: cell indices}, so name the clusters by position as the pipeline does
cluster_assignments = {str(i + 1): np.asarray(cells) for i, cells in enumerate(clusters)}

# IMPORTANT: those indices are positions into the adata that produced them. Load and subset `adata`
# exactly as the clustering run did, or the indices will point at the wrong cells.
#
# If you only have the per-cell CSV (out/clusters_after_final_merge.csv) instead of the pkl:
#   labels = pd.read_csv('out/clusters_after_final_merge.csv', index_col=0).iloc[:, -1].astype(str)
#   labels = labels.loc[adata.obs_names]                       # align to adata order
#   cluster_assignments = {c: np.flatnonzero((labels == c).to_numpy()) for c in labels.unique()}

# de_all_pairs works from cluster statistics, not the cell matrix
cluster_by_obs = np.empty(adata.n_obs, dtype=object)
for label, idx in cluster_assignments.items():
    cluster_by_obs[idx] = label
means, present, variances = tc.get_cluster_means(
    adata, cluster_assignments, cluster_by_obs, low_th=THRESHOLDS['low_thresh']   # R low.th
)
cl_size = {label: len(idx) for label, idx in cluster_assignments.items()}

# ---------------------------------------------------------------- 1. in memory
summary, detail = tc.de_all_pairs(means, variances, present, cl_size, THRESHOLDS, top_n=500)
# summary: one row per pair  [pair, P1, P2, up_num, down_num, num, up_score, down_score, score]
# detail : one row per gene  [pair, P1, P2, gene, logPval, sign, rank, lfc]; sign='up' = higher in P1
print(summary.sort_values('score').head())

# a specific pair, in memory
print(detail[detail['pair'] == '1_2'].head())

# ---------------------------------------------------------------- 2. to disk (recommended at scale)
# Define the cluster bins explicitly, exactly as post_clustering_qc.R does:
#
#   cn <- as.character(sort(unique(cl)))
#   cl.bin.size = 100
#   cl.bin = data.frame(cl = cn, bin = ceiling((1:length(cn)/cl.bin.size)))
#   save(cl.bin, file = 'cl.bin.rda')
#
# Same thing here -- clusters are sorted, then chopped into consecutive groups of bin_size:
cl_bin = tc.make_cl_bin(cluster_assignments.keys(), bin_size=100)   # {'1': 1, ..., '101': 2, ...}

# Why pass cl_bin rather than just bin_size: the binning IS the on-disk layout, so pinning it keeps
# partitions identical across runs (re-runs, added clusters, a different n_jobs) and lets you reuse
# one binning for several datasets. de_all_pairs also saves it to <dir>/_cl_bin.json, so
# tc.read_de_pairs() finds it automatically -- the equivalent of R's cl.bin.rda.
#
# Each bin-pair is written to <dir>/bin.x=X/bin.y=Y/part-0.parquet, and the bin-pair is also the unit
# of parallel work, so bin_size trades off: smaller bins -> more, smaller files and more parallel
# tasks; larger bins -> fewer, larger files. Returns (None, None) -- everything is on disk.
tc.de_all_pairs(
    means, variances, present, cl_size, THRESHOLDS,
    top_n=500,                   # R de_selected_pairs default
    out_dir='de_parquet',        # per-gene detail   (R out.dir)
    summary_dir='de_summary',    # per-pair summary  (R summary.dir)
    cl_bin=cl_bin,               # R cl.bin;  omit to have it built from bin_size=100
    n_jobs=16,                   # bin-pairs computed in parallel (R mc.cores)
)

# Read a few pairs back: only the partitions holding them are opened.
sub = tc.read_de_pairs('de_parquet', pairs=[('1', '2'), ('3', '17')])
print(sub.head())

# The whole summary is small enough to load at once, and is what the QC steps consume.
summary = tc.read_de_pairs('de_summary')
low_quality = tc.find_low_quality(summary, low_th=2)
triplets = tc.find_triplets(summary)
print(f"{len(summary):,} pairs | {len(low_quality)} low-quality clusters | {len(triplets)} triplets")


# ==============================================================================================
# NOTES: reading the output, and querying specific pairs
# ==============================================================================================
#
# ---- What the two tables mean -----------------------------------------------------------------
#
# summary -- ONE ROW PER PAIR. This is what merge/QC decisions are made from.
#   pair        "P1_P2" identifier, e.g. "3_17"
#   P1, P2      the two clusters. Ordering is lexicographic on the label string (matching R's
#               `pairs[,1] <= pairs[,2]`), so pair "19_3" is normal -- "19" sorts before "3".
#   up_num      genes significantly HIGHER in P1
#   down_num    genes significantly HIGHER in P2
#   num         up_num + down_num  -> how separable the pair is
#   up_score /  sum of -log10(padj) over those genes, each capped at 20 per gene
#   down_score  (the cap matches scrattch.bigcat de_stats_pair)
#   score       up_score + down_score. LOW score = similar clusters = merge candidates.
#
# detail -- ONE ROW PER GENE PER PAIR (top_n per direction). All 8 columns:
#
#   DIRECTION LIVES IN P1/P2, exactly as scrattch.bigcat stores it: **P1 is always the cluster the
#   gene is UP in.** So a pair's rows are split into two orientation blocks -- genes up in `a` are
#   stored (P1=a, P2=b) and genes up in `b` are stored (P1=b, P2=a). Given that, `pair` and `sign`
#   are BOTH redundant; they are kept because R keeps them and they make a row self-describing.
#
#   pair        canonical "<lo>_<hi>" (lexicographic), the SAME for both orientation blocks, so it
#               remains the single grouping/filter key for a pair. REDUNDANT: derivable from P1/P2.
#   P1          the cluster this gene is HIGHER in.
#   P2          the cluster it is compared against.
#               (P1/P2 are constant within an orientation block, not within a pair.)
#   gene        gene name
#   logPval     -log10(padj), UNCAPPED (matches R, which stores the raw value). The cap of 20 is
#               applied only when summing into a de score, so
#               summary.score == sum(min(detail.logPval, 20)) over the pair's genes.
#   sign        'up' / 'down'. REDUNDANT: sign == 'up' exactly when P1 < P2 lexicographically, i.e.
#               it records whether this row's orientation is the canonical one. Filtering on
#               `P1 == c` is the direct way to get "genes up in cluster c" -- that is what R's
#               select_top_pos_markers_ds does.
#   rank        1 = most significant within that direction (restarts at 1 for 'up' and for 'down')
#   lfc         |log fold change|, ABSOLUTE VALUE -- the direction lives in `sign`, not the number.
#               Natural log here, not log2 (see the THRESHOLDS note above).
#
# So each pair contributes up to 2*top_n rows, laid out as a block of 'up' rows ranked 1..n
# followed by a block of 'down' rows ranked 1..n:
#
#   pair P1 P2 gene  logPval sign  rank      lfc
#    1_2  1  2  g91    24.31   up     1 2.682556      <- most significant gene higher in cluster 1
#    1_2  1  2 g145    20.08   up     2 2.508728
#    1_2  2  1 g108    31.72 down     1 2.869499      <- most significant gene higher in cluster 2
#    1_2  2  1 g150    22.44 down     2 2.741372         (note P1/P2 swap: P1 is the "up" cluster)
#
# ---- Reading it -------------------------------------------------------------------------------
#
#   summary = tc.read_de_pairs('de_summary')          # whole summary; small, load it all
#   detail  = tc.read_de_pairs('de_parquet')          # WARNING: whole detail; ~300 bytes/row, so
#                                                     # all pairs of 454 clusters is ~30 GB
#
#   # most-similar pairs (merge candidates):
#   summary.nsmallest(20, 'score')[['pair', 'num', 'score']]
#
#   # pairs that would merge at R's de.score.th:
#   summary[summary.score < THRESHOLDS['score_thresh']]
#
#   # asymmetric pairs -- one cluster is a "low quality" version of the other:
#   summary[(summary.up_num < 5) & (summary.down_num > 50)]
#
# ---- Querying specific pairs (the point of the binned layout) ----------------------------------
#
# Pass `pairs=` and only the partitions containing them are opened -- the rest of the dataset is
# never read. Either orientation works; ('3','17') and ('17','3') return the same rows.
#
#   one   = tc.read_de_pairs('de_parquet', pairs=[('3', '17')])       # reads <dir>/_cl_bin.json
#   one   = tc.read_de_pairs('de_parquet', pairs=[('3', '17')], cl_bin=cl_bin)   # or pass it
#   many  = tc.read_de_pairs('de_parquet', pairs=[('3', '17'), ('1', '2'), ('40', '41')])
#   cols  = tc.read_de_pairs('de_parquet', pairs=[('3', '17')], columns=['gene', 'P1', 'logPval'])
#
#   # top 10 markers separating one pair, in each direction. Filter on P1, NOT on `sign`: `sign`
#   # is relative to the canonical lexicographic order, and '17' < '3', so for this pair 'up'
#   # actually means "higher in 17". P1 always names the cluster the gene is up in.
#   d = tc.read_de_pairs('de_parquet', pairs=[('3', '17')])
#   d[d['P1'] == '3'].nsmallest(10, 'rank')['gene'].tolist()       # higher in cluster 3
#   d[d['P1'] == '17'].nsmallest(10, 'rank')['gene'].tolist()      # higher in cluster 17
#
#   # every pair involving one cluster (build the pair list yourself, then one call):
#   others = [c for c in cluster_assignments if c != '3']
#   d3 = tc.read_de_pairs('de_parquet', pairs=[('3', c) for c in others])
#
#   # ---- top-N positive markers for a cluster -------------------------------------------------
#   # This is scrattch.bigcat's select_top_pos_markers_ds: take the genes up in the cluster against
#   # every other cluster and score each by sum(1000 - rank), so a gene that ranks high against MANY
#   # clusters wins. Filtering on P1 is what makes it a one-liner.
#   up3 = d3[d3['P1'] == '3']
#   top3 = (up3.assign(s=1000 - up3['rank'])
#              .groupby('gene')['s'].sum()
#              .nlargest(3))
#
# If you already hold a detail DataFrame in memory, DON'T filter it repeatedly with a boolean mask
# (`detail[detail.pair == p]` is a full scan, ~250 ms per lookup at 5M rows). Index it once:
#
#   by_pair = {k: v for k, v in detail.groupby('pair')}   # then by_pair['3_17'] is instant
#
# The summary's `bin.x` / `bin.y` columns come from the directory names and are only bookkeeping --
# `tc.load_cl_bin('de_parquet')` returns the cluster -> bin map if you need it.

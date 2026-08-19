"""
Post-clustering QC — Python port of scrattch.bigcat's post-clustering qc functions and scripts.

Provides:
  - create_pairs / get_pairs         : cluster-pair enumeration / parsing
  - de_all_pairs                     : exhaustive all-pairs DE (summary + per-gene detail)
  - find_doublet_by_marker           : flag clusters expressing markers of >1 cell type
  - find_low_quality                 : flag clusters that are low-quality versions of another
  - find_triplets / check_triplet /
    find_doublets                    : doublet detection by the "triplet" (A+B -> C) method

DE uses the same eBayes moderated-t as the merge step (de_ebayes), so scores are consistent.
"""
from typing import Any, Dict, List, Optional, Tuple
from collections import defaultdict
import itertools
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .de_ebayes import get_linear_fit_vals, moderate_variances
from .diff_expression import get_qdiff, filter_gene_stats

# thresholds consumed by filter_gene_stats (the "de" thresholds); score_thresh/min_genes/low_thresh
# are NOT filter args and are dropped before filtering (matching merge_clusters_by_de).
_FILTER_KEYS = ('q1_thresh', 'q2_thresh', 'cluster_size_thresh', 'qdiff_thresh', 'padj_thresh', 'lfc_thresh')

SCORE_CAP = 20.0   # per-gene -log10(padj) cap (matches scrattch de_stats_pair / calc_de_score)


def create_pairs(cluster_labels) -> List[Tuple[str, str]]:
    """All nondirectional cluster pairs (P1<P2, self excluded). Mirrors R create_pairs."""
    labs = sorted({str(c) for c in cluster_labels})
    return [(a, b) for a, b in itertools.combinations(labs, 2)]


def get_pairs(pair_strs) -> pd.DataFrame:
    """Parse 'P1_P2' strings -> DataFrame(P1, P2) indexed by the string. Mirrors R get_pairs."""
    rows = [s.split('_', 1) for s in pair_strs]
    df = pd.DataFrame(rows, columns=['P1', 'P2'], index=list(pair_strs))
    return df


def _logpval(padj: np.ndarray) -> np.ndarray:
    with np.errstate(divide='ignore'):
        lp = -np.log10(padj.astype(float))
    return np.minimum(lp, SCORE_CAP)


def make_cl_bin(clusters, bin_size: int = 100) -> Dict[str, int]:
    """
    Assign clusters to bins, mirroring scrattch.bigcat's
    `cl.bin = data.frame(cl=cn, bin=ceiling((1:length(cn)/cl.bin.size)))`.

    Clusters are sorted numerically when every label is integer-like, otherwise lexicographically,
    then chopped into consecutive groups of `bin_size` (bins are 1-based, as in R).
    """
    labels = [str(c) for c in clusters]
    try:
        order = sorted(labels, key=int)
    except ValueError:
        order = sorted(labels)
    return {c: i // bin_size + 1 for i, c in enumerate(order)}


def _bin_pair(a, b, cl_bin: Dict[str, int]) -> Tuple[int, int]:
    """Bin-pair a cluster pair belongs to, ordered bin.x <= bin.y (R writes the upper triangle only)."""
    bx, by = cl_bin[str(a)], cl_bin[str(b)]
    return (bx, by) if bx <= by else (by, bx)


def _partition_path(root: str, bin_x: int, bin_y: int) -> str:
    """Hive-style partition directory, matching R's file.path(out.dir, "bin.x=X", "bin.y=Y")."""
    return os.path.join(root, f"bin.x={bin_x}", f"bin.y={bin_y}")


def _write_partition(df: pd.DataFrame, root: str, bin_x: int, bin_y: int) -> None:
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as err:
        raise ImportError(
            "Writing binned DE results requires pyarrow (pip install pyarrow). "
            "Omit out_dir/summary_dir to keep results in memory instead."
        ) from err
    d = _partition_path(root, bin_x, bin_y)
    os.makedirs(d, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   os.path.join(d, "part-0.parquet"))


def _save_cl_bin(root: str, cl_bin: Dict[str, int]) -> None:
    """Store cl_bin beside the dataset so readers can map pairs to partitions without it being passed."""
    os.makedirs(root, exist_ok=True)
    with open(os.path.join(root, "_cl_bin.json"), "w") as f:
        json.dump({str(k): int(v) for k, v in cl_bin.items()}, f)


def load_cl_bin(root: str) -> Dict[str, int]:
    """Load the cl_bin mapping saved next to a binned DE dataset."""
    with open(os.path.join(root, "_cl_bin.json")) as f:
        return json.load(f)


def read_de_pairs(root: str,
                  pairs: Optional[List[Tuple[Any, Any]]] = None,
                  cl_bin: Optional[Dict[str, int]] = None,
                  columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Read DE rows for specific cluster pairs from a binned dataset written by `de_all_pairs`.

    Only the partitions those pairs live in are opened -- the whole point of the bin layout.
    Works for both the detail (out_dir) and summary (summary_dir) datasets.

    Parameters
    ----------
    root: dataset directory written by de_all_pairs
    pairs: cluster pairs to fetch; None reads the entire dataset
    cl_bin: cluster -> bin map; read from `root/_cl_bin.json` when omitted
    columns: subset of columns to read

    Returns
    -------
    DataFrame of the matching rows (empty if none of the pairs are present)
    """
    try:
        import pyarrow.dataset as pads
    except ImportError as err:
        raise ImportError("Reading binned DE results requires pyarrow (pip install pyarrow).") from err

    if pairs is None:
        ds = pads.dataset(root, format="parquet", partitioning="hive")
        return ds.to_table(columns=columns).to_pandas()

    if cl_bin is None:
        cl_bin = load_cl_bin(root)
    pairs = [(str(a), str(b)) for a, b in pairs]
    wanted = {f"{a}_{b}" for a, b in pairs} | {f"{b}_{a}" for a, b in pairs}

    # map the requested pairs to their partitions and open only those
    parts = sorted({_partition_path(root, *_bin_pair(a, b, cl_bin)) for a, b in pairs})
    files = [os.path.join(d, f) for d in parts if os.path.isdir(d)
             for f in sorted(os.listdir(d)) if f.endswith(".parquet")]
    if not files:
        return pd.DataFrame(columns=columns or [])
    df = pads.dataset(files, format="parquet").to_table(columns=columns).to_pandas()
    return df[df["pair"].isin(wanted)].reset_index(drop=True)


def _de_one_pair(a, b, cluster_means, present_cluster_means, cl_size, genes, sqrt_sigma,
                 stdev_unscaled, df, df_prior, df_pooled, filt, padj_alpha, top_n, want_detail):
    """eBayes DE for a single cluster pair -> (summary row, detail frame or None)."""
    means_diff = (cluster_means.loc[a] - cluster_means.loc[b]).to_frame()
    stdev_comb = np.sqrt(np.sum(stdev_unscaled.loc[[a, b]] ** 2))[0]
    df_total = min(df + df_prior, df_pooled)
    t_vals = means_diff / sqrt_sigma / stdev_comb
    p_vals = 2 * stats.t.sf(np.abs(t_vals[0]), df_total)
    _, p_adj, _, _ = multipletests(p_vals, alpha=padj_alpha, method='holm')

    s = pd.DataFrame(index=genes)
    s['p_value'] = p_vals; s['p_adj'] = p_adj; s['lfc'] = means_diff.values
    s['q1'] = present_cluster_means.loc[a].values
    s['q2'] = present_cluster_means.loc[b].values
    s['qdiff'] = get_qdiff(present_cluster_means.loc[a].values, present_cluster_means.loc[b].values)

    up = filter_gene_stats(s, 'up-regulated', cl1_size=cl_size[a], cl2_size=cl_size[b], **filt)
    down = filter_gene_stats(s, 'down-regulated', cl1_size=cl_size[a], cl2_size=cl_size[b], **filt)
    up = up.sort_values('p_adj'); down = down.sort_values('p_adj')
    up_lp = _logpval(up['p_adj'].to_numpy()); down_lp = _logpval(down['p_adj'].to_numpy())

    summ_row = {
        'pair': f"{a}_{b}", 'P1': a, 'P2': b,
        'up_num': len(up), 'down_num': len(down), 'num': len(up) + len(down),
        'up_score': float(up_lp.sum()), 'down_score': float(down_lp.sum()),
        'score': float(up_lp.sum() + down_lp.sum()),
    }
    if not want_detail:
        return summ_row, None

    uh, dh = up.head(top_n), down.head(top_n)
    uh_lp, dh_lp = up_lp[:len(uh)], down_lp[:len(dh)]
    d = pd.DataFrame({
        'pair': f"{a}_{b}", 'P1': a, 'P2': b,
        'gene': list(uh.index) + list(dh.index),
        'logPval': np.concatenate([uh_lp, dh_lp]),
        'sign': ['up'] * len(uh) + ['down'] * len(dh),
        'rank': list(range(1, len(uh) + 1)) + list(range(1, len(dh) + 1)),
        'lfc': np.abs(np.concatenate([uh['lfc'].to_numpy(), dh['lfc'].to_numpy()])),
    })
    return summ_row, (d if len(d) else None)


def de_all_pairs(
        cluster_means: pd.DataFrame,     # clusters x genes (normalized)
        cluster_variances: pd.DataFrame, # clusters x genes
        present_cluster_means: pd.DataFrame,
        cl_size: Dict[Any, int],
        thresholds: Dict[str, Any],
        pairs: Optional[List[Tuple[Any, Any]]] = None,
        top_n: int = 500,
        return_detail: bool = True,
        out_dir: Optional[str] = None,
        summary_dir: Optional[str] = None,
        cl_bin: Optional[Dict[str, int]] = None,
        bin_size: int = 100,
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Compute eBayes DE for ALL (or given) cluster pairs. Mirrors scrattch.bigcat de_all_pairs
    (+ its de_summary / de_parquet outputs).

    With `out_dir` / `summary_dir`, results stream to disk in scrattch.bigcat's binned layout --
    clusters are assigned to bins of `bin_size` and each bin-pair is written to
    `<dir>/bin.x=X/bin.y=Y/part-0.parquet` -- so an all-pairs run never holds every pair in memory
    and any pair can be read back later with `read_de_pairs` by opening only its partition.
    `_cl_bin.json` is saved alongside so readers are self-contained. Requires pyarrow.

    Parameters
    ----------
    pairs: pairs to test; None enumerates all pairs of cluster_means.index
    top_n: detail genes kept per direction per pair (R de_selected_pairs default: 500)
    return_detail: build the per-gene detail (ignored for a dataset already written to out_dir)
    out_dir / summary_dir: write detail / summary as binned parquet instead of returning them
    cl_bin: cluster -> bin map; built from `bin_size` when omitted
    bin_size: clusters per bin (R cl.bin.size default: 100)

    Returns
    -------
    summary : DataFrame [pair, P1, P2, up_num, down_num, num, up_score, down_score, score],
              or None when written to `summary_dir`
    detail  : DataFrame [pair, P1, P2, gene, logPval, sign, rank, lfc], or None when written to
              `out_dir` or when return_detail=False. `sign` = 'up' means higher in P1.
    """
    if pairs is None:
        pairs = create_pairs(cluster_means.index)
    pairs = [(str(a), str(b)) for a, b in pairs]
    filt = {k: thresholds.get(k) for k in _FILTER_KEYS}

    # eBayes fit ONCE across all clusters (same as de_pairs_ebayes)
    sigma_sq, df, stdev_unscaled = get_linear_fit_vals(cluster_variances, cl_size)
    sigma_sq_post, _var_prior, df_prior = moderate_variances(sigma_sq, df)
    df_pooled = np.sum(df)
    sqrt_sigma = np.sqrt(sigma_sq_post)
    genes = cluster_means.columns
    padj_alpha = thresholds['padj_thresh']

    def _run(a, b, want_detail):
        return _de_one_pair(a, b, cluster_means, present_cluster_means, cl_size, genes, sqrt_sigma,
                            stdev_unscaled, df, df_prior, df_pooled, filt, padj_alpha,
                            top_n, want_detail)

    to_disk = out_dir is not None or summary_dir is not None
    if not to_disk:
        summ_rows, detail_frames = [], []
        for a, b in pairs:
            row, d = _run(a, b, return_detail)
            summ_rows.append(row)
            if d is not None:
                detail_frames.append(d)
        summary = pd.DataFrame(summ_rows)
        detail = pd.concat(detail_frames, ignore_index=True) if (return_detail and detail_frames) else \
            (pd.DataFrame(columns=['pair', 'P1', 'P2', 'gene', 'logPval', 'sign', 'rank', 'lfc'])
             if return_detail else None)
        return summary, detail

    # --- binned, streamed to disk: one partition per bin-pair, nothing accumulated across bins ---
    if cl_bin is None:
        cl_bin = make_cl_bin(cluster_means.index, bin_size=bin_size)
    for root in (out_dir, summary_dir):
        if root is not None:
            _save_cl_bin(root, cl_bin)

    by_bin: Dict[Tuple[int, int], List[Tuple[str, str]]] = defaultdict(list)
    for a, b in pairs:
        by_bin[_bin_pair(a, b, cl_bin)].append((a, b))

    want_detail = out_dir is not None and return_detail
    kept_summary: List[Dict[str, Any]] = []
    for (bx, by), bin_pairs in sorted(by_bin.items()):
        rows, frames = [], []
        for a, b in bin_pairs:
            row, d = _run(a, b, want_detail)
            rows.append(row)
            if d is not None:
                frames.append(d)
        if out_dir is not None and frames:
            _write_partition(pd.concat(frames, ignore_index=True), out_dir, bx, by)
        if summary_dir is not None:
            _write_partition(pd.DataFrame(rows), summary_dir, bx, by)
        else:
            kept_summary.extend(rows)
        del rows, frames

    summary = None if summary_dir is not None else pd.DataFrame(kept_summary)
    return summary, None


def find_doublet_by_marker(cluster_means: pd.DataFrame, markers: Dict[str, List[str]], th: float = 3.5) -> pd.DataFrame:
    """
    Flag clusters whose mean expression exceeds `th` for markers of >1 distinct cell type.
    Mirrors R find_doublet_by_marker.  cluster_means: clusters x genes.

    Returns the marker means (clusters x all-markers) for the flagged doublet clusters.
    """
    groups = {k: [g for g in v if g in cluster_means.columns] for k, v in markers.items()}
    groups = {k: v for k, v in groups.items() if v}
    cl_val = pd.DataFrame({k: cluster_means[v].max(axis=1) for k, v in groups.items()})  # clusters x celltypes
    is_doublet = (cl_val > th).sum(axis=1) > 1
    all_markers = [g for v in groups.values() for g in v]
    return cluster_means.loc[is_doublet, all_markers]


def find_low_quality(summary: pd.DataFrame, low_th: int = 2) -> pd.DataFrame:
    """
    From the all-pairs summary, flag clusters that are low-quality versions of another cluster:
    a pair where one side has < low_th DE genes UP means that side is a low-quality subset.
    Mirrors R find_low_quality_big. Returns [pair, P1, P2, up_num, down_num, cl, cl_low].
    """
    df = summary[(summary['up_num'] < low_th) | (summary['down_num'] < low_th)].copy()
    # if few genes are UP in P1 (up_num small) -> P1 is the low-quality one, P2 is the reference
    up_small = df['up_num'] < low_th
    df['cl'] = np.where(up_small, df['P2'], df['P1'])
    df['cl_low'] = np.where(up_small, df['P1'], df['P2'])
    return df[['pair', 'P1', 'P2', 'up_num', 'down_num', 'cl', 'cl_low']]


def find_triplets(summary: pd.DataFrame, min_up_num: int = 30, max_down_num: int = 10,
                  min_de_num: int = 50) -> pd.DataFrame:
    """
    Enumerate candidate doublet triplets (cl_up ~ doublet of cl_down_x + cl_down_y).
    Mirrors R find_triplets_big. `summary` is the de_all_pairs summary.
    """
    s = summary
    asym = s[((s.up_num > min_up_num) & (s.down_num < max_down_num) & (s.up_num - s.down_num > min_up_num)) |
             ((s.down_num > min_up_num) & (s.up_num < max_down_num) & (s.down_num - s.up_num > min_up_num))].copy()
    up_bigger = asym.up_num > asym.down_num
    asym['cl_up'] = np.where(up_bigger, asym.P1, asym.P2)
    asym['cl_down'] = np.where(up_bigger, asym.P2, asym.P1)
    asym['up_num_o'] = np.where(up_bigger, asym.up_num, asym.down_num)
    asym['down_num_o'] = np.where(up_bigger, asym.down_num, asym.up_num)
    # keep cl_up that appear as the "bigger" side in >1 asymmetric pair
    counts = asym['cl_up'].value_counts()
    asym = asym[asym['cl_up'].isin(counts[counts > 1].index)]

    # self-join on cl_up -> pair its two down-partners
    t = asym[['cl_up', 'cl_down', 'up_num_o', 'down_num_o']]
    trip = t.merge(t, on='cl_up', suffixes=('_x', '_y'))
    trip = trip[trip.cl_down_x != trip.cl_down_y].copy()
    trip['P1'] = np.minimum(trip.cl_down_x, trip.cl_down_y)
    trip['P2'] = np.maximum(trip.cl_down_x, trip.cl_down_y)
    trip['pair'] = trip.P1.astype(str) + '_' + trip.P2.astype(str)
    # the two parents must be genuinely different from each other (both directions have many DE genes)
    diff_pairs = set(s[(s.up_num > min_de_num) & (s.down_num > min_de_num)]['pair'])
    trip = trip[trip['pair'].isin(diff_pairs)]
    trip = trip.sort_values(['cl_up']).reset_index(drop=True)
    return trip


def _dir_genes(detail_by_pair: Dict[str, pd.DataFrame], c1: str, c2: str, top_n: int) -> Dict[str, float]:
    """Genes higher in c1 than c2 (logPval), rank<=top_n. Uses stored P1<P2 + sign column."""
    a, b = (c1, c2) if c1 < c2 else (c2, c1)
    want = 'up' if c1 < c2 else 'down'   # 'up' = higher in stored P1
    d = detail_by_pair.get(f"{a}_{b}")
    if d is None:
        return {}
    d = d[(d['sign'] == want) & (d['rank'] <= top_n)]
    return dict(zip(d['gene'], d['logPval']))


def check_triplet(detail_by_pair: Dict[str, pd.DataFrame], cl_up: str, cl_x: str, cl_y: str,
                  top_n: int = 50) -> Dict[str, Any]:
    """
    Score whether cl_up is a doublet of parents cl_x + cl_y. Mirrors R check_triplet_big:
    of the genes distinguishing the two parents, what fraction does cl_up "inherit" from each.
    """
    cl, c1, c2 = str(cl_up), str(cl_x), str(cl_y)
    trunc = lambda g: float(np.minimum(np.fromiter(g.values(), float, len(g)), SCORE_CAP).sum()) if g else 0.0

    up_genes = _dir_genes(detail_by_pair, c1, c2, top_n)   # higher in parent1 than parent2
    down_genes = _dir_genes(detail_by_pair, c2, c1, top_n) # higher in parent2 than parent1
    up_s, down_s = trunc(up_genes), trunc(down_genes)

    # does the doublet (cl) inherit parent1's genes? -> overlap of (cl>c2) genes with (c1>c2) genes
    cl_vs_c2 = set(_dir_genes(detail_by_pair, cl, c2, top_n).keys())
    ou1 = {g: v for g, v in up_genes.items() if g in cl_vs_c2};  ou1_s = trunc(ou1)
    cl_vs_c1 = set(_dir_genes(detail_by_pair, cl, c1, top_n).keys())
    od1 = {g: v for g, v in down_genes.items() if g in cl_vs_c1}; od1_s = trunc(od1)

    up2 = _dir_genes(detail_by_pair, c1, cl, top_n); up2_s = trunc(up2)
    ou2 = {g: v for g, v in up2.items() if g in up_genes};        ou2_s = trunc(ou2)
    dn2 = _dir_genes(detail_by_pair, c2, cl, top_n); dn2_s = trunc(dn2)
    od2 = {g: v for g, v in dn2.items() if g in down_genes};      od2_s = trunc(od2)

    denom = up_s + down_s + up2_s + dn2_s
    score = (ou1_s + od1_s + ou2_s + od2_s) / denom if denom > 0 else 0.0
    return {
        'cl': cl, 'cl1': c1, 'cl2': c2, 'up_num': len(up_genes), 'down_num': len(down_genes),
        'score': score,
        'olap_ratio_up_1': (ou1_s / up_s) if up_s else 0.0,
        'olap_ratio_down_1': (od1_s / down_s) if down_s else 0.0,
        'olap_ratio_up_2': (ou2_s / up2_s) if up2_s else 0.0,
        'olap_ratio_down_2': (od2_s / dn2_s) if dn2_s else 0.0,
        'olap_num_up_1': len(ou1), 'olap_num_down_1': len(od1),
        'olap_num_up_2': len(ou2), 'olap_num_down_2': len(od2),
    }


def find_doublets(detail: pd.DataFrame, triplets: pd.DataFrame, top_n: int = 50,
                  score_th: float = 0.8, olap_th: float = 1.6) -> pd.DataFrame:
    """
    For each candidate cl_up, test its triplets; record the first that looks like a doublet
    (score>score_th and olap_ratio_up_1+olap_ratio_down_1 > olap_th). Mirrors find_doublets_all_big.
    Returns one row per candidate that scores best (all candidates' best result).
    """
    detail_by_pair = {p: d for p, d in detail.groupby('pair')}
    results = []
    for cl_up, tg in triplets.groupby('cl_up'):
        best = None
        for _, row in tg.iterrows():
            r = check_triplet(detail_by_pair, cl_up, row['cl_down_x'], row['cl_down_y'], top_n=top_n)
            if best is None or r['score'] > best['score']:
                best = r
            if r['score'] > score_th and (r['olap_ratio_up_1'] + r['olap_ratio_down_1']) > olap_th:
                best = r
                break
        if best is not None:
            results.append(best)
    return pd.DataFrame(results)

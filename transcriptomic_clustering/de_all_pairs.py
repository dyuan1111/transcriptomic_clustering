"""
Exhaustive all-pairs differential expression — Python port of scrattch.bigcat's `de_all_pairs`.

Provides:
  - create_pairs / get_pairs         : cluster-pair enumeration / parsing
  - make_cl_bin / load_cl_bin        : the cluster -> bin map that defines the on-disk layout
  - de_all_pairs                     : exhaustive all-pairs DE (summary + per-gene detail),
                                       in memory or streamed to a binned parquet dataset
  - read_de_pairs                    : read specific pairs back, opening only the needed partitions

Uses the same eBayes moderated-t as the merge step (de_ebayes), so scores are consistent.

Detail layout follows scrattch.bigcat: direction is encoded in P1/P2 (P1 is always the cluster the
gene is UP in), `pair` keeps the canonical lexicographic label pair, and partitions are
bin.x=<X>/bin.y=<Y> covering both orientations.
"""
from typing import Any, Dict, List, Optional, Tuple
from collections import defaultdict
import itertools
import json
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .de_ebayes import get_linear_fit_vals, moderate_variances
from .diff_expression import get_qdiff, no_gc_collect

# thresholds consumed by filter_gene_stats (the "de" thresholds); score_thresh/min_genes/low_thresh
# are NOT filter args and are dropped before filtering (matching merge_clusters_by_de).
_FILTER_KEYS = ('q1_thresh', 'q2_thresh', 'cluster_size_thresh', 'qdiff_thresh', 'padj_thresh', 'lfc_thresh')

SCORE_CAP = 20.0   # per-gene -log10(padj) cap applied when SUMMING the de score (not when storing).
                   # Mirrors scrattch de.genes.R: `up.genes = -log10(padj)` is kept uncapped, and only
                   # the copy summed into up.score/down.score is truncated at 20.


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
    """-log10(padj), UNCAPPED. Cap with SCORE_CAP only when summing into a de score."""
    with np.errstate(divide='ignore'):
        return -np.log10(padj.astype(float))


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

    # Map the requested pairs to their partitions and open only those. A pair's genes are split over
    # both orientations of its bin-pair (P1/P2 encode direction), so read (bx, by) and (by, bx).
    wanted_bins = set()
    for a, b in pairs:
        bx, by = _bin_pair(a, b, cl_bin)
        wanted_bins.add((bx, by))
        wanted_bins.add((by, bx))
    parts = sorted({_partition_path(root, bx, by) for bx, by in wanted_bins})
    files = [os.path.join(d, f) for d in parts if os.path.isdir(d)
             for f in sorted(os.listdir(d)) if f.endswith(".parquet")]
    if not files:
        return pd.DataFrame(columns=columns or [])
    df = pads.dataset(files, format="parquet").to_table(columns=columns).to_pandas()
    return df[df["pair"].isin(wanted)].reset_index(drop=True)


# Per-worker context for the parallel bin-pair loop. Populated once per process by
# _init_de_worker so the (clusters x genes) frames are shipped once per worker, not once per task.
_DE_CTX: Dict[str, Any] = {}


def _init_de_worker(ctx: Dict[str, Any]) -> None:
    _DE_CTX.clear()
    _DE_CTX.update(ctx)


def _process_bin_pair(task):
    """
    Compute DE for every pair in one bin-pair and write its partition(s).

    Returns (bin_x, bin_y, summary_rows) with summary_rows empty when the summary was written to
    disk -- so nothing large travels back to the parent.
    """
    (bin_x, bin_y), bin_pairs = task
    c = _DE_CTX
    rows, frames = [], []
    for a, b in bin_pairs:
        row, detail = _de_one_pair(
            a, b, c['cluster_means'], c['present_cluster_means'], c['cl_size'], c['genes'],
            c['sqrt_sigma'], c['stdev_unscaled'], c['df'], c['df_prior'], c['df_pooled'],
            c['filt'], c['padj_alpha'], c['top_n'], c['want_detail'],
        )
        rows.append(row)
        if detail is not None:
            frames.append(detail)
    if c['out_dir'] is not None and frames:
        # A pair's genes go to two partitions, because P1/P2 carry the direction: the genes up in `a`
        # land in (bin(a), bin(b)) and the genes up in `b` in (bin(b), bin(a)). Since tasks are keyed
        # by the NORMALISED bin-pair, this task is the only writer of either directory -- no collision.
        detail = pd.concat(frames, ignore_index=True)
        cb = c['cl_bin']
        bx = detail['P1'].map(lambda x: cb[str(x)])
        by = detail['P2'].map(lambda x: cb[str(x)])
        for (pbx, pby), grp in detail.groupby([bx, by], sort=True):
            _write_partition(grp.reset_index(drop=True), c['out_dir'], int(pbx), int(pby))
    if c['summary_dir'] is not None:
        _write_partition(pd.DataFrame(rows), c['summary_dir'], bin_x, bin_y)
        rows = []
    return bin_x, bin_y, rows


def _process_pair_chunk(chunk):
    """Compute DE for a chunk of pairs and return (summary rows, detail frame or None).

    Used by the in-memory path, where results come back to the parent rather than being written.
    """
    c = _DE_CTX
    rows, frames = [], []
    for a, b in chunk:
        row, detail = _de_one_pair(
            a, b, c['cluster_means'], c['present_cluster_means'], c['cl_size'], c['genes'],
            c['sqrt_sigma'], c['stdev_unscaled'], c['df'], c['df_prior'], c['df_pooled'],
            c['filt'], c['padj_alpha'], c['top_n'], c['want_detail'],
        )
        rows.append(row)
        if detail is not None:
            frames.append(detail)
    return rows, (pd.concat(frames, ignore_index=True) if frames else None)


def _de_one_pair(a, b, cluster_means, present_cluster_means, cl_size, genes, sqrt_sigma,
                 stdev_unscaled, df, df_prior, df_pooled, filt, padj_alpha, top_n, want_detail):
    """
    eBayes DE for a single cluster pair -> (summary row, detail frame or None).

    Uses numpy boolean masks rather than building a per-pair DataFrame and calling
    filter_gene_stats; the filter logic is replicated exactly (see diff_expression.filter_gene_stats:
    'up' compares q1 against q1_thresh / cluster_size and q2 against q2_thresh, 'down' the reverse).
    """
    means_diff = (cluster_means.loc[a] - cluster_means.loc[b]).to_frame()
    stdev_comb = np.sqrt(np.sum(stdev_unscaled.loc[[a, b]] ** 2))[0]
    df_total = min(df + df_prior, df_pooled)
    t_vals = means_diff / sqrt_sigma / stdev_comb
    p_vals = 2 * stats.t.sf(np.abs(t_vals[0]), df_total)
    # statsmodels calls gc.collect() on every multipletests() call; across an all-pairs run that is
    # the dominant cost (~90% of runtime), so neutralize it here as the merge path does.
    with no_gc_collect():
        _, p_adj, _, _ = multipletests(p_vals, alpha=padj_alpha, method='holm')

    lfc = np.asarray(means_diff.values).ravel()
    q1 = np.asarray(present_cluster_means.loc[a].values)
    q2 = np.asarray(present_cluster_means.loc[b].values)
    qdiff = get_qdiff(q1, q2)

    q1_th, q2_th = filt.get('q1_thresh'), filt.get('q2_thresh')
    cs_th, qd_th = filt.get('cluster_size_thresh'), filt.get('qdiff_thresh')
    padj_th, lfc_th = filt.get('padj_thresh'), filt.get('lfc_thresh')
    size_a, size_b = cl_size[a], cl_size[b]

    up, down = lfc > 0, lfc < 0
    if padj_th:
        m = p_adj < padj_th; up &= m; down &= m
    if lfc_th:
        m = np.abs(lfc) > lfc_th; up &= m; down &= m
    if q1_th:
        up &= q1 > q1_th; down &= q2 > q1_th
    if size_a:
        up &= q1 * size_a >= cs_th
    if size_b:
        down &= q2 * size_b >= cs_th
    if q2_th:
        up &= q2 < q2_th; down &= q1 < q2_th
    if qd_th:
        m = np.abs(qdiff) > qd_th; up &= m; down &= m

    # rank each direction by adjusted p-value (same ordering as the previous sort_values('p_adj'))
    up_idx = np.flatnonzero(up); down_idx = np.flatnonzero(down)
    up_idx = up_idx[np.argsort(p_adj[up_idx], kind='quicksort')]
    down_idx = down_idx[np.argsort(p_adj[down_idx], kind='quicksort')]
    # stored per-gene values are uncapped (as R does); the cap applies only to the score sums
    up_lp, down_lp = _logpval(p_adj[up_idx]), _logpval(p_adj[down_idx])
    up_score = float(np.minimum(up_lp, SCORE_CAP).sum())
    down_score = float(np.minimum(down_lp, SCORE_CAP).sum())

    summ_row = {
        'pair': f"{a}_{b}", 'P1': a, 'P2': b,
        'up_num': len(up_idx), 'down_num': len(down_idx), 'num': len(up_idx) + len(down_idx),
        'up_score': up_score, 'down_score': down_score,
        'score': up_score + down_score,
    }
    if not want_detail:
        return summ_row, None

    uh, dh = up_idx[:top_n], down_idx[:top_n]
    if len(uh) + len(dh) == 0:
        return summ_row, None
    gene_arr = np.asarray(genes)
    # Direction is encoded in P1/P2, as scrattch.bigcat does: P1 is always the cluster the gene is
    # UP in, so an 'up' gene is stored (P1=a, P2=b) and a 'down' gene (P1=b, P2=a). `pair` keeps the
    # canonical "a_b" for both, and `sign` is retained for readability though it is redundant
    # (sign == 'up' exactly when P1 < P2 lexicographically).
    d = pd.DataFrame({
        'pair': f"{a}_{b}",
        'P1': [a] * len(uh) + [b] * len(dh),
        'P2': [b] * len(uh) + [a] * len(dh),
        'gene': np.concatenate([gene_arr[uh], gene_arr[dh]]),
        'logPval': np.concatenate([up_lp[:len(uh)], down_lp[:len(dh)]]),
        'sign': ['up'] * len(uh) + ['down'] * len(dh),
        'rank': list(range(1, len(uh) + 1)) + list(range(1, len(dh) + 1)),
        'lfc': np.abs(np.concatenate([lfc[uh], lfc[dh]])),
    })
    return summ_row, d


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
        n_jobs: int = 1,
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
    bin_size: clusters per bin (R cl.bin.size default: 100). Independent of n_jobs, so the layout of a
              written dataset is reproducible; smaller bins give more (smaller) partitions and more
              parallel tasks.
    n_jobs: worker processes (R mc.cores). On the on-disk path the unit of work is the bin-pair (as in
            R); in memory the pairs are chunked instead, so parallelism does not require writing files.

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

    def _make_ctx(want_detail, out=None, summ=None, cl_bin_ctx=None):
        return {
            'cluster_means': cluster_means, 'present_cluster_means': present_cluster_means,
            'cl_size': cl_size, 'genes': genes, 'sqrt_sigma': sqrt_sigma,
            'stdev_unscaled': stdev_unscaled, 'df': df, 'df_prior': df_prior,
            'df_pooled': df_pooled, 'filt': filt, 'padj_alpha': padj_alpha, 'top_n': top_n,
            'want_detail': want_detail, 'out_dir': out, 'summary_dir': summ,
            'cl_bin': cl_bin_ctx,
        }

    to_disk = out_dir is not None or summary_dir is not None
    if not to_disk:
        summ_rows, detail_frames = [], []
        if n_jobs > 1 and len(pairs) > 1:
            # chunk the pairs so each worker gets several; imap keeps the output in pair order
            n_chunks = min(len(pairs), max(1, n_jobs) * 4)
            size = int(np.ceil(len(pairs) / n_chunks))
            chunks = [pairs[i:i + size] for i in range(0, len(pairs), size)]
            with mp.Pool(processes=min(n_jobs, len(chunks)),
                         initializer=_init_de_worker,
                         initargs=(_make_ctx(return_detail),)) as pool:
                for rows, d in pool.imap(_process_pair_chunk, chunks):
                    summ_rows.extend(rows)
                    if d is not None:
                        detail_frames.append(d)
        else:
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
        # NOTE: R shrinks bins with mc.cores (`cl.bin.size = min(100, length(cn)/mc.cores)`), which makes
        # the on-disk layout depend on the core count. We keep bin_size fixed so a dataset is reproducible
        # and queryable regardless of n_jobs; parallelism comes from the number of bin-pairs instead.
        cl_bin = make_cl_bin(cluster_means.index, bin_size=bin_size)
    for root in (out_dir, summary_dir):
        if root is not None:
            _save_cl_bin(root, cl_bin)

    by_bin: Dict[Tuple[int, int], List[Tuple[str, str]]] = defaultdict(list)
    for a, b in pairs:
        by_bin[_bin_pair(a, b, cl_bin)].append((a, b))
    tasks = sorted(by_bin.items())

    want_detail = out_dir is not None and return_detail
    ctx = _make_ctx(want_detail, out_dir, summary_dir, cl_bin)

    kept_summary: List[Dict[str, Any]] = []
    if n_jobs > 1 and len(tasks) > 1:
        # one task per bin-pair, exactly as R parallelizes its foreach(bin1) %:% foreach(bin2) loop
        with mp.Pool(processes=min(n_jobs, len(tasks)),
                     initializer=_init_de_worker, initargs=(ctx,)) as pool:
            for _bx, _by, rows in pool.imap_unordered(_process_bin_pair, tasks):
                kept_summary.extend(rows)
    else:
        _init_de_worker(ctx)
        try:
            for task in tasks:
                _bx, _by, rows = _process_bin_pair(task)
                kept_summary.extend(rows)
        finally:
            _DE_CTX.clear()

    summary = None if summary_dir is not None else pd.DataFrame(kept_summary).sort_values('pair').reset_index(drop=True)
    return summary, None


# ---------------------------------------------------------------------------------------------
# Marker queries over a written de_parquet dataset.
#
# Port of scrattch.bigcat's markers.parquet.R query family. All of them reduce to one idea
# (R get_gene_score_ds): score a gene by how consistently it is a TOP marker across a set of
# comparisons,
#
#     score(gene) = sum over selected pairs of (max_num - rank)          # R max.num = 1000
#
# then keep the highest-scoring genes. Because `max_num - rank` is between max_num-top_n and
# max_num-1 for every stored gene, the score is dominated by *how many* comparisons a gene wins,
# with its rank acting as the tie-break.
#
# Direction: the detail encodes it in P1 (P1 is always the cluster the gene is UP in), so
# "genes up in g" is simply the rows with P1 == g -- no `sign` filtering, exactly as R does.
# ---------------------------------------------------------------------------------------------

MAX_NUM = 1000   # R get_gene_score_ds(max.num=1000)


def _detail_files(root: str, cl_bin: Optional[Dict[str, int]] = None,
                  p1_clusters: Optional[List[str]] = None) -> List[str]:
    """
    Parquet files that can hold rows for the requested P1 clusters.

    Rows describing genes up in cluster g carry P1 == g, hence bin.x == cl_bin[g], so only the
    bin.x partitions of those clusters need to be opened. Passing p1_clusters=None reads all.
    """
    import glob
    if p1_clusters is None:
        return sorted(glob.glob(os.path.join(root, "bin.x=*", "bin.y=*", "*.parquet")))
    if cl_bin is None:
        cl_bin = load_cl_bin(root)
    bins = sorted({cl_bin[str(c)] for c in p1_clusters})
    files: List[str] = []
    for b in bins:
        files += sorted(glob.glob(os.path.join(root, f"bin.x={b}", "bin.y=*", "*.parquet")))
    return files


def _rank_sum_by_cluster(root: str,
                         p1_clusters: List[str],
                         p2_allowed: Optional[Dict[str, set]] = None,
                         genes: Optional[set] = None,
                         cl_bin: Optional[Dict[str, int]] = None,
                         max_num: int = MAX_NUM) -> pd.DataFrame:
    """
    Accumulate sum(max_num - rank) per (P1, gene) over the requested rows, streaming one partition
    file at a time so memory stays bounded on datasets with tens of millions of rows.

    p2_allowed: optional {P1: set(P2)} restricting which comparisons count for each cluster.
    """
    p1_set = {str(c) for c in p1_clusters}
    parts = []
    for f in _detail_files(root, cl_bin, p1_clusters):
        t = pd.read_parquet(f, columns=['P1', 'P2', 'gene', 'rank'])
        if t.empty:
            continue
        t['P1'] = t['P1'].astype(str)
        t = t[t['P1'].isin(p1_set)]
        if genes is not None and not t.empty:
            t = t[t['gene'].isin(genes)]
        if t.empty:
            continue
        t['rank'] = t['rank'].astype(np.int64)
        t = t[t['rank'] < max_num]                      # R: filter(rank < max.num)
        if p2_allowed is not None and not t.empty:
            t['P2'] = t['P2'].astype(str)
            keep = [p2 in p2_allowed.get(p1, ()) for p1, p2 in zip(t['P1'], t['P2'])]
            t = t[keep]
        if t.empty:
            continue
        t = t.assign(s=max_num - t['rank'])
        parts.append(t.groupby(['P1', 'gene'], sort=False)['s'].sum().reset_index())
    if not parts:
        return pd.DataFrame(columns=['P1', 'gene', 'score'])
    out = pd.concat(parts, ignore_index=True)
    out = out.groupby(['P1', 'gene'], sort=False)['s'].sum().reset_index()
    return out.rename(columns={'s': 'score'})


def _rank_sum_from_frame(de: pd.DataFrame, allowed: Dict[str, set],
                         genes: Optional[set], max_num: int) -> pd.DataFrame:
    """sum(max_num - rank) per (P1, gene) over an in-memory detail frame."""
    if de.empty:
        return pd.DataFrame(columns=['P1', 'gene', 'score'])
    t = de
    t = t[t['P1'].astype(str).isin(allowed)]
    if genes is not None:
        t = t[t['gene'].isin(genes)]
    t = t[t['rank'].astype(np.int64) < max_num]
    if t.empty:
        return pd.DataFrame(columns=['P1', 'gene', 'score'])
    keep = [str(p2) in allowed.get(str(p1), ()) for p1, p2 in zip(t['P1'], t['P2'])]
    t = t[keep]
    if t.empty:
        return pd.DataFrame(columns=['P1', 'gene', 'score'])
    t = t.assign(P1=t['P1'].astype(str), s=max_num - t['rank'].astype(np.int64))
    return t.groupby(['P1', 'gene'], sort=False)['s'].sum().reset_index().rename(columns={'s': 'score'})


def get_gene_score_ds(root: Optional[str],
                      pairs,
                      genes: Optional[set] = None,
                      cl_bin: Optional[Dict[str, int]] = None,
                      max_num: int = MAX_NUM,
                      de: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Rank-sum gene score over a set of DIRECTED comparisons. Mirrors R `get_gene_score_ds`.

    Parameters
    ----------
    root: de_parquet dataset written by de_all_pairs (may be None when `de` is supplied)
    pairs: iterable of (P1, P2) -- each means "genes UP in P1, relative to P2"
    genes: optional gene universe to restrict to (R `genes`)
    cl_bin: cluster -> bin map; read from root/_cl_bin.json when omitted
    max_num: rank ceiling, R's max.num (default 1000)
    de: optional pre-loaded detail frame [P1, P2, gene, rank] to score instead of reading from
        disk -- R's `de=` argument. Much faster when the same rows are scored repeatedly.

    Returns
    -------
    DataFrame [gene, score] sorted by descending score
    """
    pairs = [(str(a), str(b)) for a, b in pairs]
    allowed: Dict[str, set] = defaultdict(set)
    for a, b in pairs:
        allowed[a].add(b)
    if de is not None:
        df = _rank_sum_from_frame(de, allowed, genes, max_num)
    else:
        df = _rank_sum_by_cluster(root, list(allowed), allowed, genes, cl_bin, max_num)
    if df.empty:
        return pd.DataFrame(columns=['gene', 'score'])
    out = df.groupby('gene', sort=False)['score'].sum().reset_index()
    return out.sort_values(['score', 'gene'], ascending=[False, True]).reset_index(drop=True)


def select_markers_pair_group_top_ds(root: str,
                                     group1,
                                     group2,
                                     genes: Optional[set] = None,
                                     cl_bin: Optional[Dict[str, int]] = None,
                                     select_sign=('up',),
                                     n_markers: int = 20,
                                     max_num: int = MAX_NUM) -> Dict[str, List[str]]:
    """
    Top markers separating one group of clusters from another. Mirrors R
    `select_markers_pair_group_top_ds`.

    'up' scores genes higher in group1 than group2 (every g1 x g2 comparison); 'down' is the
    reverse direction. Returns {'up_genes': [...], 'down_genes': [...]} -- the key omitted from
    select_sign comes back empty.
    """
    g1 = [str(g) for g in (group1 if isinstance(group1, (list, tuple, set, pd.Index)) else [group1])]
    g2 = [str(g) for g in (group2 if isinstance(group2, (list, tuple, set, pd.Index)) else [group2])]
    result: Dict[str, List[str]] = {'up_genes': [], 'down_genes': []}
    if 'up' in select_sign:
        sc = get_gene_score_ds(root, [(a, b) for a in g1 for b in g2], genes, cl_bin, max_num)
        result['up_genes'] = list(sc['gene'].head(n_markers))
    if 'down' in select_sign:
        sc = get_gene_score_ds(root, [(a, b) for a in g2 for b in g1], genes, cl_bin, max_num)
        result['down_genes'] = list(sc['gene'].head(n_markers))
    return result


def select_top_pos_markers_ds(root: str,
                              clusters,
                              select_cl=None,
                              genes: Optional[set] = None,
                              cl_bin: Optional[Dict[str, int]] = None,
                              n_markers: int = 3,
                              max_num: int = MAX_NUM) -> Dict[str, List[str]]:
    """
    Top positive markers for each cluster against all the others. Mirrors R
    `select_top_pos_markers_ds`.

    For each cluster g in select_cl, scores the genes up in g versus every other cluster in
    `clusters` and keeps the best `n_markers`.

    Unlike R -- which loops over clusters and issues one query each -- this makes a single pass
    over the needed partitions and scores every cluster at once, which is what makes it fast
    enough to run over a whole taxonomy.

    Returns {cluster: [gene, ...]}, at most n_markers per cluster, best first.
    """
    clusters = [str(c) for c in clusters]
    select_cl = clusters if select_cl is None else [str(c) for c in select_cl]
    others = set(clusters)
    allowed = {g: others - {g} for g in select_cl}

    df = _rank_sum_by_cluster(root, select_cl, allowed, genes, cl_bin, max_num)
    if df.empty:
        return {c: [] for c in select_cl}
    # R: arrange(-score) then head(n). Gene name breaks ties so repeated runs agree.
    df = df.sort_values(['P1', 'score', 'gene'], ascending=[True, False, True])
    top = df.groupby('P1', sort=False).head(n_markers)
    out = {c: [] for c in select_cl}
    for c, g in zip(top['P1'], top['gene']):
        out[c].append(g)
    return out


# ---------------------------------------------------------------------------------------------
# Greedy ("combo") marker selection -- R select_pos_markers_ds and its helpers.
#
# select_top_pos_markers_ds above answers "what are the best N markers for this cluster overall".
# This family answers a harder question: "give me a SET of markers such that every pairwise
# comparison is covered by at least `n_markers` of them". It starts from the top-scoring markers,
# counts how many comparisons each pair still needs, and then greedily adds whichever remaining
# gene covers the most uncovered pairs, until every pair is satisfied or the budget runs out.
# ---------------------------------------------------------------------------------------------


def _load_directed_detail(root: str, pairs, genes: Optional[set], cl_bin: Optional[Dict[str, int]],
                          max_num: int = MAX_NUM) -> pd.DataFrame:
    """Detail rows [P1, P2, gene, rank] for the given DIRECTED pairs, loaded once for reuse."""
    pairs = [(str(a), str(b)) for a, b in pairs]
    allowed: Dict[str, set] = defaultdict(set)
    for a, b in pairs:
        allowed[a].add(b)
    frames = []
    for f in _detail_files(root, cl_bin, list(allowed)):
        t = pd.read_parquet(f, columns=['P1', 'P2', 'gene', 'rank'])
        if t.empty:
            continue
        t['P1'] = t['P1'].astype(str); t['P2'] = t['P2'].astype(str)
        t = t[t['P1'].isin(allowed)]
        if genes is not None and not t.empty:
            t = t[t['gene'].isin(genes)]
        if t.empty:
            continue
        t = t[t['rank'].astype(np.int64) < max_num]
        if t.empty:
            continue
        keep = [p2 in allowed.get(p1, ()) for p1, p2 in zip(t['P1'], t['P2'])]
        t = t[keep]
        if not t.empty:
            frames.append(t)
    if not frames:
        return pd.DataFrame(columns=['P1', 'P2', 'gene', 'rank'])
    return pd.concat(frames, ignore_index=True)


def check_pairs_lfc(to_add: pd.DataFrame, genes, cl_means: pd.DataFrame,
                    lfc_th: float = 2.0) -> pd.DataFrame:
    """
    How many of `genes` separate each pair by more than `lfc_th`, from cluster means alone.
    Mirrors R `check_pairs_lfc`: counts genes with `cl_means[g, P1] - cl_means[g, P2] > lfc_th`.

    Returns `to_add` with an added integer `checked` column.
    """
    genes = [g for g in ([genes] if isinstance(genes, str) else genes) if g in cl_means.columns]
    out = to_add.copy()
    out['P1'] = out['P1'].astype(str); out['P2'] = out['P2'].astype(str)
    if not genes:
        out['checked'] = 0
        return out
    m = cl_means.loc[:, genes]
    d1 = m.reindex(out['P1']).to_numpy(float)
    d2 = m.reindex(out['P2']).to_numpy(float)
    out['checked'] = ((d1 - d2) > lfc_th).sum(axis=1).astype(int)
    return out


def check_pairs_ds(root: Optional[str],
                   to_add: pd.DataFrame,
                   genes,
                   cl_bin: Optional[Dict[str, int]] = None,
                   de: Optional[pd.DataFrame] = None,
                   max_num: int = MAX_NUM) -> pd.DataFrame:
    """
    How many of `genes` are stored DE markers for each pair (rank < max_num), from the dataset.
    Mirrors R `check_pairs_ds`. Returns DataFrame [P1, P2, checked].
    """
    genes = {genes} if isinstance(genes, str) else set(genes)
    pairs = list(zip(to_add['P1'].astype(str), to_add['P2'].astype(str)))
    src = de if de is not None else _load_directed_detail(root, pairs, genes, cl_bin, max_num)
    base = pd.DataFrame({'P1': [p[0] for p in pairs], 'P2': [p[1] for p in pairs]})
    if src.empty:
        base['checked'] = 0
        return base
    t = src[src['gene'].isin(genes) & (src['rank'].astype(np.int64) < max_num)]
    if t.empty:
        base['checked'] = 0
        return base
    cnt = (t.assign(P1=t['P1'].astype(str), P2=t['P2'].astype(str))
             .groupby(['P1', 'P2'], sort=False).size().rename('checked').reset_index())
    out = base.merge(cnt, on=['P1', 'P2'], how='left')
    out['checked'] = out['checked'].fillna(0).astype(int)
    return out


def select_markers_pair_direction_ds(root: Optional[str],
                                     add_num: pd.DataFrame,
                                     genes,
                                     cl_bin: Optional[Dict[str, int]] = None,
                                     de: Optional[pd.DataFrame] = None,
                                     max_genes: int = 1000,
                                     cl_means: Optional[pd.DataFrame] = None,
                                     lfc_th: float = 2.0,
                                     max_num: int = MAX_NUM) -> Dict[str, Any]:
    """
    Greedily add markers until every pair in `add_num` has its `num` requirement met.
    Mirrors R `select_markers_pair_direction_ds`.

    add_num: DataFrame [P1, P2, num] -- how many more markers each directed pair still needs.
    Returns {'select_genes': [...], 'de': <remaining detail frame or None>}.
    """
    genes = set(genes)
    add_num = add_num.copy()
    add_num['P1'] = add_num['P1'].astype(str); add_num['P2'] = add_num['P2'].astype(str)
    select_genes: List[str] = []

    def score(pairs_df, gene_pool):
        return get_gene_score_ds(root, list(zip(pairs_df['P1'], pairs_df['P2'])),
                                 genes=gene_pool, cl_bin=cl_bin, max_num=max_num, de=de)

    gene_score = score(add_num, genes)
    while len(add_num) > 0 and genes and len(select_genes) < max_genes:
        if gene_score is None or gene_score.empty:
            break
        g = str(gene_score['gene'].iloc[0])
        if cl_means is not None:
            new_checked = check_pairs_lfc(add_num[['P1', 'P2']], g, cl_means, lfc_th)
        else:
            new_checked = check_pairs_ds(root, add_num[['P1', 'P2']], g, cl_bin, de, max_num)
        if new_checked is None or new_checked.empty:
            break
        add_num = add_num.drop(columns=['checked'], errors='ignore').merge(
            new_checked[['P1', 'P2', 'checked']], on=['P1', 'P2'], how='left')
        add_num['checked'] = add_num['checked'].fillna(0).astype(int)
        add_num['num'] = add_num['num'] - add_num['checked']
        to_remove = add_num.loc[add_num['num'] <= 0, ['P1', 'P2']]
        add_num = add_num[add_num['num'] > 0].drop(columns=['checked'])
        genes = genes - {g}
        select_genes.append(g)
        if len(add_num) == 0:
            break
        # R recomputes from scratch when most pairs dropped out, otherwise subtracts the score
        # contributed by the now-satisfied pairs. Both give the same ordering; the second is cheaper.
        if len(to_remove) > len(add_num) or len(add_num) < 1000:
            gene_score = score(add_num, genes)
        else:
            rm = get_gene_score_ds(root, list(zip(to_remove['P1'], to_remove['P2'])),
                                   genes=genes, cl_bin=cl_bin, max_num=max_num, de=de)
            if rm is None or rm.empty:
                gene_score = gene_score[gene_score['gene'] != g]
                continue
            tmp = (gene_score[gene_score['gene'] != g]
                   .merge(rm.rename(columns={'score': 'rm_score'}), on='gene', how='left'))
            tmp['rm_score'] = tmp['rm_score'].fillna(0)
            tmp['score'] = tmp['score'] - tmp['rm_score']
            gene_score = (tmp[tmp['score'] > 0][['gene', 'score']]
                          .sort_values(['score', 'gene'], ascending=[False, True])
                          .reset_index(drop=True))
        if de is not None:
            de = de[de['gene'] != g]
            if de.empty:
                break
    return {'select_genes': select_genes, 'de': de}


def select_markers_pair_group_ds(root: str,
                                 group1,
                                 group2,
                                 genes: Optional[set] = None,
                                 cl_bin: Optional[Dict[str, int]] = None,
                                 n_markers: int = 20,
                                 select_sign=('up', 'down'),
                                 max_genes: int = 50,
                                 default_markers: Optional[List[str]] = None,
                                 cl_means: Optional[pd.DataFrame] = None,
                                 lfc_th: float = 2.0,
                                 de: Optional[pd.DataFrame] = None,
                                 max_num: int = MAX_NUM) -> List[str]:
    """
    Markers covering every comparison between two cluster groups. Mirrors R
    `select_markers_pair_group_ds`.

    Starts from the top `n_markers` per direction, then greedily adds genes until every directed
    pair is covered `n_markers` times or `max_genes` is reached. Returns the marker list.
    """
    g1 = [str(g) for g in (group1 if isinstance(group1, (list, tuple, set, pd.Index)) else [group1])]
    g2 = [str(g) for g in (group2 if isinstance(group2, (list, tuple, set, pd.Index)) else [group2])]

    up_pairs = [(a, b) for a in g1 for b in g2] if 'up' in select_sign else []
    down_pairs = [(b, a) for a in g1 for b in g2] if 'down' in select_sign else []
    if default_markers is None:
        top = select_markers_pair_group_top_ds(root, g1, g2, genes, cl_bin,
                                               select_sign, n_markers, max_num)
        markers = list(top['up_genes']) + list(top['down_genes'])
    else:
        markers = list(default_markers)
    if not markers:
        return markers

    add_num = pd.DataFrame(up_pairs + down_pairs, columns=['P1', 'P2'])
    add_num['num'] = n_markers
    if cl_means is not None:
        checked = check_pairs_lfc(add_num[['P1', 'P2']], markers, cl_means, lfc_th)
    else:
        checked = check_pairs_ds(root, add_num[['P1', 'P2']], markers, cl_bin, de, max_num)
    add_num = add_num.merge(checked[['P1', 'P2', 'checked']], on=['P1', 'P2'], how='left')
    add_num['checked'] = add_num['checked'].fillna(0).astype(int)
    add_num['num'] = add_num['num'] - add_num['checked']
    add_num = add_num[add_num['num'] > 0].drop(columns=['checked'])

    budget = max_genes - len(markers)
    pool = (set(genes) if genes is not None else None)
    if len(add_num) > 0 and budget > 0 and (pool is None or len(pool) > 1):
        remaining = (pool - set(markers)) if pool is not None else None
        if remaining is None:
            remaining = set(_load_directed_detail(
                root, list(zip(add_num['P1'], add_num['P2'])), None, cl_bin, max_num)['gene'])
            remaining -= set(markers)
        more = select_markers_pair_direction_ds(root, add_num, remaining, cl_bin, de, budget,
                                                cl_means, lfc_th, max_num)
        markers = markers + more['select_genes']
    return markers


def select_pos_markers_ds(root: str,
                          clusters,
                          select_cl=None,
                          genes: Optional[set] = None,
                          cl_bin: Optional[Dict[str, int]] = None,
                          n_markers: int = 1,
                          max_genes: int = 50,
                          cl_means: Optional[pd.DataFrame] = None,
                          lfc_th: float = 2.0,
                          out_dir: Optional[str] = None,
                          overwrite: bool = True,
                          max_num: int = MAX_NUM) -> Dict[str, List[str]]:
    """
    Positive markers for each cluster, covering every comparison against the other clusters.
    Mirrors R `select_pos_markers_ds`.

    Where `select_top_pos_markers_ds` just takes the N best-scoring genes per cluster, this keeps
    adding genes until each cluster-vs-other comparison is covered `n_markers` times, so a cluster
    that is hard to separate from one particular neighbour gets extra markers for that comparison.

    out_dir/overwrite mirror R's per-cluster caching (R writes <cl>.markers.rda; here JSON), so a
    long run can be resumed without recomputing clusters that already finished.

    Returns {cluster: [gene, ...]}.
    """
    clusters = [str(c) for c in clusters]
    select_cl = clusters if select_cl is None else [str(c) for c in select_cl]
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)

    result: Dict[str, List[str]] = {}
    for cl in select_cl:
        cache = os.path.join(out_dir, f"{str(cl).replace('/', '')}.markers.json") if out_dir else None
        if cache and os.path.exists(cache) and not overwrite:
            with open(cache) as fh:
                result[cl] = json.load(fh)
            continue
        others = [c for c in clusters if c != cl]
        # load this cluster's rows once and reuse them through the greedy loop
        de = _load_directed_detail(root, [(cl, o) for o in others], genes, cl_bin, max_num)
        markers = select_markers_pair_group_ds(
            root, [cl], others, genes=genes, cl_bin=cl_bin, n_markers=n_markers,
            select_sign=('up',), max_genes=max_genes, cl_means=cl_means, lfc_th=lfc_th,
            de=de, max_num=max_num,
        )
        result[cl] = markers
        if cache:
            with open(cache, "w") as fh:
                json.dump(markers, fh)
    return result

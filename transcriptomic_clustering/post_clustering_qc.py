"""
Post-clustering QC — Python port of scrattch.bigcat's post-clustering qc functions and scripts.

Consumes the output of `de_all_pairs` (see de_all_pairs.py) to flag suspect clusters:
  - find_doublet_by_marker           : flag clusters expressing markers of >1 cell type
  - find_low_quality                 : flag clusters that are low-quality versions of another
  - find_triplets / check_triplet /
    find_doublets                    : doublet detection by the "triplet" (A+B -> C) method
"""
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .de_all_pairs import SCORE_CAP

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
                  min_de_num: int = 50, all_pairs=None, select_cl=None) -> pd.DataFrame:
    """
    Enumerate candidate doublet triplets (cl_up ~ doublet of cl_down_x + cl_down_y).
    Mirrors R `find_triplets_big`.

    A doublet looks *asymmetric* against each of its parents: many genes up in it, almost none down.
    So keep the strongly asymmetric pairs, orient each so `cl_up` is the richer side, keep clusters
    that are the richer side against more than one partner, and pair those partners up. The two
    parents must themselves be well separated (`min_de_num` DE genes in both directions), otherwise
    the "doublet" is just one cluster split in two.

    Parameters
    ----------
    summary: the de_all_pairs summary
    min_up_num / max_down_num: asymmetry thresholds (R min.up.num / max.down.num)
    min_de_num: how separated the two parents must be (R min.de.num)
    all_pairs: optional restriction to a set of pairs -- a DataFrame with a `pair` column or an
        iterable of pair strings (R `all.pairs`)
    select_cl: optional restriction to clusters (R `select.cl`)

    Returns
    -------
    DataFrame [cl_up, cl_down_x, cl_down_y, up_num_x, down_num_x, up_num_y, down_num_y, P1, P2,
    pair, pair1, pair2], sorted by cl_up then by down_num_x + down_num_y ascending -- so for each
    candidate the most convincing triplet (fewest down genes) comes first. `find_doublets` stops at
    the first triplet that passes, so this order decides which one is reported.
    """
    s = summary
    asym = s[((s.up_num > min_up_num) & (s.down_num < max_down_num) & (s.up_num - s.down_num > min_up_num)) |
             ((s.down_num > min_up_num) & (s.up_num < max_down_num) & (s.down_num - s.up_num > min_up_num))].copy()

    if all_pairs is not None:
        keep = set(all_pairs['pair']) if isinstance(all_pairs, pd.DataFrame) else set(all_pairs)
        asym = asym[asym['pair'].isin(keep)]
    if select_cl is not None:
        sel = {str(c) for c in select_cl}
        asym = asym[asym.P1.astype(str).isin(sel) & asym.P2.astype(str).isin(sel)]

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
    # R renames up.num.new.* back to up.num.* here; do the same for all four (R's rename list has a
    # typo that leaves up.num.new.y untouched -- not reproduced).
    trip = trip.rename(columns={'up_num_o_x': 'up_num_x', 'down_num_o_x': 'down_num_x',
                                'up_num_o_y': 'up_num_y', 'down_num_o_y': 'down_num_y'})

    trip['P1'] = np.minimum(trip.cl_down_x, trip.cl_down_y)
    trip['P2'] = np.maximum(trip.cl_down_x, trip.cl_down_y)
    trip['pair'] = trip.P1.astype(str) + '_' + trip.P2.astype(str)
    if all_pairs is not None:
        keep = set(all_pairs['pair']) if isinstance(all_pairs, pd.DataFrame) else set(all_pairs)
        trip = trip[trip['pair'].isin(keep)]
    # the two parents must be genuinely different from each other (both directions have many DE genes)
    diff_pairs = set(s[(s.up_num > min_de_num) & (s.down_num > min_de_num)]['pair'])
    trip = trip[trip['pair'].isin(diff_pairs)]

    # the cl_up-vs-parent pair keys, canonical (lexicographic) like every other `pair` in the dataset
    def _canon(a, b):
        a = a.astype(str); b = b.astype(str)
        return np.where(a < b, a + '_' + b, b + '_' + a)
    trip['pair1'] = _canon(trip.cl_up, trip.cl_down_x)
    trip['pair2'] = _canon(trip.cl_up, trip.cl_down_y)

    # R: arrange(cl.up, down.num.x + down.num.y) -- fewest down genes first
    trip = trip.assign(_ord=trip.down_num_x + trip.down_num_y)
    trip = trip.sort_values(['cl_up', '_ord']).drop(columns='_ord').reset_index(drop=True)
    return trip


def _dir_genes(detail_by_pair: Dict[str, pd.DataFrame], c1: str, c2: str,
               top_n: Optional[int] = None) -> Dict[str, float]:
    """
    Genes higher in c1 than c2 (logPval), rank<=top_n.

    Direction lives in P1, as in scrattch.bigcat: rows with P1 == c1 are exactly the genes up in c1.
    `pair` stays canonical, so the lookup key is still the lexicographically ordered label pair.

    top_n=None applies no rank cap -- R's check_triplet_big deliberately takes ALL genes for the
    two overlap sets while capping the four scored sets at top.n.
    """
    a, b = (c1, c2) if c1 < c2 else (c2, c1)
    d = detail_by_pair.get(f"{a}_{b}")
    if d is None:
        return {}
    m = (d['P1'].astype(str) == str(c1))
    if top_n is not None:
        m &= (d['rank'] <= top_n)
    d = d[m]
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

    # Does the doublet (cl) inherit parent1's genes? -> overlap of (cl>c2) genes with (c1>c2) genes.
    # NOTE: R does NOT cap these two sets at top.n (its `tmp.genes` has no rank filter) even though
    # the four scored sets are capped. Capping them here would shrink the overlap and understate the
    # ratios, so top_n is deliberately omitted.
    cl_vs_c2 = set(_dir_genes(detail_by_pair, cl, c2).keys())
    ou1 = {g: v for g, v in up_genes.items() if g in cl_vs_c2};  ou1_s = trunc(ou1)
    cl_vs_c1 = set(_dir_genes(detail_by_pair, cl, c1).keys())
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


def find_doublets(detail: Optional[pd.DataFrame], triplets: pd.DataFrame, top_n: int = 50,
                  score_th: float = 0.8, olap_th: float = 1.6,
                  root: Optional[str] = None, cl_bin: Optional[Dict[str, int]] = None) -> pd.DataFrame:
    """
    Score every candidate's triplets and return the results. Mirrors R `find_doublets_all_big`.

    For each `cl_up`, triplets are tested in the order `find_triplets` produced (fewest down genes
    first) and testing STOPS at the first triplet with `score > score_th` and
    `olap_ratio_up_1 + olap_ratio_down_1 > olap_th`.

    **Every tested triplet is returned**, not just the winner -- R writes them all and the caller
    filters afterwards, often at a looser threshold than the early-stop one (the reference workflow
    stops at olap 1.6 but selects at 1.4). Returning only the best row would hide those.

    Parameters
    ----------
    detail: the per-gene detail as a DataFrame; pass None to stream from `root` instead
    triplets: output of find_triplets
    root: de_parquet directory -- when given, each triplet's rows are read from the partitions of
        its three clusters instead of holding the whole detail in memory (what R does, and what
        makes this usable on a large taxonomy)
    cl_bin: cluster -> bin map, needed with `root`

    Returns
    -------
    DataFrame with one row per tested triplet: [cl, cl1, cl2, up_num, down_num, score,
    olap_ratio_*, olap_num_*].
    """
    if detail is None and root is None:
        raise ValueError("pass either `detail` or `root`")
    from .de_all_pairs import read_de_pairs

    by_pair_all = None
    if detail is not None:
        by_pair_all = {p: d for p, d in detail.groupby('pair')}

    results = []
    for cl_up, tg in triplets.groupby('cl_up', sort=False):
        for _, row in tg.iterrows():
            c1, c2 = row['cl_down_x'], row['cl_down_y']
            if by_pair_all is not None:
                by_pair = by_pair_all
            else:
                trio = [str(cl_up), str(c1), str(c2)]
                d = read_de_pairs(root, pairs=[(a, b) for i, a in enumerate(trio)
                                               for b in trio[i + 1:]], cl_bin=cl_bin)
                by_pair = {p: g for p, g in d.groupby('pair')} if len(d) else {}
            r = check_triplet(by_pair, cl_up, c1, c2, top_n=top_n)
            results.append(r)
            if r['score'] > score_th and (r['olap_ratio_up_1'] + r['olap_ratio_down_1']) > olap_th:
                break
    return pd.DataFrame(results)

import numpy as np
import pandas as pd
import pytest

from transcriptomic_clustering.post_clustering_qc import _dir_genes


def test_dir_genes_follows_p1_not_lexicographic_order():
    """
    _dir_genes must key off P1 (the cluster the gene is up in), not off `sign`. Clusters '3' and '17'
    make the trap explicit: '17' < '3' as strings, so the canonical pair is "17_3" and the rows tagged
    sign='up' are the ones higher in 17, NOT in 3.
    """
    detail = pd.DataFrame({
        'pair': ['17_3'] * 4,
        'P1':   ['17', '17', '3', '3'],
        'P2':   ['3', '3', '17', '17'],
        'gene': ['up_in_17_a', 'up_in_17_b', 'up_in_3_a', 'up_in_3_b'],
        'logPval': [30.0, 25.0, 22.0, 21.0],
        'sign': ['up', 'up', 'down', 'down'],
        'rank': [1, 2, 1, 2],
    })
    by_pair = {'17_3': detail}

    assert _dir_genes(by_pair, '3', '17', top_n=10) == {'up_in_3_a': 22.0, 'up_in_3_b': 21.0}
    assert _dir_genes(by_pair, '17', '3', top_n=10) == {'up_in_17_a': 30.0, 'up_in_17_b': 25.0}
    # top_n cuts by rank within the direction
    assert _dir_genes(by_pair, '3', '17', top_n=1) == {'up_in_3_a': 22.0}
    # unknown pair -> empty
    assert _dir_genes(by_pair, '3', '99', top_n=10) == {}


def _summary(rows):
    """Minimal de_all_pairs summary: [pair, P1, P2, up_num, down_num]."""
    df = pd.DataFrame(rows, columns=['P1', 'P2', 'up_num', 'down_num'])
    df['pair'] = df.P1 + '_' + df.P2
    return df


def test_find_triplets_orders_by_fewest_down_genes():
    """
    R ends with arrange(cl.up, down.num.x + down.num.y), and find_doublets stops at the FIRST
    triplet that passes -- so this ordering decides which doublet gets reported, and must hold.
    """
    from transcriptomic_clustering.post_clustering_qc import find_triplets
    # cluster '9' looks like a doublet of 1, 2 and 3: asymmetric against each, with different
    # numbers of down genes so the triplet ordering is unambiguous.
    rows = [
        ('1', '9', 0, 100), ('2', '9', 0, 100), ('3', '9', 0, 100),   # 9 is the "up" side of each
        ('1', '2', 200, 200), ('1', '3', 200, 200), ('2', '3', 200, 200),  # parents well separated
    ]
    s = _summary(rows)
    # give each asymmetric pair a distinct down count so the sum orders the triplets
    s.loc[s.pair == '1_9', 'up_num'] = 3
    s.loc[s.pair == '2_9', 'up_num'] = 1
    s.loc[s.pair == '3_9', 'up_num'] = 2

    trip = find_triplets(s, min_up_num=30, max_down_num=10, min_de_num=50)
    assert set(trip.cl_up) == {'9'}
    sums = list(trip.down_num_x + trip.down_num_y)
    assert sums == sorted(sums), f"triplets must be ordered by fewest down genes first, got {sums}"
    # the most convincing triplet (parents 2 and 3, down counts 1+2) comes first
    assert {trip.iloc[0].cl_down_x, trip.iloc[0].cl_down_y} == {'2', '3'}


def test_find_triplets_emits_r_columns():
    """R's output carries pair1/pair2 and up.num.x/down.num.x style columns."""
    from transcriptomic_clustering.post_clustering_qc import find_triplets
    s = _summary([('1', '9', 2, 100), ('2', '9', 2, 100), ('1', '2', 200, 200)])
    trip = find_triplets(s, min_up_num=30, max_down_num=10, min_de_num=50)
    for col in ('cl_up', 'cl_down_x', 'cl_down_y', 'up_num_x', 'down_num_x',
                'up_num_y', 'down_num_y', 'P1', 'P2', 'pair', 'pair1', 'pair2'):
        assert col in trip.columns, f"missing {col}"
    r = trip.iloc[0]
    # pair1/pair2 are canonical (lexicographic), like every other `pair` in the dataset
    assert r.pair1 == '1_9' and r.pair2 == '2_9'
    assert r.pair == '1_2'
    # orientation: cl_up is the side with more up-genes
    assert r.cl_up == '9'


def test_find_triplets_respects_all_pairs_and_select_cl():
    """R's all.pairs / select.cl restrictions."""
    from transcriptomic_clustering.post_clustering_qc import find_triplets
    s = _summary([('1', '9', 2, 100), ('2', '9', 2, 100), ('3', '9', 2, 100),
                  ('1', '2', 200, 200), ('1', '3', 200, 200), ('2', '3', 200, 200)])
    assert len(find_triplets(s, min_de_num=50)) > 0

    # restricting the clusters drops 3 entirely -> only the 1/2 triplet survives
    trip = find_triplets(s, min_de_num=50, select_cl=['1', '2', '9'])
    assert set(trip.cl_down_x) | set(trip.cl_down_y) == {'1', '2'}

    # restricting all_pairs to exclude 3's asymmetric pair does the same
    trip = find_triplets(s, min_de_num=50, all_pairs=['1_9', '2_9', '1_2'])
    assert set(trip.cl_down_x) | set(trip.cl_down_y) == {'1', '2'}


def _detail(rows):
    """Detail frame in the de_parquet schema: P1 is the cluster the gene is up in."""
    df = pd.DataFrame(rows, columns=['P1', 'P2', 'gene', 'logPval', 'rank'])
    lo = np.minimum(df.P1, df.P2); hi = np.maximum(df.P1, df.P2)
    df['pair'] = lo + '_' + hi
    df['sign'] = np.where(df.P1 == lo, 'up', 'down')
    return df


def test_check_triplet_does_not_rank_cap_the_overlap_sets():
    """
    R caps the four SCORED gene sets at top.n but deliberately leaves the two OVERLAP sets uncapped
    (its `tmp.genes` has no rank filter). The scored sets are sums, so a cap is right; the overlap
    sets are membership tests -- "is parent1's top gene also up in the doublet?" -- and a gene can
    be rank 1 for the parents while ranking far lower for the doublet.

    Capping them halves the overlap ratios and suppresses real doublets, so this pins the split.
    """
    from transcriptomic_clustering.post_clustering_qc import check_triplet
    top_n = 50
    detail = _detail([
        # the two parents differ by gA (up in 1) and gB (up in 2), both rank 1
        ('1', '2', 'gA', 10.0, 1),
        ('2', '1', 'gB', 10.0, 1),
        # the doublet '9' carries both, but they rank BEYOND top_n against each parent
        ('9', '2', 'gA', 10.0, top_n + 50),
        ('9', '1', 'gB', 10.0, top_n + 50),
    ])
    by_pair = {p: d for p, d in detail.groupby('pair')}

    r = check_triplet(by_pair, '9', '1', '2', top_n=top_n)
    # gA is one of parent1's top genes AND is up in the doublet -> full overlap, despite rank 100
    assert r['olap_ratio_up_1'] == pytest.approx(1.0), (
        "the overlap set must not be rank-capped; capping it would score this real doublet as 0")
    assert r['olap_ratio_down_1'] == pytest.approx(1.0)
    assert r['olap_num_up_1'] == 1 and r['olap_num_down_1'] == 1


def test_find_doublets_returns_every_tested_triplet():
    """
    R writes every tested triplet and the caller filters afterwards -- often at a looser threshold
    than the early-stop one (the reference workflow stops at olap 1.6 but selects at 1.4). Returning
    only the best row per candidate would hide triplets in that band.
    """
    from transcriptomic_clustering.post_clustering_qc import find_doublets
    detail = _detail([
        ('1', '2', 'gA', 10.0, 1), ('2', '1', 'gB', 10.0, 1),
        ('1', '3', 'gC', 10.0, 1), ('3', '1', 'gD', 10.0, 1),
    ])
    triplets = pd.DataFrame({'cl_up': ['9', '9'], 'cl_down_x': ['1', '1'], 'cl_down_y': ['2', '3']})
    # nothing scores high enough to early-stop, so both triplets must come back
    out = find_doublets(detail, triplets, top_n=50, score_th=0.8, olap_th=1.6)
    assert len(out) == 2, f"expected every tested triplet, got {len(out)}"
    assert list(out['cl']) == ['9', '9']


def test_find_doublets_requires_detail_or_root():
    from transcriptomic_clustering.post_clustering_qc import find_doublets
    triplets = pd.DataFrame({'cl_up': ['9'], 'cl_down_x': ['1'], 'cl_down_y': ['2']})
    try:
        find_doublets(None, triplets)
    except ValueError as e:
        assert "detail" in str(e) and "root" in str(e)
    else:
        raise AssertionError("expected ValueError when neither detail nor root is given")

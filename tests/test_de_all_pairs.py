import os

import numpy as np
import pandas as pd
import pytest

from transcriptomic_clustering.de_all_pairs import (
    SCORE_CAP,
    _logpval,
    de_all_pairs,
    load_cl_bin,
    make_cl_bin,
    read_de_pairs,
)

THRESHOLDS = {
    'q1_thresh': 0.3, 'q2_thresh': None, 'cluster_size_thresh': 10,
    'qdiff_thresh': 0.3, 'padj_thresh': 0.05, 'lfc_thresh': 0.1,
}


@pytest.fixture
def cluster_stats():
    rng = np.random.default_rng(0)
    n_cl, n_genes = 6, 40
    labels = [str(i) for i in range(1, n_cl + 1)]
    genes = [f"g{i}" for i in range(n_genes)]
    means = pd.DataFrame(rng.random((n_cl, n_genes)) * 3, index=labels, columns=genes)
    variances = pd.DataFrame(rng.random((n_cl, n_genes)) * 0.1 + 0.05, index=labels, columns=genes)
    present = pd.DataFrame(np.clip(means.values / means.values.max(), 0, 1), index=labels, columns=genes)
    return means, variances, present, {c: 100 for c in labels}


def test_make_cl_bin_matches_r_binning():
    """cl.bin = ceiling(rank / bin_size), 1-based, numerically sorted."""
    assert make_cl_bin([str(i) for i in range(1, 7)], bin_size=2) == {
        '1': 1, '2': 1, '3': 2, '4': 2, '5': 3, '6': 3
    }


def test_de_all_pairs_binned_roundtrip(cluster_stats, tmp_path):
    """
    Writing to out_dir/summary_dir must produce scrattch.bigcat's bin.x=X/bin.y=Y layout, and
    reading a pair back must give exactly what the in-memory path returns for that pair.
    """
    pytest.importorskip("pyarrow")
    means, variances, present, cl_size = cluster_stats
    out_dir, summary_dir = str(tmp_path / "de_parquet"), str(tmp_path / "de_summary")

    summary_mem, detail_mem = de_all_pairs(means, variances, present, cl_size, THRESHOLDS, top_n=5)

    summary_disk, detail_disk = de_all_pairs(
        means, variances, present, cl_size, THRESHOLDS, top_n=5,
        out_dir=out_dir, summary_dir=summary_dir, bin_size=2,
    )
    # nothing is held in memory once it is on disk
    assert summary_disk is None and detail_disk is None

    # The detail encodes direction in P1/P2 (P1 = the cluster the gene is up in), as scrattch.bigcat
    # does, so each pair contributes to both orientations of its bin-pair: 3 bins -> the full 3x3 grid.
    partitions = sorted(
        os.path.relpath(root, out_dir)
        for root, _, files in os.walk(out_dir)
        if any(f.endswith(".parquet") for f in files)
    )
    assert partitions == [
        "bin.x=1/bin.y=1", "bin.x=1/bin.y=2", "bin.x=1/bin.y=3",
        "bin.x=2/bin.y=1", "bin.x=2/bin.y=2", "bin.x=2/bin.y=3",
        "bin.x=3/bin.y=1", "bin.x=3/bin.y=2", "bin.x=3/bin.y=3",
    ]
    # the summary keeps one row per unordered pair, so it stays on the upper triangle
    summary_partitions = sorted(
        os.path.relpath(root, summary_dir)
        for root, _, files in os.walk(summary_dir)
        if any(f.endswith(".parquet") for f in files)
    )
    assert summary_partitions == [
        "bin.x=1/bin.y=1", "bin.x=1/bin.y=2", "bin.x=1/bin.y=3",
        "bin.x=2/bin.y=2", "bin.x=2/bin.y=3", "bin.x=3/bin.y=3",
    ]
    assert load_cl_bin(out_dir) == make_cl_bin(means.index, bin_size=2)

    # a pair read back from its partition matches the in-memory result exactly
    for pair in [("1", "2"), ("2", "3"), ("4", "6")]:
        key = f"{pair[0]}_{pair[1]}"
        got = read_de_pairs(out_dir, pairs=[pair]).sort_values(["gene", "sign"]).reset_index(drop=True)
        expected = detail_mem[detail_mem.pair == key].sort_values(["gene", "sign"]).reset_index(drop=True)
        pd.testing.assert_frame_equal(got[expected.columns], expected)

    # summary round-trips unchanged
    back = read_de_pairs(summary_dir).sort_values("pair").reset_index(drop=True)
    expected = summary_mem.sort_values("pair").reset_index(drop=True)
    pd.testing.assert_frame_equal(back[expected.columns], expected)


def test_logpval_is_not_capped():
    """
    R stores the raw -log10(padj) per gene and truncates only the copy it sums into the de score
    (de.genes.R: `up.genes = -log10(padj)`; `tmp = up.genes; tmp[tmp > 20] = 20`).
    """
    lp = _logpval(np.array([1e-25, 1e-30, 1e-3]))
    assert lp[0] == pytest.approx(25.0)
    assert lp[1] == pytest.approx(30.0)
    assert lp[2] == pytest.approx(3.0)
    assert (lp[:2] > SCORE_CAP).all(), "detail values must not be truncated at 20"


def test_score_is_sum_of_capped_detail(cluster_stats):
    """The cap lives in the score, not the stored value: score == sum(min(logPval, 20))."""
    means, variances, present, cl_size = cluster_stats
    summary, detail = de_all_pairs(
        means, variances, present, cl_size, THRESHOLDS, top_n=10 ** 6,
    )
    assert not detail.empty

    capped = detail.assign(c=np.minimum(detail['logPval'], SCORE_CAP))
    recomputed = capped.groupby('pair')['c'].sum()
    expected = summary.set_index('pair')['score'].reindex(recomputed.index)
    np.testing.assert_allclose(expected.values, recomputed.values, rtol=1e-9)

    # up/down split too
    for sign, col in (('up', 'up_score'), ('down', 'down_score')):
        sub = capped[capped['sign'] == sign].groupby('pair')['c'].sum()
        exp = summary.set_index('pair')[col].reindex(sub.index)
        np.testing.assert_allclose(exp.values, sub.values, rtol=1e-9)


def test_detail_encodes_direction_in_p1_p2(cluster_stats):
    """
    scrattch.bigcat stores direction in P1/P2: P1 is always the cluster the gene is UP in, `pair`
    keeps the canonical lexicographic label pair, and `sign` is retained but redundant --
    sign == 'up' exactly when P1 < P2 lexicographically.
    """
    means, variances, present, cl_size = cluster_stats
    _, detail = de_all_pairs(means, variances, present, cl_size, THRESHOLDS, top_n=5)
    assert not detail.empty

    p1 = detail['P1'].astype(str)
    p2 = detail['P2'].astype(str)
    expected_sign = np.where(p1.values < p2.values, 'up', 'down')
    assert (expected_sign == detail['sign'].values).all()

    expected_pair = [f"{a}_{b}" if a < b else f"{b}_{a}" for a, b in zip(p1, p2)]
    assert expected_pair == list(detail['pair'].astype(str))


def _write_detail(tmp_path, rows, cl_bin):
    """Write a hand-made detail dataset in the de_parquet layout (direction encoded in P1/P2)."""
    pytest.importorskip("pyarrow")
    from transcriptomic_clustering.de_all_pairs import _save_cl_bin, _write_partition
    root = str(tmp_path / "de_parquet")
    df = pd.DataFrame(rows)
    for (bx, by), grp in df.groupby([df.P1.map(cl_bin), df.P2.map(cl_bin)]):
        _write_partition(grp.reset_index(drop=True), root, int(bx), int(by))
    _save_cl_bin(root, cl_bin)
    return root


def test_get_gene_score_ds_is_rank_sum():
    """score(gene) = sum over the selected pairs of (max_num - rank), R get_gene_score_ds."""
    from transcriptomic_clustering.de_all_pairs import MAX_NUM
    assert MAX_NUM == 1000


def test_select_top_pos_markers_ds_ranks_by_breadth_then_rank(tmp_path):
    """
    A gene that is a top marker against MANY clusters must beat one that is rank-1 against a single
    cluster: each win contributes ~max_num, so breadth dominates and rank only breaks ties.
    """
    cl_bin = {'1': 1, '2': 1, '3': 1, '4': 1}
    rows = []
    # 'broad' is rank 5 against all three partners; 'narrow' is rank 1 but only against cluster 2
    for p2 in ('2', '3', '4'):
        rows.append(dict(pair=f"1_{p2}", P1='1', P2=p2, gene='broad',
                         logPval=9.0, sign='up', rank=5, lfc=1.0))
    rows.append(dict(pair="1_2", P1='1', P2='2', gene='narrow',
                     logPval=9.0, sign='up', rank=1, lfc=1.0))
    # a gene belonging to another cluster's direction, which must not leak into cluster 1's markers
    rows.append(dict(pair="1_2", P1='2', P2='1', gene='other_dir',
                     logPval=9.0, sign='down', rank=1, lfc=1.0))
    root = _write_detail(tmp_path, rows, cl_bin)

    from transcriptomic_clustering.de_all_pairs import select_top_pos_markers_ds
    got = select_top_pos_markers_ds(root, ['1', '2', '3', '4'], select_cl=['1'],
                                    cl_bin=cl_bin, n_markers=2)
    assert got['1'] == ['broad', 'narrow']

    # the gene universe is honoured
    got = select_top_pos_markers_ds(root, ['1', '2', '3', '4'], select_cl=['1'],
                                    genes={'narrow'}, cl_bin=cl_bin, n_markers=2)
    assert got['1'] == ['narrow']


def test_select_markers_pair_group_top_ds_directions(tmp_path):
    """'up' scores group1-over-group2; 'down' is the reverse direction, read off P1."""
    cl_bin = {'1': 1, '2': 1}
    rows = [
        dict(pair="1_2", P1='1', P2='2', gene='up_in_1', logPval=9.0, sign='up', rank=1, lfc=1.0),
        dict(pair="1_2", P1='2', P2='1', gene='up_in_2', logPval=9.0, sign='down', rank=1, lfc=1.0),
    ]
    root = _write_detail(tmp_path, rows, cl_bin)

    from transcriptomic_clustering.de_all_pairs import select_markers_pair_group_top_ds
    r = select_markers_pair_group_top_ds(root, ['1'], ['2'], cl_bin=cl_bin,
                                         select_sign=('up', 'down'), n_markers=5)
    assert r['up_genes'] == ['up_in_1']
    assert r['down_genes'] == ['up_in_2']
    # requesting only 'up' leaves the other side empty
    r = select_markers_pair_group_top_ds(root, ['1'], ['2'], cl_bin=cl_bin, n_markers=5)
    assert r['up_genes'] == ['up_in_1'] and r['down_genes'] == []


def test_check_pairs_lfc_counts_separating_genes():
    """Counts genes whose cluster-mean difference exceeds lfc_th, in the P1 - P2 direction."""
    from transcriptomic_clustering.de_all_pairs import check_pairs_lfc
    cl_means = pd.DataFrame({'a': [10.0, 0.0], 'b': [0.0, 0.0], 'c': [10.0, 10.0]},
                            index=['1', '2'])
    to_add = pd.DataFrame({'P1': ['1', '2'], 'P2': ['2', '1']})
    out = check_pairs_lfc(to_add, ['a', 'b', 'c'], cl_means, lfc_th=2.0)
    # 1 vs 2: only 'a' separates (10-0). 2 vs 1: nothing (a is -10, b and c are 0).
    assert list(out['checked']) == [1, 0]


def test_check_pairs_ds_counts_stored_markers(tmp_path):
    from transcriptomic_clustering.de_all_pairs import check_pairs_ds
    cl_bin = {'1': 1, '2': 1, '3': 1}
    rows = [
        dict(pair="1_2", P1='1', P2='2', gene='g1', logPval=9.0, sign='up', rank=1, lfc=1.0),
        dict(pair="1_2", P1='1', P2='2', gene='g2', logPval=8.0, sign='up', rank=2, lfc=1.0),
        dict(pair="1_3", P1='1', P2='3', gene='g1', logPval=7.0, sign='up', rank=1, lfc=1.0),
    ]
    root = _write_detail(tmp_path, rows, cl_bin)
    to_add = pd.DataFrame({'P1': ['1', '1'], 'P2': ['2', '3']})
    out = check_pairs_ds(root, to_add, ['g1', 'g2'], cl_bin=cl_bin)
    assert list(out['checked']) == [2, 1]
    # a gene that is not a marker for either pair contributes nothing
    out = check_pairs_ds(root, to_add, ['absent'], cl_bin=cl_bin)
    assert list(out['checked']) == [0, 0]


def test_select_pos_markers_ds_covers_pairs_the_top_marker_misses(tmp_path):
    """
    The greedy variant must keep adding genes until every comparison is covered -- that is what
    distinguishes it from select_top_pos_markers_ds, which just takes the best-scoring gene.
    """
    from transcriptomic_clustering.de_all_pairs import (
        check_pairs_ds, select_pos_markers_ds, select_top_pos_markers_ds)
    cl_bin = {'1': 1, '2': 1, '3': 1, '4': 1}
    rows = []
    # 'wide' marks cluster 1 against 2 and 3 but NOT 4; 'only4' marks it against 4 alone
    for p2, rk in (('2', 1), ('3', 1)):
        rows.append(dict(pair=f"1_{p2}", P1='1', P2=p2, gene='wide',
                         logPval=9.0, sign='up', rank=rk, lfc=1.0))
    rows.append(dict(pair="1_4", P1='1', P2='4', gene='only4',
                     logPval=9.0, sign='up', rank=1, lfc=1.0))
    root = _write_detail(tmp_path, rows, cl_bin)

    clusters = ['1', '2', '3', '4']
    top = select_top_pos_markers_ds(root, clusters, select_cl=['1'], cl_bin=cl_bin, n_markers=1)
    assert top['1'] == ['wide']                      # best single marker, but misses cluster 4

    combo = select_pos_markers_ds(root, clusters, select_cl=['1'], cl_bin=cl_bin,
                                  n_markers=1, max_genes=10)
    assert combo['1'][0] == 'wide'                   # starts from the top marker
    assert set(combo['1']) == {'wide', 'only4'}      # then covers the pair it missed

    to_add = pd.DataFrame({'P1': ['1', '1', '1'], 'P2': ['2', '3', '4']})
    assert (check_pairs_ds(root, to_add, combo['1'], cl_bin=cl_bin)['checked'] >= 1).all()


def test_select_pos_markers_ds_cache_roundtrip(tmp_path):
    """out_dir/overwrite mirror R's per-cluster caching so a long run can resume."""
    from transcriptomic_clustering.de_all_pairs import select_pos_markers_ds
    cl_bin = {'1': 1, '2': 1}
    rows = [dict(pair="1_2", P1='1', P2='2', gene='g1', logPval=9.0, sign='up', rank=1, lfc=1.0)]
    root = _write_detail(tmp_path, rows, cl_bin)
    cache = str(tmp_path / "cl_markers")

    first = select_pos_markers_ds(root, ['1', '2'], select_cl=['1'], cl_bin=cl_bin,
                                  n_markers=1, out_dir=cache)
    assert first['1'] == ['g1']
    assert os.path.exists(os.path.join(cache, "1.markers.json"))
    # with overwrite=False the cached value is returned even if the dataset is gone
    again = select_pos_markers_ds("/nonexistent", ['1', '2'], select_cl=['1'], cl_bin=cl_bin,
                                  n_markers=1, out_dir=cache, overwrite=False)
    assert again['1'] == ['g1']

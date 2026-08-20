import os

import numpy as np
import pandas as pd
import pytest

from transcriptomic_clustering.post_clustering_qc import (
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

    # 3 bins -> upper triangle -> 6 partitions per dataset
    partitions = sorted(
        os.path.relpath(root, out_dir)
        for root, _, files in os.walk(out_dir)
        if any(f.endswith(".parquet") for f in files)
    )
    assert partitions == [
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

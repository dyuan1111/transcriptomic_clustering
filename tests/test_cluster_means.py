import os
import pytest

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from transcriptomic_clustering import cluster_means as cm
from scipy.sparse import csr_matrix


@pytest.fixture
def adata():

    X = np.array([
        [0, 1, 3],
        [0, 6, 2],
        [1, 0, 1],
        [4, 4, 7],
        [4, 0, 0],
        [6, 0, 6],
        [8, 6, 3],
        [2, 1, 6],
        [0, 0, 7],
        [2, 1, 0],
    ])
    n_obs = X.shape[0]
    n_var = X.shape[1]
    cell_names = [f"cell_{i}" for i in range(n_obs)]
    obs = pd.DataFrame(index=cell_names)

    var_names = [f"var_{i}" for i in range(n_var)]
    var = pd.DataFrame(index=var_names)
    adata = ad.AnnData(csr_matrix(X), obs=obs, var=var)

    return adata


@pytest.fixture
def clusters(adata):
    cluster_assignments = {
        '11': [0, 3, 5, 9],
        '2': [1, 2, 6],
        '32': [4, 7],
        '4': [8]
    }

    cluster_by_obs = np.array(['11', '2', '2', '11', '32', '11', '2', '32', '4', '11'])

    cluster_means = pd.DataFrame(
        np.array([[3., 1.5, 4.],
                  [3., 4., 2.],
                  [3., 0.5, 3.],
                  [0., 0., 7.]]),
        index = ['11', '2', '32', '4'],
        columns=adata.var_names
    )

    present_cluster_means = pd.DataFrame(
        np.array([[0.5, 0.25, 0.75],
                  [(1/3), (2/3), (1/3)],
                  [0.5, 0., .5],
                  [0., 0., 1.]]),
        index = ['11', '2', '32', '4'],
        columns=adata.var_names
    )

    cluster_variances = pd.DataFrame(
        np.array([[6+2/3, 3, 10],
                  [19, 12.0, 1.0],
                  [2.0, 0.5, 18.0],
                  [0.0, 0.0, 0.0]]),
        index=['11', '2', '32', '4'],
        columns=adata.var_names
    )

    return cluster_means, present_cluster_means, cluster_variances, cluster_assignments, cluster_by_obs


def test_get_cluster_means_inmemory(adata, clusters):

    (expected_cluster_means,
     expected_present_cluster_means,
     expected_cluster_variances,
     cluster_assignments,
     cluster_by_obs) = clusters

    obtained_cluster_means, obtained_present_cluster_means, obtained_cluster_variances = \
        cm.get_cluster_means(adata, cluster_assignments, cluster_by_obs, low_th=2)

    assert obtained_cluster_means.index.equals(expected_cluster_means.index)
    assert obtained_cluster_means.columns.equals(expected_cluster_means.columns)
    assert np.allclose(obtained_cluster_means.to_numpy(), expected_cluster_means.to_numpy())
    assert obtained_present_cluster_means.index.equals(expected_present_cluster_means.index)
    assert obtained_present_cluster_means.columns.equals(expected_present_cluster_means.columns)
    assert np.allclose(obtained_present_cluster_means.to_numpy(), expected_present_cluster_means.to_numpy())
    assert obtained_cluster_variances.index.equals(expected_cluster_variances.index)
    assert obtained_cluster_variances.columns.equals(expected_cluster_variances.columns)
    assert np.allclose(obtained_cluster_variances.to_numpy(), expected_cluster_variances.to_numpy())


def test_get_cluster_means_backed(adata, clusters, tmpdir_factory):

    (expected_cluster_means,
     expected_present_cluster_means,
     expected_cluster_variances,
     cluster_assignments,
     cluster_by_obs) = clusters

    tmpdir = str(tmpdir_factory.mktemp("test_cluster_means"))
    input_file_name = os.path.join(tmpdir, "input.h5ad")

    ad.AnnData(csr_matrix(adata.X), obs=adata.obs, var=adata.var).write(input_file_name) # make tmp input file

    adata = sc.read_h5ad(input_file_name, backed='r')
    obtained_cluster_means, obtained_present_cluster_means, obtained_cluster_variances = \
        cm.get_cluster_means(adata, cluster_assignments, cluster_by_obs, low_th=2)

    assert obtained_cluster_means.index.equals(expected_cluster_means.index)
    assert obtained_cluster_means.columns.equals(expected_cluster_means.columns)
    assert np.allclose(obtained_cluster_means.to_numpy(), expected_cluster_means.to_numpy())
    assert obtained_present_cluster_means.index.equals(expected_present_cluster_means.index)
    assert obtained_present_cluster_means.columns.equals(expected_present_cluster_means.columns)
    assert np.allclose(obtained_present_cluster_means.to_numpy(), expected_present_cluster_means.to_numpy())
    assert obtained_cluster_variances.index.equals(expected_cluster_variances.index)
    assert obtained_cluster_variances.columns.equals(expected_cluster_variances.columns)
    assert np.allclose(obtained_cluster_variances.to_numpy(), expected_cluster_variances.to_numpy())

def test_present_is_exact_at_half_detection():
    """
    `present` is a detection rate k/n, so it lands on 0.5 exactly, and the q1 filter is a STRICT
    `q1 > q1_thresh` with q1_thresh routinely 0.5. It must therefore be computed as an exact rational.

    np.mean() on a sparse boolean mask is a matvec against a vector of 1/n -- it adds 1/n to itself k
    times -- so an exact half comes out one ulp ABOVE 0.5 and wrongly passes the filter. Counting and
    dividing once (as scrattch.bigcat's C++ does) keeps it exact.
    """
    n_cells, n_detect = 454, 227          # 227/454 is exactly 0.5
    col = np.zeros((n_cells, 1))
    col[:n_detect, 0] = 2.0               # above low_th=1, the rest are 0
    adata_half = ad.AnnData(csr_matrix(col), dtype=np.float64)

    assignments = {'1': np.arange(n_cells)}
    cluster_by_obs = np.array(['1'] * n_cells, dtype=object)
    _, present, _ = cm.get_cluster_means_inmemory(adata_half, assignments, low_th=1)

    value = present.to_numpy().ravel()[0]
    assert value == 0.5, f"expected exactly 0.5, got {value!r}"
    assert not (value > 0.5), "a half-detected gene must NOT pass a strict `> 0.5` filter"

    # the buggy formulation, kept here to document what is being guarded against
    naive = float(np.asarray(np.mean((adata_half.X > 1), axis=0)).ravel()[0])
    assert naive > 0.5 and naive != 0.5, "sparse np.mean should exhibit the ulp error being guarded"

from annoy import AnnoyIndex
import scanpy as sc
import numpy as np
import pandas as pd
import pytest
import csv
import os

from mock import patch
from unittest.mock import MagicMock

import transcriptomic_clustering as tc
from transcriptomic_clustering.clustering import (
    cluster_louvain_phenograph,
    _annoy_build_csr_nn_graph,
    _uniform_csr_from_nn_dict,
)

DIR_NAME = os.path.dirname(__file__)
TEST_RANDOM_SEED = 5
TEST_K = 15

@pytest.fixture
def pca_data():
    return sc.read_csv(os.path.join(DIR_NAME, "data/output_dimreduction_pca.csv"), first_column_names=True)

@pytest.fixture
def sample_graph():
    return sc.read_h5ad(os.path.join(DIR_NAME, "data/sample_graph.h5ad"))

@pytest.fixture
def sample_index_file():
    return os.path.join(DIR_NAME, "data/sample_annoy_index.ann")

@pytest.fixture
def sample_nn_dict():
    return {
        0: [2, 4],
        1: [3, 5],
        2: [0, 4],
        3: [1, 5],
        4: [0, 2],
        5: [1, 3]
    }

@pytest.fixture
def annoy_taynaud_sample_partition():
    sample_output_adata = sc.read_h5ad(os.path.join(DIR_NAME, "data/annoy_taynaud_15_nn_rs_5.h5ad"))
    sample_cluster_by_obs = sample_output_adata.obs['pheno_louvain'].values.tolist()
    return sample_cluster_by_obs

@pytest.fixture
def annoy_vtraag_sample_partition():
    sample_output_adata = sc.read_h5ad(os.path.join(DIR_NAME, "data/annoy_vtraag_15_nn_rs_5.h5ad"))
    sample_cluster_by_obs = sample_output_adata.obs['pheno_louvain'].values.tolist()
    return sample_cluster_by_obs

@pytest.fixture
def r_data():
    with open(os.path.join(DIR_NAME, "data/R_tasic_results.csv")) as csvf:
        myreader = csv.reader(csvf)
        r_list = []
        for row in myreader:
            r_list.append(row[1])
    return r_list[1:]

def test_phenograph_adata_not_inplace(pca_data, r_data):
    _, _, _, q = cluster_louvain_phenograph(pca_data, 10)
    assert q > .8

def test_phenograph_adata_inplace(pca_data, r_data):
    cluster_by_obs, _, _, q = cluster_louvain_phenograph(pca_data, 10, annotate=True)

    print(f"Cluster by obs has type {type(cluster_by_obs)}")
    print(f"pca obs has type {type(pca_data.obs['pheno_louvain'])}")
    assert cluster_by_obs.tolist() == pca_data.obs['pheno_louvain'].to_list()
    assert q > .8

def test_phenograph_outlier_label_handling():
    test_adata = MagicMock()
    test_adata.X = None
    with patch.object(tc.clustering, 'phenograph', return_value=([-1, 1, -1, 2, -1, 3], None, None)):
        cluster_by_obs, obs_by_cluster, _, _ = cluster_louvain_phenograph(test_adata, 4, False)
    assert cluster_by_obs == [4, 1, 5, 2, 6, 3]
    assert obs_by_cluster == {1: [1], 4: [0], 5: [2], 2: [3], 6: [4], 3: [5]}

def test_uniform_csr_from_nn_dict(sample_nn_dict):
    expected_matrix = np.array([
        [0., 0., 1., 0., 1., 0.],
        [0., 0., 0., 1., 0., 1.],
        [1., 0., 0., 0., 1., 0.],
        [0., 1., 0., 0., 0., 1.],
        [1., 0., 1., 0., 0., 0.],
        [0., 1., 0., 1., 0., 0.]
    ])
    result_matrix = tc.clustering._uniform_csr_from_nn_dict(sample_nn_dict).todense()
    assert np.array_equal(expected_matrix, result_matrix)

def test_jaccard_csr_from_nn_dict(sample_nn_dict):
    expected_matrix = np.array([
        [0., 0., 1/3, 0., 1/3, 0.],
        [0., 0., 0., 1/3, 0., 1/3],
        [1/3, 0., 0., 0., 1/3, 0.],
        [0., 1/3, 0., 0., 0., 1/3],
        [1/3, 0., 1/3, 0., 0., 0.],
        [0., 1/3, 0., 1/3, 0., 0.]
    ])
    result_matrix = tc.clustering._jaccard_csr_from_nn_dict(sample_nn_dict).todense()
    assert np.array_equal(expected_matrix, result_matrix)

def test_annoy_build_csr_nn_graph_reproducibility(pca_data, sample_graph, sample_index_file):
    graph_csr = tc.clustering._annoy_build_csr_nn_graph(pca_data, sample_index_file, TEST_K)
    assert np.array_equal(sample_graph.X.toarray(), graph_csr.toarray().astype('float32'))

def test_annoy_build_csr_nn_graph_reproducibility_uniform(pca_data, sample_graph, sample_index_file):
    uniform_sample_graph_arr = (sample_graph.X.toarray() > 0).astype('float32')
    graph_csr = tc.clustering._annoy_build_csr_nn_graph(pca_data, sample_index_file, TEST_K, weighting_method='uniform')
    assert np.array_equal(uniform_sample_graph_arr, graph_csr.toarray())

def test_annoy_build_csr_nn_graph_reproducibility_multithread(pca_data, sample_graph, sample_index_file):
    graph_csr = tc.clustering._annoy_build_csr_nn_graph(pca_data, sample_index_file, TEST_K, n_jobs=8)
    assert np.array_equal(sample_graph.X.toarray(), graph_csr.toarray().astype('float32'))

def test_annoy_build_csr_nn_graph_bad_weighting(pca_data, sample_index_file):
    with pytest.raises(ValueError) as err:
        graph_csr = tc.clustering._annoy_build_csr_nn_graph(pca_data, sample_index_file, TEST_K, weighting_method='fake method')
    assert str(err.value) == 'fake method is not a valid weighting option! Must use jaccard or uniform'

def test_get_annoy_knn_new_index(pca_data, sample_graph, sample_index_file):
    graph_adata = tc.clustering.get_annoy_knn(pca_data, TEST_K, random_seed=TEST_RANDOM_SEED)
    assert np.array_equal(sample_graph.X.toarray(), graph_adata.X.toarray())

def test_get_annoy_knn_reproducibility(pca_data, sample_graph, sample_index_file):
    graph_adata = tc.clustering.get_annoy_knn(pca_data, TEST_K, random_seed=TEST_RANDOM_SEED, annoy_index_filename=sample_index_file)
    assert np.array_equal(sample_graph.X.toarray(), graph_adata.X.toarray())

def test_get_annoy_knn_reproducibility_multi_thread(pca_data, sample_graph, sample_index_file):
    graph_adata = tc.clustering.get_annoy_knn(pca_data, TEST_K, n_jobs=8, random_seed=TEST_RANDOM_SEED, annoy_index_filename=sample_index_file)
    assert np.array_equal(sample_graph.X.toarray(), graph_adata.X.toarray())

def test_get_taynaud_louvain_behaviour(sample_graph):
    sample_partition = {i: i % 4 for i in range(284)}
    sample_modularity = 0.6
    expected_cluster_by_obs = [i % 4 for i in range(284)]

    with patch.object(tc.clustering.community_louvain, 'best_partition', return_value=sample_partition):
        with patch.object(tc.clustering.community_louvain, 'modularity', return_value=sample_modularity):
            cluster_by_obs, q = tc.clustering.get_taynaud_louvain(sample_graph)

    assert cluster_by_obs == expected_cluster_by_obs
    assert q == sample_modularity

def test_get_taynaud_louvain_reproducibility(sample_graph, annoy_taynaud_sample_partition):
    cluster_by_obs, _ = tc.clustering.get_taynaud_louvain(sample_graph, random_seed=TEST_RANDOM_SEED)
    assert cluster_by_obs == annoy_taynaud_sample_partition

def test_get_vtraag_leiden_reproducibility(sample_graph, annoy_vtraag_sample_partition):
    cluster_by_obs, _ = tc.clustering.get_vtraag_leiden(sample_graph, random_seed=TEST_RANDOM_SEED)
    assert cluster_by_obs == annoy_vtraag_sample_partition

def test_annoy_parallel_search_matches_reference(tmp_path):
    """
    Regression test: each pool worker loads the annoy index once (via the initializer) and
    queries a block of observations. The resulting graph must match a direct, single-index
    reference -- i.e. every observation is queried exactly once and reassembled at its own index.
    """
    rng = np.random.default_rng(0)
    data = rng.random((200, 8)).astype(np.float32)

    index_file = str(tmp_path / "index.ann")
    ai = AnnoyIndex(data.shape[1], "euclidean")
    ai.on_disk_build(index_file)
    ai.set_seed(TEST_RANDOM_SEED)
    for i, row in enumerate(data):
        ai.add_item(i, row)
    ai.build(10)

    reference = {i: ai.get_nns_by_item(i, TEST_K) for i in range(data.shape[0])}
    expected = _uniform_csr_from_nn_dict(dict(reference))

    for n_jobs in (1, 4):
        actual = _annoy_build_csr_nn_graph(
            data, index_file, k=TEST_K, n_jobs=n_jobs, weighting_method="uniform"
        )
        assert actual.shape == expected.shape
        assert (actual != expected).nnz == 0, f"graph differs for n_jobs={n_jobs}"


def test_pynndescent_handles_subsets_smaller_than_k():
    """
    The recursive pipeline drills down to groups of a few dozen cells, where k can exceed the number
    of points. annoy returns min(k, n) neighbours; pynndescent instead pads unfilled slots with -1,
    which reaches csr_matrix as a negative index and aborts graph construction. Both backends must
    survive these sizes.
    """
    pytest.importorskip("pynndescent")
    import anndata as ad
    from transcriptomic_clustering.clustering import cluster_louvain
    rng = np.random.default_rng(0)
    for n in (8, 15, 29):
        X = rng.normal(0, 1, (n, 10)).astype(np.float32)
        adata = ad.AnnData(X, obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
        for method in ('annoy', 'pynndescent'):
            _, obs_by_cluster, _, _ = cluster_louvain(
                adata, k=15, knn_method=method, weighting_method='jaccard_snn',
                louvain_method='vtraag', resolution=1.0, n_jobs=1, knn_seed=1)
            assert sum(len(v) for v in obs_by_cluster.values()) == n, (
                f"{method} lost cells at n={n}")


def test_pynndescent_clamps_n_jobs_to_numba_thread_cap():
    """
    pynndescent calls numba.set_num_threads(n_jobs), which raises above NUMBA_NUM_THREADS. The annoy
    backend uses a multiprocessing Pool with no such ceiling, so the same n_jobs that works there
    must not abort here.
    """
    pytest.importorskip("pynndescent")
    import anndata as ad
    from transcriptomic_clustering.clustering import cluster_louvain
    rng = np.random.default_rng(0)
    X = np.vstack([rng.normal(m, 0.4, (60, 8)) for m in (0, 5)]).astype(np.float32)
    adata = ad.AnnData(X, obs=pd.DataFrame(index=[f"c{i}" for i in range(120)]))
    _, obs_by_cluster, _, _ = cluster_louvain(
        adata, k=15, knn_method='pynndescent', weighting_method='jaccard_snn',
        louvain_method='vtraag', resolution=1.0, n_jobs=512, knn_seed=1)
    assert sum(len(v) for v in obs_by_cluster.values()) == 120


def test_default_knn_method_is_pynndescent():
    """
    The default backend is pynndescent: it agrees with scrattch.bigcat slightly better than annoy,
    is exactly reproducible at a fixed seed, and is more stable across seeds. Pinned here because
    changing it silently would change every downstream clustering result.
    """
    import inspect
    from transcriptomic_clustering.clustering import cluster_louvain
    assert inspect.signature(cluster_louvain).parameters['knn_method'].default == 'pynndescent'


def test_annoy_seed_is_a_working_deprecated_alias_for_knn_seed():
    """
    `annoy_seed` predates the pynndescent backend and has always seeded whichever backend is in use,
    so the accurate name is `knn_seed`. The old name must keep working -- existing pipeline configs
    pass it -- but should warn, and must produce the identical graph.
    """
    pytest.importorskip("pynndescent")
    import warnings as _w
    import anndata as ad
    from sklearn.metrics import adjusted_rand_score
    from transcriptomic_clustering.clustering import cluster_louvain
    rng = np.random.default_rng(0)
    X = np.vstack([rng.normal(m, 0.4, (120, 10)) for m in (0, 5)]).astype(np.float32)
    adata = ad.AnnData(X, obs=pd.DataFrame(index=[f"c{i}" for i in range(240)]))

    def run(**kw):
        with _w.catch_warnings(record=True) as caught:
            _w.simplefilter("always")
            cbo, _, _, _ = cluster_louvain(adata, k=15, weighting_method='jaccard_snn',
                                           louvain_method='vtraag', resolution=1.0, n_jobs=1, **kw)
            return cbo, [x for x in caught if issubclass(x.category, DeprecationWarning)]

    new, warn_new = run(knn_seed=1)
    old, warn_old = run(annoy_seed=1)
    default, warn_default = run()

    assert adjusted_rand_score(old, new) == 1.0, "the alias must give the identical graph"
    assert adjusted_rand_score(default, new) == 1.0, "the default must still be seed 1"
    assert len(warn_old) == 1, "passing annoy_seed should warn"
    assert not warn_new and not warn_default, "knn_seed and the default must not warn"

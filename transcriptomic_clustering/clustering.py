import os
import functools
import warnings
from datetime import datetime
import annoy
import community as community_louvain
import networkx as nx
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import leidenalg
from annoy import AnnoyIndex
from scanpy import AnnData
from scipy.sparse import csr_matrix, issparse
from math import log
import tempfile
from typing import List
import igraph as ig
from scanpy.external.tl import phenograph


def cluster_louvain_phenograph(
    adata: AnnData,
    k: int=15,
    annotate: bool = False,
    **kwargs):
    """
    Uses the louvain phenograph algorithm to calculate clusters in the given adata.
    Out of concern for memory the annotate behavior and behavior on AnnData of the
    scanpy pheongraph function is overridden.

    Parameters
    -----------
    adata: an AnnData of data to cluster
    k: number of nearest neighbors to use during graph construction
    annotate: if True, community labels are added to adata as observations. Otherwise
             all scanpy pheongraph outputs are returned.

    Returns
    -----------
    cluster_by_obs: an array of community labels for each cell in adata, in order
    obs_by_cluster: a map of community labels to lists of cell indices in that community in order
    graph: the calculated adjacency graph on the adata
    q: the maximum modularity of the final clustering
    """
    if 'copy' in kwargs:
        del kwargs['copy']

    cluster_by_obs, graph, q = phenograph(
        adata=adata.X,
        k=k,
        **kwargs
    )

    obs_by_cluster = _cluster_obs_list_to_dict(cluster_by_obs)

    # phenograph returns all unclustered cells with a -1 label
    # this code splits that null cluster
    if -1 in obs_by_cluster:
        max_cluster = max(list(obs_by_cluster.keys()))
        for obs in obs_by_cluster[-1]:
            obs_by_cluster[max_cluster + 1] = [obs]
            cluster_by_obs[obs] = max_cluster + 1
            max_cluster += 1

        del obs_by_cluster[-1]

    if annotate:
        adata.obs['pheno_louvain'] = cluster_by_obs

    return cluster_by_obs, obs_by_cluster, graph, q


def _cluster_obs_list_to_dict(cluster_by_obs: List[int]):
    """
    Converts observation clusters in a list to a dict
    mapping clusters to all observations in the cluster

    Parameters
    -----------
    cluster_by_obs: an array of community labels for each cell in adata, in order

    Returns
    -----------
    obs_by_cluster: a map of community labels to lists of cell indices in that community in order
    """
    obs_by_cluster = {}
    for i in range(len(cluster_by_obs)):
        cluster = cluster_by_obs[i]
        if cluster not in obs_by_cluster:
            obs_by_cluster[cluster] = []
        obs_by_cluster[cluster].append(i)

    return obs_by_cluster

def cluster_louvain(
    adata: AnnData,
    k: int=15,
    annotate: bool = False,
    nn_measure: str = 'euclidean',
    knn_method: str = 'pynndescent',
    louvain_method: str = 'taynaud',
    weighting_method: str = 'jaccard',
    annoy_trees: int = None,
    n_jobs: int = 1,
    resolution: float = 1.,
    annoy_index_filename: str = None,
    graph_filename: str = None,
    random_seed: int = None,
    jaccard_prune: float = 0.05,
    jaccard_prune_min_size: int = 50000,
    annoy_seed: int = None,
    knn_seed: int = None
):
    """
    Cluster cells by building an approximate-KNN Jaccard/SNN graph and running
    Louvain or Leiden community detection on it.

    Parameters
    -----------
    adata: an AnnData of data to cluster
    k: number of nearest neighbors to use during graph construction
    annotate: if True, community labels are added to adata as observations. Otherwise
              all scanpy pheongraph outputs are returned.
    nn_measure: metric to use for nearest neighbor evaluation. Can be "angular",
                "euclidean", "manhattan", "hamming", or "dot"
    knn_method: KNN backend.
                "pynndescent" (default) uses NN-descent, which refines an approximate graph
                iteratively. It agrees with scrattch.bigcat slightly better than annoy does
                (ARI 0.840 vs 0.829 on the 428-cluster reference), is exactly reproducible at a fixed
                seed, and is more stable across seeds (0.912 vs 0.875) -- but it builds the index
                entirely in memory.
                "annoy" matches scrattch.bigcat's BiocNeighbors::buildAnnoy: a random-projection
                forest, memory-mapped so it never has to hold the index in RAM. Prefer it when the
                dataset is too large to build in memory, or when byte-level comparability with
                scrattch.bigcat matters. See reports/9_knn_backends.md.
    louvain_method: Louvain method to use, currently can only be "taynaud"
    weighting_method: weighting method to use to use for nearest neighbors graph
                      Can currently be "jaccard" or "uniform".
    annoy_trees: Number of trees to use in annoy random forest index. More trees gives higher
                 precision but is more computationally expensive.
    n_jobs: Number of processes to use in any parallelizable steps
    resolution: Louvain resolution parameter, changes size of communities
    annoy_index_filename: File to store annoy index in, defaults to temp file
    graph_filename: File to store KNN graph AnnData in, if unset does not save graph
    random_seed: int to use as random seed
    knn_seed: seed for the KNN backend, held FIXED and decoupled from random_seed so the graph is
              seed-invariant the way scrattch.bigcat's is (BiocNeighbors::buildAnnoy uses a fixed
              internal seed that ignores set.seed). Applies to whichever backend is selected --
              annoy's index seed, or pynndescent's random_state. Defaults to 1.
    annoy_seed: deprecated alias for knn_seed, kept so existing configs keep working. It was named
                for annoy but has always seeded whichever backend is in use.

    Returns
    -----------
    cluster_by_obs: an array of community labels for each cell in adata, in order
    obs_by_cluster: a map of community labels to lists of cell indices in that community in order
    graph: the calculated adjacency graph on the adata
    q: the maximum modularity of the final clustering 
    """
    # `annoy_seed` predates the pynndescent backend; it seeds whichever backend is selected, so the
    # accurate name is `knn_seed`. The old name still works.
    if annoy_seed is not None:
        warnings.warn("annoy_seed is deprecated; use knn_seed. It seeds whichever KNN backend is "
                      "selected, not only annoy.", DeprecationWarning, stacklevel=2)
        if knn_seed is None:
            knn_seed = annoy_seed
    if knn_seed is None:
        knn_seed = 1

    if knn_method == 'annoy':
        nn_adata = get_annoy_knn(
            adata = adata,
            k = k,
            nn_measure = nn_measure,
            weighting_method = weighting_method,
            annoy_trees = annoy_trees,
            n_jobs = n_jobs,
            annoy_index_filename = annoy_index_filename,
            graph_filename = graph_filename,
            random_seed = knn_seed     # FIXED, decoupled from the clustering seed -> seed-invariant KNN graph, like R
        )
    elif knn_method == 'pynndescent':
        nn_adata = get_pynndescent_knn(
            adata = adata,
            k = k,
            nn_measure = nn_measure,
            weighting_method = weighting_method,
            n_jobs = n_jobs,
            graph_filename = graph_filename,
            random_seed = knn_seed     # same fixed seed slot, so the KNN graph stays seed-invariant
        )
    else:
        raise ValueError(f"{knn_method} is not a valid knn method! "
                         f"Choose 'annoy' (default, matches scrattch.bigcat) or 'pynndescent'.")

    # Match scrattch.bigcat jaccard_big: prune weak jaccard edges (keep weight > prune),
    # but ONLY for large graphs (n > bin_size=50000), which is where R's nbin>1 branch prunes.
    if weighting_method in ('jaccard', 'jaccard_snn') and jaccard_prune > 0 and nn_adata.n_obs > jaccard_prune_min_size:
        Xp = nn_adata.X.tocsr()
        Xp.data[Xp.data <= jaccard_prune] = 0.0
        Xp.eliminate_zeros()
        nn_adata.X = Xp

    if louvain_method == 'taynaud':
        cluster_by_obs, q = get_taynaud_louvain(
            nn_adata = nn_adata,
            resolution = resolution,
            random_seed = random_seed
        )
    elif louvain_method == 'vtraag':
        cluster_by_obs, q = get_vtraag_leiden(
            nn_adata = nn_adata,
            random_seed = random_seed
        )
    else:
        raise ValueError(f"{louvain_method} is not a valid louvain method! Only available methods are taynaud and vtraag")

    obs_by_cluster = _cluster_obs_list_to_dict(cluster_by_obs)

    if annotate:
        adata.obs['pheno_louvain'] = cluster_by_obs

    return cluster_by_obs, obs_by_cluster, nn_adata.X, q

def _uniform_csr_from_nn_dict(nn_dict):
    """
    Converts a dictionary mapping indices to lists of their nearest neighbor indices
    to a csr_matrix with ordered indices and uniform weights with value 1.

    Parameters
    -----------
    nn_dict: dictionary mapping indices to lists of their nearest neighbor indices

    Returns
    -----------
    uniform_csr: csr_matrix object with value 1 at each index (i,j) where j is a
                 nearest neighbor of i
    """
    n = len(nn_dict)
    k = len(nn_dict[0])

    csr_indices = []
    for i in range(n):
        csr_indices.extend(nn_dict.pop(i))

    csr_weights = [1 for i in range(n * k)]
    csr_indptr = [i * k for i in range(n + 1)]

    return csr_matrix((csr_weights, csr_indices, csr_indptr), shape=(n, n), dtype=float)

def _calc_jaccard(idx, max_union, nn_set_dict):
    """
    Given a dictionary mapping indices to sets of their nearest neighbor indices,
    an index, and a constant value provided for computational convenience computes the
    jaccard coefficients for the nearest neighbors of observation idx. Returns both the index
    and the coefficients for post-parallelism re-ordering.

    Parameters
    -----------
    idx: index of the observation to evaluate
    max_union: 2k, provided for computational brevity
    nn_dict: dictionary mapping indices to sets of their nearest neighbor indices

    Returns
    -----------
    idx: index of the evaluated observation
    nn_jaccard: list of jaccard coefficients for the observation's nearest neighbors
    """
    shared_neighbors = [float(len(nn_set_dict[idx] & nn_set_dict[j])) for j in nn_set_dict[idx]]
    return idx, [x / (max_union - x) for x in shared_neighbors]

def _jaccard_csr_from_nn_dict(nn_dict, n_jobs = 1):
    """
    Converts a dictionary mapping indices to lists of their nearest neighbor indices
    to a csr_matrix with ordered indices and jaccard coefficients as weights.

    Parameters
    -----------
    nn_dict: dictionary mapping indices to lists of their nearest neighbor indices
    n_joba: number of jobs to use in computing jaccard coefficients. Defaults to 1.

    Returns
    -----------
    jaccard_csr: csr_matrix object with jaccard coefficient of observations i and j at
                 index (i,j) where j is a nearest neighbor of i
    """
    n = len(nn_dict)
    k = len(nn_dict[0])

    pool = Pool(processes=n_jobs)
    weight_map = {}
    chunk_size = min(max(1, int(n / n_jobs)), 10000)
    nn_set_dict = {k: set(v) for k,v in nn_dict.items()}
    picklable_jaccard = functools.partial(_calc_jaccard, max_union = 2*k, nn_set_dict = nn_set_dict)
    for chunk_result in tqdm(pool.imap_unordered(picklable_jaccard, range(n), chunksize=chunk_size), total=n):
        weight_map[chunk_result[0]] = chunk_result[1]

    pool.close()

    csr_weights = []
    csr_indices = []
    for i in range(n):
        csr_weights.extend(weight_map[i])
        csr_indices.extend(nn_dict[i])

    csr_indptr = [i * k for i in range(n)] + [n * k]
    return csr_matrix((csr_weights, csr_indices, csr_indptr), shape=(n, n), dtype=float)


def _jaccard_snn_csr_from_nn_dict(nn_dict):
    """
    Full shared-nearest-neighbor (SNN) Jaccard graph, matching R scrattch.bigcat::jaccard_big.
    Builds the binary KNN membership matrix B (cell x cell, B[i,j]=1 if j in KNN(i)), forms the
    co-neighbor counts A = B @ B.T (A[i,j] = # neighbors shared by cells i and j), and weights every
    such pair by Jaccard = shared / (2k - shared). Unlike _jaccard_csr_from_nn_dict (edges only on
    direct KNN pairs), this creates an edge between EVERY pair of cells that shares >=1 neighbor.
    """
    n = len(nn_dict)
    k = len(nn_dict[0])
    rows = np.repeat(np.arange(n), k)
    cols = np.fromiter((j for i in range(n) for j in nn_dict[i]), dtype=np.int64, count=n * k)
    B = csr_matrix((np.ones(n * k, dtype=np.float32), (rows, cols)), shape=(n, n))
    A = (B @ B.T).tocoo()                       # shared-neighbor counts (all co-neighbor pairs)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = A.data / (2.0 * k - A.data)         # Jaccard
    return csr_matrix((w.astype(np.float64), (A.row, A.col)), shape=(n, n))


# Per-worker annoy index cache. The index is immutable and identical for every query, so each
# pool worker memory-maps it once (via _init_annoy_worker) instead of re-opening it per observation.
_ANNOY_WORKER = {}


def _init_annoy_worker(vec_len, nn_measure, annoy_index_filename):
    """
    Pool initializer: memory-map the annoy index once per worker process.

    Parameters
    -----------
    vec_len: length of the observation vectors
    nn_measure: distance metric used to evaluate nearest neighbors on this index
    annoy_index_filename: location of the file-cached annoy index
    """
    ai = AnnoyIndex(vec_len, nn_measure)
    ai.load(annoy_index_filename)
    _ANNOY_WORKER['index'] = ai


def _search_nn_chunk(idx_block, k):
    """
    Find the k nearest neighbors for a block of observations, querying the annoy index
    that _init_annoy_worker loaded once for this worker process.

    Parameters
    -----------
    idx_block: iterable of observation indices to find nearest neighbors for
    k: number of nearest neighbors to find

    Returns
    -----------
    list of (idx, nn_list) pairs, where nn_list holds the indices of the k nearest
    neighbors to the observation at idx
    """
    ai = _ANNOY_WORKER['index']
    return [(idx, ai.get_nns_by_item(idx, k)) for idx in idx_block]

def _annoy_build_csr_nn_graph(
    data_matrix: np.array,
    annoy_index_filename: str,
    k: int,
    n_jobs: int = 1,
    nn_measure: str = 'euclidean',
    weighting_method: str = 'jaccard'
):
    """
    Given a matrix of observation data and the location of a file-cached annoy index on
    that data, generates a csr_matrix of the k nearest neighbors of each observation weighted
    by the method specified in weighting.

    Parameters
    -----------
    data_matrix: a matrix of observation data
    annoy_index_filename: location of the file-cached annoy index
    k: number of nearest neighbors to find
    n_jobs: number of processes to use while computing nearest neighbors
    nn_measure: distance metric used to calculate nearest neighbors
    weighting_method: weighting method to use to use for nearest neighbors graph.
                      Can currently be "jaccard" or "uniform".

    Returns
    -----------
    nn_csr: csr_matrix of weighted neaerest neighbors
    """
    n, vec_len = data_matrix.shape
    graph_map = {}
    chunk_size = min(max(1, int(n / n_jobs)), 10000)
    # Dispatch blocks of indices; each worker loads the index once (initializer) and queries the block.
    idx_blocks = [range(start, min(start + chunk_size, n)) for start in range(0, n, chunk_size)]
    picklable_search = functools.partial(_search_nn_chunk, k=k)
    with Pool(processes=n_jobs,
              initializer=_init_annoy_worker,
              initargs=(vec_len, nn_measure, annoy_index_filename)) as pool:
        for block_result in tqdm(pool.imap_unordered(picklable_search, idx_blocks), total=len(idx_blocks)):
            for idx, nn_list in block_result:
                graph_map[idx] = nn_list

    if weighting_method == 'jaccard':
        return _jaccard_csr_from_nn_dict(graph_map, n_jobs)
    if weighting_method == 'jaccard_snn':
        return _jaccard_snn_csr_from_nn_dict(graph_map)
    if weighting_method == 'uniform':
        return _uniform_csr_from_nn_dict(graph_map)
    else:
        raise ValueError(f"{weighting_method} is not a valid weighting option! Must use jaccard, jaccard_snn or uniform")

def get_pynndescent_knn(
    adata: AnnData,
    k: int,
    nn_measure: str = 'euclidean',
    weighting_method: str = 'jaccard',
    n_jobs: int = 1,
    graph_filename: str = None,
    random_seed: int = None
):
    """
    K nearest neighbors via pynndescent (NN-descent), as an alternative to the annoy backend.

    Produces exactly the same {index: neighbour list} structure annoy's path does -- k neighbours
    per cell, self included -- so the graph weighting (jaccard / jaccard_snn / uniform) and
    everything downstream is shared and unchanged.

    Where the two differ: annoy builds a random-projection forest and never refines it, and can be
    memory-mapped from disk; NN-descent starts from a random graph and iteratively improves it, so
    it is usually more accurate at the same k but holds everything in memory. See
    reports/knn_backends.md.

    Parameters
    -----------
    adata: an AnnData of data to cluster
    k: number of nearest neighbors to find (self included, matching the annoy path)
    nn_measure: distance metric; annoy's "angular" is mapped to pynndescent's "cosine"
    weighting_method: "jaccard", "jaccard_snn" or "uniform"
    n_jobs: number of threads for index construction and querying
    graph_filename: File to store KNN graph AnnData in, if unset does not save graph
    random_seed: seed for NN-descent's random initialisation, for a reproducible graph

    Returns
    -----------
    nn_adata: An AnnData object with a weighted nearest neighbor graph of the observations
    """
    try:
        from pynndescent import NNDescent
    except ImportError as err:
        raise ImportError("pynndescent is not installed. It is the default KNN backend "
                          "(pip install pynndescent); or pass knn_method='annoy' to use the "
                          "memory-mapped annoy backend instead.") from err

    # pynndescent calls numba.set_num_threads(n_jobs), and numba raises if n_jobs exceeds
    # NUMBA_NUM_THREADS (fixed at import from the detected core count). The annoy backend uses a
    # multiprocessing Pool and has no such limit, so the same n_jobs that is fine there would abort
    # here -- clamp instead of making the caller special-case the backend.
    try:
        import numba
        max_threads = int(numba.config.NUMBA_NUM_THREADS)
        if n_jobs > max_threads:
            warnings.warn(f"pynndescent: n_jobs={n_jobs} exceeds numba's thread cap "
                          f"({max_threads}); using {max_threads}")
            n_jobs = max_threads
    except ImportError:
        pass
    n_jobs = max(1, n_jobs)

    data_matrix = adata.X
    if issparse(data_matrix):
        data_matrix = data_matrix.toarray()
    data_matrix = np.asarray(data_matrix, dtype=np.float32)

    # annoy metric names -> pynndescent metric names
    metric = {'angular': 'cosine', 'dot': 'dot'}.get(nn_measure, nn_measure)

    # The recursive pipeline drills down to very small subsets, where k can exceed the number of
    # points. annoy silently returns min(k, n) neighbours; pynndescent needs n_neighbors <= n.
    n_obs = data_matrix.shape[0]
    k_eff = max(2, min(k, n_obs))

    index = NNDescent(
        data_matrix,
        n_neighbors=k_eff,
        metric=metric,
        random_state=random_seed,
        n_jobs=n_jobs,
    )
    neighbors, _ = index.neighbor_graph          # (n, k_eff) -- includes the point itself

    # pynndescent pads unfilled slots with -1 (it happens on small or sparsely-connected subsets).
    # Those would reach csr_matrix as negative indices and abort graph construction, so fold them
    # onto the point itself -- which is what annoy effectively yields and leaves the fixed-width
    # Jaccard/SNN construction downstream valid.
    if (neighbors < 0).any():
        neighbors = np.where(neighbors < 0, np.arange(n_obs, dtype=neighbors.dtype)[:, None],
                             neighbors)
    graph_map = {i: list(map(int, neighbors[i])) for i in range(neighbors.shape[0])}

    if weighting_method == 'jaccard':
        csr_graph = _jaccard_csr_from_nn_dict(graph_map, n_jobs)
    elif weighting_method == 'jaccard_snn':
        csr_graph = _jaccard_snn_csr_from_nn_dict(graph_map)
    elif weighting_method == 'uniform':
        csr_graph = _uniform_csr_from_nn_dict(graph_map)
    else:
        raise ValueError(f"{weighting_method} is not a valid weighting option! "
                         f"Must use jaccard, jaccard_snn or uniform")

    graph_adata = AnnData(csr_graph, obs=adata.obs, var=adata.obs)
    if graph_filename:
        graph_adata.write(graph_filename)
    return graph_adata


def get_annoy_knn(
    adata: AnnData,
    k: int,
    nn_measure: str = 'euclidean',
    weighting_method: str = 'jaccard',
    annoy_trees: int = None,
    n_jobs: int = 1,
    annoy_index_filename: str = None,
    graph_filename: str = None,
    random_seed: int = None
):
    """
    Given an AnnData object of observation data, uses the annoy package to compute the k neareset neighbors
    of each observation. Returns an AnnData object with a nearest neighbor graph of the observations weighted
    by the method specified in weighting_method.

    Parameters
    -----------
    adata: an AnnData of data to cluster
    k: number of nearest neighbors to use during graph construction
    nn_measure: metric to use for nearest neighbor evaluation. Can be "angular",
                "euclidean", "manhattan", "hamming", or "dot"
    weighting_method: weighting method to use to use for nearest neighbors graph
                      Can currently be "jaccard" or "uniform".
    annoy_trees: Number of trees to use in annoy random forest index. More trees gives higher
                 precision but is more computationally expensive.
    n_jobs: Number of processes to use in any parallelizable steps.
            Note that AnnoyIndex builds differently on different n_jobs values
    annoy_index_filename: File to store annoy index in, defaults to temp file
    graph_filename: File to store KNN graph AnnData in, if unset does not save graph
    random_seed: int to use as random seed
    knn_seed: seed for the KNN backend, held FIXED and decoupled from random_seed so the graph is
              seed-invariant the way scrattch.bigcat's is (BiocNeighbors::buildAnnoy uses a fixed
              internal seed that ignores set.seed). Applies to whichever backend is selected --
              annoy's index seed, or pynndescent's random_state. Defaults to 1.
    annoy_seed: deprecated alias for knn_seed, kept so existing configs keep working. It was named
                for annoy but has always seeded whichever backend is in use.

    Returns
    -----------
    nn_adata: An AnnData object with a weighted nearest neighbor graph of the observations
    """
    data_matrix = adata.X
    ai = AnnoyIndex(data_matrix.shape[1], nn_measure)
    if annoy_index_filename and os.path.isfile(annoy_index_filename):
        ai.load(annoy_index_filename)
    else:
        if not annoy_index_filename:
            annoy_index_filename = os.path.join(tempfile.gettempdir(), f"annoy_index_{datetime.now().strftime('%Y%m%d%H%M%S')}")
        ai.on_disk_build(annoy_index_filename)
        if random_seed:
            ai.set_seed(random_seed)

        for i in range(data_matrix.shape[0]):
            ai.add_item(i, data_matrix[i])

        if not annoy_trees:
            annoy_trees = max(1, int(log(data_matrix.shape[0], 2)))

        ai.build(annoy_trees, n_jobs=n_jobs)

    csr_graph = _annoy_build_csr_nn_graph(data_matrix, annoy_index_filename, k, n_jobs, nn_measure, weighting_method)

    graph_adata = AnnData(csr_graph, obs=adata.obs, var=adata.obs)
    if graph_filename:
        graph_adata.write(graph_filename)
    return graph_adata

def get_taynaud_louvain(
    nn_adata: AnnData,
    resolution: float = 1.,
    random_seed: int = None
):
    """
    Given an AnnData object containing a weighted nearest neighbor graph, computes the
    louvain partition of the graph.

    Parameters
    -----------
    nn_adata: An AnnData object with a weighted nearest neighbor graph of the observations
    resolution: Louvain resolution parameter, changes size of communities
    random_seed: int to use as random seed 

    Returns
    -----------
    cluster_by_obs: an array of community labels for each cell in adata, in order
    q: the maximum modularity of the final clustering 
    """
    nn_coo = nn_adata.X.tocoo()
    G = nx.Graph()
    G.add_weighted_edges_from(list(zip(nn_coo.row, nn_coo.col, nn_coo.data)))

    partition = community_louvain.best_partition(G, resolution=resolution, random_state=random_seed)
    cluster_by_obs = [partition[i] for i in range(len(partition))]

    q = community_louvain.modularity(partition, G)

    return cluster_by_obs, q

def get_vtraag_leiden(
    nn_adata: AnnData,
    random_seed: int = None
):
    """
    Given an AnnData object containing a weighted nearest neighbor graph, computes the
    louvain partition of the graph. This method is faster than taynaud, but not necessarily
    as optimal.

    Parameters
    -----------
    nn_adata: An AnnData object with a weighted nearest neighbor graph of the observations
    random_seed: int to use as random seed 

    Returns
    -----------
    cluster_by_obs: an array of community labels for each cell in adata, in order
    q: the maximum modularity of the final clustering 
    """
    nn_coo = nn_adata.X.tocoo()
    nn_igraph = ig.Graph(edges=list(zip(nn_coo.row, nn_coo.col)), edge_attrs={'weights': nn_coo.data})
    optimiser = leidenalg.Optimiser()
    if random_seed:
        optimiser.set_rng_seed(random_seed)
    partition = leidenalg.ModularityVertexPartition(nn_igraph, weights='weights')
    optimiser.optimise_partition(partition)
    return partition.membership, partition.quality()
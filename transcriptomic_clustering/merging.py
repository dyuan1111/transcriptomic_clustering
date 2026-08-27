from typing import Any, Tuple, Dict, List, Optional, Set, Literal
import time
import anndata as ad
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import logging
from collections import defaultdict
import warnings
import transcriptomic_clustering as tc
from transcriptomic_clustering.markers import select_marker_genes
from transcriptomic_clustering.de_all_pairs import _de_one_pair, SCORE_CAP, _FILTER_KEYS
from transcriptomic_clustering.de_ebayes import get_linear_fit_vals, moderate_variances

logger = logging.getLogger(__name__)

DEFAULT_THRESHOLDS = {
    'q1_thresh': 0.5,
    'q2_thresh': None,
    'cluster_size_thresh': 6,
    'qdiff_thresh': 0.7,
    'padj_thresh': 0.05,
    'lfc_thresh': 1.0,
    'score_thresh': 150,
    'low_thresh': 1,
    'min_genes': 5
}

# Genes per direction per batch carried into the batch-aware conservation filter.
# scrattch.bigcat harmonize_merge.R:126 (`top.n=1000`).
CONSERVATION_TOP_N = 1000


def merge_clusters(
        adata_norm: ad.AnnData,
        adata_reduced: ad.AnnData,
        cluster_assignments: Dict[Any, List],
        cluster_by_obs: np.ndarray,
        thresholds: Dict[str, Any] = DEFAULT_THRESHOLDS,
        k: Optional[int] = 2,
        de_method: Optional[str] = 'ebayes',
        n_markers: Optional[int] = 20,
        chunk_size: Optional[int] = None,
        return_markers_df: Optional[bool] = False,
        n_jobs: Optional[int] = 1,
        batch_aware_merging: Optional[Dict[str, Any]] = None
) -> Tuple[Dict[Any, np.ndarray], Set]:
    """
    Merge clusters based on size and differential gene expression score

    1. merge small clusters
    2. merge clusters based by differential gene expression score

    Parameters
    ----------
    adata_norm:
        AnnData object of normalized data
    adata_reduced:
        AnnData object in reduced space
    cluster_assignments:
        map of cluster label to cell idx belonging to cluster
    cluster_by_obs:
        array of cells with cluster value
    thresholds:
        threshold for de calculation
    k:
        number of cluster neighbors
    de_method:
        method used for de calculation
    n_markers:
        number of marker genes to select from up and down de genes
        if 0 or None will skip
    chunk_size:
        number of observations to process in a single chunk
    return_markers_df:
        return markers as a dataframe
    n_jobs:
        number of jobs to run in parallel
    batch_aware_merging:
        None (default) merges on pooled statistics, exactly as before. A dict switches to
        batch-aware merging (`merge_clusters_by_de_batch_aware`), which tests each batch separately
        and merges only when none objects -- use it when cells come from several platforms or
        modalities and a pooled difference could be a batch artefact. Keys:

        - ``batch_obs`` (required): column of ``adata_norm.obs`` naming each cell's batch
        - ``lfc_conservation_th`` (default 0.7): fraction of measuring batches that must agree on a
          gene's fold-change direction for it to count as evidence
        - ``conservation_lfc_th`` (default: ``thresholds['lfc_thresh']``): the fold change a gene
          must clear in a batch to count as agreeing
        - ``thresholds`` (optional): ``{batch: {threshold: value}}`` overrides layered on the shared
          ``thresholds``; unlisted batches inherit them unchanged. Any threshold can be overridden.
          The validated real-world case is WMB2 production giving single-nucleus data a lower
          ``q1_thresh`` (0.3 vs 0.4) because nuclei sample less RNA and so detect any gene in fewer
          cells; a hypothetical of the same kind is a lower ``score_thresh`` for a shallow gene
          panel that cannot accumulate whole-transcriptome DE scores.
        - ``genes_by_batch`` (optional): ``{batch: genes that batch's reference measures}``, for
          multi-platform data whose references carry different gene lists. Each batch's DE runs on
          its own genes, and a gene measured by only ONE platform still counts as evidence (it is
          trivially conserved), as in R. Omit when every batch measures every gene.
        - ``pair_batch`` (optional, default None): None tests every shortlisted pair each round (the
          deterministic reference). An integer N mirrors R's ``pairBatch`` (100 in production): test
          N pairs at a time by decreasing centroid similarity and stop the round's testing as soon
          as anything is mergeable. Faster, but the merge order then depends on scan order.

        Requires ``de_method='ebayes'`` and an in-memory ``adata_norm``.

    Returns
    -------
    cluster_assignments:
        updated mapping of cluster assignments
    markers:
        if select_markers, set of genes that differentiate clusters, or empty set
    """
    if len(cluster_assignments.keys()) == 1:
        return cluster_assignments.copy()

    # Validate the batch-aware options up front -- a typo here should not surface after the
    # (expensive) cluster-mean passes, and must never fall back silently to pooled merging.
    batch_by_obs = None
    if batch_aware_merging is not None:
        allowed_keys = {'batch_obs', 'lfc_conservation_th', 'conservation_lfc_th', 'thresholds',
                        'pair_batch', 'genes_by_batch'}
        unknown_keys = set(batch_aware_merging) - allowed_keys
        if unknown_keys:
            raise ValueError(
                f"batch_aware_merging has unknown keys {sorted(unknown_keys)}; "
                f"expected {sorted(allowed_keys)}"
            )
        _pb = batch_aware_merging.get('pair_batch')
        if _pb is not None and (not isinstance(_pb, int) or _pb < 1):
            raise ValueError(f"batch_aware_merging['pair_batch'] must be a positive integer or None, got {_pb!r}")
        if 'batch_obs' not in batch_aware_merging:
            raise ValueError("batch_aware_merging requires 'batch_obs', naming a column of adata_norm.obs")
        batch_obs = batch_aware_merging['batch_obs']
        if batch_obs not in adata_norm.obs:
            raise ValueError(
                f"batch_aware_merging['batch_obs']={batch_obs!r} is not a column of adata_norm.obs. "
                f"Available: {list(adata_norm.obs.columns)}"
            )
        if de_method != 'ebayes':
            # The conservation filter matches genes by name across batches; the chisq path reports
            # positions rather than names, so it cannot supply the per-gene evidence this needs.
            raise ValueError(
                f"batch-aware merging requires de_method='ebayes', got {de_method!r}"
            )
        batch_by_obs = adata_norm.obs[batch_obs].values

    # Calculate cluster means on reduced space
    logger.info("Computing reduced cluster means")
    tic = time.perf_counter()
    cl_means_reduced, _, _ = tc.get_cluster_means(adata_reduced,
                                               cluster_assignments,
                                               cluster_by_obs,
                                               chunk_size,
                                               low_th=thresholds['low_thresh'])
    logger.info(f'Completed reduced cluster means')
    toc = time.perf_counter()
    logger.info(f'Reduced Cluster Means Elapsed Time: {toc - tic}')


    # Small-cluster handling differs by merge mode, matching the two R pipelines:
    #   pooled       -> merge_cl_big (merge_cl.R:97-125): a cluster below cluster_size_thresh POOLED
    #                   cells is merged WHOLE into its nearest reduced-space centroid.
    #   batch-aware  -> merge_cl_multiple (harmonize_merge.R:406-495): a cluster below min.cells in
    #                   EVERY batch is DISSOLVED -- its cells re-mapped individually to the nearest
    #                   big-cluster centroid (cosine, as map_cells_knn's Annoy.Cosine default), so
    #                   they may scatter. Without this, such clusters are unjudgeable by every batch
    #                   and -- since unjudgeable pairs never merge -- would survive to the output
    #                   without any evidence ever vetting them.
    cluster_assignments_merge = cluster_assignments.copy()
    tic = time.perf_counter()
    batch_thresholds = None
    if batch_aware_merging is None:
        logger.info("Merging small clusters")
        merge_small_clusters(cl_means_reduced, cluster_assignments_merge, thresholds['cluster_size_thresh'])
        logger.info(f'Completed merging small clusters')
    else:
        batch_thresholds = _resolve_batch_thresholds(
            thresholds, batch_aware_merging.get('thresholds'),
            list(pd.unique(np.asarray(batch_by_obs)))
        )
        dissolved = _dissolve_cl_small(
            cluster_assignments_merge, cl_means_reduced, adata_reduced, batch_by_obs,
            {b: t['cluster_size_thresh'] for b, t in batch_thresholds.items()},
        )
        logger.info(f'Dissolved {len(dissolved)} cl.small clusters '
                    f'(below min.cells in every batch): {sorted(map(str, dissolved))}')
    toc = time.perf_counter()
    logger.info(f'Small Clusters Elapsed Time: {toc - tic}')

    # Create new cluster_by_obs based on updated cluster assignments
    cluster_by_obs = np.zeros((adata_norm.shape[0],))
    for cl_id, idxs in cluster_assignments_merge.items():
        cluster_by_obs[idxs] = cl_id

    # The dissolve moved cells between clusters, so the reduced-space centroids that seed the
    # candidate-pair shortlist must be recomputed -- R does the same (cl.rd is rebuilt from the
    # updated cl at harmonize_merge.R:508 before the merge loop starts).
    if batch_aware_merging is not None:
        cl_means_reduced, _, _ = tc.get_cluster_means(adata_reduced,
                                                      cluster_assignments_merge,
                                                      cluster_by_obs,
                                                      chunk_size,
                                                      low_th=thresholds['low_thresh'])

    # Calculate cluster means on normalized data
    logger.info("Computing Cluster Means")
    tic = time.perf_counter()
    cl_means, present_cl_means, cl_vars = tc.get_cluster_means(adata_norm,
                                                      cluster_assignments_merge,
                                                      cluster_by_obs,
                                                      chunk_size,
                                                      low_th=thresholds['low_thresh'])
    logger.info(f'Completed Cluster Means')
    toc = time.perf_counter()
    logger.info(f'Cluster Means Elapsed Time: {toc - tic}')

    # Merge remaining clusters by differential expression
    logger.info("Merging Clusters by DE")
    tic = time.perf_counter()
    # merge_clusters_by_de mutates cluster_assignments_merge in place AND returns the merged
    # means/vars/present (rebuilt from its numpy-level updates) for downstream marker selection.
    if batch_aware_merging is None:
        cl_means, cl_vars, present_cl_means = merge_clusters_by_de(cluster_assignments_merge,
                             cl_means,
                             cl_vars,
                             present_cl_means,
                             cl_means_reduced,
                             thresholds,
                             k,
                             de_method,
                             )
    else:
        logger.info("Computing per-batch cluster means")
        tic_b = time.perf_counter()
        stats_by_batch, sizes_by_batch = tc.get_cluster_means_per_batch(
            adata_norm, cluster_assignments_merge, batch_by_obs,
            low_th=thresholds['low_thresh'],
        )
        logger.info(f'Per-batch Cluster Means Elapsed Time: {time.perf_counter() - tic_b}')
        # batch_thresholds was resolved before the cl.small dissolve; batches with no cells in any
        # cluster cannot occur here (every batch has cells, every cell has a cluster)
        cl_means, cl_vars, present_cl_means = merge_clusters_by_de_batch_aware(
            cluster_assignments_merge,
            cl_means,
            cl_vars,
            present_cl_means,
            stats_by_batch,
            sizes_by_batch,
            cl_means_reduced,
            thresholds,
            batch_thresholds,
            batch_aware_merging.get('lfc_conservation_th', 0.7),
            batch_aware_merging.get('conservation_lfc_th'),
            k,
            batch_aware_merging.get('pair_batch'),
            batch_aware_merging.get('genes_by_batch'),
        )
    logger.info(f'Completed Merging Clusters by DE')
    toc = time.perf_counter()
    logger.info(f'Merging DE Elapsed Time: {toc - tic}')

    # Select marker genes
    if return_markers_df is False and n_markers is None:
        markers = None
    else:
        logger.info('Starting Marker Selection')
        tic = time.perf_counter()
        markers = select_marker_genes(
            cluster_assignments=cluster_assignments_merge,
            cluster_means=cl_means,
            cluster_variances=cl_vars,
            present_cluster_means=present_cl_means,
            # strip merge-loop-only controls; they are not DE-filter thresholds and would be
            # forwarded to filter_gene_stats (which rejects unknown kwargs).
            thresholds={k: v for k, v in thresholds.items() if k not in ('merge_mode', 'max_cl_size')},
            n_markers=n_markers,
            de_method=de_method,
            return_markers_df=return_markers_df,
            n_jobs=n_jobs
        )
        logger.info('Completed Marker Selection')
        toc = time.perf_counter()
        logger.info(f'Marker Selection Elapsed Time: {toc - tic}')


    return cluster_assignments_merge, markers


def merge_two_clusters(
        cluster_assignments: Dict[Any, List],
        label_source: Any,
        label_dest: Any,
        cluster_means: pd.DataFrame,
        cluster_variances: Optional[pd.DataFrame]=None,
        present_cluster_means: pd.DataFrame=None
):
    """
    Merge source cluster into a destination cluster by:
    1. updating cluster means and variances
    2. updating mean of expressions present if not None
    3. updating cluster assignments

    Parameters
    ----------
    cluster_assignments:
        map of cluster label to cell idx belonging to cluster
    label_source:
        label of cluster being merged
    label_dest:
        label of cluster merged into
    cluster_means:
        dataframe of cluster means indexed by cluster label
    cluster_variances:
        dataframe of cluster variances indexed by cluster label
    present_cluster_means:
        dataframe of cluster means indexed by cluster label filtered by low_th

    Returns
    -------
    """

    merge_cluster_means_vars(cluster_assignments, label_source, label_dest, cluster_means, cluster_variances)

    if present_cluster_means is not None:
        merge_cluster_means_vars(cluster_assignments, label_source, label_dest, present_cluster_means, None)

    # merge cluster assignments
    cluster_assignments[label_dest].extend(cluster_assignments[label_source])
    cluster_assignments.pop(label_source)


def merge_cluster_means_vars(
        cluster_assignments: Dict[Any, List],
        label_source: Any,
        label_dest: Any,
        cluster_means: pd.DataFrame,
        cluster_variances: Optional[pd.DataFrame]
):
    """
    Merge source cluster into a destination cluster by:
    1. computing the updated cluster centroid (mean gene expression)
        of the destination cluster
    2. compute updated variance
    3. deleting source cluster after merged

    Parameters
    ----------
    cluster_assignments:
        map of cluster label to cell idx belonging to cluster
    label_source:
        label of cluster being merged
    label_dest:
        label of cluster merged into
    cluster_means:
        dataframe of cluster means indexed by cluster label
    cluster_variances:
        dataframe of cluster variances indexed by cluster label

    Returns
    -------
    """

    # update cluster means:
    n1 = len(cluster_assignments[label_source])
    n2 = len(cluster_assignments[label_dest])
    mean1 = cluster_means.loc[label_source]
    mean2 = cluster_means.loc[label_dest]

    mean_comb = (mean1 * n1 + mean2 * n2) / (n1 + n2)

    cluster_means.loc[label_dest] = mean_comb
    cluster_means.drop(label_source, inplace=True)

    if cluster_variances is not None:
        var1 = cluster_variances.loc[label_source]
        var2 = cluster_variances.loc[label_dest]
        
        var_comb = 1 / (n1 + n2 - 1) * (
            (n1 - 1) * var1 + n1 * (mean1 - mean_comb) ** 2 +
            (n2 - 1) * var2 + n1 * (mean2 - mean_comb) ** 2
        )

        cluster_variances.loc[label_dest] = var_comb
        cluster_variances.drop(label_source, inplace=True)


def cdist_normalized(
        X: np.ndarray,
        Y: np.ndarray,
) -> np.ndarray:
    """
    Calculate similarity metric as (1 - pairwise_distance/max_distance)

    Parameters
    ----------
    X:
        An m by n array of m clusters in an n-dimensional space
    Returns
    -------
    similarity:
        measure of similarity
    """
    similarity = cdist(X, Y, 'euclidean')
    similarity /= np.max(similarity)
    similarity *= -1
    similarity += 1

    return similarity


def calculate_similarity(
        cluster_means: pd.DataFrame,
        group_rows: List[Any],
        group_cols: List[Any]
) -> pd.DataFrame:
    """
    Calculate similarity measure between two cluster groups (group_rows, group_cols)
    based on cluster means (cluster centroids)
    If data has more than 2 dimensions use correlation coefficient as a measure of similarity
    else use normalized distance measure

    Parameters
    ----------
    cluster_means:
        Dataframe of cluster means with cluster labels as index
    group_rows:
        cluster group with clusters being merged
    group_cols:
        cluster group with destination clusters
    Returns
    -------
    similarity:
        array of similarity measure
    """
    source_means = cluster_means.loc[group_rows]
    destination_means = cluster_means.loc[group_cols]

    # #4 (bigcat alignment): use EUCLIDEAN distance on the reduced-space cluster means, matching R
    # scrattch.bigcat get_knn_pairs (method="Annoy.Euclidean", sim = 1 - dist/max(dist) == cdist_normalized).
    # The stock code used correlation distance for >2 dims, which fed the merge a DIFFERENT candidate-pair
    # set than R. Euclidean here matches R's nearest-cluster selection for both the DE-merge candidates
    # and the small-cluster pre-merge. (Ordering is invariant to the global 1-dist/max normalization.)
    similarity = cdist_normalized(source_means, destination_means)

    similarity_df = pd.DataFrame(
        similarity,
        index=group_rows,
        columns=group_cols,
        copy=False
    )
    for common in (set(group_rows) & set(group_cols)):
        similarity_df.at[common, common] = np.nan

    return similarity_df


def find_most_similar(
        similarity_df: pd.DataFrame,
) -> Tuple[Any, Any, float]:
    """

    Parameters
    ----------
    similarity_df:
        similarity metric between clusters
    Returns
    -------
    source_label, dest_label, max_similarity:
        labels of the source and destination clusters and their similarity value
    """

    similarity_df = similarity_df.transpose()

    similarity_sorted = similarity_df.unstack().sort_values(ascending=False).dropna()
    source_label, dest_label = similarity_sorted.index[0]
    max_similarity = similarity_sorted[(source_label, dest_label)]

    return source_label, dest_label, max_similarity


def find_small_clusters(
        cluster_assignments: Dict[Any, List],
        min_size: int
) -> List[Any]:

    return [k for (k, v) in cluster_assignments.items() if len(v) < min_size]


def merge_small_clusters(
        cluster_means: pd.DataFrame,
        cluster_assignments: Dict[Any, List],
        min_size: int,
):
    """
    Then merge small clusters (with size < min_size) iteratively as:

    1. calculate similarity between small and all clusters
    2. merge most-highly similar small cluster
    3. update list of small/all clusters
    4. go to 1 until all small clusters are merged

    Parameters
    ----------
    cluster_means:
        dataframe of cluster means indexed by cluster label
    cluster_assignments:
        map of cluster label to cell idx belonging to cluster
    min_size:
        smallest size that is not merged

    Returns
    -------
    cluster_assignments:
        updated mapping of cluster assignments
    """
    all_cluster_labels = list(cluster_assignments.keys())
    small_cluster_labels = find_small_clusters(cluster_assignments, min_size=min_size)

    while small_cluster_labels:
        if len(cluster_assignments.keys()) == 1:
            break

        similarity_small_to_all_df = calculate_similarity(
            cluster_means,
            group_rows=small_cluster_labels,
            group_cols=all_cluster_labels)

        source_label, dest_label, max_similarity = find_most_similar(
            similarity_small_to_all_df,
        )
        logger.debug(f"Merging small cluster {source_label} into {dest_label} -- similarity: {max_similarity}")
        merge_two_clusters(cluster_assignments, source_label, dest_label, cluster_means)

        # update labels:
        small_cluster_labels = find_small_clusters(cluster_assignments, min_size=min_size)
        all_cluster_labels = list(cluster_assignments.keys())



def merge_clusters_by_de(
    cluster_assignments: Dict[Any, List],
    cluster_means: pd.DataFrame,
    cluster_variances: pd.DataFrame,
    present_cluster_means: pd.DataFrame,
    cluster_means_rd: pd.DataFrame,
    thresholds: Dict[str, Any],
    k: Optional[int] = 2,
    de_method: Optional[Literal['ebayes', 'chisq']] = 'ebayes',
):
    """
    Merge clusters by the calculated gene differential expression score

    1. get k nearest clusters for each cluster in a reduced space
    2. calculate differential expression scores for all pairs
    3. sort scores by lowest and loop through them, merging pairs with scores lower than threshold

    Parameters
    ----------
    cluster_assignments:
        map of cluster label to cell idx belonging to cluster
    cluster_means:
        dataframe of cluster means indexed by cluster label in a normalized space
    cluster_variances:
        dataframe of cluster variances indexed by cluster label in a normalized space
    present_cluster_means:
        dataframe of cluster means indexed by cluster label filtered by low_th in a normalized space
    cluster_means_rd:
        dataframe of cluster means indexed by cluster label in a reduced space
    k:
        number of cluster neighbors
    de_method:
        method used for de calculation
    thresholds:
        threshold use de calculation

    Returns
    -------
    cluster_assignments:
        updated mapping of cluster assignments
    """
    cl_size = {k: len(v) for k, v in cluster_assignments.items()}

    thresholds = thresholds.copy()
    score_th = thresholds.pop('score_thresh')
    min_genes = thresholds.pop('min_genes')
    thresholds.pop('low_thresh')
    # max_cl_size caps the cells/cluster used in the DE test (matching scrattch.bigcat max.cl.size).
    # DISABLED by default (None) -> use ALL cells per cluster, matching R with max.cl.size=Inf.
    # Set an integer in merge_clusters_kwargs.thresholds['max_cl_size'] to re-enable the cap.
    max_cl_size = thresholds.pop('max_cl_size', None)
    # merge_mode:
    #   'aligned' (default) - R-faithful: re-examine ALL clusters' k-nearest each round, merge the single
    #      lowest-score pair + extras only < score_th/2, recompute. Highest fidelity to scrattch.bigcat,
    #      but slow (many rounds).
    #   'fast' - previous behavior: after round 1 only re-examine merged-destination clusters' neighbors,
    #      and merge EVERY non-conflicting candidate < score_th per round. Far fewer rounds -> much faster.
    merge_mode = thresholds.pop('merge_mode', 'aligned')

    # ---- numpy state for the wide (clusters x genes) mean/var/present matrices ----
    # Per-merge updates on 17k-column DataFrames via .loc[row]= / .drop are pathologically slow: each
    # assignment re-introspects every column's dtype (numpy.array construction dominated ~74% of merge
    # runtime). Instead hold means/vars/present as numpy arrays with a label->row map, combine rows
    # arithmetically in place, and rebuild small DataFrames only for the per-round DE call. Native dtype
    # is preserved so results are bit-identical to the old .loc/.drop path. Reduced-space means
    # (cluster_means_rd, ~32 cols) stay a DataFrame -- cheap, and merged with merge_cluster_means_vars.
    genes = cluster_means.columns
    labels = list(cluster_means.index)
    means_np = np.array(cluster_means.values, copy=True)
    vars_np = np.array(cluster_variances.values, copy=True)
    present_np = np.array(present_cluster_means.values, copy=True)
    row_of = {lab: i for i, lab in enumerate(labels)}   # label -> fixed row index (rows never move)
    live = set(labels)                                  # labels still active (dead rows stay but unused)

    # DE-score cache keyed by frozenset(pair) -> {'score', 'num'} (mirrors R merge_cl_big's `de.genes`
    # cache: pairs are recomputed only when new; entries touching a merged cluster are invalidated).
    de_cache: Dict[frozenset, Dict[str, float]] = {}
    half_th = score_th / 2.0
    merged_cluster_dsts = None
    while len(cluster_assignments.keys()) > 1:
        # 'aligned' re-examines the k nearest clusters of EVERY cluster each round (like R get_knn_pairs);
        # 'fast' only re-examines merged-destination clusters' neighbors after round 1 (previous behavior).
        scope = None if merge_mode == 'aligned' else merged_cluster_dsts
        logger.info(f"Getting {k} nearest clusters")
        neighbor_pairs = order_pairs(get_k_nearest_clusters(cluster_means_rd, scope, k))
        # dedupe to unique unordered pairs
        seen = set(); uniq_pairs = []
        for p in neighbor_pairs:
            key = frozenset(p)
            if key not in seen and len(key) == 2:
                seen.add(key); uniq_pairs.append(p)
        neighbor_pairs = uniq_pairs
        logger.info(f"Completed {k} nearest clusters")
        if len(neighbor_pairs) == 0:
            break

        # Step 4: compute DE only for pairs not already cached (matching R's `new.pairs`)
        to_compute = [p for p in neighbor_pairs if frozenset(p) not in de_cache]
        if to_compute:
            logger.info(f"Calculating de scores for {len(to_compute)} new pairs using {de_method}")
            # per-cluster cell counts for the DE test; if max_cl_size is None the cap is DISABLED
            # (use all cells), else cap at max_cl_size (recomputed each round so merged clusters re-cap)
            cl_size_de = cl_size if max_cl_size is None else {c: min(n, max_cl_size) for c, n in cl_size.items()}
            # rebuild small DataFrames of the LIVE clusters from the numpy state (single-block wrap, cheap;
            # matches the old path which passed the full live-cluster means/vars/present each round)
            live_labels = [lab for lab in labels if lab in live]
            live_rows = [row_of[lab] for lab in live_labels]
            means_df = pd.DataFrame(means_np[live_rows], index=live_labels, columns=genes, copy=False)
            present_df = pd.DataFrame(present_np[live_rows], index=live_labels, columns=genes, copy=False)
            if de_method == 'ebayes':
                vars_df = pd.DataFrame(vars_np[live_rows], index=live_labels, columns=genes, copy=False)
                new_scores = tc.de_pairs_ebayes(
                    to_compute, means_df, vars_df,
                    present_df, cl_size_de, thresholds,
                )
            elif de_method == 'chisq':
                new_scores = tc.de_pairs_chisq(
                    to_compute, means_df, present_df, cl_size_de, thresholds,
                )
            else:
                raise ValueError(f'Unknown de_method {de_method}, must be one of [chisq, ebayes]')
            for pr, row in new_scores.iterrows():
                de_cache[frozenset(pr)] = {'score': float(row.score), 'num': float(row.num)}

        # Candidate pairs = R test_merge: score < score_th OR num < min_genes  (#3: strict <)
        candidates = []
        for p in neighbor_pairs:
            d = de_cache[frozenset(p)]
            if d['score'] < score_th or d['num'] < min_genes:
                candidates.append((p, d['score']))
        if len(candidates) == 0:
            break
        candidates.sort(key=lambda x: x[1])   # ascending by score

        # Merge acceptance:
        #   'aligned': merge the single lowest-score pair unconditionally + extras only < score_th/2 (#1).
        #   'fast':    merge EVERY non-conflicting candidate this round (previous behavior).
        merged_clusters = set()
        merged_cluster_dsts = set()
        logger.info("Merging clusters by DE score")
        for i, (pair, score) in enumerate(candidates):
            dst_label, src_label = pair
            if merge_mode == 'aligned' and i != 0 and not (score < half_th):
                continue
            if dst_label in merged_clusters or src_label in merged_clusters:
                continue

            logger.debug(f"Merging cluster {src_label} into {dst_label} -- de score: {score}")
            # reduced-space means (small): keep the DataFrame path -- reads pre-merge sizes from
            # cluster_assignments, so must run BEFORE the assignment update below.
            merge_cluster_means_vars(cluster_assignments, src_label, dst_label, cluster_means_rd, None)
            # wide means/vars/present: numpy row combine, replicating merge_cluster_means_vars EXACTLY.
            # Subtlety: in the original, `mean2 = cluster_means.loc[dest]` is a pandas VIEW that is
            # overwritten to mean_comb (`cluster_means.loc[dest] = mean_comb`) BEFORE the variance line,
            # so the (mean2 - mean_comb)**2 term is identically zero. Only the source term survives.
            n1 = cl_size[src_label]; n2 = cl_size[dst_label]
            rs = row_of[src_label]; rd = row_of[dst_label]
            m1 = means_np[rs]; m2 = means_np[rd]
            mean_comb = (m1 * n1 + m2 * n2) / (n1 + n2)
            v1 = vars_np[rs]; v2 = vars_np[rd]
            var_comb = 1 / (n1 + n2 - 1) * (
                (n1 - 1) * v1 + n1 * (m1 - mean_comb) ** 2 +
                (n2 - 1) * v2
            )
            p1 = present_np[rs]; p2 = present_np[rd]
            present_comb = (p1 * n1 + p2 * n2) / (n1 + n2)
            means_np[rd] = mean_comb; vars_np[rd] = var_comb; present_np[rd] = present_comb
            # update assignments (mirrors merge_two_clusters) + bookkeeping
            cluster_assignments[dst_label].extend(cluster_assignments[src_label])
            cluster_assignments.pop(src_label)
            live.discard(src_label)
            merged_clusters.add(src_label)
            merged_clusters.add(dst_label)
            merged_cluster_dsts.add(dst_label)
            cl_size[dst_label] += cl_size[src_label]
            cl_size.pop(src_label)

        if not merged_clusters:
            break   # nothing merged this round -> converged
        # invalidate cached scores for any pair touching a merged cluster (means changed)
        for key in [key for key in de_cache if key & merged_clusters]:
            del de_cache[key]

    # Rebuild DataFrames of the surviving (live) clusters from the numpy state and return them.
    # cluster_assignments was mutated in place (unchanged contract); the caller reuses these merged
    # means/vars/present for marker selection, so they must reflect the merges.
    final_labels = [lab for lab in labels if lab in live]
    final_rows = [row_of[lab] for lab in final_labels]
    cluster_means = pd.DataFrame(means_np[final_rows], index=final_labels, columns=genes)
    cluster_variances = pd.DataFrame(vars_np[final_rows], index=final_labels, columns=genes)
    present_cluster_means = pd.DataFrame(present_np[final_rows], index=final_labels, columns=genes)
    return cluster_means, cluster_variances, present_cluster_means


def _resolve_batch_thresholds(
        base_thresholds: Dict[str, Any],
        overrides: Optional[Dict[Any, Dict[str, Any]]],
        batches: List[Any]
) -> Dict[Any, Dict[str, Any]]:
    """
    Per-batch thresholds: the shared thresholds, with optional per-batch overrides layered on top.

    Batches with no entry inherit the shared values unchanged, so adding a batch does not mean
    restating every threshold. Per-batch values matter when modalities differ in depth (a 300-gene
    panel cannot reach the same `score_thresh` as whole-transcriptome data); scrattch.bigcat's
    de.param.list is the same idea.
    """
    overrides = overrides or {}
    unknown = set(overrides) - set(batches)
    if unknown:
        raise ValueError(
            f"batch_aware_merging['thresholds'] has entries for batches not present in the data: "
            f"{sorted(map(str, unknown))}. Present batches: {sorted(map(str, batches))}"
        )
    resolved = {}
    for batch in batches:
        merged = base_thresholds.copy()
        merged.update(overrides.get(batch, {}))
        resolved[batch] = merged
    return resolved


def _r_pair_name(a: Any, b: Any) -> str:
    """
    The pair's name exactly as scrattch.bigcat builds it: the two labels as strings, lexicographic
    min and max, joined by '_' (paste(pmin(P1,P2), pmax(P1,P2)) on character vectors -- so labels
    6 and 22 name the pair "22_6", because "2" < "6" string-wise). Used ONLY to break score ties in
    the same order R's tapply/stable-sort does.
    """
    sa, sb = str(a), str(b)
    return f"{sa}_{sb}" if sa <= sb else f"{sb}_{sa}"


def _dissolve_cl_small(
        cluster_assignments: Dict[Any, List],
        cl_means_reduced: pd.DataFrame,
        adata_reduced: ad.AnnData,
        batch_by_obs: np.ndarray,
        batch_cs_th: Dict[Any, int],
) -> List[Any]:
    """
    Dissolve clusters that are below min.cells in EVERY batch, remapping their cells individually.

    Port of merge_cl_multiple's cl.small step (harmonize_merge.R:400-495, joint.rd.dat branch): such
    clusters can never be judged by any batch, and since unjudgeable pairs never merge they would
    otherwise survive to the output unvetted. Each cell is assigned to the nearest big-cluster
    centroid in the reduced space by COSINE similarity (map_cells_knn's Annoy.Cosine default,
    annotate.R:704), so a dissolved cluster's cells may scatter to different destinations --
    deliberately unlike the pooled path's whole-cluster merge.

    Centroids are the PRE-dissolve ones, as in R (cl.dat is computed before the remap, :427).
    Mutates cluster_assignments in place; returns the dissolved labels.
    """
    batch_by_obs = np.asarray(batch_by_obs)
    batches = list(batch_cs_th.keys())
    cl_small, cl_big = [], []
    for lab, idx in cluster_assignments.items():
        b_of = batch_by_obs[np.asarray(idx)]
        if all((b_of == b).sum() < batch_cs_th[b] for b in batches):
            cl_small.append(lab)
        else:
            cl_big.append(lab)
    if not cl_small:
        return []
    if not cl_big:
        # R returns NULL here (nothing mergeable at all); for a library call, leaving the clustering
        # untouched with a loud warning is the less destructive contract.
        warnings.warn("every cluster is below min.cells in every batch; cl.small dissolve skipped")
        return []

    small_cells = np.concatenate([np.asarray(cluster_assignments[lab]) for lab in cl_small])
    query = adata_reduced.X[small_cells, :]
    if not isinstance(query, np.ndarray):
        query = np.asarray(query.todense()) if hasattr(query, 'todense') else np.asarray(query)
    centroids = cl_means_reduced.loc[cl_big].values

    # cosine: normalize both sides, take the largest dot product
    qn = query / np.maximum(np.linalg.norm(query, axis=1, keepdims=True), 1e-12)
    cn = centroids / np.maximum(np.linalg.norm(centroids, axis=1, keepdims=True), 1e-12)
    dest = np.argmax(qn @ cn.T, axis=1)

    for lab in cl_small:
        cluster_assignments.pop(lab)
    for cell, d in zip(small_cells, dest):
        cluster_assignments[cl_big[d]].append(int(cell))
    return cl_small


def _shortlist_pairs_correlation(cluster_means_rd: pd.DataFrame, k: int) -> List[Tuple[Any, Any, float]]:
    """
    Candidate pairs exactly as merge_cl_multiple's embedding branch builds them.

    R (harmonize_merge.R:514-537): `get_cl_sim` = Pearson CORRELATION between cluster centroids
    (`:12`, dims > 2), then `sim_knn` takes each cluster's top-k *including itself* (the diagonal is
    always 1), ties resolved in column order (`harmonize.R:1082-1084`); self-pairs are then dropped,
    pairs are canonicalized by string min/max and deduped, and ordered by similarity descending.

    Two ways this differs from the pooled path's `get_k_nearest_clusters`, both deliberate there and
    wrong here: the pooled path uses normalized EUCLIDEAN distance (matching `merge_cl_big`'s
    Annoy.Euclidean), and it excludes self from the top-k, so k=4 yields 4 real neighbours where R's
    batch-aware shortlist yields 3.

    Returns [(P1, P2, sim)] sorted by sim descending, P1/P2 canonicalized as R names them.
    """
    # R's cl.rd columns come from table(cl) -- lexicographic label order; ties in sim_knn resolve in
    # that column order, so replicate it.
    labels = sorted(cluster_means_rd.index, key=str)
    X = cluster_means_rd.loc[labels].values.astype(np.float64)
    sim = np.corrcoef(X)                                   # clusters x clusters, diag 1
    k_eff = min(k, len(labels))

    seen, out = set(), []
    for i, lab in enumerate(labels):
        row = sim[i]
        th = np.sort(row)[len(labels) - k_eff]             # k-th largest (rowOrderStats)
        neighbors = [j for j in range(len(labels)) if row[j] >= th][:k_eff]   # column order, head k
        for j in neighbors:
            if j == i:
                continue
            a, b = str(lab), str(labels[j])
            key = (a, b) if a <= b else (b, a)             # pmin/pmax on strings
            if key not in seen:
                seen.add(key)
                out.append((labels[i] if a <= b else labels[j],
                            labels[j] if a <= b else labels[i],
                            float(row[j])))
    out.sort(key=lambda t: -t[2])
    return out


def _half_threshold(batch_score_th: Dict[Any, float]) -> float:
    """
    The bar for merging *extra* pairs in a round, beyond the single lowest-scoring one.

    scrattch.bigcat uses half the **mean** of the per-dataset `de.score.th` (harmonize_merge.R:399),
    not the min -- so one strict batch cannot tighten the rule for all the others. Identical to
    `score_thresh / 2` when every batch shares a threshold; only differs when they don't.
    """
    return (sum(batch_score_th.values()) / len(batch_score_th)) / 2.0


def _ebayes_fit(cluster_variances: pd.DataFrame, cl_size: Dict[Any, int]) -> Dict[str, Any]:
    """
    The eBayes variance-moderation fit that `_de_one_pair` needs, for one batch's statistics.

    Fitting per batch means each batch gets its own variance prior. That is intended: scrattch.bigcat
    likewise computes DE independently within each dataset, so a batch's verdict depends only on its
    own cells.
    """
    sigma_sq, df, stdev_unscaled = get_linear_fit_vals(cluster_variances, cl_size)
    sigma_sq_post, _var_prior, df_prior = moderate_variances(sigma_sq, df)
    return {
        'sqrt_sigma': np.sqrt(sigma_sq_post),
        'stdev_unscaled': stdev_unscaled,
        'df': df,
        'df_prior': df_prior,
        'df_pooled': np.sum(df),
    }


def _conserved_genes(
        detail_by_batch: Dict[Any, pd.DataFrame],
        batch_means: Dict[Any, Dict[Any, np.ndarray]],
        gene_pos: Dict[Any, int],
        dst_label: Any,
        src_label: Any,
        lfc_thresh: float,
        lfc_conservation_th: float,
        gene_sets_by_batch: Optional[Dict[Any, Set[str]]] = None,
) -> Set[Tuple[Any, Any]]:
    """
    Genes whose fold change points the same way in enough batches to count as evidence.

    Mirrors scrattch.bigcat's lfc-conservation filter (`de_genes_pairs_multiple`): take every
    (gene, direction) any batch called differential, then check that gene's fold change in *every*
    batch holding both clusters -- including batches too small to judge the pair, and batches that
    did not call the gene DE. Keep it only if it clears `lfc_thresh` in the stated direction in at
    least `lfc_conservation_th` of them (R: `lfc.num >= lfc.conservation.th * set.num`).

    A gene up in cluster A according to one batch and up in cluster B according to another is
    contradictory -- the signature of a batch effect -- so it is dropped for every batch, not only
    the disagreeing one. This is what makes the merge batch-aware.

    Parameters
    ----------
    detail_by_batch:
        {batch: per-gene DE detail frame with columns gene / P1, where P1 is the cluster the gene is
        up in}, as returned by `_de_one_pair`
    batch_means:
        {batch: {cluster label: gene-mean vector}} for the two clusters under test
    gene_pos:
        gene name -> column index, shared across batches (all batches index the same var axis)

    Returns
    -------
    set of (gene, up_label) pairs that survive
    """
    candidates = set()
    for detail in detail_by_batch.values():
        if detail is not None and len(detail):
            candidates.update(zip(detail['gene'], detail['P1']))

    kept = set()
    for gene, up_label in candidates:
        down_label = src_label if up_label == dst_label else dst_label
        gene_idx = gene_pos.get(gene)
        if gene_idx is None:
            continue
        measured = 0
        agreed = 0
        for batch, means in batch_means.items():
            if up_label not in means or down_label not in means:
                continue    # this batch lacks one of the clusters -- it cannot speak to this gene
            if gene_sets_by_batch is not None and gene not in gene_sets_by_batch[batch]:
                continue    # this batch's reference does not measure the gene (R: gene %in% tmp.genes,
                            # harmonize_merge.R:160) -- a platform cannot vote on a gene it cannot see.
                            # A gene measured by ONE platform therefore has set.num=1 and is trivially
                            # conserved: single-platform evidence counts, as in R.
            measured += 1
            if means[up_label][gene_idx] - means[down_label][gene_idx] > lfc_thresh:
                agreed += 1
        if measured and agreed >= lfc_conservation_th * measured:
            kept.add((gene, up_label))
    return kept


def _combine_rows(means_np, vars_np, present_np, rs, rd, n1, n2, exact=False):
    """
    Fold source row `rs` into destination row `rd`, in place.

    exact=False reproduces the pooled merge loop's arithmetic, including its deliberate asymmetry:
    the destination's `(m2 - mean_comb) ** 2` variance term is absent because in the original pandas
    path `mean2` was a view already overwritten by `mean_comb`, making that term identically zero.

    exact=True includes that term -- the true combined sample variance. The batch-aware path uses
    this, because R updates per-set sqr.means by weighted average (harmonize_merge.R:330-332), which
    is exact; matching the pooled path's asymmetry there would drift from R in every round after a
    merge.

    n1 or n2 of zero means one cluster has no cells in this batch, which the pooled path never sees;
    the surviving side is then carried over untouched.
    """
    if n1 == 0:
        return
    if n2 == 0:
        means_np[rd] = means_np[rs]
        vars_np[rd] = vars_np[rs]
        present_np[rd] = present_np[rs]
        return

    m1 = means_np[rs]
    m2 = means_np[rd]
    mean_comb = (m1 * n1 + m2 * n2) / (n1 + n2)
    if n1 + n2 > 1:
        var_comb = 1 / (n1 + n2 - 1) * (
            (n1 - 1) * vars_np[rs] + n1 * (m1 - mean_comb) ** 2 +
            (n2 - 1) * vars_np[rd] +
            (n2 * (m2 - mean_comb) ** 2 if exact else 0.0)
        )
    else:
        var_comb = np.zeros_like(vars_np[rd])
    present_comb = (present_np[rs] * n1 + present_np[rd] * n2) / (n1 + n2)
    means_np[rd] = mean_comb
    vars_np[rd] = var_comb
    present_np[rd] = present_comb


def merge_clusters_by_de_batch_aware(
    cluster_assignments: Dict[Any, List],
    cluster_means: pd.DataFrame,
    cluster_variances: pd.DataFrame,
    present_cluster_means: pd.DataFrame,
    stats_by_batch: Dict[Any, Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]],
    sizes_by_batch: Dict[Any, Dict[Any, int]],
    cluster_means_rd: pd.DataFrame,
    thresholds: Dict[str, Any],
    batch_thresholds: Dict[Any, Dict[str, Any]],
    lfc_conservation_th: float = 0.7,
    conservation_lfc_th: Optional[float] = None,
    k: Optional[int] = 2,
    pair_batch: Optional[int] = None,
    genes_by_batch: Optional[Dict[Any, List[str]]] = None,
):
    """
    Merge clusters by DE, asking each batch separately and merging only when none objects.

    The pooled merge (`merge_clusters_by_de`) computes one set of statistics over all cells and asks
    once whether two clusters are distinguishable. A difference it finds may be biology or may be one
    platform's artefact, and it cannot tell which. This variant mirrors scrattch.bigcat's
    `merge_cl_multiple`: it computes DE *within* each batch, drops genes whose fold change does not
    point the same way across batches, and merges a pair only when every batch that can judge it
    votes to merge.

    Two asymmetries carried over from R, both deliberate:

    - **A batch with no DE evidence votes to merge.** R's `test_merge` returns TRUE on an empty gene
      set. "I see nothing separating these" is a merge vote, not an abstention -- otherwise a batch
      that measures a pair poorly would silently keep it split.
    - **Objection is a veto, agreement is only permission.** One batch that can separate the pair
      keeps it split regardless of how many others would merge. Real biology visible to a single
      modality is preserved; that is the point of the exercise.

    A pair no batch can judge -- no batch holds at least `cluster_size_thresh` cells of both clusters
    -- is never merged, matching R's pre-filter. This is not a batch-artefact detector: a cluster
    seen in only one batch and distinguishable there is *kept*, because that batch vetoes.

    Clustering itself is unchanged and batch-agnostic; it runs on the integrated embedding, as does
    the nearest-neighbour search that proposes candidate pairs here. Only the merge test is split by
    batch.

    Parameters
    ----------
    cluster_assignments:
        map of cluster label to cell idx; mutated in place, as in the pooled path
    cluster_means / cluster_variances / present_cluster_means:
        pooled statistics, carried alongside the per-batch ones purely so the merged pooled
        statistics can be returned for downstream marker selection
    stats_by_batch:
        {batch: (means, present, variances)} from `get_cluster_means_per_batch`
    sizes_by_batch:
        {batch: {cluster: n_cells}}, zero-filled for clusters absent from a batch
    cluster_means_rd:
        cluster means in the reduced (integrated) space, used for the neighbour search
    thresholds:
        shared thresholds; also carries the loop controls `merge_mode` and `max_cl_size`
    batch_thresholds:
        {batch: thresholds} from `_resolve_batch_thresholds`
    lfc_conservation_th:
        fraction of measuring batches that must agree on a gene's fold-change direction
        (R `lfc.conservation.th`)
    conservation_lfc_th:
        fold change a gene must clear in a batch to count as agreeing. Defaults to
        `thresholds['lfc_thresh']`. R hardcodes `1` here (harmonize_merge.R:167) rather than reading
        `de.param$lfc.th` -- and that `1` is in log2, while this pipeline works in natural log, so it
        is exposed as a parameter instead of a literal.
    k:
        number of cluster neighbors
    genes_by_batch:
        {batch: genes that batch's reference actually measures}. Multi-platform references differ
        (WMB2: 32,285 / 32,285 / 21,899 genes), and R runs each platform's DE on its own gene list.
        With one union axis, a gene absent from a platform is indistinguishable from an unexpressed
        one, which corrupts two things: ~10^4 fake zero-variance genes distort that platform's eBayes
        prior, and the conservation denominator counts platforms that never measured the gene. Each
        batch's DE is therefore restricted to its own genes, and a gene measured by only one platform
        is trivially conserved (set.num=1) -- single-platform evidence counts, as in R. None (default)
        means every batch measures every gene.
    pair_batch:
        None (default): compute DE for every shortlisted pair each round, then merge the weakest --
        the deterministic reference behaviour. An integer N mirrors R's `pairBatch` (default 100
        there): new pairs are DE-tested N at a time in order of decreasing centroid similarity, and
        as soon as ANY cached pair is mergeable the round stops testing and merges. Cheaper, but
        which pair merges first then depends on how many pairs the round happened to compute, so the
        result is a function of scan order, not just of the input. Exposed to measure that
        cost/benefit; production R runs used 100.

    Returns
    -------
    merged pooled (cluster_means, cluster_variances, present_cluster_means)
    """
    thresholds = thresholds.copy()
    thresholds.pop('score_thresh', None)
    thresholds.pop('min_genes', None)
    thresholds.pop('low_thresh', None)
    max_cl_size = thresholds.pop('max_cl_size', None)
    merge_mode = thresholds.pop('merge_mode', 'aligned')

    batches = list(stats_by_batch.keys())
    # per-batch DE settings, unpacked once
    batch_filt = {b: {key: batch_thresholds[b].get(key) for key in _FILTER_KEYS} for b in batches}
    batch_score_th = {b: batch_thresholds[b]['score_thresh'] for b in batches}
    batch_min_genes = {b: batch_thresholds[b]['min_genes'] for b in batches}
    batch_cs_th = {b: batch_thresholds[b]['cluster_size_thresh'] for b in batches}
    batch_padj_alpha = {b: batch_thresholds[b]['padj_thresh'] for b in batches}
    if conservation_lfc_th is None:
        conservation_lfc_th = thresholds.get('lfc_thresh') or 0.0

    genes = cluster_means.columns
    gene_pos = {g: i for i, g in enumerate(genes)}
    # per-batch gene axes: column indices into the union axis, the gene names in union order, and a
    # set for membership tests. None -> every batch measures every gene.
    if genes_by_batch is not None:
        gene_col_idx, genes_by_batch_ordered, gene_sets_by_batch = {}, {}, {}
        for b, gs in genes_by_batch.items():
            gset = set(gs)
            unknown = gset - set(genes)
            if unknown:
                raise ValueError(
                    f"genes_by_batch[{b!r}] has {len(unknown)} genes not in the data, "
                    f"e.g. {sorted(unknown)[:3]}"
                )
            idx = np.flatnonzero(np.fromiter((g in gset for g in genes), dtype=bool, count=len(genes)))
            gene_col_idx[b] = idx
            genes_by_batch_ordered[b] = genes[idx]
            gene_sets_by_batch[b] = gset
    else:
        gene_col_idx = None
        genes_by_batch_ordered = None
        gene_sets_by_batch = None
    labels = list(cluster_means.index)
    row_of = {lab: i for i, lab in enumerate(labels)}    # shared by pooled and every batch
    live = set(labels)

    # Pooled state, updated by exactly the same arithmetic as the pooled path so the statistics
    # returned for marker selection match what those merges would have produced there.
    pooled = (
        np.array(cluster_means.values, copy=True),
        np.array(cluster_variances.values, copy=True),
        np.array(present_cluster_means.values, copy=True),
    )
    cl_size = {lab: len(cluster_assignments[lab]) for lab in labels}

    # Per-batch state, laid out over the FULL label set so `row_of` is shared. Clusters absent from a
    # batch keep an all-zero row that nothing reads: they are excluded from that batch's eBayes fit
    # and from every pair it judges, via sizes_by_batch.
    batch_state = {}
    cl_size_by_batch = {}
    for b in batches:
        b_means, b_present, b_vars = stats_by_batch[b]
        m = np.zeros((len(labels), len(genes)), dtype=b_means.values.dtype)
        v = np.zeros_like(m)
        p = np.zeros_like(m)
        for lab in b_means.index:
            m[row_of[lab]] = b_means.loc[lab].values
            v[row_of[lab]] = b_vars.loc[lab].values
            p[row_of[lab]] = b_present.loc[lab].values
        batch_state[b] = (m, v, p)
        cl_size_by_batch[b] = {lab: sizes_by_batch[b].get(lab, 0) for lab in labels}

    de_cache: Dict[frozenset, Dict[str, float]] = {}
    half_th = _half_threshold(batch_score_th)
    while len(cluster_assignments.keys()) > 1:
        # Shortlist exactly as merge_cl_multiple does (correlation, self-in-top-k, sim-ordered) --
        # NOT the pooled path's Euclidean get_k_nearest_clusters; see _shortlist_pairs_correlation.
        # R re-shortlists every cluster every round (there is no 'fast' scope in merge_cl_multiple).
        logger.info(f"Getting {k} nearest clusters (correlation)")
        shortlist = _shortlist_pairs_correlation(cluster_means_rd, k)
        pair_sim = {(a, b): s for a, b, s in shortlist}
        neighbor_pairs = [(a, b) for a, b, _s in shortlist]
        logger.info(f"Completed {k} nearest clusters")
        if len(neighbor_pairs) == 0:
            break

        # R drops unjudgeable pairs BEFORE batching (merge.pairs.select, harmonize_merge.R:528-533),
        # not at verdict time -- and the placement matters beyond bookkeeping: batch membership
        # determines the variance-fit population, so an unjudgeable pair left in a chunk changes
        # every other pair's score in that chunk. (Observed: an unjudgeable rank-5 pair in Python's
        # chunk 1 vs R pulling in the rank-101 pair instead.)
        neighbor_pairs = [p for p in neighbor_pairs
                          if any(cl_size_by_batch[b][p[0]] >= max(batch_cs_th[b], 1)
                                 and cl_size_by_batch[b][p[1]] >= max(batch_cs_th[b], 1)
                                 for b in batches)]
        to_compute = [p for p in neighbor_pairs if frozenset(p) not in de_cache]
        if to_compute:
            logger.info(f"Calculating batch-aware de scores for {len(to_compute)} new pairs "
                        f"across {len(batches)} batches")

            def _build_fits(chunk_pairs):
                # THE FIT POPULATION IS THE CHUNK'S OWN CLUSTERS, not all live clusters.
                # R's de_selected_pairs restricts cl to the clusters appearing in the pairs it was
                # handed (`select.cl <- unique(c(pairs$P1, pairs$P2))`, de.genes.R:1152-1154) before
                # fitting sigma/squeezeVar -- so a pair's variance model, and hence its SCORE,
                # depends on which other pairs shared its evaluation call. Early rounds (many pairs,
                # most clusters) approximate an all-cluster fit; late rounds fit over a handful of
                # clusters and score very differently. An earlier version fitted over all live
                # clusters, which reproduced round-1 scores exactly and diverged in later rounds.
                chunk_cls = set()
                for a_, b_ in chunk_pairs:
                    chunk_cls.add(a_); chunk_cls.add(b_)
                fits_, frames_ = {}, {}
                for b in batches:
                    fit_cls = [lab for lab in labels
                               if lab in chunk_cls
                               and cl_size_by_batch[b][lab] >= max(batch_cs_th[b], 1)]
                    if len(fit_cls) < 2:
                        continue
                    m, v, p_arr = batch_state[b]
                    rows = [row_of[lab] for lab in fit_cls]
                    sizes = {lab: cl_size_by_batch[b][lab] for lab in fit_cls}
                    # restrict to the genes THIS batch measures: genes absent from the batch's
                    # reference would enter the across-genes prior as fake zero-variance rows
                    if gene_col_idx is not None and b in gene_col_idx:
                        cols = gene_col_idx[b]
                        b_genes = genes_by_batch_ordered[b]
                        m_b, v_b, p_b = m[np.ix_(rows, cols)], v[np.ix_(rows, cols)], p_arr[np.ix_(rows, cols)]
                    else:
                        b_genes = genes
                        m_b, v_b, p_b = m[rows], v[rows], p_arr[rows]
                    vars_df = pd.DataFrame(v_b, index=fit_cls, columns=b_genes, copy=False)
                    fits_[b] = _ebayes_fit(vars_df, sizes)
                    frames_[b] = (
                        pd.DataFrame(m_b, index=fit_cls, columns=b_genes, copy=False),
                        pd.DataFrame(p_b, index=fit_cls, columns=b_genes, copy=False),
                        sizes if max_cl_size is None
                        else {lab: min(n, max_cl_size) for lab, n in sizes.items()},
                    )
                return fits_, frames_

            if pair_batch is not None:
                # R's pairBatch semantics (harmonize_merge.R:527-567): new pairs are tested in
                # similarity-ordered chunks, and testing stops for the round as soon as anything in
                # the cache is mergeable. Pairs left uncomputed carry no cache entry, so they are
                # re-shortlisted and tested in a later round, exactly as in R. The ordering uses the
                # shortlist's own correlation similarities.
                to_compute.sort(key=lambda p_: (-pair_sim[p_], _r_pair_name(*p_)))
                chunks = [to_compute[i:i + pair_batch] for i in range(0, len(to_compute), pair_batch)]
            else:
                chunks = [to_compute]

            b_axes = ({b: genes_by_batch_ordered.get(b, genes) for b in batches}
                      if genes_by_batch_ordered is not None else {b: genes for b in batches})
            n_done = 0
            for chunk in chunks:
                fits, frames = _build_fits(chunk)   # per-chunk fit population, as in R
                for pair in chunk:
                    verdict = _judge_pair_across_batches(
                        pair, batches, fits, frames, batch_state, cl_size_by_batch, row_of, b_axes,
                        gene_pos, batch_filt, batch_padj_alpha, batch_cs_th,
                        batch_score_th, batch_min_genes, lfc_conservation_th, conservation_lfc_th,
                        gene_sets_by_batch,
                    )
                    verdict['pair'] = pair
                    de_cache[frozenset(pair)] = verdict
                n_done += len(chunk)
                # early stop only in pair_batch mode; R checks the WHOLE cache, not just this chunk
                if pair_batch is not None and any(v['to_merge'] for v in de_cache.values()):
                    if n_done < len(to_compute):
                        logger.info(f"pair_batch={pair_batch}: mergeable pair found; "
                                    f"{len(to_compute) - n_done} shortlisted pairs left untested this round")
                    break

        # A pair is a candidate only if every judging batch voted to merge (`to_merge`); pairs no
        # batch could judge carry to_merge=False and so never merge.
        #
        # Candidates come from the WHOLE surviving cache, not just this round's shortlist: R's
        # test_merge_multiple evaluates every pair in de.genes.list (harmonize_merge.R:381), and
        # entries only leave the cache when a merge touches one of their clusters. A pair shortlisted
        # two rounds ago whose clusters are untouched is still a valid candidate -- its statistics
        # have not changed -- even if centroid drift pushed it out of the current k-NN.
        candidates = [(v['pair'], v['score']) for v in de_cache.values() if v['to_merge']]
        if len(candidates) == 0:
            break
        # Ties in score are COMMON at scale (a pair with no DE evidence anywhere scores exactly 0,
        # and the per-gene cap makes round sums like 20 recur), and a tie decides which pair gets the
        # round's unconditional merge. Breaking ties by list order would make that depend on set/hash
        # iteration order -- non-deterministic across processes for string labels. Break them exactly
        # as R does instead: R's tapply orders pairs by their name -- the two labels as strings,
        # lexicographic min and max joined by '_' (so labels 6 and 22 name the pair "22_6") -- and
        # its stable sort preserves that order among equal scores (harmonize_merge.R:394-395).
        # Adopting the same rule makes tied merges identical across the two implementations, not
        # merely deterministic within each.
        candidates.sort(key=lambda x: (x[1], _r_pair_name(*x[0])))

        merged_clusters = set()
        logger.info("Merging clusters by batch-aware DE score")
        for i, (pair, score) in enumerate(candidates):
            if merge_mode == 'aligned' and i != 0 and not (score < half_th):
                continue
            # The SURVIVING label is the larger cluster's, as in R's merge_x_y (harmonize_merge.R:
            # 317-320, pooled sizes; ties keep the table's sorted-name order, i.e. the
            # lexicographically smaller label). Label choice does not change any merged cluster's
            # membership, but it names later pairs -- and pair names break score ties -- so matching
            # it keeps tied merge ORDER identical to R in rounds after this one.
            a_label, b_label = pair
            if a_label in merged_clusters or b_label in merged_clusters:
                continue
            if cl_size[a_label] < cl_size[b_label] or \
               (cl_size[a_label] == cl_size[b_label] and str(b_label) < str(a_label)):
                dst_label, src_label = b_label, a_label
            else:
                dst_label, src_label = a_label, b_label

            # INFO, matching R's verbose=TRUE which prints every merge -- the merge sequence is
            # the primary artefact for diagnosing R/Python divergence at scale.
            logger.info(f"BATCH-AWARE MERGE {dst_label} <- {src_label} score={score} "
                         f"num={de_cache.get(frozenset(pair), {}).get('num')}")
            merge_cluster_means_vars(cluster_assignments, src_label, dst_label, cluster_means_rd, None)
            rs, rd = row_of[src_label], row_of[dst_label]
            # exact=True: R updates per-set sqr.means by weighted average (harmonize_merge.R:330-332),
            # which is an EXACT combined variance -- unlike the pooled path's deliberate asymmetry.
            _combine_rows(*pooled, rs, rd, cl_size[src_label], cl_size[dst_label], exact=True)
            for b in batches:
                _combine_rows(*batch_state[b], rs, rd,
                              cl_size_by_batch[b][src_label], cl_size_by_batch[b][dst_label],
                              exact=True)
                cl_size_by_batch[b][dst_label] += cl_size_by_batch[b][src_label]
                cl_size_by_batch[b][src_label] = 0

            cluster_assignments[dst_label].extend(cluster_assignments[src_label])
            cluster_assignments.pop(src_label)
            live.discard(src_label)
            merged_clusters.add(src_label)
            merged_clusters.add(dst_label)
            cl_size[dst_label] += cl_size[src_label]
            cl_size.pop(src_label)

        if not merged_clusters:
            break
        for key in [key for key in de_cache if key & merged_clusters]:
            del de_cache[key]

    final_labels = [lab for lab in labels if lab in live]
    final_rows = [row_of[lab] for lab in final_labels]
    means_np, vars_np, present_np = pooled
    cluster_means = pd.DataFrame(means_np[final_rows], index=final_labels, columns=genes)
    cluster_variances = pd.DataFrame(vars_np[final_rows], index=final_labels, columns=genes)
    present_cluster_means = pd.DataFrame(present_np[final_rows], index=final_labels, columns=genes)
    return cluster_means, cluster_variances, present_cluster_means


def _judge_pair_across_batches(
        pair, batches, fits, frames, batch_state, cl_size_by_batch, row_of, genes_by_batch, gene_pos,
        batch_filt, batch_padj_alpha, batch_cs_th, batch_score_th, batch_min_genes,
        lfc_conservation_th, conservation_lfc_th, gene_sets_by_batch,
) -> Dict[str, Any]:
    """
    Ask every batch that can judge `pair` whether to merge it, and combine the answers.

    Returns {'score', 'num', 'to_merge'}. `score` is the largest score any judging batch gave -- the
    strongest evidence against merging -- so ordering candidates by it merges the least separable
    pairs first, matching R's `tapply(sc, pair, max)`. Pairs no batch can judge come back
    to_merge=False with an infinite score.
    """
    dst_label, src_label = pair
    # R's pre-filter: pmin(size[P1, set], size[P2, set]) >= de.param[[set]]$min.cells, for any set
    # (harmonize_merge.R:528-533). The `max(..., 1)` guards the degenerate cs_th <= 0 case, where a
    # cluster with no cells in this batch would otherwise qualify to judge but be absent from its
    # statistics frame.
    judging = [b for b in batches
               if b in fits
               and cl_size_by_batch[b][dst_label] >= max(batch_cs_th[b], 1)
               and cl_size_by_batch[b][src_label] >= max(batch_cs_th[b], 1)]
    if not judging:
        return {'score': float('inf'), 'num': 0.0, 'to_merge': False}

    detail_by_batch = {}
    for b in judging:
        means_df, present_df, sizes_de = frames[b]
        _row, detail = _de_one_pair(
            dst_label, src_label, means_df, present_df, sizes_de, genes_by_batch[b],
            filt=batch_filt[b], padj_alpha=batch_padj_alpha[b],
            top_n=CONSERVATION_TOP_N, want_detail=True, **fits[b],
        )
        detail_by_batch[b] = detail

    # Cross-batch gene filter: a gene counts as evidence only where its fold change points the same
    # way in enough batches. Genes that flip direction between batches -- the signature of a batch
    # effect -- are dropped everywhere, which is what keeps a platform artefact from blocking a merge.
    # The denominator is EVERY batch holding both clusters, not just the judging ones: a batch too
    # small to run a DE test can still say which way a gene points, and R counts it (`set.num`).
    batch_means = {}
    for b in batches:
        if cl_size_by_batch[b][dst_label] > 0 and cl_size_by_batch[b][src_label] > 0:
            m = batch_state[b][0]
            batch_means[b] = {dst_label: m[row_of[dst_label]], src_label: m[row_of[src_label]]}
    kept = _conserved_genes(detail_by_batch, batch_means, gene_pos,
                            dst_label, src_label, conservation_lfc_th, lfc_conservation_th,
                            gene_sets_by_batch)

    # Re-score each batch over the surviving genes, then require unanimity.
    max_score = 0.0
    max_num = 0.0
    to_merge = True
    for b in judging:
        detail = detail_by_batch[b]
        if detail is None or len(detail) == 0:
            score, num = 0.0, 0
        else:
            survives = np.fromiter(
                ((g, p1) in kept for g, p1 in zip(detail['gene'], detail['P1'])),
                dtype=bool, count=len(detail),
            )
            score = float(np.minimum(detail['logPval'].values[survives], SCORE_CAP).sum())
            num = int(survives.sum())
        max_score = max(max_score, score)
        max_num = max(max_num, num)
        # R test_merge: too weak a score OR too few genes means this batch cannot tell them apart,
        # so it votes merge. A batch with no evidence at all lands here too (score 0, num 0).
        if not (score < batch_score_th[b] or num < batch_min_genes[b]):
            to_merge = False
            break

    return {'score': max_score, 'num': float(max_num), 'to_merge': to_merge}


def get_k_nearest_clusters(
        cluster_means: pd.DataFrame,
        cluster_labels: Optional[Set[Any]] = None,
        k: Optional[int] = 2
) -> List[Tuple[int, int]]:
    """
    Get k nearest neighbors for each cluster

    Parameters
    ----------
    cluster_means:
        dataframe of cluster means with cluster labels as index
    cluster_labels:
        clusters to calculate nearest neighbors for. If none will return all
    k:
        number of nearest neighbors

    Returns
    -------
    nearest_neighbors:
        list of cluster pairs
    """

    all_cluster_labels = list(cluster_means.index)
    if cluster_labels is None:
        cluster_labels = all_cluster_labels
    else:
        cluster_labels = list(cluster_labels)

    if k >= len(all_cluster_labels):
        logger.debug("k cannot be greater than or the same as the number of clusters. "
                          "Defaulting to number of clusters - 1.")
        k = len(all_cluster_labels) - 1

    similarity = calculate_similarity(
            cluster_means,
            group_rows=all_cluster_labels,
            group_cols=cluster_labels)

    similarity = similarity.unstack().dropna()

    # Get k nearest neighbors
    nearest_neighbors = set()
    for c in cluster_labels:
        # Sort similarities for a cluster
        sorted_similarities = similarity.loc[(c, )].sort_values(ascending=False)

        for i in range(k):
            neighbor_cl = sorted_similarities.index[i]

            # Make sure neighbor doesn't already exist
            if not (neighbor_cl, c) in nearest_neighbors:
                nearest_neighbors.add((c, neighbor_cl))

    return list(nearest_neighbors)


def order_pairs(
        neighbor_pairs: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """
    Order each label such that smaller label is follower by a larger
    e.g: (3,8), (6,7)

    Parameters
    ----------
    neighbor_pairs:
        list of neighbor pairs

    Returns
    -------
    ordered list of neighbor pairs
    """
    ordered_pairs = []

    for p in neighbor_pairs:
        a, b = p
        if b > a:
            ordered_pairs.append((a, b))
        else:
            ordered_pairs.append((b, a))

    return ordered_pairs


def get_cluster_assignments(
        adata: ad.AnnData,
        cluster_label_obs: str = "pheno_louvain"
) -> Dict[Any, List]:
    """

    Parameters
    ----------
    adata:
        AnnData object with with obs including cluster label
    cluster_label_obs:
        cluster label annotations in adata.obs

    Returns
    -------
    cluster_assignments:
        map of cluster label to cell idx
    """
    if cluster_label_obs not in list(adata.obs):
        raise ValueError(f"column {cluster_label_obs} is missing from obs")

    cluster_assignments = defaultdict(list)
    for i, label in enumerate(adata.obs[cluster_label_obs]):
        cluster_assignments[label].append(i)

    return cluster_assignments

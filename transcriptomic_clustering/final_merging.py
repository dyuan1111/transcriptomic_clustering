from typing import Optional, Dict, Set, List, Any, Tuple, Union
from collections import defaultdict
from dataclasses import dataclass, field

import pandas as pd
from scipy import sparse
import numpy as np
from numpy.random import default_rng
import anndata as ad

import transcriptomic_clustering as tc

import logging
import time
logger = logging.getLogger(__name__)


@dataclass
class FinalMergeKwargs:
    """Dataclass for kwargs in final_merge"""
    pca_kwargs: Dict = field(default_factory = lambda: ({}))
    filter_pcs_kwargs: Dict = field(default_factory= lambda: ({}))
    filter_known_modes_kwargs: Dict = field(default_factory = lambda: ({}))
    project_kwargs: Dict = field(default_factory = lambda: ({}))
    merge_clusters_kwargs: Dict = field(default_factory = lambda: ({}))
    latent_kwargs: Dict = field(default_factory = lambda: ({}))


def sample_clusters(
        adata: ad.AnnData,
        cluster_dict: Dict[Any, List],
        n_samples_per_clust: int = 20,
        random_seed: int = 123
) -> pd.Series:
    """
    Create a cell mask containing up to n random cells per cluster

    Parameters
    ----------
    adata
        AnnData object
    cluster_dict
        Dictionary of cluster assignments. values are lists of cell names 

    Returns
    -------
    pd.Series that is True for cells selected for sampling
    """

    rng = default_rng(random_seed)
    cell_samples_ids = []
    for k, v in cluster_dict.items():
        if len(v) > n_samples_per_clust:
            choices = rng.choice(v, size=(n_samples_per_clust,))
        else:
            choices = v
        cell_samples_ids.extend(choices)
    
    cell_samples = adata.obs_names[cell_samples_ids]

    cell_mask = pd.Series(
        index=adata.obs.index,
        dtype=bool,
    )
    cell_mask[cell_samples] = True
    return cell_mask


def _cluster_obs_dict_to_list(obs_by_cluster: Dict[int, List[int]]) -> List[int]:
    """
    Convert a dictionary of cluster assignments to a list of cluster assignments
    """
    # Find the total number of observations
    max_index = max(max(indices) for indices in obs_by_cluster.values())

    # Initialize the list of clusters with None (or some other default value)
    cluster_by_obs = [None] * (max_index + 1)

    # Fill in the list with the corresponding cluster for each observation
    for cluster, cell_ids in obs_by_cluster.items():
        for obs in cell_ids:
            cluster_by_obs[obs] = cluster

    return cluster_by_obs


def final_merge(
        adata: ad.AnnData,
        cluster_assignments: List,
        marker_genes: Optional[set] = None,
        n_samples_per_clust: int = 20,
        random_seed: Optional[int]=None,
        n_jobs: Optional[int] = 1,
        return_markers_df: Optional[bool] = False,
        final_merge_kwargs: FinalMergeKwargs = FinalMergeKwargs(),
        space: Optional[str] = None,
        rm_genes: Optional[List] = None,
) -> Tuple[List[List[int]], Union[pd.DataFrame, set]]:
    """
    Runs a final merging step on cluster assignment results.

    Step 1 builds the reduced space used for candidate-pair selection (`space`):
      * 'markers' - raw normalized marker-gene expression (cells x selected markers), no PCA.
                    Matches scrattch.bigcat `merge.R` (`rd.dat = t(norm.dat[markers, ])`).
      * 'latent'  - a precomputed embedding, via `final_merge_kwargs.latent_kwargs['latent_component']`.
      * 'pca'     - PCA on a per-cluster cell sample restricted to the marker genes, then project.
    Step 2 runs the differential-expression merge (`merge_clusters`) on that reduced space.

    Parameters
    ----------
    adata
        AnnData object (normalized).
    cluster_assignments
        List of arrays of cell ids, one array per cluster (result of onestep_clust/iter_clust).
    marker_genes
        Genes used by the 'markers' and 'pca' spaces.
    space
        Reduced space for candidate selection: 'markers', 'latent', or 'pca'. Default None keeps the
        backward-compatible behavior: 'latent' if a latent_component is set, else 'pca'.
    rm_genes
        Genes to drop from marker_genes for the 'markers' space (e.g. sex/mito), matching R `rm.genes`.

    Batch-aware merging
    -------------------
    `final_merge_kwargs.merge_clusters_kwargs` is forwarded verbatim to `merge_clusters`, so the
    final merge is made batch-aware by adding `batch_aware_merging` to it -- the same key, and the
    same meaning, as in the per-level `merge_clusters_kwargs` used by onestep/iter_clust:

        merge_clusters_kwargs = {
            'thresholds': {...},
            'de_method': 'ebayes',                  # required in batch-aware mode
            'batch_aware_merging': {
                'batch_obs': 'platform',            # column of adata.obs naming each cell's batch
                'lfc_conservation_th': 0.7,
                'thresholds': {'10X_nuclei_v3': {'q1_thresh': 0.3}},  # optional per-batch overrides
            },
        }

    Use it when cells come from several platforms or modalities and a pooled difference between two
    clusters could be a platform artefact rather than biology. Note this makes only the *merge*
    batch-aware; the reduced space that shortlists candidate pairs is unchanged, as in
    scrattch.bigcat where the embedding never votes.

    Requires an in-memory `adata` (per-batch statistics on a backed AnnData would take one full file
    pass per batch).
    """

    obs_by_cluster = defaultdict(lambda: [])
    for i, cell_ids in enumerate(cluster_assignments):
        obs_by_cluster[i] = cell_ids
    
    cluster_by_obs = _cluster_obs_dict_to_list(obs_by_cluster)

    # Validate batch-aware options BEFORE building the reduced space. merge_clusters checks them too,
    # but only once it is called -- by which point space='pca' has already paid for a full PCA and
    # projection. A typo in batch_obs should not cost that.
    _bam = final_merge_kwargs.merge_clusters_kwargs.get('batch_aware_merging')
    if _bam is not None:
        _batch_obs = _bam.get('batch_obs')
        if _batch_obs is None:
            raise ValueError(
                "batch_aware_merging requires 'batch_obs', naming a column of adata.obs"
            )
        if _batch_obs not in adata.obs:
            raise ValueError(
                f"batch_aware_merging['batch_obs']={_batch_obs!r} is not a column of adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )
        if adata.isbacked:
            raise NotImplementedError(
                "Batch-aware final merge requires an in-memory AnnData; the backed path would take "
                "one full file pass per batch. Load the data into memory first."
            )
        _counts = adata.obs[_batch_obs].value_counts().to_dict()
        logger.info(f'Final merge is BATCH-AWARE on obs[{_batch_obs!r}]: {_counts}')
        if len(_counts) < 2:
            logger.warning(
                f'Only one batch present in obs[{_batch_obs!r}]; batch-aware merging degenerates to '
                f'the pooled merge (nothing to cross-check a difference against).'
            )

    # Choose the reduced space that feeds candidate-pair selection in merge_clusters (see docstring).
    if space is None:
        space = 'latent' if final_merge_kwargs.latent_kwargs.get("latent_component") is not None else 'pca'

    if space == 'markers':
        # scrattch.bigcat merge.R: rd.dat = raw normalized marker-gene expression; candidate KNN on markers.
        if marker_genes is None:
            raise ValueError("space='markers' requires marker_genes")
        _rm = set(rm_genes) if rm_genes else set()
        markers_list = [g for g in marker_genes if g not in _rm and g in adata.var_names]
        if not markers_list:
            raise ValueError("No marker genes present in adata.var_names for the marker-gene final merge")
        logger.info(f'Final merge on marker-gene space ({len(markers_list)} genes)')
        projected_adata = adata[:, markers_list].copy()

    elif space == 'pca':

        if marker_genes is None:
            raise ValueError("Need marker genes to run PCA")

        cell_mask = sample_clusters(
            adata=adata,
            cluster_dict=obs_by_cluster,
            n_samples_per_clust=n_samples_per_clust,
            random_seed=random_seed
        )
        gene_mask = pd.Series(
            index=adata.var.index,
            dtype=bool
        )
        gene_mask[marker_genes] = True

        # Do PCA on cell samples and marker genes
        logger.info('Computing PCA on cell samples and marker genes')
        tic = time.perf_counter()
        (components, explained_variance_ratio, explained_variance, means) =  tc.pca(
            adata,
            gene_mask=gene_mask,
            cell_select=cell_mask,
            random_state=random_seed,
            **final_merge_kwargs.pca_kwargs
        )
        logger.info(f'Computed {components.shape[1]} principal components')
        toc = time.perf_counter()
        logger.info(f'PCA Elapsed Time: {toc - tic}')

        # Filter PCA
        logger.info('Filtering PCA Components')
        tic = time.perf_counter()
        components = tc.dimension_reduction.filter_components(
            components,
            explained_variance,
            explained_variance_ratio,
            **final_merge_kwargs.filter_pcs_kwargs
        )
        logger.info(f'Filtered to {components.shape[1]} principal components')
        toc = time.perf_counter()
        logger.info(f'Filter PCA Elapsed Time: {toc - tic}')

        # Project
        logger.info("Projecting into PCA space")
        tic = time.perf_counter()
        projected_adata = tc.project(
            adata, components, means,
            **final_merge_kwargs.project_kwargs
        )
        logger.info(f'Projected Adata Dimensions: {projected_adata.shape}')
        toc = time.perf_counter()
        logger.info(f'Projection Elapsed Time: {toc - tic}')

        # Filter Projection
        #Filter Known Modes
        if final_merge_kwargs.filter_known_modes_kwargs:
            logger.info('Filtering Known Modes')
            tic = time.perf_counter()

            projected_adata = tc.filter_known_modes(
                projected_adata,
                **final_merge_kwargs.filter_known_modes_kwargs
            )

            logger.info(f'Projected Adata Dimensions after Filtering Known Modes: {projected_adata.shape}')
            toc = time.perf_counter()
            logger.info(f'Filter Known Modes Elapsed Time: {toc - tic}')
        else:
            logger.info('No known modes, skipping Filter Known Modes')

    elif space == 'latent':
        logger.info('Extracting latent dims')
        tic = time.perf_counter()

        ## Extract latent dimensions
        projected_adata = tc.latent_project(adata, **final_merge_kwargs.latent_kwargs)

        toc = time.perf_counter()
        logger.info(f'Extracting latent dims Elapsed Time: {toc - tic}')

    else:
        raise ValueError(f"Unknown space {space!r}; use 'markers', 'latent', or 'pca'")

    # Merging
    logger.info('Starting Cluster Merging')
    tic = time.perf_counter()
    cluster_assignments_after_merging, markers = tc.merge_clusters(
        adata_norm=adata,
        adata_reduced=projected_adata,
        cluster_assignments=obs_by_cluster,
        cluster_by_obs=cluster_by_obs,
        return_markers_df=return_markers_df,
        n_jobs=n_jobs,
        **final_merge_kwargs.merge_clusters_kwargs
    )
    logger.info(f'Completed Merging')
    toc = time.perf_counter()
    logger.info(f'Merging Elapsed Time: {toc - tic}')

    cluster_assignments_after_merging = list(cluster_assignments_after_merging.values())

    return cluster_assignments_after_merging, markers

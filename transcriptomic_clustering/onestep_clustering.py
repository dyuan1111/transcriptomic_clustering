from typing import Dict, Optional, List, Any
import logging
from dataclasses import dataclass, field
import time

import math
import numpy as np
import scanpy as sc
import anndata as ad
import transcriptomic_clustering as tc


logger = logging.getLogger(__name__)


@dataclass
class OnestepKwargs:
    """Dataclass for kwargs in onestep_clust"""
    means_vars_kwargs: Dict = field(default_factory = lambda: ({}))
    highly_variable_kwargs: Dict = field(default_factory = lambda: ({}))
    pca_kwargs: Dict = field(default_factory = lambda: ({}))
    filter_pcs_kwargs: Dict = field(default_factory= lambda: ({}))
    filter_known_modes_kwargs: Dict = field(default_factory = lambda: ({}))
    project_kwargs: Dict = field(default_factory = lambda: ({}))
    latent_kwargs: Dict = field(default_factory = lambda: ({}))
    cluster_louvain_kwargs: Dict = field(default_factory = lambda: ({}))
    merge_clusters_kwargs: Dict = field(default_factory = lambda: ({}))


def onestep_clust(
        norm_adata: sc.AnnData,
        onestep_kwargs: OnestepKwargs=OnestepKwargs(),
        random_seed: Optional[int]=None) -> List[np.ndarray]:
    """
    Performs an entire clustering step
    * get mean and variance of each gene
    * determine highly variable genes
    * do pca on a sample of cells
    * filter known pca components
    * cluster cells using louvain clustering
    * merge clusters

    Parameters
    ----------
    norm_adata: log normalized adata (see tc.normalization for computation details)
    onestep_kwargs: Dataclass containg keyword arguements for each function (see OnestepKwargs)
    random_seed: random_seed for functions that use a random seed/random state

    Returns
    -------
    List of arrays of cell ids, one array per cluster
    Set of marker genes (if calculated)

    """
    logger.info('Starting onestep clustering')

    ## Choose either PCA or pre-defined latent space (scVI, scAlign, Seurat, etc.)
    if onestep_kwargs.latent_kwargs.get("latent_component") is None:

        # Means and Variances
        logger.info('Computing Means and Variances of genes')
        tic = time.perf_counter()

        means, variances, gene_mask = tc.get_means_vars_genes(
            adata=norm_adata,
            **onestep_kwargs.means_vars_kwargs
        )

        toc = time.perf_counter()
        logger.info(f'Means Vars Elapsed Time: {toc - tic}')

        #Highly Variable
        logger.info('Computing Highly Variable Genes')
        tic = time.perf_counter()

        highly_variable_mask = tc.highly_variable_genes(
            adata=norm_adata,
            means=means,
            variances=variances,
            gene_mask=gene_mask,
            **onestep_kwargs.highly_variable_kwargs
        )

        toc = time.perf_counter()
        logger.info(f'Highly Variable Elapsed Time: {toc - tic}')

        #PCA
        logger.info('Computing PCA')
        tic = time.perf_counter()

        (components, explained_variance_ratio, explained_variance, means) =  tc.pca(
            norm_adata,
            gene_mask=highly_variable_mask,
            random_state=random_seed,
            **onestep_kwargs.pca_kwargs
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
            **onestep_kwargs.filter_pcs_kwargs
        )

        logger.info(f'Filtered to {components.shape[1]} principal components')
        toc = time.perf_counter()
        logger.info(f'Filter PCA Elapsed Time: {toc - tic}')
        
        #Projection
        logger.info("Projecting normalized adata into PCA space")
        tic = time.perf_counter()

        projected_adata = tc.project(
            norm_adata,
            components, means,
            **onestep_kwargs.project_kwargs
        )
        logger.info(f'Projected Adata Dimensions: {projected_adata.shape}')

        toc = time.perf_counter()
        logger.info(f'Projection Elapsed Time: {toc - tic}')

        #Filter Known Modes
        if onestep_kwargs.filter_known_modes_kwargs:
            logger.info('Filtering Known Modes')
            tic = time.perf_counter()

            projected_adata = tc.filter_known_modes(projected_adata, 
                                                    **onestep_kwargs.filter_known_modes_kwargs)

            logger.info(f'Projected Adata Dimensions after Filtering Known Modes: {projected_adata.shape}')
            toc = time.perf_counter()
            logger.info(f'Filter Known Modes Elapsed Time: {toc - tic}')
        else:
            logger.info('No known modes, skipping Filter Known Modes')
    else:
        logger.info('Extracting latent dims')
        tic = time.perf_counter()

        ## Extract latent dimensions
        projected_adata = tc.latent_project(norm_adata,
                                            **onestep_kwargs.latent_kwargs)

        toc = time.perf_counter()
        logger.info(f'Extracting latent dims Elapsed Time: {toc - tic}')

    #Louvain Clustering
    logger.info('Starting Louvain Clustering')
    tic = time.perf_counter()

    cluster_louvain_kwargs = onestep_kwargs.cluster_louvain_kwargs.copy()
    k = cluster_louvain_kwargs.pop('k', 15)
    k = min(k, math.floor(projected_adata.n_vars))
    cluster_by_obs, obs_by_cluster, graph, qc = tc.cluster_louvain(
        projected_adata,
        k=k,
        random_seed=random_seed,
        **cluster_louvain_kwargs
    )

    logger.info(f'Completed Louvain Clustering, found {len(obs_by_cluster.keys())} clusters')
    toc = time.perf_counter()
    logger.info(f'Clustering Elapsed Time: {toc - tic}')

    #Merging
    cluster_sizes_before_merging = {k: len(v) for k, v in obs_by_cluster.items()}
    logger.info('Starting Cluster Merging')
    tic = time.perf_counter()

    cluster_assignments_after_merging, marker_genes = tc.merge_clusters(
        adata_norm=norm_adata,
        adata_reduced=projected_adata,
        cluster_assignments=obs_by_cluster,
        cluster_by_obs=cluster_by_obs,
        **onestep_kwargs.merge_clusters_kwargs
    )
    
    logger.info('Completed Cluster Merging')
    toc = time.perf_counter()
    logger.info(f'Merging Elapsed Time: {toc - tic}')

    logger.info('Completed One Step Clustering')
    return list(cluster_assignments_after_merging.values()), marker_genes
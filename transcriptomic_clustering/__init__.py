# -*- coding: utf-8 -*-

"""Top-level package for transcriptomic_clustering."""

import os
from logging.config import fileConfig


version_path = os.path.join(os.path.dirname(__file__), "VERSION.txt")
with open(version_path, "r") as version_file:
    __version__ = version_file.read().strip()


fileConfig(os.path.join(
    os.path.dirname(__file__), 
    'logging_config.ini')
)

# paths
from .utils.memory import memory
from .normalization import normalize
from .highly_variable_genes import highly_variable_genes
from .means_vars_genes import get_means_vars_genes
from .dimension_reduction import pca
from .projection import project, latent_project
from .clustering import cluster_louvain, cluster_louvain_phenograph
from .filter_known_modes import filter_known_modes
from .hierarchical_sorting import hclust
from .cluster_means import get_cluster_means, get_cluster_means_per_batch
from .merging import merge_clusters
from .diff_expression import de_pairs_chisq, vec_chisq_test
from .de_ebayes import de_pairs_ebayes, de_pairs_ebayes_parallel
from .final_merging import final_merge
from .de_all_pairs import (
    check_pairs_ds, check_pairs_lfc, create_pairs, de_all_pairs, get_gene_score_ds, get_pairs,
    load_cl_bin, make_cl_bin, read_de_pairs, select_markers_pair_direction_ds,
    select_markers_pair_group_ds, select_markers_pair_group_top_ds, select_pos_markers_ds,
    select_top_pos_markers_ds,
)
from .post_clustering_qc import (
    check_triplet, find_doublet_by_marker, find_doublets, find_low_quality, find_triplets,
)

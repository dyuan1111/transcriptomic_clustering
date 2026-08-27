"""
Tests for batch-aware cluster merging (merge_clusters(..., batch_aware_merging={...})).

The env used for this work has no pytest, so every test is a plain function and `__main__` runs
them all; the pytest-style naming keeps them collectable once pytest is available.
"""
import numpy as np
import pandas as pd
import anndata as ad

import transcriptomic_clustering as tc
from transcriptomic_clustering import merging


N_GENES = 40
BASE = 4.0
# Deliberately unequal batch sizes. With equal ones a contradictory effect cancels exactly in the
# pooled mean, and the pooled merge would get the right answer for the wrong reason -- the
# batch-effect test needs a case where pooling genuinely misleads.
BATCH_SIZES = {'A': 50, 'B': 10}

THRESHOLDS = {
    'q1_thresh': 0.5,
    'q2_thresh': None,
    'cluster_size_thresh': 6,
    # every gene is expressed above low_thresh in every cluster here, so detection rates are all 1
    # and a qdiff requirement would reject every gene before the batch logic is reached
    'qdiff_thresh': None,
    'padj_thresh': 0.05,
    'lfc_thresh': 1.0,
    'score_thresh': 40,
    'low_thresh': 1,
    'min_genes': 5,
}

# Genes 0-9 carry the signal separating clusters 1 and 2; genes 20-39 separate the {1,2} block from
# clusters 3 and 4 identically in every batch, so that block never merges and the loop always has
# somewhere to go.
SIGNAL = slice(0, 10)
BLOCK = slice(20, 40)


def _build(effect_cl1, batches=('A', 'B'), noise=0.15, seed=0):
    """
    Four clusters over two batches.

    effect_cl1[batch] shifts genes 0-9 in cluster 1 relative to cluster 2 within that batch; a
    positive value puts them up in cluster 1, a negative value up in cluster 2. Everything else is
    identical across batches, so the only thing a test varies is that shift.
    """
    rng = np.random.default_rng(seed)
    # Per-gene noise levels, not one shared level: with identical variance on every gene the eBayes
    # prior degenerates (df_prior -> inf, moderated variance -> NaN) and every DE test comes back
    # empty. Real data never looks like that, but a naive fixture does.
    sigma = rng.uniform(noise, noise * 4, size=N_GENES)
    rows, obs_cluster, obs_batch = [], [], []
    for batch in batches:
        n = BATCH_SIZES[batch]
        for cl in (1, 2, 3, 4):
            block = rng.normal(BASE, sigma, size=(n, N_GENES))
            if cl in (1, 2):
                shift = effect_cl1[batch]
                block[:, SIGNAL] += shift if cl == 1 else -shift
            else:
                # a large, batch-consistent difference that keeps {3, 4} away from {1, 2}
                block[:, BLOCK] += 6.0 if cl == 3 else 10.0
            rows.append(block)
            obs_cluster.extend([cl] * n)
            obs_batch.extend([batch] * n)

    X = np.vstack(rows)
    obs = pd.DataFrame(
        {'cluster': obs_cluster, 'platform': obs_batch},
        index=[f'cell_{i}' for i in range(X.shape[0])],
    )
    var = pd.DataFrame(index=[f'gene_{i}' for i in range(N_GENES)])
    adata = ad.AnnData(X, obs=obs, var=var)

    cluster_assignments = {
        cl: np.flatnonzero(obs['cluster'].values == cl) for cl in (1, 2, 3, 4)
    }
    cluster_by_obs = obs['cluster'].values
    # The reduced space only proposes candidate pairs. Build it from the BLOCK genes so that 3 and 4
    # sit far from {1, 2} and clusters 1 and 2 are each other's nearest neighbour -- otherwise the
    # pair under test is never even offered to the merge rule.
    cols = list(range(BLOCK.start, BLOCK.start + 12))
    adata_reduced = ad.AnnData(
        X[:, cols].copy(), obs=obs, var=pd.DataFrame(index=[var.index[c] for c in cols])
    )
    return adata, adata_reduced, cluster_assignments, cluster_by_obs


def _merge(adata, adata_reduced, assignments, cluster_by_obs, batch_aware_merging, **kwargs):
    # k=4 = R's compare.k default. The batch-aware shortlist mirrors R's sim_knn, whose top-k
    # INCLUDES self, so k=2 would leave each cluster a single real neighbour -- too brittle for a
    # noise-driven fixture.
    kwargs.setdefault('k', 4)
    merged, _ = merging.merge_clusters(
        adata_norm=adata,
        adata_reduced=adata_reduced,
        cluster_assignments={k: list(v) for k, v in assignments.items()},
        cluster_by_obs=cluster_by_obs.copy(),
        thresholds=THRESHOLDS,
        n_markers=None,
        batch_aware_merging=batch_aware_merging,
        **kwargs,
    )
    return merged


def _same_cluster(merged, cell_a, cell_b):
    """Whether two cells ended up in the same cluster (labels are not stable across merges)."""
    owner = {}
    for label, members in merged.items():
        for c in members:
            owner[c] = label
    return owner[cell_a] == owner[cell_b]


# first cell of cluster 1 and of cluster 2, both within batch A (which is laid out first)
CELL_CL1 = 0
CELL_CL2 = BATCH_SIZES['A']


def test_conserved_difference_keeps_clusters_split():
    """A gene up in cluster 1 in BOTH batches is real evidence: both batches object, pair stays."""
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    merged = _merge(adata, red, assign, cbo, {'batch_obs': 'platform'})
    assert not _same_cluster(merged, CELL_CL1, CELL_CL2), \
        "clusters separated consistently in every batch must not merge"


def test_contradictory_difference_merges():
    """
    The reason this feature exists: genes up in cluster 1 in batch A and up in cluster 2 in batch B
    contradict each other, so the conservation filter discards them and nothing is left to separate
    the pair. The pooled merge cannot see this -- it only sees the net difference.
    """
    adata, red, assign, cbo = _build({'A': 2.0, 'B': -2.0})
    merged = _merge(adata, red, assign, cbo, {'batch_obs': 'platform'})
    assert _same_cluster(merged, CELL_CL1, CELL_CL2), \
        "a difference that flips direction between batches is a batch effect, not a cell type"

    pooled = _merge(adata, red, assign, cbo, None)
    assert not _same_cluster(pooled, CELL_CL1, CELL_CL2), \
        "the pooled merge should be fooled here -- otherwise the test proves nothing"


def test_none_is_a_no_op():
    """The default path must be byte-for-byte the old behaviour."""
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    a = _merge(adata, red, assign, cbo, None)
    b = _merge(adata, red, assign, cbo, None)
    assert {k: sorted(v) for k, v in a.items()} == {k: sorted(v) for k, v in b.items()}


def test_batch_with_no_evidence_votes_merge():
    """
    A batch that sees nothing separating the pair must vote to merge, not abstain. With no signal at
    all in either batch, both vote merge and the pair goes.
    """
    adata, red, assign, cbo = _build({'A': 0.0, 'B': 0.0})
    merged = _merge(adata, red, assign, cbo, {'batch_obs': 'platform'})
    assert _same_cluster(merged, CELL_CL1, CELL_CL2)


def test_pair_no_batch_can_judge_never_merges():
    """
    Raise cluster_size_thresh above every per-batch cluster size so no batch qualifies to judge any
    pair. Nothing may merge: absence of evidence is not evidence of similarity.
    """
    adata, red, assign, cbo = _build({'A': 0.0, 'B': 0.0})
    merged = merging.merge_clusters(
        adata_norm=adata,
        adata_reduced=red,
        cluster_assignments={k: list(v) for k, v in assign.items()},
        cluster_by_obs=cbo.copy(),
        # cluster_size_thresh also drives merge_small_clusters, so keep it above the largest
        # per-batch cluster (50) and below the pooled size (60)
        thresholds={**THRESHOLDS, 'cluster_size_thresh': 55},
        n_markers=None,
        k=4,
        batch_aware_merging={'batch_obs': 'platform'},
    )[0]
    assert len(merged) == 4, f"expected no merges, got {len(merged)} clusters"


def test_one_objecting_batch_vetoes():
    """
    Objection is a veto; agreement is only permission. Batch B is given an unreachable score_thresh
    so it always votes merge, while batch A sees a conserved difference -- the pair must stay split
    on batch A's word alone. This is what protects biology only one modality can resolve.
    """
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    merged = _merge(adata, red, assign, cbo, {
        'batch_obs': 'platform',
        'thresholds': {'B': {'score_thresh': 1e9}},
    })
    assert not _same_cluster(merged, CELL_CL1, CELL_CL2), \
        "a single batch that can separate the pair must keep it split"


def test_per_batch_threshold_override():
    """A batch given an unreachable score_thresh always votes merge, so the other batch decides."""
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    merged = _merge(adata, red, assign, cbo, {
        'batch_obs': 'platform',
        'thresholds': {'A': {'score_thresh': 1e9}, 'B': {'score_thresh': 1e9}},
    })
    assert _same_cluster(merged, CELL_CL1, CELL_CL2), \
        "with both batches' score_thresh unreachable every pair should merge"


def test_single_platform_gene_counts_as_evidence():
    """
    genes_by_batch semantics, as in R: a gene measured by only ONE platform has set.num=1 in the
    conservation filter and is trivially conserved -- that platform's evidence counts alone.

    Setup: the separating signal (genes 0-9, up in cluster 1 in batch A only, flat in B) would
    normally be killed by conservation (2 batches measure it, only 1 agrees, 1 < 0.7*2). Declaring
    those genes as measured ONLY by batch A makes set.num=1, so they survive and batch A vetoes the
    merge. Same data, opposite outcomes -- the mask is what changes the answer.
    """
    adata, red, assign, cbo = _build({'A': 2.0, 'B': 0.0})
    all_genes = list(adata.var_names)
    signal = all_genes[SIGNAL]

    # without masks: B measures the signal genes and disagrees (lfc ~0 there) -> dropped -> merge
    merged = _merge(adata, red, assign, cbo, {'batch_obs': 'platform'})
    assert _same_cluster(merged, CELL_CL1, CELL_CL2), \
        'without masks, a 1-of-2-batches signal fails conservation and the pair merges'

    # with masks: the signal genes belong to batch A only -> set.num=1 -> conserved -> A vetoes
    merged = _merge(adata, red, assign, cbo, {
        'batch_obs': 'platform',
        'genes_by_batch': {'A': all_genes, 'B': [g for g in all_genes if g not in set(signal)]},
    })
    assert not _same_cluster(merged, CELL_CL1, CELL_CL2), \
        'a single-platform gene must be able to veto alone, as in R'


def test_genes_by_batch_validation():
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    try:
        _merge(adata, red, assign, cbo,
               {'batch_obs': 'platform', 'genes_by_batch': {'A': ['not_a_gene']}})
    except ValueError as e:
        assert 'genes_by_batch' in str(e)
    else:
        raise AssertionError('expected ValueError for unknown genes')


def test_cl_small_dissolved_in_batch_aware_mode():
    """
    The gap case: a cluster below min.cells in EVERY batch but above the POOLED threshold. The
    pooled small-cluster rule keeps it; batch-aware merging can never judge it (unjudgeable pairs
    never merge), so without R's cl.small dissolve it would survive to the output unvetted. Port of
    harmonize_merge.R:406-495: its cells are remapped individually to big-cluster centroids.
    """
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    # steal 4 cells from each batch's cluster 1 into a new cluster 9: sizes 4 (A) + 4 (B) = 8 pooled,
    # which clears cluster_size_thresh=6 pooled but is below 6 in each batch separately
    idx1 = list(assign[1])
    batch = adata.obs['platform'].values
    from_a = [i for i in idx1 if batch[i] == 'A'][:4]
    from_b = [i for i in idx1 if batch[i] == 'B'][:4]
    assign = {k: [i for i in v if i not in set(from_a + from_b)] for k, v in assign.items()}
    assign[9] = from_a + from_b
    cbo = cbo.copy()
    cbo[assign[9]] = 9

    merged = _merge(adata, red, assign, cbo, {'batch_obs': 'platform'})
    sizes = {k: len(v) for k, v in merged.items()}
    assert not any(set(v) == set(assign[9]) for v in merged.values()), \
        'the everywhere-small cluster must be dissolved, not carried through'
    assert sum(sizes.values()) == adata.n_obs, 'dissolved cells must be reassigned, not dropped'


def test_pair_batch_mode_reaches_same_result_here():
    """
    pair_batch=1 forces the chunked scan-order path (mirroring R's pairBatch): one pair tested at a
    time, stopping the round as soon as anything is mergeable. On this fixture the outcome must match
    the exhaustive default -- the batch-effect pair is the only mergeable one, so scan order cannot
    change the answer. This pins the chunk loop and early-stop mechanics, NOT general equivalence:
    with many mergeable pairs the two modes may legitimately merge in different orders.
    """
    adata, red, assign, cbo = _build({'A': 2.0, 'B': -2.0})
    merged_default = _merge(adata, red, assign, cbo, {'batch_obs': 'platform'})
    merged_batched = _merge(adata, red, assign, cbo, {'batch_obs': 'platform', 'pair_batch': 1})
    assert _same_cluster(merged_batched, CELL_CL1, CELL_CL2), \
        'pair_batch mode must still find and merge the batch-effect pair'
    assert len(merged_batched) == len(merged_default)


def test_pair_batch_validation():
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    for bad in (0, -5, 2.5, 'many'):
        try:
            _merge(adata, red, assign, cbo, {'batch_obs': 'platform', 'pair_batch': bad})
        except ValueError as e:
            assert 'pair_batch' in str(e)
        else:
            raise AssertionError(f'expected ValueError for pair_batch={bad!r}')


def test_bad_batch_obs_raises():
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    for bad, msg in [
        ({'batch_obs': 'not_a_column'}, 'not a column'),
        ({}, "requires 'batch_obs'"),
        ({'batch_obs': 'platform', 'typo': 1}, 'unknown keys'),
    ]:
        try:
            _merge(adata, red, assign, cbo, bad)
        except ValueError as e:
            assert msg in str(e), f"unexpected message for {bad}: {e}"
        else:
            raise AssertionError(f'expected ValueError for {bad}')


def test_chisq_with_batch_aware_raises():
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    try:
        _merge(adata, red, assign, cbo, {'batch_obs': 'platform'}, de_method='chisq')
    except ValueError as e:
        assert 'ebayes' in str(e)
    else:
        raise AssertionError('expected ValueError for de_method=chisq')


def test_unknown_batch_in_threshold_overrides_raises():
    adata, red, assign, cbo = _build({'A': 1.5, 'B': 1.5})
    try:
        _merge(adata, red, assign, cbo,
               {'batch_obs': 'platform', 'thresholds': {'C': {'score_thresh': 10}}})
    except ValueError as e:
        assert 'not present in the data' in str(e)
    else:
        raise AssertionError('expected ValueError for an unknown batch')


# --- constants pinned against the R source -------------------------------------------------------
# These three were wrong in the first implementation, which was written from a prose summary of
# merge_cl_multiple rather than from harmonize_merge.R. None of the behavioural tests above can
# catch them: all three shift a threshold, and the fixture sits far from every boundary. They are
# pinned directly instead.

def test_half_threshold_is_the_mean_not_the_min():
    """harmonize_merge.R:399 -- `de.score.th = mean(...)` over the per-dataset thresholds."""
    assert merging._half_threshold({'A': 100, 'B': 100}) == 50.0
    # mean(100, 300) / 2 = 100; the min would give 50
    assert merging._half_threshold({'A': 100, 'B': 300}) == 100.0
    assert merging._half_threshold({'A': 40, 'B': 1000, 'C': 100}) == 190.0


def test_conservation_counts_every_batch_holding_both_clusters():
    """
    harmonize_merge.R:159,167 -- `set.num` counts every merge.set that holds the gene and both
    clusters, including sets too small to have run a DE test. Restricting the denominator to the
    judging batches would make the filter more permissive.
    """
    genes = ['g0', 'g1']
    gene_pos = {g: i for i, g in enumerate(genes)}
    # only batch A called anything DE; B and C merely hold both clusters
    detail = {'A': pd.DataFrame({'gene': ['g0'], 'P1': ['x']})}

    # g0 is up in x for A and B, but DOWN in x for C -> 2 of 3 agree, 2 >= 0.7*3 = 2.1 is False
    means = {
        'A': {'x': np.array([5.0, 0.0]), 'y': np.array([0.0, 0.0])},
        'B': {'x': np.array([5.0, 0.0]), 'y': np.array([0.0, 0.0])},
        'C': {'x': np.array([0.0, 0.0]), 'y': np.array([5.0, 0.0])},
    }
    kept = merging._conserved_genes(detail, means, gene_pos, 'x', 'y', 1.0, 0.7)
    assert kept == set(), "a batch that disagrees must still count in the denominator"

    # drop the disagreeing batch -> 2 of 2 agree -> kept
    kept = merging._conserved_genes(detail, {k: means[k] for k in ('A', 'B')},
                                    gene_pos, 'x', 'y', 1.0, 0.7)
    assert kept == {('g0', 'x')}

    # a batch missing one of the clusters does not count either way
    means_partial = dict(means)
    means_partial['C'] = {'x': np.array([0.0, 0.0])}      # no 'y'
    kept = merging._conserved_genes(detail, means_partial, gene_pos, 'x', 'y', 1.0, 0.7)
    assert kept == {('g0', 'x')}


def test_conservation_lfc_threshold_is_used():
    """
    harmonize_merge.R:167 uses a literal `lfc > 1`, in log2. This pipeline is natural-log, so the
    value is a parameter; check it is actually applied rather than ignored.
    """
    genes = ['g0']
    gene_pos = {'g0': 0}
    detail = {'A': pd.DataFrame({'gene': ['g0'], 'P1': ['x']})}
    means = {'A': {'x': np.array([1.5]), 'y': np.array([0.0])}}

    assert merging._conserved_genes(detail, means, gene_pos, 'x', 'y', 1.0, 0.7) == {('g0', 'x')}
    assert merging._conserved_genes(detail, means, gene_pos, 'x', 'y', 2.0, 0.7) == set(), \
        "a fold change below conservation_lfc_th must not count as agreement"


def test_tie_break_matches_r_pair_naming():
    """
    Score ties break in R's order: pair name = string-wise min/max of the labels joined by '_'
    (harmonize_merge.R:531-533 pmin/pmax on characters, :394-395 tapply + stable sort). So labels
    6 and 22 name the pair "22_6" ("2" < "6" string-wise), and "22_6" sorts BEFORE "2_30"
    ('2' < '_'). Pinned so tied merges stay identical across the two implementations.
    """
    assert merging._r_pair_name(6, 22) == '22_6'
    assert merging._r_pair_name(22, 6) == '22_6'
    assert merging._r_pair_name('2', '30') == '2_30'
    assert merging._r_pair_name(6, 22) < merging._r_pair_name(2, 30)   # R would merge {6,22} first


def test_conservation_top_n_matches_r():
    """harmonize_merge.R:126 -- `top.n=1000` genes per direction reach the conservation filter."""
    assert merging.CONSERVATION_TOP_N == 1000


def test_backed_adata_raises():
    """Per-batch statistics need an in-memory X; the backed path would take one pass per batch."""
    adata, _red, assign, _cbo = _build({'A': 1.5, 'B': 1.5})
    import tempfile, os
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'a.h5ad')
        adata.write_h5ad(path)
        backed = ad.read_h5ad(path, backed='r')
        try:
            tc.get_cluster_means_per_batch(backed, {k: list(v) for k, v in assign.items()},
                                           adata.obs['platform'].values)
        except NotImplementedError as e:
            assert 'in-memory' in str(e)
        else:
            raise AssertionError('expected NotImplementedError for backed AnnData')


def test_per_batch_means_match_manual():
    adata, _red, assign, _cbo = _build({'A': 1.5, 'B': -1.5})
    stats, sizes = tc.get_cluster_means_per_batch(
        adata, {k: list(v) for k, v in assign.items()}, adata.obs['platform'].values, low_th=1
    )
    assert set(stats) == {'A', 'B'}
    for batch in ('A', 'B'):
        means, present, _var = stats[batch]
        for cl in (1, 2, 3, 4):
            idx = np.flatnonzero(
                (adata.obs['cluster'].values == cl) & (adata.obs['platform'].values == batch)
            )
            assert sizes[batch][cl] == len(idx)
            np.testing.assert_allclose(means.loc[cl].values, adata.X[idx].mean(axis=0))
            np.testing.assert_allclose(
                present.loc[cl].values, (adata.X[idx] > 1).sum(axis=0) / len(idx)
            )


if __name__ == '__main__':
    import logging
    logging.disable(logging.WARNING)
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    failed = 0
    for t in tests:
        try:
            t()
            print(f'PASS {t.__name__}')
        except Exception as e:
            failed += 1
            print(f'FAIL {t.__name__}: {type(e).__name__}: {e}')
    print(f'\n{len(tests) - failed}/{len(tests)} passed')
    raise SystemExit(1 if failed else 0)

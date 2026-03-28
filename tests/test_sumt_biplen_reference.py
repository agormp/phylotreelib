\
"""
Reference tests for majority-rule consensus + mean bipartition branch lengths.

Goal
----
Catch regressions in how sumt/phylotreelib:
  - counts bipartition frequencies (majority-rule cutoff)
  - computes mean branch length for each bipartition (conditional on presence)

These tests deliberately DO NOT call TreeSummary._add_bips / curtree.iter_bipinfo()
for the reference computation. Instead, they compute splits directly from each Tree's
(parent->child) edges using remote-children bitmasks.

Notes
-----
* We still use phylotreelib to parse trees from file. If parsing is broken, many other
  tests should already fail; here we focus on the split/length aggregation logic.
* The "reference" here is an independent implementation, not MrBayes/PAUP output.
  This avoids a brittle external dependency while still providing a strong guard
  against internal regressions.

Data
----
Uses tests/data/random_10tips_1000.trees (already referenced by other tests).
"""

from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Tuple, Iterable

import pytest
import phylotreelib as pt


# -----------------------------
# Helpers: split canonicalization
# -----------------------------

def _popcount(x):
    return x.bit_count()


def _allmask(n_leaves):
    return (1 << n_leaves) - 1


def _leaf_mask(tree: "pt.Tree", leaf):
    # Tree.cached_attributes typically includes a deterministic leaf ordering,
    # but leaf2index is already consistent with remotechildren_mask_dict usage.
    return 1 << tree.leaf2index[leaf]


def _node_mask(tree, node):
    """
    Return mask for any node (internal node id or leaf name).
    remotechildren_mask_dict usually includes masks for internal nodes; for leaves
    we fall back to leaf2index.
    """
    rem = tree.remotechildren_mask_dict  # ensures it's computed/cached
    try:
        return int(rem[node])
    except KeyError:
        # likely a leaf
        return _leaf_mask(tree, node)


def _canon_halfmask(mask, fullmask):
    """
    Canonicalize an unrooted split by taking the "smaller side" of the bipartition.
    """
    comp = fullmask ^ mask
    # choose by size first, then by integer value for deterministic tie-breaking
    a, b = mask, comp
    pa, pb = _popcount(a), _popcount(b)
    if pa < pb:
        return a
    if pb < pa:
        return b
    return min(a, b)


def iter_internal_splits_with_lengths(tree):
    """
    Yield (halfmask, length) for each INTERNAL split in the tree.

    Internal split definition here: both sides have at least 2 leaves.
    """
    n = len(tree.leaves)
    full = _allmask(n)

    # iterate all directed edges parent->child; each corresponds to one undirected split
    child_dict = tree.child_dict
    for parent in tree.intnodes:
        for child, br in child_dict[parent].items():
            m = _node_mask(tree, child)
            hm = _canon_halfmask(m, full)
            k = _popcount(hm)
            if 2 <= k <= n - 2:
                yield (hm, float(br.length))


@dataclass
class SplitAgg:
    n_trees: int
    count: Dict[int, int]        # halfmask -> number of trees where split occurs
    sumlen: Dict[int, float]     # halfmask -> sum of lengths across trees where split occurs

    def meanlen(self, hm):
        return self.sumlen[hm] / self.count[hm]

    def freq(self, hm):
        return self.count[hm] / self.n_trees


def reference_aggregate(file_path, burnin= 0):
    """
    Reference aggregation: counts/frequencies + mean lengths conditional on presence.
    """
    tf = pt.Nexustreefile(file_path)
    # discard burnin trees deterministically
    for _ in range(burnin):
        tf.readtree(returntree=False)

    count: Dict[int, int] = {}
    sumlen: Dict[int, float] = {}
    n = 0

    for t in tf:
        n += 1
        # ensure masks are available
        _ = t.remotechildren_mask_dict

        seen = set()
        for hm, blen in iter_internal_splits_with_lengths(t):
            # split should occur once per tree; guard against any accidental duplicates
            if hm in seen:
                continue
            seen.add(hm)
            count[hm] = count.get(hm, 0) + 1
            sumlen[hm] = sumlen.get(hm, 0.0) + blen

    if n == 0:
        raise RuntimeError("No trees read in reference_aggregate()")

    return SplitAgg(n_trees=n, count=count, sumlen=sumlen)


def build_sumt_sumtree(file_path, burnin= 0):
    tf = pt.Nexustreefile(file_path)
    for _ in range(burnin):
        tf.readtree(returntree=False)

    ts = pt.TreeSummary(trackbips=True, trackblen=True)
    for t in tf:
        ts.add_tree(t)

    sumtree = pt.build_sumtree(ts, treetype="con", blen="biplen", rooting=None, og=None,
                              count_burnin_filename_list=None)
    return sumtree, ts


def iter_sumtree_internal_splits_with_lengths(sumtree):
    """
    Internal splits in the *summary* tree (should be fully deterministic).
    """
    n = len(sumtree.leaves)
    full = _allmask(n)
    cd = sumtree.child_dict
    for parent in sumtree.intnodes:
        for child, br in cd[parent].items():
            m = _node_mask(sumtree, child)
            hm = _canon_halfmask(m, full)
            k = _popcount(hm)
            if 2 <= k <= n - 2:
                yield (hm, float(br.length))


# -----------------------------
# Tests
# -----------------------------

def test_biplen_means_match_reference_on_majority_splits():
    """
    For a fixed tree sample file:
      - compute reference mean lengths for each split (conditional on presence)
      - compute sumt's summary tree (--con --biplen)
      - compare branch lengths split-by-split
    """
    data_dir = Path(__file__).with_name("data")
    f = data_dir / "neanderthal.nexus.run1.t"
    ref = reference_aggregate(f, burnin=0)
    sumtree, ts = build_sumt_sumtree(f, burnin=0)
    sorted_leaves = sorted(ts.leaves)

    # Majority-rule set from reference
    majority = {hm for hm, c in ref.count.items() if (c / ref.n_trees) >= 0.5}

    # Compare only splits that appear in the summary tree AND are majority by reference
    sumtree_splits = dict(iter_sumtree_internal_splits_with_lengths(sumtree))
    common = majority.intersection(sumtree_splits.keys())
    assert common, "No common majority splits found; test data unexpected"

    # Tight tolerance: these are computed from the same underlying float lengths.
    # If this fails, it's likely a real aggregation logic regression.
    for hm in sorted(common):
        got = sumtree_splits[hm]
        exp = ref.meanlen(hm)
        assert got == pytest.approx(exp, rel=0, abs=1e-12), f"split {hm}: got {got} exp {exp}"

    # Sanity: TreeSummary counts for lengths should be conditional on presence, i.e. br.n == count[hm]
    # (TreeSummary stores leaf + internal splits; we only check internal ones here.)
    for hm in common:
        # Build a mapping from TreeSummary bipartitions -> halfmask
        # We derive halfmask from the actual bipartition's leaf sets.
        br = None
        for (b1, b2), s in ts.bipartsummary.items():
            # internal only
            if len(b1) == 1 or len(b2) == 1:
                continue
            # convert leaf set to mask using ts.leaves ordering (should match tree order)
            full = _allmask(len(ts.leaves))
            m = 0
            for leaf in b1:
                m |= 1 << sorted_leaves.index(leaf)  # ok for small n
            hm2 = _canon_halfmask(m, full)
            if hm2 == hm:
                br = s
                break
        assert br is not None, "Could not map TreeSummary bipartition to reference halfmask"
        assert br.n == ref.count[hm], f"split {hm}: TreeSummary n={br.n} ref count={ref.count[hm]}"


def test_consensus_topology_matches_reference_majority_set():
    """
    The consensus topology should contain exactly the set of majority internal splits.
    """
    data_dir = Path(__file__).with_name("data")
    f = data_dir / "neanderthal.nexus.run1.t"
    ref = reference_aggregate(f, burnin=0)
    sumtree, _ = build_sumt_sumtree(f, burnin=0)

    majority = {hm for hm, c in ref.count.items() if (c / ref.n_trees) >= 0.5}
    sumtree_splits = set(dict(iter_sumtree_internal_splits_with_lengths(sumtree)).keys())

    # In a majority-rule consensus, internal splits should be exactly the majority set.
    # (Leaves are always present but excluded here by construction.)
    assert sumtree_splits == majority

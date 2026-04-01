import copy
import csv
import itertools
import math
import random
import textwrap
from dataclasses import dataclass
from string import ascii_lowercase, digits
from types import SimpleNamespace
from typing import Dict

import numpy as np
import phylotreelib as pt
import pytest


def approx_places(value, places=7):
    return pytest.approx(value, abs=10 ** (-places))


@pytest.fixture(autouse=True)
def _inject_treedata(request, treedata):
    if request.instance is not None:
        request.instance.treedata = treedata


# Combined from test_phylotreelib_migrated.py

class TestTreeRead:
    """Tests methods for reading (not iterating) trees from treefiles"""

    def test_readtree_newick(self, tmp_path):
        """Test that readtree returns correct trees and None when Newickfile exhausted"""

        filename = tmp_path / "trees.newick"
        treelist = []
        trees = pt.TreeSet()
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            trees.addtree(mytree)
        filename.write_text(trees.newick(), encoding="utf-8")

        ntrees = len(treelist)
        treefile = pt.Newicktreefile(filename)
        for i in range(ntrees):
            tree = treefile.readtree()
            assert tree == treelist[i]

        value = treefile.readtree()
        assert value is None

    def test_readtree_nexus(self, tmp_path):
        """Test that readtree returns correct trees and None when Nexusfile exhausted"""

        filename = tmp_path / "trees.nexus"
        treelist = []
        with filename.open("w", encoding="utf-8") as filehandle:
            for treestring in self.treedata.values():
                mytree = pt.Tree.from_string(treestring)
                treelist.append(mytree)
                filehandle.write(mytree.nexus())

        treefile = pt.Nexustreefile(filename)
        for i, tree in enumerate(treefile):
            assert tree == treelist[i]

    def test_readtrees_newick(self, tmp_path):
        """Test that readtrees (plural) returns correct trees from Newickfile"""

        filename = tmp_path / "trees.newick"
        treelist = []
        trees = pt.TreeSet()
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            trees.addtree(mytree)
        filename.write_text(trees.newick(), encoding="utf-8")

        treefile = pt.Newicktreefile(filename)
        treeset_from_file = treefile.readtrees()
        for origtree, readtree in zip(treelist, treeset_from_file):
            assert origtree == readtree

    def test_readtrees_nexus(self, tmp_path):
        """Test that readtrees (plural) returns correct trees from Nexusfile"""

        filename = tmp_path / "trees.nexus"
        treelist = []
        with filename.open("w", encoding="utf-8") as filehandle:
            for treestring in self.treedata.values():
                mytree = pt.Tree.from_string(treestring)
                treelist.append(mytree)
                filehandle.write(mytree.nexus())

        treefile = pt.Nexustreefile(filename)
        newlyreadtrees = []
        for tree in treefile:
            newlyreadtrees.append(tree)
        for origtree, readtree in zip(treelist, newlyreadtrees):
            assert origtree == readtree


class TestOutGroupRooting:
    """Tests for rootout function"""

    def test_ogroot_treestrings(self):
        for treestring in self.treedata.values():
            t = pt.Tree.from_string(treestring)
            for node in t.intnodes - {t.root}:
                tc = t.copy_treeobject()
                ig = tc.remote_children(node)
                og = tc.leaves - ig
                tc.rootout(og)
                rootkid1, rootkid2 = tc.children(tc.root)
                rootremkids1 = tc.remote_children(rootkid1)
                rootremkids2 = tc.remote_children(rootkid2)
                assert ig in [rootremkids1, rootremkids2]
                assert og in [rootremkids1, rootremkids2]


class TestHas_same_root:
    """Tests for has_same_root function"""

    def test_sameroot_sameintnodenumbers(self):
        """Check returns True: two identical trees"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            assert t1.has_same_root(t2)

    def test_sameroot_diffintnodenumbers(self):
        """Check returns True: two trees with same rooted topology, different nodeids"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            delta_id = 35 + max(t1.intnodes)
            for _, node_id in enumerate(t1.intnodes):
                t2.rename_intnode(node_id, node_id + delta_id)
            assert t1.has_same_root(t2)

    def test_diffroot_sameintnodenumbers(self):
        """Check returns False: two trees with same unrooted topology, different roots"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            newroot = random.choice(tuple(t1.intnodes - {t1.root}))
            t2.deroot()
            t2.reroot(newroot, polytomy=True)
            assert not t1.has_same_root(t2)


class TestMatch_nodes:
    """Tests for match_nodes function"""

    def test_sametop_sameroot(self):
        """Test correct output from well formed examples"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            delta_id = 35 + max(t1.intnodes)
            for _, node_id in enumerate(t1.intnodes):
                t2.rename_intnode(node_id, node_id + delta_id)
            intnode1to2, unmatched_root1, unmatched_root2 = t1.match_nodes(t2)
            for id1 in intnode1to2:
                if type(id1) == int:
                    assert id1 + delta_id == intnode1to2[id1]
                else:
                    assert id1 == intnode1to2[id1]
            assert unmatched_root1 is None
            assert unmatched_root2 is None

    def test_sametop_difroot(self):
        """Test correct output from well formed examples with different root"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            if len(t1.intnodes) > 2:
                root1, root2 = random.sample(tuple(t1.intnodes - {t1.root}), 2)
                t1.reroot(root1, polytomy=True)
                t2.reroot(root2, polytomy=True)
                delta_id = 35 + max(t1.intnodes)
                for _, node_id in enumerate(t1.intnodes):
                    t2.rename_intnode(node_id, node_id + delta_id)
                intnode1to2, unmatched_root1, unmatched_root2 = t1.match_nodes(t2)
                for id1 in intnode1to2:
                    if type(id1) == int:
                        assert id1 + delta_id == intnode1to2[id1]
                    else:
                        assert id1 == intnode1to2[id1]
                assert unmatched_root1 is None
                assert unmatched_root2 is None

    def test_unmatchedroots(self):
        """Test that unmatched roots are reported correctly"""
        t1 = pt.Tree.randtree(ntips=50)
        t2 = t1.copy_treeobject(copylengths=True)
        delta_id = 35 + max(t1.intnodes)
        for _, node_id in enumerate(t1.intnodes):
            t2.rename_intnode(node_id, node_id + delta_id)
        t2.deroot()
        parent = random.choice(tuple(t2.intnodes))
        kid = next(iter(t2.children(parent)))
        blen = t2.nodedist(parent, kid)
        t2.reroot(parent, kid, dist_from_node1=blen / 2)
        t1origroot = t1.root
        t2origroot = t2.root
        intnode1to2, unmatched_root1, unmatched_root2 = t1.match_nodes(t2)
        for id1 in intnode1to2:
            if type(id1) == int:
                assert id1 + delta_id == intnode1to2[id1]
            else:
                assert id1 == intnode1to2[id1]
        assert unmatched_root1 == t1origroot
        assert unmatched_root2 == t2origroot

    def test_difleaves(self):
        """Test error raised when leaves differ"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            randleaf = random.choice(tuple(t1.leaves))
            t1.remove_leaf(randleaf)
            with pytest.raises(pt.TreeError):
                t1.match_nodes(t2)

    def test_sameleaves_diftop(self):
        """Test error raised when leaves same but topologies differ"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            while t1.topology() == t2.topology():
                t2.shuffle_leaf_names()
            with pytest.raises(pt.TreeError):
                t1.match_nodes(t2)


class TestRelationShipMethods:
    """Tests methods for determining children, parents, remote children, MRCAs, and basenode"""

    def test_adjacent(self):
        """Consistency and function of .children() and .parent() methods"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            for parent in mytree.intnodes:
                children = mytree.children(parent)
                for child in children:
                    result = mytree.parent(child)
                    assert result == parent

    def test_remote(self):
        """Consistency and function of remote_children, find_mrca, findbasenode"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            assert mytree.remote_children(mytree.root) == mytree.leaves
            for basalnode in mytree.intnodes:
                remotekids = mytree.remote_children(basalnode)
                assert basalnode == mytree.find_mrca(remotekids)
                assert basalnode == mytree.findbasenode(remotekids)
                if basalnode != mytree.root and basalnode not in mytree.children(mytree.root):
                    basal_parent = mytree.parent(basalnode)
                    other_half = mytree.leaves - remotekids
                    assert basal_parent == mytree.findbasenode(other_half)
                if basalnode in mytree.children(mytree.root):
                    basal_siblings = mytree.children(mytree.root) - {basalnode}
                    other_half = mytree.leaves - remotekids
                    if len(basal_siblings) == 1:
                        sib = next(iter(basal_siblings))
                        assert sib == mytree.findbasenode(other_half)
                    else:
                        assert mytree.root == mytree.findbasenode(other_half)
            for leaf in mytree.leaves:
                assert {leaf} == mytree.remote_children(leaf)
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        assert {"s1", "s2", "s3"} == mytree.remote_children(1)


class TestPathMethods:
    """Tests methods for determining nodepath, nodedist, and treelength"""

    def test_nodedist(self):
        """Does nodedist() report correct patristic distance?"""
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)

        assert mytree.nodedist("s1", "s2") == approx_places(0.25)
        assert mytree.nodedist("s1", "s3") == approx_places(0.5)
        assert mytree.nodedist("s1", "s4") == approx_places(0.625)
        assert mytree.nodedist("s1", "s5") == approx_places(0.625)
        assert mytree.nodedist("s2", "s3") == approx_places(mytree.nodedist("s1", "s3"))
        assert mytree.nodedist("s2", "s4") == approx_places(mytree.nodedist("s1", "s4"))
        assert mytree.nodedist("s2", "s5") == approx_places(mytree.nodedist("s1", "s5"))
        assert mytree.nodedist("s3", "s4") == approx_places(0.375)
        assert mytree.nodedist("s3", "s5") == approx_places(0.375)
        assert mytree.nodedist("s4", "s5") == approx_places(0.25)

        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)

        expecteddist = (
            0.008 + 0.03629 + 0.1025 + 0.16576 + 0.16576 + 0.00771 + 0.00386
            + 0.06371 + 0.00584 + 0.00555 + 0.00063 + 0.00616 + 0.04414 + 0.05684
        )
        assert mytree.nodedist("SIVMK", "HV1H3") == approx_places(expecteddist)

        expecteddist = 0.01652 + 0.16576 + 0.16576 + 0.19506
        anc1 = mytree.find_mrca(
            {"HV2BE", "HV2D1", "HV2SB", "HV2S2", "HV2ST", "HV2G1", "HV2RO", "HV2CA", "HV2NZ"}
        )
        assert mytree.nodedist("SIVCZ", anc1) == approx_places(expecteddist)

        expecteddist = 0.01652 + 0.16576 + 0.16576 + 0.05684 + 0.03899 + 0.02325
        anc2 = mytree.find_mrca({"HV1EL", "HV1Z2"})
        assert mytree.nodedist(anc1, anc2) == approx_places(expecteddist)

    def test_treelength(self):
        """Does treelength() return correct value?"""

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        expectedlength = 0.125 + 0.125 + 0.25 + 0.125 + 0.125 + 0.125 + 0.125
        assert mytree.length() == approx_places(expectedlength)

        treestring = self.treedata["string_with_label"]
        mytree = pt.Tree.from_string(treestring)
        expectedlength = 0.101408 + 0.071355 + 0.124263 + 0.009364 + 0.014955
        assert mytree.length() == approx_places(expectedlength)

    def test_nodepath_sanity(self):
        """Sanity check: nodepath(A,B) = reverse of nodepath(B,A)?"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            leaves = mytree.leaves
            intnodes = mytree.intnodes
            for node1 in leaves:
                for node2 in leaves:
                    path12 = mytree.nodepath(node1, node2)
                    path21 = mytree.nodepath(node2, node1)
                    path21.reverse()
                    assert path12 == path21

            for node1 in intnodes:
                for node2 in intnodes:
                    path12 = mytree.nodepath(node1, node2)
                    path21 = mytree.nodepath(node2, node1)
                    path21.reverse()
                    assert path12 == path21

            for node1 in leaves:
                for node2 in intnodes:
                    path12 = mytree.nodepath(node1, node2)
                    path21 = mytree.nodepath(node2, node1)
                    path21.reverse()
                    assert path12 == path21

    def test_nodepath_details(self):
        """Checking parent-child relationships on nodepath from root to leaves"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            leaves = mytree.leaves
            root = mytree.root
            for leaf in leaves:
                path = mytree.nodepath(root, leaf)
                for i in range(len(path) - 1):
                    parent = path[i]
                    child = path[i + 1]
                    assert child in mytree.children(parent)
                    assert parent == mytree.parent(child)


class TestTreeOutput:
    """Tests of string output"""

    def test_newick_output(self):
        """Sanity check: output from newick() can be parsed by Tree.from_string()"""
        for instring in self.treedata.values():
            mytree = pt.Tree.from_string(instring)
            outstring = mytree.newick()
            mytree2 = pt.Tree.from_string(outstring)
            assert mytree == mytree2

    def test_nexus_output(self):
        """Check that output from nexus() can be parsed by Nexustreefile()"""
        for instring in self.treedata.values():
            mytree = pt.Tree.from_string(instring)
            nexus_string = mytree.nexus()
            nexusfile = pt.Nexustreefile(filecontent=nexus_string)
            mytree2 = next(nexusfile)
            assert mytree == mytree2

    def test_contree_nexus_output(self):
        """Check that nexus output from contree result can be parsed by Nexustreefile"""
        treesummary = pt.TreeSummary(trackbips=True, trackblen=True)
        for _ in range(10):
            tree = pt.Tree.randtree(ntips=50)
            treesummary.add_tree(tree)
        stb = pt.SummaryTreeBuilder(treesummary)
        contree = stb.contree()
        nexus_string = contree.nexus()
        nexusfile = pt.Nexustreefile(filecontent=nexus_string)
        mytree = next(nexusfile)
        assert contree == mytree


class TestTreesummarytests:
    @pytest.fixture(autouse=True)
    def _setup(self, data_dir):
        self.t1_fname = data_dir / "mrbayes" / "contest.postburnin.1.t"
        self.t2_fname = data_dir / "mrbayes" / "contest.postburnin.2.t"
        con_fname = data_dir / "mrbayes" / "contest.nexus.con.tre"
        mbres_fname = data_dir / "bip_mean_var.txt"
        mb_trprobs_fname = data_dir / "mrbayes" / "contest.nexus.trprobs"
        cfile = pt.Nexustreefile(con_fname)
        self.mb_contree = cfile.readtree()
        cfile.close()
        with mbres_fname.open(encoding="utf-8") as mbfile:
            mbresults = mbfile.readlines()
        self.mbresdict = {}
        names1, names2, meanvar = mbresults[0].strip().split("|")
        bip1 = names1.strip().split()
        bip2 = names2.strip().split()
        mocktree = SimpleNamespace()
        frozenset_leaves = frozenset(bip1 + bip2)
        sorted_leaf_tup = tuple(sorted(frozenset_leaves))
        leaf2index = {leaf: i for i, leaf in enumerate(sorted_leaf_tup)}
        leaf2mask = {leaf: (1 << i) for leaf, i in leaf2index.items()}
        ntips = len(sorted_leaf_tup)
        alltips_mask = (1 << ntips) - 1
        mocktree.cached_attributes = (
            frozenset_leaves,
            sorted_leaf_tup,
            leaf2index,
            leaf2mask,
            ntips,
            alltips_mask,
        )
        for line in mbresults:
            names1, names2, meanvar = line.strip().split("|")
            bip1 = names1.strip().split()
            bipart = pt.Bipartition.from_leafset(frozenset(bip1), mocktree)
            vals = meanvar.strip().split()
            mean = float(vals[0])
            var = float(vals[1])
            self.mbresdict[bipart] = [mean, var]
        trprobfile = pt.Nexustreefile(mb_trprobs_fname)
        self.trprob_trees = trprobfile.readtrees()
        trprobfile.close()

    def test_contree(self):
        ts = pt.TreeSummary(trackbips=True, trackblen=True)
        tf1 = pt.Nexustreefile(self.t1_fname)
        for t in tf1:
            ts.add_tree(t)
        tf2 = pt.Nexustreefile(self.t2_fname)
        for t in tf2:
            ts.add_tree(t)
        stb = pt.SummaryTreeBuilder(ts)
        own_contree = stb.contree()
        own_bipdict = own_contree.bipdict()
        mb_bipdict = self.mb_contree.bipdict()

        assert self.mb_contree.topology() == own_contree.topology()

        for bip, own_branch in own_bipdict.items():
            mb_mean, mb_var = self.mbresdict[bip]
            assert mb_mean == approx_places(own_branch.length)
            assert mb_var == approx_places(own_branch.length_sd ** 2)

            bip1, bip2 = bip
            if len(bip1) != 1 and len(bip2) != 1:
                mb_branch = mb_bipdict[bip]
                own_freq = round(float(own_branch.posterior), 3)
                mb_freq = float(mb_branch.label)
                assert mb_freq == approx_places(own_freq)

    def test_treesummary_update(self):
        ts1 = pt.TreeSummary(trackbips=True, trackblen=True)
        tf1 = pt.Nexustreefile(self.t1_fname)
        for t in tf1:
            ts1.add_tree(t)
        tf1.close()
        ts2 = pt.TreeSummary(trackbips=True, trackblen=True)
        tf2 = pt.Nexustreefile(self.t2_fname)
        for t in tf2:
            ts2.add_tree(t)
        tf2.close()
        ts1.update(ts2)
        stb = pt.SummaryTreeBuilder(ts1)
        own_contree = stb.contree()
        own_bipdict = own_contree.bipdict()
        mb_bipdict = self.mb_contree.bipdict()

        assert self.mb_contree.topology() == own_contree.topology()

        for bip, own_branch in own_bipdict.items():
            mb_mean, mb_var = self.mbresdict[bip]
            assert mb_mean == approx_places(own_branch.length)
            assert mb_var == approx_places(own_branch.length_sd ** 2)

            bip1, bip2 = bip
            if len(bip1) != 1 and len(bip2) != 1:
                mb_branch = mb_bipdict[bip]
                own_freq = float(own_branch.posterior)
                mb_freq = float(mb_branch.label)
                assert mb_freq == approx_places(own_freq, places=3)

    def test_treesummary_with_topologies(self):
        ts = pt.TreeSummary(trackbips=True, tracktopo=True)
        tf1 = pt.Nexustreefile(self.t1_fname)
        for t in tf1:
            ts.add_tree(t)
        tf2 = pt.Nexustreefile(self.t2_fname)
        for t in tf2:
            ts.add_tree(t)
        for tree in self.trprob_trees:
            mb_topology = tree.topology()
            assert mb_topology in ts.biptoposummary


class TestTopologytests:
    """Tests topology related methods"""

    def test_bipdict(self):
        """Does bipdict() return correct result?"""

        treestring = self.treedata["string_with_label"]
        mytree = pt.Tree.from_string(treestring)
        total_leaves = {"KL0F07689", "KW081_13", "SBC669_26", "YAL016W"}
        mocktree = SimpleNamespace()
        frozenset_leaves = frozenset(total_leaves)
        sorted_leaf_tup = tuple(sorted(frozenset_leaves))
        ntips = len(sorted_leaf_tup)
        alltips_mask = (1 << ntips) - 1
        leaf2index = {leaf: i for i, leaf in enumerate(sorted_leaf_tup)}
        leaf2mask = {leaf: (1 << i) for leaf, i in leaf2index.items()}
        mocktree.cached_attributes = (
            frozenset_leaves,
            sorted_leaf_tup,
            leaf2index,
            leaf2mask,
            ntips,
            alltips_mask,
        )
        bipdict = mytree.bipdict()

        expectedkey = pt.Bipartition.from_leafset(frozenset(["YAL016W", "SBC669_26"]), mocktree)
        expectedlen = 0.124263
        expectedlab = "0.0507"

        assert expectedkey in bipdict
        assert bipdict[expectedkey].length == approx_places(expectedlen)
        assert bipdict[expectedkey].label == expectedlab

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        total_leaves = {"s4", "s1", "s3", "s2", "s5"}
        mocktree = SimpleNamespace()
        frozenset_leaves = frozenset(total_leaves)
        sorted_leaf_tup = tuple(sorted(frozenset_leaves))
        ntips = len(sorted_leaf_tup)
        alltips_mask = (1 << ntips) - 1
        leaf2index = {leaf: i for i, leaf in enumerate(sorted_leaf_tup)}
        leaf2mask = {leaf: (1 << i) for leaf, i in leaf2index.items()}
        mocktree.cached_attributes = (
            frozenset_leaves,
            sorted_leaf_tup,
            leaf2index,
            leaf2mask,
            ntips,
            alltips_mask,
        )
        bipdict = mytree.bipdict()
        expectedkeys = [
            pt.Bipartition.from_leafset(["s1"], mocktree),
            pt.Bipartition.from_leafset(["s2"], mocktree),
            pt.Bipartition.from_leafset(["s3"], mocktree),
            pt.Bipartition.from_leafset(["s4"], mocktree),
            pt.Bipartition.from_leafset(["s5"], mocktree),
            pt.Bipartition.from_leafset(["s4", "s5"], mocktree),
            pt.Bipartition.from_leafset(["s1", "s2"], mocktree),
        ]

        expectedvals = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.25]

        for key in bipdict:
            assert key in expectedkeys

        for i in range(len(expectedkeys)):
            key = expectedkeys[i]
            val = expectedvals[i]
            assert val == approx_places(bipdict[key].length)

    def test_topology(self):
        """Does topology() return expected result?"""

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        topology = mytree.topology()
        total_leaves = {"s4", "s1", "s3", "s2", "s5"}
        mocktree = SimpleNamespace()
        frozenset_leaves = frozenset(total_leaves)
        sorted_leaf_tup = tuple(sorted(frozenset_leaves))
        ntips = len(sorted_leaf_tup)
        alltips_mask = (1 << ntips) - 1
        leaf2index = {leaf: i for i, leaf in enumerate(sorted_leaf_tup)}
        leaf2mask = {leaf: (1 << i) for leaf, i in leaf2index.items()}
        mocktree.cached_attributes = (
            frozenset_leaves,
            sorted_leaf_tup,
            leaf2index,
            leaf2mask,
            ntips,
            alltips_mask,
        )

        expectedtop = frozenset(
            [
                pt.Bipartition.from_leafset(["s4", "s5"], mocktree),
                pt.Bipartition.from_leafset(["s4", "s5", "s3"], mocktree),
                pt.Bipartition.from_leafset(["s4", "s5", "s3", "s2"], mocktree),
                pt.Bipartition.from_leafset(["s4", "s5", "s3", "s1"], mocktree),
                pt.Bipartition.from_leafset(["s4", "s5", "s1", "s2"], mocktree),
                pt.Bipartition.from_leafset(["s1", "s5", "s3", "s2"], mocktree),
                pt.Bipartition.from_leafset(["s4", "s1", "s3", "s2"], mocktree),
            ]
        )

        assert topology == expectedtop

    def test_compatiblewith(self):
        """Does is_compatible_with() return expected result?"""

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)

        expectedbips = [
            frozenset([frozenset(["s4", "s5"]), frozenset(["s1", "s2", "s3"])]),
            frozenset([frozenset(["s4", "s5", "s3"]), frozenset(["s1", "s2"])]),
            frozenset([frozenset(["s4", "s5", "s3", "s2"]), frozenset(["s1"])]),
            frozenset([frozenset(["s4", "s5", "s3", "s1"]), frozenset(["s2"])]),
            frozenset([frozenset(["s4", "s5", "s1", "s2"]), frozenset(["s3"])]),
            frozenset([frozenset(["s1", "s5", "s3", "s2"]), frozenset(["s4"])]),
            frozenset([frozenset(["s4", "s1", "s3", "s2"]), frozenset(["s5"])]),
        ]

        badbips = [
            frozenset([frozenset(["s3", "s5"]), frozenset(["s1", "s2", "s4"])]),
            frozenset([frozenset(["s4", "s1", "s3"]), frozenset(["s5", "s2"])]),
        ]

        for bip in expectedbips:
            assert mytree.is_compatible_with(bip)

        for bip in badbips:
            assert not mytree.is_compatible_with(bip)

    def test_is_resolved(self):
        """Does is_resolved() return expected result?"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            assert mytree.is_resolved(), treestring

        mytree = pt.Tree.from_string("(A, B, C, D);")
        assert not mytree.is_resolved()

        mytree = pt.Tree.from_string("(A, (B, C, D));")
        assert not mytree.is_resolved()

        mytree = pt.Tree.from_string("(A, (B, (C, D), E));")
        assert not mytree.is_resolved()

    def test_resolve(self):
        namelist = [f"seq{i}" for i in range(50)]
        t1 = pt.Tree.from_leaves(namelist)
        t1.resolve()
        t1_string = t1.newick()
        t2 = pt.Tree.from_string(t1_string)
        assert isinstance(t2, pt.Tree)


class TestLengthLabel:
    """Tests setting and getting of branch lengths and labels"""

    def test_setlength(self):
        """Does setlength() work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"HV2BE", "HV2D1"})
        mytree.setlength(anc, "HV2BE", 0.666)
        mytree.setlength(anc, "HV2D1", 0.333)
        assert 0.666 == approx_places(mytree.nodedist(anc, "HV2BE"))
        assert 0.333 == approx_places(mytree.nodedist(anc, "HV2D1"))

    def test_labels(self):
        """Does setlabel() and getlabel() work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"HV2BE", "HV2D1"})
        mytree.setlabel(anc, "HV2BE", "Simpel_label")
        mytree.setlabel(anc, "HV2D1", "Label with space")
        assert mytree.getlabel(anc, "HV2BE") == "Simpel_label"
        assert mytree.getlabel(anc, "HV2D1") == "Label with space"


class TestTreeChanging:
    """Tests for functions that alter tree structure"""

    def test_insertnode(self):
        """Does insert_node() work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"HV2BE", "HV2D1"})
        pre_nodes = copy.copy(mytree.nodes)

        branchstruct = pt.Branchstruct(length=0.777)
        branchstruct.label = "New_branch"
        newnode = mytree.insert_node(anc, ["HV2BE", "HV2D1"], branchstruct)
        assert mytree.nodedist(newnode, anc) == approx_places(0.777)
        assert mytree.getlabel(newnode, anc) == "New_branch"
        assert mytree.parent(newnode) == anc
        assert mytree.children(newnode) == {"HV2BE", "HV2D1"}
        assert mytree.children(anc) == {newnode}
        assert len(pre_nodes) + 1 == len(mytree.nodes)
        assert pre_nodes < mytree.nodes

        anc2 = mytree.parent(anc)
        branchstruct = pt.Branchstruct(length=0.999)
        branchstruct.label = "Lower_branch"
        newnode2 = mytree.insert_node(anc2, [anc], branchstruct)
        assert mytree.nodedist(newnode2, anc2) == approx_places(0.999)
        assert mytree.getlabel(newnode2, anc2) == "Lower_branch"
        assert mytree.parent(newnode2) == anc2
        assert mytree.children(newnode2) == {anc}
        assert mytree.children(anc2) == {newnode2, "HV2SB"}

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"s1", "s2"})
        branchstruct = pt.Branchstruct(length=0.777)
        branchstruct.label = "New_branch"
        newnode = mytree.insert_node(anc, ["s1", "s2"], branchstruct)
        assert mytree.nodedist(newnode, anc) == approx_places(0.777)
        assert mytree.getlabel(newnode, anc) == "New_branch"
        assert mytree.parent(newnode) == anc
        assert mytree.children(newnode) == {"s1", "s2"}
        assert mytree.children(anc) == {newnode}

    def test_addbranch(self):
        """Does add_branch() work correctly?"""
        mytree = pt.Tree.from_leaves(["A", "B", "C", "D", "E"])
        bipart1 = frozenset([frozenset(["A", "B"]), frozenset(["C", "D", "E"])])
        bipart2 = frozenset([frozenset(["A", "B", "C"]), frozenset(["D", "E"])])
        pre_nodes = copy.copy(mytree.nodes)

        branchstruct = pt.Branchstruct(length=0.222)
        branchstruct.label = "First_branch"
        mytree.add_branch(bipart1, branchstruct)
        anc1 = mytree.parent("A")
        anc2 = mytree.parent("C")
        root = mytree.root
        assert mytree.nodedist(anc1, anc2) == approx_places(0.222)
        assert mytree.getlabel(root, anc1) == "First_branch"
        assert mytree.getlabel(root, anc2) == "First_branch"
        assert mytree.children(anc1) == {"A", "B"}
        assert mytree.children(anc2) == {"C", "D", "E"}
        assert len(pre_nodes) + 2 == len(mytree.nodes)
        assert pre_nodes < mytree.nodes

        branchstruct = pt.Branchstruct(length=0.333)
        branchstruct.label = "Second_branch"
        mytree.add_branch(bipart2, branchstruct)
        anc3 = mytree.parent("D")
        anc4 = mytree.parent(anc3)
        assert mytree.nodedist(anc3, anc4) == approx_places(0.333)
        assert mytree.getlabel(anc3, anc4) == "Second_branch"
        assert mytree.children(anc3) == {"D", "E"}
        assert mytree.children(anc4) == {anc3, "C"}

    def test_removebranch(self):
        """Does remove_branch work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        node1 = mytree.find_mrca({"HV1EL", "HV1Z2", "HV1Z8"})
        node2 = mytree.parent(node1)
        node2remotekids = mytree.remote_children(node2)
        pre_len = mytree.length()
        pre_nodes = copy.copy(mytree.nodes)

        mytree.remove_branch(node1, node2)
        assert mytree.remote_children(node2) == node2remotekids
        with pytest.raises(pt.TreeError):
            mytree.children(node1)
        assert pre_nodes - {node1} == mytree.nodes
        assert len(mytree.children(node2)) == 3
        assert pre_len == approx_places(mytree.length())

    def test_removeleaf(self):
        """Does remove_leaf() work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc1 = mytree.find_mrca({"SIVMK", "SIVML"})
        anc2 = mytree.parent(anc1)
        removed_dist = mytree.nodedist(anc1, "SIVMK")
        pre_len = mytree.length()
        pre_nodes = copy.copy(mytree.nodes)

        mytree.remove_leaf("SIVMK")
        assert mytree.children(anc2) == {"SIVML", "SIVM1"}
        assert mytree.length() == approx_places(pre_len - removed_dist)
        with pytest.raises(pt.TreeError):
            mytree.children(anc1)
        assert len(pre_nodes) - 2 == len(mytree.nodes)
        assert pre_nodes > mytree.nodes

        mytree = pt.Tree.from_string("(A, ((B, C), (D, E)));")
        anc1 = mytree.parent("A")
        anc2 = mytree.find_mrca({"B", "C", "D", "E"})
        pre_root = mytree.root

        mytree.remove_leaf("A")
        with pytest.raises(pt.TreeError):
            mytree.children(anc1)
        assert pre_root == anc1
        assert mytree.root == anc2

    def test_rename_leaf(self):
        """Does rename_leaf() work corectly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"SIVMK", "SIVML"})
        mytree.rename_leaf("SIVMK", "Monkey_Seq")
        assert mytree.children(anc) == {"SIVML", "Monkey_Seq"}
        assert mytree.parent("Monkey_Seq") == anc
        assert "SIVMK" not in mytree.remote_children(mytree.root)
        assert "SIVMK" not in mytree.nodes
        assert "SIVMK" not in mytree.leaves
        assert "Monkey_Seq" in mytree.nodes
        assert "Monkey_Seq" in mytree.leaves

        with pytest.raises(pt.TreeError):
            mytree.rename_leaf("Not_there", "new_name")


class TestRootTester:
    """Tests for rooting related methods in Tree object"""

    def test_rootmid(self):
        """Test that midpoint rooting creates expected tree on some known examples"""
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)

        mytree.rootmid()
        assert mytree.nodedist("s1") == approx_places(0.3125), f"leaves: {mytree.leaves}"
        assert mytree.nodedist("s2") == approx_places(0.3125)
        assert mytree.nodedist("s3") == approx_places(0.1875)
        assert mytree.nodedist("s4") == approx_places(0.3125)
        assert mytree.nodedist("s5") == approx_places(0.3125)

        treestring = self.treedata["HIVtree_HV1A2root"]
        mytree = pt.Tree.from_string(treestring)
        mytree.rootmid()
        assert mytree.nodedist("SIVM1") == approx_places(mytree.nodedist("HV1RH"))
        assert mytree.nodedist("SIVM1") == approx_places(0.34765)

    def test_rootminvar(self):
        """Test that rootminvar finds the minimum variance root"""
        ts1 = self.treedata["HIVtree"]
        t1 = pt.Tree.from_string(ts1)

        ts2 = self.treedata["HIVminvar"]
        t2 = pt.Tree.from_string(ts2)

        t1.rootminvar()
        assert t1 == t2


class TestDist_tree_construction:
    """Tests methods for constructing trees from distance matrix"""

    @pytest.fixture(autouse=True)
    def _setup(self, data_dir):
        self.njtreestringlist = [
            "((A:13.0,B:4.0):2.0,(C:4.0,D:10.0):2.0);",
            "((A:2,(B:1.5,C:1):3.2):4.1,(D:1.2,(E:2.7,F:3.1):5):5.2);",
        ]
        self.njdmatlist = [
            {
                "A": {"B": 17, "C": 21, "D": 27},
                "B": {"A": 17, "C": 12, "D": 18},
                "C": {"A": 21, "B": 12, "D": 14},
                "D": {"A": 27, "B": 18, "C": 14},
            },
            {
                "A": {"B": 6.7, "C": 6.2, "D": 12.5, "E": 19, "F": 19.4},
                "B": {"A": 6.7, "C": 2.5, "D": 15.2, "E": 21.7, "F": 22.1},
                "C": {"A": 6.2, "B": 2.5, "D": 14.7, "E": 21.2, "F": 21.6},
                "D": {"A": 12.5, "B": 15.2, "C": 14.7, "E": 8.9, "F": 9.3},
                "E": {"A": 19, "B": 21.7, "C": 21.2, "D": 8.9, "F": 5.8},
                "F": {"A": 19.4, "B": 22.1, "C": 21.6, "D": 9.3, "E": 5.8},
            },
        ]
        self.upgmatreelist = [
            "(((A:8.500,B:8.500):2.500,E:11.000):5.500,(C:14.000,D:14.000):2.500);",
            "((A:1.500,B:1.500):3.750,(C:2.000,D:2.000):3.250);",
        ]
        self.upgmadmatlist = [
            {
                "A": {"B": 17, "C": 21, "D": 31, "E": 23},
                "B": {"A": 17, "C": 30, "D": 34, "E": 21},
                "C": {"A": 21, "B": 30, "D": 28, "E": 39},
                "D": {"A": 31, "B": 34, "C": 28, "E": 43},
                "E": {"A": 23, "B": 21, "C": 39, "D": 43},
            },
            {
                "A": {"B": 3, "C": 11, "D": 11},
                "B": {"A": 3, "C": 10, "D": 10},
                "C": {"A": 11, "B": 10, "D": 4},
                "D": {"A": 11, "B": 10, "C": 4},
            },
        ]

        dmat_fname = data_dir / "large_distmat.tsv"
        njtree_fname = data_dir / "large_njtree.txt"
        upgmatree_fname = data_dir / "large_upgmatree.txt"

        with dmat_fname.open(encoding="utf-8") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            dictlist = []
            for rowdict in reader:
                dictlist.append(rowdict)
            self.large_distdict = dict(zip(reader.fieldnames, dictlist))
        self.largephylip_njtree = njtree_fname.read_text(encoding="utf-8")
        self.largephylip_upgmatree = upgmatree_fname.read_text(encoding="utf-8")

    def test_nj(self):
        """Verify that nj method produces correct trees for list of distance matrices"""
        for i in range(len(self.njtreestringlist)):
            inputtree = pt.Tree.from_string(self.njtreestringlist[i])
            dmat = pt.Distmatrix.from_distdict(self.njdmatlist[i])
            njtree = dmat.nj()
            assert njtree == inputtree
        inputtree = pt.Tree.from_string(self.largephylip_njtree)
        dmat = pt.Distmatrix.from_distdict(self.large_distdict)
        njtree = dmat.nj()
        assert njtree == inputtree

    def test_upgma(self):
        """Verify that upgma method produces correct trees for list of distance matrices"""
        for i in range(len(self.upgmatreelist)):
            inputtree = pt.Tree.from_string(self.upgmatreelist[i])
            dmat = pt.Distmatrix.from_distdict(self.upgmadmatlist[i])
            upgma = dmat.upgma()
            assert upgma == inputtree, f"\n  i: {i}\n  upgma:\n{upgma}\n  input:\n{inputtree}"
        inputtree = pt.Tree.from_string(self.largephylip_upgmatree)
        dmat = pt.Distmatrix.from_distdict(self.large_distdict)
        upgmatree = dmat.upgma()
        assert upgmatree == inputtree


class TestTreeDistTester:
    """Tests methods for computing tree distance and similarity"""

    def test_treedist(self):
        ts = self.treedata["HIVtree"]
        t1 = pt.Tree.from_string(ts)
        t2 = pt.Tree.from_string(ts)

        subtreenode = t1.parent("HV2BE")
        t2.spr(subtreenode, "SIVM1")
        treedist = t1.treedist(t2, normalise=False)
        assert treedist == 10


# Combined from test_phylotreelib_pytestversion.py

class Test_remove_comments:

    def test_balanced(self):
        text = 'Hello [comment] World!\n[newline]'
        assert pt.remove_comments(text) == 'Hello  World!\n'

    def test_unbalanced(self):
        text = 'Hello [comment World!\n[newline'
        with pytest.raises(pt.TreeError):
            pt.remove_comments(text)

###################################################################################################
###################################################################################################

# Tests for Interner class

###################################################################################################
###################################################################################################

# How to do this? Probably need to instantiate relevant classes and add items to interner
# and then check that different items point to same object in memory, but how?

###################################################################################################
###################################################################################################

# Tests for Branchstruct class

###################################################################################################
###################################################################################################

class Test_init_Branchstruct:

    def test_attribute_assignment(self):
        mylength = 0.1543
        mylabel = "0.95"
        branch = pt.Branchstruct(length=mylength)
        branch.label=mylabel
        assert branch.length == mylength
        assert branch.label == mylabel

class Test_copy_Branchstruct:
    def test_copy(self):
        mylength = 0.1543
        mylabel = "0.9965"
        b1 = pt.Branchstruct(length=mylength)
        b1.label=mylabel
        b1.freq = 0.9876
        b1.posterior = 0.34984
        b2 = b1.copy()
        assert b1.length == b2.length
        assert b1.label == b2.label
        assert b1.freq == b2.freq
        assert b1.posterior == b2.posterior

###################################################################################################
###################################################################################################

# Tests for Topostruct class

###################################################################################################
###################################################################################################

class Test_create_Topostruct:

    def test_attributes(self):
        ts = pt.Topostruct()
        ts.tree = "(A,(B,(C,D)));"  # not actually a tree but fine for testing...
        ts.posterior = 0.97
        ts.n = 23
        # Note: main test here is to see if class accepts (only) attribues named in __slots__
        # Assertion is done implicitly: if exception is raised above then test fails
        with pytest.raises(AttributeError):
            ts.notanatt = "not supposed to be an attribute"

###################################################################################################
###################################################################################################

# Tests for Tree class

###################################################################################################
###################################################################################################

class Test_from_string:
    def test_parse_simplestring(self, treedata):
        treestring = treedata["simplestring"]
        assert isinstance(pt.Tree.from_string(treestring), pt.Tree)

    def test_parse_blanks(self, treedata):
        treestring = treedata["string_with_blanks"]
        assert isinstance(pt.Tree.from_string(treestring), pt.Tree)

    def test_parse_newline(self, treedata):
        treestring = treedata["string_with_newlines"]
        assert isinstance(pt.Tree.from_string(treestring), pt.Tree)

    def test_parse_label(self, treedata):
        treestring = treedata["string_with_label"]
        assert isinstance(pt.Tree.from_string(treestring), pt.Tree)

    def test_parse_complexstring(self, treedata):
        treestring = treedata["complexstring"]
        assert isinstance(pt.Tree.from_string(treestring), pt.Tree)

    def test_parse_spuriousnewlines(self, treedata):
        treestring = treedata["string_with_weird_newlines"]
        assert isinstance(pt.Tree.from_string(treestring), pt.Tree)

    def test_return_correct(self, treedata):
        treestring = treedata["string_with_label"]
        mytree = pt.Tree.from_string(treestring)
        expected_leaves = {"KL0F07689", "KW081_13", "SBC669_26", "YAL016W"}
        expected_intnodes = {0, 1}
        expected_nodes = expected_leaves | expected_intnodes

        assert mytree.leaves == expected_leaves
        assert mytree.intnodes == expected_intnodes
        assert mytree.nodes == expected_nodes
        assert mytree.children(0) == {1, "KL0F07689", "KW081_13"}
        assert mytree.children(1) == {"SBC669_26", "YAL016W"}

###################################################################################################

class Test_parent_dict:

    def test_parentdict_not_set_by_default(self, treedata):
        mytree = pt.Tree.from_string(treedata["simplestring"])
        assert mytree._parent_dict == None

    def test_build_parentdict(self, treedata):
        mytree = pt.Tree.from_string(treedata["simplestring"])
        mytree.build_parent_dict()
        expected_dict = {0:None, 1:0, "s4":0, "s5":0, 2:1, "s3":1, "s1":2, "s2":2}
        assert mytree._parent_dict == expected_dict

    def test_parentdict_property(self, treedata):
        mytree = pt.Tree.from_string(treedata["simplestring"])
        pd = mytree.parent_dict
        expected_dict = {0:None, 1:0, "s4":0, "s5":0, 2:1, "s3":1, "s1":2, "s2":2}
        assert pd == expected_dict

    # Python note: add tests for whether dict is cleared in relevant situations

###################################################################################################

class Test_from_biplist:
    """Tests from_biplist() constructor"""

    def test_frombip(self):
        biplist = {frozenset([frozenset(["A"]), frozenset(["B", "C", "D"])]): pt.Branchstruct(0.1),
                   frozenset([frozenset(["B"]), frozenset(["A", "C", "D"])]): pt.Branchstruct(0.2),
                   frozenset([frozenset(["C"]), frozenset(["B", "A", "D"])]): pt.Branchstruct(0.3),
                   frozenset([frozenset(["D"]), frozenset(["B", "C", "A"])]): pt.Branchstruct(0.4),
                   frozenset([frozenset(["A", "B"]), frozenset(["C", "D"])]): pt.Branchstruct(0.5)}

        mytree = pt.Tree.from_biplist(biplist)
        assert isinstance(mytree, pt.Tree)
        assert mytree.leaves == {"A", "B", "C", "D"}
        assert mytree.intnodes == {0, 1}
        assert abs(mytree.length() - 1.5) < 1e-9  # compare floating point values with a small tolerance
        assert mytree.parent("A") == mytree.parent("B")
        assert mytree.parent("C") == mytree.parent("D")

###################################################################################################

class Test_from_topology:
    """Tests from_topology() constructor"""

    def test_frombip(self):
        topology = {frozenset([frozenset(["A"]), frozenset(["B", "C", "D"])]),
                   frozenset([frozenset(["B"]), frozenset(["A", "C", "D"])]),
                   frozenset([frozenset(["C"]), frozenset(["B", "A", "D"])]),
                   frozenset([frozenset(["D"]), frozenset(["B", "C", "A"])]),
                   frozenset([frozenset(["A", "B"]), frozenset(["C", "D"])])}

        mytree = pt.Tree.from_topology(topology)
        assert isinstance(mytree, pt.Tree)
        assert mytree.leaves == {"A", "B", "C", "D"}
        assert mytree.intnodes == {0, 1}
        assert mytree.length() == 0.0
        assert mytree.parent("A") == mytree.parent("B")
        assert mytree.parent("C") == mytree.parent("D")

###################################################################################################

class Test_from_leaf:
    """Tests from_leaf() constructor"""

    def test_fromleaf(self):
        """Does Tree.from_leaf() work correctly?"""
        # Note: this constructor relies mostly on from_string, so testing is minimal...
        leaves = ["A", "B", "C", "D", "E"]
        mytree = pt.Tree.from_leaves(leaves)
        assert isinstance(mytree, pt.Tree)
        assert mytree.children(mytree.root) == {"A", "B", "C", "D", "E"}
        assert len(mytree.intnodes) == 1
        assert len(mytree.leaves) == 5

###################################################################################################

class Test_from_branchinfo:
    """Tests from_branchinfo() constructor"""

    def test_frombranchinfo(self, treedata):
        """Does Tree.from_branchinfo() work correctly?"""
        # Testing is done by reading trees from treedata, converting them to lists,
        # recreating tree from lists, and comparing to original tree
        for treestring in treedata.values():
            origtree = pt.Tree.from_string(treestring)
            parentlist = []
            childlist = []
            lenlist = []
            lablist = []
            for parent in origtree.intnodes:
                for child in origtree.children(parent):
                    parentlist.append(parent)
                    childlist.append(child)
                    lenlist.append(origtree.nodedist(parent, child))
                    lablist.append(origtree.get_branch_attribute(parent, child, "label", default=""))
            newtree = pt.Tree.from_branchinfo(parentlist, childlist, lenlist, label=lablist)
            assert origtree == newtree

###################################################################################################

class Test_randtree:

    def test_randtree_ntips(self):
        """Does Tree.randtree() work correctly when ntips is provided?"""
        ntips = 5
        mytree = pt.Tree.randtree(ntips=ntips)

        assert isinstance(mytree, pt.Tree)
        assert len(mytree.leaves) == ntips
        assert mytree.is_resolved()
        for parent in mytree.intnodes:
            for child in mytree.child_dict[parent]:
                length = mytree.child_dict[parent][child].length
                assert length >= 0   # since lognormvariate will always produce positive values
                assert isinstance(length, float)

###################################################################################################

class Test_TreeIteration:
    """Tests iteration over file with several trees"""

    def test_iterate_newick(self, treedata, tmp_path):
        """Does iteration over tree objects from treefile work? """

        # First: construct treefile
        filename = tmp_path / "tempfile"
        trees = pt.TreeSet()
        for treestring in treedata.values():
            mytree = pt.Tree.from_string(treestring)
            trees.addtree(mytree)
        with open(filename, 'w') as file:
            file.write(trees.newick())

        # Secondly: iterate over treefile
        treefile = pt.Newicktreefile(filename)
        for tree in treefile:
            assert isinstance(tree, pt.Tree)

    def test_read_correctly_newick(self, treedata, tmp_path):
        """Do I also get the correctly read trees from treefile"""

        # First: construct treefile
        filename = tmp_path / "tempfile"
        treelist = []
        trees = pt.TreeSet()
        for treestring in treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            trees.addtree(mytree)
        with open(filename, 'w') as file:
            file.write(trees.newick())

        # Secondly: iterate over treefile, check that read trees correspond to written trees
        treefile = pt.Newicktreefile(filename)
        for i, tree in enumerate(treefile):
            assert tree == treelist[i]

    def test_iterate_nexus(self, treedata, tmp_path):
        """Does iteration over tree objects from nexus treefile work? """

        # First: construct treefile
        filename = tmp_path / "tempfile"
        with open(filename, 'w') as file:
            for treestring in treedata.values():
                mytree = pt.Tree.from_string(treestring)
                file.write(mytree.nexus())

        # Secondly: iterate over treefile
        treefile = pt.Nexustreefile(filename)
        for tree in treefile:
            assert isinstance(tree, pt.Tree)

    def test_read_correctly_nexus(self, treedata, tmp_path):
        """Do I also get the correctly read trees from nexus treefile"""

        # First: construct treefile
        filename = tmp_path / "tempfile"
        treelist = []
        with open(filename, 'w') as file:
            for treestring in treedata.values():
                mytree = pt.Tree.from_string(treestring)
                treelist.append(mytree)
                file.write(mytree.nexus())

        # Secondly: iterate over treefile, check that read trees correspond to written trees
        treefile = pt.Nexustreefile(filename)
        for i, tree in enumerate(treefile):
            assert tree == treelist[i]

###################################################################################################

class Test_str:

    def test_tree_str_representation(treedata):
        """Test the string representation of a Tree object"""

        parentlist = [0,0,0,1,1]
        childlist = ["A", "B", 1, "C", "D"]
        lenlist = [2.0, 2.0, 1.0, 1.0, 1.0]
        lablist = ["", "", "0.95", "", ""]
        my_tree = pt.Tree.from_branchinfo(parentlist, childlist, lenlist, label=lablist)

        # Expected representation
        expected_str_lines =   ["+--------------------------------------------+",
                                "|  parent  |  child  |  branchlen  |  label  |",
                                "+--------------------------------------------+",
                                "|       0  |      1  |          1  |   0.95  |",
                                "|       0  |      A  |          2  |         |",
                                "|       0  |      B  |          2  |         |",
                                "|       1  |      C  |          1  |         |",
                                "|       1  |      D  |          1  |         |",
                                "+--------------------------------------------+",
                                "Tree length: 7",
                                "",
                                "4 Leaves:",
                                "-",
                                "A",
                                "B",
                                "C",
                                "D"]
        tree_str_lines = str(my_tree).splitlines()
        # Note: intnodes are in sorted_intnode order, but order of kids not really defined??
        # therefore: test that all lines are the same, while ignoring order
        assert set(expected_str_lines) == set(tree_str_lines)

###################################################################################################

class Test_eq:
    """Tests special method __eq__() for Tree objects"""

    def test_equality(self, treedata):
        """Sanity check: can Tree.__eq__() differentiate between set of known trees?"""
        for treestring in treedata.values():
            mytree1 = pt.Tree.from_string(treestring)
            mytree2 = pt.Tree.from_string(treestring)
            assert mytree1 == mytree2

        stringlist = [
            treedata["simplestring"],
            treedata["complexstring"],
            treedata["HIVtree"]
        ]
        for i in range(len(stringlist) - 1):
            for j in range(i+1, len(stringlist)):
                mytree1 = pt.Tree.from_string(stringlist[i])
                mytree2 = pt.Tree.from_string(stringlist[j])
                assert mytree1 != mytree2

    def test_diffblens(self, treedata):
        """Check that different branch lengths results in non-equality"""
        for treestring in treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            for parent in t2.intnodes:
                for kid in t2.children(parent):
                    origblen = t2.nodedist(parent,kid)
                    t2.setlength(parent,kid,origblen*1.5)
            assert t1 != t2

###################################################################################################

class Test_hash:

    def test_same_object_same_hash(self):
        """Check that the same object produces the same hash value every time it's hashed."""
        tree = pt.Tree()
        assert hash(tree) == hash(tree)

    def test_different_objects_different_hashes(self):
        """Check that different objects produce different hash values."""
        tree1 = pt.Tree.randtree(ntips=4)
        tree2 = pt.Tree.randtree(ntips=4)

        assert tree1 is not tree2  # They are different objects
        assert hash(tree1) != hash(tree2)

    def test_use_as_dict_key(self):
        """Check that Tree objects can be used as dictionary keys."""
        tree1 = pt.Tree.randtree(ntips=4)
        tree2 = pt.Tree.randtree(ntips=4)

        my_dict = {tree1: 'value1'}

        assert tree1 in my_dict
        assert tree2 not in my_dict

        my_dict[tree2] = 'value2'
        assert tree2 in my_dict

    # If you decide to have a more complicated hash method in the future, add tests for it.

###################################################################################################

class Test_copy_treeobject:
    """Tests for copy_treeobject function"""

    def test_blen_lab(self, treedata):
        for treestring in treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=True, copyattr=True)
            assert t1 == t2
            assert t1.root == t2.root
            assert t1.has_same_root(t2)
            for parent in t1.intnodes:
                for kid in t1.children(parent):
                    t1lab = t1.getlabel(parent, kid)
                    t2lab = t2.getlabel(parent, kid)
                    assert t1lab == t2lab

    def test_blen_nolab(self, treedata):
        for treestring in treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=True, copyattr=False)
            assert t1 == t2
            assert t1.root == t2.root
            assert t1.has_same_root(t2)
            for parent in t2.intnodes:
                for kid in t2.children(parent):
                    t2lab = t2.getlabel(parent, kid)
                    assert t2lab == ""

    def test_noblen_nolab(self, treedata):
        for treestring in treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=False, copyattr=False)
            assert t1.topology() == t2.topology()
            assert t1.root == t2.root
            assert t1.has_same_root(t2)
            for parent in t2.intnodes:
                for kid in t2.children(parent):
                    t2lab = t2.getlabel(parent, kid)
                    assert t2lab == ""

###################################################################################################

class Test_build_dist_dict:

    def test_distdict(self, treedata):
        """Test that dist_dict is constructed correctly and agrees with nodedist function"""

        # First set of tests
        for treestring in treedata.values():
            mytree = pt.Tree.from_string(treestring)
            mytree.build_dist_dict()
            for n1 in mytree.nodes:
                for n2 in mytree.nodes:
                    assert pytest.approx(mytree.nodedist(n1, n2)) == mytree.dist_dict[n1][n2]

        # Specific node distance checks
        treestring = treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        mytree.build_dist_dict()
        assert pytest.approx(mytree.dist_dict["s4"]["s1"]) == 0.625
        assert pytest.approx(mytree.dist_dict["s2"]["s1"]) == 0.25
        assert pytest.approx(mytree.dist_dict["s3"][0]) == 0.25
        assert pytest.approx(mytree.dist_dict[0]["s1"]) == 0.5

        # Check after changing root
        for treestring in treedata.values():
            mytree = pt.Tree.from_string(treestring)
            mytree.rootmid()
            mytree.build_dist_dict()
            for n1 in mytree.intnodes:
                for n2 in mytree.leaves:
                    assert pytest.approx(mytree.nodedist(n1, n2)) == mytree.dist_dict[n1][n2]

        # Specific node distance checks after changing root
        treestring = treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        mytree.rootmid()
        mytree.build_dist_dict()
        assert pytest.approx(mytree.dist_dict["s4"]["s1"]) == 0.625
        assert pytest.approx(mytree.dist_dict["s2"]["s1"]) == 0.25
        assert pytest.approx(mytree.dist_dict["s3"][0]) == 0.25
        assert pytest.approx(mytree.dist_dict[0]["s1"]) == 0.5

###################################################################################################

class Test_build_path_dict:

    def test_pathdict(self, treedata):
        """Test that path_dict is constructed correctly and agrees with nodepath function"""

        # First set of tests
        for treestring in treedata.values():
            mytree = pt.Tree.from_string(treestring)
            mytree.build_path_dict()
            for n1 in mytree.intnodes:
                for n2 in mytree.leaves:
                    npath = mytree.nodepath(n1, n2)
                    npathnew = [n1]
                    n = n1
                    while n != n2:
                        n = mytree.path_dict[n][n2]
                        npathnew.append(n)
                    assert npath == npathnew

        # Specific node path checks
        treestring = treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        mytree.build_path_dict()
        assert mytree.path_dict["s4"]["s1"] == 0
        assert mytree.path_dict[0]["s1"] == 1
        assert mytree.path_dict[1]["s1"] == 2
        assert mytree.path_dict[2]["s1"] == "s1"

###################################################################################################

class Test_sorted_intnodes:

    def test_basic_functionality(self, treedata):
        for treestring in treedata.values():
            tree = pt.Tree.from_string(treestring)
            sorted_nodes = tree.sorted_intnodes()

            # All nodes in the sorted list should be internal nodes
            assert set(sorted_nodes) == tree.intnodes

            # There should be no duplicates
            assert len(sorted_nodes) == len(set(sorted_nodes))

    def test_ordering(self, treedata):
        for treestring in treedata.values():
            tree = pt.Tree.from_string(treestring)

            sorted_nodes_deepfirst = tree.sorted_intnodes(deepfirst=True)
            sorted_nodes_shallowfirst = tree.sorted_intnodes(deepfirst=False)

            # The reversed list of deepfirst should be equal to shallowfirst
            assert tuple(reversed(sorted_nodes_deepfirst)) == sorted_nodes_shallowfirst

    def test_tree_consistency(self, treedata):

        # children of node should appear before it in list when deepfirst=False
        # and after it when deepfirst=True.

        for treestring in treedata.values():
            tree = pt.Tree.from_string(treestring)

            # Check for deepfirst=True
            sorted_nodes = tree.sorted_intnodes(deepfirst=True)
            for i, node in enumerate(sorted_nodes):
                children = set(tree.children(node)) & tree.intnodes
                for child in children:
                    assert sorted_nodes.index(child) > i

            # Check for deepfirst=False
            sorted_nodes = tree.sorted_intnodes(deepfirst=False)
            for i, node in enumerate(sorted_nodes):
                children = set(tree.children(node)) & tree.intnodes
                for child in children:
                    assert sorted_nodes.index(child) < i

    def test_simplestring_sorted_intnodes(self, treedata):
        t = pt.Tree.from_string(treedata["simplestring"])
        result = t.sorted_intnodes(deepfirst=True)
        expected_result = (0,1,2)

        # Assert that the result matches the expected result
        assert result == expected_result

###################################################################################################

class Test_is_bifurcation:

    def test_bifurcations(self):
        """Check returns True: internal nodes with two descendants"""
        for i in range(10):
            t = pt.Tree.randtree(ntips=50)
            for intnode in t.intnodes:
                assert t.is_bifurcation(intnode)

    def test_trifurcations(self):
        """Check returns False: internal nodes with three descendants"""
        for i in range(10):
            t = pt.Tree.randtree(ntips=50)
            for intnode in t.intnodes:
                kid = next(iter(t.children(intnode)))
                if kid not in t.leaves:
                    t.remove_branch(intnode, kid)
                    assert not t.is_bifurcation(intnode)
                    break

###################################################################################################

class Test_n_bipartitions:

    def test_bifurcating_trees(self):
        for ntips in range(10,20):
            t = pt.Tree.randtree(ntips=ntips)
            expected_nbip = ntips - 3
            assert t.n_bipartitions() == expected_nbip

    def test_nonbifurcating_trees(self):
        for ntips in range(10,20):
            expected_nbip = ntips - 3
            t = pt.Tree.randtree(ntips=ntips)
            c = random.choice(list(t.intnodes - {t.root}))
            p = t.parent(c)
            if (p == t.root) and (t.is_bifurcation(p)):
                expected_nbip = expected_nbip
            else:
                expected_nbip = expected_nbip -1
            t.remove_branch(p,c)
            assert t.n_bipartitions() == expected_nbip

###################################################################################################

class Test_leaflist:

    def test_leaflist_sorted(self):
        chars = ascii_lowercase + digits
        namelist = [''.join(random.choices(chars, k=6)) for _ in range(25)]
        t = pt.Tree.from_leaves(leaflist=namelist)
        leafnames = t.leaflist()
        assert leafnames == sorted(namelist)

###################################################################################################

class Test_transdict:
    def test_transdict_correct_keys(self):
        ntips = random.randint(10,50)
        t = pt.Tree.randtree(ntips=ntips)
        trans_dictionary = t.transdict()
        assert set(trans_dictionary.keys()) == set(t.leaflist()), "Transdict keys do not match leaf names."

    def test_transdict_values_are_sequential(self):
        ntips = random.randint(10,50)
        t = pt.Tree.randtree(ntips=ntips)
        trans_dictionary = t.transdict()
        expected_values = [str(i) for i in range(1, len(trans_dictionary) + 1)]
        assert set(trans_dictionary.values()) == set(expected_values), "Transdict values are not sequential."

    def test_transdict_correct_key_value_pairs(self):
        ntips = random.randint(10,50)
        t = pt.Tree.randtree(ntips=ntips)
        trans_dictionary = t.transdict()
        sorted_leaves = t.leaflist()
        for idx, leaf in enumerate(sorted_leaves):
            assert trans_dictionary[leaf] == str(idx + 1), f"Leaf: {leaf} does not have expected value."

###################################################################################################

class Test_translateblock:

    def test_translateblock(self, treedata):
        for treestring in treedata.values():
            t = pt.Tree.from_string(treestring)
            trans_dictionary = t.transdict()
            result = t.translateblock(trans_dictionary)
            assert result.startswith("    translate\n"), "Output does not start with '    translate\n'"
            assert result.endswith("    ;\n"), "Output does not end with '    ;\n'"
            lines = result.split("\n")
            for line in lines[1:-3]:  # exclude the start and end lines
                assert line.strip().endswith(','), "Intermediate lines should end with a comma"
                number, name_with_comma = line.strip().split()
                assert number.isdigit(), "Number part should only contain digits"
                assert name_with_comma.endswith(",")

    def test_translateblock_key_value_formatting(self, treedata):
        for treestring in treedata.values():
            t = pt.Tree.from_string(treestring)
            trans_dictionary = t.transdict()
            result = t.translateblock(trans_dictionary)
            lines = result.split("\n")[1:-2]  # exclude the start and end lines
            nameset = set()
            for line in lines:
                number, name_with_comma = line.strip().split()
                name = name_with_comma.replace(",", "")
                nameset.add(name)
                assert name in t.leaves
                assert trans_dictionary[name] == number
            assert nameset == t.leaves

###################################################################################################

class TestSPR:
    def test_spr_raises_on_two_leaf_tree(self):
        """SPR should not work on a tree with only 2 leaves."""
        tree = pt.Tree.from_leaves(["A", "B"])
        with pytest.raises(pt.TreeError):
            tree.spr()

    def test_random_spr_no_params(self):
        """Without parameters, spr() should perform a random SPR while preserving leaf set.
        No exceptions should be raised"""
        for ntips in range(3,9):
            tree = pt.Tree.randtree(ntips=ntips)
            original_leaves = tree.leaves.copy()
            tree.spr()
            assert tree.leaves == original_leaves

    def test_random_spr_all_possible_params(self):
        """Test all possible combinations of prune and regraft nodes for range of random trees
        of different size. No exceptions should be raised"""
        for ntips in range(3,9):
            origtree = pt.Tree.randtree(ntips=ntips)
            original_leaves = origtree.leaves.copy()
            for prune_node in origtree.possible_spr_prune_nodes():
                for regraft_node in origtree.possible_spr_regraft_nodes(prune_node):
                    t = origtree.copy_treeobject()
                    t.spr(prune_node, regraft_node)
                    assert t.leaves == original_leaves

    def test_spr_with_prune_only(self):
        """When only a prune_node is given, spr() should choose a valid regraft node."""
        tree = pt.Tree.randtree(ntips=7)
        possible_prune_nodes = tree.possible_spr_prune_nodes()
        prune_node = next(iter(possible_prune_nodes))
        original_leaves = tree.leaves.copy()
        tree.spr(prune_node=prune_node)
        assert tree.leaves == original_leaves

    def test_spr_with_invalid_regraft_node(self):
        """Supplying a regraft_node that is not allowed should raise an error."""
        tree = pt.Tree.randtree(ntips=8)
        possible_prune_nodes = tree.possible_spr_prune_nodes()
        prune_node = next(iter(possible_prune_nodes))
        # Choose an invalid regraft node.
        possible_regraft_nodes = tree.possible_spr_regraft_nodes(prune_node)
        impossible_regraft_nodes = tree.nodes - possible_regraft_nodes
        invalid_regraft_node = next(iter(impossible_regraft_nodes))
        with pytest.raises(pt.TreeError):
            tree.spr(prune_node=prune_node, regraft_node=invalid_regraft_node)

    def test_spr_preserves_total_length_general_case(self):
        random.seed(0)
        tree = pt.Tree.randtree(ntips=15)
        base_len = tree.length()
        for _ in range(50):
            tree.spr()
            assert tree.length() == pytest.approx(base_len)

    def test_spr_rejects_tree_with_two_leaves(self):
        tree = pt.Tree.from_string("(A:1,B:1);")
        with pytest.raises(pt.TreeError):
            tree.spr()

    def test_spr_rejects_regraft_at_root(self):
        tree = pt.Tree.randtree(ntips=10)
        prune_node = random.choice(tuple(tree.possible_spr_prune_nodes()))
        with pytest.raises(pt.TreeError):
            tree.spr(prune_node=prune_node, regraft_node=tree.root)


###################################################################################################

class TestGraft:
    def test_graft_leaf_on_edge_preserves_edge_length_and_adds_connecting_branch(self):
        tree = pt.Tree.from_string("(A:1,B:2,C:3);")
        total_len_before = tree.length()

        parent = tree.root
        child = "A"
        connect_len = 0.7

        graftpoint = tree.graft(
            other="X",
            parent=parent,
            child=child,
            frac_from_parent=0.4,
            connect_length=connect_len,
        )

        # new leaf exists
        assert "X" in tree.leaves
        assert tree.is_parent_child_pair(graftpoint, "X")
        assert tree.child_dict[graftpoint]["X"].length == pytest.approx(connect_len)

        # edge (parent->child) length is preserved by splitting
        new_child_path_len = tree.child_dict[parent][graftpoint].length + tree.child_dict[graftpoint][child].length
        assert new_child_path_len == pytest.approx(1.0)

        # total length increased by connect_len (since we attached new leaf)
        assert tree.length() == pytest.approx(total_len_before + connect_len)

    def test_graft_tree_renames_internal_nodes_to_avoid_collision(self):
        t1 = pt.Tree.from_string("(A:1,(B:1,C:1):1);")
        t2 = pt.Tree.from_string("(D:1,(E:1,F:1):1);")
        # Force an ID collision by renaming t2 internal node(s) to overlap t1
        # (pick a non-root internal node in t2)
        t1_int = sorted([n for n in t1.intnodes if isinstance(n, int)])
        t2_int = sorted([n for n in t2.intnodes if isinstance(n, int)])
        assert t1_int and t2_int
        collide = t1_int[0]
        if collide in t2_int:
            # already colliding; keep it
            pass
        else:
            t2.rename_intnode(t2_int[0], collide)

        parent = t1.root
        child = "A"
        graftpoint = t1.graft(t2, parent=parent, child=child, frac_from_parent=0.5, connect_length=0.0)

        # All leaves from both trees present
        assert set(t1.leaves) == {"A", "B", "C", "D", "E", "F"}
        # Internal IDs should be unique
        ints = [n for n in t1.intnodes if isinstance(n, int)]
        assert len(ints) == len(set(ints))
        assert graftpoint in t1.intnodes

###################################################################################################

class TestPrune:
    """Tests for pruning related methods in Tree object"""

    def test_prune_maxlen(self):
        """Brute force: prune_maxlen finds the maximum-length subtree with nkeep leaves."""
        random.seed(0)  # keep this test deterministic
        for _ in range(100):
            ntips = random.randint(5, 11)          # combinatorial explosion beyond this
            nkeep = random.randint(3, ntips - 1)

            t1 = pt.Tree.randtree(ntips=ntips)

            lengths = []
            for discardset in itertools.combinations(t1.leaves, ntips - nkeep):
                t2 = t1.copy_treeobject()
                t2.remove_leaves(discardset)
                lengths.append(t2.length())

            t1.prune_maxlen(nkeep=nkeep)
            assert t1.length() == pytest.approx(max(lengths))

    def test_prune_maxlen_with_keeplist(self):
        """Brute force: prune_maxlen respects keeplist and still maximizes length."""
        random.seed(1)  # keep deterministic, but different sequence than test above
        for _ in range(100):
            ntips = random.randint(5, 11)
            nkeep = random.randint(3, ntips - 1)
            nkeeplist = random.randint(1, nkeep)

            t1 = pt.Tree.randtree(ntips=ntips)
            keeplist = random.sample(list(t1.leaves), nkeeplist)

            potential_to_remove = t1.leaves - set(keeplist)
            lengths = []
            for discardset in itertools.combinations(potential_to_remove, ntips - nkeep):
                t2 = t1.copy_treeobject()
                t2.remove_leaves(discardset)
                lengths.append(t2.length())

            t1.prune_maxlen(nkeep=nkeep, keeplist=keeplist)
            assert t1.length() == pytest.approx(max(lengths))

    def test_find_common_leaf(self, treedata):
        """find_common_leaf returns the leaf with minimal sum of distances to all leaves in the set."""
        for ts in treedata.values():
            t = pt.Tree.from_string(ts)
            leafset = t.leaves

            best_leaf = None
            best_sum = float("inf")
            for leaf1 in leafset:
                blensum = 0.0
                for leaf2 in leafset:
                    blensum += t.nodedist(leaf1, leaf2)
                if blensum < best_sum:
                    best_sum = blensum
                    best_leaf = leaf1

            common_function = t.find_common_leaf(leafset)
            assert common_function == best_leaf

    def test_find_central_leaf(self):
        """find_central_leaf returns the leaf closest to the tree center for this known tree."""
        t = pt.Tree.from_string("(((A:1,B:1):9,C:9):1,(D:1,E:1):9);")
        central = t.find_central_leaf(t.leaves)
        assert central == "C"

###################################################################################################

class TestSubtreeAndPruneSubtree:
    def test_subtree_leaf_returns_str_and_basalbranch(self):
        tree = pt.Tree.from_string("(A:1,(B:2,C:3):4);")
        # pick a leaf that is not root
        leaf = "A"
        sub, basal = tree.subtree(leaf)
        assert isinstance(sub, str)
        assert sub == leaf
        assert basal.length == pytest.approx(1.0)

    def test_subtree_internal_returns_tree_and_basalbranch(self):
        tree = pt.Tree.from_string("(A:1,(B:2,C:3):4);")
        # internal node is the one above B,C
        root = tree.root
        # find internal child of root
        internal = next(n for n in tree.children(root) if n in tree.intnodes)

        sub, basal = tree.subtree(internal)
        assert isinstance(sub, pt.Tree)
        assert sub.root == internal
        assert set(sub.leaves) == {"B", "C"}
        assert basal.length == pytest.approx(4.0)

    def test_prune_subtree_removes_leaves_and_returns_basalbranch(self):
        tree = pt.Tree.from_string("(A:1,(B:2,C:3):4,D:5);")
        root = tree.root
        internal = next(n for n in tree.children(root) if n in tree.intnodes)
        leaves_before = set(tree.leaves)
        total_len_before = tree.length()

        sub, basal = tree.prune_subtree(internal)
        assert isinstance(sub, pt.Tree)
        assert basal.length == pytest.approx(4.0)

        # B and C should be gone from the main tree
        assert set(tree.leaves) == leaves_before - {"B", "C"}
        # Tree length should decrease by (basal branch + internal subtree length)
        assert tree.length() == pytest.approx(total_len_before - (4.0 + 2.0 + 3.0))

###################################################################################################

class TestAddNodeOnBranch:
    def test_add_node_on_branch_splits_length_by_dist(self):
        tree = pt.Tree.from_string("(A:2,B:2);")
        # root has children A and B; pick edge root->A
        root = tree.root
        newnode = tree.add_node_on_branch(root, "A", dist_from_parent=0.5, copy_attrs="both")

        assert tree.is_parent_child_pair(root, newnode)
        assert tree.is_parent_child_pair(newnode, "A")

        upper = tree.child_dict[root][newnode].length
        lower = tree.child_dict[newnode]["A"].length
        assert upper == pytest.approx(0.5)
        assert lower == pytest.approx(1.5)

    def test_add_node_on_branch_splits_length_by_frac(self):
        tree = pt.Tree.from_string("(A:2,B:2);")
        root = tree.root
        newnode = tree.add_node_on_branch(root, "A", frac_from_parent=0.25, copy_attrs="both")

        upper = tree.child_dict[root][newnode].length
        lower = tree.child_dict[newnode]["A"].length
        assert upper == pytest.approx(0.5)
        assert lower == pytest.approx(1.5)

    def test_add_node_on_branch_random_when_no_dist_or_frac(self):
        random.seed(123)
        tree = pt.Tree.from_string("(A:2,B:2);")
        root = tree.root
        newnode = tree.add_node_on_branch(root, "A", dist_from_parent=None, frac_from_parent=None)

        upper = tree.child_dict[root][newnode].length
        lower = tree.child_dict[newnode]["A"].length
        assert 0.0 < upper < 2.0
        assert upper + lower == pytest.approx(2.0)

    def test_add_node_on_branch_rejects_non_edge(self):
        tree = pt.Tree.from_string("(A:1,(B:1,C:1):1);")
        # A is not parent of C
        with pytest.raises(pt.TreeError):
            tree.add_node_on_branch("A", "C", frac_from_parent=0.5)

###################################################################################################

class TestNextInternalNodeId:
    def test_next_internal_node_id_int_tree(self):
        tree = pt.Tree.from_string("(A:1,(B:1,C:1):1);")
        # internal nodes are integers; root is typically 0
        nid = tree.next_internal_node_id()
        assert isinstance(nid, int)
        assert nid not in tree.nodes

    def test_next_internal_node_id_string_intnodes_tree(self):
        # transmission-tree-like: internal nodes can be strings
        parentlist = ["r", "r", "x", "x"]
        childlist = ["x", "y", "A", "B"]
        lenlist = [1.0, 1.0, 2.0, 3.0]
        tree = pt.Tree.from_branchinfo(parentlist, childlist, lenlist)
        nid = tree.next_internal_node_id()
        assert isinstance(nid, int)
        assert nid not in tree.nodes

###################################################################################################

class Test_compute_sumtree:
    """Compare pre-computed summary trees (treeannotator) to own computations.
    Check various combinations of summary tree type, blen-setting, and rooting"""

    # Python note: surely this can be done with fewer lines of code, but how then
    # to get individual test PASSED messages?

    def test_hipstr_cladeheight(self, data_dir):
        """HIPSTR tree with mean node heights and rooting at best HIPSTR resolution of root clade"""
        beast_data_dir = data_dir / "beast"

        # Own computation
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.trees")
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="hip", blen="cladeheight")

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.treeannot_hipstr_mean")
        tgold = tf.readtree()

        # Comparisons
        assert town.topology_clade == tgold.topology_clade # To get separate error message
        assert town.equals(tgold, rooted=True)  # Checks rooted topology again, and blens

    def test_mrhipstr_cladeheight(self, data_dir):
        """mrHIPSTR tree with mean node heights and rooting at best HIPSTR resolution of root clade"""
        beast_data_dir = data_dir / "beast"

        # Own computation
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.trees")
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="mrhip", blen="cladeheight")

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.treeannot_mrhipstr_mean")
        tgold = tf.readtree()

        # Comparisons
        assert town.topology_clade == tgold.topology_clade # To get separate error message
        assert town.equals(tgold, rooted=True)  # Checks rooted topology again, and blens

    def test_mcc_cladeheight(self, data_dir):
        """MCC tree with mean node heights and rooting at best tree's original root"""
        beast_data_dir = data_dir / "beast"

        # Own computation
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.trees")
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, tracktopo=True,
                              track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="mcc", blen="cladeheight")

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.treeannot_mcc_mean")
        tgold = tf.readtree()

        # Comparisons
        assert town.topology_clade == tgold.topology_clade # To get separate error message
        assert town.equals(tgold, rooted=True)  # Checks rooted topology again, and blens

    def test_mcc_caheight(self, data_dir):
        """MCC tree with CA node heights and rooting at best tree's original root"""
        beast_data_dir = data_dir / "beast"

        # Own computation
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.trees")
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, tracktopo=True,
                              track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="mcc", blen="caheight",
                                    count_burnin_filename_list=[(1000, 0, (beast_data_dir / "random_10tips_1000.trees"))])

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile(beast_data_dir / "random_10tips_1000.treeannot_mcc_ca")
        tgold = tf.readtree()

        # Comparisons
        assert town.topology_clade == tgold.topology_clade # To get separate error message
        assert town.equals(tgold, rooted=True)  # Checks rooted topology again, and blens


###################################################################################################
###################################################################################################
# Tests for QuantileAccumulator
###################################################################################################
###################################################################################################

class Test_init_QuantileAccumulator:

    def test_defaults(self):
        qa = pt.QuantileAccumulator()
        assert qa.k == 7
        assert qa.shift == qa.k + 1
        assert qa.scale == (1 << qa.shift)
        assert qa.mask == qa.scale - 1
        assert isinstance(qa.counts, dict) or "defaultdict" in type(qa.counts).__name__.lower()
        assert qa.n == 0
        assert isinstance(qa.neg_bucket, int)

    def test_custom_k(self):
        qa = pt.QuantileAccumulator(k=3)
        assert qa.k == 3
        assert qa.shift == 4
        assert qa.scale == 16
        assert qa.mask == 15

###################################################################################################

class Test__bucket_QuantileAccumulator:

    def test_bucket_sentinel_for_nonpositive(self):
        qa = pt.QuantileAccumulator()
        assert qa._bucket(0.0) == qa.neg_bucket
        assert qa._bucket(-1.0) == qa.neg_bucket

    def test_bucket_sentinel_for_nonfinite(self):
        qa = pt.QuantileAccumulator()
        assert qa._bucket(float("nan")) == qa.neg_bucket
        assert qa._bucket(float("inf")) == qa.neg_bucket
        assert qa._bucket(float("-inf")) == qa.neg_bucket

    def test_bucket_key_encodes_exponent_and_mantissa_bucket(self):
        qa = pt.QuantileAccumulator(k=7)
        x = 3.141592653589793
        bkey = qa._bucket(x)
        assert bkey != qa.neg_bucket

        # decode
        e = bkey >> qa.shift
        b = bkey & qa.mask

        # re-derive via frexp
        m, e2 = math.frexp(x)
        b2 = int(m * qa.scale)

        assert e == e2
        assert b == b2

    def test_bucket_interval_contains_x(self):
        """
        For x>0 finite, bucket should correspond to interval:
            mant in [b/scale, (b+1)/scale)
            value in [ldexp(b/scale, e), ldexp((b+1)/scale, e))
        and x should be in that interval.
        """
        qa = pt.QuantileAccumulator(k=10)
        for x in [1e-12, 1e-6, 0.1, 1.0, 2.0, 10.0, 1e6, 1e12]:
            bkey = qa._bucket(x)
            e = bkey >> qa.shift
            b = bkey & qa.mask

            lo = math.ldexp(b / qa.scale, e)
            hi = math.ldexp((b + 1) / qa.scale, e)

            assert lo <= x
            assert x < hi


###################################################################################################

class Test__bucket_value_QuantileAccumulator:

    def test_bucket_value_for_sentinel_is_zero(self):
        qa = pt.QuantileAccumulator()
        assert qa._bucket_value(qa.neg_bucket) == 0.0

    def test_bucket_value_within_bucket_interval(self):
        qa = pt.QuantileAccumulator(k=8)
        xs = [0.1, 0.5, 0.999, 1.0, 1.2345, 2.0, 12345.6]
        for x in xs:
            bkey = qa._bucket(x)
            e = bkey >> qa.shift
            b = bkey & qa.mask

            lo = math.ldexp(b / qa.scale, e)
            hi = math.ldexp((b + 1) / qa.scale, e)

            v = qa._bucket_value(bkey)
            assert lo <= v
            assert v < hi

###################################################################################################

class Test_add_QuantileAccumulator:

    def test_add_increments_n_and_bucket_count(self):
        qa = pt.QuantileAccumulator(k=7)
        assert qa.n == 0
        x = 1.0
        bkey = qa._bucket(x)

        qa.add(x)
        assert qa.n == 1
        assert qa.counts[bkey] == 1

        qa.add(x)
        assert qa.n == 2
        assert qa.counts[bkey] == 2

    def test_add_sentinel_bucket_is_counted(self):
        qa = pt.QuantileAccumulator()
        qa.add(0.0)
        qa.add(-5.0)
        assert qa.n == 2
        assert qa.counts[qa.neg_bucket] == 2


###################################################################################################

class Test_merge_QuantileAccumulator:

    def test_merge_combines_counts_and_n(self):
        a = pt.QuantileAccumulator(k=7)
        b = pt.QuantileAccumulator(k=7)

        for x in [1.0, 2.0, 3.0]:
            a.add(x)
        for x in [2.0, 4.0]:
            b.add(x)

        a.merge(b)
        assert a.n == 5
        # Check some bucket counts explicitly
        assert a.counts[a._bucket(2.0)] == 2

    def test_merge_incompatible_shift_raises(self):
        a = pt.QuantileAccumulator(k=7)
        b = pt.QuantileAccumulator(k=6)  # different shift
        with pytest.raises(ValueError):
            a.merge(b)

    def test_merge_incompatible_neg_bucket_raises(self):
        a = pt.QuantileAccumulator(k=7)
        b = pt.QuantileAccumulator(k=7)
        b.neg_bucket = -123  # force incompatibility
        with pytest.raises(ValueError):
            a.merge(b)

    def test_merge_equivalent_to_single_stream(self):
        """
        If you split data across accumulators and merge, the resulting quantiles
        should match building one accumulator on all points (same bucket scheme).
        """
        data = [random.lognormvariate(0.0, 1.0) for _ in range(500)]
        probs = [0.0, 0.1, 0.5, 0.9, 1.0]

        all_in_one = pt.QuantileAccumulator(k=10)
        for x in data:
            all_in_one.add(x)

        left = pt.QuantileAccumulator(k=10)
        right = pt.QuantileAccumulator(k=10)
        for x in data[:250]:
            left.add(x)
        for x in data[250:]:
            right.add(x)

        left.merge(right)

        assert left.n == all_in_one.n
        assert left.quantiles(probs) == all_in_one.quantiles(probs)


###################################################################################################

class Test_quantiles_QuantileAccumulator:

    def test_quantiles_empty_raises(self):
        qa = pt.QuantileAccumulator()
        with pytest.raises(pt.TreeError):
            qa.quantiles([0.5])

    def test_quantiles_prob_out_of_range_raises(self):
        qa = pt.QuantileAccumulator()
        qa.add(1.0)
        with pytest.raises(pt.TreeError):
            qa.quantiles([-0.01])
        with pytest.raises(pt.TreeError):
            qa.quantiles([1.01])

    def test_quantiles_returns_same_length_and_order_as_input_probs(self):
        qa = pt.QuantileAccumulator(k=12)
        for x in [1.0, 2.0, 3.0, 4.0, 5.0]:
            qa.add(x)

        probs = [0.9, 0.1, 0.5, 0.0, 1.0]
        out = qa.quantiles(probs)

        assert len(out) == len(probs)

        # Same call but reordered probs; compare by matching probabilities
        # (not by assuming numeric values are "true quantiles")
        out2 = qa.quantiles(sorted(probs))
        mapping = dict(zip(sorted(probs), out2))
        assert out == [mapping[p] for p in probs]

    def test_quantiles_monotone_for_sorted_probs(self):
        qa = pt.QuantileAccumulator(k=10)
        data = [random.random() + 1e-9 for _ in range(200)]
        for x in data:
            qa.add(x)

        probs = [0.0, 0.25, 0.5, 0.75, 1.0]
        out = qa.quantiles(probs)
        assert out == sorted(out)

    def test_quantiles_bounds_reasonable(self):
        qa = pt.QuantileAccumulator(k=10)
        data = [random.lognormvariate(0.0, 1.0) for _ in range(500)]
        for x in data:
            qa.add(x)

        bkeys = sorted(qa.counts)
        vmin = qa._bucket_value(bkeys[0])
        vmax = qa._bucket_value(bkeys[-1])

        q0, q50, q1 = qa.quantiles([0.0, 0.5, 1.0])

        assert q0 == vmin
        assert q1 == vmax
        assert q0 <= q50 <= q1

    def test_quantiles_p0_returns_min_bucket_value(self):
        qa = pt.QuantileAccumulator(k=10)
        for x in [0.8, 1.2, 2.5]:
            qa.add(x)
        bmin = qa._bucket_value(min(qa.counts))
        assert qa.quantile(0.0) == bmin

###################################################################################################

class Test_quantile_QuantileAccumulator:

    def test_quantile_calls_quantiles_singleton(self):
        qa = pt.QuantileAccumulator(k=12)
        for x in [1.0, 10.0, 100.0]:
            qa.add(x)

        q1 = qa.quantile(0.5)
        q2 = qa.quantiles([0.5])[0]
        assert q1 == q2

    def test_quantile_empty_raises(self):
        qa = pt.QuantileAccumulator()
        with pytest.raises(pt.TreeError):
            qa.quantile(0.5)

    def test_quantile_prob_out_of_range_raises(self):
        qa = pt.QuantileAccumulator()
        qa.add(1.0)
        with pytest.raises(pt.TreeError):
            qa.quantile(-0.1)
        with pytest.raises(pt.TreeError):
            qa.quantile(1.1)

###################################################################################################

class Test_precision_quantiles_vs_exact_order_statistic:

    def ceil_rank(self, p, n):
        """Your convention: rank j in {1..n}."""
        if p == 0.0:
            return 1
        return int(math.ceil(p * n))


    def test_quantile_matches_bucket_of_target_order_statistic(self):
        """
        For your definition (rank j = ceil(p*n), p=0 -> j=1),
        qa.quantile(p) should equal the representative value of the bucket
        that contains x_(j).
        """
        rng = random.Random()
        data = [rng.lognormvariate(0.0, 1.0) for _ in range(4000)]  # all > 0
        data_sorted = sorted(data)
        probs = [0.0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1.0]

        for k in [4, 7, 10, 12]:
            qa = pt.QuantileAccumulator(k=k)
            for x in data:
                qa.add(x)

            for p in probs:
                j = self.ceil_rank(p, len(data_sorted))
                xj = data_sorted[j - 1]

                # The bucket that contains the target order statistic
                bkey = qa._bucket(xj)
                expected = qa._bucket_value(bkey)

                got = qa.quantile(p)
                assert got == expected, f"k={k}, p={p}, j={j}, xj={xj}, bkey={bkey}"

    def test_quantile_within_relative_error_bound_of_target_order_statistic(self):
        """
        With positive finite data, bucket midpoint representation guarantees
        relative error <= 2^(-(k+1)) with respect to any value in that bucket.
        In particular, the returned value should be within this distance from the target order statistic x_(j).
        """
        rng = random.Random()
        data = [rng.lognormvariate(0.0, 1.0) for _ in range(4000)]
        data_sorted = sorted(data)
        probs = [0.0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1.0]

        for k in [4, 7, 10, 12]:
            qa = pt.QuantileAccumulator(k=k)
            for x in data:
                qa.add(x)

            eps = 2 ** (-(k + 1))       # worst-case relative error bound
            tol = eps + 5e-15           # small slack for floating point

            for p in probs:
                j = self.ceil_rank(p, len(data_sorted))
                xj = data_sorted[j - 1]         # exact quantile, given my rank definition
                q_approx = qa.quantile(p)

                lo = (1.0 - tol) * xj
                hi = (1.0 + tol) * xj

                assert lo <= q_approx <= hi, (
                    f"k={k}, p={p}, j={j}, xj={xj}, q_approx={q_approx}, "
                    f"band=[{lo},{hi}]"
                )

    def test_error_bound_monotone_in_k_against_target_order_statistic(self):
        """
        Stronger + stable than comparing to numpy: compare to the target order statistic x_(j)
        under your definition. The worst-case bound shrinks as k increases.
        Observed relative error should typically not increase when k increases.
        """
        rng = random.Random()
        data = [rng.lognormvariate(0.0, 1.0) for _ in range(8000)]
        data_sorted = sorted(data)
        p = 0.99
        j = self.ceil_rank(p, len(data_sorted))
        xj = data_sorted[j - 1]

        def relerr_for_k(k):
            qa = pt.QuantileAccumulator(k=k)
            for x in data:
                qa.add(x)
            q_approx = qa.quantile(p)
            return abs(q_approx - xj) / xj

        err_k4  = relerr_for_k(4)
        err_k10 = relerr_for_k(10)
        err_k14 = relerr_for_k(14)

        assert err_k10 <= err_k4 + 1e-15
        assert err_k14 <= err_k10 + 1e-15


# Combined from test_sumt_biplen_reference.py

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

def test_biplen_means_match_reference_on_majority_splits(data_dir):
    """
    For a fixed tree sample file:
      - compute reference mean lengths for each split (conditional on presence)
      - compute sumt's summary tree (--con --biplen)
      - compare branch lengths split-by-split
    """
    f = data_dir / "mrbayes" / "neanderthal.nexus.run1.t"
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


def test_consensus_topology_matches_reference_majority_set(data_dir):
    """
    The consensus topology should contain exactly the set of majority internal splits.
    """
    f = data_dir / "mrbayes" / "neanderthal.nexus.run1.t"
    ref = reference_aggregate(f, burnin=0)
    sumtree, _ = build_sumt_sumtree(f, burnin=0)

    majority = {hm for hm, c in ref.count.items() if (c / ref.n_trees) >= 0.5}
    sumtree_splits = set(dict(iter_sumtree_internal_splits_with_lengths(sumtree)).keys())

    # In a majority-rule consensus, internal splits should be exactly the majority set.
    # (Leaves are always present but excluded here by construction.)
    assert sumtree_splits == majority


class TestParsimony:
    def _build_tree_with_states(self, treestring, leaf_states=None, internal_states=None):
        tree = pt.Tree.from_string(treestring)
        ndict = tree.nodedict

        if leaf_states:
            for leaf, state in leaf_states.items():
                ndict[leaf].state = state

        if internal_states:
            for node, state in internal_states.items():
                ndict[node].state = state

        return tree

    def _all_nodes_fit_within_optimal_set(self, fitted_tree, possible_tree):
        fitted_ndict = fitted_tree.nodedict
        possible_ndict = possible_tree.nodedict
        for node in fitted_tree.nodes:
            assert fitted_ndict[node].fit in possible_ndict[node].optimal_set

    def test_parsimony_score_all_same_state(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "X", "C": "X", "D": "X"},
        )

        possible = tree.parsimony_possible_states()
        pscore, countdict = tree.parsimony_count_changes(fitpref="X")
        fitted = tree.parsimony_assign_fits(fitpref="X")

        assert possible is not tree
        assert pscore == 0
        assert dict(countdict) == {}

        for node in possible.nodes:
            assert possible.nodedict[node].optimal_set == {"X"}
            assert fitted.nodedict[node].fit == "X"
            assert fitted.nodedict[node].wasambig is False

    def test_parsimony_score_clean_two_state_split(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "X", "C": "Y", "D": "Y"},
        )

        ab_node = tree.find_mrca({"A", "B"})
        cd_node = tree.find_mrca({"C", "D"})
        possible = tree.parsimony_possible_states()
        pscore, countdict = tree.parsimony_count_changes(fitpref="X")

        assert pscore == 1
        assert sum(countdict.values()) == pscore
        assert possible.nodedict[tree.root].optimal_set == {"X", "Y"}
        assert possible.nodedict[ab_node].optimal_set == {"X"}
        assert possible.nodedict[cd_node].optimal_set == {"Y"}
        assert set(countdict) <= {("X", "Y"), ("Y", "X")}

    def test_parsimony_score_two_states_not_matching_topology(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "Y", "C": "X", "D": "Y"},
        )

        pscore, countdict = tree.parsimony_count_changes(fitpref="X")

        assert pscore == 2
        assert sum(countdict.values()) == pscore

    def test_parsimony_score_three_states(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1,E:1);",
            leaf_states={"A": "X", "B": "Y", "C": "X", "D": "Z", "E": "X"},
        )

        pscore, countdict = tree.parsimony_count_changes(fitpref="X")

        assert pscore == 2
        assert sum(countdict.values()) == pscore

    def test_parsimony_score_multifurcation_star_tree(self):
        tree = self._build_tree_with_states(
            "(A:1,B:1,C:1,D:1);",
            leaf_states={"A": "X", "B": "X", "C": "Y", "D": "Y"},
        )

        possible = tree.parsimony_possible_states()
        pscore, countdict = tree.parsimony_count_changes(fitpref="X")

        assert possible.nodedict[tree.root].optimal_set == {"X", "Y"}
        assert pscore == 2
        assert sum(countdict.values()) == pscore

    def test_parsimony_internal_observed_state_controls_fit(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "X", "C": "Y", "D": "Y"},
        )
        ab_node = tree.find_mrca({"A", "B"})
        tree.nodedict[ab_node].state = "X"

        fitted = tree.parsimony_assign_fits(fitpref="Y")

        assert fitted.nodedict[ab_node].fit == "X"
        assert fitted.nodedict[ab_node].wasambig is False

    def test_parsimony_fitpref_steers_ambiguous_root(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "X", "C": "Y", "D": "Y"},
        )

        fitted_x = tree.parsimony_assign_fits(fitpref="X")
        fitted_y = tree.parsimony_assign_fits(fitpref="Y")
        pscore_x, countdict_x = tree.parsimony_count_changes(fitpref="X")
        pscore_y, countdict_y = tree.parsimony_count_changes(fitpref="Y")

        assert fitted_x.nodedict[fitted_x.root].fit == "X"
        assert fitted_y.nodedict[fitted_y.root].fit == "Y"
        assert pscore_x == 1
        assert pscore_y == 1
        assert sum(countdict_x.values()) == pscore_x
        assert sum(countdict_y.values()) == pscore_y

    def test_parsimony_leaf_fits_always_equal_leaf_state(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "X", "C": "Y", "D": "Y"},
        )

        fitted = tree.parsimony_assign_fits(fitpref="X")

        for leaf in fitted.leaves:
            assert fitted.nodedict[leaf].fit == fitted.nodedict[leaf].state

    def test_parsimony_all_fits_within_optimal_sets(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "X", "C": "Y", "D": "Y"},
        )

        possible = tree.parsimony_possible_states()
        fitted = tree.parsimony_assign_fits(fitpref="Y")

        self._all_nodes_fit_within_optimal_set(fitted, possible)

    def test_parsimony_possible_states_returns_copy(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": "X", "C": "X", "D": "X"},
        )

        possible = tree.parsimony_possible_states()

        assert possible is not tree

    def test_parsimony_possible_states_errors_without_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")

        with pytest.raises(pt.TreeError):
            tree.parsimony_possible_states()

    def test_parsimony_possible_states_errors_without_any_states(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        _ = tree.nodedict

        with pytest.raises(pt.TreeError):
            tree.parsimony_possible_states()

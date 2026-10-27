import copy
import csv
import itertools
import math
import pickle
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


@pytest.fixture
def annotated_print_tree():
    tree = pt.Tree.from_string("((A:0.123456789,B:2)old:3,C:4);")
    tree.set_node_attribute("A", "height", 1.23456789)
    tree.set_node_attribute("A", "height_sd", 0.1)
    tree.set_node_attribute("A", "extra_meta", {"height_80%_CI": (1.0, 1.4), "height_95%_CI": (0.8, 1.5)})
    tree.set_node_attribute(tree.root, "height", 5.0)
    tree.set_node_attribute(tree.root, "extra_meta", {"height_80%_CI": (4.5, 5.5), "height_95%_CI": (4.0, 6.0)})
    node = tree.find_mrca({"A", "B"})
    branch = tree.child_dict[tree.root][node]
    branch.posterior = 0.975
    branch.length_sd = 0.25
    branch.extra_meta = {"length_80%_CI": (2.5, 3.5), "length_95%_CI": (2.0, 4.0)}
    return tree


class TestTreeNewick:

    def test_defaults_and_call_independence(self, annotated_print_tree):
        tree = annotated_print_tree
        expected = "((A:0.123456789,B:2.0)old:3.0,C:4.0);"
        assert tree.newick() == expected
        tree.newick(printdist=False, printlabels=False, print_meta=True,
                    precision=3, node_attributes=("height", "height_95%_CI"), branch_attributes=(), )
        tree.nexus(print_meta=True, node_attributes=("height",))
        assert tree.newick() == expected

    @pytest.mark.parametrize("options, expected", [
        ({"printdist": False}, "((A,B)old,C);"),
        ({"printlabels": False}, "((A:0.1234568,B:2):3,C:4);"),
        ({"labelfield": "posterior", "precision": 3}, "((A:0.123,B:2)0.975:3,C:4);"),
        ({"transdict": {"A": "1", "B": "2", "C": "3"}},
         "((1:0.1234568,2:2)old:3,3:4);"),
    ])
    def test_explicit_options(self, annotated_print_tree, options, expected):
        options.setdefault("precision", 7)
        assert annotated_print_tree.newick(**options) == expected

    def test_selected_metadata_including_root(self, annotated_print_tree):
        result = annotated_print_tree.newick(
            print_meta=True, node_attributes=("height", "height_95%_CI"),
            branch_attributes=("length_sd", "length_95%_CI"), precision=4)
        assert "A[&height=1.235,height_95%_CI={0.8,1.5}]:0.1235" in result
        assert "old:3[&length_sd=0.25,length_95%_CI={2.0,4.0}]" in result
        assert result.endswith(")[&height=5.0,height_95%_CI={4.0,6.0}];")
        assert "height_sd" not in result
        assert "height_80%_CI" not in result
        assert "B:2" in result

    @pytest.mark.parametrize("print_meta", [False, None])
    def test_metadata_suppression(self, annotated_print_tree, print_meta):
        result = annotated_print_tree.newick(
            print_meta=print_meta, node_attributes=("height", "height_95%_CI"),
            branch_attributes=("length_sd", "length_95%_CI"), )
        assert "[&" not in result

    def test_empty_selections(self, annotated_print_tree):
        tree = annotated_print_tree
        tree.newick(print_meta=True, node_attributes=("height", "height_95%_CI"), )
        result = tree.newick(print_meta=True, node_attributes=(), branch_attributes=())
        assert "[&" not in result
        result = tree.newick(print_meta=True, node_attributes=("height",), branch_attributes=())
        assert "height=" in result
        assert "_CI=" not in result

    def test_automatic_metadata_and_all_intervals(self, annotated_print_tree):
        tree = annotated_print_tree
        branch = tree.child_dict[tree.root][tree.find_mrca({"A", "B"})]
        branch.n = 10
        branch.mean = 3.0
        branch.M2 = 2.0
        branch.parent_height = 5.0
        branch.kid_height = 2.0
        branch.quantiles = pt.QuantileAccumulator()
        result = tree.newick(print_meta=True)
        assert "height=1.23456789,height_sd=0.1" in result
        assert "posterior=0.975" in result
        assert "height_80%_CI={1.0,1.4}" in result
        assert "height_95%_CI={0.8,1.5}" in result
        assert "length_80%_CI={2.5,3.5}" in result
        assert "length_95%_CI={2.0,4.0}" in result
        assert result.endswith(")[&height=5.0,height_80%_CI={4.5,5.5},height_95%_CI={4.0,6.0}];")
        for field in ("n", "mean", "M2", "quantiles", "parent_height", "kid_height"):
            assert f",{field}=" not in result
            assert f"[&{field}=" not in result
        assert tree.newick(print_meta=True, node_attributes=None, branch_attributes=None,
                           ) == result

    def test_explicit_attributes_do_not_add_intervals(self, annotated_print_tree):
        result = annotated_print_tree.newick(
            print_meta=True, node_attributes=("height",), branch_attributes=("posterior",),
            )
        assert "[&height=1.23456789]" in result
        assert "[&posterior=0.975]" in result
        assert "_CI=" not in result
        assert "height_sd=" not in result
        assert "length_sd=" not in result

    def test_node_and_branch_selections_are_independent(self, annotated_print_tree):
        tree = annotated_print_tree
        nodes_only = tree.newick(print_meta=True, branch_attributes=())
        branches_only = tree.newick(print_meta=True, node_attributes=())
        assert "height_95%_CI" in nodes_only
        assert "length_95%_CI" not in nodes_only
        assert "length_95%_CI" in branches_only
        assert "height_95%_CI" not in branches_only
        assert "height=" not in branches_only

    def test_missing_intervals_are_skipped(self, annotated_print_tree):
        result = annotated_print_tree.newick(print_meta=True, node_attributes=("height_90%_CI", "height_95%_CI"), branch_attributes=("length_90%_CI", "length_95%_CI"))
        assert "height_90%_CI" not in result
        assert "height_95%_CI" in result
        assert "height_80%_CI" not in result

    def test_zero_false_and_collection_metadata(self, annotated_print_tree):
        tree = annotated_print_tree
        tree.set_node_attribute("B", "height", 0.0)
        tree.set_node_attribute("B", "wasambig", False)
        tree.set_node_attribute("B", "state", ('human host', 'a,"quoted"'))
        tree.set_node_attribute("B", "primary_set", {"mink", "human"})
        result = tree.newick(print_meta=True)
        assert 'B[&height=0.0,state={"human host","a,""quoted"""},primary_set={"human","mink"},wasambig=false]' in result
        assert "state=," not in result
        assert "secondary_set=" not in result

    def test_explicit_none_scalar_defaults(self, annotated_print_tree):
        tree = annotated_print_tree
        assert tree.newick(printdist=None, printlabels=None, labelfield=None,
                           precision=None, print_meta=None) == tree.newick()


class TestTreeNexus:

    @pytest.mark.parametrize("rooted, declaration", [(None, ""), (True, "[&R] "), (False, "[&U] ")])
    @pytest.mark.parametrize("translateblock", [False, True])
    @pytest.mark.parametrize("print_meta", [False, True])
    def test_rootedness_declaration(self, annotated_print_tree, rooted, declaration,
                                   translateblock, print_meta):
        tree = annotated_print_tree
        tree.rooted = rooted
        plain = tree.nexus()
        root = tree.root
        result = tree.nexus(translateblock=translateblock, print_meta=print_meta)
        expected = tree.newick(print_meta=print_meta,
                              transdict=tree.transdict() if translateblock else None)
        assert "    tree nexus_tree = " + declaration + expected in result
        assert tree.root == root
        assert tree.nexus() == plain
        if rooted is None:
            assert "[&R]" not in plain and "[&U]" not in plain
        assert tree.rooted is rooted
        assert ("height=" in result) == print_meta
        parsed = next(pt.Nexustreefile(filecontent=result))
        assert parsed.equals(tree, rooted=True)


    @pytest.mark.parametrize("sample", ["((A:1,B:1):1,C:2);", "(A:1,B:1,C:1);"])
    def test_rootedness_is_not_inferred(self, sample):
        tree = pt.Tree.from_string(sample)
        assert "[&R]" not in tree.nexus()
        assert "[&U]" not in tree.nexus()
        tree.rooted = True
        assert "[&R] " in tree.nexus()
        tree.rooted = False
        assert "[&U] " in tree.nexus()


    @pytest.mark.parametrize("rooted", [0, 1, "true", [], {}])
    def test_invalid_rootedness(self, annotated_print_tree, rooted):
        with pytest.raises(pt.TreeError, match="print_rootedness must be a boolean"):
            annotated_print_tree.nexus(print_rootedness=rooted)


    def test_automatic_metadata(self, annotated_print_tree):
        tree = annotated_print_tree
        result = tree.nexus(print_meta=True)
        assert tree.newick(print_meta=True) in result
        assert "height_95%_CI" in result
        assert "length_95%_CI" in result

    def test_forwarding_and_call_independence(self, annotated_print_tree):
        tree = annotated_print_tree
        plain = tree.nexus()
        options = dict(print_meta=True, printdist=True, labelfield="posterior", precision=3,
                       node_attributes=("height", "height_95%_CI"), branch_attributes=("length_sd", "length_95%_CI"),
                       )
        result = tree.nexus(**options)
        assert "    tree nexus_tree = " + tree.newick(**options) in result
        assert result.startswith("#NEXUS")
        assert result.endswith("\nend;\n")
        assert "height_95%_CI" in result
        assert "length_95%_CI" in result
        assert tree.nexus() == plain
        assert "[&" not in plain

    def test_translate_block_and_root_metadata(self, annotated_print_tree):
        tree = annotated_print_tree
        result = tree.nexus(translateblock=True, print_meta=True, node_attributes=("height",))
        expected = tree.newick(transdict=tree.transdict(), print_meta=True, node_attributes=("height",))
        assert "translate" in result.lower()
        assert expected in result
        assert expected.endswith(")[&height=5.0];")
        parsed = next(pt.Nexustreefile(filecontent=result))
        assert parsed.equals(tree, rooted=True)

    def test_suppression_and_none_defaults(self, annotated_print_tree):
        tree = annotated_print_tree
        assert tree.nexus(printdist=None, printlabels=None, labelfield=None,
                          precision=None, print_meta=None) == tree.nexus()
        assert "[&" not in tree.nexus(print_meta=False, node_attributes=("height",))

    @pytest.mark.parametrize('rooted,declaration', [(True, '[&R]'), (False, '[&U]'), (None, None)])
    def test_full_nexus_retains_state_while_plain_newick_omits_declaration(self, rooted, declaration):
        text = ('(A[&height=1.0000000000000002,tag={true,1,1.0,?}]'
                ':0.10000000000000002[&posterior=0.8,length_95%_CI={0.05,0.2}],'
                'B:2)Root[&height=3]:0.5[&rate=1e-10];')
        tree = pt.Tree.from_string(text, parse_meta=True)
        tree.rooted = rooted
        options = dict(print_meta=True, node_attributes=('height', 'tag'),
                       branch_attributes=('posterior', 'length_95%_CI', 'rate'))
        nexus_text = tree.nexus(**options)
        if declaration is None:
            assert '[&R]' not in nexus_text and '[&U]' not in nexus_text
        else:
            assert declaration in nexus_text
        with pt.Nexustreefile(filecontent=nexus_text, parse_meta=True) as reader:
            decoded = reader.readtree()
            assert reader.readtree() is None
        assert decoded.rooted is rooted
        newick_text = tree.newick(**options)
        assert '[&R]' not in newick_text and '[&U]' not in newick_text
        plain = pt.Tree.from_string(newick_text, parse_meta=True)
        assert plain.rooted is None
        for current in (decoded, plain):
            assert current.nodedict['A'].height == 1.0000000000000002
            tag = current.nodedict['A'].extra_meta['tag']
            assert tag == (True, 1, 1.0, pt.MISSING)
            assert tag[0] is True and type(tag[1]) is int and type(tag[2]) is float
            assert tag[3] is pt.MISSING
            branch = current.child_dict[current.root]['A']
            assert branch.length == 0.10000000000000002 and branch.posterior == 0.8
            assert branch.extra_meta['length_95%_CI'] == (0.05, 0.2)
            assert current.root_suffix.label == 'Root' and current.root_suffix.length == 0.5
            assert current.root_suffix.extra_meta['rate'] == 1e-10
        assert tree.rooted is rooted  # Export does not change interpretation.
        with pt.Nexustreefile(filecontent=tree.nexus(print_rootedness=False, **options),
                              parse_meta=True) as reader:
            assert reader.readtree().rooted is None


class TestTreeSetNewick:

    def test_metadata_and_ci_forwarding(self, annotated_print_tree):
        tree = annotated_print_tree
        trees = pt.TreeSet()
        trees.addtree(tree)
        trees.addtree(tree.copy_treeobject())
        options = dict(print_meta=True, printlabels=False, printdist=True, precision=4,
                       node_attributes=("height", "height_95%_CI"), branch_attributes=("length_sd", "length_95%_CI"),
                       transdict={"A": "1", "B": "2", "C": "3"})
        assert trees.newick(**options) == (tree.newick(**options) + "\n") * 2
        assert "[&" not in trees.newick()

    def test_default_precision(self, annotated_print_tree):
        trees = pt.TreeSet()
        trees.addtree(annotated_print_tree)
        assert trees.newick() == annotated_print_tree.newick() + "\n"
        assert trees.newick(precision=None) == trees.newick()
        assert trees.newick(precision=6) == annotated_print_tree.newick(precision=6) + "\n"
        assert trees.newick(print_meta=True) == annotated_print_tree.newick(print_meta=True) + "\n"


class TestTreeSetNexus:

    @pytest.mark.parametrize("rooted, declaration", [(None, ""), (True, "[&R] "), (False, "[&U] ")])
    @pytest.mark.parametrize("translateblock", [False, True])
    @pytest.mark.parametrize("print_meta", [False, True])
    def test_rootedness_declaration(self, annotated_print_tree, rooted, declaration,
                                   translateblock, print_meta):
        tree = annotated_print_tree
        tree.rooted = rooted
        trees = pt.TreeSet()
        trees.addtree(tree)
        trees.addtree(tree.copy_treeobject())
        plain = trees.nexus()
        result = trees.nexus(translateblock=translateblock, print_meta=print_meta)
        expected = tree.newick(print_meta=print_meta,
                              transdict=tree.transdict() if translateblock else None)
        for i in (1, 2):
            assert f"    tree t.{i} = " + declaration + expected in result
        parsed = list(pt.Nexustreefile(filecontent=result))
        assert len(parsed) == 2
        assert all(t.equals(tree, rooted=True) for t in parsed)
        assert trees.nexus() == plain
        if rooted is None:
            assert "[&R]" not in plain and "[&U]" not in plain
        assert tree.rooted is rooted


    @pytest.mark.parametrize("rooted", [0, 1, "false", [], {}])
    def test_invalid_rootedness(self, rooted):
        with pytest.raises(pt.TreeError, match="print_rootedness must be a boolean"):
            pt.TreeSet().nexus(print_rootedness=rooted)


    def test_metadata_ci_and_translate_forwarding(self, annotated_print_tree):
        tree = annotated_print_tree
        trees = pt.TreeSet()
        trees.addtree(tree)
        trees.addtree(tree.copy_treeobject())
        options = dict(print_meta=True, printdist=True, printlabels=False, precision=4,
                       node_attributes=("height", "height_95%_CI"), branch_attributes=("length_sd", "length_95%_CI"),
                       )
        result = trees.nexus(translateblock=True, **options)
        expected = tree.newick(transdict=tree.transdict(), **options)
        assert "translate" in result.lower()
        assert result.count(expected) == 2
        assert len(list(pt.Nexustreefile(filecontent=result))) == 2
        assert "[&" not in trees.nexus()

    def test_default_precision(self, annotated_print_tree):
        trees = pt.TreeSet()
        trees.addtree(annotated_print_tree)
        assert annotated_print_tree.newick() in trees.nexus()
        assert trees.nexus(precision=None) == trees.nexus()
        assert annotated_print_tree.newick(precision=6) in trees.nexus(precision=6)
        assert annotated_print_tree.newick(print_meta=True) in trees.nexus(print_meta=True)


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


class Test_rename_leaf:

    @pytest.mark.parametrize("requested,fixdups,final", [
        ("X", False, "X"), ("B", True, "B_2"),
    ])
    def test_metadata_follows_resolved_name(self, requested, fixdups, final):
        tree = pt.Tree.from_string("(A:1,B:2,B_1:3,C:4);")
        tree.set_node_attribute("A", "height", 5.0)
        tree.set_node_attribute("C", "state", "untouched")
        moved = tree.nodedict["A"]
        untouched = tree.nodedict["C"]

        tree.rename_leaf("A", requested, fixdups=fixdups)

        assert "A" not in tree.nodes
        assert "A" not in tree.nodedict
        assert tree.nodedict[final] is moved
        assert tree.get_node_attribute(final, "height") == 5.0
        assert tree.nodedict["C"] is untouched
        assert set(tree.nodedict) == {final, "C"}
        assert tree.nodedist(final, "C") == pytest.approx(5)

    def test_unannotated_tree_stays_sparse(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.rename_leaf("A", "X")
        assert tree.nodedict == {}


class Test_rename_intnode:

    @pytest.mark.parametrize("rename_root", [False, True])
    def test_metadata_follows_node_id(self, rename_root):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4);")
        old = tree.root if rename_root else tree.find_mrca({"A", "B"})
        tree.set_node_attribute(old, "height", 5.0)
        tree.set_node_attribute("A", "state", "untouched")
        moved = tree.nodedict[old]
        untouched = tree.nodedict["A"]

        tree.rename_intnode(old, 99)

        assert old not in tree.nodes
        assert old not in tree.nodedict
        assert tree.nodedict[99] is moved
        assert tree.get_node_attribute(99, "height") == 5.0
        assert tree.nodedict["A"] is untouched
        assert set(tree.nodedict) == {99, "A"}
        assert tree.nodedist("A", "C") == pytest.approx(8)
        if rename_root:
            assert tree.root == 99

    def test_unannotated_tree_stays_sparse(self):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4);")
        tree.rename_intnode(tree.root, 99)
        assert tree.nodedict == {}

    @pytest.mark.parametrize("state", [True, False, None])
    def test_root_identifier_change_preserves_interpretation_suffix_and_node_data(self, state):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        tree.rooted = state
        tree.root_suffix = pt.Branchstruct(length=0.5, label="same position")
        suffix = tree.root_suffix
        tree.set_node_attribute(tree.root, "extra_meta", {"root": [pt.MISSING]})
        record = tree.nodedict[tree.root]
        before = tree.newick()
        oldroot = tree.root
        newroot = max(tree.intnodes) + 10
        tree.rename_intnode(oldroot, newroot)
        assert tree.root == newroot and tree.rooted is state
        assert tree.root_suffix is suffix
        assert tree.nodedict[newroot] is record and oldroot not in tree.nodedict
        assert tree.newick() == before


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


    def test_supported_mutable_annotation_values_are_independent(self):
        source = pt.Branchstruct(length=3.0, label="edge")
        source.extra_meta = _copy_foundation_payload()

        copied = source.copy()

        assert copied is not source
        assert copied.length == 3.0
        assert copied.label == "edge"
        assert copied.extra_meta == source.extra_meta
        assert copied.extra_meta is not source.extra_meta
        assert isinstance(copied.extra_meta["height_95%_CI"], list)
        assert isinstance(copied.extra_meta["nested"], tuple)
        assert isinstance(copied.extra_meta["nested"][0], list)
        assert isinstance(copied.extra_meta["nested"][0][0], set)
        assert isinstance(copied.extra_meta["frozen"], frozenset)
        assert copied.extra_meta["missing"] is pt.MISSING
        assert copied.extra_meta["nested"][1][1] is pt.MISSING

        copied.extra_meta["height_95%_CI"][0] = 9.0
        copied.extra_meta["nested"][0][0].add("blue")
        copied.extra_meta["new_key"] = "copy only"

        assert source.extra_meta["height_95%_CI"] == [1.0, 2.0]
        assert source.extra_meta["nested"][0][0] == {"red"}
        assert "new_key" not in source.extra_meta

    def test_omitted_attributes_do_not_retain_annotation_mapping(self):
        source = pt.Branchstruct(length=3.0, label="edge")
        source.extra_meta = _copy_foundation_payload()

        copied = source.copy(copyattr=False)

        assert copied.length == 3.0
        assert copied.label == ""
        assert copied.extra_meta is None

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

    def test_beast_ca_fixture_native_and_generic_owners_roundtrip(self, annotation_producer_input):
        text, translation = annotation_producer_input('beast/random_10tips_1000.treeannot_mcc_ca')
        tree = pt.Tree.from_string(text, transdict=translation, parse_meta=True)
        assert len(tree.leaves) == 10 and 's00' in tree.leaves
        assert tree.rooted is None  # This actual fixture supplies no declaration.
        assert tree.root_suffix.length == 0.0
        assert tree._nodeheightdict is None
        encoded = tree.newick(
            print_meta=True,
            node_attributes=('height', 'height_median', 'height_95%_HPD',
                             'CAheight_95%_HPD', 'length', 'length_median',
                             'length_95%_HPD', 'posterior'),
            branch_attributes=(),
        )
        decoded = pt.Tree.from_string(encoded, parse_meta=True)
        for current in (tree, decoded):
            node = current.nodedict['s00']
            branch = current.child_dict[current.parent('s00')]['s00']
            assert node.height == 0.3894414960000004 and type(node.height) is float
            assert node.height_median == 0.34226349999999994
            assert node.extra_meta['height_95%_HPD'] == (0.0, 0.968169)
            assert type(node.extra_meta['height_95%_HPD']) is tuple
            assert all(type(value) is float for value in node.extra_meta['height_95%_HPD'])
            assert node.extra_meta['CAheight_95%_HPD'] == (0.0, 0.968169)
            assert node.extra_meta['length'] == 0.20585859039999968
            assert node.extra_meta['length_95%_HPD'] == (0.09111859999999994, 0.317102)
            assert node.extra_meta['length_median'] == 0.19533300000000003
            assert branch.length == 0.5695019715  # Colon length remains independent.
            assert branch.length_median is None and branch.extra_meta is None
            assert current.nodedict[current.root].height == 1.1903651784999993
            assert current.root_suffix.length == 0.0 and current.rooted is None

    def test_beast_mean_fixture_preserves_leading_declaration_and_precolon_statistics(self, annotation_producer_input):
        text, translation = annotation_producer_input('beast/random_10tips_1000.treeannot_mcc_mean')
        tree = pt.Tree.from_string(text, transdict=translation, parse_meta=True)
        assert tree.rooted is True and tree.root_suffix is None
        encoded = tree.newick(
            print_meta=True,
            node_attributes=('height_median', 'height_mean', 'height_95%_HPD',
                             'length', 'length_median', 'length_95%_HPD', 'posterior'),
            branch_attributes=(),
        )
        decoded = pt.Tree.from_string(encoded, parse_meta=True)
        assert decoded.rooted is None  # Plain Newick deliberately carries no declaration.
        for current in (tree, decoded):
            node = current.nodedict['s03']
            branch = current.child_dict[current.parent('s03')]['s03']
            assert node.height is None and node.height_median == 0.32481499999999996
            assert node.extra_meta['height_mean'] == 0.3763489003999995
            assert node.extra_meta['length'] == 0.3763489003999995
            assert branch.length == 0.13221569840952402
            assert branch.length_median is None and branch.posterior is None
            cherry = current.find_mrca({'s03', 's04'})
            assert set(current.children(cherry)) == {'s03', 's04'}
            assert current.nodedict[cherry].extra_meta['posterior'] == 0.084
            assert current.nodedict[cherry].extra_meta['height_95%_HPD'] == (0.13973200000000008, 1.092982)
            assert current.child_dict[current.parent(cherry)][cherry].posterior is None
            assert current.root_suffix is None

    def test_mrbayes_summary_fixture_translates_and_converts_branch_labels(self, annotation_producer_input):
        text, translation = annotation_producer_input('mrbayes/contest.nexus.con.tre')
        # Existing translation/format/label settings keep their positional slots.
        tree = pt.Tree.from_string(text, translation, 'newick', 'posterior', float,
                                   parse_meta=True)
        assert len(tree.leaves) == 40 and {'u16376', 'u16382'} <= tree.leaves
        assert tree.rooted is None and tree.nodedict == {}
        encoded = tree.newick(printlabels=False, print_meta=True,
                              node_attributes=(), branch_attributes=('posterior',))
        decoded = pt.Tree.from_string(encoded, parse_meta=True)
        for current in (tree, decoded):
            cherry = current.find_mrca({'u16376', 'u16382'})
            assert set(current.children(cherry)) == {'u16376', 'u16382'}
            branch = current.child_dict[current.parent(cherry)][cherry]
            assert branch.posterior == 0.914 and type(branch.posterior) is float
            assert branch.length == 3.010486e-3
            assert current.child_dict[cherry]['u16376'].length == 1.357325e-3
            assert current.rooted is None and current.nodedict == {}
            assert current.root_suffix is None
        assert '[&R]' not in tree.nexus() and '[&U]' not in tree.nexus()

    def test_annotation_input_preserves_transdict_and_root_label_converter_settings(self):
        translation = {'1': 'Alpha', '2': 'Beta', '3': 'Gamma'}
        text = ('[&U]((1[&tag=1]:1,2:2)0.75[&note="node"]:3,3:4)0.25:0;')
        tree = pt.Tree.from_string(text, translation, 'newick', 'freq', float, parse_meta=True)
        assert tree.leaves == {'Alpha', 'Beta', 'Gamma'}
        assert tree.rooted is False and tree.nodedict['Alpha'].extra_meta['tag'] == 1
        assert type(tree.nodedict['Alpha'].extra_meta['tag']) is int
        cherry = tree.find_mrca({'Alpha', 'Beta'})
        assert tree.nodedict[cherry].extra_meta['note'] == 'node'
        assert tree.child_dict[tree.root][cherry].freq == 0.75
        assert type(tree.child_dict[tree.root][cherry].freq) is float
        assert tree.root_suffix.freq == 0.25 and type(tree.root_suffix.freq) is float
        assert tree.root_suffix.length == 0.0
        encoded = tree.newick(labelfield='freq', print_meta=True,
                              node_attributes=('tag', 'note'), branch_attributes=())
        decoded = pt.Tree.from_string(encoded, label_attr_name='freq', label_type=float,
                                      parse_meta=True)
        cherry = decoded.find_mrca({'Alpha', 'Beta'})
        assert decoded.child_dict[decoded.root][cherry].freq == 0.75
        assert decoded.root_suffix.freq == 0.25 and decoded.root_suffix.length == 0.0
        assert decoded.nodedict['Alpha'].extra_meta['tag'] == 1
        assert decoded.nodedict[cherry].extra_meta['note'] == 'node'
        assert decoded.rooted is None

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

    def test_bipartition_construction_establishes_unrooted_state(self):
        source = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        source.rooted = True
        result = pt.Tree.from_biplist(source.bipdict())
        assert result.rooted is False
        assert result.topology_bipart == source.topology_bipart
        assert source.rooted is True

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

    def test_bipartition_topology_establishes_unrooted_state(self):
        source = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        result = pt.Tree.from_topology(source.topology_bipart)
        assert result.rooted is False
        assert result.topology_bipart == source.topology_bipart

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

    def test_generic_star_construction_is_unspecified(self):
        assert pt.Tree.from_leaves(["A", "B", "C"]).rooted is None

###################################################################################################

class Test_MISSING:

    def test_copy_and_pickle_preserve_singleton_identity(self):
        assert repr(pt.MISSING) == "MISSING"
        assert copy.copy(pt.MISSING) is pt.MISSING
        assert copy.deepcopy(pt.MISSING) is pt.MISSING
        assert pickle.loads(pickle.dumps(pt.MISSING)) is pt.MISSING


class Test_optional_annotation_defaults:

    def test_unset_state_and_sd_are_distinct_from_empty_and_zero(self):
        node = pt.Nodestruct()
        branch = pt.Branchstruct()
        assert node.state is None
        assert branch.length_sd is None
        assert pt.Branchstruct.default_for("length_sd") is None
        node.state = ""
        branch.length_sd = 0.0
        assert node.state == ""
        assert branch.length_sd == 0.0
        tree = pt.Tree.from_string("(A:1,B:1);")
        tree.ensure_nodedict()
        with pytest.raises(pt.TreeError, match="no nodes.*state"):
            tree.parsimony_possible_states()


class Test_set_blens_from_heights:

    @pytest.mark.parametrize("missing_node", ["root", "B"])
    def test_explicitly_missing_height_fails_before_length_writes(self, missing_node):
        tree = pt.Tree.from_string("(A:1,B:1);")
        tree.set_node_attribute(tree.root, "height", 3.0)
        tree.set_node_attribute("A", "height", 0.0)
        tree.set_node_attribute("B", "height", 0.0)
        node = tree.root if missing_node == "root" else missing_node
        tree.set_node_attribute(node, "height", pt.MISSING)
        before = tree.newick()
        with pytest.raises(pt.TreeError, match="height"):
            tree.set_blens_from_heights()
        assert tree.newick() == before
        assert tree.get_node_attribute(node, "height") is pt.MISSING

    def test_explicit_height_conversion_preserves_other_stored_data(self, length_edit_tree):
        tree = length_edit_tree
        for record in tree.nodedict.values():
            record.height = record.height * 2 + 1
        before = _length_edit_data(tree)
        _ = tree.rootdist, tree.nodeheightdict
        tree.set_blens_from_heights()
        _assert_length_edit_data(tree, before)
        assert tree.rootdist["A"] == 6
        assert tree.nodeheightdict[tree.root] == 6
        assert tree.nodedict[tree.root].height == 7


class Test_init_Nodestruct:

    def test_leaf_count_is_not_ordinary_node_metadata(self):
        assert not hasattr(pt.Nodestruct(), "nleaves")
        with pytest.raises(TypeError):
            pt.Nodestruct(nleaves=4)
        with pytest.raises(TypeError):
            pt.Nodestruct(1.0, 4)


class Test_from_cladedict:

    def test_summary_leaf_counts_do_not_leak_into_tree_metadata(self):
        original = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        clades = original.cladedict()
        root_clade = pt.Clade.from_leafset(original.leaves, original)
        assert clades[root_clade].nleaves == 4

        tree = pt.Tree.from_cladedict(clades)

        assert all(not hasattr(node, "nleaves") for node in tree.nodedict.values())
        assert "nleaves=" not in tree.newick(print_meta=True)
        tree.remove_leaf("B")
        remaining = tree.cladedict()
        assert remaining[pt.Clade.from_leafset(tree.leaves, tree)].nleaves == 3
        assert clades[root_clade].nleaves == 4

    def test_clade_hierarchy_establishes_rooted_state(self):
        source = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        source.rooted = False
        result = pt.Tree.from_cladedict(source.cladedict())
        assert result.rooted is True
        assert result.topology_clade == source.topology_clade
        assert source.rooted is False


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

    def test_directed_branch_table_does_not_imply_rootedness(self):
        tree = pt.Tree.from_branchinfo(
            ["r", "r", "x", "x"], ["x", "C", "A", "B"], [2.0, 3.0, 1.0, 1.0],
        )
        assert tree.rooted is None

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

    @pytest.mark.parametrize("model", ["yule", "coalescent"])
    @pytest.mark.parametrize("rate_sd", [0.0, 0.5])
    def test_ancestral_generation_is_rooted_despite_rate_variation(self, model, rate_sd):
        # Deterministic local draws without leaking RNG state into unrelated tests.
        before_rng = random.getstate()
        try:
            random.seed(411)
            tree = pt.Tree.randtree(ntips=6, model=model, tree_height=10.0, rate_sd=rate_sd)
        finally:
            random.setstate(before_rng)
        assert tree.rooted is True
        assert len(tree.leaves) == 6

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


    @pytest.mark.parametrize("rooted", [True, False, None])
    @pytest.mark.parametrize("copylengths,copyattr", [
        (True, True), (True, False), (False, True), (False, False),
    ])
    def test_root_interpretation_and_identity_survive_all_copy_flags(
        self, rooted, copylengths, copyattr,
    ):
        source = pt.Tree.from_string("(A:2,B:3);")
        source.rooted = rooted
        source.set_node_attribute(source.root, "height", 8.0)
        source.set_node_attribute("A", "state", "red")
        source.set_branch_attribute(source.root, "A", "label", "edge")

        copied = source.copy_treeobject(copylengths=copylengths, copyattr=copyattr)

        assert copied is not source
        assert copied.root == source.root
        assert copied.leaves == source.leaves
        assert copied.nodes == source.nodes
        assert copied.topology_clade == source.topology_clade
        assert copied.rooted is rooted
        assert copied.root_suffix is None
        assert copied.nodedist("A", "B") == (5.0 if copylengths else 0.0)
        assert copied.get_branch_attribute(copied.root, "A", "label") == (
            "edge" if copyattr else ""
        )
        if copyattr:
            assert copied.nodedict[source.root].height == 8.0
            assert copied.nodedict["A"].state == "red"
            assert copied.nodedict[source.root] is not source.nodedict[source.root]
        else:
            assert copied.nodedict == {}
        assert source.nodedict[source.root].height == 8.0
        assert source.nodedict["A"].state == "red"
        assert source.nodedist("A", "B") == 5.0

    def test_full_copy_isolates_node_branch_and_suffix_annotation_values(self):
        source = pt.Tree.from_string("(A:2,B:3);")
        source.set_node_attribute(source.root, "extra_meta", _copy_foundation_payload())
        source.set_branch_attribute(source.root, "A", "extra_meta", _copy_foundation_payload())
        source.root_suffix = pt.Branchstruct(length=0.0, label="root_label")
        source.root_suffix.extra_meta = _copy_foundation_payload()

        copied = source.copy_treeobject()
        copied.nodedict[copied.root].extra_meta["nested"][0][0].add("node_copy")
        copied.child_dict[copied.root]["A"].extra_meta["height_95%_CI"][0] = 11.0
        copied.root_suffix.extra_meta["nested"][0][0].add("stem_copy")

        assert source.nodedict[source.root].extra_meta["nested"][0][0] == {"red"}
        assert source.child_dict[source.root]["A"].extra_meta["height_95%_CI"] == [1.0, 2.0]
        assert source.root_suffix.extra_meta["nested"][0][0] == {"red"}
        assert copied.root_suffix is not source.root_suffix
        assert copied.root_suffix.length == 0.0
        assert copied.root_suffix.label == "root_label"
        assert copied.root_suffix.extra_meta["missing"] is pt.MISSING
        assert copied.nodedist("A", "B") == source.nodedist("A", "B") == 5.0

    @pytest.mark.parametrize(
        "stem_length,copylengths,copyattr,expected_length,expected_present",
        [
            (None, True, True, None, True),
            (None, False, True, None, True),
            (None, True, False, None, False),
            (None, False, False, None, False),
            (0.0, True, False, 0.0, True),
            (2.5, True, False, 2.5, True),
            (2.5, False, False, 0.0, True),
            (2.5, False, True, 0.0, True),
        ],
    )
    def test_suffix_copy_flags_distinguish_absent_and_present_stem_lengths(
        self, stem_length, copylengths, copyattr, expected_length, expected_present,
    ):
        source = pt.Tree.from_string("(A:2,B:3);")
        source.root_suffix = pt.Branchstruct(length=stem_length, label="root_label")
        source.root_suffix.extra_meta = {"rate": [1.0, 2.0]}

        copied = source.copy_treeobject(copylengths=copylengths, copyattr=copyattr)

        if not expected_present:
            assert copied.root_suffix is None
        else:
            assert copied.root_suffix is not source.root_suffix
            assert copied.root_suffix.length == expected_length
            if copyattr:
                assert copied.root_suffix.label == "root_label"
                assert copied.root_suffix.extra_meta == {"rate": [1.0, 2.0]}
                copied.root_suffix.extra_meta["rate"][0] = 99.0
            else:
                assert copied.root_suffix.label == ""
                assert copied.root_suffix.extra_meta is None
        assert source.root_suffix.length == stem_length
        assert source.root_suffix.label == "root_label"
        assert source.root_suffix.extra_meta == {"rate": [1.0, 2.0]}

    def test_suffix_with_only_present_empty_or_missing_annotations_is_retained(self):
        source = pt.Tree.from_string("(A:2,B:3);")
        source.root_suffix = pt.Branchstruct(length=None)
        source.root_suffix.extra_meta = {
            "empty_string": "", "empty_collection": (),
            "zero": 0, "false": False, "missing": pt.MISSING,
        }

        copied = source.copy_treeobject()

        assert copied.root_suffix is not None
        assert copied.root_suffix.length is None
        assert copied.root_suffix.extra_meta == source.root_suffix.extra_meta
        assert copied.root_suffix.extra_meta is not source.root_suffix.extra_meta
        assert copied.root_suffix.extra_meta["missing"] is pt.MISSING

    @pytest.mark.parametrize("field,value", [
        ("mean", 0.9), ("M2", 0.25), ("n", 3),
        ("parent_height", 0.0), ("kid_height", pt.MISSING),
    ])
    def test_suffix_retains_nondefault_native_fields(self, field, value):
        source = pt.Tree.from_string("(A:2,B:3);")
        source.root_suffix = pt.Branchstruct(length=None)
        setattr(source.root_suffix, field, value)

        copied = source.copy_treeobject()

        assert copied.root_suffix is not None
        assert copied.root_suffix is not source.root_suffix
        assert copied.root_suffix.length is None
        if value is pt.MISSING:
            assert getattr(copied.root_suffix, field) is pt.MISSING
        else:
            assert getattr(copied.root_suffix, field) == value
        assert source.copy_treeobject(copyattr=False).root_suffix is None
        assert getattr(source.root_suffix, field) is value

    def test_suffix_configured_native_label_survives_copy(self):
        source = pt.Tree.from_string("(A:2,B:3);")
        source.root_suffix = pt.Branchstruct(length=None)
        source.root_suffix.mean = 0.9

        copied = source.copy_treeobject()

        assert source.newick(labelfield="mean") == "(A:2.0,B:3.0)0.9;"
        assert copied.newick(labelfield="mean") == source.newick(labelfield="mean")

    def test_suffix_with_no_retained_information_is_omitted(self):
        source = pt.Tree.from_string("(A:2,B:3);")
        source.root_suffix = pt.Branchstruct(length=None)
        source.root_suffix.extra_meta = {}

        copied = source.copy_treeobject()

        assert copied.root_suffix is None
        assert source.root_suffix is not None
        assert source.root_suffix.length is None
        assert source.root_suffix.extra_meta == {}

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
    @pytest.mark.parametrize("annotate_survivor", [False, True])
    def test_moved_leaf_retains_independent_metadata(self, annotate_survivor):
        tree = pt.Tree.from_string("((A:1,B:2):3,(C:4,D:5):6);")
        tree.set_node_attribute("A", "height", 1.25)
        tree.set_node_attribute("A", "state", "sampled")
        tree.set_node_attribute("A", "extra_meta", {"height_95%_CI": (0.5, 2.0)})
        tree.set_node_attribute("A", "primary_set", {"human", "mink"})
        moved = tree.nodedict["A"]
        if annotate_survivor:
            tree.set_node_attribute("C", "height", 9.0)
            survivor = tree.nodedict["C"]
        total_length = tree.length()
        assert tree.nodedist("A", "C") == pytest.approx(14)
        assert tree.rootdist["A"] == pytest.approx(4)

        assert tree.spr("A", "C", frac_from_parent=0.5) == ("A", "C")

        restored = tree.nodedict["A"]
        assert restored is not moved
        assert restored.height == 1.25
        assert restored.state == "sampled"
        assert restored.extra_meta == {"height_95%_CI": (0.5, 2.0)}
        assert restored.primary_set == {"human", "mink"}
        restored.extra_meta["height_95%_CI"] = (8.0, 9.0)
        restored.primary_set.add("changed")
        assert moved.extra_meta == {"height_95%_CI": (0.5, 2.0)}
        assert moved.primary_set == {"human", "mink"}
        graftpoint = tree.parent("A")
        assert tree.nodedict[graftpoint].height is None
        assert tree.nodedict[graftpoint].state is None
        assert tree.leaves == {"A", "B", "C", "D"}
        assert set(tree.nodedict) <= tree.nodes
        assert tree.length() == pytest.approx(total_length)
        assert tree.nodedist("A", "C") == pytest.approx(3)
        assert tree.rootdist["A"] == pytest.approx(9)
        if annotate_survivor:
            assert tree.nodedict["C"] is survivor
            assert survivor.height == 9.0

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

    @pytest.mark.parametrize("state", [True, False, None])
    @pytest.mark.parametrize("promote_root", [False, True])
    def test_complete_spr_keeps_state_and_suffix_tracks_actual_root_survival(self, state, promote_root):
        text = "(A:1,((B:2,C:3):4,D:5):6);" if promote_root else "((A:1,B:2):3,(C:4,D:5):6);"
        target = "B" if promote_root else "C"
        tree = pt.Tree.from_string(text)
        tree.rooted = state
        old_root = tree.root
        tree.root_suffix = pt.Branchstruct(length=0.4, label="oldroot")
        suffix = tree.root_suffix
        tree.set_node_attribute(old_root, "state", "oldroot")
        tree.set_node_attribute("A", "state", "moved")
        tree.set_node_attribute(target, "height", 22.0)
        surviving_metadata = tree.nodedict[target]
        basal = tree.child_dict[tree.parent("A")]["A"]
        basal.label = "do not transfer with reused length"
        basal.extra_meta = {"basal_rate": 7.0}
        target_parent = tree.parent(target)
        tree.child_dict[target_parent][target].label = "split"
        tree.child_dict[target_parent][target].extra_meta = {"target_rate": 8.0}
        total = tree.length()
        before_distances = tree.rootdist
        assert tree.nodedist("A", target) == pytest.approx(13 if promote_root else 14)

        assert tree.spr("A", target, frac_from_parent=0.5) == ("A", target)

        assert tree.rooted is state
        assert tree.leaves == {"A", "B", "C", "D"}
        assert tree.length() == pytest.approx(total)
        assert tree.rootdist is not before_distances
        assert tree.rootdist["A"] == pytest.approx(12 if promote_root else 9)
        assert tree.nodedist("A", target) == pytest.approx(8 if promote_root else 3)
        assert tree.nodedict["A"].state == "moved"
        assert tree.nodedict[target] is surviving_metadata
        graftpoint = tree.parent("A")
        connecting = tree.child_dict[graftpoint]["A"]
        assert connecting.length == pytest.approx(7 if promote_root else 1)
        assert connecting.label == "" and connecting.extra_meta is None
        upper = tree.child_dict[tree.parent(graftpoint)][graftpoint]
        lower = tree.child_dict[graftpoint][target]
        assert upper.label == lower.label == ""
        assert upper.extra_meta is None and lower.extra_meta is None
        assert tree.nodedict[graftpoint].state is None
        assert tree.nodedict[graftpoint].extra_meta is None
        assert set(tree.nodedict) <= tree.nodes
        if promote_root:
            assert tree.root != old_root
            assert old_root not in tree.nodedict
            assert tree.root_suffix is None
        else:
            assert tree.root == old_root
            assert tree.root_suffix is suffix
            merged = tree.child_dict[tree.root]["B"]
            assert merged.length == pytest.approx(5)
            assert merged.label == "" and merged.extra_meta is None

    @pytest.mark.parametrize("state", [True, False, None])
    def test_internal_subtree_metadata_follows_final_ids_and_unrelated_edges_survive(self, state):
        tree = pt.Tree.from_string("(((A:1,B:2):3,C:4):5,(D:6,E:7):8);")
        # Make the detached root collide with the graftpoint ID allocated after pruning.
        original_base = tree.find_mrca({"A", "B"})
        prune_node = max(tree.intnodes) + 1
        tree.rename_intnode(original_base, prune_node)
        removed_parent = tree.parent(prune_node)
        tree.rooted = state
        old_root = tree.root
        tree.root_suffix = pt.Branchstruct(length=0.4, label="root")
        suffix = tree.root_suffix
        tree.set_node_attribute(prune_node, "height", 9.0)
        tree.set_node_attribute(prune_node, "extra_meta", {"node_note": ["moved"]})
        tree.set_node_attribute("A", "state", "leaf moved")
        tree.set_node_attribute(removed_parent, "state", "deleted")
        tree.set_node_attribute("C", "state", "survives")
        c_metadata = tree.nodedict["C"]
        tree.child_dict[prune_node]["A"].extra_meta = {"rate": [1, 2]}
        original_moved_edge = tree.child_dict[prune_node]["A"]
        tree.child_dict[removed_parent][prune_node].extra_meta = {"basal": "omit"}
        tree.child_dict[tree.root][removed_parent].label = "merged upper"
        tree.child_dict[removed_parent]["C"].extra_meta = {"merged lower": "omit"}
        de = tree.parent("D")
        tree.child_dict[de]["D"].extra_meta = {"split": "omit"}
        tree.child_dict[de]["E"].label = "unrelated"
        tree.child_dict[de]["E"].extra_meta = {"keep": [3, 4]}
        unrelated_edge = tree.child_dict[de]["E"]
        total = tree.length()
        assert tree.nodedist("A", "E") == pytest.approx(24)
        assert tree.rootdist["A"] == pytest.approx(9)

        tree.spr(prune_node, "D", frac_from_parent=0.5)

        moved_root = tree.parent("A")
        graftpoint = tree.parent(moved_root)
        assert moved_root != prune_node
        assert graftpoint == prune_node
        assert tree.nodedict[moved_root].height == 9.0
        assert tree.nodedict[moved_root].extra_meta == {"node_note": ["moved"]}
        assert tree.nodedict["A"].state == "leaf moved"
        moved_edge = tree.child_dict[moved_root]["A"]
        assert moved_edge.extra_meta == {"rate": [1, 2]}
        assert moved_edge.extra_meta is not original_moved_edge.extra_meta
        assert tree.child_dict[de]["E"] is unrelated_edge
        assert unrelated_edge.label == "unrelated"
        assert unrelated_edge.extra_meta == {"keep": [3, 4]}
        assert tree.nodedict["C"] is c_metadata
        assert removed_parent not in tree.nodes and removed_parent not in tree.nodedict
        assert tree.nodedict[graftpoint].height is None
        assert tree.nodedict[graftpoint].extra_meta is None
        for parent, child in ((tree.root, "C"), (de, graftpoint), (graftpoint, "D"),
                              (graftpoint, moved_root)):
            branch = tree.child_dict[parent][child]
            assert branch.label == "" and branch.extra_meta is None
        assert tree.child_dict[graftpoint][moved_root].length == pytest.approx(3)
        assert tree.length() == pytest.approx(total)
        assert tree.nodedist("A", "E") == pytest.approx(14)
        assert tree.rootdist["A"] == pytest.approx(15)
        assert tree.root == old_root and tree.root_suffix is suffix
        assert tree.rooted is state
        assert set(tree.nodedict) <= tree.nodes


###################################################################################################

class TestGraft:

    def test_subdivided_target_edge_starts_without_source_annotations(self):
        tree = pt.Tree.from_string("(A:8,B:2,C:3);")
        original = tree.child_dict[tree.root]["A"]
        original.label = "old segment"
        original.extra_meta = {"custom": [1]}
        node = tree.graft("D", tree.root, "A", frac_from_parent=0.25)
        for branch in (tree.child_dict[tree.root][node], tree.child_dict[node]["A"]):
            assert branch.label == "" and branch.extra_meta is None
        assert original.extra_meta == {"custom": [1]}
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

    def test_donor_metadata_survives_an_unannotated_recipient(self):
        recipient = pt.Tree.from_string("(A:2,B:4,C:6);")
        assert recipient.nodedict == {}
        donor = pt.Tree.from_string("(D:1,(E:1,F:1):1);")
        donor_internal = donor.find_mrca({"E", "F"})
        donor.set_node_attribute(donor.root, "height", 3.0)
        donor.set_node_attribute(donor_internal, "height", 2.0)
        donor.set_node_attribute("D", "state", "donor")
        root_meta = donor.nodedict[donor.root]
        internal_meta = donor.nodedict[donor_internal]
        leaf_meta = donor.nodedict["D"]

        graftpoint = recipient.graft(donor, parent=recipient.root, child="A",
                                     frac_from_parent=0.5)

        assert recipient.nodedict[recipient.parent("D")] is root_meta
        assert recipient.nodedict[recipient.find_mrca({"E", "F"})] is internal_meta
        assert recipient.nodedict["D"] is leaf_meta
        assert recipient.get_node_attribute("D", "state") == "donor"
        assert recipient.nodedict[graftpoint].height is None
        assert recipient.nodedict[graftpoint].state is None
        assert "B" not in recipient.nodedict
        assert set(recipient.nodedict) <= recipient.nodes

    def test_metadata_follows_donor_id_and_leaf_renaming(self):
        recipient = pt.Tree.from_string("(A:1,(B:1,C:1):1);")
        recipient.set_node_attribute("A", "state", "recipient")
        recipient_meta = recipient.nodedict["A"]
        donor = pt.Tree.from_string("(D:1,(E:1,F:1):1);")
        donor_internal = donor.find_mrca({"E", "F"})
        donor.set_node_attribute(donor.root, "height", 3.0)
        donor.set_node_attribute(donor_internal, "height", 2.0)
        donor.set_node_attribute("D", "state", "donor")
        root_meta = donor.nodedict[donor.root]
        internal_meta = donor.nodedict[donor_internal]
        leaf_meta = donor.nodedict["D"]

        recipient.graft(donor, parent=recipient.root, child="A",
                        frac_from_parent=0.5, graftlabel="p_")

        assert recipient.nodedict[recipient.parent("p_D")] is root_meta
        assert recipient.nodedict[recipient.find_mrca({"p_E", "p_F"})] is internal_meta
        assert recipient.nodedict["p_D"] is leaf_meta
        assert recipient.get_node_attribute("p_D", "state") == "donor"
        assert recipient.nodedict["A"] is recipient_meta
        assert "D" not in recipient.nodedict
        assert set(recipient.nodedict) <= recipient.nodes

    def test_explicit_copies_preserve_reusable_donor_and_connection_originals(self):
        recipient = pt.Tree.from_string("(A:2,B:4,C:6);")
        donor = pt.Tree.from_string("(X:1,(Y:2,Z:3):4);")
        donor.rooted = True
        donor.root_suffix = pt.Branchstruct(length=0.4, label="original root")
        donor.root_suffix.extra_meta = {"stem": ["original"]}
        original_root = donor.root
        original_suffix = donor.root_suffix
        donor.set_node_attribute(original_root, "height", 5.0)
        donor.set_node_attribute("X", "extra_meta", {"node": ["original"]})
        donor_edge = donor.child_dict[donor.parent("X")]["X"]
        donor_edge.extra_meta = {"rate": [1, 2]}
        assert donor.nodedist("X", "Y") == pytest.approx(7)
        assert donor.rootdist["X"] == pytest.approx(1)
        connection = pt.Branchstruct(length=0.7, label="connecting")
        connection.extra_meta = {"note": ["original"]}
        transferred_connection = connection.copy()

        graftpoint = recipient.graft(
            donor.copy_treeobject(), recipient.root, "A", frac_from_parent=0.5,
            connect_branchstruct=transferred_connection, graftlabel="g_")

        grafted_root = recipient.parent("g_X")
        assert recipient.child_dict[graftpoint][grafted_root] is transferred_connection
        recipient.setlength(grafted_root, "g_X", 9.0)
        recipient.nodedict["g_X"].extra_meta["node"].append("recipient edit")
        recipient.child_dict[grafted_root]["g_X"].extra_meta["rate"].append(3)
        recipient.child_dict[graftpoint][grafted_root].extra_meta["note"].append("recipient edit")

        # These are preserved originals, not the consumed copies supplied above.
        assert donor.leaves == {"X", "Y", "Z"}
        assert donor.root == original_root and donor.rooted is True
        assert donor.root_suffix is original_suffix
        assert donor.root_suffix.extra_meta == {"stem": ["original"]}
        assert donor.nodedict[original_root].height == 5.0
        assert donor.nodedict["X"].extra_meta == {"node": ["original"]}
        assert donor_edge.length == pytest.approx(1)
        assert donor_edge.extra_meta == {"rate": [1, 2]}
        assert donor.nodedist("X", "Y") == pytest.approx(7)
        assert donor.rootdist["X"] == pytest.approx(1)
        assert connection.length == pytest.approx(0.7)
        assert connection.extra_meta == {"note": ["original"]}

    @pytest.mark.parametrize("suffix_kind", ["absent", "label_only", "stem"])
    @pytest.mark.parametrize("connection_length", [None, 0.7])
    def test_donor_suffix_does_not_supply_an_implicit_connection(self, suffix_kind, connection_length):
        recipient = pt.Tree.from_string("(A:2,B:4,C:6);")
        recipient.root_suffix = pt.Branchstruct(length=0.8, label="recipient")
        recipient.root_suffix.extra_meta = {"stem_note": "retained"}
        recipient_suffix = recipient.root_suffix
        recipient_root = recipient.root
        recipient.set_node_attribute(recipient_root, "height", 10.0)
        recipient_root_metadata = recipient.nodedict[recipient_root]
        donor = pt.Tree.from_string("(D:1,(E:2,F:3):4);")
        donor.set_node_attribute(donor.root, "height", 5.0)
        donor_root_metadata = donor.nodedict[donor.root]
        if suffix_kind != "absent":
            length = None if suffix_kind == "label_only" else 0.4
            donor.root_suffix = pt.Branchstruct(length=length, label="donor")
            donor.root_suffix.extra_meta = {"donor_stem_note": "omitted"}
        kwargs = {} if connection_length is None else {"connect_length": connection_length}

        graftpoint = recipient.graft(donor, recipient_root, "A", frac_from_parent=0.5, **kwargs)

        inserted_root = recipient.parent("D")
        connection = recipient.child_dict[graftpoint][inserted_root]
        assert connection.length == pytest.approx(0.0 if connection_length is None else connection_length)
        assert connection.label == ""
        assert connection.extra_meta is None
        assert recipient.root == recipient_root
        assert recipient.root_suffix is recipient_suffix
        assert recipient.nodedict[recipient_root] is recipient_root_metadata
        assert recipient.nodedict[inserted_root] is donor_root_metadata
        assert recipient.child_dict[recipient_root][graftpoint].length + recipient.child_dict[graftpoint]["A"].length == pytest.approx(2)

    @pytest.mark.parametrize("reuse_suffix", [False, True])
    def test_explicit_connection_struct_has_precedence_and_transfers_directly(self, reuse_suffix):
        recipient = pt.Tree.from_string("(A:2,B:4,C:6);")
        donor = pt.Tree.from_string("(D:1,E:2,F:3);")
        donor.root_suffix = pt.Branchstruct(length=None, label="deliberate")
        donor.root_suffix.extra_meta = {"note": ["kept"]}
        if reuse_suffix:
            # A label-only suffix needs an explicitly assigned intended length first.
            connection = donor.root_suffix
            connection.length = 0.4
        else:
            connection = pt.Branchstruct(length=0.9, label="explicit")
            connection.extra_meta = {"note": ["kept"]}
        expected_length = connection.length

        graftpoint = recipient.graft(donor, recipient.root, "A", frac_from_parent=0.5,
                                     connect_length=99.0, connect_branchstruct=connection)

        inserted_root = recipient.parent("D")
        assert recipient.child_dict[graftpoint][inserted_root] is connection
        assert connection.length == expected_length
        assert connection.extra_meta == {"note": ["kept"]}
        assert recipient.root_suffix is None

    def test_leaf_extraction_basal_branch_can_be_reused_explicitly(self):
        source = pt.Tree.from_string("(D:1,(E:2,F:3):4);")
        source.root_suffix = pt.Branchstruct(length=5.0, label="not_the_basal_branch")
        parent = source.parent("D")
        source.child_dict[parent]["D"].label = "basal"
        source.child_dict[parent]["D"].extra_meta = {"rate": 0.7}
        leaf, basal = source.subtree("D")
        assert leaf == "D" and isinstance(leaf, str)
        assert basal is not source.child_dict[parent]["D"]
        recipient = pt.Tree.from_string("(A:2,B:4,C:6);")

        graftpoint = recipient.graft(leaf, recipient.root, "A", frac_from_parent=0.5,
                                     connect_branchstruct=basal)

        assert recipient.child_dict[graftpoint]["D"] is basal
        assert basal.length == pytest.approx(1)
        assert basal.label == "basal"
        assert basal.extra_meta == {"rate": 0.7}
        assert recipient.root_suffix is None

    @pytest.mark.parametrize("recipient_state", [True, False, None])
    @pytest.mark.parametrize("donor_state", [True, False, None])
    def test_graft_preserves_recipient_state_independently_of_donor(self, recipient_state, donor_state):
        recipient = pt.Tree.from_string("(A:2,B:4,C:6);")
        recipient.rooted = recipient_state
        recipient.root_suffix = pt.Branchstruct(length=0.8, label="recipient")
        suffix = recipient.root_suffix
        old_root = recipient.root
        donor = pt.Tree.from_string("(D:1,(E:2,F:3):4);")
        donor.rooted = donor_state
        donor.root_suffix = pt.Branchstruct(length=10.0, label="donor")

        recipient.graft(donor, recipient.root, "A", frac_from_parent=0.5)

        assert recipient.rooted is recipient_state
        assert recipient.root == old_root
        assert recipient.root_suffix is suffix

    @pytest.mark.parametrize("state", [True, False, None])
    def test_leaf_graft_preserves_recipient_state(self, state):
        tree = pt.Tree.from_string("(A:2,B:4,C:6);")
        tree.rooted = state
        tree.graft("D", tree.root, "A", frac_from_parent=0.5, connect_length=0.7)
        assert tree.rooted is state

    def test_graft_uses_state_after_prior_structural_root_collapse(self):
        recipient = pt.Tree.from_string("((A:1,B:2):3,C:4);")
        recipient.rooted = True
        recipient.collapse_bifurcating_root()
        assert recipient.rooted is None  # Separate 4v structural-collapse contract.
        donor = pt.Tree.from_string("(D:1,E:2,F:3);")
        donor.rooted = True

        recipient.graft(donor, recipient.root, "A", frac_from_parent=0.5)

        assert recipient.rooted is None

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

    @pytest.mark.parametrize("mode", [None, "none", "lower", "both"])
    def test_annotation_copy_choices_and_node_attachment(self, mode):
        tree = pt.Tree.from_string("(A:8,B:2);")
        tree.set_node_attribute("A", "state", "tip")
        original = tree.child_dict[tree.root]["A"]
        original.label = "support"
        original.posterior = 0.9
        original.extra_meta = {"length": [100], "length_95%_CI": [3, 9]}
        options = {} if mode is None else {"copy_attrs": mode}
        node = tree.add_node_on_branch(tree.root, "A", frac_from_parent=0.25, **options)
        upper = tree.child_dict[tree.root][node]
        lower = tree.child_dict[node]["A"]
        assert upper.length == 2 and lower.length == 6
        assert tree.nodedict["A"].state == "tip"
        assert tree.nodedict[node].state is tree.nodedict[node].extra_meta is None
        if mode == "both":
            assert upper.extra_meta == lower.extra_meta == original.extra_meta
            assert upper.extra_meta is not lower.extra_meta
            upper.extra_meta["length"][0] = 7
            assert lower.extra_meta["length"] == original.extra_meta["length"] == [100]
        else:
            assert upper.label == "" and upper.posterior is upper.extra_meta is None
        if mode == "lower":
            assert lower.label == "support" and lower.posterior == 0.9
            assert lower.extra_meta == original.extra_meta
            lower.extra_meta["length_95%_CI"][0] = 5
            assert original.extra_meta["length_95%_CI"] == [3, 9]
        elif mode in (None, "none"):
            assert lower.label == "" and lower.posterior is lower.extra_meta is None

    def test_explicit_upper_branch_overrides_default_without_copying(self):
        tree = pt.Tree.from_string("(A:8,B:2);")
        supplied = pt.Branchstruct(length=99, label="intentional")
        supplied.extra_meta = {"keep": [1]}
        node = tree.add_node_on_branch(tree.root, "A", frac_from_parent=0.25,
                                       upper_branchstruct=supplied)
        assert tree.child_dict[tree.root][node] is supplied
        assert supplied.length == 2 and supplied.extra_meta == {"keep": [1]}
        assert tree.child_dict[node]["A"].extra_meta is None
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
        assert qa.min_value == math.inf
        assert qa.max_value == -math.inf

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

    def test_add_tracks_observed_minimum_and_maximum(self):
        qa = pt.QuantileAccumulator()
        for x in [3.0, 1.0, 4.0, 2.0]:
            qa.add(x)

        assert qa.min_value == 1.0
        assert qa.max_value == 4.0

    def test_add_nonfinite_value_invalidates_observed_bounds(self):
        qa = pt.QuantileAccumulator()
        qa.add(2.0)
        qa.add(float("nan"))

        assert qa.min_value == -math.inf
        assert qa.max_value == math.inf


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

    def test_merge_combines_observed_bounds(self):
        a = pt.QuantileAccumulator(k=7)
        b = pt.QuantileAccumulator(k=7)

        for x in [1.0, 3.0]:
            a.add(x)
        for x in [-2.0, 5.0]:
            b.add(x)

        a.merge(b)

        assert a.min_value == -2.0
        assert a.max_value == 5.0

    def test_merge_preserves_nonfinite_bounds(self):
        a = pt.QuantileAccumulator(k=7)
        b = pt.QuantileAccumulator(k=7)
        a.add(2.0)
        b.add(float("nan"))

        a.merge(b)

        assert a.min_value == -math.inf
        assert a.max_value == math.inf
        assert not a.is_effectively_constant()

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

class Test_is_effectively_constant_QuantileAccumulator:

    def test_empty_accumulator_is_not_constant(self):
        qa = pt.QuantileAccumulator()
        assert not qa.is_effectively_constant()

    def test_identical_observations_are_constant(self):
        qa = pt.QuantileAccumulator()
        for _ in range(10):
            qa.add(2.5)

        assert qa.is_effectively_constant()

    def test_numerical_noise_is_constant(self):
        qa = pt.QuantileAccumulator()
        for x in [2.5, 2.5 + 1e-13, 2.5 - 1e-13]:
            qa.add(x)

        assert qa.is_effectively_constant()

    def test_meaningful_variation_is_not_constant(self):
        qa = pt.QuantileAccumulator()
        for x in [2.5, 2.5001]:
            qa.add(x)

        assert not qa.is_effectively_constant()

    def test_nonfinite_observation_is_not_constant(self):
        qa = pt.QuantileAccumulator()
        qa.add(2.5)
        qa.add(float("nan"))

        assert not qa.is_effectively_constant()


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


###################################################################################################

class Test_set_biplen_on_existing_tree_TreePostProcessor:

    def test_final_root_edge_writes_invalidate_warmed_distances(self):
        sample = pt.Tree.from_string("(A:3,B:7);")
        summary = pt.TreeSummary(trackbips=True, trackblen=True,
                                 trackroot=True, trackrootblen=True)
        summary.add_tree(sample)
        target = pt.Tree.from_string("(A:1,B:1);")
        parents = target.parent_dict
        assert target.nodedist("A", "B") == 2
        assert target.rootdist["A"] == target.rootdist["B"] == 1
        target.build_dist_dict()

        result = pt.TreePostProcessor(summary).set_biplen_on_existing_tree(target)

        assert result is target
        assert target.child_dict[target.root]["A"].length == pytest.approx(3)
        assert target.child_dict[target.root]["B"].length == pytest.approx(7)
        assert target.nodedist("A", "B") == pytest.approx(10)
        assert target.rootdist["A"] == pytest.approx(3)
        assert target.rootdist["B"] == pytest.approx(7)
        target.build_dist_dict()
        assert target.dist_dict["A"]["B"] == pytest.approx(10)
        assert target.parent_dict is parents


class Test_update_TreeSummary:

    def test_merge_refreshes_warmed_root_frequencies_and_order(self):
        rooted_at_a = pt.Tree.from_string("(A:1,(B:1,C:1):1);")
        rooted_at_b = pt.Tree.from_string("(B:1,(A:1,C:1):1);")
        root_a = rooted_at_a.rootbip()
        root_b = rooted_at_b.rootbip()
        recipient = pt.TreeSummary(trackroot=True)
        recipient.add_tree(rooted_at_a)
        assert recipient.rootbipsummary[root_a].freq == 1
        assert [bip for _, bip, _ in recipient.sorted_rootbips] == [root_a]
        donor = pt.TreeSummary(trackroot=True)
        donor.add_tree(rooted_at_a)
        for _ in range(3):
            donor.add_tree(rooted_at_b)

        recipient.update(donor)

        assert recipient.tree_count == 5
        assert recipient.rootbipsummary[root_a].n == 2
        assert recipient.rootbipsummary[root_a].freq == pytest.approx(0.4)
        assert recipient.rootbipsummary[root_b].n == 3
        assert recipient.rootbipsummary[root_b].freq == pytest.approx(0.6)
        assert [n for n, _, _ in recipient.sorted_rootbips] == [3, 2]
        assert [bip for _, bip, _ in recipient.sorted_rootbips] == [root_b, root_a]

    def test_merge_keeps_matching_and_transfers_unseen_records_across_all_summaries(self):
        flags = dict(trackbips=True, trackclades=True, trackroot=True, trackblen=True,
                     trackheight=True, trackrootblen=True, tracktopo=True,
                     track_subcladepairs=True, store_trees=True,
                     trackci=True, ci_probs=(0.8, 0.95))
        a = pt.Tree.from_string("((A:1,B:2):3,(C:4,D:5):6);")
        a_again = pt.Tree.from_string("((A:1,B:2):3,(C:4,D:5):6);")
        b = pt.Tree.from_string("((A:2,C:3):4,(B:5,D:6):7);")
        recipient = pt.TreeSummary(**flags)
        recipient.add_tree(a)
        donor = pt.TreeSummary(**flags)
        donor.add_tree(a_again)
        donor.add_tree(b)
        serial = pt.TreeSummary(**flags)
        for sample in (a, a_again, b):
            serial.add_tree(sample)

        families = (
            ("_bipartsummary", "bipartsummary"),
            ("_cladesummary", "cladesummary"),
            ("_rootbip_summary", "rootbipsummary"),
            ("_biptoposummary", "biptoposummary"),
            ("_cladetoposummary", "cladetoposummary"),
        )
        # Warm lazy results, then capture references BEFORE consuming donor.
        for _, public in families:
            getattr(recipient, public)
            getattr(donor, public)
        kept = {}
        adopted = {}
        histograms = {}
        interval_maps = {}
        child_pairs = {}
        for private, _ in families:
            recipient_records = getattr(recipient, private)
            donor_records = getattr(donor, private)
            common = recipient_records.keys() & donor_records.keys()
            unseen = donor_records.keys() - recipient_records.keys()
            assert common and unseen
            kept[private] = {key: recipient_records[key] for key in common}
            adopted[private] = {key: donor_records[key] for key in unseen}
            if private in ("_bipartsummary", "_cladesummary"):
                for key, record in {**kept[private], **adopted[private]}.items():
                    histograms[private, key] = record.quantiles
                for key, record in adopted[private].items():
                    interval_maps[private, key] = record.extra_meta
            if private == "_cladesummary":
                for key, record in {**kept[private], **adopted[private]}.items():
                    child_pairs[key] = record.subcladepairs
        donor_count = len(donor)

        recipient.update(donor)

        assert len(recipient) == 3
        assert len(donor) == donor_count == 2  # Permitted scalar progress read only.
        for private, _ in families:
            records = getattr(recipient, private)
            for key, reference in kept[private].items():
                assert records[key] is reference
            for key, reference in adopted[private].items():
                assert records[key] is reference
        for (private, key), histogram in histograms.items():
            assert getattr(recipient, private)[key].quantiles is histogram
        # Check transferred map identity before recipient's lazy finalization can replace it.
        for (private, key), interval_map in interval_maps.items():
            assert getattr(recipient, private)[key].extra_meta is interval_map
        for key, pair_set in child_pairs.items():
            assert recipient._cladesummary[key].subcladepairs is pair_set

        def check_recipient_against_serial():
            for _, public in families:
                actual = getattr(recipient, public)
                expected = getattr(serial, public)
                assert actual.keys() == expected.keys()
                for key, result in actual.items():
                    reference = expected[key]
                    assert result.n == reference.n
                    if public in ("bipartsummary", "cladesummary"):
                        assert result.mean == pytest.approx(reference.mean)
                        assert result.M2 == pytest.approx(reference.M2)
                        assert result.freq == pytest.approx(reference.freq)
                        assert result.quantiles.counts == reference.quantiles.counts
                        assert result.extra_meta.keys() == reference.extra_meta.keys()
                        for annotation in result.extra_meta:
                            assert result.extra_meta[annotation] == pytest.approx(reference.extra_meta[annotation])
                        if public == "cladesummary":
                            assert result.subcladepairs == reference.subcladepairs
                    elif public == "rootbipsummary":
                        assert result.freq == pytest.approx(reference.freq)
                        assert result.mean_frac_to_canon() == pytest.approx(reference.mean_frac_to_canon())
                    else:
                        assert result.posterior == pytest.approx(reference.posterior)

        check_recipient_against_serial()
        root_clade = pt.Clade.from_leafset(a.leaves, a)
        assert len(recipient.cladesummary[root_clade].subcladepairs) == 2
        a_bip = pt.Bipartition.from_leafset({"A"}, a)
        assert recipient.bipartsummary[a_bip].length == pytest.approx(4 / 3)

        # A later chunk combines into the already adopted records, without donor reuse.
        tail = pt.TreeSummary(**flags)
        tail.add_tree(b)
        tail_count = len(tail)
        recipient.update(tail)
        serial.add_tree(b)
        assert len(recipient) == 4
        assert len(tail) == tail_count == 1
        assert len(donor) == donor_count
        for private, _ in families:
            for key, reference in adopted[private].items():
                assert getattr(recipient, private)[key] is reference
        for (private, key), histogram in histograms.items():
            assert getattr(recipient, private)[key].quantiles is histogram
        check_recipient_against_serial()


class Test_cladesummary_TreeSummary:

    def test_single_observation_preserves_undefined_sd(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        treesummary = pt.TreeSummary(
            trackclades=True,
            trackheight=True,
            trackci=True,
            ci_probs=(0.95,),
        )
        treesummary.add_tree(tree)

        leaf_clade = pt.Clade.from_leafset({"A"}, tree)
        nd = treesummary.cladesummary[leaf_clade]

        assert nd.height_sd is None
        assert nd.height_median == nd.height
        assert nd.extra_meta["height_95%_CI"] == (nd.height, nd.height)

    def test_effectively_constant_leaf_height_uses_exact_mean(self):
        tree1 = pt.Tree.from_string("(A:1,B:2);")
        tree2 = pt.Tree.from_string("(A:1,B:2.0000000000005);")
        treesummary = pt.TreeSummary(
            trackclades=True,
            trackheight=True,
            trackci=True,
            ci_probs=(0.95,),
        )
        treesummary.add_tree(tree1)
        treesummary.add_tree(tree2)

        leaf_clade = pt.Clade.from_leafset({"A"}, tree1)
        nd = treesummary.cladesummary[leaf_clade]

        assert nd.height_sd == 0.0
        assert nd.height_median == nd.height
        assert nd.extra_meta["height_95%_CI"] == (nd.height, nd.height)

    def test_variable_leaf_height_keeps_quantile_estimates(self):
        tree1 = pt.Tree.from_string("(A:1,B:2);")
        tree2 = pt.Tree.from_string("(A:1,B:3);")
        treesummary = pt.TreeSummary(
            trackclades=True,
            trackheight=True,
            trackci=True,
            ci_probs=(0.95,),
        )
        treesummary.add_tree(tree1)
        treesummary.add_tree(tree2)

        leaf_clade = pt.Clade.from_leafset({"A"}, tree1)
        nd = treesummary.cladesummary[leaf_clade]

        assert nd.height_sd > 0.0
        assert nd.extra_meta["height_95%_CI"][0] != nd.extra_meta["height_95%_CI"][1]


###################################################################################################

class Test_write_into_CAHeightEstimator:

    def test_preserves_unwritten_extra_metadata(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.set_node_attribute("A", "extra_meta", {
            "height_95%_CI": (-1, -1), "height_80%_CI": (0, 10), "custom": ["keep"],
        })
        plan = pt.CAHeightEstimator.build_plan(tree, trackci=True, ci_probs=(0.95,))
        estimator = pt.CAHeightEstimator(plan, trackci=True)
        estimator.add_tree(pt.Tree.from_string("(A:1,B:2);"))
        estimator.write_into(tree)
        extra = tree.nodedict["A"].extra_meta
        assert extra["height_95%_CI"] == (tree.nodedict["A"].height,) * 2
        assert extra["height_80%_CI"] == (0, 10)
        assert extra["custom"] == ["keep"]

    def test_single_observation_preserves_undefined_sd(self):
        sumtree = pt.Tree.from_string("(A:1,B:2);")
        plan = pt.CAHeightEstimator.build_plan(
            sumtree,
            trackci=True,
            ci_probs=(0.95,),
        )
        estimator = pt.CAHeightEstimator(plan, trackci=True)
        estimator.add_tree(pt.Tree.from_string("(A:1,B:2);"))

        estimator.write_into(sumtree)
        nd = sumtree.nodedict["A"]

        assert nd.height_sd is None
        assert nd.height_median == nd.height
        assert nd.extra_meta["height_95%_CI"] == (nd.height, nd.height)

    def test_effectively_constant_leaf_height_uses_exact_mean(self):
        sumtree = pt.Tree.from_string("(A:1,B:2);")
        plan = pt.CAHeightEstimator.build_plan(
            sumtree,
            trackci=True,
            ci_probs=(0.95,),
        )
        estimator = pt.CAHeightEstimator(plan, trackci=True)
        estimator.add_tree(pt.Tree.from_string("(A:1,B:2);"))
        estimator.add_tree(pt.Tree.from_string("(A:1,B:2.0000000000005);"))

        estimator.write_into(sumtree)
        nd = sumtree.nodedict["A"]

        assert nd.height_sd == 0.0
        assert nd.height_median == nd.height
        assert nd.extra_meta["height_95%_CI"] == (nd.height, nd.height)

    def test_variable_leaf_height_keeps_quantile_estimates(self):
        sumtree = pt.Tree.from_string("(A:1,B:2);")
        plan = pt.CAHeightEstimator.build_plan(
            sumtree,
            trackci=True,
            ci_probs=(0.95,),
        )
        estimator = pt.CAHeightEstimator(plan, trackci=True)
        estimator.add_tree(pt.Tree.from_string("(A:1,B:2);"))
        estimator.add_tree(pt.Tree.from_string("(A:1,B:3);"))

        estimator.write_into(sumtree)
        nd = sumtree.nodedict["A"]

        assert nd.height_sd > 0.0
        assert nd.extra_meta["height_95%_CI"][0] != nd.extra_meta["height_95%_CI"][1]


###################################################################################################

class Test_merge_CAHeightEstimator:

    def test_merge_preserves_constant_leaf_height(self):
        sumtree = pt.Tree.from_string("(A:1,B:2);")
        plan = pt.CAHeightEstimator.build_plan(
            sumtree,
            trackci=True,
            ci_probs=(0.95,),
        )
        left = pt.CAHeightEstimator(plan, trackci=True)
        right = pt.CAHeightEstimator(plan, trackci=True)
        left.add_tree(pt.Tree.from_string("(A:1,B:2);"))
        right.add_tree(pt.Tree.from_string("(A:1,B:2.0000000000005);"))

        left.merge(right)
        left.write_into(sumtree)
        nd = sumtree.nodedict["A"]

        assert nd.height_sd == 0.0
        assert nd.height_median == nd.height
        assert nd.extra_meta["height_95%_CI"] == (nd.height, nd.height)


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

        if leaf_states:
            for leaf, state in leaf_states.items():
                tree.set_node_attribute(leaf, "state", state)

        if internal_states:
            for node, state in internal_states.items():
                tree.set_node_attribute(node, "state", state)

        return tree

    def _all_nodes_fit_within_optimal_set(self, fitted_tree, possible_tree):
        fitted_ndict = fitted_tree.nodedict
        possible_ndict = possible_tree.nodedict
        for node in fitted_tree.nodes:
            fit = fitted_ndict[node].fit
            optimal_set = possible_ndict[node].optimal_set
            if isinstance(optimal_set, list):
                assert isinstance(fit, tuple)
                assert len(fit) == len(optimal_set)
                for si, site_optimal in enumerate(optimal_set):
                    assert fit[si] in site_optimal
            else:
                assert fit in optimal_set

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
        tree.set_node_attribute(ab_node, "state", "X")

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

    def test_parsimony_possible_states_errors_without_any_states(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")

        with pytest.raises(pt.TreeError):
            tree.parsimony_possible_states()

    def test_parsimony_multisite_basic(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={
                "A": ("X", "A"),
                "B": ("X", "A"),
                "C": ("Y", "A"),
                "D": ("Y", "A"),
            },
        )

        possible = tree.parsimony_possible_states()
        fitted = tree.parsimony_assign_fits(fitpref=("X", "A"))
        pscore, countdict = tree.parsimony_count_changes(fitpref=("X", "A"))

        assert pscore == 1
        assert sum(countdict.values()) == pscore
        assert dict(countdict) == {("X", "Y"): 1} or dict(countdict) == {("Y", "X"): 1}
        assert isinstance(fitted.nodedict[fitted.root].fit, tuple)
        assert fitted.nodedict[fitted.root].fit == ("X", "A")
        assert possible.nodedict[possible.root].optimal_set == [{"X", "Y"}, {"A"}]

    def test_parsimony_multisite_wasambig_is_tuple(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={
                "A": ("X", "X"),
                "B": ("X", "Y"),
                "C": ("Y", "Y"),
                "D": ("Y", "Y"),
            },
        )

        fitted = tree.parsimony_assign_fits(fitpref=("X", "Y"))

        assert fitted.nodedict[fitted.root].fit == ("X", "Y")
        assert fitted.nodedict[fitted.root].wasambig == (True, False)

    def test_parsimony_multisite_fitpref_string_broadcast(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={
                "A": ("X", "X"),
                "B": ("X", "X"),
                "C": ("Y", "Z"),
                "D": ("Y", "Z"),
            },
        )

        fitted = tree.parsimony_assign_fits(fitpref="X")

        assert fitted.nodedict[fitted.root].fit == ("X", "X")
        assert fitted.nodedict[fitted.root].wasambig == (True, True)

    def test_parsimony_multisite_fitpref_tuple_with_none(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={
                "A": ("X", "A"),
                "B": ("X", "B"),
                "C": ("Y", "A"),
                "D": ("Y", "B"),
            },
        )

        possible = tree.parsimony_possible_states()
        fitted = tree.parsimony_assign_fits(fitpref=(None, "A"))

        assert fitted.nodedict[fitted.root].fit[1] == "A"
        assert fitted.nodedict[fitted.root].fit[0] in possible.nodedict[possible.root].optimal_set[0]
        self._all_nodes_fit_within_optimal_set(fitted, possible)

    def test_parsimony_multisite_singleton_tuple_matches_single_site_score(self):
        treestring = "((A:1,B:1):1,(C:1,D:1):1);"
        single_tree = self._build_tree_with_states(
            treestring,
            leaf_states={"A": "X", "B": "X", "C": "Y", "D": "Y"},
        )
        tuple_tree = self._build_tree_with_states(
            treestring,
            leaf_states={"A": ("X",), "B": ("X",), "C": ("Y",), "D": ("Y",)},
        )

        single_score, _ = single_tree.parsimony_count_changes(fitpref="X")
        tuple_score, _ = tuple_tree.parsimony_count_changes(fitpref=("X",))

        assert single_score == tuple_score == 1

    def test_parsimony_multisite_mixed_state_types_raise(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": "X", "B": ("X", "A"), "C": ("Y", "A"), "D": ("Y", "A")},
        )

        with pytest.raises(pt.TreeError, match="mix strings and tuples/lists"):
            tree.parsimony_possible_states()

    def test_parsimony_multisite_inconsistent_lengths_raise(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={"A": ("X",), "B": ("X", "A"), "C": ("Y", "A"), "D": ("Y", "A")},
        )

        with pytest.raises(pt.TreeError, match="inconsistent tuple/list lengths"):
            tree.parsimony_possible_states()

    def test_parsimony_multisite_fitpref_length_mismatch_raises(self):
        tree = self._build_tree_with_states(
            "((A:1,B:1):1,(C:1,D:1):1);",
            leaf_states={
                "A": ("X", "A"),
                "B": ("X", "A"),
                "C": ("Y", "B"),
                "D": ("Y", "B"),
            },
        )

        with pytest.raises(pt.TreeError, match="fitpref tuple/list length does not match"):
            tree.parsimony_assign_fits(fitpref=("X",))


class Test_set_branch_attribute:

    @pytest.mark.parametrize("attribute,value", [
        ("posterior", 0.9), ("length_sd", 0.2), ("extra_meta", {"length_95%_CI": (0.1, 0.9)}),
    ])
    def test_metadata_updates_retain_valid_caches(self, attribute, value):
        tree = pt.Tree.from_string("((A:1,B:2):3,(C:4,D:5):6);")
        parent = tree.parent("A")
        parents = tree.parent_dict
        bipartitions = tree.topology_bipart
        clades = tree.topology_clade
        order = tree.sorted_intnodes()
        distances = tree.rootdist
        heights = tree.nodeheightdict
        assert tree.nodedist("A", "B") == pytest.approx(3)

        tree.set_branch_attribute(parent, "A", attribute, value)

        assert tree.get_branch_attribute(parent, "A", attribute) == value
        assert tree.parent_dict is parents
        assert tree.topology_bipart is bipartitions
        assert tree.topology_clade is clades
        assert tree.sorted_intnodes() is order
        assert tree.rootdist is distances
        assert tree.nodeheightdict is heights
        assert tree.nodedist("A", "B") == pytest.approx(3)

    def test_length_updates_refresh_distances_without_clearing_topology(self):
        tree = pt.Tree.from_string("((A:1,B:2):3,(C:4,D:5):6);")
        parent = tree.parent("A")
        tree.set_branch_attribute(parent, "A", "posterior", 0.9)
        tree.set_node_attribute("A", "height", 42.0)
        parents = tree.parent_dict
        bipartitions = tree.topology_bipart
        clades = tree.topology_clade
        order = tree.sorted_intnodes()
        distances = tree.rootdist
        heights = tree.nodeheightdict
        assert distances["A"] == pytest.approx(4)
        assert heights["A"] == pytest.approx(7)
        assert tree.nodedist("A", "B") == pytest.approx(3)
        tree.build_dist_dict()

        tree.set_branch_attribute(parent, "A", "length", 8.0)

        assert tree.parent_dict is parents
        assert tree.topology_bipart is bipartitions
        assert tree.topology_clade is clades
        assert tree.sorted_intnodes() is order
        assert tree.rootdist is not distances
        assert tree.nodeheightdict is not heights
        assert tree.rootdist["A"] == pytest.approx(11)
        assert tree.nodeheightdict["A"] == pytest.approx(0)
        assert tree.nodedist("A", "B") == pytest.approx(10)
        tree.build_dist_dict()
        assert tree.dist_dict["A"]["B"] == pytest.approx(10)
        assert tree.get_branch_attribute(parent, "A", "posterior") == 0.9
        assert tree.get_node_attribute("A", "height") == 42.0

    def test_native_length_assignment_retains_stored_data(self, length_edit_tree):
        tree = length_edit_tree
        before = _length_edit_data(tree)
        _ = tree.rootdist
        tree.set_branch_attribute(tree.parent("A"), "A", "length", 4)
        _assert_length_edit_data(tree, before)
        assert tree.rootdist["A"] == 6


class TestNodedictRefactor:
    def test_fresh_tree_has_empty_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        assert isinstance(tree.nodedict, dict)
        assert len(tree.nodedict) == 0

    def test_set_node_attribute_creates_entry(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = tree.find_mrca({"A", "B"})
        tree.set_node_attribute(mrca, "height", 1.5)
        assert mrca in tree.nodedict
        assert tree.nodedict[mrca].height == 1.5
        assert "A" not in tree.nodedict

    def test_get_node_attribute_default(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        assert tree.get_node_attribute("A", "height", default=None) is None
        assert tree.get_node_attribute("A", "height", default=-1) == -1
        assert tree.get_node_attribute("A", "height") == ""

    def test_ensure_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.ensure_nodedict()
        assert set(tree.nodedict.keys()) == tree.nodes
        for nd in tree.nodedict.values():
            assert isinstance(nd, pt.Nodestruct)

    def test_ensure_nodedict_preserves_existing(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "state", "human")
        tree.ensure_nodedict()
        assert tree.nodedict["A"].state == "human"

    def test_remove_leaf_cleans_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.ensure_nodedict()
        assert "B" in tree.nodedict
        tree.remove_leaf("B")
        assert "B" not in tree.nodedict
        assert "A" in tree.nodedict

    def test_graft_updates_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "state", "human")
        parent_of_a = tree.parent("A")
        grandparent = tree.parent(parent_of_a)
        tree.graft("E", parent=grandparent, child=parent_of_a, connect_length=0.5)
        assert "E" in tree.nodedict
        assert "A" in tree.nodedict
        assert tree.nodedict["A"].state == "human"

    def test_reroot_updates_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.ensure_nodedict()
        tree.reroot("A", tree.parent("A"), dist_from_node1=0.5)
        assert tree.root in tree.nodedict

    def test_spr_nodedict_consistency(self):
        tree = pt.Tree.randtree(ntips=10, rate_sd=0.1)
        tree.ensure_nodedict()
        tree.spr()
        assert set(tree.nodedict.keys()) <= tree.nodes

    def test_split_off_children_updates_nodedict(self):
        tree = pt.Tree.from_string("(A:1,B:1,C:1,D:1);")
        tree.ensure_nodedict()
        root = tree.root
        new_node = tree.split_off_children(root, ["A", "B"], pt.Branchstruct(length=0.5))
        assert new_node in tree.nodedict

    def test_copy_treeobject_independent_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "state", "human")
        copy = tree.copy_treeobject()
        assert copy.nodedict["A"].state == "human"
        copy.nodedict["A"].state = "mink"
        assert tree.nodedict["A"].state == "human"

    def test_copy_treeobject_empty_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        copy = tree.copy_treeobject()
        assert isinstance(copy.nodedict, dict)
        assert len(copy.nodedict) == 0

    def test_subtree_copies_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "state", "human")
        tree.set_node_attribute("C", "state", "mink")
        mrca_ab = tree.find_mrca({"A", "B"})
        subtree, basal = tree.subtree(mrca_ab)
        assert "A" in subtree.nodedict
        assert subtree.nodedict["A"].state == "human"
        assert "C" not in subtree.nodedict
        subtree.nodedict["A"].state = "bat"
        assert tree.nodedict["A"].state == "human"

    def test_nodeheightdict_derived(self):
        tree = pt.Tree.from_string("((A:1,B:1):2,(C:3,D:3):0);")
        nhd = tree.nodeheightdict
        assert isinstance(nhd, dict)
        assert set(nhd.keys()) == tree.nodes

    def test_nodeheightdict_independent_of_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        nhd = tree.nodeheightdict
        assert len(nhd) > 0
        assert len(tree.nodedict) == 0

    def test_set_blens_from_heights_uses_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.ensure_nodedict()
        tree.set_node_attribute(tree.root, "height", 5.0)
        for kid in tree.children(tree.root):
            if kid in tree.intnodes:
                tree.set_node_attribute(kid, "height", 2.0)
            else:
                tree.set_node_attribute(kid, "height", 0.0)
        for parent in tree.intnodes - {tree.root}:
            for kid in tree.children(parent):
                tree.set_node_attribute(kid, "height", 0.0)
        tree.set_blens_from_heights()
        for kid in tree.children(tree.root):
            if kid in tree.intnodes:
                assert tree.child_dict[tree.root][kid].length == 3.0

    def test_set_blens_from_heights_empty_nodedict_raises(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        with pytest.raises(pt.TreeError):
            tree.set_blens_from_heights()

    def test_clear_caches_preserves_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "state", "human")
        tree.clear_caches()
        assert "A" in tree.nodedict
        assert tree.nodedict["A"].state == "human"

    def test_clear_length_caches_preserves_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "height", 99.0)
        _ = tree.nodeheightdict
        tree.clear_length_caches()
        assert tree.nodedict["A"].height == 99.0
        nhd = tree.nodeheightdict
        assert nhd["A"] != 99.0

    def test_parsimony_with_ensure_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "state", "human")
        tree.set_node_attribute("B", "state", "human")
        tree.set_node_attribute("C", "state", "mink")
        tree.set_node_attribute("D", "state", "mink")
        result = tree.parsimony_possible_states()
        for node in result.intnodes:
            assert result.nodedict[node].optimal_set is not None

    def test_parsimony_no_states_raises(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        with pytest.raises(pt.TreeError, match="There are no nodes with \\.state"):
            tree.parsimony_possible_states()

    def test_cladesummary_uses_summary_nodestruct(self):
        ts = pt.TreeSummary(trackclades=True, trackheight=True, tracktopo=True)
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        ts.add_tree(tree)
        for clade, sn in ts.cladesummary.items():
            assert isinstance(sn, pt.SummaryNodestruct)
            assert hasattr(sn, "mean")
            assert hasattr(sn, "M2")
            assert hasattr(sn, "n")

    def test_tree_nodedict_uses_nodestruct(self):
        ts = pt.TreeSummary(trackclades=True, trackheight=True, tracktopo=True)
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        ts.add_tree(tree)
        stb = pt.SummaryTreeBuilder(ts)
        sumtree = stb.max_clade_cred_tree()
        tpp = pt.TreePostProcessor(ts)
        tpp.annotate_sumtree(sumtree)
        for node, nd in sumtree.nodedict.items():
            assert isinstance(nd, pt.Nodestruct)
            assert not hasattr(nd, "mean")
            assert not hasattr(nd, "M2")

    def test_nodestruct_no_accumulator_slots(self):
        nd = pt.Nodestruct()
        assert not hasattr(nd, "mean")
        assert not hasattr(nd, "M2")
        assert not hasattr(nd, "subcladepairs")
        assert not hasattr(nd, "posterior")

    def test_summary_nodestruct_has_accumulator_slots(self):
        sn = pt.SummaryNodestruct()
        assert hasattr(sn, "mean")
        assert hasattr(sn, "M2")
        assert hasattr(sn, "n")
        assert hasattr(sn, "subcladepairs")
        assert hasattr(sn, "posterior")
        assert hasattr(sn, "freq")

    def test_consensus_pipeline_end_to_end(self):
        trees = []
        for _ in range(20):
            t = pt.Tree.randtree(ntips=8, rate_sd=0.1)
            trees.append(t)

        ts = pt.TreeSummary(trackbips=True, trackblen=True, trackclades=True)
        for t in trees:
            ts.add_tree(t)

        sumtree = pt.build_sumtree(ts, treetype="con", blen="biplen")
        assert len(sumtree.nodedict) > 0
        for node in sumtree.intnodes:
            if node in sumtree.nodedict:
                assert sumtree.nodedict[node].clade_cred is not None

    def test_cladeheight_pipeline_end_to_end(self):
        trees = []
        for _ in range(20):
            t = pt.Tree.randtree(ntips=8, rate_sd=0.1)
            trees.append(t)

        ts = pt.TreeSummary(trackclades=True, trackheight=True, tracktopo=True)
        for t in trees:
            ts.add_tree(t)

        sumtree = pt.build_sumtree(ts, treetype="mcc", blen="cladeheight")
        for node in sumtree.nodes:
            assert node in sumtree.nodedict
            assert sumtree.nodedict[node].height is not None
        for parent in sumtree.intnodes:
            p_height = sumtree.nodedict[parent].height
            for child in sumtree.children(parent):
                c_height = sumtree.nodedict[child].height
                blen = sumtree.child_dict[parent][child].length
                assert abs(blen - (p_height - c_height)) < 1e-10

    def test_meta_comments_sparse_nodedict(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.set_node_attribute("A", "height", 0.0)
        nwk = tree.newick(print_meta=True, node_attributes=("height",))
        assert isinstance(nwk, str)
        assert "A[&height=0.0]:1.0" in nwk

    def test_cluster_no_branchstruct_heights(self):
        tree = pt.Tree.randtree(ntips=10, rate_sd=0.1)
        clusters, basenodes = tree.cluster_n(3)
        assert len(clusters) == len(basenodes)
        for parent in tree.intnodes:
            for child, br in tree.child_dict[parent].items():
                assert br.parent_height is None
                assert br.kid_height is None

    def test_cluster_results_unchanged(self):
        tree = pt.Tree.from_string("(((A:1,B:1):2,(C:1,D:1):2):1,(E:3,F:3):1);")
        clusters, _ = tree.cluster_n(3)
        assert len(clusters) >= 3
        all_leaves = set()
        for cl in clusters:
            all_leaves.update(cl)
        assert all_leaves == tree.leaves


def _copy_foundation_payload():
    """Supported annotation containers with one mutable descendant in a tuple."""
    return {
        "height_95%_CI": [1.0, 2.0],
        "nested": ([{"red"}], ("fixed", pt.MISSING)),
        "frozen": frozenset({"north", "south"}),
        "missing": pt.MISSING,
    }


class Test_pickle_MISSING:

    def test_identity_survives_direct_and_struct_pickle_roundtrips(self):
        assert pickle.loads(pickle.dumps(pt.MISSING)) is pt.MISSING
        node = pt.Nodestruct()
        node.extra_meta = {"missing": pt.MISSING, "nested": (pt.MISSING,)}

        restored = pickle.loads(pickle.dumps(node))

        assert restored.extra_meta["missing"] is pt.MISSING
        assert restored.extra_meta["nested"][0] is pt.MISSING


class Test_init_Tree:

    def test_root_interpretation_and_suffix_start_absent(self):
        tree = pt.Tree()
        assert tree.rooted is None
        assert tree.root_suffix is None

    def test_generic_empty_construction_is_unspecified(self):
        assert pt.Tree().rooted is None


class Test_copy_Nodestruct:

    def test_annotations_and_native_nested_parsimony_values_are_independent(self):
        source = pt.Nodestruct(height=7.0)
        source.extra_meta = _copy_foundation_payload()
        source.primary_set = [{"red"}, {"green"}]
        source.optimal_set = {"red"}

        copied = source.copy()

        assert copied.height == 7.0
        assert copied.extra_meta == source.extra_meta
        assert copied.extra_meta is not source.extra_meta
        assert copied.extra_meta["missing"] is pt.MISSING
        assert isinstance(copied.extra_meta["nested"], tuple)
        assert isinstance(copied.extra_meta["nested"][0], list)
        assert isinstance(copied.primary_set, list)
        assert isinstance(copied.primary_set[0], set)
        assert isinstance(copied.optimal_set, set)

        copied.extra_meta["nested"][0][0].add("blue")
        copied.primary_set[0].add("blue")
        copied.optimal_set.add("green")

        assert source.extra_meta["nested"][0][0] == {"red"}
        assert source.primary_set == [{"red"}, {"green"}]
        assert source.optimal_set == {"red"}

class Test_metadata_pairs:

    @pytest.mark.parametrize("owner_kind", ["node", "branch"])
    @pytest.mark.parametrize("invalid", [pt.MISSING, [], "bad"])
    def test_invalid_selected_container_has_owner_diagnostic(self, owner_kind, invalid):
        owner = pt.Nodestruct() if owner_kind == "node" else pt.Branchstruct()
        owner.extra_meta = invalid
        assert pt._metadata_pairs(owner, (), owner_kind, None) == []
        for selection in (None, ("height" if owner_kind == "node" else "length_sd",)):
            with pytest.raises(pt.TreeError, match=owner_kind + r" extra_meta.*dictionary or None"):
                pt._metadata_pairs(owner, selection, owner_kind, None)

    def test_disabled_metadata_does_not_validate_container(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.set_node_attribute("A", "extra_meta", pt.MISSING)
        tree.child_dict[tree.root]["A"].extra_meta = []
        assert "[&" not in tree.newick(print_meta=False)
        assert "[&" not in tree.newick(print_meta=True, node_attributes=(), branch_attributes=())

    def test_absent_and_explicitly_empty_mapping_do_not_gain_entries(self):
        node = pt.Nodestruct()
        assert node.extra_meta is None
        assert pt._metadata_pairs(node, None, "node", None) == []
        assert node.extra_meta is None
        node.extra_meta = {}
        assert pt._metadata_pairs(node, None, "node", None) == []
        assert node.extra_meta == {}

    def test_selected_presence_types_and_nested_values(self):
        node = pt.Nodestruct()
        node.extra_meta = {
            "unset": None, "missing": pt.MISSING, "empty": "", "collection": (),
            "integer": 1, "float": 1.0, "false": False, "zero": 0,
            "nested": (1, 1.0, True, "1", (pt.MISSING,)),
        }
        assert pt._metadata_pairs(node, None, "node", None) == [
            "missing=?", 'empty=""', "collection={}", "integer=1", "float=1.0",
            "false=false", "zero=0", 'nested={1,1.0,true,"1",{?}}',
        ]
        assert "unset" in node.extra_meta
        assert node.extra_meta["unset"] is None

    def test_exact_selection_literal_keys_and_duplicate_selection(self):
        node = pt.Nodestruct()
        node.extra_meta = {"ci": "literal", "height_95%_CI": (1, 2),
                           "extra_meta": "literal container name"}
        assert pt._metadata_pairs(node, ("ci", "ci", "absent"), "node", None) == ['ci="literal"']
        assert pt._metadata_pairs(node, ("extra_meta",), "node", None) == [
            'extra_meta="literal container name"']
        assert pt._metadata_pairs(node, (), "node", None) == []
        node.extra_meta = None
        assert pt._metadata_pairs(node, ("extra_meta",), "node", None) == []

    def test_extra_precedence_and_none_native_fallback(self):
        branch = pt.Branchstruct(length=0.5)
        branch.extra_meta = {"length": 0.7, "mean": 12, "extra_meta": "literal"}
        assert pt._metadata_pairs(branch, ("length", "mean"), "branch", None) == [
            "length=0.7", "mean=12"]
        branch.extra_meta["length"] = None
        assert pt._metadata_pairs(branch, ("length",), "branch", None) == ["length=0.5"]
        automatic = pt._metadata_pairs(branch, None, "branch", None)
        assert "mean=12" in automatic
        assert 'extra_meta="literal"' in automatic
        assert not any(pair.startswith(("M2=", "n=", "quantiles=", "parent_height=", "kid_height="))
                       for pair in automatic)
        branch.extra_meta["mean"] = None
        assert not any(pair.startswith("mean=")
                       for pair in pt._metadata_pairs(branch, None, "branch", None))

    @pytest.mark.parametrize("value", [None, pt.MISSING, 1, (1, 2)])
    def test_selected_native_scalar_duplicate_storage_is_invalid(self, value):
        node = pt.Nodestruct(height=2.0)
        node.extra_meta = {"height": value, "safe": 3}
        assert pt._metadata_pairs(node, ("safe",), "node", None) == ["safe=3"]
        with pytest.raises(pt.TreeError, match="node annotation 'height'.*native field"):
            pt._metadata_pairs(node, ("height",), "node", None)

    @pytest.mark.parametrize("value", [True, "2", (), (1, 2)])
    def test_mapped_scalar_type_errors_are_selected_only(self, value):
        node = pt.Nodestruct(height=value)
        node.extra_meta = {"safe": 0}
        assert pt._metadata_pairs(node, ("safe",), "node", None) == ["safe=0"]
        with pytest.raises(pt.TreeError, match="node annotation 'height'.*number or MISSING"):
            pt._metadata_pairs(node, ("height",), "node", None)

    @pytest.mark.parametrize("value,expected", [
        (pt.MISSING, "?"), ((pt.MISSING, 2.5), "{?,2.5}"),
        ((pt.MISSING, pt.MISSING), "{?,?}"), ([1, 2.0], "{1,2.0}"),
    ])
    def test_recognized_interval_missingness_and_ordered_pairs(self, value, expected):
        node = pt.Nodestruct()
        node.extra_meta = {"height_95%_HPD": value}
        assert pt._metadata_pairs(node, None, "node", None) == ["height_95%_HPD=" + expected]

    @pytest.mark.parametrize("value", [(), (1,), (1, 2, 3), (None, 2), (True, 2),
                                             ((1,), 2), {1, 2}, "interval"])
    def test_recognized_interval_rejects_wrong_shapes(self, value):
        branch = pt.Branchstruct()
        branch.extra_meta = {"length_95%_CI": value}
        with pytest.raises(pt.TreeError, match="branch annotation 'length_95%_CI'.*endpoints"):
            pt._metadata_pairs(branch, ("length_95%_CI",), "branch", None)

    def test_interval_recognition_is_owner_aware_exact_and_not_a_range_check(self):
        node = pt.Nodestruct()
        node.extra_meta = {"height_095.0%_CI": (9, -1), "height_150%_HPD": (1, 2),
                           "length_95%_CI": "wrong owner stays generic",
                           "height_95%_hpd": (), "height_1e2%_CI": "generic"}
        assert pt._metadata_pairs(node, None, "node", None) == [
            "height_095.0%_CI={9,-1}", "height_150%_HPD={1,2}",
            'length_95%_CI="wrong owner stays generic"', "height_95%_hpd={}",
            'height_1e2%_CI="generic"',
        ]

    @pytest.mark.parametrize("value", [(None,), (1, (2, None)), [1, None]])
    def test_nested_none_is_rejected_without_changing_stored_values(self, value):
        node = pt.Nodestruct()
        node.extra_meta = {"invalid": value, "safe": "kept"}
        assert pt._metadata_pairs(node, ("safe",), "node", None) == ['safe="kept"']
        with pytest.raises(pt.TreeError, match="node annotation 'invalid'.*collection element.*MISSING"):
            pt._metadata_pairs(node, ("invalid",), "node", None)
        assert node.extra_meta["invalid"] is value

    def test_keys_and_strings_use_doubled_double_quotes_and_literal_backslashes(self):
        node = pt.Nodestruct()
        node.extra_meta = {
            'a "quoted" key': 'O\'Brien said "hello" [north [coast]]\\n',
            "&directive": "literal", "!color": "#0000FF", "1": "1",
            "a:b/c;(d)": "punctuation allowed in bare key",
        }
        assert pt._metadata_pairs(node, None, "node", None) == [
            '"a ""quoted"" key"="O\'Brien said ""hello"" [north [coast]]\\n"',
            '"&directive"="literal"', '!color="#0000FF"', '1="1"',
            'a:b/c;(d)="punctuation allowed in bare key"',
        ]

    @pytest.mark.parametrize("value", ["[", "]", "][", ("[", "]")])
    def test_each_selected_string_must_balance_its_own_square_brackets(self, value):
        node = pt.Nodestruct()
        node.extra_meta = {"bad": value, "safe": "[balanced]"}
        assert pt._metadata_pairs(node, ("safe",), "node", None) == ['safe="[balanced]"']
        with pytest.raises(pt.TreeError, match="node annotation 'bad'.*unbalanced square brackets"):
            pt._metadata_pairs(node, ("bad",), "node", None)

    @pytest.mark.parametrize("key", ["", 1, "bad[", "bad]", "bad]["])
    def test_selected_invalid_keys_raise_useful_errors(self, key):
        node = pt.Nodestruct()
        node.extra_meta = {key: 1, "safe": 2}
        assert pt._metadata_pairs(node, ("safe",), "node", None) == ["safe=2"]
        with pytest.raises(pt.TreeError, match="node annotation"):
            pt._metadata_pairs(node, (key,), "node", None)

    def test_full_and_explicit_precision_preserve_float_category(self):
        node = pt.Nodestruct()
        node.extra_meta = {"precise": 1.2345678901234567, "rounded": 1.0001,
                           "large_integer": 123456789012345678901234567890,
                           "negative_zero": -0.0}
        assert pt._metadata_pairs(node, ("precise",), "node", None) == [
            "precise=1.2345678901234567"]
        assert pt._metadata_pairs(node, ("precise", "rounded"), "node", 3) == [
            "precise=1.23", "rounded=1.0"]
        assert pt._metadata_pairs(node, ("large_integer", "negative_zero"), "node", 3) == [
            "large_integer=123456789012345678901234567890", "negative_zero=-0.0"]

    def test_sets_are_deterministic_even_when_str_values_would_tie(self):
        node = pt.Nodestruct()
        node.extra_meta = {"values": {1, "1", "z"}, "frozen": frozenset(("b", "a"))}
        assert pt._metadata_pairs(node, None, "node", None) == [
            'values={"1","z",1}', 'frozen={"a","b"}']


class Test_extra_meta_storage:

    @pytest.mark.parametrize("struct_type", [pt.Nodestruct, pt.Branchstruct, pt.SummaryNodestruct])
    def test_lazy_dictionary_replaces_public_ci(self, struct_type):
        owner = struct_type()
        assert owner.extra_meta is None
        assert not hasattr(owner, "ci")
        with pytest.raises(AttributeError):
            owner.ci = {}

    def test_native_setters_replace_the_map_and_do_not_route_unknown_names(self):
        tree = pt.Tree.from_string("(A:1,B:1);")
        supplied = {"nleaves": "imported value", "extra_meta": "literal", "ci": 3}
        tree.set_node_attribute("A", "extra_meta", supplied)
        assert tree.get_node_attribute("A", "extra_meta") is supplied
        with pytest.raises(AttributeError):
            tree.set_node_attribute("A", "unknown_annotation", 4)
        result = tree.newick(print_meta=True,
                             node_attributes=("nleaves", "extra_meta", "ci"),
                             branch_attributes=())
        assert 'A[&nleaves="imported value",extra_meta="literal",ci=3]' in result
        replacement = {"other": 5}
        tree.set_node_attribute("A", "extra_meta", replacement)
        assert tree.nodedict["A"].extra_meta is replacement
        assert replacement == {"other": 5}


class Test_writer_keyword_migration:

    @pytest.mark.parametrize("method", ["newick", "nexus"])
    @pytest.mark.parametrize("as_set", [False, True])
    def test_removed_interval_keyword_and_positional_slot_fail(self, method, as_set):
        tree = pt.Tree.from_string("(A:1,B:2);")
        target = pt.TreeSet() if as_set else tree
        if as_set:
            target.addtree(tree)
        writer = getattr(target, method)
        with pytest.raises(TypeError):
            writer(ci_labels=())
        old_args = [True, True, "label", 7, None, None, None]
        if method == "nexus":
            old_args = [True, True, "label", 7, False, None, None,
                        None, "#0000FF", "#000000"]
        with pytest.raises(TypeError):
            writer(*old_args, ())

    @pytest.mark.parametrize("method", ["newick", "nexus"])
    @pytest.mark.parametrize("as_set", [False, True])
    def test_only_emittable_branch_metadata_requires_lengths(self, method, as_set):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.set_branch_attribute(tree.root, "A", "extra_meta", {"rate": 0.5})
        target = pt.TreeSet() if as_set else tree
        if as_set:
            target.addtree(tree)
        writer = getattr(target, method)
        with pytest.raises(pt.TreeError, match="require an emitted branch length"):
            writer(printdist=False, print_meta=True, branch_attributes=("rate",))
        assert "rate=" not in writer(printdist=False, print_meta=False)
        assert "rate=" not in writer(printdist=False, print_meta=True, branch_attributes=())


class Test_full_precision_output:

    @pytest.mark.parametrize("method", ["newick", "nexus"])
    @pytest.mark.parametrize("as_set", [False, True])
    def test_default_none_and_explicit_precision_across_writers(self, method, as_set):
        value = 1.2345678901234567
        tree = pt.Tree.from_string(f"(A:{value!r},B:2);")
        tree.set_node_attribute("A", "extra_meta", {"value": value, "float": 1.0, "integer": 1})
        target = pt.TreeSet() if as_set else tree
        if as_set:
            target.addtree(tree)
        writer = getattr(target, method)
        options = dict(print_meta=True, node_attributes=("value", "float", "integer"),
                       branch_attributes=())
        full = writer(**options)
        assert full == writer(precision=None, **options)
        assert f'value={value!r},float=1.0,integer=1]' in full
        assert f':{value!r}' in full
        reduced = writer(precision=3, **options)
        assert 'value=1.23,float=1.0,integer=1]' in reduced
        assert ':1.23' in reduced
        assert writer(True, True, "label", 7) == writer(precision=7)


class Test_root_suffix_output:

    def test_label_root_node_and_stem_keep_their_owners_and_controls(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.root_suffix = pt.Branchstruct(length=0.5, label="RootName")
        tree.root_suffix.extra_meta = {"rate": 1.1, "length": 9.0}
        tree.set_node_attribute(tree.root, "height", 2.0)
        before_length = tree.length()
        result = tree.newick(print_meta=True, node_attributes=("height",),
                             branch_attributes=("rate", "length"))
        assert ')RootName[&height=2.0]:0.5[&rate=1.1,length=9.0];' in result
        assert tree.length() == before_length
        assert tree.nodedist("A", "B") == 3.0
        assert ')RootName:0.5;' in tree.newick()
        assert ')RootName;' in tree.newick(printdist=False)
        assert ')RootName' not in tree.newick(printlabels=False)
        tree.root_suffix.posterior = 0.975
        assert ')0.975:0.5;' in tree.newick(labelfield="posterior")

    @pytest.mark.parametrize("stem_length,printdist", [(None, True), (0.0, False)])
    def test_emitted_stem_metadata_needs_a_printed_stem(self, stem_length, printdist):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.root_suffix = pt.Branchstruct(length=stem_length, label="RootName")
        tree.root_suffix.extra_meta = {"rate": pt.MISSING, "unset": None}
        with pytest.raises(pt.TreeError, match="stem.*length"):
            tree.newick(printdist=printdist, print_meta=True, branch_attributes=("rate",))
        assert "RootName" in tree.newick(printdist=printdist, print_meta=False)
        assert "RootName" in tree.newick(printdist=printdist, print_meta=True,
                                         branch_attributes=("unset",))
        assert tree.root_suffix.extra_meta["rate"] is pt.MISSING

    def test_automatic_label_metadata_and_zero_stem_are_explicit(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.root_suffix = pt.Branchstruct(length=None, label="RootName")
        with pytest.raises(pt.TreeError, match="stem.*length"):
            tree.newick(print_meta=True)
        assert 'RootName;' in tree.newick(print_meta=True, branch_attributes=())
        tree.root_suffix.length = 0.0
        assert 'RootName:0.0[&' in tree.newick(print_meta=True)


class Test_print_rootedness:

    def test_mixed_tree_states_and_suppression_do_not_mutate_interpretation(self):
        trees = pt.TreeSet()
        for state in (True, False, None):
            tree = pt.Tree.from_string("(A:1,B:2);")
            tree.rooted = state
            trees.addtree(tree)
        output = trees.nexus()
        assert 'tree t.1 = [&R] ' in output
        assert 'tree t.2 = [&U] ' in output
        assert 'tree t.3 = (' in output
        assert '[&R]' not in trees.nexus(print_rootedness=False)
        assert '[&U]' not in trees.nexus(print_rootedness=False)
        assert [tree.rooted for tree in trees] == [True, False, None]
        for tree in trees:
            with pytest.raises(TypeError):
                tree.nexus(rooted=False)
        with pytest.raises(TypeError):
            trees.nexus(rooted=None)


class Test_annotate_sumtree_TreePostProcessor:

    def test_writes_only_produced_interval_keys_and_isolates_values(self):
        sample = pt.Tree.from_string("(A:1,B:2);")
        summary = pt.TreeSummary(trackclades=True, trackheight=True, trackbips=True,
                                 trackblen=True, trackci=True, ci_probs=(0.95,))
        summary.add_tree(sample)
        source_node = summary.cladesummary[pt.Clade.from_leafset({"A"}, sample)]
        source_branch = summary.bipartsummary[
            pt.Bipartition.from_leafset({"A"}, sample)]
        source_node.extra_meta = {"height_95%_CI": [1, 2], "source_only": "ignore"}
        source_branch.extra_meta = {"length_95%_CI": [pt.MISSING, 3]}
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.set_node_attribute("A", "extra_meta", {
            "height_95%_CI": (-1, -1), "height_80%_CI": (0, 10), "custom": "keep",
        })
        target_branch = tree.child_dict[tree.root]["A"]
        target_branch.extra_meta = {"length_95%_CI": (-1, -1),
                                    "length_80%_CI": (0, 10), "custom": "keep"}
        postprocessor = pt.TreePostProcessor(summary)
        postprocessor.annotate_sumtree(tree)
        target_branch = tree.child_dict[tree.root]["A"]
        node_extra = tree.nodedict["A"].extra_meta
        branch_extra = target_branch.extra_meta
        assert node_extra == {"height_95%_CI": [1, 2], "height_80%_CI": (0, 10),
                              "custom": "keep"}
        assert branch_extra == {"length_95%_CI": [pt.MISSING, 3],
                                "length_80%_CI": (0, 10), "custom": "keep"}
        assert node_extra["height_95%_CI"] is not source_node.extra_meta["height_95%_CI"]
        assert branch_extra["length_95%_CI"] is not source_branch.extra_meta["length_95%_CI"]
        assert branch_extra["length_95%_CI"][0] is pt.MISSING
        node_extra["height_95%_CI"][0] = 99
        assert source_node.extra_meta["height_95%_CI"] == [1, 2]

        source_node.extra_meta["height_95%_CI"] = None
        source_branch.extra_meta = None
        postprocessor.annotate_sumtree(tree)
        target_branch = tree.child_dict[tree.root]["A"]
        node_extra = tree.nodedict["A"].extra_meta
        branch_extra = target_branch.extra_meta
        assert node_extra["height_95%_CI"] == [99, 2]
        assert branch_extra["length_95%_CI"] == [pt.MISSING, 3]

        source_node.extra_meta["height_95%_CI"] = pt.MISSING
        source_branch.extra_meta = {"length_95%_CI": pt.MISSING}
        postprocessor.annotate_sumtree(tree)
        target_branch = tree.child_dict[tree.root]["A"]
        node_extra = tree.nodedict["A"].extra_meta
        branch_extra = target_branch.extra_meta
        assert node_extra["height_95%_CI"] is pt.MISSING
        assert branch_extra["length_95%_CI"] is pt.MISSING

        # No result must not allocate a mapping on an otherwise fresh target.
        for owner in summary.cladesummary.values():
            owner.extra_meta = {}
        source_branch.extra_meta = {}
        fresh = pt.Tree.from_string("(A:1,B:2);")
        postprocessor.annotate_sumtree(fresh)
        assert all(owner.extra_meta is None for owner in fresh.nodedict.values())
        assert all(owner.extra_meta is None for children in fresh.child_dict.values()
                   for owner in children.values())

    def test_explicit_whole_field_setters_still_replace_the_mapping(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.set_node_attribute("A", "extra_meta", {"old": 1})
        replacement = {"new": 2}
        tree.set_node_attribute("A", "extra_meta", replacement)
        assert tree.nodedict["A"].extra_meta is replacement
        assert "old" not in replacement

    def test_existing_tree_is_detached_before_summary_branch_writeback(self):
        summary = _ownership_output_summary()
        # The generic low-level constructor retains its existing ownership policy.
        tree = pt.Tree.from_biplist(summary.bipartsummary)
        before = dict(tree.iter_bipinfo())
        assert any(branch is summary.bipartsummary[bip] for bip, branch in before.items())
        source_values = {bip: (record.length, record.mean, record.M2, record.n,
                               record.quantiles)
                         for bip, record in summary.bipartsummary.items()}
        original_topology = tree.topology_bipart
        original_distance = tree.nodedist("A", "B")

        pt.TreePostProcessor(summary).annotate_sumtree(tree)

        outputs = _assert_owned_output_branches(tree, summary)
        assert tree.topology_bipart == original_topology
        assert tree.nodedist("A", "B") == original_distance
        for bip, record in summary.bipartsummary.items():
            assert (record.length, record.mean, record.M2, record.n,
                    record.quantiles) == source_values[bip]
        bip = next(iter(outputs))
        outputs[bip].extra_meta["custom"][0].add("output")
        assert summary.bipartsummary[bip].extra_meta["custom"] == [{"summary"}]
        summary.bipartsummary[bip].extra_meta["custom"][0].add("summary changed")
        assert outputs[bip].extra_meta["custom"] == [{"summary", "output"}]


class Test_merge_Branchstruct:

    def test_length_only_result_keeps_inputs_unchanged(self):
        first = pt.Branchstruct(length=2, label="same")
        second = pt.Branchstruct(length=3, label="same")
        first.posterior = second.posterior = 0.9
        first.extra_meta = {"length": 100, "empty": (), "missing": pt.MISSING}
        second.extra_meta = {"length_95%_CI": (1, 4)}
        first.n = 7
        first.mean = 2
        merged = first.merge(second)
        assert merged.length == 5
        assert merged.label == ""
        assert merged.posterior is merged.length_sd is merged.extra_meta is None
        assert merged.n == 0 and merged.mean == 0 and merged.quantiles is None
        assert first.label == second.label == "same"
        assert first.n == 7 and first.extra_meta["length"] == 100
        assert second.extra_meta == {"length_95%_CI": (1, 4)}

    def test_removed_compatibility_argument_fails(self):
        with pytest.raises(TypeError):
            pt.Branchstruct().merge(pt.Branchstruct(), check_compat=True)


class Test_iter_bipinfo:

    def test_root_merge_discards_annotations_only_on_combined_edge(self):
        tree = pt.Tree.from_string("((A:1,B:2):3,(C:4,D:5):6);")
        for child in tree.children(tree.root):
            branch = tree.child_dict[tree.root][child]
            branch.label = "root side"
            branch.extra_meta = {"side": [1]}
        tip_branch = tree.child_dict[tree.parent("A")]["A"]
        tip_branch.extra_meta = {"untouched": True}
        result = dict(tree.iter_bipinfo())
        merged = result[tree.rootbip()]
        assert merged.length == 9
        assert merged.label == "" and merged.extra_meta is None
        assert result[pt.Bipartition.from_leafset({"A"}, tree)] is tip_branch
        assert all(tree.child_dict[tree.root][child].extra_meta == {"side": [1]}
                   for child in tree.children(tree.root))


class Test_remove_leaf:

    def test_serial_join_uses_length_only_branch(self):
        tree = pt.Tree.from_string("((A:2,B:3):5,C:7,D:11);")
        parent = tree.parent("A")
        tree.child_dict[tree.root][parent].label = "upper"
        tree.child_dict[tree.root][parent].extra_meta = {"upper": pt.MISSING}
        tree.child_dict[parent]["B"].extra_meta = {"lower": [1]}
        untouched = tree.child_dict[tree.root]["C"]
        untouched.extra_meta = {"keep": "yes"}
        tree.set_node_attribute(parent, "state", "deleted")
        tree.set_node_attribute("B", "state", "survives")
        _ = tree.rootdist
        tree.remove_leaf("A")
        joined = tree.child_dict[tree.root]["B"]
        assert joined.length == 8
        assert joined.label == "" and joined.extra_meta is None
        assert joined.length_sd is None
        assert tree.child_dict[tree.root]["C"] is untouched
        assert tree.nodedict["B"].state == "survives"
        assert parent not in tree.nodedict
        assert tree.rootdist["B"] == 8

    @pytest.mark.parametrize("state", [True, False, None])
    @pytest.mark.parametrize("promote_root", [False, True])
    def test_pruning_preserves_state_and_only_clears_a_removed_root_suffix(self, state, promote_root):
        text = "(C:5,(A:1,B:2):3);" if promote_root else "((A:1,B:2):3,C:4,D:5);"
        tree = pt.Tree.from_string(text)
        tree.rooted = state
        old_root = tree.root
        tree.root_suffix = pt.Branchstruct(length=0.4, label="oldroot")
        suffix = tree.root_suffix
        removed_leaf = "C" if promote_root else "A"
        removed_parent = tree.parent(removed_leaf)
        tree.set_node_attribute(old_root, "state", "oldroot")
        tree.set_node_attribute(removed_leaf, "state", "removed")
        tree.set_node_attribute("B", "height", 12.0)
        survivor = tree.nodedict["B"]
        if promote_root:
            new_root = tree.find_mrca({"A", "B"})
            tree.set_node_attribute(new_root, "height", 20.0)
            promoted_metadata = tree.nodedict[new_root]
        parents = tree.parent_dict
        distances = tree.rootdist

        tree.remove_leaf(removed_leaf)

        assert tree.rooted is state
        assert tree.nodedict["B"] is survivor
        assert survivor.height == 12.0
        assert removed_leaf not in tree.nodedict
        assert removed_parent not in tree.nodedict
        assert set(tree.nodedict) <= tree.nodes
        assert tree.parent_dict is not parents
        assert tree.rootdist is not distances
        if promote_root:
            assert tree.root == new_root
            assert tree.root_suffix is None
            assert old_root not in tree.nodes
            assert tree.nodedict[new_root] is promoted_metadata
            assert tree.rootdist["B"] == pytest.approx(2)
        else:
            assert tree.root == old_root
            assert tree.root_suffix is suffix
            assert tree.rootdist["B"] == pytest.approx(5)


class Test_reroot:

    def test_subdivided_root_edge_starts_without_source_annotations(self):
        tree = pt.Tree.from_string("(A:8,B:2,C:3);")
        original = tree.child_dict[tree.root]["A"]
        original.label = "old segment"
        original.extra_meta = {"length_95%_CI": [1, 9]}
        tree.reroot(tree.root, "A", frac_from_node1=0.25)
        assert sorted(branch.length for branch in tree.child_dict[tree.root].values()) == [2, 6]
        assert all(branch.label == "" and branch.extra_meta is None
                   for branch in tree.child_dict[tree.root].values())
        assert original.extra_meta == {"length_95%_CI": [1, 9]}

    @pytest.mark.parametrize("state", [True, False, None])
    def test_selecting_current_root_establishes_rootedness_and_keeps_suffix(self, state):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        tree.rooted = state
        tree.root_suffix = pt.Branchstruct(length=0.0, label="same-root")
        suffix = tree.root_suffix
        root = tree.root
        before = tree.newick()
        tree.reroot(root, polytomy=True)
        assert tree.rooted is True
        assert tree.root == root
        assert tree.root_suffix is suffix
        assert tree.newick() == before

    @pytest.mark.parametrize("state", [True, False, None])
    def test_existing_node_reroot_clears_suffix_and_keeps_node_ownership(self, state):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        oldroot = tree.root
        newroot = tree.find_mrca({"A", "B"})
        tree.rooted = state
        tree.root_suffix = pt.Branchstruct(length=0.0, label="old-root")
        tree.set_node_attribute(oldroot, "extra_meta", {"old": ["own"]})
        tree.set_node_attribute(newroot, "extra_meta", {"new": ["own"]})
        old_node = tree.nodedict[oldroot]
        new_node = tree.nodedict[newroot]
        reversed_branch = tree.child_dict[oldroot][newroot]
        tree.reroot(newroot, polytomy=True)
        assert tree.rooted is True
        assert tree.root == newroot
        assert tree.root_suffix is None
        assert tree.nodedict[oldroot] is old_node
        assert tree.nodedict[newroot] is new_node
        assert tree.child_dict[newroot][oldroot] is reversed_branch

    def test_failed_selection_does_not_establish_interpretation(self):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        tree.rooted = False
        with pytest.raises(pt.TreeError, match="Need to specify node2"):
            tree.reroot("A")
        assert tree.rooted is False

    def test_reversed_existing_edges_keep_objects_values_and_attachment(self):
        tree = pt.Tree.from_string("(((A:1,B:2):3,C:4):5,D:6,E:7);")
        base = tree.find_mrca({"A", "B"})
        middle = tree.parent(base)
        oldroot = tree.root
        old_edges = [(oldroot, middle), (middle, base)]
        records = []
        for parent, child in old_edges:
            branch = tree.child_dict[parent][child]
            branch.label = "old descendant support"
            branch.posterior = pt.MISSING
            branch.extra_meta = {"directional": [child], "empty": (), "false": False}
            records.append((parent, child, branch, branch.extra_meta, branch.length))
        unrelated = tree.child_dict[oldroot]["E"]
        tree.set_node_attribute(base, "height", 9.0)
        tree.reroot(base, polytomy=True)
        for parent, child, branch, mapping, length in records:
            assert tree.child_dict[child][parent] is branch
            assert branch.length == length
            assert branch.label == "old descendant support"
            assert branch.posterior is pt.MISSING
            assert branch.extra_meta is mapping
            assert mapping == {"directional": [child], "empty": (), "false": False}
        assert tree.child_dict[oldroot]["E"] is unrelated
        assert tree.nodedict[base].height == 9.0


class Test_remove_branch:

    @pytest.mark.parametrize("removed_length", [0, 4])
    def test_redistributed_edges_are_fresh_and_deleted_node_data_is_removed(self, removed_length):
        tree = pt.Tree.from_string(f"((A:2,B:3):{removed_length},C:5,D:6);")
        node = tree.parent("A")
        tree.set_node_attribute(node, "extra_meta", {"deleted": [1]})
        tree.set_node_attribute("A", "state", "survives")
        for branch in tree.child_dict[node].values():
            branch.label = "old endpoints"
            branch.posterior = 0.9
            branch.extra_meta = {"length": pt.MISSING, "empty": ()}
        untouched = tree.child_dict[tree.root]["C"]
        untouched.extra_meta = {"keep": [1]}
        before = tree.length()
        old_distance = tree.nodedist("A", "B")
        _ = tree.rootdist
        tree.remove_branch(tree.root, node)
        assert tree.length() == before
        assert tree.nodedist("A", "B") == old_distance + removed_length
        assert tree.rootdist["A"] == 2 + removed_length / 2
        assert tree.child_dict[tree.root]["A"].length == 2 + removed_length / 2
        assert tree.child_dict[tree.root]["B"].length == 3 + removed_length / 2
        for child in ("A", "B"):
            branch = tree.child_dict[tree.root][child]
            assert branch.label == "" and branch.extra_meta is branch.posterior is None
            assert branch.n == 0 and branch.quantiles is None
        assert tree.child_dict[tree.root]["C"] is untouched
        assert untouched.extra_meta == {"keep": [1]}
        assert tree.nodedict["A"].state == "survives"
        assert node not in tree.nodedict and node not in tree.nodes


class Test__add_bips_TreeSummary:

    def test_fresh_accumulators_ignore_sample_statistics_and_leave_inputs_unchanged(self):
        first = pt.Tree.from_string("(A:1,B:2,C:3);")
        observed = first.child_dict[first.root]["A"]
        observed.label = "sample label"
        observed.posterior = 0.25
        observed.length_sd = 11.0
        observed.mean, observed.M2, observed.n = 99.0, 12.0, 17
        observed.extra_meta = {"sample_only": ["keep"]}
        observed.quantiles = pt.QuantileAccumulator(k=7)
        observed.quantiles.add(99.0)
        sample_quantiles = observed.quantiles
        sample_counts = observed.quantiles.counts.copy()
        assert first.nodedist("A", "B") == 3.0
        summary = pt.TreeSummary(trackbips=True, trackblen=True,
                                 trackci=True, ci_probs=(0.95,))
        bip = pt.Bipartition.from_leafset({"A"}, first)

        summary.add_tree(first)
        accumulator = summary._bipartsummary[bip]
        assert accumulator is not observed
        assert accumulator.n == 1
        assert accumulator.mean == 1.0
        assert accumulator.M2 == 0.0
        assert accumulator.label == ""
        assert accumulator.extra_meta is None
        assert accumulator.quantiles is not sample_quantiles
        assert accumulator.quantiles.n == 1
        accumulator_quantiles = accumulator.quantiles

        summary.add_tree(pt.Tree.from_string("(A:3,B:4,C:5);"))
        result = summary.bipartsummary[bip]

        assert result is accumulator
        assert result.quantiles is accumulator_quantiles
        assert result.n == result.quantiles.n == 2
        assert result.length == result.mean == pytest.approx(2.0)
        assert result.M2 == pytest.approx(2.0)
        assert result.length_sd == pytest.approx(math.sqrt(2.0))
        assert result.posterior == 1.0
        assert set(result.extra_meta) == {"length_95%_CI"}
        assert observed.length == 1.0
        assert observed.label == "sample label"
        assert observed.posterior == 0.25
        assert observed.length_sd == 11.0
        assert (observed.mean, observed.M2, observed.n) == (99.0, 12.0, 17)
        assert observed.extra_meta == {"sample_only": ["keep"]}
        assert observed.quantiles is sample_quantiles
        assert sample_quantiles.n == 1
        assert sample_quantiles.counts == sample_counts
        assert first.nodedist("A", "B") == 3.0

    def test_count_only_accumulation_does_not_inherit_sample_records(self):
        sample = pt.Tree.from_string("(A:1,B:2,C:3);")
        observed = sample.child_dict[sample.root]["A"]
        observed.n, observed.mean, observed.M2 = 17, 99.0, 12.0
        observed.label = "sample label"
        observed.extra_meta = {"sample_only": pt.MISSING}
        summary = pt.TreeSummary(trackbips=True)
        bip = pt.Bipartition.from_leafset({"A"}, sample)

        summary.add_tree(sample)
        result = summary.bipartsummary[bip]

        assert result is not observed
        assert result.n == 1
        assert result.posterior == 1.0
        assert result.mean == result.M2 == 0.0
        assert result.quantiles is None
        assert result.extra_meta is None
        assert result.label == ""
        assert (observed.n, observed.mean, observed.M2) == (17, 99.0, 12.0)
        assert observed.extra_meta["sample_only"] is pt.MISSING


def _ownership_output_summary():
    summary = pt.TreeSummary(trackbips=True, trackblen=True, tracktopo=True,
                             trackci=True, ci_probs=(0.95,))
    for text in ("((A:1,B:2):3,C:4,D:5);", "((A:3,B:4):5,C:6,D:7);"):
        summary.add_tree(pt.Tree.from_string(text))
    for record in summary.bipartsummary.values():
        record.extra_meta["custom"] = [{"summary"}]
        record.extra_meta["length_95%_CI"] = list(record.extra_meta["length_95%_CI"])
    return summary


def _assert_owned_output_branches(tree, summary):
    result = dict(tree.iter_bipinfo())
    for bip, branch in result.items():
        source = summary.bipartsummary[bip]
        assert branch is not source
        assert branch.length == source.length
        assert branch.length_sd == source.length_sd
        assert branch.length_median == source.length_median
        assert branch.posterior == source.posterior
        assert branch.freq == source.freq
        assert branch.bipartition_cred == source.bipartition_cred
        assert branch.extra_meta == source.extra_meta
        assert branch.extra_meta is not source.extra_meta
        assert branch.extra_meta["custom"] is not source.extra_meta["custom"]
        assert branch.extra_meta["length_95%_CI"] is not source.extra_meta["length_95%_CI"]
        assert branch.quantiles is None
        assert branch.mean == branch.M2 == 0.0
        assert branch.n == 0
        assert branch.parent_height is branch.kid_height is None
        assert source.quantiles is not None
        assert source.n == 2
    return result


class Test_contree_SummaryTreeBuilder:

    def test_repeated_consensus_outputs_are_independent_of_summary_and_each_other(self):
        summary = _ownership_output_summary()
        builder = pt.SummaryTreeBuilder(summary)

        first = builder.contree()
        second = builder.contree()

        first_bips = _assert_owned_output_branches(first, summary)
        second_bips = _assert_owned_output_branches(second, summary)
        bip = next(iter(first_bips))
        source = summary.bipartsummary[bip]
        first_bips[bip].length = 99.0
        first_bips[bip].extra_meta["custom"][0].add("first")
        assert second_bips[bip].length == source.length
        assert second_bips[bip].extra_meta["custom"] == [{"summary"}]
        assert source.extra_meta["custom"] == [{"summary"}]
        source.extra_meta["custom"][0].add("source")
        assert second_bips[bip].extra_meta["custom"] == [{"summary"}]

    def test_compatible_bipartition_insertion_owns_its_branch_and_metadata(self):
        summary = pt.TreeSummary(trackbips=True, trackblen=True,
                                 trackci=True, ci_probs=(0.95,))
        # The later EF candidate exercises the existing allcompat insertion path.
        # This test does not change the current candidate-selection algorithm.
        texts = (
            ["((A:1,B:1):1,C:1,D:1,(E:1,F:1):2,G:1);"] * 3
            + ["((A:1,B:1):1,C:1,D:1,E:1,F:1,G:1);"] * 3
            + ["(A:1,B:1,(C:1,D:1):3,E:1,F:1,G:1);"] * 4
        )
        first = None
        for text in texts:
            tree = pt.Tree.from_string(text)
            if first is None:
                first = tree
            summary.add_tree(tree)
        ef = pt.Bipartition.from_leafset({"E", "F"}, first)
        source = summary.bipartsummary[ef]
        source.extra_meta["custom"] = [{"summary"}]

        result = pt.SummaryTreeBuilder(summary).contree(allcompat=True)
        output = dict(result.iter_bipinfo())[ef]

        assert output is not source
        assert output.length == source.length == 2.0
        assert output.extra_meta is not source.extra_meta
        assert output.quantiles is None
        assert output.n == 0
        output.extra_meta["custom"][0].add("output")
        assert source.extra_meta["custom"] == [{"summary"}]

    @pytest.mark.parametrize("allcompat", [False, True])
    def test_consensus_interpretation_is_unrooted(self, rootedness_constructor_summary, allcompat):
        summary, representative = rootedness_constructor_summary
        before = representative.newick()
        tree = pt.SummaryTreeBuilder(summary).contree(allcompat=allcompat)
        assert tree.rooted is False
        assert representative.rooted is None
        assert representative.newick() == before


class Test_max_bipart_cred_tree_SummaryTreeBuilder:

    def test_mbc_output_owns_branches_without_accumulation_bookkeeping(self):
        summary = _ownership_output_summary()

        result = pt.SummaryTreeBuilder(summary).max_bipart_cred_tree()

        outputs = _assert_owned_output_branches(result, summary)
        bip = next(iter(outputs))
        outputs[bip].extra_meta["length_95%_CI"][0] = -99.0
        assert summary.bipartsummary[bip].extra_meta["length_95%_CI"][0] != -99.0

    def test_bipartition_summary_is_unrooted(self, rootedness_constructor_summary):
        summary, representative = rootedness_constructor_summary
        assert pt.SummaryTreeBuilder(summary).max_bipart_cred_tree().rooted is False
        assert representative.rooted is None


class Test_move_root:

    @pytest.mark.parametrize("state", [True, False, None])
    def test_mechanical_reorientation_does_not_establish_rootedness(self, state):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        tree.rooted = state
        newroot = tree.find_mrca({"A", "B"})
        tree.root_suffix = pt.Branchstruct(length=0.0, label="old-root")
        tree._move_root(newroot, polytomy=True)
        assert tree.rooted is (None if state is True else state)
        assert tree.root_suffix is None

    @pytest.mark.parametrize("state", [True, False, None])
    def test_mechanical_noop_preserves_state_and_suffix(self, state):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        tree.rooted = state
        tree.root_suffix = pt.Branchstruct(length=0.0, label="same-root")
        suffix = tree.root_suffix
        tree._move_root(tree.root, polytomy=True)
        assert tree.rooted is state
        assert tree.root_suffix is suffix


class Test_match_nodes:

    @pytest.mark.parametrize("state", [True, False, None])
    def test_comparison_uses_mechanical_route_and_preserves_sources(self, monkeypatch, state):
        first = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        second = pt.Tree.from_string("(A:1,(B:1,(C:1,D:1):2):1);")
        first.rooted = second.rooted = state
        first.root_suffix = pt.Branchstruct(length=0.0, label="first")
        second.root_suffix = pt.Branchstruct(length=0.0, label="second")
        first_suffix, second_suffix = first.root_suffix, second.root_suffix
        before = first.newick(), second.newick()

        def reject_deliberate_root(*args, **kwargs):
            raise AssertionError("comparison must use mechanical root movement")

        monkeypatch.setattr(pt.Tree, "reroot", reject_deliberate_root)
        mapping, _, _ = first.match_nodes(second)
        assert all(mapping[leaf] == leaf for leaf in first.leaves)
        assert (first.newick(), second.newick()) == before
        assert first.rooted is second.rooted is state
        assert first.root_suffix is first_suffix
        assert second.root_suffix is second_suffix


class Test_collapse_bifurcating_root:

    @pytest.mark.parametrize("state", [True, False, None])
    def test_actual_removal_clears_suffix_and_only_invalidates_meaningful_root(self, state):
        tree = pt.Tree.from_string("(A:1,(B:2,C:3):4);")
        tree.rooted = state
        oldroot = tree.root
        newroot = tree.find_mrca({"B", "C"})
        tree.root_suffix = pt.Branchstruct(length=0.0, label="old-root")
        tree.set_node_attribute(oldroot, "state", "deleted")
        tree.set_node_attribute(newroot, "state", "retained")
        survivor = tree.nodedict[newroot]
        tree.collapse_bifurcating_root()
        assert tree.rooted is (None if state is True else state)
        assert tree.root_suffix is None
        assert oldroot not in tree.nodedict
        assert tree.nodedict[newroot] is survivor

    @pytest.mark.parametrize("state", [True, False, None])
    def test_multifurcation_noop_preserves_state_and_suffix(self, state):
        tree = pt.Tree.from_string("(A:1,B:2,C:3);")
        tree.rooted = state
        tree.root_suffix = pt.Branchstruct(length=0.0, label="same-root")
        suffix, root = tree.root_suffix, tree.root
        tree.collapse_bifurcating_root()
        assert tree.rooted is state
        assert tree.root == root
        assert tree.root_suffix is suffix

    def test_two_leaf_structural_collapse_still_raises(self):
        tree = pt.Tree.from_string("(A:1,B:2);")
        tree.rooted = True
        with pytest.raises(pt.TreeError, match="only leaves"):
            tree.collapse_bifurcating_root()
        assert tree.rooted is True


class Test_unroot:

    @pytest.mark.parametrize("state", [True, False, None])
    @pytest.mark.parametrize("sample", ["(A:1,(B:2,C:3):4);", "(A:1,B:2,C:3);", "(A:1,B:2);"])
    def test_deliberate_unrooting_establishes_false_and_obeys_suffix_ownership(self, state, sample):
        tree = pt.Tree.from_string(sample)
        tree.rooted = state
        tree.root_suffix = pt.Branchstruct(length=0.0, label="old-root")
        suffix, oldroot = tree.root_suffix, tree.root
        leaves = tree.leaves.copy()
        tree.unroot()
        assert tree.rooted is False
        assert tree.leaves == leaves
        if tree.root == oldroot:
            assert tree.root_suffix is suffix
        else:
            assert tree.root_suffix is None


class Test_deroot:

    def test_deprecated_alias_keeps_mechanical_semantics(self):
        tree = pt.Tree.from_string("(A:1,(B:2,C:3):4);")
        tree.rooted = True
        tree.deroot()
        assert tree.rooted is None


class Test_rootminvar:

    @pytest.mark.parametrize("state", [False, None])
    def test_already_optimal_root_establishes_true_without_replacement(self, state):
        tree = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        tree.rooted = state
        tree.root_suffix = pt.Branchstruct(length=0.0, label="optimal")
        suffix, root = tree.root_suffix, tree.root
        before = tree.newick()
        tree.rootminvar()
        assert tree.rooted is True
        assert tree.root == root
        assert tree.root_suffix is suffix
        assert tree.newick() == before


class Test_rootmid:

    def test_successful_midpoint_selection_establishes_true(self, treedata):
        tree = pt.Tree.from_string(treedata["simplestring"])
        tree.rooted = False
        tree.rootmid()
        assert tree.rooted is True


class Test_rootout:

    def test_successful_outgroup_selection_establishes_true(self):
        tree = pt.Tree.from_string("((A:1,B:1):1,C:2,D:2);")
        tree.rooted = None
        tree.rootout("A")
        assert tree.rooted is True


class Test_root_maxfreq:

    @pytest.mark.parametrize("already_matching", [True, False])
    def test_success_establishes_true_even_when_current_root_matches(self, already_matching):
        sample = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        summary = pt.TreeSummary(trackroot=True)
        summary.add_tree(sample)
        target = sample.copy_treeobject() if already_matching else pt.Tree.from_string(
            "((A:1,B:1):1,C:2,D:2);"
        )
        target.rooted = False
        target.root_suffix = pt.Branchstruct(length=0.0, label="root-label")
        suffix, oldroot = target.root_suffix, target.root

        result = pt.TreePostProcessor(summary).root_maxfreq(target)

        assert result is target
        assert target.rooted is True
        assert target.rootbip() == sample.rootbip()
        if already_matching:
            assert target.root == oldroot
            assert target.root_suffix is suffix
        else:
            assert target.root_suffix is None
        assert sample.rooted is None

    def test_selects_a_different_observed_root_with_current_bipartition_api(self):
        sample = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
        summary = pt.TreeSummary(trackroot=True)
        summary.add_tree(sample)
        target = pt.Tree.from_string("((A:1,B:1):1,C:2,D:2);")
        result = pt.TreePostProcessor(summary).root_maxfreq(target)
        assert result is target
        assert target.rootbip() == sample.rootbip()
        assert target.rooted is True


class Test_subtree:

    @pytest.mark.parametrize("state", [True, False, None])
    @pytest.mark.parametrize("at_source_root", [False, True])
    def test_returned_tree_inherits_interpretation(self, state, at_source_root):
        source = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        source.rooted = state
        node = source.root if at_source_root else source.find_mrca({"A", "B"})
        before = source.newick()

        extracted, basal = source.subtree(node)

        assert isinstance(extracted, pt.Tree)
        assert extracted.rooted is state
        assert extracted.root == node
        assert source.rooted is state
        assert source.newick() == before
        assert (basal is None) is at_source_root


    @pytest.mark.parametrize("state", [True, False, None])
    def test_leaf_result_remains_name_and_independent_basal_branch(self, state):
        source = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        source.rooted = state
        original_branch = source.child_dict[source.parent("A")]["A"]
        original_branch.label = "leaf-edge"
        original_branch.extra_meta = {"colors": ["red", "blue"]}

        extracted, basal = source.subtree("A")

        assert extracted == "A"
        assert isinstance(extracted, str)
        assert basal is not original_branch
        assert basal.length == 1.0
        assert basal.label == "leaf-edge"
        basal.extra_meta["colors"].append("green")
        assert original_branch.extra_meta == {"colors": ["red", "blue"]}
        assert source.rooted is state

    def test_proper_extraction_keeps_basal_data_separate_and_copies_metadata(self):
        source = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        node = source.find_mrca({"A", "B"})
        source.root_suffix = pt.Branchstruct(length=0.0, label="source-root")
        source.root_suffix.extra_meta = {"stem": ["old"]}
        source.set_node_attribute(node, "extra_meta", {"node": ["retained"]})
        incoming = source.child_dict[source.root][node]
        incoming.label = "basal-label"
        incoming.extra_meta = {"basal": ["kept"]}
        inner = source.child_dict[node]["A"]
        inner.extra_meta = {"inner": ["kept"]}

        extracted, basal = source.subtree(node)

        assert extracted.root_suffix is None
        assert extracted.nodedict[node].extra_meta == {"node": ["retained"]}
        assert basal.label == "basal-label"
        assert basal.length == 3.0
        assert basal.extra_meta == {"basal": ["kept"]}
        extracted.nodedict[node].extra_meta["node"].append("new")
        extracted.child_dict[node]["A"].extra_meta["inner"].append("new")
        basal.extra_meta["basal"].append("new")
        assert source.nodedict[node].extra_meta == {"node": ["retained"]}
        assert inner.extra_meta == {"inner": ["kept"]}
        assert incoming.extra_meta == {"basal": ["kept"]}
        assert source.root_suffix.extra_meta == {"stem": ["old"]}

    @pytest.mark.parametrize("stem_length", [None, 0.0, 0.75])
    def test_source_root_copies_suffix_without_inventing_basal_branch(self, stem_length):
        source = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        source.root_suffix = pt.Branchstruct(length=stem_length, label="RootName")
        source.root_suffix.extra_meta = {"rate": [1.2, pt.MISSING]}
        source.set_node_attribute(source.root, "extra_meta", {"root_node": ["own"]})

        extracted, basal = source.subtree(source.root)

        assert basal is None
        assert extracted.root_suffix is not source.root_suffix
        assert extracted.root_suffix.length == stem_length
        assert extracted.root_suffix.label == "RootName"
        assert extracted.root_suffix.extra_meta["rate"][1] is pt.MISSING
        extracted.root_suffix.extra_meta["rate"][0] = 9.0
        extracted.nodedict[extracted.root].extra_meta["root_node"].append("new")
        assert source.root_suffix.extra_meta["rate"] == [1.2, pt.MISSING]
        assert source.nodedict[source.root].extra_meta == {"root_node": ["own"]}

    @pytest.mark.parametrize("empty_container", [False, True])
    def test_source_root_does_not_allocate_empty_suffix(self, empty_container):
        source = pt.Tree.from_string("(A:1,B:2,C:3);")
        if empty_container:
            source.root_suffix = pt.Branchstruct(length=None)
            source.root_suffix.extra_meta = {}
        extracted, basal = source.subtree(source.root)
        assert extracted.root_suffix is None
        assert basal is None


class Test_prune_subtree:

    @pytest.mark.parametrize("state", [True, False, None])
    def test_extracted_tree_inherits_state_before_source_mutation(self, state):
        source = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        source.rooted = state
        node = source.find_mrca({"A", "B"})
        source.set_node_attribute(node, "extra_meta", {"node": ["own"]})
        source.root_suffix = pt.Branchstruct(length=0.0, label="retained-root")
        suffix = source.root_suffix

        extracted, basal = source.prune_subtree(node)

        assert extracted.rooted is state
        assert extracted.leaves == {"A", "B"}
        assert extracted.root_suffix is None
        assert basal.length == 3.0
        assert source.rooted is state
        assert source.leaves == {"C", "D"}
        assert source.root_suffix is suffix
        source.rooted = (state is not True)
        assert extracted.rooted is state
        extracted.nodedict[node].extra_meta["node"].append("new")
        assert node not in source.nodedict

    @pytest.mark.parametrize("state", [True, False, None])
    @pytest.mark.parametrize("promote_root", [False, True])
    def test_remaining_tree_keeps_state_and_suffix_follows_root_survival(self, state, promote_root):
        text = "((A:1,B:2):3,(C:4,D:5):6);" if promote_root else "(X:7,(A:1,B:2):3,Y:8,Z:9);"
        tree = pt.Tree.from_string(text)
        tree.rooted = state
        old_root = tree.root
        basenode = tree.find_mrca({"A", "B"})
        tree.root_suffix = pt.Branchstruct(length=0.4, label="oldroot")
        suffix = tree.root_suffix
        tree.set_node_attribute(old_root, "state", "root")
        tree.set_node_attribute(basenode, "height", 10.0)
        survivor = "C" if promote_root else "X"
        tree.set_node_attribute(survivor, "state", "retained")
        metadata = tree.nodedict[survivor]

        extracted, basal = tree.prune_subtree(basenode)

        assert tree.rooted is state
        assert isinstance(extracted, pt.Tree)
        assert extracted.leaves == {"A", "B"}
        assert basal.length == pytest.approx(3)
        assert basenode not in tree.nodes
        assert basenode not in tree.nodedict
        assert tree.nodedict[survivor] is metadata
        assert set(tree.nodedict) <= tree.nodes
        if promote_root:
            assert tree.root != old_root
            assert old_root not in tree.nodedict
            assert tree.root_suffix is None
        else:
            assert tree.root == old_root
            assert tree.root_suffix is suffix


class Test_remove_leaves:
    @pytest.mark.parametrize("state", [True, False, None])
    def test_batch_removals_preserve_interpretation_through_root_promotion(self, state):
        tree = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        tree.rooted = state
        old_root = tree.root
        new_root = tree.find_mrca({"A", "B"})
        tree.root_suffix = pt.Branchstruct(length=0.4, label="oldroot")
        tree.set_node_attribute("A", "state", "retained")
        retained = tree.nodedict["A"]

        tree.remove_leaves(["C", "D"])

        assert tree.rooted is state
        assert tree.root == new_root
        assert old_root not in tree.nodes
        assert tree.root_suffix is None
        assert tree.leaves == {"A", "B"}
        assert tree.nodedict["A"] is retained


class Test_prune_maxlen:
    @pytest.mark.parametrize("state", [True, False, None])
    @pytest.mark.parametrize("select_only", [False, True])
    def test_wrapper_preserves_state_and_selection_only_preserves_suffix(self, state, select_only):
        tree = pt.Tree.from_string("(C:5,(A:1,B:2):3);")
        tree.rooted = state
        old_root = tree.root
        tree.root_suffix = pt.Branchstruct(length=0.4, label="oldroot")
        suffix = tree.root_suffix

        result = tree.prune_maxlen(nkeep=2, keeplist=["A", "B"], return_leaves=select_only)

        assert tree.rooted is state
        if select_only:
            assert result == {"A", "B"}
            assert tree.root == old_root
            assert tree.root_suffix is suffix
            assert tree.leaves == {"A", "B", "C"}
        else:
            assert result is None
            assert tree.root != old_root
            assert tree.root_suffix is None
            assert tree.leaves == {"A", "B"}


class Test_nj:

    def test_neighbor_joining_establishes_unrooted_state(self):
        matrix = pt.Distmatrix.from_numpy_array(
            np.array([[0., 2., 4., 4.], [2., 0., 4., 4.],
                      [4., 4., 0., 2.], [4., 4., 2., 0.]]),
            ["A", "B", "C", "D"],
        )
        assert matrix.nj().rooted is False


class Test_upgma:

    def test_upgma_establishes_rooted_hierarchy(self):
        matrix = pt.Distmatrix.from_numpy_array(
            np.array([[0., 2., 4., 4.], [2., 0., 4., 4.],
                      [4., 4., 0., 2.], [4., 4., 2., 0.]]),
            ["A", "B", "C", "D"],
        )
        assert matrix.upgma().rooted is True


@pytest.fixture
def rootedness_constructor_summary():
    summary = pt.TreeSummary(
        trackbips=True, trackclades=True, trackblen=True, trackheight=True,
        tracktopo=True, track_subcladepairs=True, store_trees=True,
    )
    sample = pt.Tree.from_string("((A:1,B:1):1,(C:1,D:1):1);")
    sample.rooted = None
    summary.add_tree(sample)
    return summary, sample


class Test_max_clade_cred_tree_SummaryTreeBuilder:

    def test_clade_summary_is_rooted(self, rootedness_constructor_summary):
        summary, representative = rootedness_constructor_summary
        assert pt.SummaryTreeBuilder(summary).max_clade_cred_tree().rooted is True
        assert representative.rooted is None


class Test_hipstr_tree_SummaryTreeBuilder:

    @pytest.mark.parametrize("majrule", [False, True])
    def test_ancestral_summary_is_rooted(self, rootedness_constructor_summary, majrule):
        summary, representative = rootedness_constructor_summary
        assert pt.SummaryTreeBuilder(summary).hipstr_tree(majrule=majrule).rooted is True
        assert representative.rooted is None


class Test_add_tree_TreeSummary:

    def test_first_representative_is_retained_by_reference_without_data_mutation(self):
        first = pt.Tree.from_string("((A:1,B:2):3,C:4,D:5);")
        first.rooted = True
        first.set_node_attribute("A", "extra_meta", {"sample_only": ["keep"]})
        snapshot = first.newick(print_meta=True)
        summary = pt.TreeSummary(trackbips=True, trackclades=True, trackblen=True,
                                 trackheight=True, tracktopo=True, store_trees=True)

        summary.add_tree(first)
        summary.add_tree(pt.Tree.from_string("((A:3,B:4):5,C:6,D:7);"))
        # Finalizing computed statistics must not rewrite the retained sample.
        _ = summary.bipartsummary, summary.cladesummary

        assert next(iter(summary.biptoposummary.values())).tree is first
        assert next(iter(summary.cladetoposummary.values())).tree is first
        assert first.rooted is True
        assert first.newick(print_meta=True) == snapshot
        assert first.nodedict["A"].extra_meta == {"sample_only": ["keep"]}


@pytest.fixture
def length_edit_tree():
    tree = pt.Tree.from_string("((A:1,B:1):2,C:3,D:3);")
    tree.rooted = True
    tree.root_suffix = pt.Branchstruct(length=0.5, label="root")
    tree.root_suffix.extra_meta = {"stem": [pt.MISSING]}
    heights = tree.nodeheightdict.copy()
    for node in tree.nodes:
        tree.set_node_attribute(node, "height", heights[node])
        tree.set_node_attribute(node, "state", "")
        tree.set_node_attribute(node, "extra_meta", {"height_80%_CI": [0, pt.MISSING]})
    for children in tree.child_dict.values():
        for branch in children.values():
            branch.label = "support"
            branch.length_sd = branch.length_median = 0.0
            branch.extra_meta = {"length": [17], "length_80%_CI": [0, pt.MISSING]}
    return tree


def _length_edit_data(tree):
    return (tree.root, tree.nodes.copy(), tree.root_suffix,
            {node: (nd, nd.height, nd.state, copy.deepcopy(nd.extra_meta))
             for node, nd in tree.nodedict.items()},
            [(branch, branch.label, branch.length_sd, branch.length_median,
              copy.deepcopy(branch.extra_meta))
             for children in tree.child_dict.values() for branch in children.values()])


def _assert_length_edit_data(tree, before):
    root, nodes, suffix, node_data, branch_data = before
    assert tree.root == root and tree.nodes == nodes and tree.rooted is True
    assert tree.root_suffix is suffix
    assert suffix.length == 0.5 and suffix.label == "root"
    assert suffix.extra_meta == {"stem": [pt.MISSING]}
    for node, (record, height, state, extra) in node_data.items():
        assert tree.nodedict[node] is record
        assert (record.height, record.state, record.extra_meta) == (height, state, extra)
    for record, label, sd, median, extra in branch_data:
        assert (record.label, record.length_sd, record.length_median, record.extra_meta) == (label, sd, median, extra)


class Test_setlength:

    def test_length_assignment_retains_annotations_and_refreshes_length_queries(self, length_edit_tree):
        tree = length_edit_tree
        before = _length_edit_data(tree)
        parents, distances = tree.parent_dict, tree.rootdist
        tree.setlength(tree.parent("A"), "A", 4)
        _assert_length_edit_data(tree, before)
        assert tree.rootdist is not distances and tree.rootdist["A"] == 6
        assert tree.parent_dict is parents

    def test_caller_scaling_loop_changes_only_structural_lengths(self, length_edit_tree):
        tree = length_edit_tree
        before = _length_edit_data(tree)
        total = tree.length()
        _ = tree.rootdist
        for parent, children in list(tree.child_dict.items()):
            for child in list(children):
                tree.setlength(parent, child, tree.child_dict[parent][child].length * 2)
        _assert_length_edit_data(tree, before)
        assert tree.length() == total * 2 and tree.rootdist["A"] == 6


class Test_clear_length_caches:

    def test_raw_length_assignment_requires_explicit_invalidation(self, length_edit_tree):
        tree = length_edit_tree
        before = _length_edit_data(tree)
        cached = tree.rootdist
        tree.child_dict[tree.parent("A")]["A"].length = 4
        assert tree.rootdist is cached and cached["A"] == 3
        tree.clear_length_caches()
        _assert_length_edit_data(tree, before)
        assert tree.rootdist is not cached and tree.rootdist["A"] == 6


class Test_NewickStringParser_parse:
    @pytest.mark.parametrize('text', ['', '(A,B)', '(A:1,B:2', '(A,B));', '(A,B),(C,D);', '(A,B);(C,D);', '((A,B);', '(A,B),C;', '(A,B))R;', '(A,B)):1;', '(A,B)),;'])
    def test_rejects_incomplete_or_multiple_topologies(self, text):
        with pytest.raises(pt.TreeError):
            pt.NewickStringParser().parse(pt.Tree(), text)

    def test_reuse_releases_tree_after_success_and_failure(self):
        parser = pt.NewickStringParser()
        first = parser.parse(pt.Tree(), '(A:1,B:2);')
        with pytest.raises(pt.TreeError):
            parser.parse(pt.Tree(), '(X,');
        second = parser.parse(pt.Tree(), '(C:3,D:4);')
        assert first.leaves == {'A', 'B'} and second.leaves == {'C', 'D'}
        assert parser.treeobj is None and parser.node_stack == []

    def test_duplicate_translated_taxa(self):
        parser = pt.NewickStringParser({'1': 'A', '2': 'A'})
        with pytest.raises(pt.TreeError, match='Duplicated leafnames'):
            parser.parse(pt.Tree(), '(1,2);')

    def test_root_suffix_is_structured_and_not_an_edge(self):
        parser = pt.NewickStringParser(label_attr_name='posterior', label_type=float)
        tree = parser.parse(pt.Tree(), '(A:1,B:2)0.95:0;')
        assert tree.root_suffix.posterior == 0.95
        assert tree.root_suffix.length == 0.0
        assert tree.rooted is None and tree.parent_dict[tree.root] is None
        assert tree.nodedict == {}

    def test_label_only_root_suffix_has_no_invented_length(self):
        tree = pt.Tree.from_string('(A,B)Root;')
        assert tree.root_suffix.label == 'Root' and tree.root_suffix.length is None


class Test_MetacommentParser_parse:
    """Typed annotation payload grammar, independently of positional ownership."""

    @pytest.mark.parametrize("text", [
        '[&x=1,y={2,3}]',
        'x=1,y={2,3}',
        ' \n[& x = 1 , y = { 2, 3 } ]\t',
    ])
    def test_comment_or_payload(self, text):
        assert pt.MetacommentParser().parse(text) == {'x': 1, 'y': (2, 3)}

    def test_scalar_types_and_recursive_collections(self):
        values = pt.MetacommentParser().parse(
            '[&integer=-12,positive=+12,zero=0,decimal=1.0,exponent=-2e-3,'
            'yes=true,no=false,missing=?,empty={},'
            'nested={1,1,?,{},{{2},3.0},""}]'
        )
        assert values['integer'] == -12 and type(values['integer']) is int
        assert values['positive'] == 12 and type(values['positive']) is int
        assert values['zero'] == 0 and type(values['zero']) is int
        assert values['decimal'] == 1.0 and type(values['decimal']) is float
        assert values['exponent'] == -0.002 and type(values['exponent']) is float
        assert values['yes'] is True and values['no'] is False
        assert values['missing'] is pt.MISSING
        assert values['empty'] == ()
        assert values['nested'] == (1, 1, pt.MISSING, (), ((2,), 3.0), '')
        assert type(values['nested'][4][1]) is float

    def test_complete_bare_tokens_and_nonfinite_strings(self):
        tokens = ('north_coast', '#0000FF', '1.2.3', '1e-', '123ABC',
                  'True', 'TRUE', 'NA', 'None', 'null', 'nan', 'NaN',
                  'inf', '+inf', '-Infinity', '??', '?suffix', '0x10',
                  'a:b/c;(d)!', '&R')
        text = '[&' + ','.join('k' + str(i) + '=' + token
                               for i, token in enumerate(tokens)) + ']'
        values = pt.MetacommentParser().parse(text)
        assert tuple(values.values()) == tokens
        assert all(type(value) is str for value in values.values())

    def test_exact_nonempty_key_identity(self):
        values = pt.MetacommentParser().parse(
            '[&height=1,"Height"=2,01=3,!color=#0000FF,'
            'p:q/r;(s)=-1," meta key "=4,"&W"=5,'
            "'O''Brien'=6,\"Unicode [nordøst]\"=7]"
        )
        assert tuple(values) == ('height', 'Height', '01', '!color',
                                 'p:q/r;(s)', ' meta key ', '&W',
                                 "O'Brien", 'Unicode [nordøst]')
        assert values['01'] == 3 and values['!color'] == '#0000FF'

    @pytest.mark.parametrize('text,expected', [
        ('[&x=""]', ''),
        ("[&x='']", ''),
        ("[&x='O''Brien']", "O'Brien"),
        ('[&x="said ""hello"""]', 'said "hello"'),
        ('[&x=""""]', '"'),
        ("[&x='''']", "'"),
        ("[&x='said \"hello\"']", 'said "hello"'),
        ('[&x="O\'Brien"]', "O'Brien"),
        (r'[&x="\n\u0041\\"]', r'\n\u0041\\'),
        ('[&x=" north\ncoast ø "]', ' north\ncoast ø '),
        ('[&x=",{}=:;()"]', ',{}=:;()'),
        ('[&x="sample [north [coastal]]"]', 'sample [north [coastal]]'),
    ])
    def test_quoted_strings(self, text, expected):
        value = pt.MetacommentParser().parse(text)['x']
        assert value == expected and type(value) is str

    def test_quoted_scalars_stay_strings_recursively(self):
        value = pt.MetacommentParser().parse('[&x={"1",\'true\',"?","None",""}]')['x']
        assert value == ('1', 'true', '?', 'None', '')
        assert all(type(item) is str for item in value)

    @pytest.mark.parametrize('text', [
        '[&x=1,x=2]', '[&x=1,x=1]', '[&x=?,x=?]',
        '[&height=1,"height"=1]', "[&'O''Brien'=1,\"O'Brien\"=2]",
    ])
    def test_duplicate_decoded_keys(self, text):
        with pytest.raises(pt.TreeError, match='duplicate annotation key'):
            pt.MetacommentParser().parse(text)

    @pytest.mark.parametrize('text', [
        '[&x={1 2}]', '[&x=1 y=2]', '[&x={,1}]', '[&,x=1]',
        '[&x={1,,2}]', '[&x=1,,y=2]', '[&x={1,}]', '[&x=1,y=2,]',
        '[&x={{1,}}]', '[&x=]', '[&x={1,}]', '[&x={1,',
        '[&x=1,,]', '[&x={1],y=2]',
    ])
    def test_strict_value_and_separator_positions(self, text):
        with pytest.raises(pt.TreeError):
            pt.MetacommentParser().parse(text)

    @pytest.mark.parametrize('text', [
        '[&]', '[& ]', '', '  \n ', '[&""=1]', "[&''=1]",
        '[&=1]', '[&x 1.2]', '[&x=north coast]', '[&x=ab[note]cd]',
        '[&ab[note]cd=1]', '[&x="unterminated]', "[&x='unterminated]",
        '[&x=1]extra', '[&x=1][&y=2]',
    ])
    def test_malformed_payloads(self, text):
        with pytest.raises(pt.TreeError):
            pt.MetacommentParser().parse(text)

    @pytest.mark.parametrize('text', [
        '[&x="["]', '[&x="]"]', '[&x="][north"]',
        '[&left="[",right="]"]', '[&x={1,{"["}},y="]"]',
        '[&"["=1,"]"=2]',
    ])
    def test_each_quoted_key_and_value_balances_brackets(self, text):
        with pytest.raises(pt.TreeError, match='unbalanced square brackets'):
            pt.MetacommentParser().parse(text)

    @pytest.mark.parametrize('text', [
        '[&W 0.5]', '[&&NHX:S=human]', '[&height 1.2]',
        '[&W "text=does-not-make-a-pair"]', '[&R]', '[& U ]',
    ])
    def test_directives_and_flags_are_not_key_value_lists(self, text):
        with pytest.raises(pt.TreeError):
            pt.MetacommentParser().parse(text)

    def test_flag_named_and_quoted_ampersand_assignments(self):
        assert pt.MetacommentParser().parse(
            '[&R=1,U=2,r=3,"&W"=4]'
        ) == {'R': 1, 'U': 2, 'r': 3, '&W': 4}

    def test_mapping_is_not_applied_by_the_payload_helper(self):
        assert pt.MetacommentParser().parse(
            '[&height="?",height_95%_CI={1,"?"},length=7]'
        ) == {'height': '?', 'height_95%_CI': (1, '?'), 'length': 7}

    def test_reuse_after_success_and_failure(self):
        parser = pt.MetacommentParser()
        assert parser.parse('[&x=1]') == {'x': 1}
        with pytest.raises(pt.TreeError, match='duplicate'):
            parser.parse('[&x=2,x=3]')
        assert parser.parse('[&x=?,y={}]') == {'x': pt.MISSING, 'y': ()}
        with pytest.raises(pt.TreeError, match='unterminated quoted'):
            parser.parse('[&x="unfinished]')
        assert parser.parse('x=4') == {'x': 4}

    def test_writer_compatible_typed_nested_values(self):
        original = (pt.MISSING, True, False, 0, 1, 1.0, -0.0, 1e-3,
                    '1', '?', 'nordøst [coast]', 'said "yes"', (), ('x',))
        key = 'quoted key [north]'
        text = ('[&' + pt._format_metadata_key(key) + '='
                + pt._format_metadata_value(original, None) + ']')
        decoded = pt.MetacommentParser().parse(text)[key]
        assert decoded == original
        assert decoded[0] is pt.MISSING
        assert type(decoded[3]) is int and type(decoded[5]) is float
        assert math.copysign(1.0, decoded[6]) == -1.0

    def test_error_has_local_offset_key_and_excerpt(self):
        with pytest.raises(pt.TreeError) as caught:
            pt.MetacommentParser().parse('[&height=1,"height"=2]')
        message = str(caught.value)
        assert 'offset' in message and "'height'" in message and 'near' in message


class Test_AnnotatedNewickParser_parse:
    def test_positional_ownership_native_fields_and_suffix(self):
        text = '[&R] ((A[&height=2,location="north coast"]:[&length=7]0.5[&rate=0.7],B:1)Inner[&height_95%_HPD={?,2.5}]:2,C:3)Root[&height=4]:0[&posterior=0.9];'
        tree = pt.Tree.from_string(text, parse_meta=True)
        assert tree.rooted is True and tree.root_suffix.label == 'Root'
        assert tree.root_suffix.length == 0.0 and tree.root_suffix.posterior == 0.9
        assert tree.nodedict['A'].height == 2
        assert tree.nodedict['A'].extra_meta == {'location': 'north coast'}
        branch = tree.child_dict[tree.parent('A')]['A']
        assert branch.length == 0.5 and branch.extra_meta == {'length': 7, 'rate': 0.7}
        inner = tree.parent('A')
        assert tree.nodedict[inner].extra_meta == {'height_95%_HPD': (pt.MISSING, 2.5)}
        assert tree.nodedict[tree.root].height == 4
        assert tree._nodeheightdict is None

    @pytest.mark.parametrize('parse_meta', [False, True])
    def test_declarations_are_structural_in_both_modes(self, parse_meta):
        tree = pt.Tree.from_string('[ordinary [nested];] [& r ] (A,B)Root:1;', parse_meta=parse_meta)
        assert tree.rooted is True and tree.root_suffix.label == 'Root'
        assert tree.root_suffix.length == 1
        for text in ('[&R][&R](A,B);', '[&R][&U](A,B);', '(A[&U],B);', '(A,B):1[&R];'):
            with pytest.raises(pt.TreeError):
                pt.Tree.from_string(text, parse_meta=parse_meta)

    def test_disabled_payloads_are_not_decoded(self, monkeypatch):
        def forbidden(*args):
            raise AssertionError('metadata decoder called')
        monkeypatch.setattr(pt.MetacommentParser, 'parse', forbidden)
        tree = pt.Tree.from_string('[&U][&W 0.5](A[&x=bad grammar],B)[&R=1]:2[&&NHX:x];')
        assert tree.rooted is False and tree.nodedict == {}
        assert tree.root_suffix.length == 2 and tree.root_suffix.extra_meta is None

    @pytest.mark.parametrize('text', ['(A[&x=1][&x=1],B);', '(A[&height=1,"height"=2],B);', '(A:[&x=1]1[&x=2],B);', '(A[&height="2"],B);', '(A[&height_95%_HPD={1}],B);', '(A[&W 0.5],B);', '[&lnP=-12](A,B);'])
    def test_metadata_errors_are_not_silently_dropped(self, text):
        with pytest.raises(pt.TreeError):
            pt.Tree.from_string(text, parse_meta=True)

    def test_disjoint_keys_merge_and_node_branch_names_remain_distinct(self):
        tree = pt.Tree.from_string('(A[&x=1][ordinary][&y=2]:1[&x=3],B);', parse_meta=True)
        assert tree.nodedict['A'].extra_meta == {'x': 1, 'y': 2}
        assert tree.child_dict[tree.root]['A'].extra_meta == {'x': 3}

    def test_node_branch_mapping_is_owner_specific(self):
        tree = pt.Tree.from_string('(A[&length_median=2,height=1]:3[&height=4,length_median=5],B);', parse_meta=True)
        assert tree.nodedict['A'].height == 1
        assert tree.nodedict['A'].extra_meta == {'length_median': 2}
        branch = tree.child_dict[tree.root]['A']
        assert branch.length_median == 5 and branch.extra_meta == {'height': 4}

    def test_reusable_parser_clears_declaration_and_duplicate_state(self):
        parser = pt.AnnotatedNewickParser(parse_meta=True)
        first = parser.parse(pt.Tree(), '[&U](A[&x=1],B);')
        with pytest.raises(pt.TreeError):
            parser.parse(pt.Tree(), '[&R](A[&x=1,x=2],B);')
        last = parser.parse(pt.Tree(), '(A[&x=2],B);')
        assert first.rooted is False and last.rooted is None
        assert parser.treeobj is None and parser._meta_seen is None
        assert parser.meta_parser._text == ''

    def test_balanced_comment_brackets_are_not_quote_shielded(self):
        text = '(A[&left="[",right="]"],B);'
        assert pt.Tree.from_string(text).nodedict == {}
        with pytest.raises(pt.TreeError, match='brackets'):
            pt.Tree.from_string(text, parse_meta=True)
        for text in ('(A[nested [unfinished],B);', '(A],B);'):
            with pytest.raises(pt.TreeError):
                pt.Tree.from_string(text)

    def test_full_precision_roundtrip_preserves_types_owners_and_state(self):
        tree = pt.Tree.from_string('[&U](A[&i=1,f=1.0,b=false,s="true",empty="",seq={1,?,{}}]:0.10000000000000002[&length=9],B:2)Root[&height=2.0000000000000004]:0.5[&rate=1e-10];', parse_meta=True)
        encoded = tree.newick(printlabels=True, print_meta=True,
                             node_attributes=('i','f','b','s','empty','seq','height'),
                             branch_attributes=('length','rate'))
        decoded = pt.Tree.from_string("[&U]" + encoded, parse_meta=True)
        assert decoded.rooted is False and decoded.root_suffix.label == 'Root'
        assert decoded.root_suffix.extra_meta['rate'] == 1e-10
        assert decoded.nodedict['A'].extra_meta == tree.nodedict['A'].extra_meta
        assert type(decoded.nodedict['A'].extra_meta['f']) is float
        assert type(decoded.nodedict['A'].extra_meta['i']) is int
        assert decoded.nodedict[decoded.root].height == 2.0000000000000004
        assert decoded.child_dict[decoded.root]['A'].length == 0.10000000000000002


class Test_AnnotatedNewickParser_tokenize:
    @pytest.mark.parametrize('parse_meta', [False, True])
    def test_ordinary_comments_within_words_do_not_split_tokens(self, parse_meta):
        tree = pt.Tree.from_string('(A[note; [nested]]B:1[note].5,C:2)Ro[note]ot;', parse_meta=parse_meta)
        assert tree.leaves == {'AB', 'C'}
        assert tree.child_dict[tree.root]['AB'].length == 1.5
        assert tree.root_suffix.label == 'Root'

    def test_duplicate_block_error_reports_the_second_occurrence(self):
        text = '(A[&x=1][&x=1],B);'
        with pytest.raises(pt.TreeError) as caught:
            pt.Tree.from_string(text, parse_meta=True)
        assert f'Tree offset {text.rindex("[&x=1]")}' in str(caught.value)


class Test_comment_end:

    def test_nested_comments_use_brackets_without_quote_interpretation(self):
        text = '[ordinary "[" [nested] "]"; end]rest'
        assert pt._comment_end(text, 0) == text.index("rest")
        assert pt._comment_end("[ordinary [nested]", 0) is None
        quoted_close = '[x="]";rest'
        assert pt._comment_end(quoted_close, 0) == quoted_close.index("]") + 1


class Test_read_tree_text_Newicktreefile:

    @pytest.mark.parametrize("buffer_size", [7, 65536])
    def test_raw_preserves_comments_semicolons_root_suffix_and_adjacent_trees(self, monkeypatch, buffer_size):
        monkeypatch.setattr(pt.TreefileBase, "_BUFFER_SIZE", buffer_size)
        first = '[&R](A[&note="x;y"],B)Root[&height=3]:0.0[&rate=2];'
        second = '(C[ordinary "[" [nested] "]"; note],D);'
        with pt.Newicktreefile(filecontent=first + second + " [trailing comment]") as reader:
            assert reader.read_tree_text() == first
            assert reader.get_treestring() == second
            assert reader.readtree(False) is None
            assert reader.read_tree_text() is None
            assert not hasattr(reader, "below_root")

    @pytest.mark.parametrize("text", ["(A,B)", "(A[unfinished;B);", "(A],B);"])
    def test_unfinished_and_stray_comment_input_fails_lexically(self, text):
        with pt.Newicktreefile(filecontent=text) as reader:
            with pytest.raises(pt.TreeError):
                reader.read_tree_text()

    def test_empty_content_and_comment_only_tail_are_clean_eof(self):
        for text in ("", " \n[complete [nested] comment] \n"):
            with pt.Newicktreefile(filecontent=text) as reader:
                assert reader.read_tree_text() is None


class Test_read_tree_text_Nexustreefile:

    def test_concatenated_complete_nexus_documents_remain_supported(self):
        first = pt.Tree.from_string("(A:1,B:2);")
        first.rooted = True
        second = pt.Tree.from_string("(C:3,D:4);")
        second.rooted = False
        content = first.nexus(translateblock=True) + second.nexus(translateblock=True)
        with pt.Nexustreefile(filecontent=content) as reader:
            trees = list(reader)
        assert [tree.leaves for tree in trees] == [{"A", "B"}, {"C", "D"}]
        assert [tree.rooted for tree in trees] == [True, False]
        assert trees[0].nodedist("A", "B") == 3.0
        assert trees[1].nodedist("C", "D") == 7.0

    @pytest.mark.parametrize("middle", ["#NEXUSgarbage\n", "#INVALID\n"])
    def test_repeated_signature_is_validated(self, middle):
        content = '#NEXUS\nbegin trees;tree x=(A,B);end;' + middle + 'begin trees;tree y=(C,D);end;'
        with pt.Nexustreefile(filecontent=content) as reader:
            assert reader.readtree().leaves == {"A", "B"}
            with pytest.raises(pt.TreeError, match="Expected #NEXUS"):
                reader.readtree()

    def test_repeated_signature_inside_an_unclosed_block_is_not_discarded(self):
        content = '#NEXUS\nbegin trees;tree x=(A,B);#NEXUS\nbegin trees;tree y=(C,D);end;'
        with pt.Nexustreefile(filecontent=content) as reader:
            assert reader.readtree().leaves == {"A", "B"}
            with pytest.raises(pt.TreeError):
                reader.readtree()

    @pytest.mark.parametrize("buffer_size", [7, 65536])
    def test_headers_unrelated_blocks_and_original_post_equals_text(self, monkeypatch, buffer_size):
        monkeypatch.setattr(pt.TreefileBase, "_BUFFER_SIZE", buffer_size)
        first = ' [&R] (1[&height=1],2)Root:0.0[&rate=2];'
        second = '\n(1,2);'
        content = ('[leading [comment]]\n#NEXUS\n'
                   'begin taxa; dimensions ntax=2; taxlabels A B; end;'
                   'begin figtree; set value=[ordinary; [nested]]; tree fake=(X,Y); end;'
                   'begin trees; translate 1 A, 2 B;'
                   'tree first [p=0.2] =' + first
                   + 'utree second=' + second + 'end;')
        with pt.Nexustreefile(filecontent=content) as reader:
            assert reader.read_tree_text() == first
            assert reader._last_command_kind == "TREE"
            parser = reader.parser_obj
            assert reader.readtree(False) == second
            assert reader._last_command_kind == "UTREE"
            assert reader.parser_obj is parser
            assert reader.transdict == {"1": "A", "2": "B"}
            assert reader.read_tree_text() is None
            assert not hasattr(reader, "below_root")

    def test_translation_and_command_context_remain_stable_after_later_block(self):
        content = ('#NEXUS\nbegin trees;translate 1 A,2 B;'
                   'tree first = (1,2);end;'
                   'begin trees;translate 1 C,2 D;'
                   'utree second = (1,2);end;')
        with pt.Nexustreefile(filecontent=content) as reader:
            first = reader.read_tree_text()
            first_context = reader.parser_obj, reader._last_command_kind
            first_mapping = reader.transdict
            second = reader.read_tree_text()
            second_context = reader.parser_obj, reader._last_command_kind
            assert reader.transdict is not first_mapping
            assert first_mapping == {"1": "A", "2": "B"}
            assert first_context[0] is not second_context[0]
            # Workers receive only raw text and stable parser/command context.
            old_parser, old_kind = pickle.loads(pickle.dumps(first_context))
            one = pt.Tree._from_string_private(old_parser, first, command_kind=old_kind)
            two = pt.Tree._from_string_private(second_context[0], second,
                                               command_kind=second_context[1])
            assert one.leaves == {"A", "B"}
            assert one.rooted is None
            assert two.leaves == {"C", "D"}
            assert two.rooted is False

    def test_later_tree_block_without_translate_does_not_reuse_old_mapping(self):
        content = ('#NEXUS\nbegin trees;translate 1 A,2 B;tree x=(1,2);end;'
                   'begin trees;tree y=(1,2);end;')
        with pt.Nexustreefile(filecontent=content) as reader:
            assert reader.readtree().leaves == {"A", "B"}
            assert reader.readtree().leaves == {"1", "2"}
            assert reader.transdict is None


class Test_readtree_TreefileBase:

    @pytest.mark.parametrize("parse_meta", [False, True])
    @pytest.mark.parametrize("nexus", [False, True])
    def test_root_flags_native_values_and_suffix_use_same_parser(self, parse_meta, nexus):
        description = '[&R](A[&height=1.0],B)Root[&height=3]:0.0[&rate=2];'
        reader_type = pt.Nexustreefile if nexus else pt.Newicktreefile
        content = '#NEXUS\nbegin trees;tree x=' + description + 'end;' if nexus else description
        with reader_type(filecontent=content, parse_meta=parse_meta) as reader:
            tree = next(reader)
            assert tree.rooted is True
            assert tree.root_suffix.label == "Root"
            assert tree.root_suffix.length == 0.0
            assert not hasattr(tree, "below_root")
            if parse_meta:
                assert tree.nodedict["A"].height == 1.0
                assert tree.nodedict[tree.root].height == 3
                assert tree.root_suffix.extra_meta == {"rate": 2}
            else:
                assert tree.nodedict == {}
                assert tree.root_suffix.extra_meta is None
            assert reader.readtree() is None
            with pytest.raises(TypeError):
                reader.__next__(False)


class Test_skip_trees_TreefileBase:

    @pytest.mark.parametrize("nexus", [False, True])
    def test_skip_ignores_payload_and_declaration_semantics_without_parsing(self, monkeypatch, nexus):
        first = '(A[&height={}],B);'
        second = '[&R][&U](C,D);'
        third = '(E[ordinary; [nested]],F);'
        reader_type = pt.Nexustreefile if nexus else pt.Newicktreefile
        content = ('#NEXUS\nbegin trees;tree a=' + first + 'utree b=' + second
                   + 'tree c=' + third + 'end;') if nexus else first + second + third
        with reader_type(filecontent=content, parse_meta=True) as reader:
            with monkeypatch.context() as patch:
                def forbidden(*args, **kwargs):
                    raise AssertionError("skip constructed or parsed a Tree")
                patch.setattr(pt.Tree, "__init__", forbidden)
                patch.setattr(pt.AnnotatedNewickParser, "parse", forbidden)
                patch.setattr(reader, "read_tree_text", forbidden)
                assert reader.skip_trees(2) == 2
            assert reader.readtree().leaves == {"E", "F"}
            assert reader.skip_trees(5) == 0


class Test_fast_treecount:

    @pytest.mark.parametrize("fileformat", ["newick", "nexus"])
    def test_count_is_lexical_without_tree_or_description_construction(self, monkeypatch, tmp_path, fileformat):
        descriptions = ('(A[&height={}],B);', '[&R][&U](C,D);', '(E[&note="x;y"],F);')
        if fileformat == "nexus":
            text = '#NEXUS\nbegin trees;' + ''.join(
                f'tree t{index}={value}' for index, value in enumerate(descriptions)) + 'end;'
        else:
            text = ''.join(descriptions)
        filename = tmp_path / "trees.txt"
        filename.write_text(text)
        def forbidden(*args, **kwargs):
            raise AssertionError("count constructed topology or a full tree description")
        monkeypatch.setattr(pt.Tree, "__init__", forbidden)
        monkeypatch.setattr(pt.AnnotatedNewickParser, "parse", forbidden)
        monkeypatch.setattr(pt.Newicktreefile, "read_tree_text", forbidden)
        monkeypatch.setattr(pt.Nexustreefile, "read_tree_text", forbidden)
        assert pt.fast_treecount(filename, fileformat=fileformat) == 3


class Test_count_trees_by_parsing:

    def test_historical_namespace_and_string_fallback_call_shapes(self, tmp_path):
        filename = tmp_path / "trees.newick"
        filename.write_text('(A[&x="a;b"],B);(C,D);')
        assert pt.count_trees_by_parsing(filename, SimpleNamespace(informat="newick")) == 2
        assert pt.count_trees_by_parsing(filename, "newick") == 2
        assert pt.count_trees_by_parsing(filename, args=SimpleNamespace(informat="newick")) == 2


class Test_new_Treefile:

    @pytest.mark.parametrize("nexus", [False, True])
    def test_signature_detection_and_removed_flag_arguments(self, tmp_path, nexus):
        filename = tmp_path / "trees.txt"
        content = '[leading [nested]]\n#nexus\nbegin trees;tree x=(A,B);end;' if nexus else '(A,B);'
        filename.write_text(content)
        with pt.Treefile(filename, parse_meta=True) as reader:
            assert isinstance(reader, pt.Nexustreefile if nexus else pt.Newicktreefile)
            assert reader.readtree().leaves == {"A", "B"}
        with pytest.raises(TypeError):
            pt.Treefile(filename, False)
        with pytest.raises(TypeError):
            pt.Treefile(filename, strip_comments=False)
        with pytest.raises(TypeError):
            pt.Newicktreefile(None, '(A,B);', False)
        with pytest.raises(TypeError):
            pt.Newicktreefile(filecontent='(A,B);', strip_comments=False)
        with pytest.raises(TypeError):
            pt.Nexustreefile(filecontent='#NEXUS\nbegin trees;end;', strip_comments=False)


@pytest.fixture
def annotation_producer_input(data_dir):
    """Extract the first one-line description from these retained producer fixtures.

    This deliberately bounded test helper is not an alternative Nexus reader.
    """
    import re

    def first(relative_path):
        text = (data_dir / relative_path).read_text()
        table = re.search(r'\btranslate\b(.*?);', text, re.IGNORECASE | re.DOTALL)
        assert table is not None
        translation = dict(entry.split() for entry in table.group(1).split(','))
        tree = re.search(r'^\s*tree\b[^=]*=\s*(.*)$', text, re.IGNORECASE | re.MULTILINE)
        assert tree is not None and tree.group(1).rstrip().endswith(';')
        return tree.group(1).strip(), translation

    return first

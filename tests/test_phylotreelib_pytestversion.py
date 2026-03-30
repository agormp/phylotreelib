import phylotreelib as pt
import copy
import itertools
import math
import numpy as np
import pytest
import random
import textwrap
from string import ascii_lowercase, digits

###################################################################################################
###################################################################################################

# Tests for loose functions

###################################################################################################
###################################################################################################

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

    def test_randtree_leaflist(self):
        """Does Tree.randtree() work correctly when leaflist is provided?"""
        leaflist = ['leaf1', 'leaf2', 'leaf3', 'leaf4', 'leaf5']
        mytree = pt.Tree.randtree(leaflist=leaflist)

        assert isinstance(mytree, pt.Tree)
        assert set(mytree.leaves) == set(leaflist)
        assert mytree.is_resolved()

    def test_randtree_ntips(self):
        """Does Tree.randtree() work correctly when ntips is provided?"""
        ntips = 5
        mytree = pt.Tree.randtree(ntips=ntips)

        assert isinstance(mytree, pt.Tree)
        assert len(mytree.leaves) == ntips
        assert mytree.is_resolved()

    def test_randtree_randomlen(self):
        """Does Tree.randtree() correctly assign random lengths when randomlen is True?"""
        ntips = 5
        mytree = pt.Tree.randtree(ntips=ntips, randomlen=True)

        assert isinstance(mytree, pt.Tree)
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
        t = pt.Tree.randtree(leaflist=namelist)
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
            tree = pt.Tree.randtree(ntips=ntips, randomlen=True)
            original_leaves = tree.leaves.copy()
            tree.spr()
            assert tree.leaves == original_leaves

    def test_random_spr_all_possible_params(self):
        """Test all possible combinations of prune and regraft nodes for range of random trees
        of different size. No exceptions should be raised"""
        for ntips in range(3,9):
            origtree = pt.Tree.randtree(ntips=ntips, randomlen=True)
            original_leaves = origtree.leaves.copy()
            for prune_node in origtree.possible_spr_prune_nodes():
                for regraft_node in origtree.possible_spr_regraft_nodes(prune_node):
                    t = origtree.copy_treeobject()
                    t.spr(prune_node, regraft_node)
                    assert t.leaves == original_leaves

    def test_spr_with_prune_only(self):
        """When only a prune_node is given, spr() should choose a valid regraft node."""
        tree = pt.Tree.randtree(ntips=7, randomlen=True)
        possible_prune_nodes = tree.possible_spr_prune_nodes()
        prune_node = next(iter(possible_prune_nodes))
        original_leaves = tree.leaves.copy()
        tree.spr(prune_node=prune_node)
        assert tree.leaves == original_leaves

    def test_spr_with_invalid_regraft_node(self):
        """Supplying a regraft_node that is not allowed should raise an error."""
        tree = pt.Tree.randtree(ntips=8, randomlen=True)
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
        tree = pt.Tree.randtree(ntips=15, randomlen=True)
        base_len = tree.length()
        for _ in range(50):
            tree.spr()
            assert tree.length() == pytest.approx(base_len)

    def test_spr_rejects_tree_with_two_leaves(self):
        tree = pt.Tree.from_string("(A:1,B:1);")
        with pytest.raises(pt.TreeError):
            tree.spr()

    def test_spr_rejects_regraft_at_root(self):
        tree = pt.Tree.randtree(ntips=10, randomlen=True)
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

            t1 = pt.Tree.randtree(ntips=ntips, randomlen=True)

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

            t1 = pt.Tree.randtree(ntips=ntips, randomlen=True)
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

        # Own computation
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.trees" )
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="hip", blen="cladeheight")

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.treeannot_hipstr_mean")
        tgold = tf.readtree()

        # Comparisons
        assert town.topology_clade == tgold.topology_clade # To get separate error message
        assert town.equals(tgold, rooted=True)  # Checks rooted topology again, and blens

    def test_mrhipstr_cladeheight(self, data_dir):
        """mrHIPSTR tree with mean node heights and rooting at best HIPSTR resolution of root clade"""

        # Own computation
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.trees" )
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="mrhip", blen="cladeheight")

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.treeannot_mrhipstr_mean")
        tgold = tf.readtree()

        # Comparisons
        assert town.topology_clade == tgold.topology_clade # To get separate error message
        assert town.equals(tgold, rooted=True)  # Checks rooted topology again, and blens

    def test_mcc_cladeheight(self, data_dir):
        """MCC tree with mean node heights and rooting at best tree's original root"""

        # Own computation
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.trees" )
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, tracktopo=True,
                              track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="mcc", blen="cladeheight")

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.treeannot_mcc_mean")
        tgold = tf.readtree()

        # Comparisons
        assert town.topology_clade == tgold.topology_clade # To get separate error message
        assert town.equals(tgold, rooted=True)  # Checks rooted topology again, and blens

    def test_mcc_caheight(self, data_dir):
        """MCC tree with CA node heights and rooting at best tree's original root"""

        # Own computation
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.trees" )
        tsum = pt.TreeSummary(trackclades=True, trackheight=True, tracktopo=True,
                              track_subcladepairs=True)
        for t in tf:
            tsum.add_tree(t)
        town = tsum.compute_sumtree(treetype="mcc", blen="caheight",
                                    count_burnin_filename_list=[(1000, 0, ( data_dir / "random_10tips_1000.trees" ))])

        # Gold standard computed by treeannotator
        tf = pt.Nexustreefile( data_dir / "random_10tips_1000.treeannot_mcc_ca")
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

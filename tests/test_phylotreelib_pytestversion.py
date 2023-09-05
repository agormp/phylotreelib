import phylotreelib as pt
import math
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

    def test_nested_comments(self):
        text = 'Hello [outer [inner] comment] World!'
        assert pt.remove_comments(text) == 'Hello  World!'

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
        branch = pt.Branchstruct(length=mylength, label=mylabel)
        assert branch.length == mylength
        assert branch.label == mylabel

class Test_copy_Branchstruct:
    def test_copy(self):
        mylength = 0.1543
        mylabel = "0.9965"
        b1 = pt.Branchstruct(length=mylength, label=mylabel)
        b1.myatt1 = "howdy"
        b1.myatt2 = 0.34984
        b2 = b1.copy()
        assert b1.length == b2.length
        assert b1.label == b2.label
        assert b1.myatt1 == b2.myatt1
        assert b1.myatt2 == b2.myatt2

###################################################################################################
###################################################################################################

# Tests for Topostruct class

###################################################################################################
###################################################################################################

class Test_create_Topostruct:

    def test_attributes(self):
        ts = pt.Topostruct()
        w = 0.345
        t = "(A,(B,(C,D)));"  # not actually a tree but fine for testing...
        f = 0.97
        x = "not supposed to be an attribute"
        # Note: main test here is to see if class accepts (only) attribues named in __slots__
        # Assertion is done implicitly: if exception is raised here then test fails
        ts.weight = w
        ts.tree = t
        ts.freq = f
        with pytest.raises(AttributeError):
            ts.notanatt = x

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
                    lablist.append(origtree.getlabel(parent, child))
            newtree = pt.Tree.from_branchinfo(parentlist, childlist, lenlist, lablist)
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
        my_tree = pt.Tree.from_branchinfo(parentlist, childlist, lenlist, lablist)

        # Expected representation
        expected_str_lines =   ["|-----------------------------------------|",
                                "|  Node  |  Child  |  Distance  |  Label  |",
                                "|-----------------------------------------|",
                                "|     0  |      1  |         1  |   0.95  |",
                                "|     0  |      A  |         2  |         |",
                                "|     0  |      B  |         2  |         |",
                                "|     1  |      C  |         1  |         |",
                                "|     1  |      D  |         1  |         |",
                                "|-----------------------------------------|",
                                "",
                                "4 Leafs:",
                                "-----",
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
        """Test copy_treeobject with copylengths=True, copylabels=True"""
        for treestring in treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=True, copylabels=True)
            assert t1 == t2
            assert t1.root == t2.root
            assert t1.has_same_root(t2)
            for parent in t1.intnodes:
                for kid in t1.children(parent):
                    t1lab = t1.getlabel(parent, kid)
                    t2lab = t2.getlabel(parent, kid)
                    assert t1lab == t2lab

    def test_blen_nolab(self, treedata):
        """Test copy_treeobject with copylengths=True, copylabels=False"""
        for treestring in treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=True, copylabels=False)
            assert t1 == t2
            assert t1.root == t2.root
            assert t1.has_same_root(t2)
            for parent in t2.intnodes:
                for kid in t2.children(parent):
                    t2lab = t2.getlabel(parent, kid)
                    assert t2lab == ""

    def test_noblen_nolab(self, treedata):
        """Test copy_treeobject with copylengths=False, copylabels=False"""
        for treestring in treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=False, copylabels=False)
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
            assert list(reversed(sorted_nodes_deepfirst)) == sorted_nodes_shallowfirst

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
        expected_result = [0,1,2]

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
                kid = t.children(intnode).pop()
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


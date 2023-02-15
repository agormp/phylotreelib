# -*- coding: utf-8 -*-
# Unit tests for phylotreelib.py
# UNDER CONSTRUCTION (I know - should have been written prior to the code...)
# Simple usage: python3 test_phylotreelib.py
# verbose output: python3 test_phylotreelib.py -v

import phylotreelib as pt
import unittest
import copy
import csv
import tempfile
import os
import itertools
import random
from io import StringIO
import pytest

###################################################################################################
###################################################################################################

class TreeTestBase(unittest.TestCase):
    """Base class for test cases: contains tree data used by many tests"""

    treedata = {"simplestring" :
                        "(((s1:0.12500,s2:0.12500):0.25000,s3:0.12500):0.12500,s4:0.12500,S5:0.12500);",

               "string_with_blanks" :
                        "(((s1:0.12500, s2:0.12500): 0.25000,s3: 0.12500): 0.12500, s4:0.12500 ,S5:0.12500);",

               "string_with_newlines" :
                         """(((s1:0.12500,s2:0.12500):0.25000,
                                            s3:0.12500):0.12500,s4:0.12500,S5:0.12500);""",

               "string_with_label" :
                         """(KL0F07689: 0.101408, KW081_13: 0.071355,
                            (SBC669_26: 0.009364, YAL016W: 0.014955)0.0507 : 0.124263);""",

                "string_with_weird_newlines":
                         """((SIVCZ:0.19506,(((HV1EL:0.07434000000000002,HV1Z2:0.05618000000000001)
                            :0.023249999999999993, HV1Z8:0.09005000000000002):0.03898999999999997,(HV1MN:
                            0.08418999999999993,((HV1C4:0.09445999999999999,
                            HV1KB:0.08659):0.007830000000000004,(HV1W1:0.08096999999999999,(((HV1A2:0.07867000000000002,
                            HV1JR:0.05537000000000003):0.006249999999999978,((HV1BR:0.006839999999999957,
                            HV1H3:0.007709999999999995):0.0038600000000000856,
                            HV1B1:0.008610000000000007):0.06370999999999993):0.005840000000000067,
                            (HV1RH:0.09736999999999996,HV1S1:0.05908000000000002):0.0027200000000000557)
                            :0.005549999999999944):6.300000000000194E-4):0.006160000000000054):0.0441399999
                            9999996):0.05684):0.16576,
                            ((((((HV2BE:0.08623000000000003,HV2D1:0.08114999999999994):0.024009999999999976,
                            HV2SB:0.09655000000000002):0.003249999999999975,(((HV2S2:0.005730000000000013,
                            HV2ST:0.006879999999999997):0.06057000000000001,HV2G1:0.060250000000000026):0.012259999999999993,
                            HV2RO:0.08071):0.011629999999999974):0.004300000000000026,HV2CA:0.09954000000000002):0.0015800000000000258,
                            HV2NZ:0.10410999999999998):0.01651999999999998,((SIVMK:0.008000000000000007,
                            SIVML:0.007839999999999958):0.036290000000000044,SIVM1:0.04786999999999997):0.10249999999999998):0.16576);""",

                "complexstring" :
                         """((gi|13272696|gb|AF346973.1|:0.001343,(gi|13272836|gb|AF346983.1|:0.000510,
                            gi|13272822|gb|AF346982.1|:0.000993):0.000161):0.000026,(gi|13273158|gb|AF347006.1|:0.000499,
                            (gi|13272766|gb|AF346978.1|:0.000332,((gi|13272710|gb|AF346974.1|:0.000222,
                            gi|13272724|gb|AF346975.1|:0.000279):0.000071,
                            gi|13272808|gb|AF346981.1|:0.000242):0.000107):0.000034):0.000081,
                            (((gi|13272668|gb|AF346971.1|:0.000930,((((gi|13273074|gb|AF347000.1|:0.001594,
                            ((gi|13272584|gb|AF346965.1|:0.001348,gi|13273116|gb|AF347003.1|:0.001283):0.000043,
                            ((gi|13272598|gb|AF346966.1|:0.000743,gi|13272682|gb|AF346972.1|:0.000760):0.000043,
                            ((gi|13272850|gb|AF346984.1|:0.000657,(gi|13272920|gb|AF346989.1|:0.000530,
                            (gi|13272934|gb|AF346990.1|:0.000322,
                            gi|13273214|gb|AF347010.1|:0.000680):0.000033):0.000063):0.000279,
                            (((gi|13272654|gb|AF346970.1|:0.000192,gi|13272948|gb|AF346991.1|:0.000121):0.000068,
                            gi|13272780|gb|AF346979.1|:0.000559):0.000244,(gi|13273242|gb|AF347012.1|:-0.000000,
                            gi|13273256|gb|AF347013.1|:0.000063):0.000445):0.000598):0.000043):0.000025):0.000143):0.000016,
                            (((gi|13273284|gb|AF347015.1|:0.000508,gi|13272990|gb|AF346994.1|:0.000244):0.000298,
                            gi|13272612|gb|AF346967.1|:0.000955):0.000034,
                            gi|13272794|gb|AF346980.1|:0.000593):0.000151):0.000027,
                            gi|13273270|gb|AF347014.1|:0.001198):0.000008,((((gi|13272962|gb|AF346992.1|:0.001676,
                            (((gi|13272640|gb|AF346969.1|:0.000000,gi|13273018|gb|AF346996.1|:0.000062):0.000951,
                            (gi|13272626|gb|AF346968.1|:0.000075,gi|13273032|gb|AF346997.1|:0.000239):0.000803):0.000182,
                            gi|13272892|gb|AF346987.1|:0.000852):0.000250):0.000765,
                            gi|13272878|gb|AF346986.1|:0.001647):0.000504,(((gi|5835121|ref|NC_001643.1|:0.020174,
                            gi|5835135|ref|NC_001644.1|:0.019162):0.059964,gi|195972535|emb|AM948965.1|:0.005335):0.004078,
                            ((gi|13273200|gb|AF347009.1|:0.000242,gi|13273186|gb|AF347008.1|:0.000197):0.001923,
                            (gi|13272864|gb|AF346985.1|:0.000389,(gi|13273060|gb|AF346999.1|:0.000051,
                            gi|13273046|gb|AF346998.1|:0.000012):0.000801):0.000959):0.000823):0.000241):0.000610,
                            ((gi|13272752|gb|AF346977.1|:0.000502,gi|13272738|gb|AF346976.1|:0.000375):0.000581,
                            gi|13273004|gb|AF346995.1|:0.001111):0.000438):0.000482):0.000276):0.000025,
                            (gi|13272556|gb|AF346963.1|:0.000776,gi|13272570|gb|AF346964.1|:0.000978):0.000039):0.000093,
                            (((gi|13273088|gb|AF347001.1|:0.000969,((gi|13273172|gb|AF347007.1|:0.000449,
                            gi|13272976|gb|AF346993.1|:0.000490):0.000313,
                            gi|13273228|gb|AF347011.1|:0.000627):0.000033):0.000112,
                            gi|13272906|gb|AF346988.1|:0.000890):0.000045,((gi|13273102|gb|AF347002.1|:0.000553,
                            gi|13273144|gb|AF347005.1|:0.000825):0.000155,
                            gi|13273130|gb|AF347004.1|:0.000597):0.000268):0.000023):0.000019);""",

                "HIVtree" :
                         """((SIVCZ:0.19506,(((HV1EL:0.07434000000000002,HV1Z2:0.05618000000000001):0.023249999999999993,
                            HV1Z8:0.09005000000000002):0.03898999999999997,(HV1MN:0.08418999999999993,((HV1C4:0.09445999999999999,
                            HV1KB:0.08659):0.007830000000000004,(HV1W1:0.08096999999999999,(((HV1A2:0.07867000000000002,
                            HV1JR:0.05537000000000003):0.006249999999999978,((HV1BR:0.006839999999999957,
                            HV1H3:0.007709999999999995):0.0038600000000000856,
                            HV1B1:0.008610000000000007):0.06370999999999993):0.005840000000000067,
                            (HV1RH:0.09736999999999996,HV1S1:0.05908000000000002):0.0027200000000000557)
                            :0.005549999999999944):6.300000000000194E-4):0.006160000000000054):0.04413999999999996):0.05684):0.16576,
                            ((((((HV2BE:0.08623000000000003,HV2D1:0.08114999999999994):0.024009999999999976,
                            HV2SB:0.09655000000000002):0.003249999999999975,(((HV2S2:0.005730000000000013,
                            HV2ST:0.006879999999999997):0.06057000000000001,HV2G1:0.060250000000000026):0.012259999999999993,
                            HV2RO:0.08071):0.011629999999999974):0.004300000000000026,HV2CA:0.09954000000000002):0.0015800000000000258,
                            HV2NZ:0.10410999999999998):0.01651999999999998,((SIVMK:0.008000000000000007,
                            SIVML:0.007839999999999958):0.036290000000000044,SIVM1:0.04786999999999997):0.10249999999999998):0.16576);""",
                "HIVtree_HV1A2root" :
                            """((HV1JR:0.05537000000000003,(((HV1BR:0.006839999999999957,HV1H3:0.007709999999999995):0.0038600000000000856,
                            HV1B1:0.008610000000000007):0.06370999999999993,((HV1RH:0.09736999999999996,HV1S1:0.05908000000000002):0.0027200000000000557,
                            (HV1W1:0.08096999999999999,((HV1C4:0.09445999999999999,HV1KB:0.08659):0.007830000000000004,(HV1MN:0.08418999999999993,
                            (((HV1EL:0.07434000000000002,HV1Z2:0.05618000000000001):0.023249999999999993,HV1Z8:0.09005000000000002):0.03899,
                            (SIVCZ:0.19506000000000004,((((((HV2BE:0.08623,HV2D1:0.08114999999999992):0.024009999999999976,
                            HV2SB:0.09655):0.003249999999999975,(((HV2S2:0.005730000000000013,HV2ST:0.006879999999999997):0.060569999999999985,
                            HV2G1:0.06025):0.012259999999999993,HV2RO:0.08070999999999998):0.011629999999999974):0.004300000000000026,
                            HV2CA:0.09953999999999999):0.0015800000000000258,HV2NZ:0.10410999999999995):0.01651999999999998,
                            ((SIVMK:0.008000000000000007,SIVML:0.007839999999999958):0.036290000000000044,
                            SIVM1:0.04786999999999997):0.10249999999999995):0.33152):0.05684):0.044139999999999985):0.006160000000000054):
                            6.300000000000194E-4):0.005549999999999944):0.005840000000000067):0.006249999999999978):0.03933500000000001,
                            HV1A2:0.03933500000000001);""",
                "HIVminvar":
                            """((((HV1Z8:0.09005,(HV1EL:0.07434,HV1Z2:0.05618):0.02325):0.03899,(HV1MN:0.08419,((HV1W1:0.08097,
                                (((HV1A2:0.07867,HV1JR:0.05537):0.00625,((HV1BR:0.00684,HV1H3:0.00771):0.00386,HV1B1:0.00861):0.06371):0.00584,
                                (HV1S1:0.05908,HV1RH:0.09737):0.00272):0.00555):0.00063,(HV1C4:0.09446,
                                HV1KB:0.08659):0.00783):0.00616):0.04414):0.05684,SIVCZ:0.19506):0.132787,
                                ((((((HV2D1:0.08115,HV2BE:0.08623):0.02401,HV2SB:0.09655):0.00325,(HV2RO:0.08071,
                                (HV2G1:0.06025,(HV2ST:0.00688,HV2S2:0.00573):0.06057):0.01226):0.01163):0.0043,HV2CA:0.09954):0.00158,
                                HV2NZ:0.10411):0.01652,(SIVM1:0.04787,(SIVML:0.00784,SIVMK:0.008):0.03629):0.1025):0.198733);"""
                            }

########################################################################################
########################################################################################

class StringConstruction(TreeTestBase):
    """Tests whether Tree.from_string can parse various valid newick strings and return corresponding Tree object"""

    def test_parse_simplestring(self):
        """Can Tree.from_string parse simple newick string and return Tree object?"""
        treestring = self.treedata["simplestring"]
        self.assertTrue(isinstance(pt.Tree.from_string(treestring), pt.Tree))

    def test_parse_blanks(self):
        """Can Tree.from_string parse newick string with blanks and return Tree object?"""
        treestring = self.treedata["string_with_blanks"]
        self.assertTrue(isinstance(pt.Tree.from_string(treestring), pt.Tree))

    def test_parse_newline(self):
        """Can Tree.from_string parse newick string with newlines and return Tree object?"""
        treestring = self.treedata["string_with_newlines"]
        self.assertTrue(isinstance(pt.Tree.from_string(treestring), pt.Tree))

    def test_parse_label(self):
        """Can Tree.from_string parse newick string with branch labels and return Tree object?"""
        treestring = self.treedata["string_with_label"]
        self.assertTrue(isinstance(pt.Tree.from_string(treestring), pt.Tree))

    def test_parse_complexstring(self):
        """Can Tree.from_string parse complex newick string and return Tree object?"""
        treestring = self.treedata["complexstring"]
        self.assertTrue(isinstance(pt.Tree.from_string(treestring), pt.Tree))

    def test_parse_spuriousnewlines(self):
        """Can Tree.from_string parse newick string with newlines in unexpected places?"""
        treestring = self.treedata["string_with_weird_newlines"]
        self.assertTrue(isinstance(pt.Tree.from_string(treestring), pt.Tree))

    def test_return_correct(self):
        """Is Tree object correct?"""
        # Note: only tested on one treestring so far. Assumed to be representative
        treestring = self.treedata["string_with_label"]
        mytree = pt.Tree.from_string(treestring)
        expected_leaves = {"KL0F07689", "KW081_13", "SBC669_26", "YAL016W"}
        expected_intnodes = {0, 1}
        expected_nodes = expected_leaves | expected_intnodes

        self.assertEqual(mytree.leaves, expected_leaves)
        self.assertEqual(mytree.intnodes, expected_intnodes)
        self.assertEqual(mytree.nodes, expected_nodes)
        self.assertEqual(mytree.children(0), {1, "KL0F07689", "KW081_13"})
        self.assertEqual(mytree.children(1), {"SBC669_26", "YAL016W"})

########################################################################################
########################################################################################

class TreeIteration(TreeTestBase):
    """Tests iteration over file with several trees"""

    def test_iterate_newick(self):
        """Does iteration over tree objects from treefile work? """

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="wt", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        trees = pt.TreeSet()
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            trees.addtree(mytree)
        filehandle.write(trees.newick())
        filehandle.close()

        # Secondly: iterate over treefile
        treefile = pt.Newicktreefile(filename)
        for tree in treefile:
            treestring = tree.newick()

        # Clean up
        os.remove(filename)

    def test_read_correctly_newick(self):
        """Do I also get the correctly read trees from treefile"""

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="wt", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        treelist = []
        trees = pt.TreeSet()
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            trees.addtree(mytree)
        filehandle.write(trees.newick())
        filehandle.close()

        # Secondly: iterate over treefile, check that read trees correspond to written trees
        treefile = pt.Newicktreefile(filename)
        for i, tree in enumerate(treefile):
            self.assertEqual(tree, treelist[i])

        # Clean up
        os.remove(filename)

    def test_iterate_nexus(self):
        """Does iteration over tree objects from nexus treefile work? """

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="a", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            filehandle.write(mytree.nexus())
        filehandle.close()

        # Secondly: iterate over treefile
        treefile = pt.Nexustreefile(filename)
        for tree in treefile:
            treestring = tree.nexus()

        # Clean up
        os.remove(filename)

    def test_read_correctly_nexus(self):
        """Do I also get the correctly read trees from nexus treefile"""

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="a", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        treelist = []
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            filehandle.write(mytree.nexus())
        filehandle.close()

        # Secondly: iterate over treefile, check that read trees correspond to written trees
        treefile = pt.Nexustreefile(filename)
        i = 0
        for tree in treefile:
            self.assertEqual(tree, treelist[i])
            i += 1

        # Clean up
        os.remove(filename)

########################################################################################
########################################################################################

class TreeRead(TreeTestBase):
    """Tests methods for reading (not iterating) trees from treefiles"""

    def test_readtree_newick(self):
        """Test that readtree returns correct trees and None when Newickfile exhausted"""

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="wt", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        treelist = []
        trees = pt.TreeSet()
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            trees.addtree(mytree)
        filehandle.write(trees.newick())
        filehandle.close()

        # Secondly: read all trees from treefile
        # check that read trees correspond to written trees
        ntrees = len(treelist)
        treefile = pt.Newicktreefile(filename)
        for i in range(ntrees):
            tree = treefile.readtree()
            self.assertEqual(tree, treelist[i])

        # Check that None is returned when file is exhausted
        value = treefile.readtree()
        self.assertIsNone(value)

        # Clean up
        os.remove(filename)

    def test_readtree_nexus(self):
        """Test that readtree returns correct trees and None when Nexusfile exhausted"""

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="a", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        treelist = []
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            filehandle.write(mytree.nexus())
        filehandle.close()

        # Secondly: read all trees from treefile
        # check that read trees correspond to written trees
        ntrees = len(treelist)
        treefile = pt.Nexustreefile(filename)
        for i,tree in enumerate(treefile):
            self.assertEqual(tree, treelist[i])

        # Clean up
        os.remove(filename)

    def test_readtrees_newick(self):
        """Test that readtrees (plural) returns correct trees from Newickfile"""

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="wt", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        treelist = []
        trees = pt.TreeSet()
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            trees.addtree(mytree)
        filehandle.write(trees.newick())
        filehandle.close()

        # Secondly: read all trees from treefile in one go using readtrees()
        # check that read trees correspond to written trees
        treefile = pt.Newicktreefile(filename)
        treeset_from_file = treefile.readtrees()
        for origtree, readtree in zip(treelist, treeset_from_file):
            self.assertEqual(origtree, readtree)

        # Clean up
        os.remove(filename)

    def test_readtrees_nexus(self):
        """Test that readtrees (plural) returns correct trees from Nexusfile"""

        # First: construct treefile
        fileobject = tempfile.NamedTemporaryFile(mode="a", encoding="UTF-8", delete=False)
        filename = fileobject.name
        filehandle = fileobject.file
        treelist = []
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            treelist.append(mytree)
            filehandle.write(mytree.nexus())
        filehandle.close()

        # Secondly: read all trees from treefile in one go using readtrees()
        # check that read trees correspond to written trees
        treefile = pt.Nexustreefile(filename)
        newlyreadtrees = []
        for tree in treefile:
            newlyreadtrees.append(tree)
        for origtree, readtree in zip(treelist, newlyreadtrees):
            self.assertEqual(origtree, readtree)

        # Clean up
        os.remove(filename)


########################################################################################
########################################################################################

class LeafConstruction(unittest.TestCase):
    """Tests from_leaf() constructor"""

    def test_fromleaf(self):
        """Does Tree.from_leaf() work correctly?"""
        # Note: this constructor relies mostly on from_string, so testing is minimal...
        leaves = ["A", "B", "C", "D", "E"]
        self.assertTrue(isinstance(pt.Tree.from_leaves(leaves), pt.Tree))
        mytree = pt.Tree.from_leaves(leaves)
        self.assertEqual(mytree.children(mytree.root), {"A", "B", "C", "D", "E"})
        self.assertEqual(len(mytree.intnodes), 1)
        self.assertEqual(len(mytree.leaves), 5)


########################################################################################
########################################################################################

class BiplistConstructor(unittest.TestCase):
    """Tests from_biplist() constructor"""

    def test_frombip(self):
        """Does Tree.from_biplist() work correctly?"""

        biplist = {frozenset([frozenset(["A"]), frozenset(["B", "C", "D"])]): pt.Branchstruct(0.1),
                   frozenset([frozenset(["B"]), frozenset(["A", "C", "D"])]): pt.Branchstruct(0.2),
                   frozenset([frozenset(["C"]), frozenset(["B", "A", "D"])]): pt.Branchstruct(0.3),
                   frozenset([frozenset(["D"]), frozenset(["B", "C", "A"])]): pt.Branchstruct(0.4),
                   frozenset([frozenset(["A", "B"]), frozenset(["C", "D"])]): pt.Branchstruct(0.5)}

        self.assertTrue(isinstance(pt.Tree.from_biplist(biplist), pt.Tree))
        mytree = pt.Tree.from_biplist(biplist)
        # print("\n##############################################\nIn test_frombip\n\n")
        # print(mytree)
        self.assertEqual(mytree.leaves, {"A", "B", "C", "D"})
        self.assertEqual(mytree.intnodes, {0, 1})
        self.assertAlmostEqual(mytree.length(), 1.5)

########################################################################################
########################################################################################

class BranchinfoConstruction(TreeTestBase):
    """Tests from_branchinfo() constructor"""

    def test_frombranchinfo(self):
        """Does Tree.from_branchinfo() work correctly?"""
        # Testing is done by reading trees from treedata, converting them to lists,
        # recreating tree from lists, and comparing to original tree
        for treestring in self.treedata.values():
            origtree = pt.Tree.from_string(treestring)
            parentlist = []
            childlist = []
            lenlist = []
            lablist = []
            for parent in origtree.intnodes:
                for child in origtree.children(parent):
                    parentlist.append(parent)
                    childlist.append(child)
                    lenlist.append(origtree.nodedist(parent,child))
                    lablist.append(origtree.getlabel(parent,child))
            newtree = pt.Tree.from_branchinfo(parentlist, childlist, lenlist, lablist)
            self.assertTrue(origtree == newtree)

########################################################################################
########################################################################################

class EqualityTest(TreeTestBase):
    """Tests special method __eq__() for Tree objects"""

    def test_equality(self):
        """Sanity check: can Tree.__eq__() differentiate between set of known trees?"""
        for treestring in self.treedata.values():
            mytree1 = pt.Tree.from_string(treestring)
            mytree2 = pt.Tree.from_string(treestring)
            self.assertEqual(mytree1, mytree2)

        stringlist = []
        stringlist.append(self.treedata["simplestring"])
        stringlist.append(self.treedata["complexstring"])
        stringlist.append(self.treedata["HIVtree"])
        for i in range(len(stringlist) - 1):
            for j in range(i+1, len(stringlist) - 1):
                mytree1 = pt.Tree.from_string(stringlist[i])
                mytree2 = pt.Tree.from_string(stringlist[j])
                self.assertNotEqual(mytree1, mytree2)

    def test_diffblens(self):
        """Check that different branch lengths results in non-equality"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            for parent in t2.intnodes:
                for kid in t2.children(parent):
                    origblen = t2.nodedist(parent,kid)
                    t2.setlength(parent,kid,origblen*1.5)
            self.assertNotEqual(t1, t2)

########################################################################################
########################################################################################

# To be implemented: test of sorted_intnodes

########################################################################################
########################################################################################

class Has_same_root(TreeTestBase):
    """Tests for has_same_root function"""
    def test_sameroot_sameintnodenumbers(self):
        """Check returns True: two identical trees"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            self.assertTrue(t1.has_same_root(t2))

    def test_sameroot_diffintnodenumbers(self):
        """Check returns True: two trees with same rooted topology, different nodeids"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            delta_id = 35 + max(t1.intnodes)
            for i,id in enumerate(t1.intnodes):
                t2.rename_intnode(id, id+delta_id)
            self.assertTrue(t1.has_same_root(t2))

    def test_diffroot_sameintnodenumbers(self):
        """Check returns False: two trees with same unrooted topology, different roots"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            newroot = random.choice(tuple(t1.intnodes - {t1.root}))
            t2.deroot()
            t2.reroot(newroot, polytomy=True)
            self.assertFalse(t1.has_same_root(t2))


########################################################################################
########################################################################################

class Is_bifurcation(TreeTestBase):
    """Tests for is_bifurcation function"""
    def test_bifurcations(self):
        """Check returns True: internal nodes with two descendants"""
        for i in range(10):
            t = pt.Tree.randtree(ntips=50)
            for intnode in t.intnodes:
                self.assertTrue(t.is_bifurcation(intnode))

    def test_trifurcations(self):
        """Check returns False: internal nodes with three descendants"""
        for i in range(10):
            t = pt.Tree.randtree(ntips=50)
            for intnode in t.intnodes:
                kid = t.children(intnode).pop()
                if kid not in t.leaves:
                    t.remove_branch(intnode, kid)
                    self.assertFalse(t.is_bifurcation(intnode))
                    break


########################################################################################
########################################################################################

class Baseml_rstfile(TreeTestBase):
    """Tests for class Baseml_rstfile"""
    pass

########################################################################################
########################################################################################

class Copy_treeobject(TreeTestBase):
    """Tests for copy_treeobject function"""

    def test_blen_lab(self):
        """Test copy_treeobject with copylengths=True, copylabels=True"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=True, copylabels=True)
            self.assertEqual(t1,t2)
            self.assertEqual(t1.root, t2.root)
            self.assertTrue(t1.has_same_root(t2))
            for parent in t1.intnodes:
                for kid in t1.children(parent):
                    t1lab = t1.getlabel(parent, kid)
                    t2lab = t2.getlabel(parent, kid)
                    self.assertEqual(t1lab,t2lab)

    def test_blen_nolab(self):
        """Test copy_treeobject with copylengths=True, copylabels=False"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=True, copylabels=False)
            self.assertEqual(t1,t2)
            self.assertEqual(t1.root, t2.root)
            self.assertTrue(t1.has_same_root(t2))
            for parent in t2.intnodes:
                for kid in t2.children(parent):
                    t2lab = t2.getlabel(parent, kid)
                    self.assertEqual("",t2lab)

    def test_noblen_nolab(self):
        """Test copy_treeobject with copylengths=True, copylabels=False"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = t1.copy_treeobject(copylengths=False, copylabels=False)
            self.assertEqual(t1.topology(),t2.topology())
            self.assertEqual(t1.root, t2.root)
            self.assertTrue(t1.has_same_root(t2))
            for parent in t2.intnodes:
                for kid in t2.children(parent):
                    t2lab = t2.getlabel(parent, kid)
                    self.assertEqual("",t2lab)


########################################################################################
########################################################################################

class Match_nodes(TreeTestBase):
    """Tests for match_nodes function"""

    def test_sametop_sameroot(self):
        """Test correct output from well formed examples"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            delta_id = 35 + max(t1.intnodes)
            for i,id in enumerate(t1.intnodes):
                t2.rename_intnode(id, id+delta_id)
            intnode1to2,unmatched_root1,unmatched_root2 = t1.match_nodes(t2)
            for id1 in intnode1to2:
                if type(id1) == int:
                    self.assertEqual(id1+delta_id, intnode1to2[id1])
                else:
                    self.assertEqual(id1, intnode1to2[id1])
            self.assertEqual(unmatched_root1, None)
            self.assertEqual(unmatched_root2, None)

    def test_sametop_difroot(self):
        """Test correct output from well formed examples with different root"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            if len(t1.intnodes) > 2:
                root1,root2 = random.sample(tuple(t1.intnodes - {t1.root}), 2)
                t1.reroot(root1, polytomy=True)
                t2.reroot(root2, polytomy=True)
                delta_id = 35 + max(t1.intnodes)
                for i,id in enumerate(t1.intnodes):
                    t2.rename_intnode(id, id+delta_id)
                intnode1to2,unmatched_root1,unmatched_root2 = t1.match_nodes(t2)
                for id1 in intnode1to2:
                    if type(id1) == int:
                        self.assertEqual(id1+delta_id, intnode1to2[id1])
                    else:
                        self.assertEqual(id1, intnode1to2[id1])
                self.assertEqual(unmatched_root1, None) # Both rooted at polytomies
                self.assertEqual(unmatched_root2, None)

    def test_unmatchedroots(self):
        """Test that unmatched roots are reported correctly"""
        t1 = pt.Tree.randtree(ntips=50)
        t2 = t1.copy_treeobject(False,False)
        delta_id = 35 + max(t1.intnodes)
        for i,id in enumerate(t1.intnodes):
            t2.rename_intnode(id, id+delta_id)
        t2.deroot()
        parent = random.choice(tuple(t2.intnodes))
        kid = t2.children(parent).pop()
        blen = t2.nodedist(parent,kid)
        t2.reroot(parent,kid,node1dist=blen/2)
        t1origroot = t1.root
        t2origroot = t2.root
        intnode1to2,unmatched_root1,unmatched_root2 = t1.match_nodes(t2)
        for id1 in intnode1to2:
            if type(id1) == int:
                self.assertEqual(id1+delta_id, intnode1to2[id1])
            else:
                self.assertEqual(id1, intnode1to2[id1])
        self.assertEqual(unmatched_root1, t1origroot)
        self.assertEqual(unmatched_root2, t2origroot)

    def test_difleaves(self):
        """Test error raised when leaves differ"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            randleaf = random.choice(tuple(t1.leaves))
            t1.remove_leaf(randleaf)
            self.assertRaises(pt.TreeError, t1.match_nodes, t2)

    def test_sameleaves_diftop(self):
        """Test error raised when leaves same but topologies differ"""
        for treestring in self.treedata.values():
            t1 = pt.Tree.from_string(treestring)
            t2 = pt.Tree.from_string(treestring)
            while t1.topology() == t2.topology():
                t2.shuffle_leaf_names()
            self.assertRaises(pt.TreeError, t1.match_nodes, t2)


########################################################################################
########################################################################################

class RelationShipMethods(TreeTestBase):
    """Tests methods for determining children, parents, remote children, MRCAs, and basenode in Tree object"""

    def test_adjacent(self):
        """Consistency and function of .children() and .parent() methods"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            for parent in mytree.intnodes:
                children = mytree.children(parent)
                for child in children:
                    result = mytree.parent(child)
                    self.assertEqual(result, parent)

    def test_remote(self):
        """Consistency and function of .remote_children(), .find_mrca(), findbasenode(), and build_remchild_dict methods"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            self.assertEqual(mytree.remote_children(mytree.root), mytree.leaves)
            for basalnode in mytree.intnodes:
                remotekids = mytree.remote_children(basalnode)
                self.assertEqual(basalnode, mytree.find_mrca(remotekids))
                self.assertEqual(basalnode, mytree.findbasenode(remotekids))
            for leaf in mytree.leaves:
                self.assertEqual({leaf}, mytree.remote_children(leaf))
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        self.assertEqual({"s1", "s2", "s3"}, mytree.remote_children(1))  # I hope internal nodes get numbered as I expect...

########################################################################################
########################################################################################

class PathMethods(TreeTestBase):
    """Tests methods for determining nodepath, nodedist, and treelength in Tree object"""

    def test_nodedist(self):
        """Does nodedist() report correct patristic distance?"""
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)

        # Check all pairwise distances in simplestring tree
        # Note: using assertAlmostEqual since I'm comparing floats (computer representation issues)
        # Default behavior is to compare to 7 decimal places
        self.assertAlmostEqual(mytree.nodedist("s1", "s2"), 0.25)
        self.assertAlmostEqual(mytree.nodedist("s1", "s3"), 0.5)
        self.assertAlmostEqual(mytree.nodedist("s1", "s4"), 0.625)
        self.assertAlmostEqual(mytree.nodedist("s1", "S5"), 0.625)
        self.assertAlmostEqual(mytree.nodedist("s2", "s3"), mytree.nodedist("s1", "s3"))
        self.assertAlmostEqual(mytree.nodedist("s2", "s4"), mytree.nodedist("s1", "s4"))
        self.assertAlmostEqual(mytree.nodedist("s2", "S5"), mytree.nodedist("s1", "S5"))
        self.assertAlmostEqual(mytree.nodedist("s3", "s4"), 0.375)
        self.assertAlmostEqual(mytree.nodedist("s3", "S5"), 0.375)
        self.assertAlmostEqual(mytree.nodedist("s4", "S5"), 0.25)

        # Check a few pairwise distances in HIVtree
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)

        # Leaf to leaf distance between SIVMK and HV1H3
        expecteddist = (0.008 + 0.03629 + 0.1025 + 0.16576 + 0.16576 + 0.00771 + 0.00386 +
                        0.06371 + 0.00584 + 0.00555 + 0.00063 + 0.00616 + 0.04414 + 0.05684)
        self.assertAlmostEqual(mytree.nodedist("SIVMK", "HV1H3"), expecteddist)

        # Distance from internal node (ancestor of all HV2s) to leaf (SIVCZ)
        expecteddist = 0.01652 + 0.16576 + 0.16576 + 0.19506
        anc1 = mytree.find_mrca({"HV2BE", "HV2D1", "HV2SB", "HV2S2", "HV2ST", "HV2G1",
                                        "HV2RO", "HV2CA", "HV2NZ"})
        self.assertAlmostEqual(mytree.nodedist("SIVCZ", anc1), expecteddist)

        #  Internal node (ancestor of all HV2s) to internal node (ancestor of HV1EL and HV1Z2)
        expecteddist = 0.01652 + 0.16576 + 0.16576 + 0.05684 + 0.03899 + 0.02325
        anc2 = mytree.find_mrca({"HV1EL", "HV1Z2"})
        self.assertAlmostEqual(mytree.nodedist(anc1, anc2), expecteddist)

    def test_treelength(self):
        """Does treelength() return correct value?"""

        # NOTE: just checking for a few simple cases
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        expectedlength = 0.125 + 0.125 + 0.25 + 0.125 + 0.125 + 0.125 + 0.125
        self.assertAlmostEqual(mytree.length(), expectedlength)

        treestring = self.treedata["string_with_label"]
        mytree = pt.Tree.from_string(treestring)
        expectedlength = 0.101408 + 0.071355 + 0.124263 + 0.009364 + 0.014955
        self.assertAlmostEqual(mytree.length(), expectedlength)

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
                    self.assertEqual(path12, path21)

            for node1 in intnodes:
                for node2 in intnodes:
                    path12 = mytree.nodepath(node1, node2)
                    path21 = mytree.nodepath(node2, node1)
                    path21.reverse()
                    self.assertEqual(path12, path21)

            for node1 in leaves:
                for node2 in intnodes:
                    path12 = mytree.nodepath(node1, node2)
                    path21 = mytree.nodepath(node2, node1)
                    path21.reverse()
                    self.assertEqual(path12, path21)

    def test_nodepath_details(self):
        """Checking parent-child relationships of all node-pairs on nodepath from root to leaves"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            leaves = mytree.leaves
            root = mytree.root
            for leaf in leaves:
                path = mytree.nodepath(root, leaf)
                for i in range(len(path)-1):
                    parent = path[i]
                    child = path[i + 1]
                    self.assertTrue(child in mytree.children(parent))
                    self.assertEqual(parent, mytree.parent(child))

########################################################################################
########################################################################################

# To be implemented: test of cladegrep and shuffle_leaf_names

########################################################################################
########################################################################################

class TreeOutput(TreeTestBase):
    """Tests of stringoutput"""

    # NOTE: I am also testing string output in the iteration class.

    def test_newick_output(self):
        """Sanity check: check that output from newick() can be parsed by Tree.from_string()"""
        for instring in self.treedata.values():
            mytree = pt.Tree.from_string(instring)
            outstring = mytree.newick()
            mytree2 = pt.Tree.from_string(outstring)
            self.assertEqual(mytree, mytree2)

    def test_nexus_output(self):
        """Check that output from nexus() can be parsed by Nexustreefile()"""
        for instring in self.treedata.values():
            mytree = pt.Tree.from_string(instring)
            nexus_string = mytree.nexus()
            nexusfile = pt.Nexustreefile(filecontent=nexus_string)
            mytree2 = next(nexusfile)
            self.assertEqual(mytree, mytree2)

    def test_contree_nexus_output(self):
        """Check that nexus output from result of contree can be parsed by Nexustreefile"""

        # Generate a random tree, create 10 variants using random spr, compute consensus tree
        # This function only checks that nexus string can be parsed, not that contree is correct
        treesummary = pt.TreeSummary()
        rand_tree = pt.Tree.randtree(ntips=50, randomlen=True, name_prefix="test_")
        for i in range(10):
            tree = rand_tree.copy_treeobject()
            most_distant, maxdist = tree.find_most_distant(tree.root, tree.leaves)
            parent = tree.parent(most_distant)
            grandparent = tree.parent(parent)
            subtree_leaves = tree.remote_children(grandparent)
            rest_leaves = tree.leaves - subtree_leaves
            try:
                regraft_node = random.choice(list(rest_leaves))
            except:
                print(rest_leaves)
            tree.spr(grandparent, regraft_node)
            treesummary.add_tree(tree)
        contree = treesummary.contree()
        nexus_string = contree.nexus(print_leaflabels=False)
        nexusfile = pt.Nexustreefile(filecontent=nexus_string)
        mytree = next(nexusfile)

########################################################################################
########################################################################################

class Treesummarytests(TreeTestBase):

    def setUp(self):
        testdir_path = os.path.dirname(__file__)
        self.t1_fname = os.path.join(testdir_path, 'contest.postburnin.1.t')
        self.t2_fname = os.path.join(testdir_path, 'contest.postburnin.2.t')
        con_fname = os.path.join(testdir_path, 'contest.nexus.con.tre')
        mbres_fname = os.path.join(testdir_path, 'bip_mean_var.txt')
        mb_trprobs_fname = os.path.join(testdir_path, 'contest.nexus.trprobs')
        cfile = pt.Nexustreefile(con_fname)
        self.mb_contree = cfile.readtree()
        cfile.close()
        with open(mbres_fname) as mbfile:
            mbresults = mbfile.readlines()
        self.mbresdict = {}
        for line in mbresults:
            names1, names2, meanvar = line.strip().split("|")
            bip1 = frozenset(names1.strip().split())
            bip2 = frozenset(names2.strip().split())
            bipart = frozenset([bip1,bip2])
            vals = meanvar.strip().split()
            mean = float(vals[0])
            var = float(vals[1])
            self.mbresdict[bipart] = [mean,var]
        trprobfile = pt.Nexustreefile(mb_trprobs_fname)
        self.trprob_trees = trprobfile.readtrees()
        trprobfile.close()

    @pytest.mark.slow
    def test_contree(self):
        ts = pt.TreeSummary()
        tf1 = pt.Nexustreefile(self.t1_fname)
        for t in tf1:
            ts.add_tree(t)
        tf2 = pt.Nexustreefile(self.t2_fname)
        for t in tf2:
            ts.add_tree(t)
        own_contree = ts.contree()
        own_bipdict = own_contree.bipdict()
        mb_bipdict = self.mb_contree.bipdict()

        # test topology of own contree is same as topology found by mrbayes
        self.assertEqual(self.mb_contree.topology(), own_contree.topology())

        # test freq and branch length (mean and var) for each branch
        for bip, own_branch in own_bipdict.items():
            mb_mean,mb_var = self.mbresdict[bip]
            self.assertAlmostEqual(mb_mean, own_branch.length)
            # self.assertAlmostEqual(mb_var, own_branch.var) #implement when changed contree to one-pass version

            bip1, bip2 = bip
            if len(bip1)!=1 and len(bip2)!=1:
                mb_branch = mb_bipdict[bip]
                own_freq = float(own_branch.label)
                mb_freq = float(mb_branch.label)
                self.assertAlmostEqual(mb_freq, own_freq)

    @pytest.mark.slow
    def test_treesummary_update(self):
        ts1 = pt.TreeSummary()
        tf1 = pt.Nexustreefile(self.t1_fname)
        for t in tf1:
            ts1.add_tree(t)
        tf1.close()
        ts2 = pt.TreeSummary()
        tf2 = pt.Nexustreefile(self.t2_fname)
        for t in tf2:
            ts2.add_tree(t)
        tf2.close()
        ts1.update(ts2)
        own_contree = ts1.contree()
        own_bipdict = own_contree.bipdict()
        mb_bipdict = self.mb_contree.bipdict()

        # test topology of own contree is same as topology found by mrbayes
        self.assertEqual(self.mb_contree.topology(), own_contree.topology())

        # test freq and branch length (mean and var) for each branch
        for bip, own_branch in own_bipdict.items():
            mb_mean,mb_var = self.mbresdict[bip]
            self.assertAlmostEqual(mb_mean, own_branch.length)
            # self.assertAlmostEqual(mb_var, own_branch.var) #implement when changed contree to one-pass version

            bip1, bip2 = bip
            if len(bip1)!=1 and len(bip2)!=1:
                mb_branch = mb_bipdict[bip]
                own_freq = float(own_branch.freq)
                mb_freq = float(mb_branch.label)
                self.assertAlmostEqual(mb_freq, own_freq, places=3)

    @pytest.mark.slow
    def test_bigtreesummary(self):
        # Note: I am checking that mrbayes topologies are the same I found.
        # should also check topology frequencies
        ts = pt.BigTreeSummary()
        tf1 = pt.Nexustreefile(self.t1_fname)
        for t in tf1:
            ts.add_tree(t)
        tf2 = pt.Nexustreefile(self.t2_fname)
        for t in tf2:
            ts.add_tree(t)
        for tree in self.trprob_trees:
            mb_topology = tree.topology()
            self.assertIn(mb_topology, ts.toposummary)

########################################################################################
########################################################################################

class Topologytests(TreeTestBase):
    """Tests topology related methods"""

    def test_bipdict(self):
        """Does bipdict() return correct result?"""

        # Note: I am exploring just two cases in detail
        treestring = self.treedata["string_with_label"]
        mytree = pt.Tree.from_string(treestring)
        bipdict = mytree.bipdict()

        # Test that key for central branch is present, and that length and label are as expected
        expectedkey = frozenset([frozenset(["YAL016W", "SBC669_26"]),
                               frozenset(["KW081_13", "KL0F07689"])])
        expectedlen = 0.124263
        expectedlab = "0.0507"

        self.assertTrue(expectedkey in list(bipdict.keys()))
        self.assertAlmostEqual(bipdict[expectedkey].length, expectedlen)
        self.assertEqual(bipdict[expectedkey].label, expectedlab)

        # Case #2 - check that all keys are present, and point to expected branch length
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        bipdict = mytree.bipdict()

        expectedkeys = [frozenset([frozenset(["s4", "S5"]), frozenset(["s1", "s2", "s3"])]),
                        frozenset([frozenset(["s4", "S5", "s3"]), frozenset(["s1", "s2"])]),
                        frozenset([frozenset(["s4", "S5", "s3", "s2"]), frozenset(["s1"])]),
                        frozenset([frozenset(["s4", "S5", "s3", "s1"]), frozenset(["s2"])]),
                        frozenset([frozenset(["s4", "S5", "s1", "s2"]), frozenset(["s3"])]),
                        frozenset([frozenset(["s1", "S5", "s3", "s2"]), frozenset(["s4"])]),
                        frozenset([frozenset(["s4", "s1", "s3", "s2"]), frozenset(["S5"])])]

        expectedvals = [0.125, 0.25, 0.125, 0.125, 0.125, 0.125, 0.125]

        for key in bipdict:
            self.assertTrue(key in expectedkeys)

        for i in range(len(expectedkeys)):
            key = expectedkeys[i]
            val = expectedvals[i]
            self.assertAlmostEqual(val, bipdict[key].length)



    def test_topology(self):
        """Does topology() return expected result?"""

        # Explore one case
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        topology = mytree.topology()

        expectedtop = frozenset([frozenset([frozenset(["s4", "S5"]), frozenset(["s1", "s2", "s3"])]),
                        frozenset([frozenset(["s4", "S5", "s3"]), frozenset(["s1", "s2"])]),
                        frozenset([frozenset(["s4", "S5", "s3", "s2"]), frozenset(["s1"])]),
                        frozenset([frozenset(["s4", "S5", "s3", "s1"]), frozenset(["s2"])]),
                        frozenset([frozenset(["s4", "S5", "s1", "s2"]), frozenset(["s3"])]),
                        frozenset([frozenset(["s1", "S5", "s3", "s2"]), frozenset(["s4"])]),
                        frozenset([frozenset(["s4", "s1", "s3", "s2"]), frozenset(["S5"])])])

        self.assertEqual(topology, expectedtop)

    def test_compatiblewith(self):
        """Does is_compatible_with() return expected result?"""

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)

        expectedbips = [frozenset([frozenset(["s4", "S5"]), frozenset(["s1", "s2", "s3"])]),
                        frozenset([frozenset(["s4", "S5", "s3"]), frozenset(["s1", "s2"])]),
                        frozenset([frozenset(["s4", "S5", "s3", "s2"]), frozenset(["s1"])]),
                        frozenset([frozenset(["s4", "S5", "s3", "s1"]), frozenset(["s2"])]),
                        frozenset([frozenset(["s4", "S5", "s1", "s2"]), frozenset(["s3"])]),
                        frozenset([frozenset(["s1", "S5", "s3", "s2"]), frozenset(["s4"])]),
                        frozenset([frozenset(["s4", "s1", "s3", "s2"]), frozenset(["S5"])])]

        badbips  = [frozenset([frozenset(["s3", "S5"]), frozenset(["s1", "s2", "s4"])]),
                        frozenset([frozenset(["s4", "s1", "s3"]), frozenset(["S5", "s2"])])]

        for bip in expectedbips:
            assert mytree.is_compatible_with(bip)

        for bip in badbips:
            self.assertFalse(mytree.is_compatible_with(bip))


    def test_is_resolved(self):
        """Does is_resolved() return expected result?"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            self.assertTrue(mytree.is_resolved(), treestring)

        mytree = pt.Tree.from_string("(A, B, C, D);")
        self.assertFalse(mytree.is_resolved(), )

        mytree = pt.Tree.from_string("(A, (B, C, D));")
        self.assertFalse(mytree.is_resolved(), )

        mytree = pt.Tree.from_string("(A, (B, (C, D), E));")
        self.assertFalse(mytree.is_resolved(), )


    def test_resolve(self):
        # UNDER CONSTRUCTION: use assert method that checks for exception being raised
        namelist = ["seq{}".format(i) for i in range(50)]
        t1 = pt.Tree.from_leaves(namelist)
        t1.resolve()
        t1_string = t1.newick()
        t2 = pt.Tree.from_string(t1_string)


########################################################################################
########################################################################################

class LengthLabel(TreeTestBase):
    """Tests setting and getting of branch lengths and labels"""

    def test_setlength(self):
        """Does setlength() work correctly?"""
        # Check a few cases in HIVtree
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"HV2BE", "HV2D1"})
        mytree.setlength(anc, "HV2BE", 0.666)
        mytree.setlength(anc, "HV2D1", 0.333)
        self.assertAlmostEqual(0.666, mytree.nodedist(anc, "HV2BE"))
        self.assertAlmostEqual(0.333, mytree.nodedist(anc, "HV2D1"))

    def test_labels(self):
        """Does setlabel() and getlabel() work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"HV2BE", "HV2D1"})
        mytree.setlabel(anc, "HV2BE", "Simpel_label")
        mytree.setlabel(anc, "HV2D1", "Label with space")
        self.assertEqual(mytree.getlabel(anc, "HV2BE"), "Simpel_label")
        self.assertEqual(mytree.getlabel(anc, "HV2D1"), "Label with space")

########################################################################################
########################################################################################

class TreeChanging(TreeTestBase):
    """Tests for functions that alter tree structure"""

    def test_insertnode(self):
        """Does insert_node() work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"HV2BE", "HV2D1"})
        pre_nodes = copy.copy(mytree.nodes)

        branchstruct = pt.Branchstruct(length=0.777, label="New_branch")
        newnode = mytree.insert_node(anc, ["HV2BE", "HV2D1"], branchstruct)
        self.assertAlmostEqual(mytree.nodedist(newnode, anc), 0.777)
        self.assertEqual(mytree.getlabel(newnode, anc), "New_branch")
        self.assertEqual(mytree.parent(newnode), anc)
        self.assertEqual(mytree.children(newnode), {"HV2BE", "HV2D1"})
        self.assertEqual(mytree.children(anc), {newnode})
        self.assertEqual(len(pre_nodes) + 1, len(mytree.nodes))
        self.assertTrue(pre_nodes < mytree.nodes)

        anc2 = mytree.parent(anc)
        branchstruct = pt.Branchstruct(length= 0.999, label="Lower_branch")
        newnode2 = mytree.insert_node(anc2, [anc], branchstruct)
        self.assertAlmostEqual(mytree.nodedist(newnode2, anc2), 0.999)
        self.assertEqual(mytree.getlabel(newnode2, anc2), "Lower_branch")
        self.assertEqual(mytree.parent(newnode2), anc2)
        self.assertEqual(mytree.children(newnode2), {anc})
        self.assertEqual(mytree.children(anc2), {newnode2, "HV2SB"})

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"s1", "s2"})
        branchstruct = pt.Branchstruct(length=0.777, label="New_branch")
        newnode = mytree.insert_node(anc, ["s1", "s2"], branchstruct)
        self.assertAlmostEqual(mytree.nodedist(newnode, anc), 0.777)
        self.assertEqual(mytree.getlabel(newnode, anc), "New_branch")
        self.assertEqual(mytree.parent(newnode), anc)
        self.assertEqual(mytree.children(newnode), {"s1", "s2"})
        self.assertEqual(mytree.children(anc), {newnode})

    def test_addbranch(self):
        """Does add_branch() work correctly?"""
        mytree = pt.Tree.from_leaves(["A", "B", "C", "D", "E"])
        bipart1 = frozenset([frozenset(["A", "B"]),
                            frozenset(["C", "D", "E"])])
        bipart2 = frozenset([frozenset(["A", "B", "C"]),
                            frozenset(["D", "E"])])
        pre_nodes = copy.copy(mytree.nodes)           # Copy mutable

        branchstruct = pt.Branchstruct(length=0.222, label="First_branch")
        mytree.add_branch(bipart1, branchstruct)
        anc1 = mytree.parent("A")
        anc2 = mytree.parent("C")
        root = mytree.root
        # print("\n##############################################\nIn test_addbranch\n\n")
        # print(mytree)
        # print("anc1: {}     anc2: {}".format(anc1, anc2))
        self.assertAlmostEqual(mytree.nodedist(anc1, anc2), 0.222)
        self.assertEqual(mytree.getlabel(root, anc1), "First_branch")
        self.assertEqual(mytree.getlabel(root, anc2), "First_branch")
        self.assertEqual(mytree.children(anc1), {"A", "B"})
        self.assertEqual(mytree.children(anc2), {"C", "D", "E"})
        self.assertEqual(len(pre_nodes) + 2, len(mytree.nodes))
        self.assertTrue(pre_nodes < mytree.nodes)

        branchstruct = pt.Branchstruct(length=0.333, label="Second_branch")
        mytree.add_branch(bipart2, branchstruct)
        anc3 = mytree.parent("D")
        anc4 = mytree.parent(anc3)
        # print("\n##############################################\nIn test_addbranch\n\n")
        # print(mytree)
        # print("anc3: {}     anc4: {}".format(anc3, anc4))
        self.assertAlmostEqual(mytree.nodedist(anc3, anc4), 0.333)
        self.assertEqual(mytree.getlabel(anc3, anc4), "Second_branch")
        self.assertEqual(mytree.children(anc3), {"D", "E"})
        self.assertEqual(mytree.children(anc4), {anc3, "C"})


    def test_removebranch(self):
        """Does remove_branch work correctly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        node1 = mytree.find_mrca({"HV1EL", "HV1Z2", "HV1Z8"})
        node2 = mytree.parent(node1)
        node2remotekids = mytree.remote_children(node2)
        removed_dist = mytree.nodedist(node1, node2)
        pre_len = mytree.length()
        pre_nodes = copy.copy(mytree.nodes)

        mytree.remove_branch(node1, node2)
        self.assertEqual(mytree.remote_children(node2), node2remotekids)
        self.assertRaises(pt.TreeError, mytree.children, node1)
        self.assertEqual(pre_nodes - {node1}, mytree.nodes)
        self.assertEqual(len(mytree.children(node2)), 3)
        self.assertAlmostEqual(pre_len - removed_dist, mytree.length())

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
        self.assertEqual(mytree.children(anc2), {"SIVML", "SIVM1"})
        self.assertAlmostEqual(mytree.length(), pre_len - removed_dist)
        self.assertRaises(pt.TreeError, mytree.children, anc1)
        self.assertEqual(len(pre_nodes) - 2, len(mytree.nodes))
        self.assertTrue(pre_nodes > mytree.nodes)


        # Check behavior when leaf directly attached to root in bifurcation
        mytree = pt.Tree.from_string("(A, ((B, C), (D, E)));")
        anc1 = mytree.parent("A")
        anc2 = mytree.find_mrca({"B", "C", "D", "E"})
        pre_root = mytree.root

        mytree.remove_leaf("A")
        self.assertRaises(pt.TreeError, mytree.children, anc1)
        self.assertEqual(pre_root, anc1)
        self.assertEqual(mytree.root, anc2)

    def test_rename_leaf(self):
        """Does rename_leaf() work corectly?"""
        treestring = self.treedata["HIVtree"]
        mytree = pt.Tree.from_string(treestring)
        anc = mytree.find_mrca({"SIVMK", "SIVML"})
        mytree.rename_leaf("SIVMK", "Monkey_Seq")
        self.assertEqual(mytree.children(anc), {"SIVML", "Monkey_Seq"})
        self.assertEqual(mytree.parent("Monkey_Seq"), anc)
        self.assertTrue("SIVMK" not in mytree.remote_children(mytree.root))
        self.assertTrue("SIVMK" not in mytree.nodes)
        self.assertTrue("SIVMK" not in mytree.leaves)
        self.assertTrue("Monkey_Seq" in mytree.nodes)
        self.assertTrue("Monkey_Seq" in mytree.leaves)

        self.assertRaises(pt.TreeError, mytree.rename_leaf, "Not_there", "new_name")

########################################################################################
########################################################################################

class DistPathTester(TreeTestBase):
    """Tests for constructing dist and path dictionaries"""

    def test_distdict(self):
        """Test that dist_dict is constructed correctly and agrees with nodedist function"""
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            mytree.build_dist_dict()
            for n1 in mytree.nodes:
                for n2 in mytree.nodes:
                    self.assertAlmostEqual(mytree.nodedist(n1, n2), mytree.dist_dict[n1][n2])

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        mytree.build_dist_dict()
        self.assertAlmostEqual(mytree.dist_dict["s4"]["s1"], 0.625)
        self.assertAlmostEqual(mytree.dist_dict["s2"]["s1"], 0.25)
        self.assertAlmostEqual(mytree.dist_dict["s3"][0], 0.25)
        self.assertAlmostEqual(mytree.dist_dict[0]["s1"], 0.5)

        # Check that things still work after changing root
        for treestring in self.treedata.values():
            mytree = pt.Tree.from_string(treestring)
            mytree.rootmid()
            mytree.build_dist_dict()
            for n1 in mytree.intnodes:
                for n2 in mytree.leaves:
                    self.assertAlmostEqual(mytree.nodedist(n1, n2), mytree.dist_dict[n1][n2])

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        mytree.rootmid()
        mytree.build_dist_dict()
        self.assertAlmostEqual(mytree.dist_dict["s4"]["s1"], 0.625)
        self.assertAlmostEqual(mytree.dist_dict["s2"]["s1"], 0.25)
        self.assertAlmostEqual(mytree.dist_dict["s3"][0], 0.25)
        self.assertAlmostEqual(mytree.dist_dict[0]["s1"], 0.5)



    def test_pathdict(self):
        """Test that path_dict is constructed correctly and agrees with nodepath function"""
        for treestring in self.treedata.values():
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
                    self.assertEqual(npath, npathnew)

        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)
        mytree.build_path_dict()
        self.assertEqual(mytree.path_dict["s4"]["s1"], 0)
        self.assertEqual(mytree.path_dict[0]["s1"], 1)
        self.assertEqual(mytree.path_dict[1]["s1"], 2)
        self.assertEqual(mytree.path_dict[2]["s1"], "s1")

########################################################################################
########################################################################################

class RootTester(TreeTestBase):
    """Tests for rooting related methods in Tree object"""

    def test_rootmid(self):
        """Test that midpoint rooting creates expected tree on some known examples"""
        treestring = self.treedata["simplestring"]
        mytree = pt.Tree.from_string(treestring)

        mytree.rootmid()
        self.assertAlmostEqual(mytree.nodedist("s1"), 0.3125, msg="leaves: {}".format(mytree.leaves))
        self.assertAlmostEqual(mytree.nodedist("s2"), 0.3125)
        self.assertAlmostEqual(mytree.nodedist("s3"), 0.1875)
        self.assertAlmostEqual(mytree.nodedist("s4"), 0.3125)
        self.assertAlmostEqual(mytree.nodedist("S5"), 0.3125)

        treestring = self.treedata["HIVtree_HV1A2root"]
        mytree = pt.Tree.from_string(treestring)
        mytree.rootmid()
        self.assertAlmostEqual(mytree.nodedist("SIVM1"), mytree.nodedist("HV1RH"))
        self.assertAlmostEqual(mytree.nodedist("SIVM1"), 0.34765)

    def test_rootminvar(self):
        """Test that rootminvar finds the minimum variance root"""

        ts1 = self.treedata["HIVtree"]
        t1 = pt.Tree.from_string(ts1)

        ts2 = self.treedata["HIVminvar"]
        t2 = pt.Tree.from_string(ts2)

        t1.rootminvar()
        self.assertEqual(t1, t2)

########################################################################################
########################################################################################

class PruneTester(TreeTestBase):
    """Tests for pruning related methods in Tree object"""



    def test_prune_maxlen(self):
        """Use brute force to test that prune_maxlen finds the longest tree with given number of leaves"""
        ntips = 14    # Note: combinatorial explosion means execution time rapidly increases
        nkeep = 7
        t1 = pt.Tree.randtree(ntips=ntips, randomlen=True)
        lengths = []
        for discardset in itertools.combinations(t1.leaves, (ntips - nkeep)):
            t2 = copy.deepcopy(t1)
            t2.remove_leaves(discardset)
            lengths.append(t2.length())
        t1.prune_maxlen(nkeep = nkeep)
        self.assertAlmostEqual( t1.length(), max(lengths))

    def test_find_common_leaf(self):
        """Test that find_common_leaf() function returns the typical leaf from clade"""
        # Loop over all treestrings in treedata, and for each check consistency
        for ts in self.treedata.values():
            t = pt.Tree.from_string(ts)
            leaflist = t.leaves
            minsum = t.length() * len(leaflist)  # A value that is larger than any blensum
            # This is really just an independent way of doing what function does. Maybe change test
            for leaf1 in leaflist:
                blensum = 0.0
                for leaf2 in leaflist:
                    blensum += t.nodedist(leaf1, leaf2)
                if blensum < minsum:
                    common_test = leaf1
                    minsum = blensum
            common_function = t.find_common_leaf(leaflist)
            self.assertEqual(common_function, common_test)

    def test_find_central_leaf(self):
        """Test that find_central_leaf() function returns leaf closest to center in clade"""
        t = pt.Tree.from_string("(((A:1,B:1):9,C:9):1,(D:1,E:1):9);")
        central = t.find_central_leaf(t.leaves)
        self.assertEqual( central, "C" )


########################################################################################
########################################################################################

class dist_tree_construction(TreeTestBase):
    """Tests methods for constructing trees from distance matrix"""

    def setUp(self):
        """Set up distance matrices and corresponding tree strings for testing"""
        self.njtreestringlist = ["((A:13.0,B:4.0):2.0,(C:4.0,D:10.0):2.0);",
                                "((A:2,(B:1.5,C:1):3.2):4.1,(D:1.2,(E:2.7,F:3.1):5):5.2);"]
        self.njdmatlist = [
                        {
                            "A": {"B":17, "C":21, "D":27},
                            "B": {"A":17, "C":12, "D":18},
                            "C": {"A":21, "B":12, "D":14},
                            "D": {"A":27, "B":18, "C":14}
                        },
                        {
                            "A": {"B":6.7,  "C":6.2,  "D":12.5, "E":19,   "F":19.4},
                            "B": {"A":6.7,  "C":2.5,  "D":15.2, "E":21.7, "F":22.1},
                            "C": {"A":6.2,  "B":2.5,  "D":14.7, "E":21.2, "F":21.6},
                            "D": {"A":12.5, "B":15.2, "C":14.7, "E":8.9,  "F":9.3},
                            "E": {"A":19,   "B":21.7, "C":21.2, "D":8.9,  "F":5.8},
                            "F": {"A":19.4, "B":22.1, "C":21.6, "D":9.3,  "E":5.8}
                         }
                     ]
        self.upgmatreelist = ["(((A:8.500,B:8.500):2.500,E:11.000):5.500,(C:14.000,D:14.000):2.500);",
                              "((A:1.500,B:1.500):3.750,(C:2.000,D:2.000):3.250);"]
        self.upgmadmatlist = [
                                {
                                    "A":{"B":17, "C":21, "D":31, "E":23},
                                    "B":{"A":17, "C":30, "D":34, "E":21},
                                    "C":{"A":21, "B":30, "D":28, "E":39},
                                    "D":{"A":31, "B":34, "C":28, "E":43},
                                    "E":{"A":23, "B":21, "C":39, "D":43}
                                },
                                {
                                    "A":{"B":3,  "C":11, "D":11},
                                    "B":{"A":3,  "C":10, "D":10},
                                    "C":{"A":11, "B":10, "D":4},
                                    "D":{"A":11, "B":10, "C":4}
                                }
                            ]

        # Apparently necessary to set basepath like this when reading data files for test
        testdir_path = os.path.dirname(__file__)
        dmat_fname = os.path.join(testdir_path, 'large_distmat.tsv')
        njtree_fname = os.path.join(testdir_path, 'large_njtree.txt')
        upgmatree_fname = os.path.join(testdir_path, 'large_upgmatree.txt')

        with open(dmat_fname) as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            dictlist = []
            for rowdict in reader:
                dictlist.append(rowdict)
            self.large_distdict = dict(zip(reader.fieldnames, dictlist))
        with open(njtree_fname) as infile:
            self.largephylip_njtree = infile.read()
        with open(upgmatree_fname) as infile:
            self.largephylip_upgmatree = infile.read()

    def test_nj(self):
        """Verify that nj method produces correct trees for list of distance matrices"""
        for i in range(len(self.njtreestringlist)):
            inputtree = pt.Tree.from_string(self.njtreestringlist[i])
            dmat = pt.Distmatrix.from_distdict(self.njdmatlist[i])
            njtree = dmat.nj()
            self.assertEqual(njtree, inputtree)
        inputtree = pt.Tree.from_string(self.largephylip_njtree)
        dmat = pt.Distmatrix.from_distdict(self.large_distdict)
        njtree = dmat.nj()
        self.assertEqual(njtree, inputtree)

    def test_upgma(self):
        """Verify that upgma method produces correct trees for list of distance matrices"""
        for i in range(len(self.upgmatreelist)):
            inputtree = pt.Tree.from_string(self.upgmatreelist[i])
            dmat = pt.Distmatrix.from_distdict(self.upgmadmatlist[i])
            upgma = dmat.upgma()
            self.assertEqual(upgma, inputtree,  msg='\n  i: {}\n  upgma:\n{}\n  input:\n{}'.format(i, upgma, inputtree))
        inputtree = pt.Tree.from_string(self.largephylip_upgmatree)
        dmat = pt.Distmatrix.from_distdict(self.large_distdict)
        upgmatree = dmat.upgma()
        self.assertEqual(upgmatree, inputtree)

########################################################################################
########################################################################################

class TreeDistTester(TreeTestBase):
    """Tests methods for computing tree distance and similarity"""

    def test_treedist(self):
        ts = self.treedata["HIVtree"]
        t1 = pt.Tree.from_string(ts)
        t2 = pt.Tree.from_string(ts)

        # Move one leaf to a new position, resulting in known treedistance
        subtreenode = t1.parent("HV2BE")
        t2.spr(subtreenode, "SIVM1")
        treedist = t1.treedist(t2, normalise=False)
        self.assertEqual(treedist, 10)     # Note: I really should double check this is correct...

########################################################################################
########################################################################################

if __name__ == "__main__":
    unittest.main()

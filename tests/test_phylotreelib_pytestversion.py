import phylotreelib as pt
import pytest
import random

###################################################################################################
###################################################################################################

# Tests for loose functions

###################################################################################################
###################################################################################################

class Test_remove_comments:

    def test_balanced(self, balanced_input):
        assert pt.remove_comments(balanced_input) == 'Hello  World!\n'

    def test_unbalanced(self, unbalanced_input):
        with pytest.raises(pt.TreeError):
            pt.remove_comments(unbalanced_input)

    def test_nested_comments(self, nested_comments):
        assert pt.remove_comments(nested_comments) == 'Hello  World!'

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
        b1.myatt3 =  {1:1, 2:4, 3:9}
        b2 = b1.copy()
        assert b1.length == b2.length
        assert b1.label == b2.label
        assert b1.myatt1 == b2.myatt1
        assert b1.myatt2 == b2.myatt2
        assert b1.myatt3 == b2.myatt3
        assert b1.myatt3 is not b2.myatt3  # Deepcopy

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
        ts.weight = w    # Note: if exception is raised here then test fails
        ts.tree = t
        ts.freq = f
        with pytest.raises(AttributeError):
            ts.notanatt = x

###################################################################################################
###################################################################################################

# Tests for Tree class

###################################################################################################
###################################################################################################


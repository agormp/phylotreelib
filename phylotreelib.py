"""Classes and methods for analyzing, manipulating, and building phylogenetic trees"""
# Anders Gorm Pedersen
# Section for Bioinformatics, DTU Health Technology, Technical University of Denmark
# agpe@dtu.dk
import copy
import functools
import itertools
import math
import random
import re
import statistics
import sys
from io import StringIO
from operator import itemgetter
import numpy as np

###################################################################################################
###################################################################################################
##
##  Implementation notes:
##        (1) In principle all leafnames are duplicated many times (once per bipartition in
##            bipartsummary and toposummary). However, python appears to automagically
##            map all identical unbroken strings (=no whitespace) to the same id:
##
##            >>> a='tetM_X90939_Streptococcus_pneu'
##            >>> b='tetM_X90939_Streptococcus_pneu'
##            >>> a is b
##            True
##
##            >>> a='does not work if string contains blanks'
##            >>> b='does not work if string contains blanks'
##            >>> a is b
##            False
##
##            However, since this is dependent on implementation of python, I should probably
##            code around it explicitly
##
##            NOTE: I have now added explicit sys.intern() for all leafnames.
##            check that this behaves as expected wrt performance (and remove Globals perhaps)
##
##        (2) Although there is some extra overhead in using classes to emulate structs,
##            using dicts instead does not make a big difference performancewise (I tried).
##
##        (3) Treestrings could be parsed recursively in fewer lines, but this works a lot
##            less efficiently than the iterated version currently in Tree.from_string()
##            (I tested it).
##
###################################################################################################
###################################################################################################

###################################################################################################
###################################################################################################
#
# Various functions used by methods, that do not fit neatly in any class
###################################################################################################

def remove_comments(text):
    """Takes input string and strips away commented text, delimited by '[' and ']'.
        Also deals with nested comments."""

    # Python note: could be simplified

    # Before spending any time:
    # bail if there are no comment delimiters in string
    # raise exception if comment delimiters not balanced
    if "[" not in text:
        return text
    elif text.count("[") != text.count("]"):
        raise TreeError("String contains different number of left and right comment delimiters")

    # Preprocess delims for use in re etc
    leftdelim = re.escape("[")
    rightdelim = re.escape("]")

    # Construct sorted list of tuples of the form [(0, 'start'), (5, 'stop'), (7, 'start'), ...]
    delimlist = [(match.start(), match.end(), "start") for match in re.finditer(leftdelim, text)]
    delimlist.extend([(match.start(), match.end(), "stop") for match in re.finditer(rightdelim, text)])
    delimlist.sort()

    # Traverse text; along the way copy text not inside comment-delimiter pairs.
    # Use stack ("unmatched_starts") to keep track of nesting
    unmatched_starts = 0
    prevpos = 0
    processed_text = []
    for (match_start, match_end, match_type) in delimlist:
        if match_type == "start":
            unmatched_starts += 1
            if unmatched_starts == 1:                               # Beginning of new comment region
                processed_text.append(text[prevpos:match_start])
        elif match_type == "stop":
            unmatched_starts -= 1
            if unmatched_starts == 0:                               # End of comment region
                prevpos = match_end
            elif unmatched_starts == -1:                            # Error: more right delims than left delims
                raise TreeError("Unmatched end-comment delimiter. Context: '{}'".format(text[prevpos-10:prevpos+10]))

    # Add final block of text if relevant (i.e., if text does not stop with rightdelim), return processed text
    if prevpos < len(text):
        processed_text.append(text[prevpos:])
    return "".join(processed_text)

###################################################################################################
###################################################################################################

class Interner():
    """Class used for interning various objects."""

    # Stores dictionaries of objects to be interned: {obj:obj}
    # Interner methods returns *pointer* to object
    # Python note: is there a potential issue with hash collisions
    def __init__(self):
        self.leafset_interndict = {}
        self.clade_interndict = {}
        self.bip_interndict = {}
        self.unhashable_dict = {}

    def intern_leafset(self, leafset):
        return self.leafset_interndict.setdefault(leafset, leafset)

    def intern_clade(self, clade):
        return self.clade_interndict.setdefault(clade, clade)

    def intern_bipart(self, bipart):
        return self.bip_interndict.setdefault(bipart, bipart)

    def store_unhashable(self, name, obj):
        """
        Stores the unhashable object only if it hasn't been stored before.
        Always returns the stored object (whether newly stored or previously stored).
        """
        return self.unhashable_dict.setdefault(name, obj)

###################################################################################################
###################################################################################################

class Branchstruct:
    """Class that emulates a struct. Keeps branch-related info"""

    def __init__(self, length=0.0, label=""):
        self.length = length
        self.label = label

    def __str__(self):
        return f"{str(self.length)}\t´{self.label}´\n"

    def __repr__(self):
        return self.__str__()

    ###############################################################################################

    def copy(self):
        """Returns copy of Branchstruct object, with all attributes included"""
        obj = Branchstruct()
        for attrname, value in vars(self).items():
            setattr(obj, attrname, value)
        return obj

###################################################################################################
###################################################################################################

class Nodestruct:
    """Class that emulates a struct. Keeps node-related info"""

    def __init__(self, depth=0.0):
        self.depth = depth

    def __str__(self):
        return f"{str(self.depth)}\n"

    def __repr__(self):
        return self.__str__()

    ###############################################################################################

    def copy(self):
        """Returns copy of Nodestruct object, with all attributes included"""
        obj = Nodestruct()
        for attrname, value in vars(self).items():
            setattr(obj, attrname, value)
        return obj

###################################################################################################
###################################################################################################

class Topostruct:
    """Class that emulates a struct. Keeps topology-related info"""

    __slots__ = ["weight", "tree", "freq"]

    # Python note: perhaps replace with dataclass, available since python 3.7
    pass

###################################################################################################
###################################################################################################

class Bipartition:
    """Class that represents a bipartition: the split of leaves in two sets that corresponds
    to a branch in a tree. Class has a fast method for hashing and therefore useful when
    comparing Bipartitions"""

    __slots__ = ["leaf_set", "leaf_list", "leaf2index",
                 'indices', '_hash_value']

    def __init__(self, leafset1, all_leaves_set, sorted_leaf_list, leaf2index):
        """Initialise bipartition objects based on only one half of bipartition.
        If the complement of leafset1 is smaller, than that will be stored instead.
        If they have same size: store leafset with smaller hash.

        leafset1: one half of bipartition (to save time in caller)
        all_leaves_set: set of all leaf names
        sorted_leaf_list: sorted list of all leaf names
        leaf2index: dict of {leaf:index in sorted_leaf_list}

        Last 3 params will typically point to attributes in parent Tree object"""

        self.leaf_set = all_leaves_set
        self.leaf_list = sorted_leaf_list
        self.leaf2index = leaf2index
        leafset1 = frozenset(leafset1)

        nleaves = len(self.leaf_set)
        l1 = len(leafset1)
        l2 = nleaves - l1
        if l1 < l2:
            self.indices = sorted([leaf2index[leaf] for leaf in leafset1])
            self._hash_value = hash(leafset1)
        elif l2 < l1:
            leafset2 = frozenset(self.leaf_set) - leafset1
            self.indices = sorted([leaf2index[leaf] for leaf in leafset2])
            self._hash_value = hash(leafset2)
        else:
            leafset2 = frozenset(self.leaf_set) - frozenset(leafset1)
            h1, h2 = hash(leafset1), hash(leafset2)
            if h1 < h2:
                indices = sorted([leaf2index[leaf] for leaf in leafset1])
                self.indices = indices
                self._hash_value = h1
            else:
                indices = sorted([leaf2index[leaf] for leaf in leafset2])
                self.indices = indices
                self._hash_value = h2

    def __hash__(self):
        # Return the precomputed hash value
        return self._hash_value

    def __eq__(self, other):
        # Python note: assumes hash collisions will never occur!
        # Empirically this seems to be the case for realistic data, but is not guaranteed
        return self._hash_value == other._hash_value

    # Python note: this allows unpacking as if the class was a tuple: bip1, bip2 = bipartition
    def __iter__(self):
        # Convert indices to leaf values using set comprehension
        primary_set = frozenset({self.leaf_list[i] for i in self.indices})
        complement_set = frozenset(self.leaf_set) - primary_set
        return iter((primary_set, complement_set))

    def __str__(self):
        bip1,bip2 = self.get_bipartitions()
        if len(bip1) < len(bip2):
            return f"\n{str(bip1)}\n"

    def __repr__(self):
        return self.__str__()

    def get_bipartitions(self):
        # Convert to sets before returning
        primary_set = frozenset({self.leaf_list[i] for i in self.indices})
        complement_set = frozenset(self.leaf_set) - primary_set
        return (primary_set, complement_set)

###################################################################################################
###################################################################################################

class Clade:
    """Class that represents a clade: the set of leaves descended from an internal node"""

    __slots__ = ["all_leaves_set", "leaf_list", "leaf2index",
                 'indices', '_hash_value']

    def __init__(self, leafset, all_leaves_set, sorted_leaf_list, leaf2index):
        """Initialise clade objects.

        leafset: set of leaves descending from an internal node
        all_leaves_set: set of all leaf names in tree
        sorted_leaf_list: sorted list of all leaf names in tree
        leaf2index: dict of {leaf:index in sorted_leaf_list}

        Last 3 params will typically point to attributes in parent Tree object"""

        self.all_leaves_set = all_leaves_set
        self.leaf_list = sorted_leaf_list
        self.leaf2index = leaf2index
        leafset = frozenset(leafset)

        self.indices = sorted([leaf2index[leaf] for leaf in leafset])
        self._hash_value = hash(leafset)

    def __hash__(self):
        # Return the precomputed hash value
        return self._hash_value

    def __eq__(self, other):
        # Python note: assumes hash collisions will never occur!
        # Empirically this seems to be the case for realistic data, but is not guaranteed
        return self._hash_value == other._hash_value

    # Python note: this allows unpacking as if the class was a tuple: c1 = myclade
    def __iter__(self):
        # Convert indices to leaf values using set comprehension
        leaves = frozenset({self.leaf_list[i] for i in self.indices})
        return iter((leaves))

    def __str__(self):
        myclade = self.get_clade()
        return f"\n{str(myclade)}\n"

    def __repr__(self):
        return self.__str__()

    def get_clade(self):
        leaves = frozenset({self.leaf_list[i] for i in self.indices})
        return (leaves)

###################################################################################################
###################################################################################################

class TreeError(Exception):
    pass

###################################################################################################
###################################################################################################

class NewickStringParser:
    """Class creating parser for specific Newick tree string. Used in Tree.from_string"""

    def __init__(self, treeobj, treestring, transdict=None):
        # Construct Tree object that will be filled out by parser
        # Tree is represented as a dictionary of dictionaries. The keys in the top dictionary
        # are the internal nodes which are numbered consecutively. Each key has an
        # associated value that is itself a dictionary listing the children: keys are
        # child nodes, values are Branchstructs containing "length" and "label" fields.
        # Leafs are identified by a string instead of a number

        # NOTE: interprets non-leaf labels as belonging to an internal branch (not to
        # an internal node). The label is attached to the same branch as the branch length
        treeobj.child_dict = {}
        treeobj.leaves = set()
        treeobj.intnodes = set()
        treeobj.root = 0
        self.tree = treeobj

        # Hack to remove whitespace from treestring
        self.treestring = "".join(treestring.split())

        # For keeping track of tree state during parsing (different from parser state)
        self.node_stack = []

        # For checking if there are duplicated leaf names after parsing
        self.ntips = 0

        # If transdict was supplied: use corresponding handler when adding leaves
        if transdict:
            self.transdict = transdict
            self._handle_add_leaf = NewickStringParser._handle_add_leaf_with_transdict
        else:
            self._handle_add_leaf = NewickStringParser._handle_add_leaf

        # Dispatch dictionary: specifies which function to use for any given combination of
        # state and token-type. Functions use token-value as input an return next state.
        self.dispatch = {
            "TREE_START": {
                "(": NewickStringParser._handle_add_root_intnode
            },
            "INTNODE_START": {
                "(": NewickStringParser._handle_add_intnode,
                "NUM_NAME": self._handle_add_leaf
            },
            "LEAF": {
                ":": NewickStringParser._handle_transition_brlen,
                ",": NewickStringParser._handle_transition_child,
                ")": NewickStringParser._handle_intnode_end
            },
            "EXPECTING_BRLEN": {
                "NUM_NAME": NewickStringParser._handle_add_brlen
            },
            "BRLEN": {
                ",": NewickStringParser._handle_transition_child,
                ")": NewickStringParser._handle_intnode_end
            },
            "EXPECTING_CHILD": {
                "(": NewickStringParser._handle_add_intnode,
                "NUM_NAME": self._handle_add_leaf
            },
            "INTNODE_END": {
                ")": NewickStringParser._handle_intnode_end,
                ",": NewickStringParser._handle_transition_child,
                ":": NewickStringParser._handle_transition_brlen,
                "NUM_NAME": NewickStringParser._handle_label,
                ";": NewickStringParser._handle_transition_tree_end
            },
            "LABEL": {
                ")": NewickStringParser._handle_intnode_end,
                ",": NewickStringParser._handle_transition_child,
                ":": NewickStringParser._handle_transition_brlen
            }
        }

    ###############################################################################################

    def parse(self):
        dispatch = self.dispatch
        delimset = set(",;:()")
        state = "TREE_START"
        tree_parts_list = re.split(r'([(),;:])', self.treestring)
        tree_parts_list = list(filter(None, tree_parts_list))  # Remove empty strings from split()
        for token_value in tree_parts_list:
            if token_value in delimset:
                token_type = token_value
            else:
                token_type = "NUM_NAME"
            try:
                handler = dispatch[state][token_type]
                state = handler(self, token_value)
            except KeyError:
                self._handle_parse_error(state, token_value, token_type, self.treestring)

        self.sanitychecks()
        self.tree.nodes = self.tree.leaves | self.tree.intnodes

        return self.tree

    ###############################################################################################

    def sanitychecks(self):
        if self.node_stack:
            msg = f"Nodes remaining on stack after parsing treestring: {self.node_stack}"
            raise TreeError(msg)
        if self.ntips > len(self.tree.leaves):
            raise TreeError(f"Duplicated leafnames in treestring:\n{self.treestring}")

    ###############################################################################################

    def _handle_parse_error(self, state, token_value, token_type, treestring):
        # If unexpected token-type was encountered: first check if parentheses are balanced:
        if treestring.count("(") != treestring.count(")"):
            msg = "Imbalance in tree-string: different number of left- and right-parentheses\n"
            msg += f"Left: ({treestring.count('(')}  Right: {treestring.count(')')})"
            raise TreeError(msg)
        else:
            # If that was not the problem: report current parser-state, token-type, and token-value
            msg = ("Parsing error: unexpected token-type for this state:\n"
                   f"Parser state: {state}\n"
                   f"Token_type:   {token_type}\n"
                   f"Token-value:  {token_value}\n"
                   f"Tree-string:  {treestring}\n")
            raise TreeError(msg)

    ###############################################################################################

    @staticmethod
    def _handle_add_root_intnode(instance, token_value):
        instance.nodeno = 0
        instance.tree.child_dict[instance.nodeno] = {}
        instance.node_stack.append(instance.nodeno)
        instance.tree.intnodes.add(instance.nodeno)
        return "INTNODE_START"

    ###############################################################################################

    @staticmethod
    def _handle_add_intnode(instance, token_value):
        instance.nodeno += 1
        instance.tree.child_dict[instance.nodeno] = {}
        parent = instance.node_stack[-1]
        instance.tree.child_dict[parent][instance.nodeno] = Branchstruct()
        instance.node_stack.append(instance.nodeno)
        instance.tree.intnodes.add(instance.nodeno)
        return "INTNODE_START"

    ###############################################################################################

    @staticmethod
    def _handle_add_leaf(instance, name):
        child = sys.intern(name)
        parent = instance.node_stack[-1]
        instance.tree.child_dict[parent][child] = Branchstruct()
        instance.node_stack.append(child)
        instance.tree.leaves.add(child)
        instance.ntips += 1
        return "LEAF"

    ###############################################################################################

    @staticmethod
    def _handle_add_leaf_with_transdict(instance, name):
        child = sys.intern(instance.transdict[name])
        parent = instance.node_stack[-1]
        instance.tree.child_dict[parent][child] = Branchstruct()
        instance.node_stack.append(child)
        instance.tree.leaves.add(child)
        return "LEAF"

    ###############################################################################################

    @staticmethod
    def _handle_transition_child(instance, token_value):
        instance.node_stack.pop()
        return "EXPECTING_CHILD"

    ###############################################################################################

    @staticmethod
    def _handle_transition_brlen(instance, token_value):
        return "EXPECTING_BRLEN"

    ###############################################################################################

    @staticmethod
    def _handle_add_brlen(instance, brlen_string):
        try:
            brlen = float(brlen_string)
            child = instance.node_stack[-1]
            parent = instance.node_stack[-2]
            instance.tree.child_dict[parent][child].length = brlen
            return "BRLEN"
        except ValueError:
            raise TreeError(f"Expected branch length: {brlen_string}")

    ###############################################################################################

    @staticmethod
    def _handle_intnode_end(instance, token_value):
        instance.node_stack.pop()
        return "INTNODE_END"

    ###############################################################################################

    @staticmethod
    def _handle_label(instance, label):
        child = instance.node_stack[-1]
        parent = instance.node_stack[-2]
        instance.tree.child_dict[parent][child].label = label
        return "LABEL"

    ###############################################################################################

    @staticmethod
    def _handle_transition_tree_end(instance, token_value):
        instance.node_stack.pop()
        return "TREE"

###################################################################################################
###################################################################################################

class Tree:
    """Class representing basic phylogenetic tree object."""

    # Implementation note: Tree objects can be constructed from several different kinds of things:
    # including Newick tree strings, Bipartition lists, and a list of leaves.
    # The Tree class therefore has several alternate constructors implemented as classmethods
    # The main constructor "__init__" is therefore mostly empty

    def __init__(self):
        self._parent_dict = None
        self._remotechildren_dict = None
        self._frozenset_leaves = None
        self._sorted_leaf_list = None
        self._leaf2index = None
        self.dist_dict = None
        self.path_dict = None
        self._remotechildren_dict = None
        self.interner = None
        self._sorted_intnodes_deep = None
        self._sorted_intnodes_shallow = None
        self._rootdist = None
        self._nodedepthdict = None
        self._topology_bipart = None
        self._topology_clade = None

    ###############################################################################################

    def clear_caches(self):
        self._parent_dict = None
        self._remotechildren_dict = None
        self._frozenset_leaves = None
        self._sorted_leaf_list = None
        self._leaf2index = None
        self.dist_dict = None
        self.path_dict = None
        self._remotechildren_dict = None
        self.interner = None
        self._sorted_intnodes_deep = None
        self._sorted_intnodes_shallow = None
        self.nodedist.cache_clear()
        self._rootdist = None
        self._nodedepthdict = None
        self._topology_bipart = None
        self._topology_clade = None

    ###############################################################################################

    @classmethod
    def from_string(cls, orig_treestring, transdict=None, interner=None):
        """Constructor: Tree object from tree-string in Newick format"""

        # All action is in NewickStringParser
        obj = cls()
        parser = NewickStringParser(obj, orig_treestring, transdict)
        obj = parser.parse()
        obj.interner = interner
        if obj.interner:
            obj.leaves = obj.interner.store_unhashable("leaves", obj.leaves)
            obj.intnodes = obj.interner.store_unhashable("intnodes", obj.intnodes)
            obj.nodes = obj.interner.store_unhashable("nodes", obj.nodes)
        del parser
        return obj

    ###############################################################################################

    @classmethod
    def from_biplist(cls, biplist, interner=None):
        """Constructor: Tree object from bipartition list"""

        # Input is a bipartitionlist (actually a dictionary of bipartition:Branchstruct pairs):
        # Names of leaves on one side of a branch are represented as an immutable set of leaves
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree is represented as a dictionary of bipartition:Branchstruct pairs
        obj = cls()

        # Extract set of leaves
        part1, part2 = next(iter(biplist))  # First key=set of two leaf name sets
        obj.leaves = part1 | part2          # Concatenate them to get all leafnames
        obj.intnodes = {0}                  # Will be built as we go along
        obj.root = 0                        # Root is node zero at start
                                            # This may change after re-rooting etc

        # Tree is represented as a dictionary of dictionaries. The keys in the top dictionary
        # are the internal nodes, which are numbered consecutively. Each key has an associated
        # value that is itself a dictionary listing the children: keys are child nodes, values are
        # Branchstructs containing "length" (float) and "label" (str) fields.
        # Leafs are identified by a string instead of a number

        # Construct tree dictionary (main data structure, essentially a child list)
        obj.child_dict = {}
        obj.child_dict[0]={}
        maxnode = 0

        # Start by building star-tree, this is resolved branch-by-branch later on
        # Python note: use startree constructor?
        for leaf in obj.leaves:
            obj.child_dict[0][leaf]= None

        # Iterate over all bipartitions, for each: add extra branch and/or update Branchstruct
        for (bip1, bip2), branchstruct in biplist.items():

            # If bipartition represents external branch: update relevant Branchstruct
            if len(bip1) == 1 or len(bip2) == 1:
                if len(bip1) == 1:
                    (leaf, ) = bip1      # A one-member tuple used for value unpacking (pop)
                else:
                    (leaf, ) = bip2

                # Find childdict containing leaf, and update branchstruct
                for childdict in obj.child_dict.values():
                    if leaf in childdict:
                        childdict[leaf] = branchstruct
                        break

            # If bipartition represents internal branch: add branch to tree, transfer Branchstruct
            else:
                if len(bip1) > len(bip2):
                    bip1,bip2 = bip2,bip1

                # Determine which group of leaves to move
                # Note: one part of bipartition will necessarily have root as its MRCA
                # since the members are present on both sides of root. It is the other part of
                # the bipartition (where all members are on same side of root) that should be moved
                # For a star-tree the resolution will be random (both have root as their MRCA)
                mrca1 = obj.find_mrca(bip1)
                if mrca1 != 0:
                    insertpoint = mrca1
                    active_bip = bip1
                else:
                    mrca2 = obj.find_mrca(bip2)
                    insertpoint = mrca2
                    active_bip = bip2

                # Determine which of insertpoints children to move, namely all children
                # that are either in the active bipartition OR children whose descendants are
                movelist = []
                for child in obj.children(insertpoint):
                    if child in active_bip:
                        movelist.append(child)
                    elif obj.remotechildren_dict[child] <= active_bip:
                        movelist.append(child)

                # Construct new internal node (and therefore branch), transfer Branchstruct
                maxnode += 1
                obj.child_dict[maxnode] = {}
                obj.child_dict[insertpoint][maxnode] = branchstruct
                obj.intnodes.add(maxnode)       # Add new node to list of internal nodes

                # Move relevant children to new node, transfer Branchstructs
                for child in movelist:
                    obj.child_dict[maxnode][child] = obj.child_dict[insertpoint][child]
                    del obj.child_dict[insertpoint][child]

                # reset obj caches which are now obsolete
                obj.clear_caches()

        obj.nodes = set(obj.leaves | obj.intnodes)
        obj.interner = interner
        if obj.interner:
            obj.leaves = obj.interner.store_unhashable("leaves", obj.leaves)
            obj.intnodes = obj.interner.store_unhashable("intnodes", obj.intnodes)
            obj.nodes = obj.interner.store_unhashable("nodes", obj.nodes)
        return obj

    ###############################################################################################

    @classmethod
    def from_topology(cls, topology, interner=None):
        """Constructor: Tree object from topology"""

        # Input is a topology, i.e., a set of bipartitions
        # Names of leaves on one side of a branch are represented as an immutable set of leaves
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # Huge overlap with from_biplist: difference is this doesnt have Branchstructs to begin with
        # Think about mergin code, with flag for branchstruct presence
        obj = cls()

        # Extract set of leaves
        part1, part2 = next(iter(topology)) # First item=set of two leaf name sets
        obj.leaves = part1 | part2          # Concatenate them to get all leafnames
        obj.intnodes = {0}                  # Will be built as we go along
        obj.root = 0                        # Root is node zero at start
                                            # This may change after re-rooting etc

        # Tree is represented as a dictionary of dictionaries. The keys in the top dictionary
        # are the internal nodes, which are numbered consecutively. Each key has an associated
        # value that is itself a dictionary listing the children: keys are child nodes, values are
        # Branchstructs containing "length" (float) and "label" (str) fields.
        # Leafs are identified by a string instead of a number

        # Construct tree dictionary (main data structure, essentially a child list)
        obj.child_dict = {}
        obj.child_dict[0]={}
        maxnode = 0

        # Start by building star-tree, this is resolved branch-by-branch later on
        for leaf in obj.leaves:
            obj.child_dict[0][leaf]= Branchstruct()

        # Iterate over all bipartitions, for each: add extra branch and/or update Branchstruct
        for bip1, bip2 in topology:

            # If bipartition represents internal branch: add branch to tree
            if len(bip1) > 1 and len(bip2) > 1:
                if len(bip1) > len(bip2):
                    bip1,bip2 = bip2,bip1

                # Determine which group of leaves to move
                # Note: one part of bipartition will necessarily have root as its MRCA
                # since the members are present on both sides of root. It is the other part of
                # the bipartition (where all members are on same side of root) that should be moved
                # For a star-tree the resolution will be random (both have root as their MRCA)
                mrca1 = obj.find_mrca(bip1)
                if mrca1 != 0:
                    insertpoint = mrca1
                    active_bip = bip1
                else:
                    mrca2 = obj.find_mrca(bip2)
                    insertpoint = mrca2
                    active_bip = bip2

                # Determine which of insertpoints children to move, namely all children
                # that are either in the active bipartition OR children whose descendants are
                movelist = []
                for child in obj.children(insertpoint):
                    if child in active_bip:
                        movelist.append(child)
                    elif obj.remotechildren_dict[child] <= active_bip:
                        movelist.append(child)

                # Construct new internal node (and therefore branch)
                maxnode += 1
                obj.child_dict[maxnode] = {}
                obj.child_dict[insertpoint][maxnode] = Branchstruct()
                obj.intnodes.add(maxnode)       # Add new node to list of internal nodes

                # Move relevant children to new node, transfer Branchstructs
                for child in movelist:
                    obj.child_dict[maxnode][child] = obj.child_dict[insertpoint][child]
                    del obj.child_dict[insertpoint][child]

            # Reset obj caches which are now obsolete
            obj.clear_caches()

        obj.nodes = set(obj.leaves | obj.intnodes)

        obj.interner = interner
        if obj.interner:
            obj.leaves = obj.interner.store_unhashable("leaves", obj.leaves)
            obj.intnodes = obj.interner.store_unhashable("intnodes", obj.intnodes)
            obj.nodes = obj.interner.store_unhashable("nodes", obj.nodes)

        return obj

    ###############################################################################################

    @classmethod
    def from_cladedict(cls, cladedict, interner=None):

        clade = next(iter(cladedict))
        leaves = clade.all_leaves_set
        obj = cls.from_leaves(leaves)

        for clade,node in cladedict.items():
            clade_leaves = clade.get_clade()
            if len(clade_leaves)>1 and (clade_leaves != leaves):
                mrca = obj.find_mrca(clade_leaves)
                movelist = []
                for child in obj.children(mrca):
                    if child in clade_leaves:
                        movelist.append(child)
                    elif obj.remotechildren_dict[child] <= clade_leaves:
                        movelist.append(child)
                obj.insert_node(mrca, movelist, Branchstruct(label=node.label))

                # Reset obj caches which are now obsolete
                obj.clear_caches()

        obj.nodes = set(obj.leaves | obj.intnodes)
        obj.interner = interner
        if obj.interner:
            obj.leaves = obj.interner.store_unhashable("leaves", obj.leaves)
            obj.intnodes = obj.interner.store_unhashable("intnodes", obj.intnodes)
            obj.nodes = obj.interner.store_unhashable("nodes", obj.nodes)

        return obj

    ###############################################################################################

    @classmethod
    def from_leaves(cls, leaflist, interner=None):
        """Constructor: star-tree object from list of leaves"""

        treelist = ["("]
        for name in leaflist:
            treelist.append(name)
            treelist.append(",")
        del treelist[-1]
        treelist.append(");")
        return cls.from_string("".join(treelist), interner)

    ###############################################################################################

    @classmethod
    def from_branchinfo(cls, parentlist, childlist, lenlist=None, lablist=None, interner=None):
        """Constructor: Tree object from information about all branches in tree

        Information about one branch is conceptually given as:
            parentnodeID, childnodeID, [length], [label]

        The function takes as input 2 to 4 separate lists containing:
            IDs of parents (internal nodes, so integer values)
            ID of children (internal or leaf nodes, so integer or string)
            Length of branches (optional)
            Label of branches (optional)

        The four lists are assumed to have same length and be in same order (so index n in
        each list corresponds to same branch).

        Note: most IDs appear multiple times in lists
        Note 2: can be used as workaround so user can specify IDs for internal nodes"""

        nbranches = len(parentlist)
        if lenlist is None:
            lenlist = [0.0]*nbranches
        if lablist is None:
            lablist = [""]*nbranches
        for lst in [childlist, lenlist, lablist]:
            if len(lst) != nbranches:
                msg = "All lists provided to from_branchinfo() must have same length:\n"
                msg += str(lst)
                raise TreeError(msg)

        obj = cls()                    # Ensures class will be correct also for subclasses of Tree
        obj.child_dict = {}
        obj.leaves = set()
        obj.intnodes = set()

        for i in range(nbranches):
            parent = parentlist[i]     # Perhaps check types are OK?
            child = childlist[i]
            blen = lenlist[i]
            lab = lablist[i]
            if parent in obj.child_dict:
                obj.child_dict[parent][child] = Branchstruct(blen, lab)
            else:
                obj.child_dict[parent] = { child:Branchstruct(blen, lab) }
            obj.intnodes.add(parent)
            if isinstance(child, str):
                obj.leaves.add(child)

        # Root node is the parent node that is not also in childlist
        diffset = set(parentlist) - set(childlist)
        obj.root = diffset.pop()

        obj.nodes = obj.leaves | obj.intnodes

        obj.interner = interner
        if obj.interner:
            obj.leaves = obj.interner.store_unhashable("leaves", obj.leaves)
            obj.intnodes = obj.interner.store_unhashable("intnodes", obj.intnodes)
            obj.nodes = obj.interner.store_unhashable("nodes", obj.nodes)

        return obj

    ###############################################################################################
    @classmethod
    def randtree(cls, leaflist=None, ntips=None, randomlen=False, name_prefix="s", interner=None):
        """Constructor: tree with random topology from list of leaf names OR number of tips"""

        # Implementation note: random trees are constructed by randomly resolving star-tree
        # Should perhaps use actual bifurcating process to generate random trees instead?
        # At least when adding branch lengths (otherwise distribution of brlens on tree will
        # be quite different from real trees, thus biasing statistical inference

        if leaflist is None and ntips is None:
            msg = "Must specify either list of leafnames or number of tips to create random tree"
            raise TreeError(msg)
        if leaflist is not None and ntips is not None:
            msg = "Only specify either list of leafnames or number of tips to create random tree"\
                  " (not both)"
            raise TreeError(msg)

        # If leaflist given:
        #   construct startree using names, then resolve to random bifurcating topology
        if leaflist is not None and ntips is None:
            tree = cls.from_leaves(leaflist, interner)

        # If ntips given:
        #   construct list of zeropadded, numbered, names,
        #   then construct startree using names, finally resolve to random topology
        else:
            ndigits = len(str(ntips))          # Number of digits required to write max taxon number
            namelist = []
            for i in range(ntips):
                name = "{prefix}{num:0{width}d}".format(prefix=name_prefix, num=i, width=ndigits)
                namelist.append( name )
            tree = cls.from_leaves(namelist, interner)   # Star tree with given number of leaves

        tree.resolve()                         # Randomly resolve to bifurcating tree

        if randomlen:
            for parent in tree.intnodes:
                for child in tree.child_dict[parent]:
                    tree.child_dict[parent][child].length = random.lognormvariate(math.log(0.2), 0.3)

        return tree

    ###############################################################################################

    def __iter__(self):
        """Returns iterator object for Tree object. Yields subtrees with .basalbranch attribute"""

        # Basal branch struct may contain useful information: e.g., label and length below subtree
        # Single node trees consisting of leaves are also considered subtrees (REMOVE????)
        class SubtreeIterator():
            def __init__(self, fulltree):
                self.basenodes = fulltree.sorted_intnodes()
                self.basenodes.extend(fulltree.leaflist())
                self.basenodes.remove(fulltree.root)
                self.i = 0
                self.fulltree = fulltree

            def __iter__(self):
                return self

            def __next__(self):
                if self.i >= len(self.basenodes):
                    raise StopIteration
                self.i += 1
                basenode = self.basenodes[self.i - 1]
                (subtree, basalbranch) = self.fulltree.subtree(basenode, return_basalbranch=True)
                subtree.basalbranch = basalbranch
                return subtree

        return SubtreeIterator(self)

    ###############################################################################################

    def __str__(self):
        """Prints table of parent-child relationships including branch lengths and labels"""

        # Starts by building a table of strings containing all information
        # Then formats table into a string

        # Headers
        table = []
        table.append(["Node", "Child", "Distance", "Label"])

        # Build table of parent-child relationships
        for node in self.sorted_intnodes():
            for kid in self.children(node):
                nodstr = str(node)
                kidstr = str(kid)
                dist = "{num:.6g}".format(num=self.child_dict[node][kid].length)
                label = self.child_dict[node][kid].label
                table.append([nodstr, kidstr, dist, label])

        # Find widest string in each column
        maxwidth = [0]*4
        for row in table:
            for i, word in enumerate(row):
                if len(word) > maxwidth[i]:
                    maxwidth[i] = len(word)

        totwidth = maxwidth[0]+maxwidth[1]+maxwidth[2]+maxwidth[3] + 19     # Flanking space

        # Build string from table:
        # Header line
        tabstring = "|" + "-" * (totwidth)  + "|" + "\n"
        for j in range(4):
            tabstring += "|  "+ table[0][j].center(maxwidth[j]) + "  "
        tabstring += "|\n"
        tabstring += "|" + "-" * (totwidth)  + "|" + "\n"

        # Rest of table
        for i in range(1, len(table)):
            for j in range(4):
                tabstring += "|  "+ table[i][j].rjust(maxwidth[j]) + "  "
            tabstring += "|\n"
        tabstring += "|" + "-" * (totwidth)  + "|" + "\n"

        # Add list of leaves to tablestring
        tabstring += "\n%d Leafs:\n" % len(self.leaves)
        tabstring += "-" * maxwidth[1] + "\n"
        for leaf in sorted(self.leaves):
            tabstring += "%s\n" % leaf

        return tabstring

    ###############################################################################################

    def __eq__(self, other, blenprecision=0.005):
        """Implements equality testing for Tree objects"""

        # Two trees are identical if they have the same leaves, the same topology
        # and the same branchlengths. Branch labels are ignored. Rooting is ignored
        # NB: floating point comparison of relative difference.
        # Precision chosen based on empirical comparison between own and PHYLIP tree (...)
        if self.leaves != other.leaves:
            return False
        if self.topology() != other.topology():
            return False
        bipself = self.bipdict()
        bipother = other.bipdict()
        for bipart in bipself:
            len1 = bipself[bipart].length
            len2 = bipother[bipart].length
            if len1 != 0 and len2 != 0:
                if (abs(len1 - len2) / len1) > blenprecision:   # Floating point comparison of relative diff
                    return False
            if (len1 == 0 and len2 > 0) or (len1 > 0 and len2 == 0):
                return False

        # If we made it this far without returning, then Tree objects must be identical
        return True

    ###############################################################################################

    def __hash__(self):
        """Implements hashing for Tree objects, so they can be used as keys in dicts"""

        # Using hash = id of self.
        # NOTE: this does NOT live up to reasonable hash-criteria... Change at some point.
        # NOTE2: Also unsure about effect on performance
        return id(self)

    ###############################################################################################

    @staticmethod
    def is_property(cls, attr):
        return isinstance(getattr(cls, attr, None), property)

    def clear_attributes(self, remlist=None, keeplist=None):
        if remlist and keeplist:
            raise ValueError("Only one of 'remlist' and 'keeplist' can be provided.")

        if not remlist and not keeplist:
            raise ValueError("One of 'remlist' or 'keeplist' must be provided.")

        all_attributes = set(dir(self))

        if keeplist:
            attrs_to_remove = all_attributes - set(keeplist)
        elif remlist:
            attrs_to_remove = set(remlist)

        for attr in attrs_to_remove:
            # Avoid properties, methods, and special attributes (those that start with '__')
            if not attr.startswith("__") and not self.is_property(Tree, attr) and not callable(getattr(self, attr)) :
                setattr(self, attr, None)

    ###############################################################################################

    @property
    def parent_dict(self):
        """Lazy evaluation of _parent_dict when needed"""
        if self._parent_dict == None:
            self.build_parent_dict()
        return self._parent_dict

    ###############################################################################################

    def build_parent_dict(self):
        """Constructs _parent_dict enabling faster lookups, when needed"""

        self._parent_dict = {}
        for parent in self.intnodes:
            for child in self.child_dict[parent]:
                self._parent_dict[child] = parent
        self._parent_dict[self.root] = None     # Add special value "None" as parent of root

    ###############################################################################################

    @property
    def remotechildren_dict(self):
        """Lazy evaluation of _remotechildren_dict when needed"""
        if self._remotechildren_dict == None:
            self.build_remotechildren_dict()
        return self._remotechildren_dict

    ###############################################################################################

    def build_remotechildren_dict(self):
        """Constructs dict of all {parent:{remotechildren}} pairs in efficient manner.
        This dict can then be used directly or by remote_children() to speed up enquiries."""

        remdict = self._remotechildren_dict = {}
        for parent in self.sorted_intnodes(deepfirst=False):
            remdict[parent] = set()
            kidstack = list(self.child_dict[parent])
            while kidstack:
                curnode = kidstack.pop()
                if curnode in self.leaves:
                    remdict[parent].add(curnode)
                elif curnode in remdict:
                    remdict[parent].update(remdict[curnode])
                else:
                    kidstack.extend(self.child_dict[curnode])
        for node in self.leaves:
            remdict[node] = {node}

    ###############################################################################################

    @property
    def frozenset_leaves(self):
        if self._frozenset_leaves == None:
            if self.interner:
                self._frozenset_leaves = self.interner.intern_leafset(frozenset(self.leaves))
            else:
                self._frozenset_leaves = frozenset(self.leaves)
        return self._frozenset_leaves

    ###############################################################################################

    @property
    def sorted_leaf_list(self):
        if self._sorted_leaf_list == None:
            sl = sorted(self.leaves)
            if self.interner:
                self._sorted_leaf_list = self.interner.store_unhashable("sorted_leaf_list", sl)
            else:
                self._sorted_leaf_list = sl
        return self._sorted_leaf_list

    ###############################################################################################

    @property
    def leaf2index(self):
        if self._leaf2index == None:
            self._leaf2index = {}
            for i,leaf in enumerate(self.sorted_leaf_list):
                self._leaf2index[leaf] = i
            if self.interner:
                self._leaf2index = self.interner.store_unhashable("leaf2index", self._leaf2index)
        return self._leaf2index

    ###############################################################################################

    @property
    def rootdist(self):
        """Property (dictionary) giving the distance from each node in tree to the root"""

        if self._rootdist is None:
            self._rootdist = {self.root: 0.0}
            children = list(self.children(self.root))
            while children:
                child = children.pop()
                parent = self.parent_dict[child]
                self._rootdist[child] = (self._rootdist[parent] +
                                         self.child_dict[parent][child].length)
                if child in self.intnodes:
                    children.extend(self.children(child))
        return self._rootdist

    ###############################################################################################

    @property
    def nodedepthdict(self):
        if self._nodedepthdict is None:
            self._nodedepthdict = {}
            maxdist = 0
            for leaf in self.leaves:
                if self.rootdist[leaf] > maxdist:
                    maxdist = self.rootdist[leaf]
            rootdepth = self._nodedepthdict[self.root] = maxdist
            for node in self.nodes - {self.root}:
                self._nodedepthdict[node] = rootdepth - self.rootdist[node]
        return self._nodedepthdict

    ###############################################################################################

    @property
    def topology_bipart(self):
        """Returns set of Bipartitions representation of topology"""

        # Names of leaves on one side of a branch are represented as an immutable set.
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree topology is represented as a set of bipartitions
        # This is essentially a naked version of a bipdict

        if self._topology_bipart == None:
            bipdict = self.bipdict()
            self._topology_bipart = frozenset(bipdict.keys())

        return self._topology_bipart

    ###############################################################################################

    @property
    def topology_clade(self):
        """Returns set of Clades representation of topology"""

        # Names of leaves on one side of a branch are represented as an immutable set.
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree topology is represented as a set of bipartitions
        # This is essentially a naked version of a bipdict

        if self._topology_clade == None:
            cladedict = self.cladedict()
            self._topology_clade = frozenset(cladedict.keys())

        return self._topology_clade

    ###############################################################################################

    def copy_treeobject(self, copylengths=True, copylabels=True, interner=None):
        """Returns copy of Tree object. Copies structure and branch lengths.
        Caches and any user-added attributes are not copied.
        Similar to effect of copy.deepcopy but customized and much faster"""

        obj = Tree()
        obj.root = self.root
        obj.leaves = self.leaves.copy()
        obj.intnodes = self.intnodes.copy()
        obj.nodes = self.nodes.copy()
        obj.interner = interner
        obj.child_dict = {}
        origtree = self.child_dict
        newtree = obj.child_dict
        for parent in origtree:
            newtree[parent] = {}
            for child in origtree[parent]:
                if copylengths:
                    blen = origtree[parent][child].length
                else:
                    blen = 0.0
                if copylabels:
                    lab = origtree[parent][child].label
                else:
                    lab = ""
                newtree[parent][child] = Branchstruct(blen,lab)
        return obj

    ###############################################################################################

    def build_dist_dict(self):
        """Construct dictionary keeping track of all pairwise distances between nodes"""

        # Data structures and algorithm inspired by the Floyd-Warshall algorithm, but modified and
        # faster than O(n^3) since it is on a tree (unique paths)

        # Python note: maybe I could check for existence before recomputing
        # (but then important to clear after changes!)

        # Python note 2: dict.fromkeys does something clever about presizing dict so there is less
        # of a performance hit when it is later added to, hence the slightly odd initialisation
        # (25% faster than dict comprehension)
        dist = self.dist_dict = dict.fromkeys(self.nodes)
        for key in dist:
            dist[key] = {}
        tree = self.child_dict
        combinations = itertools.combinations

        # Traverse tree starting from root, breadth-first (sorted_intnodes)
        # This is required for below algorithm to work
        intnodes = self.sorted_intnodes()
        for parent in intnodes:
            children = tree[parent].keys()
            if dist[parent]:
                prev_contacts = dist[parent].keys()
                for child in children:
                    childlen = tree[parent][child].length
                    for prev_contact in prev_contacts:
                        totlen = dist[prev_contact][parent] + childlen
                        dist[prev_contact][child] = totlen
                        dist[child][prev_contact] = totlen
            for (child1, child2) in combinations(children, 2):
                totlen = tree[parent][child1].length + tree[parent][child2].length
                dist[child1][child2] = totlen
                dist[child2][child1] = totlen
            for child in children:
                totlen = tree[parent][child].length
                dist[parent][child] = totlen
                dist[child][parent] = totlen

        # Fill in diagonal (zero entries), just in case
        for node in self.nodes:
            dist[node][node] = 0

    ###############################################################################################

    def build_path_dict(self):
        """Construct dictionary keeping track of all pairwise paths between nodes"""

        # Data structures and algorithm inspired by the Floyd-Warshall algorithm,
        # but modified and faster than O(n^3) since it is on a tree (unique paths)
        # Python note: dict.fromkeys does something clever about presizing dict so there is less
        # of a performance hit when it is later added to, hence the slightly odd initialisation
        path = self.path_dict = dict.fromkeys(self.nodes)
        for key in path:
            path[key] = {}
        tree = self.child_dict
        combinations = itertools.combinations

        # Traverse tree starting from root, breadth-first (sorted_intnodes)
        # This is required for algorithm to work
        intnodes = self.sorted_intnodes()
        for parent in intnodes:
            children = tree[parent].keys()
            if path[parent]:
                prev_contacts = path[parent].keys()
                for child in children:
                    for prev_contact in prev_contacts:
                        path[prev_contact][child] = path[prev_contact][parent]
                        path[child][prev_contact] = parent
            for (child1, child2) in combinations(children, 2):
                path[child1][child2] = path[child2][child1] = parent
            for child in children:
                path[parent][child] = child
                path[child][parent] = parent

    ###############################################################################################

    def sorted_intnodes(self, deepfirst=True):
        """Returns sorted intnode list for breadth-first traversal of tree"""

        # "intnodes" is a set, meaning iteration occurs in no defined order.
        # This function returns a list sorted such that deep nodes generally go before
        # shallow nodes (deepfirst=False reverses this)
        if deepfirst and self._sorted_intnodes_deep:
            return self._sorted_intnodes_deep
        if not deepfirst and self._sorted_intnodes_shallow:
            return self._sorted_intnodes_shallow

        # Add nodes one tree-level at a time.
        # First root, then children of root, then children of those, etc
        sorted_nodes = []
        curlevel = {self.root}
        while curlevel:
            sorted_nodes.extend(curlevel)
            nextlevel = []

            # For each node in current level: add those children that are also internal nodes
            for node in curlevel:
                nextlevel.extend(self.children(node) & self.intnodes)

            curlevel = nextlevel

        self._sorted_intnodes_deep = sorted_nodes.copy()
        if not deepfirst:
            sorted_nodes.reverse()
            self._sorted_intnodes_shallow = sorted_nodes.copy()

        return sorted_nodes

    ###############################################################################################

    def is_bifurcation(self, node):
        """Checks if internal node is at bifurcation (has two children)"""
        try:
            nkids = len(self.children(node))
            return (nkids == 2)
        except:
            raise TreeError("Node is leaf. Can't check for bifurcation when no children")

    ###############################################################################################

    def n_bipartitions(self):
        """Returns the number of bipartitions (= number of internal branches) in tree
        Note: if root is at bifurcation, then those 2 branches = 1 bipartition"""

        nbip = 0
        for n1 in self.intnodes:
            for n2 in self.children(n1):
                if n2 in self.intnodes:
                    nbip +=1
        if self.is_bifurcation(self.root):
            nbip -= 1
        return nbip

    ###############################################################################################

    def n_branches(self):
        """Returns the number of branches in tree"""
        nbr = 0
        for node in self.intnodes:
            nbr += len(self.child_dict[node])
        return nbr

    ###############################################################################################

    def leaflist(self):
        """Returns list of leaf names sorted alphabetically"""

        leafnamelist = list(self.leaves)
        leafnamelist.sort()
        return leafnamelist

    ###############################################################################################

    def transdict(self):
        """Returns dictionary of {name:number_as_string} for use in translateblocks"""

        leafnamelist = self.leaflist()
        transdict = {}
        for i,leafname in enumerate(leafnamelist):
            transdict[leafname] = f"{i+1}"
        return transdict

    ###############################################################################################

    def translateblock(self, transdict):
        translist = ["    translate\n"]
        for number,name in transdict.items():
            translist.append(f"        {name:<4s}  {number}")
            translist.append(",\n")
        translist[-1] = "\n    ;\n"
        translateblock = "".join(translist)
        return translateblock

    ###############################################################################################

    def children(self, parent):
        """Returns set containing parent's immediate descendants"""

        # Python note: does not seem to benefit from lru_caching, and leads to multiple problems

        try:
            return set(self.child_dict[parent].keys())
        except KeyError as err:
            msg = "Node %s is not an internal node" % parent
            raise TreeError(msg) from err

    ###############################################################################################

    # @functools.lru_cache(maxsize=None)
    def remote_children(self, parent):
        """Returns set containing all leaves that are descendants of parent"""

        # If "parent" is a leaf, then return a set consisting of only itself
        if parent in self.leaves:
            return {parent}

        # Traverse the tree iteratively to find remote children:
        #       if kid is leaf: add it to list of remote children.
        #       if kid is intnode: push its children on stack
        kidstack = set( self.child_dict[parent] )
        remotechildren = set()

        while kidstack:
            curnode = kidstack.pop()
            if curnode in self.leaves:
                remotechildren.add(curnode)
            else:
                kidstack.update( self.child_dict[curnode] )

        return remotechildren

    ###############################################################################################

    def remote_nodes(self, parent):
        """Returns set containing all nodes (intnodes and leaves) that are descendants of parent.
        This set includes parent itself"""

        # If "parent" is a leaf, then return a set consisting of only itself
        if parent in self.leaves:
            return {parent}

        # Traverse the tree iteratively to find remote nodes:
        kidstack = set( self.child_dict[parent] )
        remotenodes = {parent}

        while kidstack:
            curnode = kidstack.pop()
            remotenodes.add(curnode)
            if curnode in self.intnodes:
                kidstack.update( self.child_dict[curnode] )

        return remotenodes

    ###############################################################################################

    def parent(self, node):
        """Returns parent of node"""

        # Python note: should this be distinct from .parent_dict?
        # Perhaps to be used for one-offs when I dont want to trigger construction of full dict
        # Then in code .parent_dict can be used explicitly for when repeated usage is expected
        # Would then perhaps have to be iteration through intnodes until match found.
        try:
            return self.parent_dict[node]
        except KeyError as err:
            raise TreeError("Node {} does not exist (as a key in parent_dict)".format(node)) from err

    ###############################################################################################

    def get_branchstruct(self, node1, node2):
        """Returns Branchstruct object from branch between node1 and node2"""

        if node1 == self.parent(node2):
            parent = node1
            child = node2
        elif node2 == self.parent(node1):
            parent = node1
            child = node2
        else:
            msg = f"Nodes {node1} and {node2} are not adjacent in tree. There is no branch between them"
            raise TreeError(msg)

        return self.child_dict[parent][child]

    ###############################################################################################

    def match_nodes(self, other):
        """Compares two identical trees with potentially different internal node IDs.
        Returns tuple containing following:
            Dictionary giving mapping from nodeid in self to nodeid in other (also leaves)
            unmatched_root1: "None" or id of unmatched root in self if root at bifurcation
            unmatched_root2: "None" or id of unmatched root in other if root at bifurcation

        Note: The last two are only different from None if the trees dont have the same
        exact rooting
        """

        unmatched_root1 = None
        unmatched_root2 = None

        if self.leaves != other.leaves:
            raise TreeError("Trees have different leaves. Can't match intnodes")
        elif self.topology() != other.topology():  # Also checked in sameroot. Could use try
            raise TreeError("Trees have different topologies. Can't match intnodes")
        elif not self.has_same_root(other):

            # Point "self" and "other" to copies of objects so originals are unchanged
            self = self.copy_treeobject(copylengths=False, copylabels=False)
            other = other.copy_treeobject(copylengths=False, copylabels=False)

            # Keep track of original roots if bifurcations
            if self.is_bifurcation(self.root):
                unmatched_root1 = self.root
                self.deroot()
            if other.is_bifurcation(other.root):
                unmatched_root2 = other.root
                other.deroot()

            # Pick arbitrary internal node to root two trees on.
            arbitrarykid = random.choice(tuple(self.leaves))
            newroot1 = self.parent(arbitrarykid)
            newroot2 = other.parent(arbitrarykid)
            self.reroot(newroot1, polytomy=True)
            other.reroot(newroot2, polytomy=True)

        # Now possible to match internal nodes based on their offspring
        node1to2 = dict()
        childpairset = {(leaf,leaf) for leaf in self.leaves}
        while childpairset:
            kid1,kid2 = childpairset.pop()
            parent1 = self.parent(kid1)
            parent2 = other.parent(kid2)
            node1to2[parent1] = parent2
            if parent1 != self.root:
                childpairset.add((parent1,parent2))

        # Add leaves
        # These are identical in two trees, but often useful to have in dict for downstream use
        for leafnode in self.leaves:
            node1to2[leafnode] = leafnode      # Leaf names are same. Add to dict

        return (node1to2, unmatched_root1, unmatched_root2)

    ###############################################################################################

    def rename_intnodes_to_match(self, other):
        """Takes as input a tree (other) with the same topology as self, but with potentially
        different internal nodeIDs. Renames the internal nodeIDs in self so they are the
        same as those in other. Returns copy of self with new nodeIDs"""

        # Create copy of self, and rename all intnodes so there wont be clashes during renaming
        maxselfid = max(self.intnodes)
        minnewid = maxselfid + 1
        newtree = self.copy_treeobject()
        for origid in list(self.intnodes):
            newtree.rename_intnode(origid, origid+minnewid)

        # Now rename
        self2other, unmatched_root1, unmatched_root2 = newtree.match_nodes(other)
        if unmatched_root1 or unmatched_root2:
            raise TreeError("The two trees have different topology - can not match intnodes")
        for intnode in list(newtree.intnodes):
            newtree.rename_intnode(intnode, self2other[intnode])
        return newtree

    ###############################################################################################

    def nearleafs(self, leaf1, maxdist):
        """Returns set of leaves that are less than maxdist from leaf, measured along branches"""

        otherleaves = self.leaves - {leaf1}
        neighbors = set()
        for leaf2 in otherleaves:
            if self.nodedist(leaf1, leaf2) < maxdist:
                neighbors.add(leaf2)

        return neighbors

    ###############################################################################################

    def nearest_n_leaves(self, leaf1, n_neighbors):
        """Returns set of N leaves closest to leaf along tree (patristic distance)"""

        # Python note: numpy.argsort may be faster, but difficult to include ties (n)

        leaflist = self.leaves.copy() - {leaf1}            # Set of all leaves except leaf1
        leaflist = list(leaflist)
        distlist = self.nodedistlist(leaf1, leaflist)
        dist_leaf_list = list(zip(distlist, leaflist))     # List of (dist, leaf2) tuples
        dist_leaf_list.sort()                              # Sort on distance (first item in tuple)
        maxdist = dist_leaf_list[n_neighbors - 1][0]       # Maximum distance to include
        nearest_leaves = set()
        for dist, leaf in dist_leaf_list:
            if dist <= maxdist:
                nearest_leaves.add(leaf)
            else:
                break

        return nearest_leaves

    ###############################################################################################

    def find_mrca(self, leaves):
        """Finds Most Recent Common Ancestor for the provided set of leaves.
        MRCA for a leaf is leaf itself"""

        leafset = set(leaves)
        if not leafset <= self.leaves:
            leafstring = ", ".join(map(str, leafset))
            msg = "Some nodes in set are not part of tree: %s" % leafstring
            raise TreeError(msg)

        # special case: mrca for leaf is leaf itself
        if len(leafset) == 1:
            return next(iter(leafset))

        # pick random starting node among leafset, and find its parent node
        random_leaf = next(iter(leafset))
        parent = self.parent(random_leaf)

        # Walk down the tree from the initially picked node, until remkids include all of "leafset"
        while not leafset <= self.remotechildren_dict[parent]:
            parent = self.parent(parent)

        return parent

    ###############################################################################################

    def find_central_leaf(self, leaflist):
        """Finds central leaf for the provided list of leaves.
        Defined as having approximately equal distance to the two farthest leaves in leaflist"""

        nleaves = len(leaflist)
        #  Cluster has one member: return it (yeah, well...)
        if nleaves == 1:
            return leaflist[0]

        #  Cluster has only two members: Pick the leaf farthest from root
        if nleaves == 2:
            if self.nodedist(leaflist[0]) > self.nodedist(leaflist[1]):
                return leaflist[0]
            else:
                return leaflist[1]

        #  Cluster has more than two members:
        # Find two most distant leafs in leaflist (which "spread out" subtree).
        # Then find the leaf that has approximately the same distance to these two
        # (interpreted as being in a sense halfway between them...)
        basenode = self.find_mrca(leaflist)
        sub = self.subtree(basenode)
        (dist, leaf1, leaf2) = sub.diameter(return_leaves=True)
        smallest_diff = dist       # Pick value certain to be larger than all dist differences
        for leaf in leaflist:
            diff = abs(sub.nodedist(leaf, leaf1) - sub.nodedist(leaf, leaf2))
            if diff < smallest_diff:
                central_leaf = leaf
                smallest_diff = diff
        return central_leaf

    ###############################################################################################

    def find_common_leaf(self, leaflist):
        """Finds common leaf for the provided list of leaves.
        Defined as having the smallest average distance to remaining leaves"""

        nleaves = len(leaflist)
        #  Cluster has one member: return it (yeah, well...)
        if nleaves == 1:
            return leaflist[0]

        #  Cluster has two members: Pick the leaf farthest from root
        if nleaves == 2:
            if self.nodedist(leaflist[0]) > self.nodedist(leaflist[1]):
                return leaflist[0]
            else:
                return leaflist[1]

        # Cluster has more than two members:
        # Find leaf having the smallest average distance to remaining leaves
        dlist = []
        for leaf1 in leaflist:
            d = 0.0
            for leaf2 in leaflist:
                if leaf1 != leaf2:
                    d += self.nodedist(leaf1, leaf2)
            avedist = d / ( nleaves - 1 )
            dlist.append((avedist, leaf1))
        dlist.sort()
        typical_leaf = dlist[0][1]
        return typical_leaf

    ###############################################################################################

    def findbasenode(self, leafset):
        """Finds node that is at the base of all leaves in leafset."""

        # Specifically, the provided leafset is assumed to form one part of a bipartition
        # that is compatible with the tree. That bipartition corresponds to an internal
        # branch delimited by two internal nodes. This function determines which of these
        # nodes is nearer to the provided leafset. This could perhaps be described as
        # the unrooted variety of an MRCA

        # NOTE: this approach does not work when leafset is part of polytomy (i.e., basenode
        # has other descendants in addition to leafset. Think hard about how to fix this now...

        # Conversion to ensure comparison with sets can be done.
        # Could I assume set() without conversion?
        leafset = set(leafset)

        # Take care of special case where leafset equals all leafs in tree
        if leafset == self.leaves:
            return self.root

        # If root is at bifurcation: start by checking whether one of its kids match
        # If not: remove root from further investigation to avoid trouble
        rootkids = list(self.children(self.root))
        if len(rootkids) == 2:
            if leafset == self.remote_children(rootkids[0]):
                return rootkids[0]
            elif leafset == self.remote_children(rootkids[1]):
                return rootkids[1]
            else:
                intnodes = self.intnodes - {self.root}
        else:
            intnodes = self.intnodes

        # Iterate over remaining possible pairs of internal nodes
        for node1 in intnodes:
            for node2 in self.children(node1):
                bipart2 = self.remote_children(node2)    # returns set(node) if node is a leaf
                bipart1 = self.leaves - bipart2

                # Return basenode when bipartition has been found
                if leafset == bipart1:
                    return node1
                if leafset==bipart2:
                    return node2

        # If we fell off for loops, then leafset is incompatible with tree: panic
        # Construct string listing all entries in leafset (for error message)
        leafstring = ""
        for leaf in leafset:
            leafstring += str(leaf)
            leafstring += ", "
        leafstring = leafstring[:-2]    # Remove trailing comma and blank

        msg = f"The following is not a monophyletic group in this tree:\n{leafstring}"
        raise TreeError(msg)

    ###############################################################################################

    def find_bipart_nodes(self, bipartition):
        """Given a Bipartition as input: return the two nodes on the current tree that delimit
        the branch corresponding to the bipartition (if present in tree)"""

        (bip1, bip2) = bipartition
        mrca = self.find_mrca(bip1)
        bip = bip1
        if mrca == self.root:
            mrca = self.find_mrca(bip2)
            bip = bip2
        if self.remotechildren_dict[mrca] != bip:
            raise TreeError(f"Bipartition not present in tree: \n{bipartition}")

        return (self.parent(mrca), mrca)

    ###############################################################################################

    @functools.lru_cache(maxsize=None)
    def nodedist(self,node1,node2=None):
        """Returns distance between node1 and node2 along tree (patristic distance)"""

        # Python note: make recursive to gain advantage of keeping intermediate nodedists in cache?

        # Special case when node1 == node2:
        if node1 == node2:
            return 0.0

        # Node2 defaults to root if not given
        if node2 is None:
            node2 = self.root

        # Local copies to speed up access
        root = self.root
        pdict = self.parent_dict
        tree = self.child_dict

        # Find path from node1 back to root Keep track of cumulated distances along the way
        child1 = node1
        node1_ancdist = {}
        node1_ancdist[node1] = 0.0
        while child1 != root:
            parent1 = pdict[child1]
            node1_ancdist[parent1] = node1_ancdist[child1] + tree[parent1][child1].length
            child1 = parent1

        # Find path from node2 back to node on node1's path. Keep track of cumulated distances
        child2 = node2
        cumdist2 = 0.0
        while child2 not in node1_ancdist:
            parent2 = pdict[child2]
            cumdist2 += tree[parent2][child2].length
            child2 = parent2

        # Compute combined distance
        nodedist = cumdist2 + node1_ancdist[child2]
        return nodedist

    ###############################################################################################

    def nodedistlist(self, node1, nodelist):
        """Returns list of distances from node1 to nodes in nodelist (same order as nodelist)"""

        self.build_dist_dict()
        distlist = []
        for node2 in nodelist:
            distlist.append(self.dist_dict[node1][node2])
        return distlist

    ###############################################################################################

    def nodedepth(self,node):
        """Returns depth of node: distance from furthest leaf-level to node"""

        # The depth of node N can be found from the depth of the root as:
        # depth_N = depth_root - dist(root, N)
        # First: Find rootdepth. Use cached value if present
        try:
            rootdepth = self.rootdepth
        except AttributeError:
            maxdist = 0.0
            for leaf in self.leaflist():
                rootdist = self.nodedist(leaf)
                if rootdist > maxdist:
                    maxdist = rootdist
            rootdepth = maxdist
            self.rootdepth = rootdepth

        # Second: Find root-dist of specified node, and compute and return depth of node
        depth = rootdepth - self.nodedist(node)
        return depth

    ###############################################################################################

    def nodepath_fromdict(self, node1, node2):
        """Returns path between node1 and node2 along tree, from preconstructed path_dict"""

        try:
            path = [node1]
            n = node1
            while n != node2:
                n = self.path_dict[n][node2]
                path.append(n)
            return path
        except AttributeError as err:
            msg = "The path dictionary has not been constructed. Cannot use nodepath_fromdict"
            raise TreeError(msg) from err

    ###############################################################################################

    def nodepath(self, node1, node2):
        """Returns path between node1 and node2 along tree."""

        # If path_dict exists, use method tuned for that data structure
        if self.path_dict is not None:
            return self.nodepath_fromdict(node1, node2)

        # Find path from node1 to root
        root = self.root
        node1path = [node1]
        child = node1
        while child != root:
            parent = self.parent(child)
            node1path.append(parent)
            child = parent

        # Find path from node2 to root (or to first node that is also on node1's path)
        node2path = [node2]
        child = node2
        while child not in node1path:
            parent = self.parent(child)
            node2path.append(parent)
            child = parent
        intersect = child        # child is now the intersection between the two paths

        # Clean up paths
        node2path = node2path[:-1]                  # Remove intersection from node2path
        intersect_index = node1path.index(intersect)
        node1path = node1path[:intersect_index+1]   # Remove all upstream of intersect in node1path

        # Merge paths to form full path
        while node2path:
            lastnode = node2path.pop()
            node1path.append(lastnode)

        # node1path now contains entire path between nodes and is returned
        return node1path

    ###############################################################################################

    def length(self):
        """Returns tree length (sum of all branch lengths)"""

        # Python note: Maybe add lru cache, but then clear cache if leaves are added or removed

        treelength = 0.0
        for node in self.intnodes:
            for child in self.children(node):
                treelength += self.child_dict[node][child].length

        return treelength

    ###############################################################################################

    def height(self):
        """Returns height of tree: Largest root-to-tip distance"""
        nodeset = self.leaves
        node1 = self.root
        maxdist = self.find_most_distant(node1, nodeset)[1]

        return maxdist

    ###############################################################################################

    def diameter(self, return_leaves=False):
        """Return diameter: longest leaf-leaf distance along tree.
        If return_leaves is True: Return tuple with (maxdist, Leaf1, Leaf2)"""

        # Find the two leaves having the largest pairwise distance. Neat, 2-step algorithm:
        # (1) Pick random leaf, leaf1, find longest path to other leaf, leaf2
        # (2) Starting at leaf2, find longest path, this is longest path in tree! (It's true...)

        # Local copy for faster lookup (necessary?)
        nodedist = self.nodedist

        # Step 1: pick random leaf (leaf1) and find longest path to other leaf (leaf2)
        leaf1 = next(iter(self.leaves))       # Get random leaf from set without popping it
        maxdist = 0.0
        for node2 in self.leaves:
            dist = nodedist(leaf1, node2)
            if dist > maxdist:
                maxdist = dist
                leaf2 = node2

        # Step 2: Find longest path starting at leaf2, this is longest path in tree
        maxdist = 0.0
        for node2 in self.leaves:
            dist = nodedist(leaf2, node2)
            if dist > maxdist:
                maxdist = dist
                leaf3 = node2

        # Return requested result
        # Python note: Bad idea to have varying return values. Decide on output format
        if return_leaves:
            return (maxdist, leaf2, leaf3)
        else:
            return maxdist

    ###############################################################################################

    def find_most_distant(self, node1, nodeset):
        """Finds node in nodeset that is most distant from node1"""

        # Local copy for faster lookups
        nodedist = self.nodedist

        maxdist = 0.0
        for node2 in nodeset:
            dist = nodedist(node1, node2)
            if dist > maxdist:
                maxdist = dist
                most_distant = node2

        return (most_distant, maxdist)

    ###############################################################################################

    def patristic_distdict(self):
        """Return nested dictionary giving all pairwise, patristic distances:
        dict[node1][node2] = patristic distance"""

        # Initialize empty nested 2D dictionary with sequence names as keys.
        distdict = dict.fromkeys(self.leaves)
        for name in self.leaves:
            distdict[name] = dict.fromkeys(self.leaves)

        # Fill dictionary with values
        for leaf1, leaf2 in itertools.combinations(self.leaves, 2):
            dist = self.nodedist(leaf1,leaf2)
            distdict[leaf1][leaf2] = dist
            distdict[leaf2][leaf1] = dist

        return distdict

    ###############################################################################################

    def average_pairdist(self, leaflist, return_median=False):
        """Return average or median pairwise, patristic distance between leaves in leaflist"""
        # Python note: better to return list of pair distances, which can then be averaged etc.

        distlist = []
        for leaf1, leaf2 in itertools.combinations(leaflist, 2):
            distlist.append(self.nodedist(leaf1, leaf2))
        if return_median:
            return statistics.median(distlist)
        else:
            return statistics.mean(distlist)

    ###############################################################################################

    def average_ancdist(self, leaflist, return_median=False):
        """Return average or median patristic distance from leaves to their MRCA"""
        # Python note: better to return list of pair distances, which can then be averaged etc.

        ancnode = self.find_mrca(leaflist)
        distlist = []
        for leaf in leaflist:
            distlist.append(self.nodedist(leaf, ancnode))
        if return_median:
            return statistics.median(distlist)
        else:
            return statistics.mean(distlist)

    ###############################################################################################

    def shuffle_leaf_names(self):
        """Shuffles the names of all leaves"""

        # NOTE: this is mostly useful for statistical (resampling-style) analysis of whether
        # patterns of clustered leaves are significant

        # Construct list of old names
        oldnames = list(self.leaves)

        # Construct list of new names by shuffling old names
        newnames = list(self.leaves)
        random.shuffle(newnames)

        # Construct temporary list of names, to avoid duplicate names during name-changing process
        tmpnames = ["t" + str(i) for i in range(len(oldnames))]

        # Change old names to temporary names
        for (oldname, newname) in zip(oldnames, tmpnames):
            self.rename_leaf(oldname, newname)

        # Change temporary names to new names
        for (oldname, newname) in zip(tmpnames, newnames):
            self.rename_leaf(oldname, newname)

        # Clear lru_caches (which cannot be edited manually)
        # self.remote_children.cache_clear()
        self.clear_caches()

    ###############################################################################################

    def cladegrep(self, pattern, minsize = 2):
        """Finds clades (monophyletic groups) where all leaves contain specified pattern"""

        # Note: Assumes meaningful placement of root (all clades taken as starting at intnode
        # and progressing downstream, away from root.)

        # Step 1: Find all clades > minsize where all members match pattern
        matching_clades = []

        # self.remote_children() is expensive computationally. Ensure that intnodes are visited
        # in reasonably rational order such that its cache is used as much as possible
        # This means: shallow internal nodes (closer to leaves) should go first
        sortedintnodes = self.sorted_intnodes(deepfirst=False)
        for parent in sortedintnodes:
            clade = self.remote_children(parent)
            matching_clades.append(clade)

            # Drop clade if too small
            if len(clade) < minsize:
                del matching_clades[-1]

            # Drop clade if any members do not match pattern
            else:
                for leaf in clade:
                    if leaf.find(pattern) == -1:
                        del matching_clades[-1]
                        break

        # Step #2: Remove clades that are subsets of other clades
        maximal_clades = matching_clades[:]  # Copy of list (can't change list while iterating)
        for clade1 in matching_clades:
            for clade2 in matching_clades:
                if clade1 < clade2:
                    maximal_clades.remove(clade1)
                    break   # No reason to compare clade1 to any more

        return maximal_clades

    ###############################################################################################

    def newick(self, printdist=True, printlabels=True, print_leaflabels=False,
                precision=6, labelfield="label", transdict=None):
        """Returns Newick format tree string representation of tree object"""

        # Distances and (internal branch) labels are printed unless user explicitly request no printing
        # Transdict is meant for printing Nexus files with translate blocks: In these cases
        # the leaf names should be replaced by the leaf number instead, so transdict is reverse
        # of info in translate block: {name:number} instead of {number:name}
        # NOTE: This could probably be done slightly faster by iteration (instead of recursion)
        # for instance by using a stack, but unlikely that this function will be heavily used...
        # NOTE 2: I am using getattr() on "labelfield" to allow run-time specification of what
        # field in Branchstruct to use as label. Perhaps this is an indication that Branchstruct
        # should be a dict in the first place (instead of a class, that only contains data)

        def append_children(parentnode, labelfield):
            """Recursive function that has main responsibility for building Newick tree string"""

            for child in self.children(parentnode):

                branchstruct = self.child_dict[parentnode][child]
                dist = branchstruct.length
                label = getattr(branchstruct, labelfield)

                if child in self.leaves:
                    if transdict:
                        treelist.append(transdict[child])     # Should transdict errors be caught here?
                    else:
                        treelist.append(child)
                    if label != "" and print_leaflabels:
                        treelist.append(" ")
                        treelist.append("{}".format(label))
                    if printdist:
                        treelist.append(":{num:.{prec}g}".format(num=dist, prec=precision))
                else:
                    treelist.append("(")
                    append_children(child, labelfield)
                    treelist.append(")")

                    if label != "" and printlabels:
                        treelist.append("{}".format(label))
                    if printdist:
                        treelist.append(":{num:.{prec}g}".format(num=dist, prec=precision))

                treelist.append(",")

            del treelist[-1]            # Remove last comma when no more siblings

        # EXECUTION STARTS HERE!
        # Newick tree is built left-to-right as list, and finally converted to string
        # Root requires special treatment, rest of tree managed by recursion
        root = self.root
        treelist = ["("]
        append_children(root, labelfield)
        treelist.append(");")
        treestring = "".join(treelist)
        return treestring

    ###############################################################################################

    def nexus(self, printdist=True, printlabels=True, print_leaflabels=False,
                precision=6, labelfield="label", translateblock=False):
        """Returns nexus format tree as a string"""

        # Construct header
        stringlist = ["#NEXUS\n\nbegin trees;\n"]

        # If translateblock is requested: add translateblock to stringlist
        if translateblock:
            transdict = self.transdict()
            stringlist.append(self.translateblock(transdict))

        # Add newick tree string
        stringlist.append("   tree nexus_tree = ")
        if translateblock:
            stringlist.append(self.newick(printdist, printlabels, print_leaflabels, precision, labelfield, transdict))
        else:
            stringlist.append(self.newick(printdist, printlabels, print_leaflabels, precision, labelfield))

        # Add footer
        stringlist.append("\nend;\n")

        return "".join(stringlist)

    ###############################################################################################

    def figtree(self, printdist=True, printlabels=True, print_leaflabels=False, precision=6,
                colorlist=None, color="0000FF"):
        """Returns figtree format tree as a string"""

        # Implementation note: Rudimentary and mostly for coloring.
        # Should I add figtree section at end with various settings?

        # Construct header
        header = "#NEXUS\nbegin taxa\n\tdimensions ntax={};\n".format(len(self.leaflist()))
        header += "\ttaxlabels\n\t"
        stringlist = [header]
        for leaf in self.leaflist():
            stringlist.append("\t{}".format(leaf))
            if leaf in colorlist:
                col = color
            else:
                col = "000000"   # default color is black
            stringlist.append("[&!color=#{}]\n".format(col))
        stringlist.append(";\nend;\n")
        stringlist.append("\nbegin trees;\n\ttree nexus_tree = ")

        # Add newick tree string
        stringlist.append(self.newick(printdist, printlabels, print_leaflabels, precision))

        # Add footer
        stringlist.append("\nend;\n")

        return "".join(stringlist)

    ###############################################################################################

    def bipdict(self, keep_remchild_dict = False):
        """Returns tree in the form of a "bipartition dictionary" """

        # Names of leaves on one side of a branch are represented as an immutable set
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree is represented as a dictionary where the keys are bipartitions
        # The values are Branchstructs
        bipartition_dict = {}

        # For each branch: find bipartition representation, add this and Branchstruct to list.
        # Remote kids of node most distant from root (or node itself) forms one part of bipartition
        # Note: if root has two kids, then the root bipartition is added twice
        # This will be dealt with below
        for child, child_remkids in self.remotechildren_dict.items():
            if child != self.root:
                parent = self.parent(child)
                bipartition = Bipartition(child_remkids, self.frozenset_leaves,
                                          self.sorted_leaf_list, self.leaf2index)
                if self.interner:
                    bipartition = self.interner.intern_bipart(bipartition)
                origbranch = self.child_dict[parent][child]
                bipartition_dict[bipartition] = origbranch.copy()

        # If root is attached to exactly two nodes, then two branches correspond to the same
        # bipartition. Clean up by collapsing two branches (add lengths, compare labels)
        rootkids = self.children(self.root)
        if len(rootkids) == 2:
            kid1, kid2 = rootkids
            bipart1 = self.remotechildren_dict[kid1]
            rootbip = Bipartition(bipart1, self.frozenset_leaves,
                                      self.sorted_leaf_list, self.leaf2index)
            if self.interner:
                rootbip = self.interner.intern_bipart(rootbip)

            # Create new collapsed branch, sum distances from both kids
            combined_len = (self.child_dict[self.root][kid1].length +
                                                self.child_dict[self.root][kid2].length)
            bipartition_dict[rootbip] = Branchstruct(combined_len)

            # Deal with labels intelligently
            lab1 = self.child_dict[self.root][kid1].label
            lab2 = self.child_dict[self.root][kid2].label
            if (lab1 is not None) and (lab2 is None):
                lab = lab1
            elif (lab1 is None) and (lab2 is not None):
                lab = lab2
            else:
                lab = lab1
            bipartition_dict[rootbip].label = lab

        # Python note: to save memory. Maybe this should be dealt with centrally?
        if not keep_remchild_dict:
            self._remotechildren_dict = None

        return bipartition_dict

    ###############################################################################################

    def cladedict(self, keep_remchild_dict = False):
        """Returns tree in the form of a "clade dictionary" """

        # Names of leaves in clade are represented as an immutable set
        # The entire tree is represented as a dictionary where the keys are clades
        # The values are Nodestructs
        clade_dict = {}

        # For each node: find clade representation, add this and Nodestruct to list.
        for node, node_remkids in self.remotechildren_dict.items():
            clade = Clade(node_remkids, self.frozenset_leaves,
                          self.sorted_leaf_list, self.leaf2index)
            if self.interner:
                clade = self.interner.intern_clade(clade)
            nodedepth = self.nodedepthdict[node]
            clade_dict[clade] = Nodestruct(nodedepth)

        # Python note: to save memory. Maybe this should be dealt with centrally?
        if not keep_remchild_dict:
            self._remotechildren_dict = None

        return clade_dict

    ###############################################################################################

    def topology(self):
        """Returns set of sets of sets representation of topology ("naked bipdict")"""

        # Names of leaves on one side of a branch are represented as an immutable set.
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree topology is represented as a set of bipartitions
        # This is essentially a naked version of a bipdict

        # DEPRECATED: use either topology_bipart or topology_clade
        return self.topology_bipart

    ###############################################################################################

    def rootbip(self):
        """For a tree rooted at a bifurcation: returns a tuple giving the following information
        about the bipartition on which the root is located:
                (Bipartition, leafset1, blen1, leafset2, blen2)
        where leafset1 and leafset2 are the two halves of the bipartition, and blen1 and blen2
        are the lengths of the branches leading from the root to their two basal nodes"""

        rootkids = list(self.children(self.root))

        # If rooted at multifurcation: raise error
        # Python note: rethink logic here...
        if len(rootkids) > 2:
            msg = "Input tree rooted at multifurcation - not possible to assign root to bipartition"
            raise TreeError(msg)

        leafset1 = self.remote_children(rootkids[0])
        leafset2 = self.remote_children(rootkids[1])
        blen1 = self.nodedist(self.root, rootkids[0])
        blen2 = self.nodedist(self.root, rootkids[1])
        rootbip = Bipartition(leafset1, self.frozenset_leaves, self.sorted_leaf_list, self.leaf2index)

        return rootbip, leafset1, blen1, leafset2, blen2

    ###############################################################################################

    def check_bip_compatibility(self, bipart):
        """Checks the compatibility between bipartition and tree.
        Returns tuple of: is_present, is_compatible, insert_tuple
                where insert_tuple = None or (parentnode, childmovelist)
        is_present:
            True if bipartition is already present in tree. Implies "is_compatible = True"
        is_compatible:
            True if bipartition is compatible with tree. "is_present" can be True or False
        insert_tuple:
            If is_compatible: Tuple of (parentnode, childmovelist) parameters for insert_node
            If not is_compatible: None
        """
        bip1, bip2 = bipart
        if (bip1 | bip2) != self.leaves:
            raise TreeError("Bipartition has different set of leaves than tree")

        # Special case of startree: bip not present, but compatible
        # Return right away to avoid rest of function being indented...
        if len(self.intnodes) == 1:
            is_present = False
            is_compatible = True
            insert_tuple = (self.root, bip1)
            return is_present, is_compatible, insert_tuple

        # Special case where bipartition corresponds to leaf branch: always compatible
        if len(bip1) == 1 or len(bip2) == 1:
            if len(bip1) == 1:
                leaf = next(iter(bip1))
            else:
                leaf = next(iter(bip2))
            parent = self.parent(leaf)
            is_present = True
            is_compatible = True
            insert_tuple = (parent, leaf)
            return is_present, is_compatible, insert_tuple

        # In all other cases: at least one part of bipartition will necessarily have root as its MRCA
        # (because its members are present on both sides of root).
        # It is the other part of the bipartition that can potentially be moved
        # (Occasionally both parts have root as MRCA, but this will be covered by the more general
        # solution below)

        # First: set bip1 to smallest bip (increases probability that we wont have to compute mrca(bip2))
        if len(bip1) > len(bip2):
            bip1,bip2 = bip2,bip1
        mrca1 = self.find_mrca(bip1)
        if mrca1 != self.root:
            insertpoint = mrca1
            active_bip = bip1
        else:
            mrca2 = self.find_mrca(bip2)
            insertpoint = mrca2
            active_bip = bip2

        # Determine which of insertpoint's children to move (namely all children
        # that are either in the active bipartition OR children whose descendants are)
        movelist = []
        moveable_descendants = []
        for child in self.children(insertpoint):
            if child in active_bip:
                movelist.append(child)
                moveable_descendants.append(child)
            else:
                remkids = self.remotechildren_dict[child]
                if remkids <= active_bip:
                    movelist.append(child)
                    moveable_descendants.extend(remkids)

        # If moveable_descendants != active_bip, then bipartition is not compatible with tree
        # (and not present)
        # Otherwise: if insertpoint is at bifurcation, then bipart is already present in tree
        # (and compatible)
        # Final possibility: bipart is compatible but not present
        if set(moveable_descendants) != active_bip:
            is_present = False
            is_compatible = False
            insert_tuple = None
        elif self.is_bifurcation(insertpoint):
            is_present = True
            is_compatible = True
            insert_tuple = (insertpoint, movelist)
        else:
            is_present = False
            is_compatible = True
            insert_tuple = (insertpoint, movelist)

        return is_present, is_compatible, insert_tuple

    ###############################################################################################

    def is_compatible_with(self, bipart):
        """Checks whether a given bipartition is compatible with the tree.
        Note: also returns True if bipartition is already in tree"""

        is_present, is_compatible, insert_tuple = self.check_bip_compatibility(bipart)
        return is_compatible

    ###############################################################################################

    def bipart_is_present(self, bipart):
        """Checks whether a given bipartition is present in tree"""

        is_present, is_compatible, insert_tuple = self.check_bip_compatibility(bipart)
        return is_present

    ###############################################################################################

    def is_resolved(self):
        """Checks whether tree is fully resolved (no polytomies)"""

        n_leaves = len(self.leaves)
        n_nodes = n_leaves + len(self.intnodes)

        if self.is_bifurcation(self.root):
            maxnodes = 2 * n_leaves - 1
        else:
            maxnodes = 2 * n_leaves - 2

        return(n_nodes == maxnodes)

    ###############################################################################################

    def resolve(self):
        """Randomly resolves multifurcating tree by by adding zero-length internal branches."""

        # Find nodes with > 2 children, add to list of nodes needing to be resolved
        unresolved_nodes = []
        for node in self.intnodes:
            numkids = len(self.child_dict[node])   # Note: not safe to use .children() method while
                                             # changing tree (cache will break)
            if numkids > 2:
                unresolved_nodes.append(node)

        # Keep adding extra internal nodes until there are no unresolved nodes left
        while unresolved_nodes:
            intnode1 = unresolved_nodes.pop()
            kids = self.child_dict[intnode1].keys()
            kids = list(kids)

            # Divide children into two random subsets
            # Note: subsets can contain mix of leaves and intnodes
            # Python note: is it just as random, and simpler, to do one at a time?
            subset1_size = random.randint(1, len(kids) - 1)
            subset1 = random.sample(kids, subset1_size)
            subset2 = set(kids) - set(subset1)
            subset2_size = len(subset2)

            # For each subset:
            #   if more than 1 member: insert branch and intnode2 between intnode1 and subset
            #   if more than 2 members: also add intnode2 to list of unresolved nodes
            # After this, intnode1 is resolved (two branches emanating from it)
            if subset1_size > 1:
                branchstruct = Branchstruct()
                intnode2 = self.insert_node(intnode1, subset1, branchstruct)
                if subset1_size > 2:
                    unresolved_nodes.append(intnode2)

            if subset2_size > 1:
                branchstruct = Branchstruct()
                intnode2 = self.insert_node(intnode1, subset2, branchstruct)
                if subset2_size > 2:
                    unresolved_nodes.append(intnode2)

    ###############################################################################################

    def set_branch_attribute(self, node1, node2, attrname, attrvalue):
        """Set the value of any branch attribute.
        attrname: Name of attribute (e.g., "length")
        attrvalue: Value of attribute (e.g. 0.153)"""

        branch = self.child_dict[node1][node2]
        setattr(branch, attrname, attrvalue)

        self.clear_caches()    # Python note: only lengthrelated caches actually - refactor

    ###############################################################################################

    def setlength(self, node1, node2, length):
        """Sets length of branch connecting node1 and node2"""

        # Python note: maybe deprecate and use set_branch_attribute for all such problems?

        if node1 == self.parent(node2):
            parent = node1
            child = node2
        elif node2 == self.parent(node1):
            parent = node2
            child = node1
        else:
            msg = "There is no branch connecting node {} and {}".format(node1, node2)
            raise TreeError(msg)

        self.child_dict[parent][child].length = length

        self.clear_caches()    # Python note: only lengthrelated caches actually - refactor

    ###############################################################################################

    def setlabel(self, node1, node2, label):
        """Sets label on branch connecting node1 and node2"""

        if node1 == self.parent(node2):
            parent = node1
            child = node2
        elif node2 == self.parent(node1):
            parent = node2
            child = node1
        else:
            msg = "There is no branch connecting node %s and %s" % (node1, node2)
            raise TreeError(msg)

        self.child_dict[parent][child].label = label

    ###############################################################################################

    def getlabel(self, node1, node2):
        """Gets label on branch connecting node1 and node2"""

        if node1 == self.parent(node2):
            parent = node1
            child = node2
        elif node2 == self.parent(node1):
            parent = node2
            child = node1
        else:
            msg = "There is no branch connecting node %s and %s" % (node1, node2)
            raise TreeError(msg)

        return self.child_dict[parent][child].label

    ###############################################################################################

    def set_nodeid_labels(self):
        """Sets labels to be the same as the child node ID
        Allows use of e.g. Figtree to show nodeIDs as nodelabels"""

        for parent in self.intnodes:
            for kid in self.children(parent):
                self.setlabel(parent, kid, str(kid))

    ###############################################################################################

    def subtree(self, basenode, return_basalbranch=False):
        """Returns subtree rooted at basenode as Tree object"""

        # Note: rooting matters!
        # Note 2: basenode may be leaf!

        if return_basalbranch:
            if basenode == self.root:
                msg = "Can not return branch below root node"
                raise TreeError(msg)
            parent = self.parent(basenode)
            basalbranch = self.child_dict[parent][basenode]
            basalbranchcopy = basalbranch.copy()

        # Special case: basenode is leaf => subtree is minimal tree with two nodes (root and leaf)
        # Note: branchlength below leaf is used for the single branch in this tree
        if basenode in self.leaflist():
            other = Tree.from_string(f"({basenode});")
            blen = self.nodedist(self.parent(basenode), basenode)
            other.setlength(other.root, basenode, blen)

        # If basenode is internal: subtree has more than one leaf
        else:
            # Create empty Tree object. Transfer relevant subset of self's data structure to other
            other = Tree()
            other.child_dict = {}
            other.intnodes = {basenode}
            other.leaves = set()
            other.root = basenode
            curlevel = [basenode]
            while curlevel:
                nextlevel = []
                for parent in curlevel:
                    other.child_dict[parent] = {}
                    kids = self.children(parent)
                    for kid in kids:
                        other.child_dict[parent][kid] = copy.deepcopy(self.child_dict[parent][kid])
                    intnode_kids = kids & self.intnodes
                    other.intnodes.update(intnode_kids)
                    nextlevel.extend(intnode_kids)
                    leaf_kids = kids & self.leaves
                    other.leaves.update(leaf_kids)
                curlevel = nextlevel
            other.nodes = other.leaves | other.intnodes

        self.clear_caches()

        # Python note: possibly bad idea to have different possible returnvalues.
        # Simplify and deal with it at consumer end
        if return_basalbranch:
            return (other, basalbranchcopy)
        else:
            return other

    ###############################################################################################

    def graft(self, other, node1, node2=None, blen1=0, blen2=0, graftlabel=None,
                graft_with_other_root=False):
        """Graft other tree to self

        tree2 (other) intnodes will be renamed if names clash with those in tree1.
        node1: node in tree1 (self) below which tree2 (other) will be grafted. Cannot be root1
        node2: node in tree2 (other) below which tree2 will be attached (default is root of tree2)
        blen1: length of branch added to tree1 below graftpoint (lower of two newly created branches)
        blen2: length of branch above graft point and below tree2 (upper of two newly created branches)
        graftlabel: prepend value of "label" to leaf names on t2 (e.g: "graft_s1")
        graft_with_other_root: use root of other as graftpoint (i.e., do not add extra basal
                               branch between other.root and self.graftpoint)"""

        # Check that node1 is not root in tree1 (self)
        if node1 == self.root:
            raise TreeError("It is not possible to graft other tree below root-node: {}".format(node1))

        # Add new intnode on branch in tree1 where tree2 is to be grafted
        parent1 = self.parent(node1)
        branchstruct = Branchstruct(length=blen1)
        graftpoint = self.insert_node(parent1, [node1], branchstruct)

        # If node2 is not given: set to root node of tree2
        if node2 is None:
            node2 = other.root

        # If node2 is not root of tree2 then re-root on branch below node2.
        # After this, grafting can happen below root2
        elif node2 != other.root:
            other.deroot()
            node2parent = other.parent(node2)
            other.reroot(node2parent, node2)

        # Rename other's internal nodes if names clash with self's
        renameset = self.intnodes & other.intnodes
        if renameset:
            newnum = max(self.intnodes | other.intnodes)
            for oldnum in renameset:
                newnum += 1
                other.rename_intnode(oldnum, newnum)

        # prepend label to leaf names on grafted subtree if requested
        if graftlabel is not None:
            for oldname in other.leaflist():
                newname = "{}{}".format(graftlabel, oldname)
                other.rename_leaf(oldname, newname)

        # Update main data structure (self.child_dict dictionary) by merging with dict from other
        self.child_dict.update(other.child_dict)
        # Link subtree to graftpoint in self.child_dict
        self.child_dict[graftpoint][other.root] = Branchstruct(length=blen2)

        # Update look-up lists and caches
        self.nodes.update( other.nodes )
        self.intnodes.update( other.intnodes )
        self.leaves.update( other.leaves )

        # If requested: use other.root as graftpoint
        # (remove branch between graftpoint and other.root)
        # This is particularly useful if other consists of only root and single leaf
        if graft_with_other_root:
            self.remove_branch(graftpoint, other.root)

        # Reset caches and lists
        self.clear_caches()

    ###############################################################################################

    def cluster_n(self, nclust):
        """Divides tree into 'nclust' clusters based on distance from root.

           Returns tuple containing: list with sets of leafnames (one set per cluster)
                                     list of basenodes of clusters"""

        # Finds clusters by conceptually cutting across branches at height where there are nclust
        # groups downstream of cutpoint. This essentially finds equally spaced clusters.
        # Leafs that are upstream of cutpoint (closer to root) also form clusters of their own

        # Parses tree-dict, and for each "from" and "to" node, computes distance from root (height)
        # Saves this as two new atributes on Branchstruct: "parent_height" and "kid_height"

        # Also creates sorted list of (pheight, nkids) tuples containing internal node heights
        # and number of branches emanating from node
        # This information can be used to infer number of clusters when cutting above a given node
        # in the following way: If two branches emanate from root, then two clusters will be formed
        # by cutting above it. BUT, further out the tree: For each N branches emanating from an
        # internal node, N - 1 additional clusters will be formed by cutting above it
        # This is because one branch goes INTO the node (it was already part of one cluster)
        # (e.g., if two branches emanate from next node futher out, then ONE additional cluster
        # will be formed by cutting downstream of it)

        # Python note: Lazy implementation where I use .nodedist() function for each node in tree.
        # Could probably be sped up by using information already in .child_dict dict (thereby essentially
        # looking at lower branches only once)

        # PYTHON NOTE: simplify. Maybe always return list of leaf-sets.
        # Make sure that root dists are computed from bottom of tree first

        # Sanity check: There can not be more clusters than leaves
        nleaves = len(self.leaves)
        if nclust > nleaves:
            msg = "Requested {} clusters but tree only has {} leaves".format(nclust, nleaves)
            raise TreeError(msg)

        # Find nodeheights and number of emanating branches for all internal nodes
        nodeheightlist = []
        for parent, kid_dict in self.child_dict.items():
            pheight = self.nodedist(parent)
            nkids = len(kid_dict)
            nodeheightlist.append((pheight, nkids))
            for kid in kid_dict:
                kheight = self.nodedist(kid)
                kid_dict[kid].parent_height = pheight
                kid_dict[kid].kid_height = kheight
        nodeheightlist.sort()
        self.nodeheightlist = nodeheightlist

        # Find height cutoff corresponding to nclust groups:
        # walk along list of internal nodes until nclust groups reached
        cutno = nodeheightlist[0][1]   # Initialize cutno to number of branches emanating from root
        level = 0
        cutoff = 0.0
        while cutno < nclust:          # Increase level until >= required number of clusters
            level += 1
            cutoff = nodeheightlist[level][0]          # Height above which we want to cut branches
            cutno += (nodeheightlist[level][1] - 1)    # For each n additional branches cut,
                                                       # there will be n - 1 additional clusters

        # We now know cutoff:
        # height above which we can cut all branches in order to obtain >= nclust groups

        # For each branch in tree:
        # Check if cut at cutoff will cross this branch.
        # If so: add node above it to list of base nodes
        clusterlist = []                # list of sets of leaves (each set is one cluster)
        cluster_basenodes = []          # List of basenodes of clusters
        for parent, kid_dict in self.child_dict.items():
            for kid in kid_dict:
                if kid_dict[kid].parent_height <= cutoff < kid_dict[kid].kid_height:
                    cluster = self.remote_children(kid)
                    clusterlist.append(cluster)
                    cluster_basenodes.append(kid)       # NOTE: some basenodes may be leaves
                elif (kid in self.leaves) and kid_dict[kid].kid_height < cutoff:
                    cluster = set([kid])
                    clusterlist.append(cluster)
                    cluster_basenodes.append(kid)       # NOTE: some basenodes may be leaves

        return (clusterlist, cluster_basenodes)

    ###############################################################################################

    def cluster_cut(self, cutoff):
        """Divides tree into clusters by cutting across tree "cutoff" distance from root.
           Returns list containing sets of leafnames"""

        clusterlist = []            # List of leaf sets (each set is one cluster)
        cluster_basenodes = []      # List of basenodes of clusters
        cluster_leaves = set()       # set containing all leaves that are put in clusters
        # For each branch in tree: Find out if cut at cutoff will cross branch.
        # If branch will be cut: remote_children of kid form one cluster. Add set to clusterlist
        for parent in self.sorted_intnodes():
            parent_rootdist = self.nodedist(parent, self.root)
            for kid in self.children(parent):
                kid_rootdist = self.nodedist(kid, self.root)
                if parent_rootdist <= cutoff <= kid_rootdist:
                    cluster = self.remote_children(kid)
                    clusterlist.append(cluster)
                    cluster_basenodes.append(kid)
                    cluster_leaves.update(cluster)
        unclassified = self.leaves - cluster_leaves      # Set containing unclassified leaves
                                                         # (leaves that are below cutpoint)

        return (clusterlist, cluster_basenodes, unclassified)

    ###############################################################################################

    def insert_node(self, parent, childnodes, branchstruct):
        """Inserts an extra node between parent and children listed in childnodes list
        (so childnodes are now attached to newnode instead of parent).
        The branchstruct will be attached to the branch between parent and newnode.
        Branches to childnodes retain their original branchstructs.
        The node number of the new node is returned"""

        self._parent_dict = None

        if parent not in self.nodes:
            msg = f"Node {parent} does not exist"
            raise TreeError(msg)

        # Local copies for faster access
        tree = self.child_dict

        newnode = max(self.intnodes) + 1
        tree[newnode] = {}

        # Add new internal node as child of "parent"
        tree[parent][newnode] = branchstruct

        # Move childnodes from previous parent to new node
        for child in childnodes:
            tree[newnode][child] = tree[parent][child]
            del tree[parent][child]

        # Update self.intnodes and self.nodes to include new node.
        self.intnodes.add(newnode)
        self.nodes.add(newnode)

        self.clear_caches()

        return newnode

    ###############################################################################################

    def add_branch(self, bipart, branchstruct):
        """Adds branch represented by bipartition to unresolved tree."""

        # NOTE: huge overlap with TreeFromBiplist - should use this function there as well!

        # Sanity check: is bipartition compatible with tree?
        if not self.is_compatible_with(bipart):
            raise TreeError("Bipartition is not compatible with tree: %s" % bipart)

        part1, part2 = bipart
        if len(part1) > len(part2):
            part1,part2 = part2,part1    # Increases probability we wont have to compute mrca(part2)

        # Determine where to insert new node
        # In the special case of a star tree: add two new internal nodes, and move each half of
        # bipartition away from root (branch length will be divided equally between two branches)
        if len(self.intnodes) == 1:
            branchstruct1 = Branchstruct(length=branchstruct.length/2, label=branchstruct.label)
            branchstruct2 = Branchstruct(length=branchstruct.length/2, label=branchstruct.label)
            self.insert_node(self.root, part1, branchstruct1)
            self.insert_node(self.root, part2, branchstruct2)

        # In all other cases: one part of bipartition will necessarily have root as its MRCA
        #       (because its members are present on both sides of root).
        #       It is the other part of the bipartition that should be moved
        #       (the one where all members are on same side of root)
        else:
            mrca1 = self.find_mrca(part1)
            if mrca1 != self.root:      # If mrca2 is root, insert at mrca1, and move part1
                insertpoint = mrca1
                active_bip = part1
            else:                       # If mrca1 is root, insert at mrca2, and move part2
                mrca2 = self.find_mrca(part2)
                insertpoint = mrca2
                active_bip = part2

            # Determine which of insertpoint's children to move (namely all children
            # that are either in the active bipartition OR children whose descendants are)
            movelist = []
            for child in self.children(insertpoint):
                if child in active_bip:
                    movelist.append(child)
                elif (child not in self.leaves) and (self.remotechildren_dict[child] <= active_bip):
                    movelist.append(child)

            # Add branch at determined position
            self.insert_node(insertpoint, movelist, branchstruct)

        # Clear lru_caches (which cannot be edited manually)
        # self.remote_children.cache_clear()
        self.clear_caches()

    ###############################################################################################

    def remove_branch(self, node1, node2):
        """Removes branch connecting node1 and node2 (thereby possibly creating polytomy)
        Length of removed branch is distributed among descendant branches.
        This means tree length is conserved.
        Descendant nodes will be farther apart from each other, but closer to outside nodes."""

        if node1 == self.parent(node2):
            parent = node1
            child = node2
        elif node2 == self.parent(node1):
            parent = node2
            child = node1
        else:
            msg = "There is no branch connecting node %s and %s" % (node1, node2)
            raise TreeError(msg)

        if node1 in self.leaves or node2 in self.leaves:
            msg = "Attempting to remove external branch"
            raise TreeError(msg)

        # Move children of "child" so they are attached directly to "parent"
        # Keep track of length of removed branch, so it can be distributed among descendants
        lostlen = self.nodedist(parent,child)
        grandchildren = self.children(child)
        addlen = lostlen / len(grandchildren)
        for grandchild in grandchildren :
            self.child_dict[parent][grandchild] = Branchstruct(self.child_dict[child][grandchild].length + addlen,
                                                         self.child_dict[child][grandchild].label)

        # Delete "child" node and link from parent to child. Update intnodes and nodes
        del self.child_dict[child]
        del self.child_dict[parent][child]
        self.intnodes.remove(child)
        self.nodes = self.leaves | self.intnodes

        # Update _parent_dict
        del self._parent_dict[child]
        for grandchild in grandchildren:
            self._parent_dict[grandchild] = parent

        self.clear_caches()

    ###############################################################################################

    def remove_leaves(self, leaflist):
        """Removes leaves in list from tree, cleans up so remaining tree structure is sane"""

        for leaf in leaflist:
            self.remove_leaf(leaf)

    ###############################################################################################

    def remove_leaf(self, leaf):
        """Removes named leaf from tree, cleans up so remaining tree structure is sane"""

        parent = self.parent(leaf)
        childset = self.children(parent)
        root = self.root

        # If leaf is part of bifurcation AND is directly attached to root, then
        # the "other child" of the root must become the new root
        if (len(childset) == 2) and (leaf in self.children(root)):
            [child2] = childset - {leaf}                    # Remaining item is other child
            del self.child_dict[root]                             # Remove entry for old root
            self.intnodes.remove(root)
            self.nodes.remove(root)
            self.root = child2                              # child2 is new root
            del self._parent_dict[child2]                   # clean up parent_dict: Note: not lazy? change?

        # If leaf is part of bifurcation but NOT attached directly to root, then parent
        # must also be removed from tree, and the remaining child needs to be grafted
        # onto grandparent with proper cumulated distance
        # Keep branch label (if any) of the branch from grandparent to parent
        elif len(childset) == 2:
            [child2] = childset - {leaf}                    # Remaining item is other child
            child2dist = self.child_dict[parent][child2].length   # Remember dist to other child
            grandparent = self.parent(parent)

            # Add remaining child to grandparent
            # self.child_dict[grandparent][leaf2] = Branchstruct(self.child_dict[grandparent][parent].length,
            #                                              self.child_dict[grandparent][parent].label)
            self.child_dict[grandparent][child2] = self.child_dict[grandparent][parent]
            self.child_dict[grandparent][child2].length += child2dist   # Cumulated distance
            del self.child_dict[parent]                           # Delete parent and leaf
            del self.child_dict[grandparent][parent]              # Also remove pointer from gp to p
            del self._parent_dict[leaf]                      # Remove unused entries in parent_dict
            del self._parent_dict[parent]
            self._parent_dict[child2] = grandparent          # Update parent_dict for leaf2
            self.intnodes.remove(parent)
            self.nodes.remove(parent)

        # If leaf is part of multifurcation, then no special cleanup needed
        else:
            del self.child_dict[parent][leaf]
            del self._parent_dict[leaf]

        # Remove leaf entry from global leaflist. Update intnodeslist
        self.leaves.remove(leaf)
        self.nodes.remove(leaf)

        self.clear_caches()

    ###############################################################################################

    def add_leaf(self, parent, newleafname, branchstruct):
        """Adds new leaf to existing intnode ´parent´"""

        if parent not in self.intnodes:
            raise TreeError(f"Parent is not an existing internal node: {parent}")
        if newleafname in self.leaves:
            raise TreeError(f"Leaf already exists: {newleafname}")
        self.child_dict[parent][newleafname] = branchstruct
        self.nodes.add(newleafname)
        self.leaves.add(newleafname)

        self.clear_caches()

    ###############################################################################################

    def collapse_clade(self, leaflist, newname="clade"):
        """Replaces clade (leaves in leaflist) with single leaf.
        Branch length is set to average dist from basenode parent to leaves"""

        if len(leaflist) == 1:
            # Special case where there is only one leaf in leaflist:
            # Do not collapse anything, but change name to newname (?)
            oldname = leaflist.pop()
            self.rename_leaf(oldname, newname)
        else:
            # Find average distance from parent of basenode to leaves (median - use mean instead?)
            mrca = self.find_mrca(leaflist)
            mrca_parent = self.parent(mrca)
            avdist = self.average_ancdist(leaflist, return_median=True)
            avdist += self.nodedist(mrca, mrca_parent)

            # Remove all but one of the leaves in leaflist
            # (hackish way of keeping leaf node for subsequent renaming...)
            leaflist = list(leaflist)
            subleaflist = leaflist[1:]
            self.remove_leaves(subleaflist)

            # Rename remaining leaf to newname and set branch length to average dist
            self.rename_leaf(leaflist[0], newname)
            self.setlength(mrca_parent, newname, avdist)

        self.clear_caches()

    ###############################################################################################

    def nameprune(self, sep="_", keep_pattern=None):
        """Prune leaves based on name redundancy:
        Find subtrees where all leaves have same start of name (up to first "_")"""

        # A bit of a hack for very specific project...
        # I should find way to generalize (e.g. ask for pattern/regular expression to match)

        # Find namestarts (up to first occurrence of sep) that occur more than once among leaves
        seen = set()
        dups = set()
        for name in self.leaflist():
            if not keep_pattern or keep_pattern not in name:
                namestart = name.split(sep)[0]
                if namestart in seen:
                    dups.add(namestart)
                else:
                    seen.add(namestart)

        # Find clades where all members contain one of the duplicated name starts
        remlist = []
        for dupname in dups:
            remlist.extend(self.cladegrep(dupname))

        # From each clade: remove all but one of the leaves with matching name starts
        for nameset in remlist:
            nameset.pop()                   # Randomly remove one member (will be kept)
            for name in nameset:
                self.remove_leaf(name)

    ###############################################################################################

    def numberprune(self, nkeep, keeplist=None, keep_common_leaves=False,
                    keep_most_distant=False, return_leaves=False, enforce_n = False):
        """Prune tree so 'nkeep' leaves remain, approximately evenly spaced over tree.

        "keeplist" can be used to specify leaves that _must_ be retained.
        'keep_common_leaves' requests preferential retainment of leaves with many neighbors
        (default is to keep leaves that are as equally spaced as possible)
        'keep_most_distant' requests that the two most distant leaves in tree
                    (which spread out the diameter) should be kept
        'return_leaves': return selected leaves, but do not actually prune tree
        'enforce_n' enforce exactly N leaves in pruned tree
                    (normally leaves in includelist and most distant are additional to N)"""

        keepset = set()
        if keeplist:
            keepset.update(keeplist)

        # If requested: Find and retain two most distant leaves in tree
        if keep_most_distant:
            (_, leaf1, leaf2) = self.diameter(return_leaves=True)
            keepset.update((leaf1, leaf2))

        # If enforce_n has not been requested:
        # Find N clusters in addition to any members of keeplist or the two most distant leaves
        if not enforce_n:
            clusters = self.cluster_n(nclust=nkeep)[0]   # list containing sets of leafnames

        # If enforce_n: Iteratively find N, N-1, ... clusters until total retained number of leave
        # (including those in keeplist etc) is == N
        else:
            n_clusters = nkeep
            found_n = False
            while not found_n:
                clusters = self.cluster_n(n_clusters)[0]
                coveredclusters = []
                for cluster in clusters:
                    intersection = cluster & keepset
                    if intersection:
                        coveredclusters.append(cluster)
                n_retained = len(keepset) - len(coveredclusters) + len(clusters)
                if n_retained == nkeep:
                    found_n = True
                else:
                    n_clusters -= 1

            # We have now identified:
            #    (1) The value of N that will result in the requested number of retained leaves
            #    (2) The clusters in which members of keepset are located.
            # Remove these clusters before proceeding
            for cluster in coveredclusters:
                clusters.remove(cluster)

        # For each cluster: Find one representative member, and add this to keepset
        for cluster in clusters:
            if keep_common_leaves:
                keep_leaf = self.find_common_leaf(list(cluster))
            else:
                keep_leaf = self.find_central_leaf(list(cluster))
            keepset.add(keep_leaf)

        self.clear_caches()

        # If requested: return selected leaves without pruning tree
        if return_leaves:
            return keepset

        # Otherwise: prune tree so only leaves in keepset are retained, by removing all others
        else:
            discardset = self.leaves - keepset
            self.remove_leaves(discardset)
            return None

    ###############################################################################################

    def prune_maxlen(self, nkeep, return_leaves=False):
        """Prune tree so remaining nkeep leaves spread out maximal percentage of branch length"""

        possible_branches = set()     # Possible basal branches for starting next path to leaf
                                      # (node1 is on path, and node2 is not)
        used_branches = set()         # Branches that are on the path
        keep_leaves = set()           # Leaves to keep in tree

        # Midpoint root to make initialisation simpler
        # (costly - should rewrite algorithm to start anywhere.
        # On the other hand I would need diameter anyway...)
        self.rootmid()

        # Place central data structures and functions in local namespace for faster lookup
        nodedist = self.nodedist
        nodepath = self.nodepath
        remote_children = self.remote_children
        children = self.children

        # Initialise: add two branches emanating from root to list of possible starting branches
        # Note: this only works due to midpoint rooting
        # (which ensures root will be on the first, longest path between two leafs)
        rootkids = self.children( self.root )
        for child in rootkids:
            possible_branches.add( (self.root, child) )

        # Until we have added nkeep leaves to path:
        # find longest newpath from existing path to leaf, add to path
        while len(keep_leaves) < nkeep:

            # Among possible starting branches:
            # find the one having the max possible distance to a remote child
            maxdist = 0.0
            for (parent, child) in possible_branches:
                for leaf in remote_children(child):
                    dist = nodedist(parent,leaf)
                    if dist > maxdist:
                        maxdist = dist
                        node1, node2, keepleaf = parent, child, leaf

            # Add the found leaf to list of leaves.
            # Remove the basal branch that was used from possible starting branches
            keep_leaves.add( keepleaf )
            possible_branches = possible_branches - { (node1, node2) }

            # Update possible_branches and used branches based on newly added path
            newpath = nodepath( node1, keepleaf )
            for i in range( len(newpath) - 1 ):
                parent, child1 = newpath[i], newpath[i+1]
                used_branches.add( ( parent, child1 ) )
                otherkids = children(parent) - {child1}
                for child2 in otherkids:
                    if (parent, child2) not in used_branches:
                        possible_branches.add( (parent, child2) )

        self.clear_caches()

        # If requested: return selected leaves without pruning tree
        # Otherwise: prune tree so only leaves in keepset are retained
        if return_leaves:
            return keep_leaves
        else:
            discard_leaves = self.leaves - keep_leaves
            self.remove_leaves(discard_leaves)

    ###############################################################################################

    def transname(self, namefile):
        """Translate all leaf names using oldname/newname pairs in namefile"""

        with open(namefile, "r") as transfile:
            transdict = {}
            for line in transfile:
                words = line.split()
                oldname = words[0]
                newname = words[1]
                transdict[oldname] = newname

        # Create copy of original set of leaf names to avoid iterating over set while changing it
        orignames = copy.copy(self.leaves)
        for oldname in orignames:
            newname = transdict[oldname]
            self.rename_leaf(oldname, newname)

    ###############################################################################################

    def rename_leaf(self, oldname, newname, fixdups=False):
        """Changes name of one leaf. Automatically fixes duplicates if requested"""

        if oldname not in self.leaves:
            msg = "Leaf %s does not exist" % oldname
            raise TreeError(msg)

        if newname in self.leaves:
            if not fixdups:
                msg = "Attempted to create duplicate leafname: %s" % newname
                raise TreeError(msg)
            i = 1
            fixedname = newname + "_" + str(i)
            while fixedname in self.leaves:
                i += 1
                fixedname = newname + "_" + str(i)
            newname = fixedname

        parent = self.parent(oldname)
        self._parent_dict = None
        self.child_dict[parent][newname] = self.child_dict[parent][oldname]
        del self.child_dict[parent][oldname]

        # Update self.leaves and self.nodes
        self.leaves.add(newname)
        self.leaves.remove(oldname)
        self.nodes.add(newname)
        self.nodes.remove(oldname)

        # # Update self.parent_dict if it exists:
        # if self.parent_dict is not None:
        #     self.parent_dict[newname] = self.parent_dict[oldname]
        #     del self.parent_dict[oldname]

        self.clear_caches()

    ###############################################################################################

    def rename_intnode(self, oldnum, newnum):
        """Changes number of one internal node"""

        if oldnum not in self.intnodes:
            msg = "Internal node {} does not exist".format(oldnum)
            raise TreeError(msg)

        if newnum in self.intnodes:
            msg = f"There is already an internal node with the number {newnum}"
            raise TreeError(msg)

        # Make a note of original's parent and children
        kidlist = self.children(oldnum)
        parent = self.parent(oldnum)            # Will be None if oldnum is root

        # Update main data structure (child list in self.child_dict)
        self.child_dict[newnum] = {}
        for child in kidlist:
            self.child_dict[newnum][child] = self.child_dict[oldnum][child]
        del self.child_dict[oldnum]
        if oldnum == self.root:
            self.root = newnum
        else:
            self.child_dict[parent][newnum] = self.child_dict[parent][oldnum]
            del self.child_dict[parent][oldnum]

        # Update look-up lists, caches, and root-marker if relevant
        self.nodes.add(newnum)
        self.nodes.remove(oldnum)
        self.intnodes.add(newnum)
        self.intnodes.remove(oldnum)
        self._parent_dict = None
        # if parent is not None:
        #     self.parent_dict[newnum] = self.parent_dict[oldnum]
        #     del self.parent_dict[oldnum]
        # for child in kidlist:
        #     self.parent_dict[child] = newnum

        self.clear_caches()

    ###############################################################################################

    def treedist_RF(self, other, normalise=False, rooted=False, return_details=False):
        """Compute symmetric tree distance (Robinson Foulds) between self and other tree.
        normalise: divide RF distance by the total number of bipartitions in the two trees
        rooted: take position of root into account.
        return_details: return list including intermediate values also:
            [treedist, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2,
             tree1_unique_biparts, tree2_unique_biparts, shared_biparts]
            where the last three are sets of bipartitions"""

        # Python note: the rooted measure is found by adding an extra leaf to the root node
        # before computing distance. This is the same as counting clades for internal branches
        # on rooted tree

        # Check that trees are comparable (have identical leaf sets)
        # Python note: maybe allow automatically discarding different leaves
        if self.leaves != other.leaves:
            raise TreeError("Can't compute treedist: two trees have different leaf sets")

        # Find set of bipartitions in each tree
        # If rooted measure requested: first copy trees, then add extra leaf to root nodes
        if rooted:
            tree1 = self.copy_treeobject()
            tree1.add_leaf(tree1.root, "root", Branchstruct())
            tree2 = other.copy_treeobject()
            tree2.add_leaf(tree2.root, "root", Branchstruct())
        else:
            tree1 = self
            tree2 = other
        tree1_biparts = tree1.topology()
        tree2_biparts = tree2.topology()

        # Find bipartitions unique to tree1 and tree2
        tree1_unique_biparts = tree1_biparts - tree2_biparts
        tree2_unique_biparts = tree2_biparts - tree1_biparts

        # Find shared bipartitions
        shared_biparts = tree1_biparts & tree2_biparts

        # Compute distance
        n_shared = len(shared_biparts) - len(tree1.leaves)     # Only internal branches counts!!!
        n_uniq1 = len(tree1_unique_biparts)
        n_uniq2 = len(tree2_unique_biparts)
        n_bip1 = len(tree1_biparts) - len(tree1.leaves)        # Only internal branches counts!!!
        n_bip2 = len(tree2_biparts) - len(tree2.leaves)       # Only internal branches counts!!!
        treedist = n_uniq1 + n_uniq2
        if normalise:
            treedist = treedist / (n_bip1 + n_bip2)

        if return_details:
            return [treedist, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2,
                    tree1_unique_biparts, tree2_unique_biparts, shared_biparts]
        else:
            return treedist

    ###############################################################################################

    def _pathdist_as_ndarray(self):
        """Utillity function for treedist_pathdiff: return pathdiff matrix as vector:
        For each pair of leaves: count number of edges on path between them,
        return flattened version of pairwise matrix (in alphabetical leaf-pair order)
        Note: traversing root counts as one edge (so I deroot treecopy before counting)"""

        def dist_iterator(tree, namepairs):
            for n1, n2 in namepairs:
                yield len(tree.nodepath(n1,n2)) - 1

        tree = self.copy_treeobject()
        tree.deroot()
        leafnames = tree.leaflist()     # Sorted alphabetically
        namepairs = itertools.combinations(leafnames, 2)
        nleaves = len(leafnames)
        npairs = nleaves * (nleaves - 1) // 2
        distarray = np.fromiter(dist_iterator(tree, namepairs), dtype=float, count=npairs)
        return distarray

    ###############################################################################################

    def treedist_pathdiff(self, other):
        """Compute path difference tree-distance between self and other:
        Euclidean distance between nodepath-dist matrices considered as vectors.
        Measure described in M.A. Steel, D. Penny, Syst. Biol. 42 (1993) 126–141
        """

        self_distvec = self._pathdist_as_ndarray()
        other_distvec = other._pathdist_as_ndarray()

        return np.linalg.norm(self_distvec - other_distvec)

    ###############################################################################################

    def treedist(self, other, normalise=True, verbose=False):
        """Deprecated: Use treedist_RF instead.
        Compute symmetric tree distance (Robinson Foulds) between self and other tree.
        Normalised measure returned by default"""

        # Check that trees are comparable (have identical leaf sets)
        if self.leaves != other.leaves:
            raise TreeError("Can't compute treedist: two trees have different leaf sets")

        # Find set of bipartitions in each tree
        # Recall that: Names of leafs on one side of a branch are represented as an immutable set.
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree topology is represented as a set of bipartitions
        tree1_biparts = self.topology()
        tree2_biparts = other.topology()

        # Find bipartitions unique to tree1 and tree2 (using set arithemtic)
        tree1_unique_biparts = tree1_biparts - tree2_biparts
        tree2_unique_biparts = tree2_biparts - tree1_biparts

        # Find shared bipartitions
        shared_biparts = tree1_biparts & tree2_biparts

        # Compute distance
        n_shared = len(shared_biparts) - len(self.leaves)     # Only internal branches counts!!!
        n_uniq1 = len(tree1_unique_biparts)
        n_uniq2 = len(tree2_unique_biparts)
        n_bip1 = len(tree1_biparts) - len(self.leaves)        # Only internal branches counts!!!
        n_bip2 = len(tree2_biparts) - len(other.leaves)       # Only internal branches counts!!!
        symdif = n_uniq1 + n_uniq2
        symdif_norm = symdif / (n_bip1 + n_bip2)

        # Return requested values
        # Python note: really should not have different possible returns... Refactor
        if verbose:
            return symdif, symdif_norm, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2
        elif normalise:
            return symdif_norm
        else:
            return symdif

    ###############################################################################################

    def treesim(self, other, verbose=False):
        """Compute normalised symmetric similarity between self and other tree"""

        symdif, symdif_norm, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2 = self.child_dictdist(other, verbose=True)
        symsim_norm = 1.0 - symdif_norm

        if verbose:
            return symdif, symdif_norm, symsim_norm, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2
        else:
            return symsim_norm

    ###############################################################################################

    def has_same_root(self, other):
        """Compares two trees. Returns True if topologies are same and rooted in same place"""

        if self.topology() != other.topology():
            raise TreeError("Tree topologies are different: Rootings can not be compared")
        else:
            t1partition = set()   # Set of sets of kid's remotechildren
            for kid in self.children(self.root):
                t1partition.add(frozenset(self.remote_children(kid)))
            t2partition = set()
            for kid in other.children(other.root):
                t2partition.add(frozenset(other.remote_children(kid)))
            if t1partition == t2partition:
                return True
            else:
                return False

    ###############################################################################################

    def deroot(self):
        """If root is at bifurcation: remove root node, connect adjacent nodes"""
        root = self.root
        rootkids = self.children(root)
        if len(rootkids) == 2:      # If root is at bifurcation
            kid1, kid2 = rootkids

            # Length of new, combined branch is sum of distances
            distsum = self.child_dict[root][kid1].length + self.child_dict[root][kid2].length

            # Deal with labels semi-intelligently
            # If only one label set: use that. Otherwise pick lab1
            lab1 = self.child_dict[root][kid1].label
            lab2 = self.child_dict[root][kid2].label
            if lab1 != "" and lab2 == "":
                lab = lab1
            elif lab1 == ""  and lab2 != "":
                lab = lab2
            else:
                lab = lab1                  # If agree pick 1, if disagree: pick 1 randomly

            # Python note: should handle situation where Branchstruct has additional attributes
            # Use introspection to find attributes and combine intelligently?
            if kid1 in self.intnodes:
                self.child_dict[kid1][kid2] = Branchstruct(length = distsum, label = lab)
                self.root = kid1
            elif kid2 in self.intnodes:
                self.child_dict[kid2][kid1] = Branchstruct(length = distsum, label = lab)
                self.root = kid2
            else:
                raise TreeError("Cannot deroot tree: only leaves are present")

            # remove previous root
            del self.child_dict[root]

            # update intnode and node attributes and clear caches, which are now unreliable
            self.intnodes = set(self.child_dict.keys())
            self.nodes = self.intnodes | self.leaves
            self.clear_caches()

    ###############################################################################################

    def reroot(self, node1, node2=None, polytomy=False, node1dist=0.0):
        """Places new root on branch between node1 and node2, node1dist from node1"""

        # If tree is to be rooted at basal polytomy, then node1 (base of outgroup) is the new root
        if polytomy:
            newroot = node1

        # If polytomy not requested: new root must be inserted on branch between node1 and node2:
        # Determine which node is parent of other, determine branchlength and label for branch
        # between nodes 1 and 2, figure out what the new branch lengths will be after splitting
        # the branch, and finally insert new node. Bail out if the nodes are not neighbors
        else:
            if node2 is None:
                msg = "Need to specify node2 to reroot() method when rooting at bifurcation"
                raise TreeError(msg)
            if node1dist > self.nodedist(node1,node2):
                msg = ("Parameter node1dist too large:\n"
                       + f"    node1dist = {node1dist} > dist(node1, node2) = {self.nodedist(node1,node2)}")
                raise TreeError(msg)
            if node1 == self.parent(node2):
                parent = node1
                child = node2
                parent_to_root_dist = node1dist
                root_to_child_dist = self.nodedist(node1,node2) - node1dist
            elif node2 == self.parent(node1):
                parent = node2
                child = node1
                parent_to_root_dist = self.nodedist(node1,node2) - node1dist
                root_to_child_dist = node1dist
            else:
                msg = "Node {} and {} are not neighbors in tree".format(node1, node2)
                raise TreeError(msg)

            # New branch (from parent to newroot) will inherit all attributes from original branch
            # (from parent to child), except for length, which is split between the two branches
            # Special case: if child is leaf, then original branch will have no branch support
            # In this case set label to "1.0" (maybe only do this if rest of labels are branch support??)
            newbranch = self.child_dict[parent][child].copy()
            newbranch.length = parent_to_root_dist
            if newbranch.label == "":
                newbranch.label = "1.0"
            newroot = self.insert_node(parent, [child], newbranch)
            self.child_dict[newroot][child].length = root_to_child_dist

        # Things that were already downstream of newroot do not need to be moved, but things that
        # were previously upstream need to be moved downstream, which is done by reversing the
        # links on the direct path going back from newroot to oldroot
        oldroot = self.root
        path_to_old = self.nodepath(newroot, oldroot)
        for i in range( len(path_to_old) - 1 ):
            newparent, oldparent = path_to_old[i], path_to_old[i+1]
            self.child_dict[newparent][oldparent] = self.child_dict[oldparent][newparent]
            del self.child_dict[oldparent][newparent]
            # self.parent_dict[oldparent] = newparent

        self.clear_caches()

        # Update root info:
        self.root = newroot

    ###############################################################################################

    def rootmid(self):
        """Performs midpoint rooting of tree"""

        # Sanity check: if tree has zero length, then midpoint rooting is not possible
        if self.length() == 0.0:
            raise TreeError("All branch lengths are zero - midpoint rooting not possible")

        # Remove previous root if present at bifurcation
        self.deroot()

        # Find the two leaves having the largest pairwise distance.
        (maxdist, leaf1, leaf2) = self.diameter(return_leaves = True)
        midway = maxdist/2.0

        # Get path between leaf1 and leaf2
        path = self.nodepath(leaf1, leaf2)

        # Find the branch that contains the midpoint of the tree:
        # Work backwards through path, stop when cumulated branch length exceeds midway
        cumdist = 0.0
        while cumdist < midway:
            (node1, node2) = (path[-1], path[-2])
            cumdist += self.nodedist(node1, node2)
            path.pop()

        # Place root on current branch, correct distance from one end
        node2dist = cumdist - midway
        self.reroot(node2, node1, False, node2dist)

    ###############################################################################################

    def rootminvar(self):
        """Performs minimum variance rooting of tree"""

        # Based on results in the paper:
        # "Minimum variance rooting of phylogenetic trees and implications for
        # species tree reconstruction", Uyen Mai, Erfan Sayyari, Siavash Mirarab

        # Sanity check: if tree has zero length, then minimum variance rooting is not possible
        if self.length() == 0.0:
            raise TreeError("All branch lengths are zero - minimum variance rooting not possible")

        # Find all pairwise distances between nodes
        # Store all intnode-to-leaf distances in numpy 2D array (matrix)
        # Each row has leaf-dists for one intnode
        # Python note: This can be optimized
        self.build_dist_dict()
        intnodes = self.sorted_intnodes()
        leaves = self.leaflist()
        nodes = intnodes + leaves
        distmat = np.empty(shape=(len(nodes), len(leaves)))
        for i,node in enumerate(nodes):
            for j,leaf in enumerate(leaves):
                distmat[i,j] = self.dist_dict[node][leaf]

        # Compute variance of root-to-tip distance for all nodes
        # Python note: perhaps do this just-in-time when needed
        # (could I omit some computations for leaves for instance?)
        distvar = distmat.var(axis=1)             # axis=1: variance of each row in matrix

        # For each branch in tree:
        # find local minimum variance point using equation 7 in Mai paper (mimimum of parabola)
        # Location of point is stored as triplet:
        # parentnode, childnode, distance of point from parent on branch
        # Keep track of overall minimum variance and corresponding location

        # Initialise overall minimum to be distance zero from first intnode on one of its branches
        minvar = distvar[0]
        minparent = nodes[0]
        minchild = self.children(minparent).pop()
        minpardist = 0.0

        # Iterate over branches in tree,
        # compute minimum variance point, and update overall minimum if relevant
        # Names of intermediate variables have been set to match those in Mai paper
        for u,parent in enumerate(intnodes):
            for child in self.children(parent):
                v = nodes.index(child)
                STu = distmat[u].sum()           # Sum of distances to all leaves from parent
                SIv = 0.0                        # Sum of distances to remote_children from child
                for remchild in self.remote_children(child):
                    SIv += self.dist_dict[child][remchild]
                v_size = len(self.remote_children(child))
                ev = self.dist_dict[parent][child]
                n = len(leaves)
                alpha = ( 2*STu - 4*(SIv + v_size*ev) ) / n
                beta = 1 - 2*v_size/n
                a = 1 - beta**2
                b = alpha - (2*STu*beta / n)
                c = distvar[u]
                if a != 0:
                    x = -b / (2*a)                  # Equation 4. x = distance from parent on branch
                else:
                    x = -math.inf                   # Indicator value to show no root in polynomium
                if 0 < x < ev:
                    xvar = a*x**2 + b*x + c
                elif distvar[u] < distvar[v]:
                    x = 0.0
                    xvar = distvar[u]
                else:
                    x = ev
                    xvar = distvar[v]
                if xvar < minvar:
                    minvar = xvar
                    minparent = parent
                    minchild = child
                    minpardist = x

        # Reroot on global minimum variance point

        # If minparent is not root:
        # remove old root and reroot using the minparent and minpardist already found:
        if minparent != self.root:
            self.deroot()
            self.reroot(node1=minparent, node2=minchild, node1dist=minpardist)

        # If minparent IS root:
        # remove old root and reroot using the minparent and minpardist already found:
        else:
            # Minpardist == 0: Do nothing. Current root is minimal variance root
            if minpardist == 0.0:
                return
            # Minpardist != 0: Do something, depending on whether root is at bifurcation or not
            else:
                rootkids = self.children(self.root)
                # Bifurcation:
                # old root will be removed. Find new minparent (one of root's other children)
                if len(rootkids) == 2:
                    rootkids.remove(minchild)
                    minparent = rootkids.pop()    # Chose one of root's other children as minparent
                    minpardist = minpardist + self.dist_dict[self.root][minparent]
                    self.deroot()
                    self.reroot(node1=minparent, node2=minchild, node1dist=minpardist)
                    return

                # Multifurcation: old root will not be removed. Do not need to find new minparent
                else:
                    self.reroot(node1=minparent, node2=minchild, node1dist=minpardist)
                    return


    ###############################################################################################

    def rootout(self,outgroup, polytomy=False):
        """Roots tree on outgroup"""

        # Remove previous root if present at bifurcation
        self.deroot()

        # Find pair of internal nodes corresponding to ingroup:outgroup bipartition
        if type(outgroup) is str:       # if just a single name is given:
                                        # protect string from iteration to single letters
            outgroup = [outgroup]
        outgroup = frozenset(outgroup)
        ingroup = self.leaves - outgroup
        outbase = self.findbasenode(outgroup)
        inbase = self.findbasenode(ingroup)

        # If outgroup should form basal polytomy with ingroup: root on outbase
        if polytomy:
            self.reroot(outbase, inbase, polytomy=True)

        # If tree does not have branch lengths: skip computation of lengths
        elif self.length() == 0.0:
            self.reroot(outbase, inbase)

        # Else: compute where on branch to place root
        else:

            # Find longest base-leaf distance in outgroup:
            distances = set()
            for leaf in outgroup:
                leafdist = self.nodedist(leaf, outbase)
                distances.add(leafdist)
            max_out_dist = max(distances)

            # Find longest base-leaf distance in ingroup:
            distances = set()
            for leaf in ingroup:
                leafdist = self.nodedist(leaf, inbase)
                distances.add(leafdist)
            max_in_dist = max(distances)

            # Find length of branch separating ingroup and outgroup
            inoutdist = self.nodedist(outbase, inbase)

            # If possible, then put root exactly midway between extremes
            diameter = max_in_dist + max_out_dist + inoutdist
            midway = diameter/2
            really_close = 0.1

            # Case 1: outgroup has only one taxon
            if len(outgroup) == 1:
                # if possible place root at midpoint
                if midway < inoutdist:
                    outbasedist = midway
                # If not, place it real close to midpoint
                else:
                    outbasedist = (1 - really_close) * inoutdist

            # Case 2: outgroup has more than one taxon
            # If possible place root at midpoint
            elif (midway > max_in_dist) and (midway > max_out_dist):
                outbasedist = midway - max_out_dist     # This is distance from outbase to new root
            # If not then put root really close to relevant internal node
            else:
                if midway < max_out_dist:   # Root should be close to outbase
                    outbasedist = really_close*self.nodedist(outbase, inbase)
                else:                       # Root should be close to inbase
                    outbasedist = (1-really_close)*self.nodedist(outbase, inbase)

            # Root tree at the computed position on the ingroup:outgroup branch
            self.reroot(outbase, inbase, polytomy, outbasedist)

    ###############################################################################################
    def possible_spr_prune_nodes(self):
        """Utililty function when using spr function: where is it possible to prune"""

        possible_prune_nodes = self.nodes - {self.root}
        rootkids = list(self.children(self.root))
        if len(rootkids) == 2:
            if rootkids[0] in self.leaves and rootkids[1] in self.intnodes:
                possible_prune_nodes = possible_prune_nodes - {rootkids[1]}
            elif rootkids[1] in self.leaves and rootkids[0] in self.intnodes:
                possible_prune_nodes = possible_prune_nodes - {rootkids[0]}
        return possible_prune_nodes

    ###############################################################################################

    def possible_spr_regraft_nodes(self, prune_node):
        """Utility function when using spr function: where is it possible to regraft
        prune_node: the node below which pruning will take place (before regrafting)"""

        # Python note: could return preprocessed subtree and treecopy to save computation
        # but this is pretty fast and interface would be less clear.
        # But in case of bottleneclk: also return subtree and treecopy perhaps

        treecopy = self.copy_treeobject()
        subtree = treecopy.subtree(prune_node)
        for leaf in subtree.leaves:
            treecopy.remove_leaf(leaf)
        possible_regraft_nodes = treecopy.leaves - {treecopy.root}
        return possible_regraft_nodes

    ###############################################################################################

    def spr(self, prune_node=None, regraft_node=None):
        """Subtree Pruning and Regrafting.

        prune_node: basenode of subtree that will be pruned.
        regraft_node: node in remaining treestump below which subtree will be grafted

        If no parameters are specified (both are None): perform random SPR
        If only prune_node is specified: choose random regraft_node

        Must specify either both parameters, no parameters, or only prune_node
        """

        # Sanity check: can't perfom SPR on tree with only two leaves
        if len(self.leaves) == 2:
            raise TreeError("Can not perform SPR on tree with only 2 leaves")

        # Invalid argument check
        if prune_node is None and regraft_node is not None:
            msg = ("You only specified regraft_node. "
                   "Must specify either both parameters, no parameters, or only prune_node")
            raise TreeError(msg)

        # Select random prune_node or check the one provided
        possible_prune_nodes = self.possible_spr_prune_nodes()
        if prune_node == None:
            prune_node = random.choice(list(possible_prune_nodes))
        else:
            if prune_node not in possible_prune_nodes:
                raise TreeError(f"Can not prune below {prune_node}")

        # Choose random regraft_node or check the one provided
        possible_regraft_nodes = self.possible_spr_regraft_nodes(prune_node)
        if regraft_node is None:
            regraft_node = random.choice(list(possible_regraft_nodes))
        elif regraft_node not in possible_regraft_nodes:
            msg = f"Specified regraft_node {regraft_node} is not compatible with prune_node"
            raise TreeError(msg)

        # Pruning: Remove subtree
        isleaf = prune_node in self.leaves      # Has to be set before pruning!
        subtree = self.subtree(prune_node)
        for leaf in self.remote_children(prune_node):
            self.remove_leaf(leaf)

        # Regraft: Add subtree back onto remaining tree
        # Special treatment when pruning single leaf (to avoid superfluous internal node)
        self.graft(subtree, regraft_node, graft_with_other_root=isleaf)

###################################################################################################
###################################################################################################
###################################################################################################

class TreeSet():
    """Class for storing and manipulating a number of trees, which all have the same leafs"""

    def __init__(self):
        self.treelist = []

    ###############################################################################################

    def __getitem__(self, index):
        """Implements indexing of treeset.

        Simple index returns single tree.
        Slice returns TreeSet object with selected subset of trees"""
        if isinstance(index, slice):
            newtreeset = TreeSet()
            newtreeset.treelist = self.treelist[index]
            return newtreeset
        else:
            return self.treelist[index]

    ###############################################################################################

    def __len__(self):
        return len(self.treelist)

    ###############################################################################################

    def __iter__(self):
        """Returns fresh iterator object allowing iteration over Treeset (which is itself an iterable)"""
        return self.TreeSetIterator(self)

    ###############################################################################################

    class TreeSetIterator():

        def __init__(self, treeset):
            self.treelist = treeset.treelist
            self.i = 0
            self.max = len(treeset.treelist)

        def __iter__(self):
            return self

        def __next__(self):
            if self.i < self.max:
                next_tree = self.treelist[self.i]
                self.i += 1
                return next_tree
            else:
                raise StopIteration

    ###############################################################################################

    def addtree(self, tree):
        """Adds Tree object to Treeset object"""
        self.treelist.append(tree)

    ###############################################################################################

    def addtreeset(self, treeset):
        """Adds all trees in TreeSet object to this TreeSet object"""
        for tree in treeset:
            self.addtree(tree)

    ###############################################################################################

    def rootmid(self):
        """Performs midpoint rooting on all trees in TreeSet"""
        for tree in self.treelist:
            tree.rootmid()

    ###############################################################################################

    def nexus(self, printdist=True, printlabels=True, translateblock=True):
        """Returns nexus format tree as a string"""

        # Construct header
        stringlist = ["#NEXUS\n\nbegin trees;\n"]

        # If translateblock is requested: add translateblock to stringlist
        if translateblock:
            transdict = self[0].transdict()
            stringlist.append(self[0].translateblock(transdict))

        # Add newick tree strings
        for i, tree in enumerate(self.treelist):
            stringlist.append("    tree t.{} = ".format(i+1))
            if translateblock:
                stringlist.append(tree.newick(printdist=printdist, printlabels=printlabels, transdict=transdict))
            else:
                stringlist.append(tree.newick(printdist=printdist, printlabels=printlabels))
            stringlist.append("\n")

        # Add footer
        stringlist.append("end;\n")

        return "".join(stringlist)

    ###############################################################################################

    def newick(self, printdist=True, printlabels=True):
        """Returns newick format tree as a string"""

        stringlist = []
        for tree in self.treelist:
            stringlist.append(tree.newick(printdist, printlabels))
            stringlist.append("\n")
        return "".join(stringlist)


###################################################################################################
###################################################################################################
###################################################################################################

class RootBipStruct:
    def __init__(self, leafset1, blen1, leafset2, blen2):
        combined_length = blen1 + blen2
        self.count = 1
        self.leafset1 = leafset1
        self.fraction1 = blen1 / combined_length
        self.leafset2 = leafset2
        self.fraction2 = blen2 / combined_length

    def add(self, leafset1, blen1, leafset2, blen2):
        """Adds branch length fractions to the current sum and increments the count."""
        # Make sure blen1 and blen2 refer to the correct leafsets
        if leafset1 != self.leafset1:
            blen2,blen1 = blen1,blen2
        combined_length = blen1 + blen2
        self.count += 1
        self.fraction1 += blen1 / combined_length
        self.fraction2 += blen2 / combined_length

    def merge(self, other):
        """Merges this RootBipStruct with another (for same bipartition)"""
        fraction1, fraction2 = other.fraction1, other.fraction2
        if other.leafset1 != self.leafset1:
            fraction1, fraction2 = fraction2, fraction1
        self.count += other.count
        self.fraction1 += fraction1
        self.fraction2 += fraction2

    def avg_frac(self, leafset):
        if leafset == self.leafset1:
            return self.fraction1 / self.count
        else:
            return self.fraction2 / self.count

###################################################################################################
###################################################################################################

class TreeSummary():
    """Class summarizing bipartitions and branch lengths (but not topologies) from many trees"""

    def __init__(self, trackbips=True, trackclades=False, trackroot=False):
        """TreeSummary constructor. Initializes relevant data structures"""
        self.transdict = None
        self.translateblock = None
        self.tree_count = 0
        self.tree_weight_sum = 0.0
        self._bipartsummary = {}         # Dict: {bipartition:branchstruct with extra fields}
        self._bipartsummary_processed = False
        self._cladesummary = {}         # Dict: {clade:nodestruct with extra fields}
        self._cladesummary_processed = False
        self._sorted_biplist = None
        self.trackroot = trackroot
        self.trackbips = trackbips
        self.trackclades = trackclades
        self._sorted_rootbips = None
        if trackroot:
            self._rootbip_summary = {}
        else:
            self._rootbip_summary = None

    ###############################################################################################

    def __len__(self):
        return self.tree_count

    ###############################################################################################

    @property
    def bipartsummary(self):
        """Property method for lazy evaluation of freq, var, and sem for branches"""
        if not self._bipartsummary_processed:
            for branch in self._bipartsummary.values():
                branch.freq = branch.SUMW / self.tree_weight_sum
                n = branch.bip_count
                if n > 1:
                    branch.var = branch.T * n / ((n - 1) * branch.SUMW)
                    branch.sem = math.sqrt(branch.var)/math.sqrt(n)
                else:
                    branch.var = "NA"
                    branch.sem = "NA"
            self._bipartsummary_processed = True

        return self._bipartsummary

    ###############################################################################################

    @property
    def cladesummary(self):
        """Property method for lazy evaluation of freq, var, and sem for node depths"""
        if not self._cladesummary_processed:
            for node in self._cladesummary.values():
                node.freq = node.SUMW / self.tree_weight_sum
                n = node.clade_count
                if n > 1:
                    node.var = node.T * n / ((n - 1) * node.SUMW)
                    node.sem = math.sqrt(node.var)/math.sqrt(n)
                else:
                    node.var = "NA"
                    node.sem = "NA"
            self._cladesummary_processed = True

        return self._cladesummary

    ###############################################################################################

    @property
    def sorted_biplist(self):
        """Return list of bipartitions.
        First external (leaf) bipartitions sorted by leafname.
        Then internal bipartitions sorted by freq"""

        if self._sorted_biplist == None:
            leafbips = []
            internalbips = []

            for bip, branch in self.bipartsummary.items():
                (bip1,bip2) = bip
                if len(bip1) == 1:
                    leafname = next(iter(bip1))
                    leafbips.append((leafname, bip))
                elif len(bip2) == 1:
                    leafname = next(iter(bip2))
                    leafbips.append((leafname, bip))
                else:
                    internalbips.append((branch.freq,bip))

            leafbips = sorted(leafbips, key=itemgetter(0))
            internalbips = sorted(internalbips, key=itemgetter(0), reverse=True)
            self._sorted_biplist = leafbips + internalbips

        return self._sorted_biplist

    ###############################################################################################

    @property
    def sorted_rootbips(self):
        """Return list of root-bipartitions (branches where root has been seen), sorted by
        occurrence (count on tree samples added to TreeSummary)"""

        if self._sorted_rootbips == None:
            self._sorted_rootbips = []
            maxcount = 0
            for bip,rootbipstruct in self._rootbip_summary.items():
                self._sorted_rootbips.append((rootbipstruct.count, bip, rootbipstruct))
            self._sorted_rootbips.sort(key=itemgetter(0), reverse=True)
        return self._sorted_rootbips

    ###############################################################################################

    def add_branchid(self):
        """Adds attribute .branchID to all bipartitions in .bipartsummary
        External bipartitions are labeled with the leafname.
        Internal bipartitions are labeled with consecutive numbers by decreasing frequency"""

        # Python note: Start enumeration at 1
        for freqrank, (sortkey,bipart) in enumerate(self.sorted_biplist, 1):
            if type(sortkey) == str:
                self.bipartsummary[bipart].branchID = sortkey
            else:
                self.bipartsummary[bipart].branchID = freqrank

    ###############################################################################################

    def add_tree(self, curtree, weight=1.0):
        """Add tree object to treesummary, update all relevant bipartition summaries"""

        # Main interface to TreeSummary.
        # Takes tree object, updates relevant measures
        # First time entered: build set of leaves for consistency checking.
        # Also compute transdict and translateblock for tree reporting (and storage?)
        if self.tree_count == 0:
            self.leaves = curtree.leaves
            self.transdict = curtree.transdict()
            self.translateblock = curtree.translateblock(self.transdict)
        elif curtree.leaves != self.leaves:
            msg = f"Leaves on tree number {self.tree_count} are different than previous trees"
            raise TreeError(msg)

        self.tree_count += 1
        self.tree_weight_sum += weight       # The weighted equivalent of tree_count

        if self.trackroot:
            self._add_root(curtree)

        cladedict = None
        if self.trackclades:
            cladedict = self._add_clade(curtree, weight)

        bipdict = None
        if self.trackbips:
            bipdict = self._add_bip(curtree, weight)

        return bipdict, cladedict

    ###############################################################################################

    def _add_root(self, curtree):
        """Helper method for add_tree: handles roots"""

        bipartition, leafset1, blen1, leafset2, blen2 = curtree.rootbip()
        if bipartition in self._rootbip_summary:
            self._rootbip_summary[bipartition].add(leafset1, blen1, leafset2, blen2)
        else:
            self._rootbip_summary[bipartition] = RootBipStruct(leafset1, blen1, leafset2, blen2)

    ###############################################################################################

    def _add_clade(self, curtree, weight):
        """Helper method to add_tree: handles clades"""

        self._cladesummary_processed = False

        cladedict = curtree.cladedict()
        for clade,nodestruct in cladedict.items():
            depth = nodestruct.depth

            # If clade has been seen before: update existing info
            if clade in self._cladesummary:
                Q = depth - self._cladesummary[clade].depth
                TEMP = self._cladesummary[clade].SUMW + weight
                R = Q*weight/TEMP
                self._cladesummary[clade].depth += R
                self._cladesummary[clade].T += R * self._cladesummary[clade].SUMW * Q
                self._cladesummary[clade].SUMW = TEMP
                self._cladesummary[clade].clade_count += 1

            # If bipartition has never been seen before: add to dict and add online attributes
            else:
                self._cladesummary[clade]=nodestruct
                self._cladesummary[clade].clade_count = 1
                self._cladesummary[clade].SUMW = weight
                self._cladesummary[clade].depth = depth
                self._cladesummary[clade].T = 0.0

        return cladedict

    ###############################################################################################

    def _add_bip(self, curtree, weight):
        """Helper method to add_tree: handles bipartitions"""

        self._bipartsummary_processed = False
        self._sorted_biplist = None

        # I am interested in being able to compute weighted frequency of a bipartition as well as
        # the weighted mean and weighted variance of the branch length for that bipartition.
        # In order to do this I follow the robust (= no underflow/overflow problems), one-pass
        # approach described in D.H.D. West, "Updating Mean and Variance Estimates: An Improved
        # Method", Communications of the ACM, 22(9), 1979.
        # A number of variables are used, some of which are stored in the summary dictionary.
        # These variables have mostly been named according to the original paper.
        # Exceptions are "bip_count" which was "N", "weight" which was "W", "mean" which was "M",
        # and "brlen" which was "X".
        # Note: mean branch length is stored in two attributes: mean and length
        # This is due to other functions that expect the attribute .length to be present
        bipdict = curtree.bipdict()
        for bipart,branchstruct in bipdict.items():
            brlen = branchstruct.length

            # If bipartition has been seen before: update existing info
            if bipart in self._bipartsummary:
                Q = brlen - self._bipartsummary[bipart].length
                TEMP = self._bipartsummary[bipart].SUMW + weight
                R = Q*weight/TEMP
                self._bipartsummary[bipart].length += R
                self._bipartsummary[bipart].T += R * self._bipartsummary[bipart].SUMW * Q
                self._bipartsummary[bipart].SUMW = TEMP
                self._bipartsummary[bipart].bip_count += 1

            # If bipartition has never been seen before: add to dict and add online attributes
            else:
                self._bipartsummary[bipart]=branchstruct
                self._bipartsummary[bipart].bip_count = 1
                self._bipartsummary[bipart].SUMW = weight
                self._bipartsummary[bipart].length = brlen
                self._bipartsummary[bipart].T = 0.0

        return bipdict

    ###############################################################################################

    def update(self, other):
        """Merge this object with external treesummary"""

        # Sanity check: do two treesummaries refer to same set of leaves?
        if self.leaves != other.leaves:
            msg = "Not all trees have same set of leaves."
            raise TreeError(msg)

        # Update treecount and weight
        self.tree_count += other.tree_count
        self.tree_weight_sum += other.tree_weight_sum

        if self.trackbips:
            self._updatebip(other)

        if self.trackclades:
            self._updateclade(other)

        if self.trackroot:
            self._updateroot(other)


    ###############################################################################################

    def _updatebip(self, other):

        # Merge "treesummary.bipartsummary" with "self.bipartsummary"
        other_bipsum = other.bipartsummary
        self_bipsum = self.bipartsummary

        for bipart in other_bipsum:
            # If bipart already in self.bipartsummary, update fields
            if bipart in self_bipsum:

                sumw1 = self_bipsum[bipart].SUMW
                sumw2 = other_bipsum[bipart].SUMW
                mean1 = self_bipsum[bipart].length
                mean2 = other_bipsum[bipart].length
                t1 = self_bipsum[bipart].T
                t2 = other_bipsum[bipart].T

                self_bipsum[bipart].bip_count += other_bipsum[bipart].bip_count
                self_bipsum[bipart].length = (mean1*sumw1 + mean2*sumw2)/(sumw1+sumw2)
                self_bipsum[bipart].SUMW += other_bipsum[bipart].SUMW

                # Note: the following expression was arrived at empirically!
                # I have not proven this is correct, but it does seem to work...
                self_bipsum[bipart].T = t1+t2+sumw1*sumw2*(mean2-mean1)*(mean2-mean1)/(sumw1+sumw2)

            # If bipartition has never been seen before: transfer Branchstruct from other_bipsum:
            else:
                self_bipsum[bipart] = other_bipsum[bipart]

        self._bipartsummary_processed = False
        self._sorted_biplist = None

    ###############################################################################################

    def _updateclade(self, other):

        # Merge "treesummary.bipartsummary" with "self.bipartsummary"
        other_cladesum = other.cladesummary
        self_cladesum = self.cladesummary

        for clade in other_cladesum:
            # If bipart already in self.bipartsummary, update fields
            if clade in self_cladesum:

                sumw1 = self_cladesum[clade].SUMW
                sumw2 = other_cladesum[clade].SUMW
                mean1 = self_cladesum[clade].depth
                mean2 = other_cladesum[clade].depth
                t1 = self_cladesum[clade].T
                t2 = other_cladesum[clade].T

                self_cladesum[clade].clade_count += other_cladesum[clade].clade_count
                self_cladesum[clade].depth = (mean1*sumw1 + mean2*sumw2)/(sumw1+sumw2)
                self_cladesum[clade].SUMW += other_cladesum[clade].SUMW

                # Note: the following expression was arrived at empirically!
                # I have not proven this is correct, but it does seem to work...
                self_cladesum[clade].T = t1+t2+sumw1*sumw2*(mean2-mean1)*(mean2-mean1)/(sumw1+sumw2)

            # If bipartition has never been seen before: transfer Branchstruct from other_bipsum:
            else:
                self_cladesum[clade] = other_cladesum[clade]

        self._cladesummary_processed = False
        self._sorted_biplist = None


    ###############################################################################################

    def _updateroot(self, other):

        # Python note: am i missing something here???

        other_rootbipsum = other._rootbip_summary
        self_rootbipsum = self._rootbip_summary
        for rootbip in other_rootbipsum:
            if rootbip in self_rootbipsum:
                self_rootbipsum[rootbip].merge(other_rootbipsum[rootbip])
            else:
                self_rootbipsum[rootbip] = other_rootbipsum[rootbip]

    ###############################################################################################

    def log_bipart_credibility(self, biptopology):
        """Compute log bipartition credibility for topology (sum of log(freq) for all branches)"""

        bipartsummary = self.bipartsummary
        logsum = 0.0
        for bipartition in biptopology:
            logsum += math.log(bipartsummary[bipartition].freq)
        return logsum

    ###############################################################################################

    def log_clade_credibility(self, cladetopology):
        """Compute log clade credibility for topology (sum of log(freq) for all clades)"""

        cladesummary = self.cladesummary
        logsum = 0.0
        for clade in cladetopology:
            logsum += math.log(cladesummary[clade].freq)
        return logsum

    ###############################################################################################

    def contree(self, cutoff=0.5, allcompat=False, labeldigits=3):
        """Returns a consensus tree built from selected bipartitions."""

        if cutoff < 0.5:
            msg = "Consensus tree cutoff has to be at least 0.5"
            raise TreeError(msg)

        # Transfer biparts and branches with freq>cutoff to new bipdict, create tree
        conbipdict = {}
        i = 0
        for _, bip in self.sorted_biplist:
            i += 1
            branch = self.bipartsummary[bip]
            if branch.freq < cutoff:
                break
            branch.label = f"{branch.freq:.{labeldigits}f}"
            conbipdict[bip] = branch
        contree = Tree.from_biplist(conbipdict)

        # If allcompat has been requested: add remaining, compatible bipartitions to contree
        if allcompat:
            for j in range(i, len(self.sorted_biplist)):
                if contree.is_resolved():
                    break
                _,bip = self.sorted_biplist[j]
                branch = self.bipartsummary[bip]
                branch.label= f"{branch.freq:.{labeldigits}f}"
                is_present, is_compatible, insert_tuple = contree.check_bip_compatibility(bip)
                if is_compatible and (not is_present):
                    parentnode, childnodes = insert_tuple
                    contree.insert_node(parentnode, childnodes, branch)
                    contree._remotechildren_dict = None

        return contree

    ###############################################################################################

    def root_maxfreq(self, summary_tree):
        """Uses info about root bipartitions in TreeSummary to place root on summary tree.
        Divides length of root bipartition among two branches in accordance with average
        fraction of lengths seen for this rootbip across all trees."""

        # Starting with most frequent root location: find one that is compatible
        # Python note: should i just try number 1 on sorted list?
        if summary_tree.is_bifurcation(summary_tree.root):
            cur_rootbip, _, _, _, _ = summary_tree.rootbip()
        else:
            cur_rootbip = None
        for count, bip, summary_rootbipstruct in self.sorted_rootbips:
            if summary_tree.bipart_is_present(bip):
                # If tree already rooted correctly: do not reroot
                if (cur_rootbip is None) or (bip != cur_rootbip):
                    parent,child = summary_tree.find_bipart_nodes(bip)
                    summary_tree.deroot()  # Python note: necessary?
                                           # reroot seems to assume not rooted at birfurcation
                                           # rethink reroot function and others depending on it!
                    summary_tree.reroot(child, parent)
                summary_tree.rootcred = count / self.tree_count

                # Divide branch lengths for two rootkids according to fractions
                # seen for this rootbip across trees in ._rootbip_summary
                kid1,kid2 = summary_tree.children(summary_tree.root)
                biplen = summary_tree.nodedist(kid1, kid2)
                kid1_remkids = summary_tree.remotechildren_dict[kid1]
                dist_to_kid1 = biplen * summary_rootbipstruct.avg_frac(kid1_remkids)
                dist_to_kid2 = biplen - dist_to_kid1
                summary_tree.child_dict[summary_tree.root][kid1].length = dist_to_kid1
                summary_tree.child_dict[summary_tree.root][kid2].length = dist_to_kid2

                return summary_tree

        # If we did not return by now, then bipart not in contree
        raise TreeError(f"Summary_tree tree not compatible with any observed root locations")

    ###############################################################################################

    def set_mean_node_depths(self, summary_tree):
        """Set branch lengths on summary tree based on mean node depth for clades corresponding
        to parent and child nodes (blen = depth_parent - depth_child).

        NOTE 1: only meaningful if input trees are based on a clock model.
        NOTE 2: only works if all clades in tree have been observed at least once. The option
                will therefore not work with all rootings, and may also fail for majority rule
                consensus trees
        NOTE 3: only uses node depths from monopyletic clades (so some values may be set
        based on very few trees)"""

        all_leaves = summary_tree.frozenset_leaves
        sorted_leafs = summary_tree.sorted_leaf_list
        leaf2index = summary_tree.leaf2index

        try:
            for parent in summary_tree.sorted_intnodes(deepfirst=True):
                p_remkids = summary_tree.remotechildren_dict[parent]
                p_clade = Clade(p_remkids, all_leaves, sorted_leafs, leaf2index)
                p_depth = self.cladesummary[p_clade].depth
                for child in summary_tree.children(parent):
                    c_remkids = summary_tree.remotechildren_dict[child]
                    c_clade = Clade(c_remkids, all_leaves, sorted_leafs, leaf2index)
                    c_depth = self.cladesummary[c_clade].depth
                    blen = p_depth - c_depth
                    summary_tree.setlength(parent, child, blen)
        except KeyError as e:
            raise TreeError("Problem while setting mean node depths: the following clade has not been "
                            + "observed among input trees: check rooting of tree."
                            + f"\n{e.args[0]}")

        return summary_tree

    ###############################################################################################

    def compute_rootcred(self, tree):
        """Returns root credibility (frequency of tree's root among observed trees) based
        on current root of summary_tree and information in self._rootbip_summary"""

        # Python note: mostly relevant for MCC trees. Other tree types will typically
        # get the .rootcred attribute set when calling .root_max_freq()

        tree_rootbip, _, _, _, _ = tree.rootbip()
        summary_rootbipstruct = self._rootbip_summary[tree_rootbip]
        return summary_rootbipstruct.count / self.tree_count

###################################################################################################
###################################################################################################

class BigTreeSummary(TreeSummary):
    """Class summarizing bipartitions, branch lengths, and topologies from many trees"""

    # Does everything TreeSummary does and also keeps track of topologies
    # (topology list is potentially quite big, which is the reason for not including it in TS)

    def __init__(self, store_trees=False, trackbips=True, trackclades=False, trackroot=False):
        TreeSummary.__init__(self, trackbips, trackclades, trackroot)
        self._biptoposummary = {}
        self._biptoposummary_processed = False
        self._cladetoposummary = {}
        self._cladetoposummary_processed = False
        self.store_trees = store_trees

    ###############################################################################################

    @property
    def biptoposummary(self):
        """Property method for lazy evaluation of topostruct.freq"""
        if not self._biptoposummary_processed:
            for topostruct in self._biptoposummary.values():
                topostruct.freq = topostruct.weight / self.tree_weight_sum
            self._biptoposummary_processed = True

        return self._biptoposummary

    ###############################################################################################

    @property
    def cladetoposummary(self):
        """Property method for lazy evaluation of topostruct.freq"""
        if not self._cladetoposummary_processed:
            for topostruct in self._cladetoposummary.values():
                topostruct.freq = topostruct.weight / self.tree_weight_sum
            self._cladetoposummary_processed = True

        return self._cladetoposummary

    ###############################################################################################

    def add_tree(self, curtree, weight=1.0):
        """Add tree to treesummary, update all summaries"""

        self._biptoposummary_processed = False
        self._cladetoposummary_processed = False

        # Superclass method takes care of updating n_trees and all bipart/clade-related info
        # Python note: bipdict or cladedict are None when bip or clade not tracked
        bipdict, cladedict = TreeSummary.add_tree(self, curtree, weight)

        if self.trackbips:
            self._addbiptopo(bipdict, curtree, weight)

        if self.trackclades:
            self._addcladetopo(cladedict, curtree, weight)

    ###############################################################################################

    def _addbiptopo(self, bipdict, curtree, weight):

        # If biptopology has never been seen before, then add it and initialize count
        # If topology HAS been seen before then update count
        topology = frozenset(bipdict.keys())
        if topology in self._biptoposummary:
            self._biptoposummary[topology].weight += weight
        else:
            self._biptoposummary[topology]=Topostruct()
            self._biptoposummary[topology].weight = weight
            if self.store_trees:
                curtree.clear_caches()
                self._biptoposummary[topology].tree = curtree

    ###############################################################################################

    def _addcladetopo(self, cladedict, curtree, weight):

        # If biptopology has never been seen before, then add it and initialize count
        # If topology HAS been seen before then update count
        topology = frozenset(cladedict.keys())
        if topology in self._biptoposummary:
            self._cladetoposummary[topology].weight += weight
        else:
            self._cladetoposummary[topology]=Topostruct()
            self._cladetoposummary[topology].weight = weight
            if self.store_trees:
                curtree.clear_caches()
                self._cladetoposummary[topology].tree = curtree

    ###############################################################################################

    def update(self, other):
        """Merge this object with other treesummary"""

        self._biptoposummary_processed = False
        self._cladetoposummary_processed = False

        # Superclass method takes care of updating:
        # tree_count, tree_weight_sum, and bipart/clade-summary
        TreeSummary.update(self, other)

        for biptopology in other._biptoposummary:
            # If topology already in self.toposummary, update count
            if biptopology in self._biptoposummary:
                self._biptoposummary[biptopology].weight += other._biptoposummary[biptopology].weight
            # If topology has never been seen before, simply transfer entry
            else:
                self._biptoposummary[biptopology]=other._biptoposummary[biptopology]

        for cladetopology in other._cladetoposummary:
            # If topology already in self.toposummary, update count
            if cladetopology in self._cladetoposummary:
                self._cladetoposummary[cladetopology].weight += other._cladetoposummary[cladetopology].weight
            # If topology has never been seen before, simply transfer entry
            else:
                self._cladetoposummary[cladetopology]=other._cladetoposummary[cladetopology]

    ###############################################################################################

    def max_bipart_cred_tree(self, labeldigits=3):
        """Find maximum bipartition credibility tree. Return tuple of (maxcredtree, maxlogcred)"""

        maxlogcred = -math.inf
        for biptopology in self.biptoposummary:
            logcred = self.log_bipart_credibility(biptopology)
            if logcred > maxlogcred:
                maxlogcred = logcred
                maxlogcredbiptopo = biptopology

        maxcredbipdict = {}
        for bipartition in maxlogcredbiptopo:
            branch = self.bipartsummary[bipartition]
            branch.label = f"{round(branch.freq, labeldigits)}"
            maxcredbipdict[bipartition] = branch

        # Build tree from bipartitions  in new bipdict
        maxcredtree = Tree.from_biplist(maxcredbipdict)

        return maxcredtree, maxlogcred

    ###############################################################################################

    def max_clade_cred_tree(self, labeldigits=3):
        """Find maximum clade credibility tree. Return tuple of (maxcredtree, maxlogcred)"""

        maxlogcred = -math.inf
        for clade_topology in self.cladetoposummary:
            logcred = self.log_clade_credibility(clade_topology)
            if logcred > maxlogcred:
                maxlogcred = logcred
                maxlogcred_cladetopo = clade_topology

        maxcred_bipdict = {}
        for clade in maxlogcred_cladetopo:
            nodestruct = self.cladesummary[clade]
            nodestruct.label = f"{round(nodestruct.freq, labeldigits)}"
            maxcred_bipdict[clade] = nodestruct
        maxcredtree = Tree.from_cladedict(maxcred_bipdict)

        return maxcredtree, maxlogcred

###################################################################################################
###################################################################################################
###################################################################################################

class TreefileBase():
    """Abstract base-class for representing tree file objects."""

    # Classes for specific formats inherit from this class and add extra stuff as needed.
    # NOTE: i am opening files in "read text" with encoding UTF-8. Will this work across platforms?

    def __init__(self, filename, filecontent, interner):

        num_args = (filename is not None) + (filecontent is not None)
        if num_args != 1:
            raise TreeError("TreefileBase __init__ requires either filename or filecontent (not both)")
        elif filecontent:
            self.treefile = StringIO(filecontent)
        else:
            self.treefile = open(filename, mode="rt", encoding="UTF-8")

        self.interner = interner
        self.buffer = ""                # Used for keeping leftovers after reading whole line
        self.below_root = None

    ###############################################################################################

    def __enter__(self):
        """Implements context manager behaviour for TreefileBase types.
        Usage example:
            with pt.Newicktreefile(filename) as tf:
                mytree = tf.readtree()
            mytree.rootminvar()
        """
        return self

    ###############################################################################################

    def __exit__(self, type, value, traceback):
        """Implements context manager behaviour for TreefileBase types.
        Usage example:
            with pt.Newicktreefile(filename) as tf:
                mytree = tf.readtree()
            mytree.rootminvar()
        """
        # Python note: consider adding code to handle exceptions (but maybe ok to just pass on?)
        self.close()

    ###############################################################################################

    def get_treestring(self):
        """Return next tree-string"""

        # We are now at beginning of treestring, read until semi-colon encountered
        # Raise StopIteration when EOF has been reached
        stringlist = [self.buffer]

        if ";" not in self.buffer:           # Only read on if end of treestring not reached
            for line in self.treefile:
                stringlist.append(line)
                if ";" in line:
                    break

        treestring = "".join(stringlist)

        # If we got this far and still haven't found ";" then EOF must have been reached
        # Signal EOF to caller, which can then clean up and stop iteration
        if ";" not in treestring:
            return None

        stringparts = treestring.split(";")
        treestring = "".join([stringparts[0], ";"])
        self.buffer = "".join(stringparts[1:])

        # Finally: check if there is information below the root node,
        # i.e., between last right parenthesis and the semicolon.
        # If so then save as self.below_root
        # Python note: can possibly be done smarter. And what should happen with below root part?
        last_right_paren = treestring.rfind(')')
        if last_right_paren != -1:
            self.below_root = treestring[last_right_paren + 1:-1]
            if self.below_root:
                treestring = treestring[:last_right_paren + 1] + ";"
        return treestring

    ###############################################################################################

    def readtree(self, returntree=True):
        """Reads one treestring from file and returns as Tree object if requested.
        Returns None when exhausted file"""

        try:
            tree = self.__next__(returntree)
            return tree
        except StopIteration:
            self.treefile.close()
            return None

    ###############################################################################################

    def readtrees(self, discardprop=0.0):
        """Reads trees from file and returns as TreeSet object. Can discard fraction of trees"""

        # Avoid erroneous error message when running pylint ("self" is OK for iteration)
        # pylint: disable=not-an-iterable

        treeset = TreeSet()
        for tree in self:
            treeset.addtree(tree)
        if discardprop > 0:
            ndiscard = int(discardprop * len(treeset))
            treeset = treeset[ndiscard:]

        return treeset

    ###############################################################################################

    def close(self):
        """For explicit closing of Treefile before content exhausted"""
        self.treefile.close()

###################################################################################################
###################################################################################################


class Newicktreefile(TreefileBase):
    """Class representing Newick tree file. Iteration returns tree-objects"""

    def __init__(self, filename=None, filecontent=None, interner=None):
        TreefileBase.__init__(self, filename, filecontent, interner)
        # HACK!!! Minimal file format check:
        # Read first three lines in file, check whether any of them contains "#NEXUS".
        # If so then this is presumably NOT a Newick file (but a nexus file...) => exit.
        # If not, then hope that this IS a Newick file. Reset file pointer and proceed.
        filestart = self.treefile.readline()
        filestart += self.treefile.readline()
        filestart += self.treefile.readline()
        if filestart.find("#NEXUS") != -1:
            msg = "File does not appear to be in Newick format"
            raise TreeError(msg)
        else:
            self.treefile.seek(0)       # Reset pointer to start of file

    ###############################################################################################

    def __iter__(self):
        return self

    ###############################################################################################

    def __next__(self, returntree=True):
        treestring = self.get_treestring()
        if treestring is None:
            self.treefile.close()
            raise StopIteration
        else:
            treestring = remove_comments(treestring)
            if returntree:
                tree = Tree.from_string(treestring, self.interner)
                tree.below_root = self.below_root
                return tree

###################################################################################################
###################################################################################################
###################################################################################################


class Nexustreefile(TreefileBase):
    """Class representing Nexus tree file. Iteration returns tree object or None"""

    ###############################################################################################

    def __init__(self, filename=None, filecontent=None, interner=None):
        """Read past NEXUS file header, parse translate block if present"""

        TreefileBase.__init__(self, filename, filecontent, interner)

        ###########################################################################################

        def skip_comment(line):
            """Reads past a NEXUS comment"""

            if line.count("[") > line.count("]"):
                for extraline in self.treefile:
                    line += extraline
                    if line.count("[") == line.count("]"):
                        break

            # Now we should have an equal number of "[" and "]"
            # Remove comments
            return remove_comments(line)

        ###########################################################################################

        def read_until(pattern, skipcomment=False):
            """Read and return filecontent up to and including 'pattern'"""

            readsofar = ""

            # Read one line until text pattern encountered. Skip comments, chek for EOF
            for line in self.treefile:
                if "[" in line and skipcomment:
                    line = skip_comment(line)
                elif line == "":
                    msg = "Bug: EOF reached before read_until pattern was encountered"
                    raise TreeError(msg)

                readsofar += line
                if pattern.search(readsofar):
                    break

            return readsofar

        ###########################################################################################

        # EXECUTION of __init__ STARTS HERE
        self.transdict = None

        # Regular expressions that are used frequently during execution
        # DOTALL flag forces "." to also match newlines
        # "*?" means non-greedy matching - match as little text as possible
        self.begin_trees_pattern = re.compile(r"""
                                    begin\s+trees           # "begin trees"
                                    (\s+[\w\-\/\.]+\s*)?    # possible block name
                                    ;                       # semicolon after "begin trees"
                                    """, re.IGNORECASE | re.VERBOSE)

        self.tree_header_pattern = re.compile(r"""
                                    ^.*?u?tree              # Anything up to the first "(u)tree"
                                    \s+(\*\s)?\s*           # tree name may be preceeded by "* "
                                    [\w\-\/\.]+             # Tree name
                                    \s*=                    # Whitespace and "="
                                    """, re.IGNORECASE | re.DOTALL | re.VERBOSE)

        self.end_pattern = re.compile(r"\send(block)?;", re.IGNORECASE)

        # Read past "begin trees;" (take all allowed variants into account)
        filestart = read_until(self.begin_trees_pattern, skipcomment=True)

        # Minimal fileformat check: start of NEXUS file should contain "#NEXUS". If not then exit:
        if filestart.find("#NEXUS") == -1:      # The string "#NEXUS" was not found
            msg = "File does not appear to be in NEXUS format"
            raise TreeError(msg)

        # Read past first "tree <NAME> =" statement. Keep text that was read in buffer
        self.buffer = read_until(self.tree_header_pattern, skipcomment=True)

        # If buffer contains a "translate" block: parse it
        if re.search("translate", self.buffer, re.IGNORECASE):

            self.transdict = {}

            # Remove start of buffer
            pattern = re.compile(r"^.*translate\s*", re.IGNORECASE | re.DOTALL)
            self.buffer = pattern.sub("", self.buffer)

            # Remove end of buffer, only "translate" block core remains after this
            # NOTE: I originally had following complicated pattern
            # This was found to break on some trees.
            # Why did I use that? Are there some cases that I am no longer covering?
            # re.compile(";\s+u?tree\s+(\*\s)?\s*[\w\-\/\.]+\s*=.*", re.IGNORECASE | re.DOTALL)
            pattern = re.compile(r";\s+.*tree.*", re.IGNORECASE | re.DOTALL)
            transblock = pattern.sub("", self.buffer)

            # Split on commas, generating list of "code-whitespace-origname" pairs
            # then split each of these on whitespace, and build translation dictionary
            translist = transblock.split(",")
            for transpair in translist:
                words = transpair.split()
                code = words[0]
                realname = sys.intern(words[1])
                self.transdict[code] = realname

        # Keep everything after "tree <NAME> =" in self.buffer: may contain tree!
        # Non-greedy mathching *?: find first = and then every non-start parenthesis
        pattern = re.compile("^.*?=[^(]*", re.DOTALL)
        self.buffer = pattern.sub("", self.buffer)

    ###############################################################################################

    def __iter__(self):
        return self

    ###############################################################################################

    def __next__(self, returntree=True):

        treestring = self.get_treestring()
        if treestring is None:
            self.treefile.close()
            raise StopIteration

        # remove comments in brackets if present
        # remove leading "tree NAME =" (compiled regexp "tree_header_pattern")
        # Implementation note: does not deal with figtree comments!
        treestring = remove_comments(treestring)
        treestring = self.tree_header_pattern.sub("", treestring)

        # If "end;" statement has been reached: terminate for loop, do NOT return tree object
        if self.end_pattern.search(treestring):
            self.treefile.close()
            raise StopIteration

        # Return tree object if requested
        if returntree:
            tree = Tree.from_string(treestring, self.transdict, self.interner)
            tree.below_root = self.below_root
            return tree

###################################################################################################
###################################################################################################
###################################################################################################

class Treefile:
    """Factory for making Newick or Nexus treefile objects. Autodetects fileformat"""

    def __new__(klass, filename):

        def read_until_non_comment(filename):
            with open(filename, 'r') as file:
                content = []
                comment_level = 0  # Counter to track the nesting level of comments
                for line in file:
                    stripped_line = ""
                    i = 0
                    while i < len(line):
                        char = line[i]
                        if char == '[':  # Start of a comment or deeper nesting
                            comment_level += 1
                        elif char == ']':  # End of a comment or moving up a nesting level
                            comment_level -= 1
                            if comment_level < 0:  # Malformed comment structure
                                raise ValueError("Mismatched comment brackets in the file.")
                            i += 1  # Skip the closing bracket
                            continue
                        elif comment_level == 0:  # Not inside a comment
                            stripped_line += char
                        i += 1
                    content.append(line)
                    if stripped_line.strip():  # If there's any non-commented text
                        break
                return ''.join(content)

        headertext = read_until_non_comment(filename)
        if "#nexus" in headertext.lower():
            return Nexustreefile(filename)
        else:
            return Newicktreefile(filename)

###################################################################################################
###################################################################################################
###################################################################################################
##
##          Treebuilding section: classes and methods for building trees
##
###################################################################################################
###################################################################################################

class Distmatrix(object):
    """Class representing distance matrix for set of taxa. Knows how to compute trees"""

    def __init__(self):
        self.dmat = None
        self.namelist = None
        self.n = None
        self.name2index = None
        self.index2name = None

    #######################################################################################

    @classmethod
    def from_distdict(cls, distdict):
        """Construct Distmatrix object from nested dictionary of dists: distdict[name1][name2] = dist"""

        # Note: Should add more careful parsing and error checking at some point: check all pairs are accounted for!
        self = cls()
        self.namelist = sorted(list(distdict.keys()))
        self.n = len(self.namelist)
        self.dmat = np.zeros((self.n, self.n))
        self.name2index = dict(zip(self.namelist, range(self.n)))
        self.index2name = dict(zip(range(self.n), self.namelist))

        for name1, name2 in itertools.combinations(self.namelist,2):
            dist = distdict[name1][name2]
            self.setdist(name1, name2, dist)

        return self

    #######################################################################################

    @classmethod
    def from_numpy_array(cls, nparray, namelist):
        """
        Construct Distmatrix object from numpy array and corresponding list of names
        Names in namelist must be in same order as indices in numpy 2D array
        """
        self = cls()
        self.namelist = namelist
        self.n = len(namelist)
        self.dmat = nparray
        self.name2index = dict(zip(self.namelist, range(self.n)))
        self.index2name = dict(zip(range(self.n), self.namelist))

        return self

    #######################################################################################

    @classmethod
    def from_distfile(cls, distfilename):
        """
        Construct Distmatrix object from file containing rows of: name1 name2 distance
        """
        self = cls()
        name1list = []
        name2list = []
        distlist = []
        with open(distfilename, "r") as distfile:
            for line in distfile:
                words = line.split()
                name1list.append(words[0])
                name2list.append(words[1])
                distlist.append(float(words[2]))
        uniqnames = set(name1list) | set(name2list)
        self.namelist = sorted(uniqnames)
        self.n = len(self.namelist)
        self.name2index = dict(zip(self.namelist, range(self.n)))
        self.index2name = dict(zip(range(self.n), self.namelist))

        self.dmat = np.zeros((self.n, self.n))
        for name1,name2,dist in zip(name1list, name2list, distlist):
            self.setdist(name1, name2, dist)

        return self

    #######################################################################################

    def __str__(self):
        """Returns distance matrix as string"""

        # Format of output like this:
        ##        S5    s1  s2  s3  s4
        ##        0.0   0.625   0.625   0.5 0.5
        ##        0.625 0.0 0.25    0.75    0.75
        ##        0.625 0.25    0.0 0.75    0.75
        ##        0.5   0.75    0.75    0.0 0.25
        ##        0.5   0.75    0.75    0.25    0.0
        namelist = sorted(list(self.namelist))
        tmplist = []

        # Header line: all names in order, tab separated
        tmplist.append("\t")
        for name in namelist:
            tmplist.append(str(name))
            tmplist.append("\t")
        tmplist.pop()
        tmplist.append("\n")

        # Body of matrix. Same order as namelist
        for name1 in namelist:
            tmplist.append(name1)
            tmplist.append("\t")
            for name2 in namelist:
                tmplist.append(str(self.getdist(name1,name2)))
                tmplist.append("\t")
            tmplist.pop()
            tmplist.append("\n")

        return "".join(tmplist)

    #######################################################################################

    def clean_names(self, illegal=",:;()[]", rep="_"):
        """Rename items to avoid characters that are problematic in Newick tree strings:
        Replaces all occurrences of chars in 'illegal' by 'rep'"""

        if self.namelist is None:
            raise TreeError("There are no items in Distmatrix object. Can not run clean_names()")
        illegal_esc_list = [re.escape(char) for char in illegal]
        illegal_esc_list.append("_") # To avoid two underscores in a row
        illegal_esc_string = "".join(illegal_esc_list)
        regex = f"[{illegal_esc_string}]+"
        for old in self.namelist:
            new = re.sub(regex,"_",old)
            old_new_tuples.append((old, new))
        for old,new in old_new_tuples:
            self.rename(old, new)

    #######################################################################################

    def rename(self, oldname, newname):
        """Changes name of one item from oldname to newname"""

        self.namelist.remove(oldname)
        self.namelist.append(newname)
        i = self.name2index[oldname]
        del self.name2index[oldname]
        self.name2index[newname] = i
        self.index2name[i] = newname

    #######################################################################################

    def getdist(self, name1, name2):
        """Returns distance between named entries"""

        i1 = self.name2index[name1]
        i2 = self.name2index[name2]
        dist = self.dmat[i1, i2]
        return dist

    #######################################################################################

    def setdist(self, name1, name2, dist):
        """Sets distance between named entries"""

        i1 = self.name2index[name1]
        i2 = self.name2index[name2]
        self.dmat[i1, i2] = self.dmat[i2, i1] = dist

    #######################################################################################

    def avdist(self):
        """Returns average dist in matrix (not including diagonal)"""
        ndist = self.n * (self.n - 1)           # diagonal entries are zero
        avdist = np.sum(self.dmat) / ndist
        return avdist

    ###############################################################################################

    def nj(self):
        """Computes neighbor joining tree, returns Tree object"""

        # Construct star-tree. This will be resolved node by node during algorithm
        njtree = Tree.from_leaves(self.namelist)
        rootnode = njtree.root

        # Local copies of object attributes that can be changed during algorithm
        # Note: unresolved_nodes maps indices from dmat to actual names
        # dmat and unresolved_nodes are updated during algorithm to always be in sync
        dmat = self.dmat.copy()
        unresolved_nodes = self.namelist.copy()

        # Main loop: continue merging nearest neighbors until only two nodes left
        while len(unresolved_nodes) > 2:
            n = len(unresolved_nodes)

            # Compute njdist = (n-2) * d(n1, n2) - udist(n1) - udist(n2)
            udist = dmat.sum(axis=0)          # udist = summed dist to all other nodes.
            udist_col = udist.reshape(n, 1)   # column vector version of udist
            njdist = (n - 2) * dmat
            njdist -= udist                   # subtract udist vector from all rows, in-place
            njdist -= udist_col               # subtract udist vector from all cols, in-place

            # Find closest neighbors according to njdist
            np.fill_diagonal(njdist, np.inf)  # set diagonal to inf (so it wont be min)
            flat_ix = np.argmin(njdist)
            i1, i2 = np.unravel_index(flat_ix, njdist.shape)
            if i1 > i2:                       # Ensure i1 < i2 (used below)
                i1, i2 = i2, i1
            nb1 = unresolved_nodes[i1]
            nb2 = unresolved_nodes[i2]

            # Connect two nearest nodes on tree (insert new node below them)
            dist_12 = dmat[i1,i2]
            udist_1 = udist[i1]
            udist_2 = udist[i2]
            branchstruct = Branchstruct()
            newnode = njtree.insert_node(rootnode, [nb1, nb2], branchstruct)
            dist1 = 0.5 * dist_12 + 0.5 * (udist_1 - udist_2) / (n - 2)
            dist2 = 0.5 * dist_12 + 0.5 * (udist_2 - udist_1) / (n - 2)
            njtree.setlength(newnode, nb1, dist1)
            njtree.setlength(newnode, nb2, dist2)

            # Update distance matrix.
            # Remove 2 merged entries (rows and cols), add 1 new row and col for newnode
            # Removal is done by deleting row and column i2, and overwriting i1 by new node
            dist_new = 0.5 * (dmat[i1] + dmat[i2] - dist_12) # distvector from new
            dmat[i1] = dmat[:,i1] = dist_new
            dmat = np.delete(dmat, i2, axis=0)
            dmat = np.delete(dmat, i2, axis=1)

            # Update list of unresolved nodes
            unresolved_nodes[i1] = newnode
            del unresolved_nodes[i2]         # Note: all indices >=i2 are shifted same as array

        # After loop. Set length of branch conecting final two nodes
        nb1, nb2 = unresolved_nodes[0], unresolved_nodes[1]
        dist = dmat[0, 1]
        njtree.deroot()
        njtree.setlength(nb1, nb2, dist)

        return njtree

    ###############################################################################################

    def upgma(self):
        """Computes UPGMA tree, returns Tree object"""

        # Construct star-tree. This will be resolved node by node during algorithm
        upgmatree = Tree.from_leaves(self.namelist)
        rootnode = upgmatree.root

        # Local copies of Distmatrix object attributes that can be changed during algorithm
        dmat = self.dmat.copy()
        remaining_nodes = self.namelist.copy()
        name2ix = self.name2index.copy()
        ix2name = self.index2name.copy()
        n = len(remaining_nodes)

        # Dicts keeping track of node depths and cluster sizes
        clus_size = dict.fromkeys(range(n), 1)      # Number of leaves in cluster (initially 1)
        depth = dict.fromkeys(range(n), 0.0)        # Depth of node (leaves = 0)

        # List keeping track of minimum value (and its index) in each row
        # Should lead to O(n^2) instead of O(n^3) time complexity (see Felsenstein)
        np.fill_diagonal(dmat, np.inf)              # set diagonal to inf (so it wont be min)
        min_colix = np.argmin(dmat, axis=1)         # col-index of smallest value in each row
        min_rowval = dmat[np.arange(n), min_colix]  # smallest value in each row

        # Main loop: continue merging nearest nodes until only two nodes left
        while len(remaining_nodes) > 2:

            # Find closest neighbors
            i1 = np.argmin(min_rowval)          # row index for overall smallest value
            i2 = min_colix[i1]                  # corresponding column index
            n1 = ix2name[i1]
            n2 = ix2name[i2]

            # Connect two nearest nodes on tree (insert new node below them)
            branchstruct = Branchstruct()
            newnode = upgmatree.insert_node(rootnode, [n1, n2], branchstruct)
            depth_12 = dmat[i1,i2]/2
            dist_new1 = depth_12 - depth[i1]
            dist_new2 = depth_12 - depth[i2]
            upgmatree.setlength(newnode, n1, dist_new1)
            upgmatree.setlength(newnode, n2, dist_new2)

            # Update distance matrix
            # Remove two merged entries (row and col), replace by row and col for newnode
            # Removal is here done by setting values to inf (i2),
            # or by overwriting previous values with new node (i1)
            dist_new = (clus_size[i1] * dmat[i1] + clus_size[i2] * dmat[i2]) / (clus_size[i1] + clus_size[i2])
            dmat[i1] = dmat[:,i1] = dist_new
            dmat[i2] = dmat[:,i2] = np.inf

            # Update lists keeping track of smallest value:
            # For each row in dmat:
            #       if smallest value was previously in one of removed cols: find new min
            #       else if new value (in dist_new) smaller than previous: replace
            #       set i2 rowinfo to inf so it will not be minimum (= delete info)
            # This complicated approach is O(n) where checking entire dmat would be O(n^2), so worth it
            merged_indexes = {i1,i2}
            for i, (prev_mincol, prev_val, new_val) in enumerate( zip( min_colix, min_rowval, dist_new ) ):
                if  prev_mincol in merged_indexes:          # Previous minimum has been changed by merge:
                    new_min_colix = np.argmin(dmat[i])      # find new minimum on this row of dmat
                    min_colix[i] = new_min_colix
                    min_rowval[i] = dmat[i,new_min_colix]
                else:
                    if prev_val > new_val:
                        min_colix[i] = i1
                        min_rowval[i] = new_val

            # Update lists and dicts keeping track of current nodes and their names
            remaining_nodes.append(newnode)
            remaining_nodes.remove(n1)
            remaining_nodes.remove(n2)
            del name2ix[n1]
            del name2ix[n2]
            del ix2name[i2]
            ix2name[i1] = newnode
            name2ix[newnode] = i1
            clus_size[i1] = clus_size[i1] + clus_size[i2]
            del clus_size[i2]
            depth[i1] = depth_12
            del depth[i2]

        # After loop. Set length of branch conecting final two nodes
        n1, n2 = remaining_nodes[0], remaining_nodes[1]
        i1 = name2ix[n1]
        i2 = name2ix[n2]
        depth_12 = dmat[i1,i2]/2
        dist_new1 = depth_12 - depth[i1]
        dist_new2 = depth_12 - depth[i2]
        upgmatree.setlength(rootnode, n1, dist_new1)
        upgmatree.setlength(rootnode, n2, dist_new2)

        return upgmatree


###################################################################################################
###################################################################################################

# Placeholder: Insert test code here and run module in standalone mode
def main():
    pass

###################################################################################################


if __name__ == "__main__":
    main()

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

def remove_comments(text, leftdelim, rightdelim=None):
    """Takes input string and strips away commented text, delimited by 'leftdelim' and 'rightdelim'.
        Also deals with nested comments."""

    # NOTE: only deals with block comments at present
    # Python note: maybe this is too general. Will I ever use multichar delims?
    def wordsoverlap(w1, w2):
        for i in range(1, len(w2)):
            if w1.startswith(w2[i:]):
                return True
        for i in range(1, len(w1)):
            if w2.startswith(w1[i:]):
                return True
        return False

    if leftdelim == rightdelim:
        raise SeqError("Left and right delimiters are identical")
    elif leftdelim in rightdelim:
        raise SeqError("Left delimiter is substring of right delimiters")
    elif rightdelim in leftdelim:
        raise ExcepSeqErrortion("Right delimiter is substring of left delimiters")
    elif wordsoverlap(leftdelim, rightdelim):
        raise SeqError("Right and left delimiters overlap")

    # Preprocess delims for use in re etc
    leftdelim = re.escape(leftdelim)
    rightdelim = re.escape(rightdelim)

    # Construct sorted list of tuples of the form [(0, 'start'), (5, 'stop'), (7, 'start'), ...]
    delimlist = [(match.start(), match.end(), "start") for match in re.finditer(leftdelim, text)]
    # If text contains no starts (=> no comments): return un-altered text
    if not delimlist:
        return text
    else:
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
                raise Exception("Unmatched end-comment delimiter. Context: '{}'".format(text[prevpos-10:prevpos+10]))

    # Add final block of text if relevant (i.e., if text does not stop with rightdelim), return processed text
    if prevpos < len(text):
        processed_text.append(text[prevpos:])
    return "".join(processed_text)

###################################################################################################
###################################################################################################

class Globals():
    """Class containing globally used functions and labels."""

    # I'm not convinced this is the way to go. Module instead?"""
    # Global repository for bipartitions, to avoid redundant saving in topologies etc.
    biparts = {}

###################################################################################################
###################################################################################################

class Interner():
    """Class used for interning various objects."""

    # Stores dictionaries of leafset, bipartitions, and topologies
    # Interner methods returns *pointer* to leafset or bipartition
    # Could perhaps just use one dict for interning *anything*, but worry about hash collisions?

    def __init__(self):
        self.leafsets = {}
        self.biparts = {}
        self.topo = {}

    def intern_leafset(self, leafset):
        if leafset not in self.leafsets:
            self.leafsets[leafset]=leafset
        return self.leafsets[leafset]

    def intern_bipart(self, bipart):
        if bipart not in self.biparts:
            self.biparts[bipart]=bipart
        return self.biparts[bipart]

    def intern_topology(self, topology):
        if topology not in self.topo:
            self.topo[topology]=topology
        return self.topo[topology]


###################################################################################################
###################################################################################################

class Branchstruct():
    """Class that emulates a struct. Keeps branch-related info"""

    def __init__(self, length=0.0, label=""):
        self.length = length
        self.label = label

    ###############################################################################################

    def copy(self):
        """Returns copy of Branchstruct object, with all attributes included"""

        # Python note: reconsider copying attributes besides blen and lab?
        obj = Branchstruct()
        for attrname, value in vars(self).items():
            setattr(obj, attrname, value)
        return obj

###################################################################################################
###################################################################################################

class Topostruct():
    """Class that emulates a struct. Keeps topology-related info"""

    __slots__ = ["weight", "tree", "freq"]

    # Python note: perhaps replace with dataclass, available since python 3.7
    pass

###################################################################################################
###################################################################################################

class TreeError(Exception):
    pass

###################################################################################################
###################################################################################################

class Tree():
    """Class representing basic phylogenetic tree object."""

    # Implementation note: Tree objects can be constructed from several different kinds of things:
    # including Newick tree strings, Bipartition lists, and a list of leaves.
    # The Tree class therefore has several alternate constructors implemented as classmethods
    # The main constructor "__init__" is therefore mostly empty

    def __init__(self):
        self._parent_dict = None         # Dict node:parent relationships (only built if required)
        self.dist_dict = None
        self.path_dict = None
        self.sorted_intnode_cache = None

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
            for child in self.tree[parent]:
                self._parent_dict[child] = parent
        self._parent_dict[self.root] = None     # Add special value "None" as parent of root

    ###############################################################################################

    @classmethod
    def from_string(cls, orig_treestring, transdict=None):
        """Constructor: Tree object from tree-string in Newick format"""

        obj = cls()

        # NOTE: interprets non-leaf labels as belonging to an internal branch (not to
        # an internal node). The label is attached to the same branch as the branch length

        # Remove whitespace (string methods are much faster than regexp.sub)
        treestring = orig_treestring.replace(" ", "")
        treestring = treestring.replace("\t", "")
        treestring = treestring.replace("\n", "")
        treestring = treestring.replace("\r", "")
        treestring = treestring.replace("\f", "")

        # Sanity check: number of left- and right-parentheses should match
        if treestring.count("(") != treestring.count(")"):
            msg = "Imbalance in tree-string: different number of left- and right-parentheses\n"
            msg += "Left: ({0}  Right: {1})".format(treestring.count("("), treestring.count(")"))
            raise TreeError(msg)

        # Break treestring up into a list of parentheses, commas, names, and branch lengths
        # ("tokenize")
        # Python note: characters that are not one of the below, will be quietly discarded
        # this includes things such as quotes, ampersands, dollarsigns etc.
        # possibly this is a good idea, but may cause trouble (figtree format for instance)
        tree_parts = re.compile(r"""                # Save the following sub-patterns:
                \(          |                       # (1) a left parenthesis
                \)          |                       # (2) a right parenthesis
                ,           |                       # (3) a comma
                ;           |                       # (4) a semicolon
                :[+-]?\d*\.?\d+(?:[eE][-+]?\d+)? |  # (5) a colon followed by a branch length
                                                    #     possibly negative,
                                                    #     possibly using exponential notation
                [\w\-\/\.\*\|]+                  |  # (6) a name/label (one or more alphanumeric)
                \[.*?\]                             # (7) a bracketted comment
                                                    #     typically from FigTree
                                                    #     placed as branch label
                """, re.VERBOSE)
        tree_parts_list = tree_parts.findall(treestring)
        # Tree is represented as a dictionary of dictionaries. The keys in the top dictionary
        # are the internal nodes which are numbered consecutively. Each key has an
        # associated value that is itself a dictionary listing the children: keys are
        # child nodes, values are Branchstructs containing "length" and "label" fields.
        # Leafs are identified by a string instead of a number
        # NOTE: self.tree should ONLY be used by this object. Make pseudo-private?
        obj.tree = {}                  # Essentially a "child-list"
        obj.leaves = set()             # Set of leaf names. For speedy lookups
        obj.intnodes = set()           # Set of internal node IDs. For speedy lookups
        obj.root = 0                   # Root is node zero at start. (May change)
        obj.belowroot = None

        # Preprocess parts list to remove and store any labels or lengths below root node
        # Any text here (between the semicolon and the last right parenthesis)
        # will be stores as _one_ string in "self.belowroot"
        i = -2
        extraparts = []
        while tree_parts_list[i] != ")":
            extraparts.append(tree_parts_list[i])
            i -= 1
        if extraparts:
            extraparts.reverse()
            obj.belowroot = "".join(extraparts)
            del tree_parts_list[(i+1):-1]

        # Parse tree_parts_list left to right. Use a stack to keep track of current node
        nodeno = -1
        node_stack = []

        for part in tree_parts_list:

            # A left parenthesis indicates that we are about to examine a new internal node.
            if part == "(" :
                nodeno += 1
                obj.tree[nodeno] = {}          # Create child-list for new node

                if nodeno != 0:                 # If we are not at the root then add new node
                    parent = node_stack[-1]     # to previous node's list of children
                    obj.tree[parent][nodeno] = Branchstruct()

                node_stack.append(nodeno)       # Push new node onto stack
                obj.intnodes.add(nodeno)       # Add node to list of internal nodes

            # A right parenthesis indicates that we have finished examining a node
            # I should maybe catch the "IndexError: pop from empty list" somewhere
            elif part == ")":
                del node_stack[-1]              # Remove last examined node from stack

            # A colon indicates this is a distance
            elif part[0] == ":":
                dist = float(part[1:])                  # Remove colon and convert to float
                child = node_stack[-1]
                parent = node_stack[-2]
                obj.tree[parent][child].length =  dist # Add dist to relevant child-list

            # A comma indicates that we have finished examining a node
            elif part == ",":
                del node_stack[-1]                      # Remove last examined node from stack

            # A semicolon indicates the end of the treestring has been reached
            elif part == ";":
                del node_stack[-1]                      # Clean up stack

            # If nothing else matched then this must be a name or label
            # If previous part was simultaneously a right parenthesis or a leaf name,
            # then the token must be a branch label
            elif prevpart in [")", "leafname"]:
                child = node_stack[-1]
                parent = node_stack[-2]
                obj.tree[parent][child].label = part

            # Last possibility (I hope): the name is a leafname
            else:
                # If translation dictionary was supplied: change name accordingly
                if transdict:
                    child = sys.intern(transdict[part])  # Already interned. Not sure this is needed...
                else:
                    child = sys.intern(part)    # "Intern" strings so only one copy in memory
                                                # (may happen automatically anyway?...)

                # Check for duplicated leaf names
                if child in obj.leaves:
                    msg = "Leaf name present more than once: {}".format(child)
                    raise TreeError(msg)

                parent=node_stack[-1]                      # Add new leaf node to previous node's
                obj.tree[parent][child] = Branchstruct()   # list of children

                node_stack.append(child)                   # Push new leaf node onto stack
                obj.leaves.add(child)                      # Also update list of leaves

                part = "leafname"

            prevpart = part

        # If nodes remain on stack after parsing, then something was wrong with tree-string
        if node_stack:
            msg = "Imbalance in tree-string: '%s'" % orig_treestring
            raise TreeError(msg)


        obj.nodes = obj.leaves | obj.intnodes
        return obj

    ###############################################################################################
    @classmethod
    def from_biplist(cls, biplist):
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
        obj.tree = {}
        obj.tree[0]={}
        maxnode = 0

        # Start by building star-tree, this is resolved branch-by-branch later on
        # Python note: use startree constructor?
        for leaf in obj.leaves:
            obj.tree[0][leaf]= None

        # Iterate over all bipartitions, for each: add extra branch and/or update Branchstruct
        for (bip1, bip2), branchstruct in biplist.items():

            # If bipartition represents external branch: update relevant Branchstruct
            if len(bip1) == 1 or len(bip2) == 1:
                if len(bip1) == 1:
                    (leaf, ) = bip1      # A one-member tuple used for value unpacking (pop)
                else:
                    (leaf, ) = bip2

                # Find childdict containing leaf, and update branchstruct
                for childdict in obj.tree.values():
                    if leaf in childdict:
                        childdict[leaf] = branchstruct
                        break

            # If bipartition represents internal branch: add branch to tree, transfer Branchstruct
            else:
                mrca1 = obj.find_mrca(bip1)
                mrca2 = obj.find_mrca(bip2)

                # Determine which group of leaves to move
                # Note: one part of bipartition will necessarily have root as its MRCA
                # since the members are present on both sides of root. It is the other part of
                # the bipartition (where all members are on same side of root) that should be moved
                # For a star-tree the resolution will be random (both have root as their MRCA)
                if mrca1 == 0:                  # If mrca is root, move other group
                    insertpoint = mrca2
                    active_bip = bip2
                else:
                    insertpoint = mrca1
                    active_bip = bip1

                # Determine which of insertpoints children to move, namely all children
                # that are either in the active bipartition OR children whose descendants are
                moveset = []
                for child in obj.children(insertpoint):
                    if child in active_bip:
                        moveset.append(child)
                    elif (child not in obj.leaves) and (obj.remote_children(child) <= active_bip):
                        moveset.append(child)

                # Construct new internal node (and therefore branch), transfer Branchstruct
                maxnode += 1
                obj.tree[maxnode] = {}
                obj.tree[insertpoint][maxnode] = branchstruct
                obj.intnodes.add(maxnode)       # Add new node to list of internal nodes

                # Move relevant children to new node, transfer Branchstructs
                for child in moveset:
                    obj.tree[maxnode][child] = obj.tree[insertpoint][child]
                    del obj.tree[insertpoint][child]

        obj.nodes = set(obj.leaves | obj.intnodes)
        return obj

    ###############################################################################################
    @classmethod
    def from_topology(cls, topology):
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
        obj.tree = {}
        obj.tree[0]={}
        maxnode = 0

        # Start by building star-tree, this is resolved branch-by-branch later on
        for leaf in obj.leaves:
            obj.tree[0][leaf]= Branchstruct()

        # Iterate over all bipartitions, for each: add extra branch and/or update Branchstruct
        for bip1, bip2 in topology:

            # If bipartition represents internal branch: add branch to tree
            if len(bip1) > 1 and len(bip2) > 1:
                mrca1 = obj.find_mrca(bip1)
                mrca2 = obj.find_mrca(bip2)

                # Determine which group of leaves to move
                # Note: one part of bipartition will necessarily have root as its MRCA
                # since the members are present on both sides of root. It is the other part of
                # the bipartition (where all members are on same side of root) that should be moved
                # For a star-tree the resolution will be random (both have root as their MRCA)
                if mrca1 == 0:                  # If mrca is root, move other group
                    insertpoint = mrca2
                    active_bip = bip2
                else:
                    insertpoint = mrca1
                    active_bip = bip1

                # Determine which of insertpoints children to move, namely all children
                # that are either in the active bipartition OR children whose descendants are
                moveset = []
                for child in obj.children(insertpoint):
                    if child in active_bip:
                        moveset.append(child)
                    elif (child not in obj.leaves) and (obj.remote_children(child) <= active_bip):
                        moveset.append(child)

                # Construct new internal node (and therefore branch)
                maxnode += 1
                obj.tree[maxnode] = {}
                obj.tree[insertpoint][maxnode] = Branchstruct()
                obj.intnodes.add(maxnode)       # Add new node to list of internal nodes

                # Move relevant children to new node, transfer Branchstructs
                for child in moveset:
                    obj.tree[maxnode][child] = obj.tree[insertpoint][child]
                    del obj.tree[insertpoint][child]

        obj.nodes = set(obj.leaves | obj.intnodes)
        return obj

    ###############################################################################################
    @classmethod
    def from_leaves(cls, leaflist):
        """Constructor: star-tree object from list of leaves"""

        treelist = ["("]
        for name in leaflist:
            treelist.append(name)
            treelist.append(",")
        del treelist[-1]
        treelist.append(");")
        return cls.from_string("".join(treelist))

    ###############################################################################################
    @classmethod
    def from_branchinfo(cls, parentlist, childlist, lenlist=None, lablist=None):
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
        obj.tree = {}
        obj.leaves = set()
        obj.intnodes = set()

        for i in range(nbranches):
            parent = parentlist[i]     # Perhaps check types are OK?
            child = childlist[i]
            blen = lenlist[i]
            lab = lablist[i]
            if parent in obj.tree:
                obj.tree[parent][child] = Branchstruct(blen, lab)
            else:
                obj.tree[parent] = { child:Branchstruct(blen, lab) }
            obj.intnodes.add(parent)
            if isinstance(child, str):
                obj.leaves.add(child)

        # Root node is the parent node that is not also in childlist
        diffset = set(parentlist) - set(childlist)
        obj.root = diffset.pop()

        obj.nodes = obj.leaves | obj.intnodes
        return obj

    ###############################################################################################
    @classmethod
    def randtree(cls, leaflist=None, ntips=None, randomlen=False, name_prefix="s"):
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
            tree = cls.from_leaves(leaflist)

        # If ntips given:
        #   construct list of zeropadded, numbered, names,
        #   then construct startree using names, finally resolve to random topology
        else:
            ndigits = len(str(ntips))          # Number of digits required to write max taxon number
            namelist = []
            for i in range(ntips):
                name = "{prefix}{num:0{width}d}".format(prefix=name_prefix, num=i, width=ndigits)
                namelist.append( name )
            tree = cls.from_leaves(namelist)   # Star tree with given number of leaves

        tree.resolve()                         # Randomly resolve to bifurcating tree

        if randomlen:
            for parent in tree.intnodes:
                for child in tree.tree[parent]:
                    tree.tree[parent][child].length = random.lognormvariate(math.log(0.2), 0.3)

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
                dist = "{num:.6g}".format(num=self.tree[node][kid].length)
                label = self.tree[node][kid].label
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

    def copy_treeobject(self, copylengths=True, copylabels=True):
        """Returns copy of Tree object. Copies structure and branch lengths.
        Caches and any user-added attributes are not copied.
        Similar to effect of copy.deepcopy but customized and much faster"""

        # Python note: reconsider copying attributes besides blen and lab?
        obj = Tree()
        obj.root = self.root
        obj.leaves = self.leaves.copy()
        obj.intnodes = self.intnodes.copy()
        obj.nodes = self.nodes.copy()
        obj.tree = {}
        origtree = self.tree
        newtree = obj.tree
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
        tree = self.tree
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
        tree = self.tree
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

        self.sorted_intnode_cache = sorted_nodes

        if not deepfirst:
            sorted_nodes.reverse()

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
            return set(self.tree[parent].keys())
        except KeyError as err:
            msg = "Node %s is not an internal node" % parent
            raise TreeError(msg) from err

    ###############################################################################################
    #@functools.lru_cache(maxsize=None)
    def remote_children(self, parent):
        """Returns set containing all leaves that are descendants of parent"""

        # If "parent" is a leaf, then return a set consisting of only itself
        if parent in self.leaves:
            return {parent}

        # Traverse the tree iteratively to find remote children:
        #       if kid is leaf: add it to list of remote children.
        #       if kid is intnode: push its children on stack
        kidstack = set( self.tree[parent] )
        remotechildren = set()

        while kidstack:
            curnode = kidstack.pop()
            if curnode in self.leaves:
                remotechildren.add(curnode)
            else:
                kidstack.update( self.tree[curnode] )

        return remotechildren

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
        """Finds Most Recent Common Ancestor for the provided set of leaves"""

        leafset = set(leaves)
        if not leafset <= self.leaves:
            # Construct string listing all entries in leafset (for error message)
            stringlist = []
            for leaf in leafset:
                stringlist.append(str(leaf))
                stringlist.append(", ")
            leafstring = "".join(stringlist)
            leafstring = leafstring[:-2]    # Remove trailing comma and blank

            msg = "Some nodes in set are not part of tree: %s" % leafstring
            raise TreeError(msg)

        min_numkids = len(self.leaves)
        mrca = self.root
        for node in self.sorted_intnodes(deepfirst=False):
            remkids = self.remote_children(node)
            if remkids == leafset:
                mrca = node
                break
            if leafset <= remkids:
                numkids = len(remkids)
                if numkids < min_numkids:
                    mrca = node
                    min_numkids = numkids

        return mrca

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
        tree = self.tree

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
                treelength += self.tree[node][child].length

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

                branchstruct = self.tree[parentnode][child]
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

    def bipdict(self, interner=None):
        """Returns tree in the form of a "bipartition dictionary" """

        # Names of leaves on one side of a branch are represented as an immutable set
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree is represented as a dictionary where the keys are bipartitions
        # The values are Branchstructs
        # Interning: store bipartitions in global dict to avoid duplicating object
        # This is probably mostly useful when creating toposummaries in BigTreeSummary

        bipartition_dict = {}
        if interner:
            leaves = interner.intern_leafset(frozenset(self.leaves))
        else:
            leaves = frozenset(self.leaves)

        # For each branch: find bipartition representation, add this and Branchstruct to list.
        # Remote kids of node most distant from root (or node itself) forms one part of bipartition
        # Other part is then found as diff between all leaves and bipart1
        # Python note: Sorting pays off because remote_children cache is built in rational order
        sortedintnodes = self.sorted_intnodes(deepfirst=False)

        # DEBUG: I think root should be removed from sortedintnodes (dealt with separately below)
        # DEBUG 2: interning in global dict does help with toposummary, but may cause problems
        # because there is a reference to all bipartitions regardless of what happens to them
        # which means they cant be garbage collected.
        for node1 in sortedintnodes:
            for node2 in self.children(node1):
                if interner:
                    bipart1 = interner.intern_leafset(frozenset(self.remote_children(node2)))
                    bipart2 = interner.intern_leafset(leaves - bipart1)
                    bipartition = interner.intern_bipart(frozenset([bipart1, bipart2]))
                else:
                    bipart1 = frozenset(self.remote_children(node2))
                    bipart2 = leaves - bipart1
                    bipartition = frozenset([bipart1, bipart2])

                # Interning:
                # Bipartitions can be kept in global dict to avoid redundant saving (they can be huge)
                # This is especially useful in connection with toposummaries, where potentially
                # tens of thousands of trees all contain tens of bipartitions, but where a large
                # fraction of the bipartitions are used in most trees).
                bipartition_dict[bipartition] = self.tree[node1][node2]

        # If root is attached to exactly two nodes, then two branches correspond to the same
        # bipartition. Clean up by collapsing two branches (add lengths, compare labels)
        root = self.root
        rootkids = self.children(root)

        if len(rootkids) == 2:

            # First: find out what bipartition root is involved in.
            kid1, kid2 = rootkids
            if interner:
                bipart1 = interner.intern_leafset(frozenset(self.remote_children(kid1)))
                bipart2 = interner.intern_leafset(leaves - bipart1)
                bipartition = interner.intern_bipart(frozenset([bipart1, bipart2]))
            else:
                bipart1 = frozenset(self.remote_children(kid1))
                bipart2 = leaves - bipart1
                bipartition = frozenset([bipart1, bipart2])

            # Create new branch
            bipartition_dict[bipartition] = Branchstruct(self.tree[root][kid1].length,
                                                         self.tree[root][kid1].label)

            # Now, add distance to other kid
            bipartition_dict[bipartition].length += self.tree[root][kid2].length

            # Deal with labels intelligently
            lab1 = self.tree[root][kid1].label
            lab2 = self.tree[root][kid2].label

            # If only one label is set use  that.
            if (lab1 is not None) and (lab2 is None):
                lab = lab1
            elif (lab1 is None) and (lab2 is not None):
                lab = lab2
            # If both or none of the labels are set: pick label1 randomly
            else:
                lab = lab1

            bipartition_dict[bipartition].label = lab

        return bipartition_dict

    ###############################################################################################

    def topology(self):
        """Returns set of sets of sets representation of topology ("naked bipdict")"""

        # Names of leaves on one side of a branch are represented as an immutable set.
        # A bipartition is represented as an immutable set of two such (complementary) sets
        # The entire tree topology is represented as a set of bipartitions
        # This is essentially a naked version of a bipdict

        bipdict = self.bipdict()
        topology = frozenset(bipdict.keys())

        return topology

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

        # In all other cases: at least one part of bipartition will necessarily have root as its MRCA
        # (because its members are present on both sides of root).
        # It is the other part of the bipartition that can potentially be moved
        # (Occasionally both parts have root as MRCA, but this will be covered by the more general
        # solution below)
        mrca1 = self.find_mrca(bip1)
        mrca2 = self.find_mrca(bip2)

        if mrca1 == self.root:
            insertpoint = mrca2
            active_bip = bip2
        else:
            insertpoint = mrca1
            active_bip = bip1

        # Determine which of insertpoint's children to move (namely all children
        # that are either in the active bipartition OR children whose descendants are)
        movelist = []
        moveable_descendants = []
        for child in self.children(insertpoint):
            if child in active_bip:
                movelist.append(child)
                moveable_descendants.append(child)
            else:
                remkids = self.remote_children(child)
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
            numkids = len(self.tree[node])   # Note: not safe to use .children() method while
                                             # changing tree (cache will break)
            if numkids > 2:
                unresolved_nodes.append(node)

        # Keep adding extra internal nodes until there are no unresolved nodes left
        while unresolved_nodes:
            intnode1 = unresolved_nodes.pop()
            kids = self.tree[intnode1].keys()
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

    def setlength(self, node1, node2, length):
        """Sets length of branch connecting node1 and node2"""

        if node1 == self.parent(node2):
            parent = node1
            child = node2
        elif node2 == self.parent(node1):
            parent = node2
            child = node1
        else:
            msg = "There is no branch connecting node {} and {}".format(node1, node2)
            raise TreeError(msg)

        self.tree[parent][child].length = length

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

        self.tree[parent][child].label = label

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

        return self.tree[parent][child].label

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
            basalbranch = self.tree[parent][basenode]

        # Special case: basenode is leaf => subtree is minimal tree with two nodes (root and leaf)
        if basenode in self.leaflist():
            other = Tree.from_string("({});".format(basenode))

        # If basenode is internal: subtree has more than one leaf
        else:
            # Create empty Tree object. Transfer relevant subset of self's data structure to other
            other = Tree()
            other.tree = {}
            other.intnodes = {basenode}
            other.leaves = set()
            other.root = basenode
            curlevel = [basenode]
            while curlevel:
                nextlevel = []
                for parent in curlevel:
                    other.tree[parent] = {}
                    kids = self.children(parent)
                    for kid in kids:
                        other.tree[parent][kid] = copy.deepcopy(self.tree[parent][kid])
                    intnode_kids = kids & self.intnodes
                    other.intnodes.update(intnode_kids)
                    nextlevel.extend(intnode_kids)
                    leaf_kids = kids & self.leaves
                    other.leaves.update(leaf_kids)
                curlevel = nextlevel
            other.nodes = other.leaves | other.intnodes

        # Python note: possibly bad idea to have different possible returnvalues.
        # Simplify and deal with it at consumer end
        if return_basalbranch:
            return (other, basalbranch)
        else:
            return other

    ###############################################################################################

    def graft(self, other, node1, node2=None, blen1=0, blen2=0, graftlabel=None):
        """Graft other tree to self

        tree2 (other) intnodes will be renamed if names clash with those in tree1.
        node1: node in tree1 (self) below which tree2 (other) will be grafted. Cannot be root1
        node2: node in tree2 (other) below which tree2 will be attached (default is root of tree2)
        blen1: length of branch added to tree1 below graftpoint (lower of two newly created branches)
        blen2: length of branch above graft point and below tree2 (upper of two newly created branches)
        graftlabel: prepend value of "label" to leaf names on t2 (e.g: "graft_s1")"""

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

        # Update main data structure (self.tree dictionary) by merging with dict from other
        self.tree.update(other.tree)
        # Link subtree to graftpoint in self.tree
        self.tree[graftpoint][other.root] = Branchstruct(length=blen2)

        # Update look-up lists and caches
        self.nodes.update( other.nodes )
        self.intnodes.update( other.intnodes )
        self.leaves.update( other.leaves )

        # Reset caches and lists
        self.sorted_intnode_cache = None
        self._parent_dict = None
        self.sorted_intnode_cache = None
        self.nodedist.cache_clear()


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
        # Could probably be sped up by using information already in .tree dict (thereby essentially
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
        for parent, kid_dict in self.tree.items():
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
        for parent, kid_dict in self.tree.items():
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
        tree = self.tree

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

        # Erase caches
        self.sorted_intnode_cache = None

        # Clear lru_caches (which cannot be edited manually)
        #self.remote_children.cache_clear()
        self.nodedist.cache_clear()

        return newnode

    ###############################################################################################

    def add_branch(self, bipart, branchstruct):
        """Adds branch represented by bipartition to unresolved tree."""

        # NOTE: huge overlap with TreeFromBiplist - should use this function there as well!

        # Sanity check: is bipartition compatible with tree?
        if not self.is_compatible_with(bipart):
            raise TreeError("Bipartition is not compatible with tree: %s" % bipart)

        part1, part2 = bipart
        mrca1 = self.find_mrca(part1)
        mrca2 = self.find_mrca(part2)

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
            if mrca1 == self.root:      # If mrca1 is root, insert at mrca2, and move part2
                insertpoint = mrca2
                active_bip = part2
            else:                       # If mrca2 is root, insert at mrca1, and move part1
                insertpoint = mrca1
                active_bip = part1

            # Determine which of insertpoint's children to move (namely all children
            # that are either in the active bipartition OR children whose descendants are)
            movelist = []
            for child in self.children(insertpoint):
                if child in active_bip:
                    movelist.append(child)
                elif (child not in self.leaves) and (self.remote_children(child) <= active_bip):
                    movelist.append(child)

            # Add branch at determined position
            self.insert_node(insertpoint, movelist, branchstruct)

        # Clear lru_caches (which cannot be edited manually)
        #self.remote_children.cache_clear()
        self.nodedist.cache_clear()

    ###############################################################################################

    def remove_branch(self, node1, node2):
        """Removes branch connecting node1 and node2 (thereby creating polytomy)"""

        # Python note: Length of removed branch is lost: treelength and patristic distances change
        # Should I distribute lost branch length to retain treelength?
        # Eg, retain treelength and patristic distances within subtree (+ blen/n to each child)
        # Alternatively retain patristic distances to outside (+ blen to each child),
        # but then treelength would increase

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
        for grandchild in self.children(child):
            self.tree[parent][grandchild] = Branchstruct(self.tree[child][grandchild].length,
                                                         self.tree[child][grandchild].label)

        # Delete "child" node and link from parent to child. Update intnodes and nodes
        del self.tree[child]
        del self.tree[parent][child]
        self.intnodes.remove(child)
        self.nodes = self.leaves | self.intnodes

        # Python note: would be simple to just update _parent_dict instead of setting to None
        self._parent_dict = None

        # Update self.sorted_intnode_cache if it exists
        if self.sorted_intnode_cache is not None:
            self.sorted_intnode_cache.remove(child)

        # Clear lru_caches (which cannot be edited manually)
        #self.remote_children.cache_clear()
        self.nodedist.cache_clear()

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
            del self.tree[root]                             # Remove entry for old root
            self.intnodes.remove(root)
            self.nodes.remove(root)
            self.root = child2                              # child2 is new root
            del self._parent_dict[child2]                   # clean up parent_dict: Note: not lazy? change?

        # If leaf is part of bifurcation but NOT attached directly to root, then parent
        # must also be removed from tree, and the remaining child needs to be grafted
        # onto grandparent with proper cumulated distance
        # Only the branch label (if any) of the internal branch is kept
        elif len(childset) == 2:
            [child2] = childset - {leaf}                      # Remaining item is other child
            child2dist = self.tree[parent][child2].length   # Remember dist to other child
            grandparent = self.parent(parent)

            # Add remaining child to grandparent
            # self.tree[grandparent][leaf2] = Branchstruct(self.tree[grandparent][parent].length,
            #                                              self.tree[grandparent][parent].label)
            self.tree[grandparent][child2] = self.tree[grandparent][parent]
            self.tree[grandparent][child2].length += child2dist   # Cumulated distance
            del self.tree[parent]                           # Delete parent and leaf
            del self.tree[grandparent][parent]              # Also remove pointer from gp to p
            del self.parent_dict[leaf]                      # Remove unused entries in parent_dict
            del self.parent_dict[parent]
            self.parent_dict[child2] = grandparent           # Update parent_dict for leaf2
            self.intnodes.remove(parent)
            self.nodes.remove(parent)

        # If leaf is part of multifurcation, then no special cleanup needed
        else:
            del self.tree[parent][leaf]
            del self.parent_dict[leaf]

        # Remove leaf entry from global leaflist. Update intnodeslist
        self.leaves.remove(leaf)
        self.nodes.remove(leaf)

        # Erase self.sorted_intnode_cache if it exists (could I salvage something?)
        self.sorted_intnode_cache = None

        # Clear lru_caches (which cannot be edited manually)
        #self.remote_children.cache_clear()
        self.nodedist.cache_clear()

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
                    if nodedist(parent,leaf) > maxdist:
                        maxdist = nodedist(parent,leaf)
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
        self.tree[parent][newname] = self.tree[parent][oldname]
        del self.tree[parent][oldname]

        # Update self.leaves and self.nodes
        self.leaves.add(newname)
        self.leaves.remove(oldname)
        self.nodes.add(newname)
        self.nodes.remove(oldname)

        # # Update self.parent_dict if it exists:
        # if self.parent_dict is not None:
        #     self.parent_dict[newname] = self.parent_dict[oldname]
        #     del self.parent_dict[oldname]

        # Clear lru_caches (which cannot be edited manually)
        #self.remote_children.cache_clear()
        self.nodedist.cache_clear()

    ###############################################################################################

    def rename_intnode(self, oldnum, newnum):
        """Changes number of one internal node"""

        if oldnum not in self.intnodes:
            msg = "Internal node {} does not exist".format(oldnum)
            raise TreeError(msg)

        if newnum in self.intnodes:
            msg = "There is already an internal node with the number {}".format(oldnum)
            raise TreeError(msg)

        # Make a note of original's parent and children
        kidlist = self.children(oldnum)
        parent = self.parent(oldnum)            # Will be None if oldnum is root

        # Update main data structure (child list in self.tree)
        self.tree[newnum] = {}
        for child in kidlist:
            self.tree[newnum][child] = self.tree[oldnum][child]
        del self.tree[oldnum]
        if oldnum == self.root:
            self.root = newnum
        else:
            self.tree[parent][newnum] = self.tree[parent][oldnum]
            del self.tree[parent][oldnum]

        # Update look-up lists, caches, and root-marker if relevant
        self.nodes.add(newnum)
        self.nodes.remove(oldnum)
        self.intnodes.add(newnum)
        self.intnodes.remove(oldnum)
        if self.sorted_intnode_cache is not None:
            i = self.sorted_intnode_cache.index(oldnum)
            self.sorted_intnode_cache[i] = newnum
        self._parent_dict = None
        # if parent is not None:
        #     self.parent_dict[newnum] = self.parent_dict[oldnum]
        #     del self.parent_dict[oldnum]
        # for child in kidlist:
        #     self.parent_dict[child] = newnum

        # Clear lru_caches (which cannot be edited manually)
        #self.remote_children.cache_clear()
        self.nodedist.cache_clear()

    ###############################################################################################

    def treedist(self, other, normalise=True, verbose=False):
        """Compute symmetric tree distance (Robinson Foulds) between self and other tree.
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

        symdif, symdif_norm, n_shared, n_uniq1, n_uniq2, n_bip1, n_bip2 = self.treedist(other, verbose=True)
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
            distsum = self.tree[root][kid1].length + self.tree[root][kid2].length

            # Deal with labels semi-intelligently
            # If only one label set: use that. Otherwise pick lab1
            lab1 = self.tree[root][kid1].label
            lab2 = self.tree[root][kid2].label
            if lab1 != "" and lab2 == "":
                lab = lab1
            elif lab1 == ""  and lab2 != "":
                lab = lab2
            else:
                lab = lab1                  # If agree pick 1, if disagree: pick 1 randomly

            # Python note: should handle situation where Branchstruct has additional attributes
            # Use introspection to find attributes and combine intelligently?
            if kid1 in self.intnodes:
                self.tree[kid1][kid2] = Branchstruct(length = distsum, label = lab)
                self.root = kid1
            elif kid2 in self.intnodes:
                self.tree[kid2][kid1] = Branchstruct(length = distsum, label = lab)
                self.root = kid2
            else:
                raise TreeError("Cannot deroot tree: only leaves are present")

            # remove previous root
            del self.tree[root]

            # update intnodelist and delete caches, which are now unreliable
            self.intnodes = set(self.tree.keys())
            #self.remote_children.cache_clear()
            self.nodedist.cache_clear()
            self.sorted_intnode_cache = None
            self.dist_dict = None
            self._parent_dict = None

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
            newbranch = self.tree[parent][child].copy()
            newbranch.length = parent_to_root_dist
            newroot = self.insert_node(parent, [child], newbranch)
            self.tree[newroot][child].length = root_to_child_dist

        # Things that were already downstream of newroot do not need to be moved, but things that
        # were previously upstream need to be moved downstream, which is done by reversing the
        # links on the direct path going back from newroot to oldroot
        oldroot = self.root
        path_to_old = self.nodepath(newroot, oldroot)
        for i in range( len(path_to_old) - 1 ):
            newparent, oldparent = path_to_old[i], path_to_old[i+1]
            self.tree[newparent][oldparent] = self.tree[oldparent][newparent]
            del self.tree[oldparent][newparent]
            # self.parent_dict[oldparent] = newparent

        # Clean up: clear caches, which are now unreliable
        # Python note: replace with lazy evaluation, using @property?
        #self.remote_children.cache_clear()
        self.nodedist.cache_clear()
        self.sorted_intnode_cache = None
        self._parent_dict = None

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

    def spr(self,subtree_node, regraft_node):
        """Subtree Pruning and Regrafting.

            subtree_node: basenode of subtree that will be pruned.
            regraft_node: node in tree below which subtree will be grafted"""

        # Subtree Prune:
        # Create Tree object corresponding to subtree. Remove subtree from self one leaf at a time
        # Python note: should check that regraft node is in self after removing subtree
        subtree = self.subtree(subtree_node)
        for leaf in self.remote_children(subtree_node):
            self.remove_leaf(leaf)

        # Regraft: Add subtree back onto remaining tree
        self.graft(subtree, regraft_node)

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

class TreeSummary():
    """Class summarizing bipartitions and branch lengths (but not topologies) from many trees"""

    def __init__(self, interner=None):
        """TreeSummary constructor. Initializes relevant data structures"""
        self.transdict = None
        self.translateblock = None
        self.tree_count = 0
        self.tree_weight_sum = 0.0
        self._bipartsummary = {}         # Dict: {bipartition:branchstruct with extra fields}
        self._bipartsummary_processed = False
        self._sorted_biplist = None
        if not interner:
            self.interner = Interner()
        else:
            self.interner = interner    # To enable sharing interned bipartition info between
                                        # different TreeSummary instances

    ###############################################################################################

    def __len__(self):
        return self.tree_count

    ###############################################################################################

    @property
    def bipartsummary(self):
        """Property method for lazy evaluation of freq, var, and sem for branches"""
        if not self._bipartsummary_processed:
            for bipartition, branch in self._bipartsummary.items():
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
    def sorted_biplist(self):
        """Return list of bipartitions.
        First external (leaf) bipartitions sorted by leafname.
        Then internal bipartitions sorted by freq"""

        if self._sorted_biplist == None:
            leafbips = []
            internalbips = []

            for bip, branch in self.bipartsummary.items():
                bip1,bip2 = bip
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

        self._bipartsummary_processed = False
        self._sorted_biplist = None

        # Main interface to TreeSummary.
        # Takes tree object, updates relevant measures
        # First time entered: build set of leaves for consistency checking.
        # Also compute transdict and translateblock for tree reporting (and storage?)
        if self.tree_count == 0:
            self.leaves = curtree.leaves
            self.transdict = curtree.transdict()
            self.translateblock = curtree.translateblock(self.transdict)
        elif curtree.leaves != self.leaves:
            msg = "Leaves on tree number {} are different than previous trees".format(self.tree_count)
            raise TreeError(msg)

        self.tree_count += 1
        self.tree_weight_sum += weight       # The weighted equivalent of tree_count
        bipdict = curtree.bipdict(interner=self.interner)

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

        return bipdict  # Experimental: so bigtreesummary can reuse and avoid topology computation

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

    def log_clade_credibility(self, topology):
        """Compute log clade credibility for topology (sum of log(freq) for all branches)"""

        bipartsummary = self.bipartsummary
        logsum = 0.0
        for bipartition in topology:
            logsum += math.log(bipartsummary[bipartition].freq)
        return logsum

    ###############################################################################################

    def contree(self, cutoff=0.5, allcompat=False, labeldigits=3):
        """Returns a consensus tree built from selected bipartitions"""

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

        return contree

###################################################################################################
###################################################################################################
###################################################################################################

class BigTreeSummary(TreeSummary):
    """Class summarizing bipartitions, branch lengths, and topologies from many trees"""

    # Does everything TreeSummary does and also keeps track of topologies
    # (topology list is potentially quite big, which is the reason for not including it in TS)

    def __init__(self, interner=None, store_trees=False):
        TreeSummary.__init__(self, interner=None)
        self._toposummary = {}
        self._toposummary_processed = False
        self.store_trees = store_trees

    ###############################################################################################

    @property
    def toposummary(self):
        """Property method for lazy evaluation of topostruct.freq"""
        if not self._toposummary_processed:
            for topostruct in self._toposummary.values():
                topostruct.freq = topostruct.weight / self.tree_weight_sum
            self._toposummary_processed = True

        return self._toposummary

    ###############################################################################################

    def add_tree(self,curtree, weight=1.0):
        """Add tree to treesummary, update all summaries"""

        self._toposummary_processed = False

        # Superclass method takes care of updating n_trees and all bipart-related info
        # Also returns bipdict so we wont have to recompute here
        bipdict = TreeSummary.add_tree(self, curtree, weight)

        # If topology has never been seen before, then add it and initialize count
        # If topology HAS been seen before then update count
        # Python: A bit messy that interning is here done in add_tree() and not in topology()
        # This is to avoid recomputing bipdict in topology function (or storing bipdictcache)
        topology = self.interner.intern_topology(frozenset(bipdict.keys()))
        if topology in self._toposummary:
            self._toposummary[topology].weight += weight
        else:
            self._toposummary[topology]=Topostruct()
            self._toposummary[topology].weight = weight
            if self.store_trees:
                self._toposummary[topology].tree = curtree

    ###############################################################################################

    def update(self, other):
        """Merge this object with other treesummary"""

        self._toposummary_processed = False

        # Superclass method takes care of updating:
        # tree_count, tree_weight_sum, and bipartsummary
        TreeSummary.update(self, other)

        # Merge other.toposummary with self.toposummary
        for topology in other._toposummary:
            # If topology already in self.toposummary, update count
            if topology in self._toposummary:
                self._toposummary[topology].weight += other._toposummary[topology].weight
            # If topology has never been seen before, simply transfer entry
            else:
                self._toposummary[topology]=other._toposummary[topology]

    ###############################################################################################

    def max_clade_cred_tree(self, labeldigits=3):
        """Find and return maximum clade credibility tree"""

        maxlogcred = -math.inf
        for topology in self.toposummary:
            logcred = self.log_clade_credibility(topology)
            if logcred > maxlogcred:
                maxlogcred = logcred
                maxlogcredtopo = topology

        maxcredbipdict = {}
        for bipartition in maxlogcredtopo:
            branch = self.bipartsummary[bipartition]
            branch.label = f"{round(branch.freq, labeldigits)}"
            maxcredbipdict[bipartition] = branch

        # Build tree from bipartitions  in new bipdict
        maxcredtree = Tree.from_biplist(maxcredbipdict)
        return maxcredtree, maxlogcred

###################################################################################################
###################################################################################################
###################################################################################################


class Treefile():
    """Abstract base-class for representing tree file objects."""

    # Classes for specific formats inherit from this class and add extra stuff as needed.
    # NOTE: i am opening files in "read text" with encoding UTF-8. Will this work across platforms?

    def __init__(self, filename=None, filecontent=None):

        num_args = (filename is not None) + (filecontent is not None)
        if num_args != 1:
            raise TreeError("Treefile __init__ requires either filename or filecontent (not both)")
        elif filecontent:
            self.treefile = StringIO(filecontent)
        else:
            self.treefile = open(filename, mode="rt", encoding="UTF-8")

        self.buffer = ""                # Used for keeping leftovers after reading whole line

    ###############################################################################################

    def __enter__(self):
        """Implements context manager behaviour for Treefile types.
        Usage example:
            with pt.Newicktreefile(filename) as tf:
                mytree = tf.readtree()
            mytree.rootminvar()
        """
        return self

    ###############################################################################################

    def __exit__(self, type, value, traceback):
        """Implements context manager behaviour for Treefile types.
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

        return treestring

    ###############################################################################################

    def readtree(self):
        """Reads one tree from file and returns as Tree object. Returns None when exhausted file"""

        try:
            tree = next(self)
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
###################################################################################################


class Newicktreefile(Treefile):
    """Class representing Newick tree file. Iteration returns tree-objects"""

    def __init__(self, filename=None, filecontent=None):
        Treefile.__init__(self, filename, filecontent)
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

    def __next__(self):
        treestring = self.get_treestring()
        if treestring is None:
            self.treefile.close()
            raise StopIteration
        else:
            treestring =  remove_comments(treestring, leftdelim = "[", rightdelim="]")
            return Tree.from_string(treestring)

###################################################################################################
###################################################################################################
###################################################################################################


class Nexustreefile(Treefile):
    """Class representing Nexus tree file. Iteration returns tree object or None"""

    ###############################################################################################

    def __init__(self, filename=None, filecontent=None):
        """Read past NEXUS file header, parse translate block if present"""

        Treefile.__init__(self, filename, filecontent)

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
            return remove_comments(line, leftdelim = "[", rightdelim="]")

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
        self.buffer = read_until(self.tree_header_pattern)

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

    def __next__(self, noreturn=False):

        treestring = self.get_treestring()
        if treestring is None:
            self.treefile.close()
            raise StopIteration

        # remove comments in brackets if present
        # remove leading "tree NAME =" (compiled regexp "tree_header_pattern")
        # Implementation note: does not deal with figtree comments!
        treestring = remove_comments(treestring, leftdelim = "[", rightdelim="]")
        treestring = self.tree_header_pattern.sub("", treestring)

        # If "end;" statement has been reached: terminate for loop, do NOT return tree object
        if self.end_pattern.search(treestring):
            self.treefile.close()
            raise StopIteration

        # Return tree object if requested
        if noreturn:
            return None
        else:
            return Tree.from_string(treestring, self.transdict)

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


# # Placeholder: Insert test code here and run module in standalone mode
def main():
    pass

###################################################################################################


if __name__ == "__main__":
    main()
